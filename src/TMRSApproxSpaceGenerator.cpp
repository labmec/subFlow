//
//  TMRSApproxSpaceGenerator.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#include "TMRSApproxSpaceGenerator.h"
#include "TPZMHMixedMeshWithTransportControl.h"
#include "TPZCompMeshTools.h"
#include "TPZMixedDarcyWithFourSpaces.h"
#include "TPZMHMixedMesh4SpacesControl.h"
#include "TPZFastCondensedElement.h"
#include "TPZReservoirTools.h"
#include "TPZPostProcessResProp.h"
#include "TPZDarcyMemory.h"
#include "TPZDarcyFlowWithMem.h"
#include "TMRSDarcyFractureFlowWithMem.h"
#include "TMRSDarcyFractureGlueFlowWithMem.h"
#include "TPZLagrangeMultiplierCS.h"
#include "TPZCompElHDivCollapsed.h"
#include "TMRSDarcyMemory.h"
#include "TMRSTransportMemory.h"
#include "TPZNullMaterialCS.h"
#include "pzcompelwithmem.h"
#include "pzsmanal.h"
#include "TPZRefPatternTools.h"
#include "TPZRefPattern.h"
#include <pzshapequad.h>
#include <pzshapelinear.h>
#include "pzvec_extras.h"
#include "TPZSimpleTimer.h"
#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger logger("imrs");
#endif

using namespace std;

void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSApproxSpaceGenerator::TMRSApproxSpaceGenerator()
{
    mGeometry = nullptr;
    mMixedOperator = nullptr;
    mTransportOperator = nullptr;
    mHybridizer = nullptr;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSApproxSpaceGenerator &TMRSApproxSpaceGenerator::operator=(const TMRSApproxSpaceGenerator &other){
    DebugStop();
    return *this;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSApproxSpaceGenerator::~TMRSApproxSpaceGenerator(){
    if (mHybridizer){
        delete mHybridizer;
    }
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::Read(TPZStream &buf, void *context){
    DebugStop();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

int TMRSApproxSpaceGenerator::ClassId() const{
    DebugStop();
    return -1;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetGeometry(TPZGeoMesh * geometry){
    mGeometry = geometry;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetGeometry(TPZGeoMesh * gmeshfine,TPZGeoMesh * gmeshcoarse){
    if (InitMatIdForMergeMeshes() == -1000) {
        DebugStop(); // MatId for MHM should be set in main!
    }
    bool useMHM = mSimData.mTNumerics.m_mhm_mixed_Q;
    const bool needsMergeMeshes = mSimData.mTNumerics.m_need_merge_meshes_Q;
    if(needsMergeMeshes){
        MergeMeshes(gmeshfine, gmeshcoarse); // this fills mSubdomainIndexGel that is used for MHM
    }
    
    mGeometry = gmeshfine;
#ifdef PZDEBUG
    {
        int64_t nelem = mGeometry->NElements();
        for (int64_t el = 0; el<nelem; el++) {
            TPZGeoEl *gel = mGeometry->Element(el);
            if(gel->Dimension() != mGeometry->Dimension()) continue;
            if(gel->HasSubElement()) continue;
            int geldomain =-1;
            if(useMHM){
                geldomain = mSubdomainIndexGel[el];
                if(geldomain == -1) DebugStop();
            }
           
            int firstside = gel->FirstSide(2);
            int lastside = gel->NSides()-1;
            for(int side = firstside; side < lastside; side++)
            {
                TPZGeoElSide gelside(gel,side);
                for (auto neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
                    int neighdim = neigh.Element()->Dimension();
                    if(neighdim != 3) continue;
                    int neighind = neigh.Element()->Index();
                    
                    if (useMHM) {
                        int neighdomain = mSubdomainIndexGel[neighind];
                        if(neighdomain != geldomain)
                        {
                            if(!neigh.HasNeighbour(this->mSimData.mTGeometry.m_skeletonMatId))
                            {
                                DebugStop();
                            }
                        }
                    }
                   
                }
            }
        }
    }
#endif
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::ApplyUniformRefinement(int nelref){
    
    
    for (int iref=0; iref<nelref; iref++) {
        int nel = mGeometry->NElements();
        for (int iel = 0; iel <nel; iel++) {
            TPZGeoEl *gel = mGeometry->Element(iel);
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            TPZVec<TPZGeoEl *> sons;
            gel->Divide(sons);
        }
    }
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::ApplyUniformRefinement(){
    std::cout << "Applying uniform refinement numref = " << mSimData.mTGeometry.mnref << "\n";
    
    ApplyUniformRefinement(mSimData.mTGeometry.mnref);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::PrintGeometry(std::string name, bool vtkFile, bool textfile)
{
    if (!mGeometry) {
        DebugStop();
    }
    std::stringstream text_name;
    std::stringstream vtk_name;
    if (vtkFile) {
        vtk_name   << name << "_geometry"  << ".vtk";
        std::ofstream vtkfile(vtk_name.str().c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, vtkfile, true);
    }
    if (textfile) {
        text_name  << name << "_geometry" << ".txt";
        std::ofstream textfile(text_name.str().c_str());
        mGeometry->Print(textfile);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::AddAtomicMaterials(const int dim, TPZCompMesh* cmesh,
                                                  std::set<int>& matids,
                                                  std::set<int>& bcmatids,
                                                  const bool isInsertBCs) {
    GetMaterialIds(dim, matids, bcmatids);
    
    // Domain
    for (auto it:matids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dim);
        cmesh->InsertMaterialObject(nullmat);
    }
    
    // BCs
    if(isInsertBCs){
        for (auto it:bcmatids) {
            TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dim-1);
            cmesh->InsertMaterialObject(nullmat);
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SplitConnectsAtInterface(TPZCompElSide& compside) {
    
    TPZCompMesh* fluxmesh = compside.Element()->Mesh();
    const int gmeshdim = fluxmesh->Reference()->Dimension();
    
    // ===> Find 3D element neighbor
    TPZCompElSide compsideright;
    TPZGeoElSide gleft(compside.Reference());
    TPZGeoElSide neigh = gleft.Neighbour();
    int icon = 0;
    while(neigh != gleft){
        TPZGeoEl* gelnei = neigh.Element();
        if (!gelnei || gelnei->Dimension() != gmeshdim) {
            neigh++;
            continue;
        }
        compsideright = TPZCompElSide(gelnei->Reference(), neigh.Side());
        break;
    }
    if (!compsideright.Element()) {
        DebugStop(); // Could not find neighbor element!
    }
    
    compside.SplitConnect(compsideright);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// return the connect index of the element referenced by geoelside
static int64_t GeoElSideConnectIndex(const TPZGeoElSide &gelside)
{
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(gelside.Element()->Reference());
    if(!intel) DebugStop();
    if(intel->NSideConnects(gelside.Side()) != 1) DebugStop();
    return intel->SideConnectIndex(0, gelside.Side());
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh){
    bool useMHM = mSimData.mTNumerics.m_mhm_mixed_Q;
    // ===> Create the fracture hdivcollapsed elements
    // These should be unconnected with the 3D elements so we reset references in the gmesh
    // They are, however, connected (with connects) between themselves
    TPZGeoMesh* gmesh = cmesh->Reference();
    gmesh->ResetReference(); // So it does not try to use connects of neighbor elements
    gmesh->SetReference(cmesh);
    const int gmeshdim = gmesh->Dimension();
    cmesh->SetDimModel(gmeshdim-1);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(gmeshdim-1);
    // organize the elements by matid
    std::map<int,TPZStack<int64_t>> matTogelindex;
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el < nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        const int hassubel = gel->HasSubElement();
        if (hassubel) { // the mesh can be uniformly refined
            continue;
        }
        const int matid = gel->MaterialId();
        matTogelindex[matid].Push(el);
    }
    int nfrac = mSimData.mTFracProperties.m_fracprops.size();
    for(auto frac : mSimData.mTFracProperties.m_fracprops)
    {
        int fracmatid = frac.first;
        std::set<int>& frabcids = frac.second.m_fracbc;
//        int fracbcid = frac.second.m_fracbc;
        int fracintersect = frac.second.m_fracIntersectMatID;
        for(auto gelindex : matTogelindex[fracmatid])
        {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            TPZInterpolationSpace* hdivcollapsed = nullptr;
            if (gmeshdim == 2){
                hdivcollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>(*cmesh,gel);
            }
            else {
                if(gel->Type() == MElementType::ETriangle){
                    hdivcollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeTriang>(*cmesh,gel);
                }
                else if (gel->Type() == MElementType::EQuadrilateral){
                    hdivcollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,gel);
                }
                else{
                    DebugStop();
                }
            }
        }
		
		// Here we are creating the boundary compels of the fracture
        TPZStack<int64_t> allbcindexes;
        for(const int bcid : frabcids){
            for(int64_t& ind : matTogelindex[bcid]){
                allbcindexes.Push(ind);
            }
        }
        cmesh->AutoBuild(allbcindexes);
		
		// Here, we fix the side orient of neighboring fracture elements
		for(auto gelindex : matTogelindex[fracmatid]) {
			TPZGeoEl *gel = gmesh->Element(gelindex);
			TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(gel->Reference());
			if(!intel) DebugStop();
			for (int iside = gel->FirstSide(1); iside < gel->FirstSide(2); iside++) {
				TPZGeoElSide gelside(gel,iside);
                // this the side orient of the current element
				const int gelSideOrient = intel->GetSideOrient(iside);
#ifdef PZ_LOG
                if (logger.isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "\ninitial side orient\n";
                    sout << "gel index " << gel->Index() << " side " << iside << " side orient " << gelSideOrient << std::endl;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
				for(TPZGeoElSide neigh = gelside.Neighbour() ; neigh != gelside ; neigh++){
					TPZGeoEl* neighgel = neigh.Element();
					if(neighgel->HasSubElement()) continue;
                    if(neighgel->Dimension() > gmeshdim-1) continue;
                    // if the neighbour matid is neither a fracture or a fracture boundary no action
                    int neighmatid = neighgel->MaterialId();
                    // we don;t adjust if the neighbour is not a fracture or "my" boundary
//                    auto neighisfrac = IsFracMatId(neighmatid);
                    
                    bool neighisfrac=mSimData.mTFracProperties.m_fracprops.find(neighmatid) != mSimData.mTFracProperties.m_fracprops.end();
                    bool neighisbound= false;
                    int intersecMatId=mSimData.mTFracProperties.m_fracprops[fracmatid].m_fracIntersectMatID;
                    for (auto bc: mSimData.mTFracProperties.m_fracprops[fracmatid].m_fracbc) {
                        if(bc==neighmatid){
                            neighisbound=true;
                            break;
                        }
                    }

                    if(!neighisfrac && !neighisbound) {
                        continue;
                    }
					// we dont adjust if the neighmatid is different and if there is not intersection element
                    TPZGeoElSide gelintersect = gelside.HasNeighbour(intersecMatId);
                    
                    if(neighisbound && gelintersect)
                    {
                        DebugStop();
                    }
                    // the neighbour is a different fracture but there is no intersection, we dont adjust
                    if(neighisfrac && neighmatid != fracmatid && !gelintersect) {
                        continue;
                    }
                    
                    if(gelintersect) {
                        intel->SetSideOrient(iside, 1);
                        {
                            int sideorient = intel->GetSideOrient(iside);
#ifdef PZ_LOG
                            if (logger.isDebugEnabled()) {
                                std::stringstream sout;
                                sout << "\nAdjusting orientation because of intersection\n";
                                sout << "gel index " << gel->Index() << " side " << iside << " side orient " << sideorient << std::endl;
                                LOGPZ_DEBUG(logger, sout.str())
                            }
#endif
                        }
                        break;
                    }
					TPZInterpolatedElement* neighintel = dynamic_cast<TPZInterpolatedElement*>(neighgel->Reference());
					if(!neighintel) continue;
					const int neighside = neigh.Side();
					
                    // we will hybridize, therefore all sideorients should be 1
                    // neighintel->NConnects()==1 means neighintel is typy HDivbound
					if(neighisbound){
						intel->SetSideOrient(iside, 1);
						neighintel->SetSideOrient(neighside, 1);
                        int sideorient = intel->GetSideOrient(iside);
                        int neighsideorient = neighintel->GetSideOrient(neighside);
#ifdef PZ_LOG
                            if (logger.isDebugEnabled()) {
                                std::stringstream sout;
                                sout << "\nAdjusting orientation of boundary\n";
                                sout << "gel index " << gel->Index() << " side " << iside << " side orient " << sideorient << std::endl;
                                sout << "gel index " << neighintel->Reference()->Index() << " side " << neighside << " side orient " << neighsideorient << std::endl;
                                LOGPZ_DEBUG(logger, sout.str())
                            }
#endif
						continue;
					}
                    // here we assume there is no fracture intersection
                    // if there is a computational element, it has to be of the same fracmatid
                    // no other elements have been loaded into geomesh
					const int neighSideOrient = neighintel->GetSideOrient(neighside);
					
					if(gelSideOrient*neighSideOrient > 0){
#ifdef PZ_LOG
                            if (logger.isDebugEnabled()) {
                                std::stringstream sout;
                                sout << "\nAdjusting wrong neighbouring orientations gelSideOrient " << gelSideOrient << " neighSideOrient " << neighSideOrient << std::endl;
                                LOGPZ_DEBUG(logger, sout.str())
                            }
#endif
						const int newNeighSideOrient = -gelSideOrient;
						neighintel->SetSideOrient(neighside, newNeighSideOrient);
                        {
                            int sideorient = intel->GetSideOrient(iside);
                            int neighsideorient = neighintel->GetSideOrient(neighside);
#ifdef PZ_LOG
                            if (logger.isDebugEnabled()) {
                                std::stringstream sout;
                                sout << "\ngel index " << gel->Index() << " side " << iside << " side orient " << sideorient << std::endl;
                                sout << "gel index " << neighintel->Reference()->Index() << " side " << neighside << " side orient " << neighsideorient << std::endl;
                                LOGPZ_DEBUG(logger, sout.str())
                            }
#endif
                        }
					}
				}
			}
		}
		
		// Now, we split the connects of sides of fracture elements that have intersections
        for(auto gelindex : matTogelindex[fracintersect])
        {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            // SplitConnect will put an element with fracintersect matid as first neighbour of the fracture elements
            // The only elements visible to the geometric mesh are the elements of the current frac matid
            TPZReservoirTools::SplitConnect(gel);
            int64_t nel = gmesh->NElements();
            mSubdomainIndexGel.Resize(nel, -1);
            // create a geometric element for interface elements
            int interfacematid = mSimData.mTGeometry.mInterface_material_id;

            std::set<int> domains;
            TPZGeoElSide gelside(gel);
            TPZGeoElSide previous(gelside);
            // this will look for a geoelside that points to me
            previous--;
            if(previous.Element()->MaterialId() != fracmatid) DebugStop();
            if(previous.Neighbour() != gelside) DebugStop();
            // setting the domain of the newly created elements
            // creating an interface geometric element
            if (useMHM) {
                int domain = mSubdomainIndexGel[previous.Element()->Index()];
                if(domain == -1) DebugStop();
                mSubdomainIndexGel[gel->Index()] = domain;
                domains.insert(domain);
                TPZGeoElBC gbc(gelside,interfacematid);
                nel = mGeometry->NElements();
                mSubdomainIndexGel.Resize(nel, -1);
                mSubdomainIndexGel[gbc.CreatedElement()->Index()] = domain;
            }
            
            // look for the other element with the fractintersect matid
            // this element was created in splitconnect
            for(auto neigh = gelside.Neighbour(); neigh != gelside; neigh++)
            {
                int domain =-1;
                int neighmatid = neigh.Element()->MaterialId();
                if (useMHM) {
                    domain = mSubdomainIndexGel[neigh.Element()->Index()];
                    // here we check for any fracture matid
                    if(IsFracMatId(neighmatid)) {
                        domains.insert(domain);
                        
                    };
                }
                
                // here we check for the current fracture being processed
                if(neighmatid == fracmatid)
                {
                    // the fracture element must belong to some subdomain
                   
                    if(domain == -1 && useMHM) DebugStop();
                    auto bound = neigh.Neighbour();
                    // the first neighbour of the fracture element needs to be a fracintersect
                    if(bound.Element()->MaterialId() != fracintersect) DebugStop();
                    
                    TPZGeoElBC gbc(bound,interfacematid);
                    if (useMHM) {
                        mSubdomainIndexGel[bound.Element()->Index()] = domain;
                        nel = mGeometry->NElements();
                        mSubdomainIndexGel.Resize(nel, -1);
                        mSubdomainIndexGel[gbc.CreatedElement()->Index()] = domain;
                        // the pressure geometric element was created, if needed, in split connect??
                        auto pressure = neigh.HasNeighbour(mSimData.mTGeometry.m_pressureMatId);
                        if(!pressure) DebugStop();
                        mSubdomainIndexGel[pressure.Element()->Index()] = domain;
                    }
                    
                }
            }
            // if the neighbouring elements belong to different macro domains then the pressure will stay in the father mesh
            if(domains.size()>1 && useMHM)
            {
                auto pressure = gelside.HasNeighbour(mSimData.mTGeometry.m_pressureMatId);
                if(!pressure) DebugStop();
                mSubdomainIndexGel[pressure.Element()->Index()] = -1;
            }
        }

        // reset the references of the elements of the fracture
        for(auto gelindex : matTogelindex[fracmatid])
        {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            if(!gel->Reference()) DebugStop();
            gel->ResetReference();
        }
        for(auto gelindex : allbcindexes)
        {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            if(!gel->Reference()) DebugStop();
            gel->ResetReference();
        }
        for(auto gelindex : matTogelindex[fracintersect])
        {
            TPZGeoEl *gel = gmesh->Element(gelindex);
            if(!gel->Reference()) DebugStop();
            gel->ResetReference();
            TPZGeoElSide gelside(gel);
            // look for the other fracintersect elements we have created
            for(TPZGeoElSide neigh = gelside.Neighbour(); neigh != gelside; neigh++)
            {
                if(neigh.Element()->MaterialId() == fracintersect)
                {
                    if(!neigh.Element()->Reference()) DebugStop();
                    neigh.Element()->ResetReference();
                }
            }
        }
    }
    
    // ===> Set HdivCollapsed elements connection with adjacent 3d elements
    // This is done through top and bottom connects of the hdivcollapsed el.
    // Since bottom connect's default sideorient is negative one (-1),
    // we should connect it with the 3D element side that has sideorient positive one (+1)
    cmesh->LoadReferences(); // So the hdivcollapsed element can see the 3D the flux mesh adjacent elements to it
    
    // if there are superposed fracture elements then first connects need to be created to link the superposed
    // fracture elements
    // then the extremes need to link to the 3D elements
    for(auto cel : cmesh->ElementVec()) {
        if (!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        const int matid = gel->MaterialId();
        //if (matid != FractureMatId()) continue;
        if (!IsFracMatId(matid)) continue;
        int64_t index;
        TPZInterpolatedElement* hdivcollapsed = nullptr;
        if (gmeshdim == 2) {
            hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>* >(cel);
        }
        else{
            if(gel->Type() == MElementType::ETriangle){
                hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang>* >(cel);
            }
            else if (gel->Type() == MElementType::EQuadrilateral){
                hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* >(cel);
            }
        }
        
        if (!hdivcollapsed) {
            DebugStop();
        }
        
        int nconnects = hdivcollapsed->NConnects();
        // for overlapping fractures the bottom and top connects are initialized all at once
        try {
            auto bottomconnectindex = hdivcollapsed->ConnectIndex(nconnects);
        } catch (...) {
            PZError << "This warning was not a problem since I am at a try catch";
            continue;
        }
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neigh = gelside.Neighbour();
        // left and right volumetric elements
        // left has side orient 1 , right has side orient -1
        std::pair<TPZCompElSide,TPZCompElSide> leftright;
        std::pair<int,int> leftrightdomain;
        std::pair<int64_t,int64_t> leftrightcindex;
        // hdivdomain : domain of the fracture element
        int hdivdomain =-1;
        if(useMHM){
            hdivdomain = mSubdomainIndexGel[gel->Index()];
        }
        
        int icon = 0;
        // loop over all neighbours of the fracture element
        while(neigh != gelside){
            TPZGeoEl* gelnei = neigh.Element();
            if (!gelnei || gelnei->Dimension() != gmeshdim) {
                neigh++;
                continue;
            }
            // we found a volumetric element
            TPZCompElSide compside = TPZCompElSide(gelnei->Reference(), neigh.Side());
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(compside.Element());
            const int sideorient = intel->GetSideOrient(compside.Side());

            // split the connect shared by two adjacent 3d elements if it is the first time
            if(icon==0) SplitConnectsAtInterface(compside);
            if (sideorient == 1){
                leftright.first = compside;
                if(useMHM){
                  leftrightdomain.first = mSubdomainIndexGel[gelnei->Index()];
                }
               
                leftrightcindex.first = compside.ConnectIndex();
                icon++;
            } else if(sideorient == -1){
                leftright.second = compside;
                if(useMHM){
                leftrightdomain.second = mSubdomainIndexGel[gelnei->Index()];
                }
                leftrightcindex.second = compside.ConnectIndex();
                icon++;
            }
            neigh++;
        }
        // this would imply more than 2 3D elements linked to the fracture
        if(icon != 2) DebugStop();
        std::pair<TPZConnect,TPZConnect> cleftright = {cmesh->ConnectVec()[leftrightcindex.first],cmesh->ConnectVec()[leftrightcindex.second]};
        if(cleftright.first.HasDependency() && leftrightdomain.first == leftrightdomain.second && useMHM) DebugStop();
        if(cleftright.second.HasDependency() && leftrightdomain.first == leftrightdomain.second && useMHM) DebugStop();
        // verify if we need to swap the connects
		// NOTE Jul 2022: This code was need to run case 2 of flemisch benchmark with mhm
		//                It seems, it is not needed anymore since it works without it. If something happens, please check.
//        bool needswap = false;
//
//        if((hdivdomain == leftrightdomain.first && cleftright.first.HasDependency()) ||
//           (hdivdomain == leftrightdomain.second && cleftright.second.HasDependency()))
//        {
//            needswap = true;
//        }
//        if(!cleftright.first.HasDependency() && !cleftright.second.HasDependency())
//        {
//            auto skel = gelside.HasNeighbour(18);
//            if(!skel) DebugStop();
//            TPZCompEl *cskel = skel.Element()->Reference();
//            if(!cskel) DebugStop();
//            int skelcindex = cskel->ConnectIndex(0);
//            if(hdivdomain == leftrightdomain.first && leftrightcindex.first == skelcindex) needswap = true;
//            if(hdivdomain == leftrightdomain.second && leftrightcindex.second == skelcindex) needswap = true;
//        }
//        if(needswap)
//        {
            TPZInterpolatedElement *intelL = dynamic_cast<TPZInterpolatedElement*>(leftright.first.Element());
            TPZInterpolatedElement *intelR = dynamic_cast<TPZInterpolatedElement*>(leftright.second.Element());
            int Llocindex = intelL->SideConnectLocId(0, leftright.first.Side());
            int Rlocindex = intelR->SideConnectLocId(0, leftright.second.Side());
            intelL->SetConnectIndex(Llocindex, leftrightcindex.second);
            intelR->SetConnectIndex(Rlocindex, leftrightcindex.first);
            leftrightcindex.first = leftright.first.ConnectIndex();
            leftrightcindex.second = leftright.second.ConnectIndex();
//        }
        // verify if the order is OK
        TPZConnect &cL = cmesh->ConnectVec()[leftrightcindex.first];
        // verify if we need to swap the connects
        if(hdivdomain == leftrightdomain.first && cL.HasDependency() && useMHM)
        {
            DebugStop();
        }
        TPZConnect &cR = cmesh->ConnectVec()[leftrightcindex.second];
        // verify if we need to swap the connects
        if(hdivdomain == leftrightdomain.second && cR.HasDependency() && useMHM)
        {
            DebugStop();
        }
        int fracgluematid = mSimData.mTFracIntersectProperties.m_FractureGlueId;
        if((fracgluematid != -10000) && gelside.HasNeighbour(fracgluematid))
        {
            // there are overlapping fracture elements, they have been sorted
            TPZGeoElSide first3D(leftright.first.Reference());
            TPZGeoElSide last3D(leftright.second.Reference());
            int64_t bottomcindex = leftrightcindex.first;
            neigh = first3D.Neighbour();
            // create a new connect for the HDivBound elements
            while(neigh != last3D)
            {
                TPZCompEl *cel = neigh.Element()->Reference();
                int nc = cel->NConnects();
                if(nc == 1)
                {
                    TPZConnect &c = cel->Connect(0);
                    int64_t cindex = cmesh->AllocateNewConnect(c.NShape(), c.NState(), c.Order());
                    cel->SetConnectIndex(0, cindex);
                }
                neigh++;
            }
            TPZGeoElSide prev = first3D;
            neigh = prev.Neighbour();
            TPZGeoElSide next = neigh.Neighbour();
            // establish the connect indices for the HDivCollapsed elements
            while(neigh != last3D)
            {
                TPZCompEl *cel = neigh.Element()->Reference();
                int nc = cel->NConnects();
                if(nc != 1)
                {
                    int64_t cindexprev = GeoElSideConnectIndex(prev);
                    int64_t cindexnext = GeoElSideConnectIndex(next);
                    cel->SetConnectIndex(nc, cindexprev);
                    cel->SetConnectIndex(nc+1, cindexnext);
                }
                prev++;
                neigh++;
                next++;
            }
        } else
        {
            hdivcollapsed->SetConnectIndex(nconnects, leftrightcindex.first);
            hdivcollapsed->SetConnectIndex(nconnects+1, leftrightcindex.second);
        }
        cleftright.first.RemoveDepend();
        cleftright.second.RemoveDepend();
    }
    
    cmesh->ExpandSolution();
//    std::ofstream out("cmeshaftercollapsed.txt");
//    cmesh->Print(out);
    
    cmesh->SetDimModel(gmeshdim);
#ifdef PZDEBUG
//    {
//        cmesh->ComputeNodElCon();
//        int64_t ncon = cmesh->ConnectVec().NElements();
//        for(int64_t ic = 0; ic<ncon; ic++)
//        {
//            TPZConnect &c = cmesh->ConnectVec()[ic];
//            if(c.NElConnected() > 2 && c.NShape() != 4) DebugStop();
//        }
//    }
#endif
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::AdjustOrientBoundaryEls(TPZCompMesh* cmesh, std::set<int>& buildmatids) {
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(buildmatids.find(matid) == buildmatids.end()) continue;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        intel->SetSideOrient(gel->NSides()-1, 1);
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour(gelside.Neighbour());
        for( ;neighbour != gelside; neighbour++)
        {
            TPZGeoEl *neighgel = neighbour.Element();
            TPZCompEl *neighcel = neighgel->Reference();
            if(!neighcel) continue;
            int nc = neighcel->NConnects();
            if(nc>1 && neighgel->Dimension() == gel->Dimension()+1)
            {
                TPZInterpolationSpace *neighintel = dynamic_cast<TPZInterpolationSpace *>(neighcel);
                neighintel->SetSideOrient(neighbour.Side(), 1);
            }
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateFractureHDivCompMesh(TPZCompMesh* cmesh,
                                                          std::set<int>& matids, std::set<int>& bcids,
                                                          std::set<int>& matids_dim2, std::set<int>& bcids_dim2){
    
    TPZGeoMesh* gmesh = cmesh->Reference();
    const int dim = gmesh->Dimension();
    
    // ===> Creating space for 3D elements
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false); // we need to disconnect by hand at the fracture location later
    std::set<int> _3dmatids(matids);
    _3dmatids.insert(bcids.begin(),bcids.end());
    cmesh->AutoBuild(_3dmatids);
    
#ifdef PZDEBUG
    {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() != dim) continue;
            int firstside = gel->FirstSide(2);
            for(int is = firstside; is< gel->NSides()-1; is++)
            {
                TPZGeoElSide gelside(gel,is);
                auto neigh = gelside.Neighbour();
                if(neigh == gelside) DebugStop();
                while(neigh != gelside){
                    if(neigh.Element()->Dimension() == dim) break;
                    neigh = neigh.Neighbour();
                }
                if(neigh.Dimension() != dim) continue;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                int sideorient = intel->GetSideOrient(is);
                TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(neigh.Element()->Reference());
                int sideorientneigh = intelneigh->GetSideOrient(neigh.Side());
                if(sideorient*sideorientneigh != -1) DebugStop();
            }
        }
//        std::ofstream out("HdivMesh.txt");
//        cmesh->Print(out);
    }
#endif
    // ===> Creating fracture element
    // properly order overlapping fracture elements
    // create fracture glue elements between the fractures
    OrderOverlappingFractures();
    // create HDivBound elements corresponding to the fracture glue elements
    _3dmatids.clear();
    
    // split the volumetric elements in the flux mesh
    CreateFractureHDivCollapsedEl(cmesh);
    cmesh->CleanUpUnconnectedNodes();
    
    if(0)
    {
        std::ofstream out("EuNaoSouLouco.txt");
        cmesh->Print(out);
    }
    // ===> Create BCs for fracture
    // make sure all boundary elements of fractures have positive sideorient
    cmesh->LoadReferences();
//    AdjustOrientBoundaryEls(cmesh,bcids_dim2);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(3);
	
#ifdef PZDEBUG
//    {
//        std::ofstream out("gmesh.txt");
//        gmesh->Print(out);
//    }
#endif
    
    // ===> Fix blocks
    cmesh->InitializeBlock();
    
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * TMRSApproxSpaceGenerator::HdivFluxCmesh(int order){
    
    const bool hasFrac = isFracSim();
    
    // -----------> Problem dimension
    if (!mGeometry)
        DebugStop();
    const int dim = mGeometry->Dimension();
    
    // -----------> Creating hdiv comp mesh
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    cmesh->SetDefaultOrder(order);
    
    // -----------> Inserting atomic matrix materials
    std::set<int> matids, bcids;
//    bcids.insert(mSimData.mTGeometry.m_skeletonMatId);
    if(this->mSimData.mTNumerics.m_mhm_mixed_Q)
    {
        bcids.insert(mSimData.mTGeometry.m_skeletonMatId-1);
    }
    int gluematid = mSimData.mTFracIntersectProperties.m_FractureGlueId;
    if(gluematid > 0) bcids.insert(gluematid);
    AddAtomicMaterials(dim,cmesh,matids,bcids);


    // -----------> Setting space and building mesh
    if (hasFrac) {
        std::set<int> matids_dim2, bcids_dim2;
        const int dimfrac = dim-1;
        AddAtomicMaterials(dimfrac,cmesh,matids_dim2,bcids_dim2); // Inserting atomic fracture materials
        CreateFractureHDivCompMesh(cmesh,matids,bcids,matids_dim2,bcids_dim2);
    }
    else{
        // simple autobuild suffice
        cmesh->SetAllCreateFunctionsHDiv();
        cmesh->AutoBuild();
        cmesh->InitializeBlock();
    }

#ifdef PZDEBUG
    CheckSideOrientOfFractureEls();
#endif
    
#ifdef PZDEBUG
//    std::ofstream sout("q_mesh.txt");
//    cmesh->Print(sout);
#endif
#ifdef PZDEBUG
    // This method checks the number of elements for each matid
//    {
//        TPZGeoMesh *gmesh = mGeometry;
//        int64_t nel = gmesh->NElements();
//        std::map<int,int> numels;
//        for (int64_t el = 0; el<nel; el++) {
//            TPZGeoEl *gel = gmesh->Element(el);
//            int matid = gel->MaterialId();
//            numels[matid]++;
//        }
//        std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
//        for(auto it: numels){
//            std::cout << "For matid " << it.first << " number of elements " << it.second << std::endl;
//        }
//    }
#endif

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CheckSideOrientOfFractureEls() {
    if(mSimData.mTNumerics.m_SpaceType != TMRSDataTransfer::TNumerics::E4Space) {
        // Hybridization is made differently in other space types. Thus, this check does not make sense unless it is E4Space
        return;
    }
    TPZGeoMesh* gmesh = mGeometry;
    for(TPZGeoEl* gel : gmesh->ElementVec()){
        const int gelmatid = gel->MaterialId();
        const int firstedge = gel->FirstSide(1);
        const int lastedge = gel->FirstSide(2);
        
        if(!mSimData.mTFracProperties.isFracMatId(gelmatid)){
            continue;
        }
        
        auto fracprop =  mSimData.mTFracProperties.m_fracprops[gelmatid];
        const int matidintersect = fracprop.m_fracIntersectMatID;
        std::set<int>& bcids = fracprop.m_fracbc;
        
        TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(gel->Reference());
        for (int is = firstedge; is < lastedge; is++) {
            const int sideorientgel = intel->GetSideOrient(is);
            TPZGeoElSide gside(gel,is);
            TPZGeoElSide neig = gside.Neighbour();
            bool isIntersect = false;
            for(; gside != neig ; neig++){
                TPZGeoEl* neiggel = neig.Element();
                const int neigmatid = neiggel->MaterialId();
                if (neigmatid == matidintersect) {
                    if(neiggel->Dimension() != 1) DebugStop();
                    // there is a neighbour that is an intersection
                    // note: we assume intersections are neighbors that are always before in the list of neighbors
                    isIntersect = true;
                }
                if(bcids.find(neigmatid) != bcids.end()){
                    if(neiggel->Dimension() != 1) DebugStop();
                    // This edge is at at bc of this fracture, then its sideorient should be positive
                    if(sideorientgel != 1) DebugStop();
                }
                if (neiggel->MaterialId() == gelmatid) {
                    if(neiggel->Dimension() != 2) DebugStop();
                    TPZInterpolatedElement* intelneig = dynamic_cast<TPZInterpolatedElement*>(neiggel->Reference());
                    const int sideorientneig = intelneig->GetSideOrient(neig.Side());
                    if(isIntersect){
                        if(sideorientgel != 1 || sideorientneig != 1) DebugStop();
                    }
                    else{
                        if(sideorientgel*sideorientneig != -1) DebugStop();
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

#include "pzshapequad.h"

/// create an HDiv mesh for mortar approximation
TPZCompMesh *TMRSApproxSpaceGenerator::HDivMortarFluxCmesh(char fluxmortarlagrange)
{
    std::set<int> matids, bcmatids;
    int dimension = mGeometry->Dimension();
    int nstate = 1;
    GetMaterialIds(dimension, matids, bcmatids);
    //    bcmatids.insert(-11);
    
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    cmesh->SetName("FluxMortarMesh");
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    for (auto it:matids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dimension,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    for (auto it:bcmatids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dimension-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId,dimension-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    
    
    std::set<int> buildmatids(matids);
    // create all flux elements as discontinuous elements
    cmesh->AutoBuild(buildmatids);
  {
    std::ofstream myfile("FluxMeshGood.txt");
    cmesh->Print(myfile);
  }
    
//  {
//    std::ofstream out("cmeshbef.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
//  }

    
    // add the Hdiv wrapper elements as boundary of the volumetric elements
    mGeometry->ResetReference();
    {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(mSimData.mTGeometry.m_HdivWrapMatId,dimension-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
#ifdef PZDEBUG
    {
        {
            std::ofstream out("FluxMortarMesh.txt");
            //            cmesh->Print(out);
        }
    }
#endif
    int64_t nelflux = cmesh->NElements();
    for (int64_t el = 0; el<nelflux; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != dimension) continue;
        int firstside = 0;
        for(int dim = 0; dim < dimension-1; dim++) firstside+=gel->NSides(dim);
        for (int side=firstside; side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if(neighbour.Element()->MaterialId() != mSimData.mTGeometry.m_HdivWrapMatId) DebugStop();
            //            int hdiv_orient = gel->NormalOrientation(side);
            gel->SetReference(cel);
            TPZCompEl *celwrap = cmesh->CreateCompEl(neighbour.Element());
            //            TPZInterpolationSpace *space = dynamic_cast<TPZInterpolationSpace *>(celwrap);
            //            space->SetSideOrient(neighside,hdiv_orient);
            gel->ResetReference();
            neighbour.Element()->ResetReference();
        }
    }
    
//  {
//    std::ofstream out("cmeshafter.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
//  }
//
    buildmatids.clear();
    
    // WHO CARES WHAT MATID 19 MEANS!!!
    // BEAUTIFUL PROGRAMMING STYLE TO PUT THIS HARDCODED!!
    // inserting the skeleton material ids
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        buildmatids.insert(19);
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(19,2,1);
        cmesh->InsertMaterialObject(nullmat);
    }
    buildmatids.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    buildmatids.insert(bcmatids.begin(),bcmatids.end());
    int borderOrder =mSimData.mTNumerics.m_MortarBorderElementFluxOrder;
    cmesh->SetDefaultOrder(borderOrder);
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    // create all flux elements as discontinuous elements
    cmesh->AutoBuild(buildmatids);
    
  {
    std::ofstream myfile("FluxMeshGood.txt");
    cmesh->Print(myfile);
  }

    
    // vector with computational element pointers for creating fracture elements
    TPZVec<TPZCompEl *> fracsupport(mGeometry->NElements(),0);
    
    
    // set the lagrange level
    {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() < dimension && gel->MaterialId() == mSimData.mTGeometry.m_zeroOrderHdivFluxMatId)
            {
                fracsupport[gel->Index()] = cel;
                cel->Connect(0).SetLagrangeMultiplier(fluxmortarlagrange);
            }
        }
    }
    cmesh->SetDimModel(dimension-1);
    matids.clear();
    bcmatids.clear();
    GetMaterialIds(dimension-1, matids, bcmatids);
    for(auto matid : matids)
    {
        cmesh->InsertMaterialObject(new TPZNullMaterial(matid,dimension-1,1));
    }
    for(auto matid : bcmatids)
    {
        cmesh->InsertMaterialObject(new TPZNullMaterial(matid,dimension-2,1));
    }
    bool complicated = false;
    if(complicated)
    {
        // data structure that stores for each geometric fracture element the two element/side volumetric elements
        // it is facing
        typedef std::pair<TPZGeoElSide,TPZGeoElSide> sidepair;
        // first geometric fracture element
        // second : left/right geoelsides of volumetric elements
        typedef std::pair<TPZGeoEl*,sidepair> gelsideandpair;
        std::list<gelsideandpair> fracElements;
        {
            int64_t nel = mGeometry->NElements();
            for(int64_t el = 0; el<nel; el++)
            {
                TPZGeoEl *gel = mGeometry->Element(el);
                if(!gel || gel->Dimension() != dimension-1) continue;
                int matid = gel->MaterialId();
                // we are looking for fracture elements
                if(matids.find(matid) == matids.end()) continue;
                TPZGeoElSide gelside(gel);
                TPZGeoElSide neighbour = gelside.Neighbour();
                sidepair leftright;
                int nfound = 0;
                // look for volumetric neighbours of the fracture element
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension() == dimension)
                    {
                        if(neighbour.Element()->NormalOrientation(neighbour.Side()) == 1)
                        {
                            leftright.first = neighbour;
                            nfound++;
                        }
                        else
                        {
                            leftright.second = neighbour;
                            nfound++;
                        }
                    }
                    neighbour=neighbour.Neighbour();
                }
                if(nfound != 2) DebugStop();
                gelsideandpair gs;
                gs.first = gel;
                gs.second = leftright;
                fracElements.push_back(gs);
            }
            // now, we have the elements and the H(div) elements connected
            // we will now create the HDiv collapsed elements, link the connects and correct the interface matid
            for(auto it : fracElements)
            {
                TPZGeoEl *fracgel = it.first;
                TPZGeoElSide fracgelside(fracgel);
                TPZGeoElSide leftgelside = it.second.first;
                TPZGeoElSide rightgelside = it.second.second;
                // find the element with material id m_zeroOrderHdivFluxMatId
                TPZGeoElSide zeroflux = fracgelside.HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                TPZCompElSide zerofluxcomp = zeroflux.Reference();
                zerofluxcomp.SetElement(fracsupport[zeroflux.Element()->Index()]);
                if(!zerofluxcomp) DebugStop();
                int zerofluxorder = zerofluxcomp.Element()->Connect(0).Order();
                if(!zeroflux) DebugStop();
                if(!zerofluxcomp) DebugStop();
                TPZGeoElSide intfaceleft;
                TPZGeoElSide intfaceright;
                // based on the order of creation of the geometric wrap elements we reach the interface elements
                {
                    TPZGeoElSide hdivwrap = leftgelside.Neighbour();
                    TPZGeoElSide intface1 = hdivwrap.Neighbour();
                    TPZGeoElSide pressmortar = intface1.Neighbour();
                    intfaceleft = pressmortar.Neighbour();
                }
                {
                    TPZGeoElSide hdivwrap = rightgelside.Neighbour();
                    TPZGeoElSide intface1 = hdivwrap.Neighbour();
                    TPZGeoElSide pressmortar = intface1.Neighbour();
                    intfaceright = pressmortar.Neighbour();
                }
                if(intfaceleft.Element()->MaterialId() != mSimData.mTGeometry.m_negLagrangeMatId) DebugStop();
                if(intfaceright.Element()->MaterialId() != mSimData.mTGeometry.m_posLagrangeMatId) DebugStop();
                // create a second zero flux element and second connect (hybridizing the mesh
                TPZGeoElBC gbc(zeroflux,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                TPZGeoElSide zeroflux2 = gbc.CreatedElement();
                zeroflux.Element()->ResetReference();
                cmesh->SetDefaultOrder(zerofluxorder);
                TPZCompEl *celflux2 = cmesh->CreateCompEl(zeroflux2.Element());
                TPZCompElSide zerofluxcomp2 = zeroflux2.Reference();
                
                zeroflux2.Element()->ResetReference();
                cmesh->SetDefaultOrder(1);
                TPZInterpolationSpace *HDivCollapsed = 0;
                if(zeroflux.Element()->Type() == ETriangle)
                {
                    HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeTriang>(*cmesh,fracgel);
                }
                else if(zeroflux.Element()->Type() == EQuadrilateral)
                {
                    HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,fracgel);
                }
                int nconnects = HDivCollapsed->NConnects();
                // bottom connect index
                int64_t cindex1 = HDivCollapsed->ConnectIndex(nconnects-2);
                // top connect index
                int64_t cindex2 = HDivCollapsed->ConnectIndex(nconnects-1);
                int nsides = fracgel->NSides();
                // changing orientation of top
                HDivCollapsed->SetSideOrient(nsides, -1);
                // linking bottom to first zero flux element - interface left - zeroflux - frac pressure
                HDivCollapsed->SetConnectIndex(nconnects-2, zerofluxcomp.Element()->ConnectIndex(0));
                // linking top to second zero flux element - interface right - zeroflux - frac pressure
                HDivCollapsed->SetConnectIndex(nconnects-1, zerofluxcomp2.Element()->ConnectIndex(0));
                cmesh->ConnectVec()[cindex1].DecrementElConnected();
                cmesh->ConnectVec()[cindex2].DecrementElConnected();
                cmesh->ConnectVec()[cindex1].SetSequenceNumber(-1);
                cmesh->ConnectVec()[cindex2].SetSequenceNumber(-1);
                cmesh->ConnectVec().SetFree(cindex1);
                cmesh->ConnectVec().SetFree(cindex2);
            }
        }
    }
    else // not so complicated
    {
        // create the fracture elements and connect them to the volumetric elements
        int64_t nel = mGeometry->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZGeoEl *gel = mGeometry->Element(el);
            if(!gel || gel->Dimension() != dimension-1) continue;
            int matid = gel->MaterialId();
            // we are looking for fracture elements
            if(matids.find(matid) == matids.end()) continue;
            // now we have a fracture element
            TPZGeoElSide gelside(gel);
            // find the element with material id m_zeroOrderHdivFluxMatId
            TPZGeoElSide zeroflux = gelside.HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
            if(!zeroflux) DebugStop();
            TPZCompElSide zerofluxcomp = zeroflux.Reference();
            zerofluxcomp.SetElement(fracsupport[zeroflux.Element()->Index()]);
            if(!zerofluxcomp) DebugStop();
            int zerofluxorder = zerofluxcomp.Element()->Connect(0).Order();
            
            int domain=0;
            TPZGeoElSide zeroflux2;
            if(mSimData.mTNumerics.m_mhm_mixed_Q){
                if(mSubdomainIndexGel.size() != mGeometry->NElements()) DebugStop();
                domain = mSubdomainIndexGel[zeroflux.Element()->Index()];
                // create a second zero flux element and second connect (hybridizing the mesh
                TPZGeoElBC gbc(zeroflux,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                zeroflux2 = gbc.CreatedElement();
                {
                    auto nel = mGeometry->NElements();
                    mSubdomainIndexGel.Resize(nel,-1);
                    mSubdomainIndexGel[gbc.CreatedElement()->Index()] = domain;
                }
            }
            else{
                // create a second zero flux element and second connect (hybridizing the mesh
                TPZGeoElBC gbc(zeroflux,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                zeroflux2 = gbc.CreatedElement();
            }
                      
             zeroflux.Element()->ResetReference();
             cmesh->SetDefaultOrder(zerofluxorder);
             TPZCompEl *celflux2 = cmesh->CreateCompEl(zeroflux2.Element());
             TPZCompElSide zerofluxcomp2 = zeroflux2.Reference();
             zeroflux2.Element()->ResetReference();
             cmesh->SetDefaultOrder(1);
             TPZInterpolationSpace *HDivCollapsed = 0;
             if(zeroflux.Element()->Type() == ETriangle)
             {
                 HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeTriang>(*cmesh,gel);
             }
             else if(zeroflux.Element()->Type() == EQuadrilateral)
             {
                 HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,gel);
             }
             int nconnects = HDivCollapsed->NConnects();
             // bottom connect index (which will be substituted)
             int64_t cindex1 = HDivCollapsed->ConnectIndex(nconnects);
             // top connect index (which will be substituted)
             int64_t cindex2 = HDivCollapsed->ConnectIndex(nconnects+1);
             int nsides = gel->NSides();
             // changing orientation of top
//             HDivCollapsed->SetSideOrient(nsides, -1);
//             HDivCollapsed->SetSideOrient(nsides-1, -1);
             // linking bottom to first zero flux element - interface left - zeroflux - frac pressure
//             HDivCollapsed->SetConnectIndex(nconnects-2, zerofluxcomp.Element()->ConnectIndex(0));
             HDivCollapsed->SetConnectIndex(nconnects, zerofluxcomp.Element()->ConnectIndex(0));
             // linking top to second zero flux element - interface right - zeroflux - frac pressure
//             HDivCollapsed->SetConnectIndex(nconnects-1, zerofluxcomp2.Element()->ConnectIndex(0));
             HDivCollapsed->SetConnectIndex(nconnects+1, zerofluxcomp2.Element()->ConnectIndex(0));
//             cmesh->ConnectVec()[cindex1].DecrementElConnected();
//             cmesh->ConnectVec()[cindex2].DecrementElConnected();
//             cmesh->ConnectVec()[cindex1].SetSequenceNumber(-1);
//             cmesh->ConnectVec()[cindex2].SetSequenceNumber(-1);
//             cmesh->ConnectVec().SetFree(cindex1);
//             cmesh->ConnectVec().SetFree(cindex2);
        }
    }
    
//    TPZNullMaterial* volume = new TPZNullMaterial(-11,1,1);
//    cmesh->InsertMaterialObject(volume);
//    bcmatids.insert(-11);

    
    cmesh->CleanUpUnconnectedNodes();
    // insert the fracture hdiv elements
    cmesh->SetDefaultOrder(1);
    cmesh->SetDimModel(dimension-1);
    cmesh->SetAllCreateFunctionsHDiv();
    // add the boundary HDiv elements
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild(bcmatids);
    cmesh->SetDimModel(dimension);
    cmesh->ExpandSolution();
#ifdef PZDEBUG
    {
        {
            //            std::ofstream out("FluxMortarMesh.txt");
            //            cmesh->Print(out);
        }
        {
            //            std::ofstream out("FluxMortar.vtk");
            //            TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
        }
        {
            //            std::ofstream out("GMeshAfterFluxMortar.vtk");
            //            TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out);
        }
        
        cmesh->LoadReferences();
        int64_t nelgeo = mGeometry->NElements();
        std::set<int> matidsloc(matids);
        matidsloc.insert(bcmatids.begin(),bcmatids.end());
        for(int64_t el = 0; el<nelgeo; el++)
        {
            TPZGeoEl *gel = mGeometry->Element(el);
            int matid = gel->MaterialId();
            if(matidsloc.find(matid) != matidsloc.end())
            {
                TPZCompEl *cel = gel->Reference();
                if(!cel) DebugStop();
            }
        }
    }
#endif
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *TMRSApproxSpaceGenerator::PressureMortarCmesh(char firstlagrangepressure,char lagrangepressure, char lagrangemortar)
{
    std::set<int> matids, bcmatids;
    int dimension = mGeometry->Dimension();
    int nstate = 1;
    GetMaterialIds(dimension, matids, bcmatids);
    // GetMaterialIds will accumulate in the matids set!
    GetMaterialIds(dimension-1,matids,bcmatids);
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    cmesh->SetName("PressureWithMortar.txt");
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    for (auto it:matids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dimension,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    cmesh->AutoBuild(matids);
    int64_t nconnects = cmesh->NConnects();
    for(int64_t ic = 0; ic<nconnects; ic++)
    {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(lagrangepressure);
    }
    // one connect for each element should be set at lower level
    {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            cel->Connect(0).SetLagrangeMultiplier(firstlagrangepressure);
        }
    }
    {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(mSimData.mTGeometry.m_MortarMatId,dimension,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    // create discontinous elements of dimension-1
    cmesh->SetDimModel(dimension-1);
    std::set<int> mortarids = {mSimData.mTGeometry.m_MortarMatId};
    int borderOrder = mSimData.mTNumerics.m_MortarBorderElementPresOrder;
    cmesh->SetDefaultOrder(borderOrder);
    if(borderOrder == 0)
    {
        cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    cmesh->AutoBuild(mortarids);
    int64_t ncon_new = cmesh->NConnects();
    for (int64_t ic = nconnects; ic<ncon_new; ic++) {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(lagrangemortar);
    }
    cmesh->SetDimModel(dimension);
    
    // insert the fracture hdiv elements
    cmesh->SetDimModel(dimension-1);
    GetMaterialIds(dimension-1, matids, bcmatids);
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild(matids);
    cmesh->SetDimModel(dimension);
    
#ifdef PZDEBUG
    {
        std::ofstream out("PressureMortarCMesh.txt");
        cmesh->Print(out);
    }
    {
        std::ofstream out("PressureMortarCMesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    }
    
#endif
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * TMRSApproxSpaceGenerator::DiscontinuousCmesh(int order, char lagrange){
    
    // -----------> Problem dimension
    if (!mGeometry)
        DebugStop();
    int dimension = mGeometry->Dimension();
    
    // -----------> Creating discontinuous comp mesh
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    cmesh->SetDefaultOrder(order);

    // -----------> Setting space and building mesh for each dimension
    const bool isInsertBCs = false; // Pressure mesh does not require BCs
    for (int idim = 0; idim <= dimension; idim++) {
        std::set<int> matids, bcids;
        // the pressure elements are only created in the original pressure mesh
        if(idim == 1 && order > 0)
        {
            matids.insert(this->mSimData.mTGeometry.m_pressureMatId);
        }
        AddAtomicMaterials(idim,cmesh,matids,bcids,isInsertBCs);
        if (matids.size() == 0) continue;
        
        cmesh->SetDimModel(idim);
        if(order > 0)
        {
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh->SetAllCreateFunctionsDiscontinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
        cmesh->AutoBuild(matids);
        cmesh->InitializeBlock();
    }
    
    if(lagrange > 0){
        int64_t ncon = cmesh->NConnects();
        for(int64_t i=0; i<ncon; i++){
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(lagrange);
        }
    }
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::GroupConnectsBySubdomain(TPZCompMesh *cmesh)
{
    std::map<int,int> domaintoconnect;
    int64_t nel = cmesh->NElements();
	const int cmeshdim = cmesh->Dimension();
    int count = 0;
    for(int64_t el = 0; el<nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel || cel->NConnects() != 1) DebugStop();
        TPZGeoEl *gel = cel->Reference();
		if(gel->Dimension() != cmeshdim) {
			delete cel;
			continue;
		}
        int64_t gelindex = gel->Index();
        int mhm_domain = mSubdomainIndexGel[gelindex];
        int condensed = -1;
        if(domaintoconnect.find(mhm_domain) == domaintoconnect.end())
        {
            condensed = count;
            domaintoconnect[mhm_domain] = count++;
        }
        else
        {
            condensed = domaintoconnect[mhm_domain];
        }
        cel->SetConnectIndex(0, condensed);
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * TMRSApproxSpaceGenerator::TransportCmesh(){
    
    if (!mGeometry) {
        DebugStop();
    }
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZNullMaterial<STATE> * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(0);
    
    //    for(int dim = 0; dim <= dimension; dim++)
    //    {
    //        std::set<int> matids, bcmatids;
    //        GetMaterialIds(dim, matids, bcmatids);
    //
    //        int nstate = 1;
    //        for(auto material_id : matids)
    //        {
    //            volume = new TPZNullMaterial(material_id,dim,nstate);
    //            cmesh->InsertMaterialObject(volume);
    //        }
    //        for(auto material_id : bcmatids)
    //        {
    //            volume = new TPZNullMaterial(material_id,dim-1,nstate);
    //            cmesh->InsertMaterialObject(volume);
    //        }
    //        std::set<int> allmat(matids);
    //        allmat.insert(bcmatids.begin(),bcmatids.end());
    ////        if(allmat.size() == 0) continue;
    //        cmesh->SetDimModel(dim);
    //        cmesh->SetAllCreateFunctionsDiscontinuous();
    //        cmesh->AutoBuild(allmat);
    //    }
    //
    //    TPZTracerFlow * volume = nullptr;
    //    cmesh->SetDefaultOrder(0);
    std::set<int> volIds;
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		int material_id = chunk.second;
		//            volume = new TMRSMultiphaseFlow<TMRSMemory>(material_id,d);
		//            volume->SetDataTransfer(mSimData);
		
		volume = new TPZNullMaterial(material_id,dimension,1);
		volIds.insert(material_id);
		//            volume->SetDataTransfer(mSimData);
		
		cmesh->InsertMaterialObject(volume);
	}
	
    //            volume->SetDataTransfer(mSimData);
    volume = new TPZNullMaterial(11,2,1);
    cmesh->InsertMaterialObject(volume);
    volIds.insert(11);
    
    if (!volume) {
        DebugStop();
    }
    
    std::set<int> boundaryId;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
	for (auto& chunk : mSimData.mTBoundaryConditions.mBCTransportMatIdToTypeValue) {
		int bc_id   = chunk.first;
		std::pair<int,REAL>& typeAndVal = chunk.second;
		int bc_type = typeAndVal.first;
		val2[0]   = typeAndVal.second;
		TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
		boundaryId.insert(bc_id);
		cmesh->InsertMaterialObject(face);
	}
	
    //
    cmesh->InitializeBlock();
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    for (auto gel:mGeometry->ElementVec()) {
        int gelId = gel->MaterialId();
        for (auto Mat_id: volIds) {
            if (gelId == Mat_id) {
                CreateTransportElement(0, cmesh, gel, false);
            }
        }
        for (auto Mat_id: boundaryId) {
            if (gelId == Mat_id) {
                CreateTransportElement(0, cmesh, gel, true);
            }
        }
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int s_order = 0;
    cmesh->SetDefaultOrder(s_order);
    cmesh->ExpandSolution();
    //    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "s_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    std::ofstream file("transportAtomic.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file);
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void  TMRSApproxSpaceGenerator::BuildAuxTransportCmesh(){
    
    // Basically, this mesh creates interface elements between volumes, fractures, intersections.
    // These interfaces will handle the tranport for a given flux from the pressure/flux problem
    
    if (!mGeometry)
        DebugStop();

    int dimension = mGeometry->Dimension();

    mTransportOperator = new TPZCompMesh(mGeometry); // This mesh is only used to set the interface elements
    TPZTracerFlow * volume = nullptr;
    mTransportOperator->SetDefaultOrder(0);
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    
    // ---------------> Adding volume materials
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    std::set<int> volIds;
    std::set<int> boundaryId;
	
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		std::cout << "physical name = " << material_name << std::endl;
		int material_id = chunk.second;
		volIds.insert(material_id);
		volume = new TPZTracerFlow(material_id,dimension); // NOTE: this material is not used for computations
		mTransportOperator->InsertMaterialObject(volume);
	}
    
    // ---------------> Adding volume boundary condition materials
	for (auto& chunk : mSimData.mTBoundaryConditions.mBCTransportMatIdToTypeValue) {
		int bc_id   = chunk.first;
		std::pair<int,REAL>& typeAndVal = chunk.second;
		int bc_type = typeAndVal.first;
		val2[0]   = typeAndVal.second;
		boundaryId.insert(bc_id);
		TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
		mTransportOperator->InsertMaterialObject(face);
	}

    
    // ---------------> Adding fracture materials
    if (isFracSim()){
        int val = 0;
        if(isThereFracIntersection()){
            val = 1;
        }
		for(auto chunk : mSimData.mTFracProperties.m_fracprops){
            auto fracmatid= chunk.first;
            auto fracprop= chunk.second;
            volume = new TPZTracerFlow(fracmatid,dimension-1);
            mTransportOperator->InsertMaterialObject(volume);
            volIds.insert(fracmatid);
            
            for(auto& bcsIds : fracprop.m_fracbc){
                volume = new TPZTracerFlow(bcsIds,dimension-2);
                mTransportOperator->InsertMaterialObject(volume);
                boundaryId.insert(bcsIds);
            }
            int oka=0;
//
//			std::string material_name = chunk.first;
////            sim_data.mTFracProperties.m_fracprops[matid] = fracprop;
//			int material_id = chunk.second;
//			volume = new TPZTracerFlow(material_id,dimension-1);
//			mTransportOperator->InsertMaterialObject(volume);
//			volIds.insert(material_id);
            

            
		}
        
		for(auto chunk : mSimData.mTGeometry.mDomainFracIntersectionNameAndMatId){
			std::string material_name = chunk.first;
			int material_id = chunk.second;
			volume = new TPZTracerFlow(material_id,dimension-2);
			mTransportOperator->InsertMaterialObject(volume);
			volIds.insert(material_id);
		}
             
        if (!volume)
            DebugStop();
        
        // Adding interface materials
        // mInterface_material_idFracSup =It is the interface associated with the upper side of the fracture
        // mInterface_material_idFracInf =It is the interface associated with the lower side of the fracture
        // mInterface_material_idFracFrac =It is the interface between two fracture elements.
        // mIterface_material_idFracBound =It is the interface between the boundaries and fracture elements.
        int fracvol1ID = mSimData.mTGeometry.mInterface_material_idFracInf;
        int fracvol2ID = mSimData.mTGeometry.mInterface_material_idFracSup;
        int fracFracID = mSimData.mTGeometry.mInterface_material_idFracFrac;
        int fracbounId = mSimData.mTGeometry.mInterface_material_idFracBound;
        
//         fracbounId =301;
        // Here, the interface materials are created between fracture/volume, fracture/fracture, and fracture/boundary
        TPZBndCond * face = volume->CreateBC(volume,fracvol1ID,dimension-1,val1,val2);
        mTransportOperator->InsertMaterialObject(face);
        TPZBndCond * face4 = volume->CreateBC(volume,fracvol2ID,dimension-1,val1,val2);
        mTransportOperator->InsertMaterialObject(face4);
        TPZBndCond * face2 = volume->CreateBC(volume,fracFracID,dimension-2,val1,val2);
        mTransportOperator->InsertMaterialObject(face2);
        TPZBndCond * face5 = volume->CreateBC(volume,fracbounId,dimension-2,val1,val2);
        mTransportOperator->InsertMaterialObject(face5);
        boundaryId.insert(fracbounId);
        
        if(mSimData.mTGeometry.m_pressureMatId){
           
            volume = new TPZTracerFlow(mSimData.mTGeometry.m_pressureMatId,dimension-2);
            
//            TPZBndCond * face6 = volume->CreateBC(volume,mSimData.mTGeometry.m_pressureMatId,dimension-2,val1,val2);
            mTransportOperator->InsertMaterialObject(volume);
            volIds.insert(mSimData.mTGeometry.m_pressureMatId);
        }
       
        
    }
  
    // ---------------> Adding Interface materials between volume/volume
    // mInterface_material_id = It is the interface between to volumetric elements
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    {
        TPZTracerFlow * interface = new TPZTracerFlow (transport_matid,dimension-1);
        mTransportOperator->InsertMaterialObject(interface);
    }
    
    // ---------------> Creating Transport Volumetric and Boundary elements
    mTransportOperator->SetAllCreateFunctionsDiscontinuous();
    for (auto gel:mGeometry->ElementVec()) {
        int gelId = gel->MaterialId();
        for (auto Mat_id: volIds) {
            if (gelId == Mat_id) {
                CreateTransportElement(0, mTransportOperator, gel, false);
            }
        }
        for (auto Mat_id: boundaryId) {
            if (gelId == Mat_id) {
                CreateTransportElement(0, mTransportOperator, gel, true);
            }
        }
    }
    
    mTransportOperator->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int s_order = 0;
    mTransportOperator->SetDefaultOrder(s_order);
    mTransportOperator->ExpandSolution();
    mTransportOperator->SetDimModel(dimension);
    {
        mTransportOperator->Reference()->ResetReference();
        mTransportOperator->LoadReferences();
        std::ofstream file2("TransportWithOutInterfaces.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(mTransportOperator, file2);
        CreateInterfaces(mTransportOperator); // Here is where the interface elements are actually created!
        std::ofstream file("TransportWithInterfaces.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(mTransportOperator, file);
        
    }
 
#ifdef PZDEBUG2
    std::ofstream transport("transport_cmesh.txt");
    mTransportOperator->Print(transport);
#endif

}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh * TMRSApproxSpaceGenerator::DiscontinuousCmesh(TPZAlgebraicDataTransfer &Atransfer){
    
    if (!mGeometry) {
        DebugStop();
    }
    int order = 0;
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZL2Projection<STATE> * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(order);
    
    std::map<std::string,int> &DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    int nstate = 1;
    TPZVec<STATE> sol;
    
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		std::cout << "physical name = " << material_name << std::endl;
		int materia_id = chunk.second;
		volume = new TPZL2Projection(materia_id,dimension,nstate, sol);
		cmesh->InsertMaterialObject(volume);
	}
	
    if (!volume) {
        DebugStop();
    }
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // The index of the computational volume elements in the transport mesh identified by material id
    //    std::map<int,TPZVec<int64_t>> fVolumeElements;
    int64_t vol_index = 0;
    TPZVec<int64_t> &elindices = Atransfer.fVolumeElements;
    for (int64_t i=0; i<elindices.size(); i++) {
        int64_t el = elindices[i];
        TPZCompEl *cel = mTransportOperator->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        int matid = gel->MaterialId();
        if(cmesh->FindMaterial(matid) == 0) continue;
        cmesh->CreateCompEl(gel);
        vol_index++;
    }
    cmesh->InitializeBlock();
    
    
    
#ifdef PZDEBUG
    //    std::stringstream file_name;
    //    if (order == 0) {
    //        file_name << "s_cmesh" << ".txt";
    //    }
    //    else{
    //        file_name << "p_cmesh" << ".txt";
    //    }
    //    std::ofstream sout(file_name.str().c_str());
    //    cmesh->Print(sout);
#endif
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildMixedMultiPhysicsCompMesh(int order){
    
    cout << "\n---------------------- Building multiphysics cmesh ----------------------" << endl;
    TPZSimpleTimer timer_mp;
    
    bool cond1 = mSimData.mTNumerics.m_four_approx_spaces_Q;
    bool cond2 = mSimData.mTNumerics.m_mhm_mixed_Q;
    
    // Sanity checks
    if (isFracSim() && !mSimData.mTFracProperties.m_fracprops.size()) {
		DebugStop();
    }
    if (isThereFracIntersection() && mSimData.mTFracIntersectProperties.m_IntersectionId == -10000){
        DebugStop(); // if simulation has frac/frac intersections, this matid should have been set
    }
    
    switch(mSimData.mTNumerics.m_SpaceType)
    {
        case TMRSDataTransfer::TNumerics::E2Space:
            if(cond1) DebugStop();
            if(cond2) DebugStop();
            BuildMixed2SpacesMultiPhysicsCompMesh(order);
            break;
        case TMRSDataTransfer::TNumerics::E4Space:
            //            if(!cond1) DebugStop();
            //            if(!cond2) DebugStop();
            BuildMixed4SpacesMultiPhysicsCompMesh(order);
            break;
        case TMRSDataTransfer::TNumerics::E4SpaceMortar:
            if(!cond1) DebugStop();
            //            if(cond2) DebugStop();
            BuildMixed4SpacesMortarMesh();
            break;
        case TMRSDataTransfer::TNumerics::E4Space1Hybridization:
            if(!cond1) DebugStop();
            //            if(cond2) DebugStop();
            BuildMixed4SpacesHybridized(order);
            break;

        default:
            DebugStop();
    }
    
    std::cout << "\nTotal time build multiphysics cmesh: " << timer_mp.ReturnTimeDouble()/1000 << " seconds" << std::endl;
    cout << "\n---------------------- Finished Building multiphysics cmesh ----------------------" << endl;
    //    std::string name_ref = "mhm_geo";
    //    PrintGeometry(name_ref);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildMixed2SpacesMultiPhysicsCompMesh(int order){
		
	// This function is currently only being used to test 2-D domains living in 3-D.
	// We use it mainly to check if we can reproduce cte pressure for two 2-D domains intersecting each other
	// perpendicularly. This generates an intersection that needs to be hybridized to reproduce cte pressure.
	// Note that even though the domains are 2-D, we do not treat them as fractures. They thus receive a volume material TMRSDarcyFlowWithMem
	
    int dimension = mGeometry->Dimension();
	if(dimension != 2) DebugStop(); // This code is only being tested for 2D. It should work for 3D without any fractures. If you do run with 3D domain, be aware that is not tested for that!
	if(isFracSim()) DebugStop(); // Not implemented to have fractures (i.e. elements of one lower dimension than the mesh's dimension)
	
	
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
    
    mMixedOperator->SetDefaultOrder(order);
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		std::cout << "physical name = " << material_name << std::endl;
		int materia_id = chunk.second;
		volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,dimension);
		//            volume = new TPZMixedDarcyFlow(materia_id, d);
		//            volume->SetPermeability(1.0);
		volume->SetDataTransfer(mSimData);
        volume->SetAxisymmetry(mSimData.mTNumerics.m_is_axisymmetric);
		mMixedOperator->InsertMaterialObject(volume);
	}
	
    if (!volume) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
	for(auto &chunk : mSimData.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
		int bc_id   = chunk.first;
		std::pair<int,REAL>& typeAndVal = chunk.second;
		int bc_type = typeAndVal.first;
		val2[0]  = typeAndVal.second;
		TPZBndCondT<REAL> * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
		if (HasForcingFunctionBC()) {
			face->SetForcingFunctionBC(mForcingFunctionBC,1);
		}
		mMixedOperator->InsertMaterialObject(face);
	}
    
    TPZManVector<TPZCompMesh *, 3> mesh_vec(3);
    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order,1);
    mesh_vec[2] = TransportCmesh();
    TPZManVector<int,5> active_approx_spaces(3);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 0;
    mMixedOperator->SetDimModel(dimension);
    
    if (isThereFracIntersection()) {
        mHybridizer = new TPZHybridizeHDiv(mesh_vec);
        HybridizeIntersections(mesh_vec);
        mHybridizer->InsertPeriferalMaterialObjects(mMixedOperator);
    }
    //    mMixedOperator->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec);
    
    // Creates interface elements in case there is hybridization for fracture intersection
    if (mHybridizer){
        CreateIntersectionInterfaceElements(mesh_vec);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "mixed_cmesh_two_spaces" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    mMixedOperator->Print(sout);
#endif
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildMixed4SpacesMortarMesh(){
    int dimension = mGeometry->Dimension();
    std::cout << __PRETTY_FUNCTION__ << " on input nel geom " << mGeometry->NElements() << std::endl;
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    std::set<int> matsWithMem, matsWithOutMem;
    GetMaterialIds(dimension, matsWithMem, matsWithOutMem);
    
    InsertGeoWrappersForMortar();
    TPZManVector<TPZCompMesh *> meshvec(4,0);
    // hdiv mesh
    char fluxmortar = 5;
    meshvec[0] = HDivMortarFluxCmesh(fluxmortar);
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> fatherAndSons;
        TPZReservoirTools::TakeFatherSonsCorrespondence(meshvec[0], this->mSubdomainIndexGel, fatherAndSons);
        TPZReservoirTools::PutFluxElementsinSubdomain(meshvec[0], this->mSubdomainIndexGel, fatherAndSons);
        TPZReservoirTools::AddDependency(fatherAndSons);
    }
    // pressure mesh
    char firstpressurelagrange = 1;
    char pressurelagrange = 3;
    char pressuremortar = 4;
    meshvec[1] = PressureMortarCmesh(firstpressurelagrange,pressurelagrange,pressuremortar);
    
    if (isThereFracIntersection()) {
        mHybridizer = new TPZHybridizeHDiv(meshvec);
		set<int> fracmatIds;
		if(mSimData.mTFracProperties.m_fracprops.size()){
			for (auto el : mSimData.mTFracProperties.m_fracprops) fracmatIds.insert(el.first);
		}
		else{
            std::map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
			for (it = mSimData.mTFracProperties.m_fracprops.begin(); it != mSimData.mTFracProperties.m_fracprops.end(); it++){
				int matfracid = it->first;
				fracmatIds.insert(matfracid);
			}
		}
		mHybridizer->IdsToHybridize() = fracmatIds;
        const int intersectionPressureLossId = mSimData.mTFracIntersectProperties.m_IntersectionPressureLossId;
        if (intersectionPressureLossId > -10000) {
            mHybridizer->fHDivWrapMatid = intersectionPressureLossId;
        }
        HybridizeIntersections(meshvec);
    }
    
    
    std::ofstream file2("PressureCmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(meshvec[1], file2);
    
    // distributed flux mesh
    int porder = 0;
    char distfluxlagrange = 2;
    meshvec[2] = DiscontinuousCmesh(porder,distfluxlagrange);
    // constant pressure mesh
    char avpressurelagrange = 6;
    meshvec[3] = DiscontinuousCmesh(porder,avpressurelagrange);
    // transport mesh
    // I believe we don't need a transport mesh anymore
    // meshvec[4] = TransportCmesh();
    
    // create the multiphysics mesh
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
    // TPZMixedDarcyFlow *volume = nullptr;
    mMixedOperator->SetDefaultOrder(1);
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    std::cout << "Creating material objects\n";
    
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		std::cout << "physical name = " << material_name <<
		" material id " << chunk.second << " dimension " << dimension << std::endl;
		int materia_id = chunk.second;
		volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,dimension);
		TMRSMemory defaultmem;
		// neste ponto podemos inserir as propriedades de permeabilidade absoluta
		volume->SetDefaultMem(defaultmem);
		//            volume = new TPZMixedDarcyFlow(materia_id, d);
		//             volume->SetPermeability(1.0);
		volume->SetDataTransfer(mSimData);
        volume->SetAxisymmetry(mSimData.mTNumerics.m_is_axisymmetric);
		mMixedOperator->InsertMaterialObject(volume);
	}
	
    TMRSDarcyFractureFlowWithMem<TMRSMemory> * fracmat = nullptr;
    for(auto chunk : mSimData.mTGeometry.mDomainFracNameAndMatId){
        std::string material_name = chunk.first;
        std::cout << "physical name = " << material_name <<
        " material id " << chunk.second << " dimension " << dimension-1 << std::endl;
        int materia_id = chunk.second;
        matsWithMem.insert(materia_id);
        fracmat = new TMRSDarcyFractureFlowWithMem<TMRSMemory>(materia_id,dimension-1);
        TMRSMemory defaultmem;
        // neste ponto podemos inserir as propriedades de permeabilidade absoluta
        //        defaultmem.m_kappa = 0.00000001;
        //        defaultmem.m_kappa_inv = 1.0/0.00000001;
        //        defaultmem.SetPermeability(0.00000001);
        //        defaultmem.m_p=100000000.0;
        //        fracmat->SetDefaultMem(defaultmem);
        //            volume = new TPZMixedDarcyFlow(materia_id, d);
        //             volume->SetPermeability(1.0);
        fracmat->SetDataTransfer(mSimData);
        mMixedOperator->InsertMaterialObject(fracmat);
    }
    if (!volume) {
        DebugStop();
    }
    
    {
        TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
		for(auto &chunk : mSimData.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
			int bc_id   = chunk.first;
			std::pair<int,REAL>& typeAndVal = chunk.second;
			int bc_type = typeAndVal.first;
			val2[0]  = typeAndVal.second;
			TPZBndCondT<REAL> * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
			if (HasForcingFunctionBC()) {
				face->SetForcingFunctionBC(mForcingFunctionBC,1);
			}
			mMixedOperator->InsertMaterialObject(face);
		}
    }
    
    {
        TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
		
		if(!fracmat) DebugStop();
		for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
			int bc_id   = chunk.first;
			std::pair<int,REAL>& typeAndVal = chunk.second;
			int bc_type = typeAndVal.first;
			matsWithOutMem.insert(bc_id);
			val2[0] = typeAndVal.second;
			if(bc_type == 2){
				val2[0] = 0.0;
				val1(0,0) = typeAndVal.second;
			}
			TPZBndCondT<REAL>* face = fracmat->CreateBC(volume,bc_id,bc_type,val1,val2);
			if (HasForcingFunctionBC()) {
				face->SetForcingFunctionBC(mForcingFunctionBC,1);
			}
			mMixedOperator->InsertMaterialObject(face);
		}
    }
    
    if(mHybridizer) {
        mHybridizer->InsertPeriferalMaterialObjects(mMixedOperator);
    }
    
    {
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<STATE> *nullmat = new TPZNullMaterialCS(mSimData.mTGeometry.m_HdivWrapMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
    {
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<STATE> *nullmat = new TPZNullMaterialCS(mSimData.mTGeometry.m_MortarMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
    {
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<STATE> *nullmat = new TPZNullMaterialCS(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
    //    mGeometry->ResetReference();
    //    mMixedOperator->LoadReferences();
    
    //    TMRSDarcyFlowWithMem<TMRSMemory>*volume2 = new TMRSDarcyFlowWithMem<TMRSMemory>(19,2);
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        TPZNullMaterialCS<STATE> *volume2 = new TPZNullMaterialCS(mSimData.mTGeometry.m_skeletonMatId,2,1);
        mMixedOperator->InsertMaterialObject(volume2);
    }
    
    
    matsWithOutMem.insert(mSimData.mTGeometry.m_HdivWrapMatId);
    matsWithOutMem.insert(mSimData.mTGeometry.m_MortarMatId);
    matsWithOutMem.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    matsWithOutMem.insert(mSimData.mTGeometry.m_skeletonMatId);
    
    TPZManVector<int> active_approx_spaces(4,1);
    //    active_approx_spaces[4] = 0;
    mMixedOperator->SetDimModel(3);
    gSinglePointMemory = true;
    
    // NS to Jose: Should this method with these arguments be commited in PZ???

//    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,meshvec,matsWithMem, matsWithOutMem);
//    DebugStop(); // look up
    // NS to Jose: Using this for now. Erase later
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,meshvec);
    
    // Creates interface elements in case there is hybridization for fracture intersection
    if (mHybridizer){
        CreateIntersectionInterfaceElements(meshvec);
    }


    //Insert fractures properties
    InitializeFracProperties(mMixedOperator);
    // insert the interface elements
    InsertInterfaceElements();
#ifdef PZDEBUG
    {
        std::stringstream file_name;
        file_name  << "mixed_cmesh_four_space_mortar" << ".txt";
        std::ofstream sout(file_name.str().c_str());
//        mMixedOperator->Print(sout);
    }
#endif
    mMixedOperator->ComputeNodElCon();
    
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        if(mSimData.mTNumerics.m_UseSubstructures_Q)
        {
            std::cout<<"Num Eq Mixed: "<<mMixedOperator->NEquations()<<std::endl;
            HideTheElements(mMixedOperator); // creating subcompmeshes
            int nels = mMixedOperator->NElements();
            for(int iel =0; iel<nels; iel++){
                TPZCompEl *cel = mMixedOperator->Element(iel);
                if(!cel){continue;}
                TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
                if(subcmesh){
                    
                    
                    std::set<int64_t> seed, groups;
                    subcmesh->ComputeNodElCon();
                    
                    
                    
                    TPZReservoirTools::TakeSeedElements(subcmesh, seed);
                    TPZReservoirTools::GroupNeighbourElements(subcmesh,seed,groups );
                    subcmesh->ComputeNodElCon();
                    std::set<int> volmatId;
                    //                    std::ofstream subcm("PrintSubm.txt");
                    //                    subcmesh->Print(subcm);
                    int nel = subcmesh->NElements();
                    TPZReservoirTools::CondenseElements(subcmesh, pressuremortar, false,volmatId);
                    std::set<int64_t> groups2;
                    // this will act only on volumetric elements
                    TPZReservoirTools::GroupNeighbourElements(subcmesh, groups, groups2);
                    subcmesh->ComputeNodElCon();
                    // this shouldn't affect the fracture elements as they won't have condensable connects
                    TPZReservoirTools::CondenseElements(subcmesh, fluxmortar, false);
                    subcmesh->ComputeNodElCon();
                    //                    int numThreads =0;
                    //                    int preconditioned =0;
                    //                    TPZAutoPointer<TPZGuiInterface> guiInterface;
                    //                    subcmesh->SetAnalysisSkyline(numThreads, preconditioned, guiInterface);
//                    subcmesh->SetAnalysisSkyline(0, 0, guiInterface);
                    //                    std::ofstream subcm2("PrintSubm2.txt");
                    //                    subcmesh->Print(subcm2);
                    
                }
            }
            mMixedOperator->ExpandSolution();
            std::cout<<"Num Eq Mixed MHM: "<<mMixedOperator->NEquations()<<std::endl;
            
        }
        else{
            std::cout<<"Num Eq Mixed: "<<mMixedOperator->NEquations()<<std::endl;

            // group and condense the H(div) space (only dimension of the mesh)
            std::set<int64_t> seed, groups;
            int64_t nel = mMixedOperator->NElements();
            int dim = mMixedOperator->Dimension();
            for (int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = mMixedOperator->Element(el);
                if(!cel){
                    continue;
                }
                TPZGeoEl *gel = cel->Reference();
                if(gel->Dimension() == dim) seed.insert(el);
            }
            // this will only group volumetric elements
            TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, seed, groups);
            mMixedOperator->ComputeNodElCon();
            
            std::set<int> volmatId;
            volmatId.insert(10);
            TPZReservoirTools::CondenseElements(mMixedOperator, pressuremortar, false,volmatId);
            
            
            std::set<int64_t> groups2;
            
            // this will act only on volumetric elements
            TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, groups, groups2);
            mMixedOperator->ComputeNodElCon();
            // this shouldn't affect the fracture elements as they won't have condensable connects
            TPZReservoirTools::CondenseElements(mMixedOperator, fluxmortar, false);
            mMixedOperator->ComputeNodElCon();
            std::cout<<"Num Eq Mixed: "<<mMixedOperator->NEquations()<<std::endl;

        }
        
    }
    else
    {
        
        // group and condense the H(div) space (only dimension of the mesh)
        std::cout<<"Num Eq Mixed: "<<mMixedOperator->NEquations()<<std::endl;
        std::set<int64_t> seed, groups;
        int64_t nel = mMixedOperator->NElements();
        int dim = mMixedOperator->Dimension();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = mMixedOperator->Element(el);
            if(!cel){
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() == dim) seed.insert(el);
        }
        // this will only group volumetric elements
        TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, seed, groups);
        mMixedOperator->ComputeNodElCon();
        
        std::set<int> volmatId;
        volmatId.insert(10);
        TPZReservoirTools::CondenseElements(mMixedOperator, pressuremortar, false,volmatId);
        
        
        std::set<int64_t> groups2;
        
        // this will act only on volumetric elements
        
        
        TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, groups, groups2);
        mMixedOperator->ComputeNodElCon();
        // this shouldn't affect the fracture elements as they won't have condensable connects
        TPZReservoirTools::CondenseElements(mMixedOperator, fluxmortar, false);
        std::cout<<"Num Eq Mixed: "<<mMixedOperator->NEquations()<<std::endl;
    }
    
    
    
#ifdef PZDEBUG
    {
        std::stringstream file_name;
        file_name  << "mixed_cmesh_four_space_mortar_two_condense" << ".txt";
        std::ofstream sout(file_name.str().c_str());
        //        mMixedOperator->Print(sout);
    }
#endif
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::InsertInterfaceElements()
{
    int dim = mGeometry->Dimension();
    std::set<int> matids, bcmatids, fracmatids, fracbcmatids;
    GetMaterialIds(dim, matids, bcmatids);
    GetMaterialIds(dim-1, fracmatids, fracbcmatids);
    //    bcmatids.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    TPZLagrangeMultiplierCS<STATE> *mat1 = new TPZLagrangeMultiplierCS<STATE>(mSimData.mTGeometry.m_posLagrangeMatId,dim,1);
    mMixedOperator->InsertMaterialObject(mat1);
    TPZLagrangeMultiplierCS<STATE> *mat2 = new TPZLagrangeMultiplierCS<STATE>(mSimData.mTGeometry.m_negLagrangeMatId,dim,1);
    mat2->SetMultiplier(-1.);
    mMixedOperator->InsertMaterialObject(mat2);
    mGeometry->ResetReference();
    mMixedOperator->LoadReferences();
    int64_t nel = mMixedOperator->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = mMixedOperator->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        // the interface will be generated between HDiv wrapper and the mortar pressure
        if(matid != mSimData.mTGeometry.m_HdivWrapMatId && matid != mSimData.mTGeometry.m_MortarMatId)
            continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide leftgel(gelside), rightgel;
        TPZGeoElSide IntfaceSide = leftgel.Neighbour();
        if(matid == mSimData.mTGeometry.m_HdivWrapMatId)
        {
            rightgel = IntfaceSide.Neighbour();
        }
        else if(matid == mSimData.mTGeometry.m_MortarMatId)
        {
            // from a mortar element we need to add an interface to either a boundary flux
            // or a zero order flux
            // if the mortar element has a fracture neighbour, then there are two neighbouring
            // zero order flux elements. Depending on the interface matid we need to connect to either
            //  zero order flux element
            int domain =0;
            if(mSimData.mTNumerics.m_mhm_mixed_Q){
                domain = mSubdomainIndexGel[gel->Index()];
                if(domain == -1) DebugStop();
            }
            
            TPZGeoElSide BCGelside = gelside.HasNeighbour(bcmatids);
            TPZGeoElSide FracGelside = gelside.HasNeighbour(fracmatids);
            TPZGeoElSide Zerofluxside = gelside.HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
            
            int intfacematid = IntfaceSide.Element()->MaterialId();
            if(BCGelside)
            {
                rightgel = BCGelside;
            }
            else if(Zerofluxside)
            {
                if(mSimData.mTNumerics.m_mhm_mixed_Q){
                    int64_t zerodomain = mSubdomainIndexGel[Zerofluxside.Element()->Index()];
                    if(zerodomain != domain)
                    {
                        Zerofluxside = Zerofluxside.Neighbour().HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                        if(!Zerofluxside) DebugStop();
                        zerodomain = mSubdomainIndexGel[Zerofluxside.Element()->Index()];
                        if(zerodomain != domain) DebugStop();
                    }
                    // if the lagrange multiplier is positive, the this is the top lagrange multiplier
                    // the top flux element is the second zero flux element
                    if(FracGelside && intfacematid == mSimData.mTGeometry.m_posLagrangeMatId)
                    {
                        Zerofluxside = Zerofluxside.Neighbour();
                        if(Zerofluxside.Element()->MaterialId() != mSimData.mTGeometry.m_zeroOrderHdivFluxMatId)
                        {
                            DebugStop();
                        }
                    }
                    rightgel = Zerofluxside;
                }
                else{
                    
                    if(FracGelside && intfacematid == mSimData.mTGeometry.m_posLagrangeMatId)
                    {
                        Zerofluxside = Zerofluxside.Neighbour();
                        if(Zerofluxside.Element()->MaterialId() != mSimData.mTGeometry.m_zeroOrderHdivFluxMatId)
                        {
                            DebugStop();
                        }
                    }
                    rightgel = Zerofluxside;
                }
                
            }
            else
            {
                DebugStop();
            }
        }
        TPZCompElSide leftcel = leftgel.Reference();
        TPZCompElSide rightcel = rightgel.Reference();
        TPZMultiphysicsInterfaceElement *intface =
        new TPZMultiphysicsInterfaceElement(*mMixedOperator,IntfaceSide.Element(),leftcel,rightcel);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::InsertGeoWrappersForMortar()
{
    int dimension = mGeometry->Dimension();
    std::set<int> matids, bcmatids;
    GetMaterialIds(dimension, matids, bcmatids);
    bcmatids.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    // around each volumetric element create an
    // - hdiv wrap element
    // - interface geo element
    // - pressure element
    // - second interface element
    // - a flux hdiv boundary element (conditionally)
    int64_t nel = mGeometry->NElements();
    int nElSubDomain = mSubdomainIndexGel.size();
//    if(mSimData.mTNumerics.m_mhm_mixed_Q){
//        if(nel != nElSubDomain){
//            DebugStop();
//        }
//    }
    
    // For each side of a 3D element, create geowrappers for it
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = mGeometry->Element(el);
        if(!gel || gel->Dimension()!= dimension || gel->HasSubElement()) continue;
        int firstside = 0;
        for(int dim = 0; dim < dimension-1; dim++)
        {
            firstside+=gel->NSides(dim);
        }
        for (int side=firstside; side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            GeoWrappersForMortarGelSide(gelside,bcmatids);
        }
    }
    
    nel = mGeometry->NElements();
    nElSubDomain = mSubdomainIndexGel.size();
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        if(nel != nElSubDomain){
            DebugStop();
        }
        std::ofstream out("gmesh_withwrap.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out, mSubdomainIndexGel);
    }
    else{
        std::ofstream out("gmesh_withwrap.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::GeoWrappersForMortarGelSide(TPZGeoElSide &gelside, std::set<int> bcmatids){
    
    TPZGeoEl *gel = gelside.Element();
    int side = gelside.Side();
    int gelindex = gel->Index();
    int hdiv_orient = gel->NormalOrientation(side);
    int subDomainIndex =-1;
    int subDomainIndexNeig =-1;
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        subDomainIndex = mSubdomainIndexGel[gelindex];
        subDomainIndexNeig = FindNeighSubDomain(gelside);
        if(subDomainIndex<0){
            DebugStop();
        }
    }
    int first_lagrange = mSimData.mTGeometry.m_posLagrangeMatId;
    int second_lagrange = mSimData.mTGeometry.m_negLagrangeMatId;
    bool cond1 = hdiv_orient < 0 ;
    bool cond2 = false;
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        cond2 =subDomainIndexNeig!=-1 && (subDomainIndexNeig > subDomainIndex);
    }
    if(cond1 )
    {
        first_lagrange = mSimData.mTGeometry.m_negLagrangeMatId;
        second_lagrange = mSimData.mTGeometry.m_posLagrangeMatId;
    }
    if(cond2){
        first_lagrange = mSimData.mTGeometry.m_negLagrangeMatId;
        second_lagrange = mSimData.mTGeometry.m_posLagrangeMatId;
    }
    
    
    TPZGeoElBC gbc1(gelside,mSimData.mTGeometry.m_HdivWrapMatId);
    TPZGeoElSide gelwrapside(gbc1.CreatedElement());
    int nBCCreated = 4;
    int index1 =gbc1.CreatedElement()->Index();
    
    TPZGeoElBC gbc2(gelwrapside,first_lagrange);
    TPZGeoElSide gelintface1(gbc2.CreatedElement());
    int index2 =gbc2.CreatedElement()->Index();
    
    
    TPZGeoElBC gbc3(gelintface1,mSimData.mTGeometry.m_MortarMatId);
    TPZGeoElSide gelmortar(gbc3.CreatedElement());
    int index3 =gbc3.CreatedElement()->Index();
    
    TPZGeoElBC gbc4(gelmortar,second_lagrange);
    TPZGeoElSide gelinterface2(gbc4.CreatedElement());
    int index4 =gbc4.CreatedElement()->Index();
    
    int size =mSubdomainIndexGel.size();
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        mSubdomainIndexGel.Resize(size + nBCCreated,-1);
        if(index1 != size || index2!=(size+1) || index3!=(size+2) || index4!=(size+3)){
            DebugStop();
        }
        mSubdomainIndexGel[index1]=subDomainIndex;
        mSubdomainIndexGel[index2]=subDomainIndex;
        mSubdomainIndexGel[index3]=subDomainIndex;
        mSubdomainIndexGel[index4]=subDomainIndex;
    }
    
    
    // this method has to be adjusted if we create MHM meshes
    if(!gelinterface2.HasNeighbour(bcmatids)){
        
        TPZGeoElBC gbc5(gelinterface2,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
        if(mSimData.mTNumerics.m_mhm_mixed_Q){
            mSubdomainIndexGel.Resize(size + nBCCreated +1,-1);
            mSubdomainIndexGel[size + nBCCreated]=subDomainIndex;
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

int TMRSApproxSpaceGenerator::FindNeighSubDomain(TPZGeoElSide &gelside){
    TPZGeoEl *gel = gelside.Element();
    int indexGel = gel->Index();
    int subdomainIndex =-1;
    TPZGeoElSide neigh = gelside.Neighbour();
    while(neigh!= gelside){
        int neighIndex=neigh.Element()->Index();
        if(neigh.Element()->Dimension() != gel->Dimension()){
            neigh = neigh.Neighbour();
            continue;
        }
        subdomainIndex = mSubdomainIndexGel[neighIndex];
        break;
    }
    
    return subdomainIndex;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::GetMaterialIds(int dim, std::set<int> &matids, std::set<int> &bcmatids)
{
#ifdef PZDEBUG
//    std::cout << "\n===> GetMaterialIds - Identifying material objects for dimension " << dim << std::endl;
#endif
    if(dim == mGeometry->Dimension())
    {
        // @TODO why not use the reservoir property datastructure?
		// NS: Do you think we should we save the data in two different data structures? Or just erase this one and keep everything in the reservoir property datastructure
        std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
        for (auto chunk : DomainDimNameAndPhysicalTag) {
#ifdef PZDEBUG
            std::string material_name = chunk.first;
//            std::cout << "physical name = " << material_name << " matid " << chunk.second<< std::endl;
#endif
            matids.insert(chunk.second);
        }
        
        for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
            int bc_id   = chunk.first;
#ifdef PZDEBUG
//            std::cout << "boundary condition matid " << bc_id << std::endl;
#endif
            bcmatids.insert(bc_id);
        }
    }
    else if (dim == mGeometry->Dimension()-1) {
        for (auto chunk : mSimData.mTGeometry.mDomainFracNameAndMatId) {
#ifdef PZDEBUG
            std::string material_name = chunk.first;
//            std::cout << "physical name = " << material_name << " matid " << chunk.second<< std::endl;
#endif
            matids.insert(chunk.second);
        }
		
		for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
			int bc_id   = chunk.first;
#ifdef PZDEBUG
			//                std::cout << "boundary condition matid " << bc_id << std::endl;
#endif
			bcmatids.insert(bc_id);
		}
        
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::InsertGeoWrappers()
{
    int dimension = mGeometry->Dimension();
    std::set<int> matids, bcmatids;
    GetMaterialIds(dimension, matids, bcmatids);

    const int64_t nel = mGeometry->NElements();
    for(int64_t el = 0; el<nel; el++)    {
        TPZGeoEl *gel = mGeometry->Element(el);
        if(!gel || gel->Dimension()!= dimension || gel->HasSubElement()) continue;
        int firstside = 0;
        for(int dim = 0; dim < dimension-1; dim++)
        {
            firstside+=gel->NSides(dim);
        }
        for (int side=firstside; side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neigh = gelside.HasNeighbour(bcmatids);
            if (neigh && neigh != gelside.Neighbour()) {
                DebugStop(); // We require that the bc is the first neighbor of a volumetric element
            }
            if(!neigh){
                const int sideorient = gel->NormalOrientation(side);
                if (!gelside.HasNeighbour(mSimData.mTGeometry.m_posLagrangeMatId)) {
                    TPZGeoElBC gbcpres(gelside,mSimData.mTGeometry.m_pressureMatId);
                }
                if (sideorient == 1) {
                    TPZGeoElBC gbcinterface(gelside,mSimData.mTGeometry.m_posLagrangeMatId);
                }
                else if(sideorient == -1){
                    TPZGeoElBC gbcinterface(gelside,mSimData.mTGeometry.m_negLagrangeMatId);
                }
                else
                    DebugStop();
                
                TPZGeoElBC gbcwrap(gelside,mSimData.mTGeometry.m_HdivWrapMatId);
            }
        }
    }
    
    std::ofstream out("gmesh_withwrap.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out);
    ofstream outtxt("gmesh.txt");
    mGeometry->Print(outtxt);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildMixed4SpacesHybridized(int order) {
    
    // Jun 2022: The idea for this method was conceived but never implemented. At this date,
	// we are using no hibridization between the 3d elements
    DebugStop();
    
    // ========================================================
    // 1) Create GeoEls used for hybridization
    // Two Geoels are created at the faces of volumetric elements (except where there is bcs):
    // Wraps and interface (positive or negative). It is REQUIRED that a volume face has that
    // neighborhood order: wrap and interface
    // In between the interface, a pressure lagrange element is created
    // TODO: NATHAN CHANGE THE NAME to posInterfaceLagrangeMatID
    InsertGeoWrappers();
    
    // ========================================================
    // 2) Create Flux HDiv computational mesh
//    HdivFluxMeshHybridized();
    
    // ========================================================
    // 3) Create Pressure L2 computational mesh
    
    // ========================================================
    // 4) Create Distributed Flux computational mesh

    // ========================================================
    // 5) Create Constant Pressure computational mesh

    // ========================================================
    // 6) Create Multiphysics computational mesh
//    const int dimension = mGeometry->Dimension();
//    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::AddMultiphysicsMaterialsToCompMesh(const int order, std::set<int> &MatsWithmem, std::set<int> &MatsWitOuthmem) {
    
    const int dimension = mGeometry->Dimension();
    mMixedOperator->SetDefaultOrder(order);
    
    // ---------------> Adding volume materials
    cout << "\n---------------------- Inserting materials in multiphysics cmesh ----------------------" << endl;
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		int material_id = chunk.second;
		volume = new TMRSDarcyFlowWithMem<TMRSMemory>(material_id,dimension);
		volume->SetDataTransfer(mSimData);
        volume->SetAxisymmetry(mSimData.mTNumerics.m_is_axisymmetric);
		mMixedOperator->InsertMaterialObject(volume);
		std::cout << "Added volume material w/ physical name = " << material_name << " and id = " << material_id << std::endl;
		MatsWithmem.insert(material_id);
	}
	
    
    if(!volume) DebugStop();
        
    // ---------------> Adding volume boundary condition materials
    auto& FlowFunctionBCmap = mSimData.mTBoundaryConditions.mBCFlowMatIdToFunctionId;
	for(auto &chunk : mSimData.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
		TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
		int bc_id   = chunk.first;
		std::pair<int,REAL>& typeAndVal = chunk.second;
		int bc_type = typeAndVal.first;
		val2[0]  = typeAndVal.second;
		TPZBndCondT<REAL> * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        int functionBCId = FlowFunctionBCmap[bc_id].first;
        if (functionBCId != 0) {
            auto functionBC = FlowFunctionBCmap[bc_id].second;
            face->SetForcingFunctionBC(functionBC,1);
        }
		mMixedOperator->InsertMaterialObject(face);
		MatsWitOuthmem.insert(bc_id);
		std::cout << "Added volume BC material w/ id = " << bc_id << " and type = " << bc_type << std::endl;
	}
    
    // ---------------> Adding fracture materials
    if (isFracSim()) {
        TMRSDarcyFractureFlowWithMem<TMRSMemory> * fracmat = nullptr;
        for(auto chunk : mSimData.mTGeometry.mDomainFracNameAndMatId){
            std::string material_name = chunk.first;
            int material_id = chunk.second;
            fracmat = new TMRSDarcyFractureFlowWithMem<TMRSMemory>(material_id,dimension-1);
            fracmat->SetDataTransfer(mSimData);
            mMixedOperator->InsertMaterialObject(fracmat);
            MatsWithmem.insert(material_id);
            std::cout << "Added frac material w/ physical name = " << material_name << " and id = " << material_id << std::endl;
        }
        if (!fracmat) {
            DebugStop();
        }
        
        // ---------------> Adding fracture boundary condition materials
		if(!fracmat) DebugStop();
		for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
			TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
			int bc_id   = chunk.first;
			std::pair<int,REAL>& typeAndVal = chunk.second;
			int bc_type = typeAndVal.first;
            val2[0] = typeAndVal.second;
			if(bc_type == 2){
				val2[0] = 0.0;
				val1(0,0) = typeAndVal.second;
			}
			TPZBndCondT<REAL>* face = fracmat->CreateBC(fracmat,bc_id,bc_type,val1,val2);
			if (HasForcingFunctionBC()) {
				face->SetForcingFunctionBC(mForcingFunctionBC,1);
			}
			mMixedOperator->InsertMaterialObject(face);
//			std::cout << "Added frac BC material w/ id = " << bc_id << " and type = " << bc_type << std::endl;
			MatsWitOuthmem.insert(bc_id);

		}
        // add the fracture glue material object
        TMRSDarcyFractureGlueFlowWithMem *mat = new TMRSDarcyFractureGlueFlowWithMem(mSimData.mTFracIntersectProperties.m_FractureGlueId,
                                                mSimData.mTFracIntersectProperties.m_FractureGluePerm);
        mMixedOperator->InsertMaterialObject(mat);
        MatsWithmem.insert(mSimData.mTFracIntersectProperties.m_FractureGlueId);
    }
    {
        // material of pressure lagrange multiplier
        int pressureMatId = mSimData.mTGeometry.m_pressureMatId;
        MatsWitOuthmem.insert(pressureMatId);
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<> *nullmat = new TPZNullMaterialCS<>(pressureMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
    if(mHybridizer){
        mHybridizer->InsertPeriferalMaterialObjects(mMixedOperator);
        MatsWitOuthmem.insert(mHybridizer->fLagrangeInterface);
        MatsWitOuthmem.insert(mHybridizer->fHDivWrapMatid);
        MatsWitOuthmem.insert(mHybridizer->fInterfaceMatid.first);
        MatsWitOuthmem.insert(mHybridizer->fInterfaceMatid.second);
    }
	
	if(mSimData.mTNumerics.m_mhm_mixed_Q){
//		TPZNullMaterialCS<STATE> *volume2 = new TPZNullMaterialCS<>(mSimData.mTGeometry.m_skeletonMatId,2,1);
//		mMixedOperator->InsertMaterialObject(volume2);
//		MatsWitOuthmem.insert(mSimData.mTGeometry.m_skeletonMatId);
        TPZNullMaterialCS<STATE> *volume3 = new TPZNullMaterialCS<>(mSimData.mTGeometry.m_skeletonMatId-1,2,1);
        mMixedOperator->InsertMaterialObject(volume3);
        MatsWitOuthmem.insert(mSimData.mTGeometry.m_skeletonMatId-1);
	}
    
    mMixedOperator->SetDimModel(dimension);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void TMRSApproxSpaceGenerator::GetTransportMaterials(std::set<int> &MatsWithmem, std::set<int> &MatsWitOuthmem){
    
    const int dimension = mGeometry->Dimension();
    // ---------------> Adding volume materials
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;

	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		int material_id = chunk.second;
		MatsWithmem.insert(material_id);
	}
	
    
    if(MatsWithmem.size()==0) DebugStop();
        
    // ---------------> Adding volume boundary condition materials
    for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
        int bc_id   = chunk.first;
        MatsWitOuthmem.insert(bc_id);
    }
    
    // ---------------> Adding fracture materials
    if (isFracSim()) {
        for(auto chunk : mSimData.mTGeometry.mDomainFracNameAndMatId){
            std::string material_name = chunk.first;
            int material_id = chunk.second;
            MatsWithmem.insert(material_id);
        }
//        if (!fracmat) {
//            DebugStop();
//        }
        
        // ---------------> Adding fracture boundary condition material
        for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
            int bc_id   = chunk.first;
            MatsWitOuthmem.insert(bc_id);
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetLagrangeMultiplier4Spaces(TPZVec<TPZCompMesh *>& mesh_vec) {
    int nmesh = mesh_vec.size();
    if(nmesh > 6) nmesh = 6;
    for (int imesh = 0; imesh < nmesh; imesh++)
    {
        char lagrange = imesh;
        if(imesh == 5) lagrange = 6;
        // First is pressure mesh
        int64_t ncon = mesh_vec[imesh]->NConnects();
        for(int64_t i=0; i<ncon; i++){
            TPZConnect &newnod = mesh_vec[imesh]->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(lagrange);
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildMixed4SpacesMultiPhysicsCompMesh(int order){
         
	TPZGeoMesh *gmesh = mGeometry;
	const bool isMHM = mSimData.mTNumerics.m_mhm_mixed_Q;
	
    
    // ========================================================
    // Creating atomic comp meshes
    TPZManVector<TPZCompMesh *, 7> mesh_vec(7);

    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order,0);
    mesh_vec[2] = DiscontinuousCmesh(0,0);
    mesh_vec[3] = DiscontinuousCmesh(0,0);
	
	if(isMHM){
		mesh_vec[4] = DiscontinuousCmesh(0,0);
		GroupConnectsBySubdomain(mesh_vec[4]);
		mesh_vec[5] = DiscontinuousCmesh(0,0);
		GroupConnectsBySubdomain(mesh_vec[5]);
		mesh_vec[6] = DiscontinuousCmesh(0,0);
		if(mSubdomainIndexGel.size() != gmesh->NElements()) DebugStop();
	}
	else{
		mesh_vec.resize(5);
		mesh_vec[4] = DiscontinuousCmesh(0,0);
	}

	
    // assign a subdomain to the lower level elements
	if (isMHM){
		IdentifySubdomainForLowdimensionElements(mesh_vec[0]);
	}

    if(isMHM && mSubdomainIndexGel.size() != gmesh->NElements()) DebugStop();

	if(isMHM && mSubdomainIndexGel.size()){
		VerifySubdomainIntegrity();
	}
    // ========================================================
    // Setting lagrange multiplier order
    SetLagrangeMultiplier4Spaces(mesh_vec);
    
    // ========================================================
    // Setting active spaces
    TPZManVector<int,7> active_approx_spaces(7);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
	active_approx_spaces[3] = 1;
	if (isMHM) {
		active_approx_spaces[4] = 1;
		active_approx_spaces[5] = 1;
		active_approx_spaces[6] = 0;
	}
	else{
		active_approx_spaces.resize(5);
		active_approx_spaces[4] = 0;
	}
    
    
    // ========================================================
    // Build multiphysics comp mesh
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    std::set<int> matsWithMem, matsWithOutMem;
    AddMultiphysicsMaterialsToCompMesh(order,matsWithMem, matsWithOutMem);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec,matsWithMem, matsWithOutMem );
    if(isFracSim()){
//        mMixedOperator->LoadReferences();
        this->InitializeMemoryFractureGlue();
    }
    
    {
        std::ofstream file("MixOpew.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(mMixedOperator, file);
    }
    
    if(isMHM && mSubdomainIndexGel.size() != gmesh->NElements()) DebugStop();
    
	if(isMHM && mSubdomainIndexGel.size()){
		VerifySubdomainIntegrity();
	}

	// ========================================================
    // Creates interface elements for hybridized intersections
    CreateIntersectionInterfaceElements();
    
    if(isMHM && mSubdomainIndexGel.size() != gmesh->NElements()) DebugStop();

    // ========================================================
    // Initialize all fractures properties (TMRSMemory)
    InitializeFracProperties(mMixedOperator);
    
#ifdef PZDEBUG
    {
//		ofstream out("mphysics.vtk");
//		TPZVTKGeoMesh::PrintCMeshVTK(mMixedOperator, out);
//        std::ofstream sout("mixed_cmesh_four_spaces.txt");
//        mMixedOperator->Print(sout);
    }
#endif

	cout << "\nNequations before hiding the elements = " << mMixedOperator->NEquations() << endl;
#ifdef PZDEBUG
	if (isMHM) {
        // visualize the subdomain association of the geometric element
//        ofstream out("gmeshdomain.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out,mSubdomainIndexGel);
    }
#endif
	// ========================================================
	// In case MHM, put the elements in submeshes
    // Verify the integrity of the subdomain indices
	if (isMHM) HideTheElements(mMixedOperator);
	
    // ========================================================
    // Condensing elements
	TPZReservoirTools::CondenseElements(mMixedOperator, 3, true);
	// TPZReservoirTools::PushConnectBackward(mMixedOperator, 3, 7);
	
#ifdef PZDEBUG
//    {
//        std::ofstream out("SubStructuredMesh.txt");
//        mMixedOperator->Print(out);
//    }
#endif
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildTransportMultiPhysicsCompMesh(){
    
    if (mSimData.mTNumerics.m_four_approx_spaces_Q) {
        BuildTransport4SpacesMultiPhysicsCompMesh();
    }else{
        BuildTransport2SpacesMultiPhysicsCompMesh();
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildTransport2SpacesMultiPhysicsCompMesh(){
    
    if (!mMixedOperator || !mGeometry) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,3> mixed_meshvec = mMixedOperator->MeshVector();
    TPZManVector<TPZCompMesh *,3> transport_meshvec(3);
    
    transport_meshvec[0] = mixed_meshvec[0];
    transport_meshvec[1] = mixed_meshvec[1];
    transport_meshvec[2] = TransportCmesh();
    
    
    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
//    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
        TPZTracerFlow * volume = nullptr;
    
    mTransportOperator->SetDefaultOrder(0);
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;

	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		std::cout << "physical name = " << material_name << std::endl;
		int materia_id = chunk.second;
		//            volume = new TMRSMultiphaseFlow<TMRSMemory>(materia_id,d);
		//            volume->SetDataTransfer(mSimData);
		volume = new  TPZTracerFlow(materia_id,dimension);
		mTransportOperator->InsertMaterialObject(volume);
	}
    
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
	for (auto& chunk : mSimData.mTBoundaryConditions.mBCTransportMatIdToTypeValue) {
		int bc_id   = chunk.first;
		std::pair<int,REAL>& typeAndVal = chunk.second;
		int bc_type = typeAndVal.first;
		val2[0]   = typeAndVal.second;
		TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
		mTransportOperator->InsertMaterialObject(face);
	}
    
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    {
//        TMRSMultiphaseFlow<TMRSMemory> * interface = new TMRSMultiphaseFlow<TMRSMemory>(transport_matid,dimension-1);
//        interface->SetDataTransfer(mSimData);
//        mTransportOperator->InsertMaterialObject(interface);
        
                TPZTracerFlow * interface = new TPZTracerFlow(transport_matid,dimension-1);
        //        interface->SetDataTransfer(mSimData);
                mTransportOperator->InsertMaterialObject(interface);
    }
    
    mTransportOperator->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 1;
    
    TPZMultiphysicsCompMesh *mult = dynamic_cast<TPZMultiphysicsCompMesh *>(mTransportOperator);
    if(mult){
        mult->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,transport_meshvec);
    }
    else{
        DebugStop();
    }
    //    mTransportOperator->BuildMultiphysicsSpace(active_approx_spaces,transport_meshvec);
    
    
    
    {
        mTransportOperator->Reference()->ResetReference();
        mTransportOperator->LoadReferences();
        
        TPZManVector<std::vector<int64_t>,4> cel_indexes(4);
        
        TPZManVector<int64_t,3> left_mesh_indexes(2,0);
        left_mesh_indexes[0] = 0;
        left_mesh_indexes[1] = 2;
        TPZManVector<int64_t,3> right_mesh_indexes(1,0);
        right_mesh_indexes[0] = 2;
        
        int64_t nel = mTransportOperator->NElements();
        for (int64_t el = 0; el < nel; el++) {
            
            TPZCompEl *cel = mTransportOperator->Element(el);
            if(!cel) DebugStop();
            TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!celmp) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            
            int gel_dim = gel->Dimension();
            cel_indexes[gel_dim].push_back(el);
            
        }
        
        for (auto cel_index: cel_indexes[dimension]) { // Higher dimension case
            TPZCompEl *cel = mTransportOperator->Element(cel_index);
            TPZMultiphysicsElement * celmult = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!celmult) {
                DebugStop();
            }
            
            if (!cel){continue;};
            
            TPZGeoEl *gel = cel->Reference();
            if (!gel){continue;};
            int nsides = gel->NSides();
            
            for (int iside = gel->NNodes(); iside < nsides; iside++) {
                
                TPZGeoElSide gelside(gel,iside);
                TPZCompElSide celside_l(cel,iside);
                TPZGeoElSide neig = gelside.Neighbour();
                bool condition =false;
                while (neig !=gelside && condition != true) {
                    
                    
					for (auto chunk : DomainDimNameAndPhysicalTag) {
						int material_id = chunk.second;
						if (neig.Element()->MaterialId() == material_id) {
							condition = true;
						}
						
					}
					
                    if (condition == false) {
						for(auto& chunk : mSimData.mTBoundaryConditions.mBCTransportMatIdToTypeValue){
							int bc_id   = chunk.first;
							if (neig.Element()->MaterialId() == bc_id) {
								condition = true;
							}
						}
                    }
                    
                    
                    if (condition == false) {
                        neig=neig.Neighbour();
                    }
                }
                
                TPZGeoEl *neihel = neig.Element();
                TPZCompElSide celside_r = neig.Reference();
                
                if ((neihel->Dimension() == gel->Dimension()) && (gel->Id() < neihel->Id()) ) {
                    TPZGeoElBC gbc(gelside,transport_matid);
                    
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), celside_l,celside_r);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                }
                if ((neihel->Dimension() == dimension - 1)) { // BC cases
                    
                    TPZGeoElBC gbc(gelside,neihel->MaterialId());
                    
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), celside_l,celside_r);
                    
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    
                }
                
            }
            
        }
    }
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    mTransportOperator->Print(transport);
#endif
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::BuildTransport4SpacesMultiPhysicsCompMesh(){
    
    if (!mMixedOperator || !mGeometry) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,5> mixed_meshvec = mMixedOperator->MeshVector();
    TPZManVector<TPZCompMesh *,5> transport_meshvec(5);
    
    transport_meshvec[0] = mixed_meshvec[0];
    transport_meshvec[1] = mixed_meshvec[1];
    transport_meshvec[2] = mixed_meshvec[2];
    transport_meshvec[3] = mixed_meshvec[3];
    transport_meshvec[4] = TransportCmesh();
    
    
    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    //    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
    TPZTracerFlow * volume = nullptr;
    mTransportOperator->SetDefaultOrder(0);
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    
	for (auto chunk : DomainDimNameAndPhysicalTag) {
		std::string material_name = chunk.first;
		std::cout << "physical name = " << material_name << std::endl;
		int material_id = chunk.second;
		//            volume = new TMRSMultiphaseFlow<TMRSMemory>(material_id,d);
		//            volume->SetDataTransfer(mSimData);
		
		volume = new TPZTracerFlow(material_id,dimension);
		//            volume->SetDataTransfer(mSimData);
		
		mTransportOperator->InsertMaterialObject(volume);
	}
    
    volume = new TPZTracerFlow(10,2);
    mTransportOperator->InsertMaterialObject(volume);
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
	for (auto& chunk : mSimData.mTBoundaryConditions.mBCTransportMatIdToTypeValue) {
		int bc_id   = chunk.first;
		std::pair<int,REAL>& typeAndVal = chunk.second;
		int bc_type = typeAndVal.first;
		val2[0]   = typeAndVal.second;
		TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
		mTransportOperator->InsertMaterialObject(face);
	}
    //crear controle no simdata
    //    TPZBndCond * face3 = volume->CreateBC(volume,-11,1,val1,val2);
    //    mTransportOperator->InsertMaterialObject(face3);
    
    int fracvol1ID = mSimData.mTGeometry.mInterface_material_idFracInf;
    int fracvol2ID = mSimData.mTGeometry.mInterface_material_idFracSup;
    int fracFracID = mSimData.mTGeometry.mInterface_material_idFracFrac;
    int fracbounId = mSimData.mTGeometry.mInterface_material_idFracBound;
    
    TPZBndCond * face = volume->CreateBC(volume,fracvol1ID,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face);
    
    TPZBndCond * face4 = volume->CreateBC(volume,fracvol2ID,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face4);
    
    TPZBndCond * face2 = volume->CreateBC(volume,fracFracID,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face2);
    //
    TPZBndCond * face5 = volume->CreateBC(volume,fracbounId,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face5);
    
    
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    {
        //        TMRSMultiphaseFlow<TMRSMemory> * interface = new TMRSMultiphaseFlow<TMRSMemory>(transport_matid,dimension-1);
        TPZTracerFlow * interface = new TPZTracerFlow (transport_matid,dimension-1);
        //        interface->SetDataTransfer(mSimData);
        mTransportOperator->InsertMaterialObject(interface);
        
    }
    
    mTransportOperator->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(5); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 0;
    active_approx_spaces[3] = 0;
    active_approx_spaces[4] = 1;
    //    mTransportOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,transport_meshvec);
    
    TPZMultiphysicsCompMesh *mult = dynamic_cast<TPZMultiphysicsCompMesh *>(mTransportOperator);
    if(mult){
        mult->BuildMultiphysicsSpace(active_approx_spaces,transport_meshvec);
    }
    else{
        DebugStop();
    }
    
    
    
    {
        mTransportOperator->Reference()->ResetReference();
        mTransportOperator->LoadReferences();
        CreateInterfaces(mTransportOperator);
    }
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    mTransportOperator->Print(transport);
#endif
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetDataTransfer(TMRSDataTransfer & SimData){
    mSimData = SimData;
    
    if (mSimData.mTGeometry.mGmeshFileName=="") {
        std::string geoname = "PreProcess/meshes/"+mSimData.mSimulationName + "_nLayers_"+ std::to_string(mSimData.mTGeometry.mnLayers)  +"_nRef_"+std::to_string(mSimData.mTGeometry.mnref)+".txt" ;
        mSimData.mTGeometry.mGmeshFileName = geoname;
    }
    
    std::ifstream file(mSimData.mTGeometry.mGmeshFileName);
    
    
//    if(!mGeometry){
//        if (file) {
//            std::cout<<"The geometric mesh will be loaded from the " + mSimData.mTGeometry.mGmeshFileName + " file."<<std::endl;
//            std::string filename = mSimData.mTGeometry.mGmeshFileName;
//            TPZPersistenceManager::OpenRead(filename);
//            TPZSavable *restore = TPZPersistenceManager::ReadFromFile();
//            mGeometry = dynamic_cast<TPZGeoMesh *>(restore);
//            PrintGeometry(mSimData.mSimulationName,1,0);
//            std::cout<<"The geometric mesh has been loaded successfully."<<std::endl;
//        }
//        else{
//            std::cout<<" Geometric mesh information has not been entered. Please enter the mesh in a text file or set it in the approximation space object."<<std::endl;
//            DebugStop();
//        }
//        
//    }
    
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer & TMRSApproxSpaceGenerator::GetDataTransfer(){
    return mSimData;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::LinkMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZCompMesh * TransportOperator){
    
	DebugStop(); // Is this function deprecated? If not, it needs serious fixing
	
	// Commented on May 2022
//    TPZMultiphysicsCompMesh *mult = dynamic_cast<TPZMultiphysicsCompMesh *>(TransportOperator);
//    if(!mult){
//        DebugStop();
//    }
//    AdjustMemory(MixedOperator, mult);
//    for (auto item : mSimData.mTGeometry.mDomainNameAndMatId[2]) {
//        int material_id = item.second;
//        UnifyMaterialMemory(material_id, MixedOperator, mult);
//        FillMaterialMemory(material_id, MixedOperator, mult);
//    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::AdjustMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator){
    
    if (!MixedOperator || !TransportOperator) {
        DebugStop();
    }
    
    /// Adjust integration rule
    /// o Stands for reservoir
    /// d Stands for pressure
    
    TPZCompMesh * cmesh_res = MixedOperator;
    TPZCompMesh * cmesh_tra = TransportOperator;
    
    cmesh_tra->LoadReferences();
    int nel_res = cmesh_res->NElements();
    int gmesh_dim = cmesh_tra->Reference()->Dimension();
    
    // Scanning structure
    std::vector<std::pair<int64_t, int64_t>> cel_pairs;
    for (long el = 0; el < nel_res; el++) {
        TPZCompEl *cel_res = cmesh_res->Element(el);
        if (!cel_res) {
            continue;
        }
        
        TPZGeoEl * gel = cel_res->Reference();
        if (!gel) {
            continue;
        }
        
        if (gel->Dimension() != gmesh_dim || gel->HasSubElement()) {
            continue;
        }
        
        /// Finding the other computational element
        TPZCompEl * cel_tra = gel->Reference();
        if (!cel_tra) {
            continue;
        }
        
        int64_t cel_res_index = cel_res->Index();
        int64_t cel_tra_index = cel_tra->Index();
        cel_pairs.push_back(std::make_pair(cel_res_index, cel_tra_index));
        cel_tra->SetFreeIntPtIndices();  // operation involving resize. It is not thread safe.
    }
    
    int nel = cel_pairs.size();
#ifdef USING_TBB
    tbb::parallel_for(size_t(0), size_t(nel), size_t(1), [&cel_pairs,&cmesh_res,&cmesh_tra] (size_t & i)
                      {
                          int64_t cel_res_index = cel_pairs[i].first;
                          int64_t cel_geo_index = cel_pairs[i].second;
                          TPZCompEl *cel_res = cmesh_res->Element(cel_res_index);
                          TPZCompEl * cel_tra = cmesh_tra->Element(cel_geo_index);
                          
                          const TPZIntPoints & rule = cel_res->GetIntegrationRule();
                          TPZIntPoints * cloned_rule = rule.Clone();
                          TPZManVector<int64_t,20> indices;
                          cel_res->GetMemoryIndices(indices);
                          cel_tra->SetMemoryIndices(indices);
                          cel_tra->SetIntegrationRule(cloned_rule);
                      }
                      );
    
#else
    for (long i = 0; i < nel; i++) {
        
        int64_t cel_res_index = cel_pairs[i].first;
        int64_t cel_tra_index = cel_pairs[i].second;
        TPZCompEl *cel_res = cmesh_res->Element(cel_res_index);
        TPZCompEl * cel_tra = cmesh_tra->Element(cel_tra_index);
        
        const TPZIntPoints & rule = cel_res->GetIntegrationRule();
        TPZIntPoints * cloned_rule = rule.Clone();
        TPZManVector<int64_t,20> indices;
        cel_res->GetMemoryIndices(indices);
        cel_tra->SetMemoryIndices(indices);
        cel_tra->SetIntegrationRule(cloned_rule);
    }
#endif
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::UnifyMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator) {
    
    if (!MixedOperator || !TransportOperator) {
        DebugStop();
    }
    
    TPZCompMesh * cmesh_o = MixedOperator;
    TPZCompMesh * cmesh_d = TransportOperator;
    
    /// Apply memory link
    TPZMaterial * material_o = cmesh_o->FindMaterial(material_id);
    TPZMaterial * material_d = cmesh_d->FindMaterial(material_id);
    if (!material_o || !material_d) {
        DebugStop();
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory_o = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material_o);
    TPZMatWithMem<TMRSMemory> * mat_with_memory_d = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material_d);
    if (!mat_with_memory_o || !mat_with_memory_d) {
        DebugStop();
    }
    
    mat_with_memory_d->SetMemory(mat_with_memory_o->GetMemory());
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::FillMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator){
    
    if (!MixedOperator || !TransportOperator) {
        DebugStop();
    }
    
    TPZCompMesh * cmesh = MixedOperator;
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    if (!material) {
        DebugStop();
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
    if (!mat_with_memory) {
        DebugStop();
    }
    
    
    // Set initial porosity, permeability, saturations, etc ...
    {
        std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        int ndata = memory_vector->NElements();
        
#ifdef USING_TBB
        tbb::parallel_for(size_t(0), size_t(ndata), size_t(1), [&memory_vector] (size_t & i) {
            TMRSMemory &mem = memory_vector.get()->operator [](i);
            
            
            mem.m_sw = 0.0;
            mem.m_phi = 1.0;
            
            REAL kappa = 1.0;
            mem.m_kappa.Resize(3, 3);
            mem.m_kappa.Zero();
            mem.m_kappa_inv.Resize(3, 3);
            mem.m_kappa_inv.Zero();
            for (int i = 0; i < 3; i++) {
                mem.m_kappa(i,i) = kappa;
                mem.m_kappa_inv(i,i) = 1.0/kappa;
            }
            
        }
                          );
        
#else
        for (int i = 0; i < ndata; i++) {
            TMRSMemory &mem = memory_vector.get()->operator [](i);
            mem.m_sw = 0.0;
            mem.m_phi = 1.0;
            REAL kappa = 1.0;
            mem.m_kappa.Resize(3, 3);
            mem.m_kappa.Zero();
            mem.m_kappa_inv.Resize(3, 3);
            mem.m_kappa_inv.Zero();
            //            kappa *= rand() % 40 +1;
            for (int j = 0; j < 3; j++) {
                mem.m_kappa(j,j) = kappa;
                mem.m_kappa_inv(j,j) = 1.0/kappa;
            }
        }
        
        
        //        int nels = cmesh->NElements();
        //        for (int iel = 0; iel< nels; iel++) {
        //            TPZCompEl *cel = cmesh->Element(iel);
        //            if(!cel){
        //                continue;
        //
        //            }
        //            TPZGeoEl *gel = cel->Reference();
        //            if (!gel || gel->HasSubElement()) {
        //                continue;
        //            }
        //
        //            if (cel->Dimension()!= cmesh->Dimension()) {
        //                continue;
        //            }
        //            if (!MixedOperator->Element(iel)) {
        //                continue;
        //            }
        //            TPZVec<int64_t> indices;
        //            cel->GetMemoryIndices(indices);
        //            TPZVec<REAL> qsi(3,0.25);
        //            qsi[2]=0.0;
        //            TPZVec<REAL> point(3,0.0);
        //            gel->X(qsi, point);
        ////            if (gel->MaterialId()!=2){
        ////                continue;
        ////            }
        //
        ////            int val = rand() % 100;
        //
        //
        //            REAL kappa =  1000*(sin(point[0])*sin(point[1]) + 2);
        //
        //
        //
        //
        ////            REAL kappa = 100000.0 + 1*(sin(point[0])*sin(point[1])+2);
        //
        //
        //            for (auto memIndex: indices) {
        //                if (memIndex<=0) {
        //                    continue;
        //                }
        //
        //
        //                TMRSMemory &mem = memory_vector.get()->operator [](memIndex);
        //                mem.m_sw = 0.0;
        //                mem.m_phi = 0.1;
        //                mem.m_kappa.Resize(3, 3);
        //                mem.m_kappa.Zero();
        //                mem.m_kappa_inv.Resize(3, 3);
        //                mem.m_kappa_inv.Zero();
        //                for (int j = 0; j < 3; j++) {
        //                    mem.m_kappa(j,j) = kappa;
        //                    mem.m_kappa_inv(j,j) = 1.0/kappa;
        //                }
        //            }
        //
        //
        //        }
        //
        
        
        
#endif
        
        
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::FillMaterialMemoryDarcy(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZAlgebraicTransport *algebraicTranspor){
    
    if (!MixedOperator || !algebraicTranspor) {
        DebugStop();
    }
    
    TPZCompMesh * cmesh = MixedOperator;
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    if (!material) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZDarcyMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPZDarcyMemory> * >(material);
    if (!mat_with_memory) {
        DebugStop();
    }
    int nels = cmesh->NElements();
    for (int iel =0; iel<=nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if (!cel) {
            continue;
        }
        TPZManVector<int64_t,20> indices;
        cel->GetMemoryIndices(indices);
        
        
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetUpdateMaterialMemory(int material_id, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q){
    
    if (!cmesh) {
        DebugStop();
    }
    
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    if (!material) {
        DebugStop();
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
    if (mat_with_memory) {
        mat_with_memory->SetUpdateMem(update_memory_Q);
        return;
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory_trans = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
    if (mat_with_memory_trans) {
        mat_with_memory_trans->SetUpdateMem(update_memory_Q);
        return;
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetUpdateMemory(int dimension, TMRSDataTransfer & sim_data, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q){
	if (dimension != 3) {
		DebugStop(); //mDomainNameAndMatId only knows about 3D materials
	}
    for (auto item : sim_data.mTGeometry.mDomainNameAndMatId) {
        int material_id = item.second;
        SetUpdateMaterialMemory(material_id, cmesh, update_memory_Q);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
{
    //    {
    //        std::ofstream out("gmeshref.txt");
    //        gmesh->Print(out);
    //    }
    coarseindices.Resize(gmesh->NElements());
    int count = 0;
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->Father()) continue;
        coarseindices[count] = el;
        count++;
    }
    coarseindices.Resize(count);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::InitializeFracProperties(TPZMultiphysicsCompMesh * MixedOperator)
{
    //
    TPZCompMesh * cmesh = MixedOperator;
    if (!MixedOperator) {
        DebugStop();
    }
    //
    int n_mats = mSimData.mTReservoirProperties.m_permeabilitiesbyId.size();
    for (auto iter = mSimData.mTReservoirProperties.m_permeabilitiesbyId.begin(); iter != mSimData.mTReservoirProperties.m_permeabilitiesbyId.end(); ++iter){
     
        int matId = iter->first;
        REAL kappa = iter->second;
        
        TPZMaterial * material1 = cmesh->FindMaterial(matId); ;
        TPZMatWithMem<TMRSMemory> * mat_with_memory1 = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material1);
        
        
        // Set initial porosity, permeability, saturations, etc ...
        {
            std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector1 = mat_with_memory1->GetMemory();
            int ndata1 = memory_vector1->NElements();
            
            for (int i = 0; i < ndata1; i++) {
                TMRSMemory &mem = memory_vector1.get()->operator [](i);
                mem.m_kappa.Resize(3, 3);
                mem.m_kappa.Zero();
                mem.m_kappa_inv.Resize(3, 3);
                mem.m_kappa_inv.Zero();
                for (int j = 0; j < 3; j++) {
                    mem.m_kappa(j,j) = kappa;
                    mem.m_kappa_inv(j,j) = 1.0/kappa;
                }
            }
        }
    }
    
    
//    TPZMaterial * material = cmesh->FindMaterial(FractureMatId()); //matIdFractures;
//    if (!material) {
//        return;
//    }
    TPZMaterial * material = nullptr;
	map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
    for (it = mSimData.mTFracProperties.m_fracprops.begin(); it != mSimData.mTFracProperties.m_fracprops.end(); it++)
    {
        int matFrac = it->first;
        material = cmesh->FindMaterial(matFrac); //matIdFractures;
        if (!material) {
            DebugStop();
        }
        TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
        if (!mat_with_memory) {
            DebugStop();
        }
        
        // Set initial porosity, permeability, saturations, etc ...
        {
            std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
            int ndata = memory_vector->NElements();
            
            for (int i = 0; i < ndata; i++) {
                TMRSMemory &mem = memory_vector.get()->operator [](i);
                mem.m_sw = 0.0;
				mem.m_phi = it->second.m_porosity;
				const REAL fracwidth = it->second.m_width;
                mem.m_ad = fracwidth;
                REAL kappa = it->second.m_perm;
                mem.m_kappa.Resize(3, 3);
                mem.m_kappa.Zero();
                mem.m_kappa_inv.Resize(3, 3);
                mem.m_kappa_inv.Zero();
                for (int j = 0; j < 3; j++) {
                    mem.m_kappa(j,j) = kappa;
                    mem.m_kappa_inv(j,j) = 1.0/kappa;
                }
            }
        }
        
    }
    
//    TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
//    if (!mat_with_memory) {
//        DebugStop();
//    }
//
//    // Set initial porosity, permeability, saturations, etc ...
//    {
//        std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
//        int ndata = memory_vector->NElements();
//
//        for (int i = 0; i < ndata; i++) {
//            TMRSMemory &mem = memory_vector.get()->operator [](i);
//            mem.m_sw = 0.0;
//            mem.m_phi = 0.1;
//            REAL kappa = mSimData.mTFracProperties.m_Permeability;
//            mem.m_kappa.Resize(3, 3);
//            mem.m_kappa.Zero();
//            mem.m_kappa_inv.Resize(3, 3);
//            mem.m_kappa_inv.Zero();
//            for (int j = 0; j < 3; j++) {
//                mem.m_kappa(j,j) = kappa;
//                mem.m_kappa_inv(j,j) = 1.0/kappa;
//            }
//        }
//    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::findNeighElementbyMatId(TPZGeoElSide &gelside, std::vector<TPZGeoElSide> &neihside, std::set<int> VolMatIds){
    
    int verif =0;
    TPZGeoElSide NeihSideAux = gelside.Neighbour();
    while (gelside != NeihSideAux) {
        for (auto Mat_id: VolMatIds) {
            if (NeihSideAux.Element()->MaterialId() == Mat_id) {
                neihside.push_back(NeihSideAux);
                verif=1;
                break;
            }
        }
        NeihSideAux = NeihSideAux.Neighbour();
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateElementInterfaces(TPZGeoEl *gel){
    
    const int domaindim = mGeometry->Dimension();
    int dimension = gel->Dimension();
    int nsides = gel->NSides();
    int ncoord = gel->NCornerNodes();
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
//    auto frac = mSimData.mTGeometry.mDomainFracNameAndMatId;
//    int fracmat = frac["Fractures"];
    std::set<int> FracMatID, fracbcs;
    GetMaterialIds(dimension-1, FracMatID, fracbcs);
    
    int val =0;
    if (dimension==3) {
        val = gel->NSides(1);
    }
    std::set<int> VolMatIds;
    std::set<int> Boundaries;
    GetMaterialIds(dimension, VolMatIds, Boundaries);
    
    
    for (int iside = ncoord+val; iside < nsides-1; iside++) {
        TPZGeoElSide gelside(gel,iside);
        std::vector<TPZGeoElSide> gelneigVec;
        //Volumetric Elements
        TPZCompElSide celside_l=gelside.Reference();
        findNeighElementbyMatId(gelside,gelneigVec,VolMatIds);
        if(gelneigVec.size()==1){
            TPZGeoElSide gelneig =gelneigVec[0];
            std::vector <TPZGeoElSide >gelfracVec;
            findNeighElementbyMatId(gelside, gelfracVec, FracMatID);
            if (gelfracVec.size()==0) {
                if (gel->Id() < gelneig.Element()->Id()) {
                    TPZCompElSide celside_r = gelneig.Reference();
                    TPZGeoElBC gbc(gelside,transport_matid);
                    TPZInterfaceElement *mp_interface_el = new TPZInterfaceElement(*mTransportOperator, gbc.CreatedElement(), celside_l, celside_r);
                }
            }
            else{
                //Interface Vol-Frac
//                std::cout<<"Not created Vol-Vol (Frac-Vol -> ok)"<<std::endl;
            }
            
        }
        else if(gelneigVec.size()==0){ //BoundaryElements
            findNeighElementbyMatId(gelside,gelneigVec,Boundaries);
            if (gelneigVec.size() != 1)
                DebugStop(); // if has no volume neighbors, it HAS to have a bc neighbor
            TPZGeoElSide gelneig = gelneigVec[0];
            TPZCompElSide celside_r = gelneig.Reference();
            int matId = gelneig.Element()->MaterialId();
            TPZGeoElBC gbc(gelside,matId);
            TPZInterfaceElement *mp_interface_el = new TPZInterfaceElement(*mTransportOperator, gbc.CreatedElement(), celside_l, celside_r);
        }
        else
            DebugStop(); // Can only have one or zero volume neighbor
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateInterfaces(TPZCompMesh *cmesh){
    
    // Here the interfaces are created based on the previously set material ids
    
    int nels = cmesh->NElements();
    int dim = cmesh->Dimension();
    std::set<int> matidsfrac;
    std::set<int> bcmatids;
    GetMaterialIds( dim-1,matidsfrac, bcmatids);
    
    // cel_indexes is a vector of vectors that stores the elements by dimension. It also skips
    // elements such as BC elements, so it only selects fracture elements for dim-1. For dim, it selects everything.
    TPZManVector<std::vector<int64_t>,4> cel_indexes(4);
    for (int64_t el = 0; el < nels; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        
        int gel_dim = gel->Dimension();
        if (gel_dim == dim-1) {
            int gelid = gel->MaterialId();
            if (!IsFracMatId(gelid)) {
                continue;
            }
        }
        cel_indexes[gel_dim].push_back(el);
    }
    //Creating Interfaces vol-vol and vol-bound without vol-fractures.
    for(auto elindex: cel_indexes[dim]){ // loop over volume elements
        TPZCompEl *cel = cmesh->Element(elindex);
        TPZGeoEl *gel = cel->Reference();
        CreateElementInterfaces(gel);
    }
    //Creating Interfaces frac-vol frac-frac and frac-boundary
    for(auto elindex: cel_indexes[dim-1]){ // loop over fracture elements
        TPZCompEl *cel = cmesh->Element(elindex);
        TPZGeoEl *gel = cel->Reference();
        CreateFracInterfaces(gel);
    }
    
    //    std::ofstream file("NewInterfaces.vtk");
    //    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void TMRSApproxSpaceGenerator::CreateFracInterfaces(TPZGeoEl *gel){
    
    int dimension = gel->Dimension();
    int nsides = gel->NSides();
    int ncoord = gel->NCornerNodes();
    int transport_interface_matid = mSimData.mTGeometry.mInterface_material_id;
    int matIdFracSup = mSimData.mTGeometry.mInterface_material_idFracSup;
    int matIdFracInf = mSimData.mTGeometry.mInterface_material_idFracInf;
    int matIdFracFrac = mSimData.mTGeometry.mInterface_material_idFracFrac;
    int matIdFracBound = mSimData.mTGeometry.mInterface_material_idFracBound;

    std::set<int> FracNeihVec;
    int dim = gel->Mesh()->Dimension();
    std::set<int> matidsbcfrac;
    GetMaterialIds( dim-1, FracNeihVec, matidsbcfrac);
    
    
    for (auto& chunk : mSimData.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
        int bc_id   = chunk.first;
        matidsbcfrac.insert(bc_id);
    }


    //intersecmatids
    std::set<int> FracNeihVec2;
    for(auto ifrac:FracNeihVec ){
        int matfracintersec =ifrac+2;
        FracNeihVec2.insert(matfracintersec);
    }
    
    std::set<int> FracIntersec;
    FracIntersec.insert(mSimData.mTGeometry.m_pressureMatId);
    
    
    // ---------------> Adding volume materials
    std::map<std::string,int> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainNameAndMatId;
    std::set<int> VolMatIds;
    int dimen = gel->Mesh()->Dimension();
    
    for (auto chunk : DomainDimNameAndPhysicalTag) {
        int material_id = chunk.second;
        VolMatIds.insert(material_id);
    }
        
   
    for (int iside = ncoord; iside < nsides; iside++) {
        TPZGeoElSide gelside(gel,iside);
        std::vector<TPZGeoElSide> gelneigVec;
        TPZCompElSide celside_l=gelside.Reference();
        // Create interface between fracture/volume elements
        // One has to be sup and the other inf. This decision is done randomly based on the order of the neighbors
        if(iside == nsides-1 ){
            findNeighElementbyMatId(gelside,gelneigVec,VolMatIds);
            if(gelneigVec.size()==2){
                for (int ivol = 0; ivol<2; ivol++) {
                    TPZGeoElSide gelneig =gelneigVec[ivol];
                    int matid = mSimData.mTGeometry.mInterface_material_idFracSup;
                    if (ivol == 0) {
                        matid = mSimData.mTGeometry.mInterface_material_idFracInf;
                    }
                    CreateInterfaceElements(gelneig, gelside, matid);
                   
                }
            }
        }
        else{ // Create interface between fracture/fracture and fracture/boundary
//            std::set<int> FracNeihVec;
            //Verify if is frac-frac frac-bound or frac-intersec
            int matfrac = gelside.Element()->MaterialId();
            auto fracprop = mSimData.mTFracProperties.m_fracprops[matfrac];
            std::set<int> matsbcfrac =fracprop.m_fracbc;
            std::set<int> myFracId;
            std::set<int> myFracIntersectId;
            std::set<int> intersecId;
            intersecId.insert( mSimData.mTGeometry.m_pressureMatId );
            myFracId.insert(matfrac);
            int matIntersec = fracprop.m_fracIntersectMatID;
            myFracIntersectId.insert(matIntersec);
            
            bool is_fracfrac=false;
            bool is_fracbound=false;
            bool is_fracintersect=false;
            
            
            std::vector<TPZGeoElSide> gelneigMyFrac;
            std::vector<TPZGeoElSide> gelneigMyBCFrac;
            std::vector<TPZGeoElSide> gelneigMyIntersecFrac;
            
            findNeighElementbyMatId(gelside,gelneigMyFrac,myFracId);
            findNeighElementbyMatId(gelside,gelneigMyBCFrac,matsbcfrac);
            findNeighElementbyMatId(gelside,gelneigMyIntersecFrac,myFracIntersectId);
            
            
            
            if (gelneigMyFrac.size()==1 && gelneigMyIntersecFrac.size() ==0) {
                TPZGeoElSide gelneig =gelneigMyFrac[0];
                if (gel->Id() < gelneig.Element()->Id()) {
                    int matid = matIdFracFrac;
                    CreateInterfaceElements(gelside, gelneig, matid);
                }
                is_fracfrac=true;
            }
            if (gelneigMyBCFrac.size()==1) {
                int matid = -1;
                for(auto id:matsbcfrac){
                    matid=id;
                }
                CreateInterfaceElements(gelside, gelneigMyBCFrac[0], matid);
                is_fracbound=true;
            }
                
            if (gelneigMyIntersecFrac.size()>=1) {
                std::vector<TPZGeoElSide> realIntersect;
                findNeighElementbyMatId(gelside,realIntersect,intersecId);
                if(!realIntersect.size()){
                    DebugStop();
                }
                else{
                    int matid=100;
                    CreateInterfaceElements(gelside, realIntersect[0], matid);
                }
                is_fracintersect=true;
            }
            if (!is_fracfrac && !is_fracbound && !is_fracintersect) {
                std::vector<TPZGeoElSide> realBC;
                int dim = gel->Dimension()+1;
                std::set<int>matids;
                std::set<int> bcmatids;
                GetMaterialIds(dim, matids, bcmatids);
                findNeighElementbyMatId(gelside,realBC,bcmatids);
                TPZGeoElSide gelsidecorrect;
                bool found = false;
                for(auto gelside2: realBC){
                    int dimel=gelside2.Element()->Dimension();
                    if(dimel== gel->Dimension()-1){
                        found=true;
                        CreateInterfaceElements(gelside, gelside2, gelside2.Element()->MaterialId());
                        break;
                    }
                }
                if(!found){
                    DebugStop();
                }
               
            }

        }
    }
}
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateInterfaceElements(TPZGeoElSide &gelside, TPZGeoElSide &gelneig, int matid){
    
    TPZGeoElBC gbc(gelside,matid);
    TPZGeoEl *gel = gelside.Element();
    TPZCompElSide celside_l = gelside.Reference();
    TPZCompElSide celside_r = gelneig.Reference();
    
    if(gel->Id() < gelneig.Element()->Id()){
        TPZInterfaceElement *mp_interface_el = new TPZInterfaceElement(*mTransportOperator, gbc.CreatedElement(), celside_l,celside_r);
    }
    else{
        TPZInterfaceElement *mp_interface_el = new TPZInterfaceElement(*mTransportOperator, gbc.CreatedElement(),celside_r, celside_l);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC){
    int dimension = gel->Dimension();
    int matid = gel->MaterialId();
 
    cmesh->SetDimModel(dimension);
    if (is_BC) {
        cmesh->SetDimModel(dimension+1);
    }
    TPZCompEl * cel = cmesh->ApproxSpace().CreateCompEl(gel, *cmesh);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
    TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
    if (intel){
        intel->PRefine(p_order);
    } else if (intelDisc) {
        intelDisc->SetDegree(p_order);
        intelDisc->SetTrueUseQsiEta();
    } else {
        DebugStop();
    }
    gel->ResetReference();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::HideTheElements(TPZCompMesh *cmesh){
    int KeepOneLagrangian = 6;
    cmesh->ComputeNodElCon();

    
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    //    TPZGeoMesh *gmesh = fCMesh->Reference();
    //    gmesh->ResetReference();
    if(mSubdomainIndexGel.size() != gmesh->NElements()) DebugStop();
    
    int64_t nel = mSubdomainIndexGel.size();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel){
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if(!cel){continue;}
        int indexel = cel->Index();
        
        //        TPZCompEl *cel2 = cmesh->ElementVec()[indexel];
        //        if(!cel2){
        //            DebugStop();
        //        }
        int64_t domain = mSubdomainIndexGel[el];
        if (domain == -1) {
//            std::cout<<"matId of element not condensed "<<gel->MaterialId()<<std::endl;
            continue;
        }
        ElementGroups[domain].insert(indexel);
    }
    
    if (ElementGroups.size() <= 0) // Change this 0 to another number to debug
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<int64_t,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << cmesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
//    {
//        std::map<int64_t,std::set<int64_t>> connecttoel;
//        int64_t nel = cmesh->NElements();
//        for (int64_t el = 0; el < nel; el++) {
//            TPZCompEl *cel = cmesh->Element(el);
//            int nc = cel->NConnects();
//            for (int ic = 0; ic<nc; ic++) {
//                int64_t cindex = cel->ConnectIndex(ic);
//                connecttoel[cindex].insert(el);
//            }
//        }
//        auto els = connecttoel[667];
//        for(auto it:els)
//        {
//            TPZCompEl *cel = cmesh->Element(it);
//            cel->Print();
//        }
//        TPZConnect &c = cmesh->ConnectVec()[667];
//        c.Print(*cmesh);
//    }
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(cmesh, ElementGroups, submeshindices, KeepOneLagrangian);
    
//    std::cout << "After putting in substructures\n";
    //    fMHMtoSubCMesh = submeshindices;
    cmesh->ComputeNodElCon();
    {
        int64_t nc = cmesh->NConnects();
        for (int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = cmesh->ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency())
            {
                c.RemoveDepend();
            }
        }
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    
    for(auto it : submeshindices)
    {
        TPZCompEl *cel = cmesh->Element(it.second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);

        const int nthreads = mSimData.mTNumerics.m_nThreadsMixedProblem;
        const bool isUseSparse = true;
        if(isUseSparse){
            // Pardiso does pivoting which is necessary in some problems
            subcmesh->SetAnalysisSparse(nthreads);
        }
        else{
            // Skyline does not do pivoting and unfortunately, it is necessary for certain problems.
            // So, using it might lead to zero pivot and therefore it is strongly not recommended.
            subcmesh->SetAnalysisSkyline(nthreads, 0);
        }
    }
    //    GroupandCondenseElements();
    //    GroupandCondenseElementsEigen();
    
    std::cout << "\n\t=======> Finished substructuring\n";
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

const bool TMRSApproxSpaceGenerator::isThereFracIntersection() const {
    const int matIdIntersection = mSimData.mTFracIntersectProperties.m_IntersectionId;
    for (auto gel : mGeometry->ElementVec()){
        const int matid = gel->MaterialId();
        if (matid == matIdIntersection) { // intersection matid
            return true;
        }
    }
    
    return false;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::HybridizeIntersections(TPZVec<TPZCompMesh *>& meshvec_Hybrid) {
    
    if (!mHybridizer){
        DebugStop();
    }
    cout << "\n==> Starting hybridizing intersections..." << endl;
        
    const int matIdIntersection = mSimData.mTFracIntersectProperties.m_IntersectionId;
    TPZCompMesh* fluxmesh = meshvec_Hybrid[0];

    // Need to decrease dimension otherwise 1d elements are not created
    const int dim = fluxmesh->Reference()->Dimension();
    fluxmesh->SetDimModel(2);
    fluxmesh->SetAllCreateFunctionsHDiv();

    
    TPZGeoMesh* gmesh = fluxmesh->Reference();
    fluxmesh->LoadReferences();
    // what does this do??? In how many spaces do you insert material objects???
    mHybridizer->InsertPeriferalMaterialObjects(meshvec_Hybrid);
    
    int dimfrac = 2;
    for (auto gel : gmesh->ElementVec()) {
        if(!gel || gel->HasSubElement()) continue;
        const int gelmatid = gel->MaterialId();
        if (gelmatid != matIdIntersection) {
            continue;
        }
        if (gel->Dimension() != 1) {
            DebugStop();
        }
        const bool isIntersectEnd = false; // this was used to set pressure at an intersection end. TODO: Delete?
        // Search for first neighbor that that is domain
        TPZGeoElSide gelside(gel);
#ifdef PZDEBUG
        TPZGeoElSide test = gelside.Neighbour().HasNeighbour(matIdIntersection);
        if(test && test != gelside){
            // Why are there two intersection at the same location?!?!
            DebugStop();
            
            // This could be called instead of DebugStop to erase the duplicates
            TPZGeoEl* geltest = test.Element();
            const int64_t duplicateIndex = geltest->Index();
            geltest->RemoveConnectivities();
            delete geltest;
            gmesh->ElementVec()[duplicateIndex] = nullptr;
        }
#endif
        TPZGeoElSide neigh = gelside.Neighbour();
        
        // loop over the neighbours of the 1d intersection element
        while(neigh != gelside){
            TPZGeoEl* gelneigh = neigh.Element();
            int neighmatid = gelneigh->MaterialId();
            int neighdim = gelneigh->Dimension();
            // is there only one matidfrac????
//            if (neighmatid == matidfrac && neighdim == dimfrac) {
            if ( IsFracMatId(neighmatid) && neighdim == dimfrac) {
#ifdef PZDEBUG
//                cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " is an intersection element" << endl;
//                cout << "===> Trying to split the connects of the flux mesh and create pressure element..." << endl;
#endif
                // the neighbour has to be a fracture element
                TPZCompEl* celneigh = gelneigh->Reference();
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (celneigh);
                if (!intel)
                    DebugStop();
                
                const int side = neigh.Side();
                // same as neigh.Reference()
                TPZCompElSide celsideleft(intel, side);
                // why pass three redundant parameters : celsideleft, intel, side?
                bool isNewInterface = mHybridizer->HybridizeInterface(celsideleft,intel,side,meshvec_Hybrid,isIntersectEnd);
                if (isNewInterface) {
#ifdef PZDEBUG
//                    cout << "=====> Connects splitted succesfuly!" << endl;
#endif
                    // break out of the while loop
                    break;
                }
                else{
                    // just make sure an interface is split only once??
                    DebugStop();
                }
            }
            neigh = neigh.Neighbour();
        } // while
    }
    
    // Set createApproxSpace createCompEl functions back to domain max dimension
    fluxmesh->SetDimModel(dim);
    fluxmesh->SetAllCreateFunctionsHDiv();
    
    cout << "==> Finished hybridizing intersections..." << endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateIntersectionInterfaceElements(TPZVec<TPZCompMesh *>& meshvec_Hybrid) {
    TPZCompMesh* cmeshpressure = mMixedOperator->MeshVector()[1];
    mMixedOperator->Reference()->ResetReference();
    cmeshpressure->LoadReferences();
    const int lagrangematid = mHybridizer->lagrangeInterfaceMatId();
    const int lagrangematidend = mHybridizer->lagrangeInterfaceEndMatId();
    TPZStack<int64_t> pressureindices;
    for (auto cel : cmeshpressure->ElementVec()) {
        if(!cel) continue;
        const int celmatid = cel->Material()->Id();
        if (celmatid != lagrangematid && celmatid != lagrangematidend) {
            continue;
        }
        TPZGeoEl* gel = cel->Reference();
        pressureindices.Push(gel->Index());
        mHybridizer->CreateInterfaceElementsForGeoEl(mMixedOperator, meshvec_Hybrid, gel);
    }
    // identify the domain indices of the interface elements
    SetInterfaceDomains(pressureindices,mHybridizer->fInterfaceMatid);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CreateIntersectionInterfaceElements() {
    // identify the domain indices of the interface elements
    int interfacematid = mSimData.mTGeometry.mInterface_material_id;
    TPZLagrangeMultiplierCS<> * mat = new TPZLagrangeMultiplierCS<>(interfacematid, 1);
    mMixedOperator->InsertMaterialObject(mat);
    std::set<int> fracintersectmatid;
    for(auto &frac : mSimData.mTFracProperties.m_fracprops)
    {
        fracintersectmatid.insert(frac.second.m_fracIntersectMatID);
    }
    int pressurematid = mSimData.mTGeometry.m_pressureMatId;
    int64_t nelem = mGeometry->NElements();
    for (int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = mGeometry->Element(el);
        if(!gel) continue;
        int matid = gel->MaterialId();
        if(fracintersectmatid.find(matid) == fracintersectmatid.end()) continue;
        TPZGeoElSide gelside(gel);
        TPZCompElSide celside(gelside.Reference());
        if(!celside) DebugStop();
        TPZGeoElSide interface(gelside.Neighbour());
        if(interface.Element()->MaterialId() != interfacematid) DebugStop();
        TPZGeoElSide pressure = gelside.HasNeighbour(pressurematid);
        if(!pressure) DebugStop();
        TPZCompElSide cpressure(pressure.Reference());
        if(!cpressure) DebugStop();
        new TPZMultiphysicsInterfaceElement(*this->mMixedOperator,interface.Element(),celside,cpressure);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh) {
    
    cout << "\n---------------------- Starting MergeMeshes for MHM data structure ----------------------\n" << endl;
    TPZSimpleTimer timer_mm;
    
    if(fInitMatIdForMergeMeshes == -1000) {
        cout << "ERROR! Please, set TMRSApproxSpaceGenerator::::fInitMatIdForMergeMeshes" << endl;
        DebugStop();
    }
    
    int fine_skeleton_matid = mSimData.mTGeometry.m_skeletonMatId;
    int coarse_skeleton_matid = 18;
    std::map<int,int64_t> MatFinetoCoarseElIndex;
//  std::map<int64_t,int64_t> NodeCoarseToNodeFine;
    std::vector<int64_t> NodeCoarseToNodeFine(coarsemesh->NNodes(),-1);
    std::map<int64_t,int64_t> ElCoarseToElFine;
    
//    int temp_bc_mat = -10;
//    // create boundary elements for elements without neighbour of dimension dim
//    { // OBS: why do we need this!!?? We dont since gmsh contains the boundary elements
//        std::map<int, int> num_created;
//        int64_t nel_fine = finemesh->NElements();
//        for (int64_t el = 0; el<nel_fine; el++) {
//            TPZGeoEl *gel = finemesh->Element(el);
//            int dim = gel->Dimension();
//            int nsides = gel->NSides();
//            int firstside = nsides-gel->NSides(dim-1)-1;
//            for (int side = firstside; side<nsides-1; side++) {
//                TPZGeoElSide gelside(gel,side);
//                TPZGeoElSide neighbour = gelside.Neighbour();
//                while(neighbour.Element()->Dimension() != dim){
//                    neighbour = neighbour.Neighbour();
//                }
//                if(neighbour == gelside){
//                    TPZGeoElBC gelbc(gelside,temp_bc_mat);
//                    num_created[gel->MaterialId()]++;
//                }
//            }
//        }
//#ifdef PZDEBUG
//        for (auto it : num_created) {
//            std::cout << "For matid " << it.first << " number of elements created " << it.second << std::endl;
//        }
//#endif
//    }
    // Find the correspondence between coarse nodes and fine nodes (MOST EXPENSIVE OPERATION)
    // For each coarse node, find its correspondence in the fine mesh and store the information in NodeCoarseToNodeFine
    // This is used later to create the skeleton mesh
    { // Why do we need this?
        int64_t nnode_coarse = coarsemesh->NNodes();
        for (int64_t n = 0; n<nnode_coarse; n++) {
            TPZGeoNode &no = coarsemesh->NodeVec()[n];
            if(no.Id() == -1) continue;
            TPZManVector<REAL,3> co(3);
            no.GetCoordinates(co);
            int64_t fineindex;
            TPZGeoNode *finenode = finemesh->FindNode(co,fineindex);
            if(!finenode) DebugStop();
            NodeCoarseToNodeFine[n] = fineindex;
        }
    }
    // identify the correspondence between the material id of the fine mesh and the coarse element index
    // of the coarse mesh
    // this also defines the subdomain of the elements
    {
    int64_t first3DCoarse = 0;
    int dim = coarsemesh->Dimension();
    {
        int64_t nelcoarse = coarsemesh->NElements();
        for (int64_t el=0; el<nelcoarse; el++) {
            TPZGeoEl *gel = coarsemesh->Element(el);
            if(gel->Dimension() == dim){
                first3DCoarse = el;
                break;
            }
        }
    }
    
    // Create Map between matid of fine elements and coarse element index
    int64_t nel_fine = finemesh->NElements();
    mSubdomainIndexGel.Resize(nel_fine);
    mSubdomainIndexGel.Fill(-1);
    for (int64_t el = 0; el<nel_fine; el++) {
        auto *gel = finemesh->Element(el);
        if(!gel || gel->HasSubElement()) continue;
        if(gel->Dimension() != dim) continue;
        int matid = gel->MaterialId();
        mSubdomainIndexGel[el] = matid-fInitMatIdForMergeMeshes+first3DCoarse;
#ifdef PZDEBUG
        if(MatFinetoCoarseElIndex.find(matid) == MatFinetoCoarseElIndex.end()){
            TPZManVector<REAL,3> xcenter(3);
            TPZGeoElSide gelside(gel);
            gelside.CenterX(xcenter);
            TPZManVector<REAL,3> qsi(dim,0.);
            int64_t coarse_index = 0;
            
            TPZGeoEl *gelcoarse = coarsemesh->FindElementCaju(xcenter, qsi, coarse_index, dim);
            if(gelcoarse->IsInParametricDomain(qsi,0)) {
                if(coarse_index-first3DCoarse != matid - fInitMatIdForMergeMeshes) DebugStop();
                MatFinetoCoarseElIndex[matid] = coarse_index;
            } else {
                std::cout << "FindElementCaju failed qsi  = " << qsi << std::endl;
                MatFinetoCoarseElIndex[matid] = matid-fInitMatIdForMergeMeshes+first3DCoarse;
            }
        }
#else
        MatFinetoCoarseElIndex[matid] = matid-fInitMatIdForMergeMeshes+first3DCoarse;
#endif
    }
#ifdef PZDEBUG
//    for(auto it : MatFinetoCoarseElIndex){
//        std::cout << "Fine mat id " << it.first << " coarse element index " << it.second << std::endl;
//    }
#endif
    }
    // modify the material id of the boundary elements of the fine mesh (EXPENSIVE OPERATION)
//    { // Not needed since gmsh contains the boundary elements
//    int64_t nel_fine = finemesh->NElements();
//    int meshdim = finemesh->Dimension();
//    std::map<int,int> created_by_mat;
//    for (int64_t el = 0; el < nel_fine; el++) {
//        TPZGeoEl *gel = finemesh->Element(el);
//        if(gel->MaterialId() == temp_bc_mat){
//            int dim = gel->Dimension();
//            if(dim != meshdim-1) continue;
//            TPZGeoElSide gelside(gel);
//            TPZManVector<REAL,3> xcenter(3);
//            gelside.CenterX(xcenter);
//            int64_t elindex3D = 0;
//            TPZManVector<REAL, 3> qsi3D(3,0.);
//            coarsemesh->FindElementCaju(xcenter, qsi3D, elindex3D, meshdim);
//            TPZGeoEl *coarsegel3D = coarsemesh->Element(elindex3D);
//            int coarseside3D = coarsegel3D->WhichSide(qsi3D);
//            if(coarsegel3D->SideDimension(coarseside3D) != dim) DebugStop();
//            TPZGeoElSide BCSide(coarsegel3D,coarseside3D);
//            TPZGeoElSide neighbour = BCSide.Neighbour();
//            while(neighbour != BCSide){
//                if(neighbour.Element()->Dimension() == dim) break;
//                neighbour = neighbour.Neighbour();
//            }
//            if(neighbour == BCSide) DebugStop();
//            int bc_id = neighbour.Element()->MaterialId();
//            created_by_mat[bc_id]++;
//            gel->SetMaterialId(bc_id);
//        }
//    }
//    for (auto it : created_by_mat) {
//        std::cout << "For matid " << it.first << " number of elements created " << it.second << std::endl;
//    }
//    auto subsize = mSubdomainIndexGel.size();
//    }
    auto finesize = finemesh->NElements();
    mSubdomainIndexGel.Resize(finesize, -1);
    // create a Skeleton element between the large elements of the coarse mesh
    std::map<std::pair<int64_t,int64_t>, int64_t> CoarseFaceEl;
    {
    int64_t nel = coarsemesh->NElements();
    int dim = coarsemesh->Dimension();
    for(int64_t el = 0; el<nel; el++){
        TPZGeoEl *gel = coarsemesh->Element(el);
        int geldim = gel->Dimension();
        if(!gel || gel->HasSubElement()) continue;
        if(geldim != 3) continue;
        int firstside = gel->NSides()-gel->NSides(dim-1)-1; // first face side
        for (int side = firstside; side < gel->NSides()-1; side++) { // skip volume
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside){ // look for volume neighbor
                if(neighbour.Element()->Dimension() == dim) break;
                neighbour = neighbour.Neighbour();
            }
            if(neighbour == gelside) continue;
            int64_t neighindex = neighbour.Element()->Index();
            TPZGeoElBC gelbc(gelside,coarse_skeleton_matid);
            std::pair<int64_t, int64_t> leftright(el,neighindex);
            if(neighindex < el) leftright = std::pair<int64_t, int64_t>(neighindex,el);
            CoarseFaceEl[leftright] = gelbc.CreatedElement()->Index();
            // will create two skeleton elements for each interface. Why?
        }
    }
    }
    
    // duplicate the skeleton elements of the coarse mesh within the fine mesh
    {
    int64_t nelcoarse = coarsemesh->NElements();
    int meshdim = coarsemesh->Dimension();
    for (int64_t el = 0; el<nelcoarse; el++) {
        auto gel = coarsemesh->Element(el);
        if(!gel || gel->HasSubElement()) continue;
        int matid = gel->MaterialId();
        if(matid != coarse_skeleton_matid) continue; // only skeleton elements
        int nnode = gel->NNodes();
        TPZManVector<int64_t, 8> nodeindices(nnode);
        for(int n=0; n<nnode; n++){
            int64_t node_index_coarse = gel->NodeIndex(n);
            int64_t node_index_fine = NodeCoarseToNodeFine[node_index_coarse];
            nodeindices[n] = node_index_fine;
        }
        auto eltype = gel->Type();
        int64_t fine_index;
        finemesh->CreateGeoElement(eltype, nodeindices, matid, fine_index);
        ElCoarseToElFine[el] = fine_index;
    }
    }
    // the pair represents the subdomain indices of the elements
    // the integer is the element index of the skeleton element in the fine mesh
    std::map<std::pair<int64_t,int64_t>, int64_t> FineFaceEl;
    
    for(auto it : CoarseFaceEl){ // loop over skeleton coarse face elements of the coarse mesh
        auto coarsepair = it.first;
        int64_t coarseface = it.second;
        //        if(ElCoarseToElFine.find(coarsepair.first) == ElCoarseToElFine.end()) DebugStop();
        //        if(ElCoarseToElFine.find(coarsepair.second) == ElCoarseToElFine.end()) DebugStop();
        if(ElCoarseToElFine.find(coarseface) == ElCoarseToElFine.end()) DebugStop();
        FineFaceEl[coarsepair] = ElCoarseToElFine[coarseface];
    }
    finemesh->BuildConnectivity();
    // modify the material id of the volumetric elements of the fine mesh
    // Why not do this earlier?
    {
        int64_t nel_fine = finemesh->NElements();
        int dim = finemesh->Dimension();
        for (int64_t el = 0; el < nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() != dim) continue;
            int matid = gel->MaterialId();
            if(MatFinetoCoarseElIndex.find(matid) == MatFinetoCoarseElIndex.end()){
                continue;
            }
            int64_t coarse_index = MatFinetoCoarseElIndex[matid];
//            int64_t fine_index = ElCoarseToElFine[coarse_index]; // Do we need this?
            TPZGeoEl *father = coarsemesh->Element(coarse_index);
            int fathermatid = father->MaterialId();
            gel->SetMaterialId(fathermatid);
        }
    }
    
    // create face elements along the small elements as sons of macroscopic faces
    {
    // identify lists of element/sides that connect two subdomains (in the fine mesh)
    std::map<std::pair<int64_t,int64_t>, std::list<int64_t>> facelist;
    {
        int64_t nel = finemesh->NElements();
        int dim = finemesh->Dimension();
        for(int64_t el = 0; el<nel; el++){
            TPZGeoEl *gel = finemesh->Element(el);
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() != dim) continue;
            int64_t domain = mSubdomainIndexGel[el];
            if(domain == -1) DebugStop();
            int firstside = gel->FirstSide(dim-1); // face sides
            for (int side = firstside; side < gel->NSides()-1; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside){
                    if(neighbour.Element()->Dimension() == dim) break;
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == gelside) continue;
                int64_t neighdomain = mSubdomainIndexGel[neighbour.Element()->Index()];
                if(neighdomain == -1) DebugStop();
                if(neighdomain < domain){
                    TPZGeoElBC gbc(neighbour,fine_skeleton_matid);
                    std::pair<int64_t,int64_t> leftright(neighdomain,domain);
                    facelist[leftright].push_back(gbc.CreatedElement()->Index());
                }
            }
        }
    }
    // create the refinement patterns between small element/side and skeleton elements
    // facelist : key : left/right domain
    // second : list of geometric element indexes of (dim-1) face elements
    for (auto it : facelist) {
        if(FineFaceEl.find(it.first) == FineFaceEl.end()) DebugStop();
        int64_t fine_skel = FineFaceEl[it.first]; // change name
        int nelmesh = it.second.size()+1;
        TPZManVector<TPZGeoEl *> gelvec(nelmesh);
        gelvec[0] = finemesh->Element(fine_skel);
        int64_t count = 1;
        for(auto itel : it.second) gelvec[count++] = finemesh->Element(itel);
#ifdef PZDEBUG
        REAL Area = gelvec[0]->Volume();
        REAL Sum = 0.;
        for(int i=1; i<gelvec.size(); i++) Sum += gelvec[i]->Volume();
        REAL diff = Area-Sum;
//        std::cout << "Skeleton area of el " << fine_skel << " area " << Area << " sum of small " << Sum << std::endl;
#endif
        TPZGeoMesh gmeshrefpattern;
        TPZRefPatternTools::GenerateGMeshFromElementVec(gelvec,gmeshrefpattern);
        TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmeshrefpattern);
        TPZGeoEl *gelcoarse = finemesh->Element(fine_skel);
        gelcoarse->SetRefPattern(refpat);
        for(int i=1; i<gelvec.size(); i++){
            gelvec[i]->SetFather(gelvec[0]);
            gelvec[0]->SetSubElement(i-1, gelvec[i]);
        }
    }
    auto subsize = mSubdomainIndexGel.size();
    auto finesize = finemesh->NElements();
    mSubdomainIndexGel.Resize(finesize, -1);
    }
    // complement the domain of the lower dimensional elements. If all volumetric neighbours share
    // the same subdomain, the element belongs to "that" domain
    {
		int64_t nel = finemesh->NElements();
        int dim = finemesh->Dimension();
		for (int64_t el = 0; el<nel; el++) {
			TPZGeoEl *gel = finemesh->Element(el);
			if(!gel || gel->HasSubElement()) continue;
			int geldim = gel->Dimension();
			int64_t domain = mSubdomainIndexGel[el];
			if(geldim == dim && domain == -1) DebugStop();
			if(geldim < dim && domain != -1) continue;
			TPZGeoElSide gelside(gel);
			TPZGeoElSide neighbour = gelside.Neighbour();
			std::set<int64_t> neighdomains;
			while(neighbour != gelside){
				int64_t locdomain = mSubdomainIndexGel[neighbour.Element()->Index()];
				if(locdomain != -1) neighdomains.insert(locdomain);
				neighbour = neighbour.Neighbour();
			}
			if(neighdomains.size() == 1){
				mSubdomainIndexGel[el] = *neighdomains.begin();
			}
		}
    }
    // set the boundary of the fractures to no flow **** WATCH OUT FOR THIS **** TO BE ADJUSTED
//    {
//        int64_t nel = finemesh->NElements();
//        int dim = finemesh->Dimension();
//        for (int64_t el = 0; el<nel; el++) {
//            TPZGeoEl *gel = finemesh->Element(el);
//            int matid = gel->MaterialId();
//            if(matid == temp_bc_mat && gel->Dimension() != dim-2){
//                std::cout << "gel index " << gel->Index() << " dim " << gel->Dimension() << " matid " << matid << std::endl;
//            }
//            if(gel->Dimension() != dim-2) continue;
//            if(matid == temp_bc_mat) matid = -11;
//            gel->SetMaterialId(matid);
//        }
//    }
	
	// TODO: Find a better way to set the mSubDomainIndex of the fracture and intersection elements
	// For now, just to make it work, we are assigning the subdomaindex of the first
	// higher dimensional element it finds.
    // @TODO verify if the fracture elements havent been assigned a subdomain in another part of the code
	for (auto gel: finemesh->ElementVec()) {
		if(!gel || gel->HasSubElement()) continue;
//		if(gel->MaterialId() != FractureMatId()) continue;
        if(!IsFracMatId(gel->MaterialId())) continue;
		TPZGeoElSide gelside(gel);
        int geldomain = mSubdomainIndexGel[gel->Index()];
		TPZGeoElSide neigh = gelside.Neighbour();
		while (neigh != gelside && geldomain == -1) {
			TPZGeoEl* neighel = neigh.Element();
			if (neighel->Dimension() == 3) {
				const int64_t neighelindex = neighel->Index();
				const int64_t domainSubIndex = mSubdomainIndexGel[neighelindex];
				if (mSubdomainIndexGel[gel->Index()] != -1)
					break; // this guy is already good
				mSubdomainIndexGel[gel->Index()] = domainSubIndex;
                geldomain = domainSubIndex;
				break;
			}
			neigh++;
		}
		if (neigh == gelside)
			DebugStop();
	}
	
	for (auto gel: finemesh->ElementVec()) {
		if(!gel || gel->HasSubElement()) continue;
		if(gel->MaterialId() != mSimData.mTFracIntersectProperties.m_IntersectionId) continue;
        auto geldomain = mSubdomainIndexGel[gel->Index()];
		TPZGeoElSide gelside(gel);
		TPZGeoElSide neigh = gelside.Neighbour();
		while (neigh != gelside && geldomain == -1) {
			TPZGeoEl* neighel = neigh.Element();
//			if (neighel->Dimension() == 2 && neighel->MaterialId() == FractureMatId()) {
            if (neighel->Dimension() == 2 &&  IsFracMatId(neighel->MaterialId())) {
				const int64_t neighelindex = neighel->Index();
				const int64_t domainSubIndex = mSubdomainIndexGel[neighelindex];
				if (mSubdomainIndexGel[gel->Index()] != -1)
					break; // this guy is already good
				mSubdomainIndexGel[gel->Index()] = domainSubIndex;
                geldomain = domainSubIndex;
				break;
			}
			neigh++;
		}
		if (neigh == gelside)
			DebugStop();
	}
	
	
#ifdef PZDEBUG
    if(0)
    {
        int64_t nel = finemesh->NElements();
        std::map<int,int> numels;
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            int matid = gel->MaterialId();
            numels[matid]++;
        }
        std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
        for(auto it: numels){
            std::cout << "For matid " << it.first << " number of elements " << it.second << std::endl;
        }
    }
	{
//		ofstream out("meshsubdomains.vtk");
		if(mSubdomainIndexGel.size() != finemesh->NElements()) DebugStop();
//		TPZVTKGeoMesh::PrintGMeshVTK(finemesh, out, mSubdomainIndexGel);
	}
#endif
	
#ifdef PZDEBUG
	CheckMeshIntegrity(finemesh);
#endif
    
    std::cout << "Total time: " << timer_mm.ReturnTimeDouble()/1000 << " seconds" << std::endl;
    cout << "\n---------------------- Finished MergeMeshes for MHM data structure ----------------------\n" << endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::CheckMeshIntegrity(TPZGeoMesh* gmesh) {
	for(auto gel : gmesh->ElementVec()) {
		const int geldim = gel->Dimension();
		if(geldim != 3) continue;
		
		const int64_t mdomainindex = mSubdomainIndexGel[gel->Index()];
		const int firstsideFace = gel->FirstSide(2);
		for (int iside = firstsideFace; iside < gel->NSides() - 1; iside++) {
			std::set<int64_t> neighdomains = {mdomainindex};
			TPZGeoElSide gelside(gel,iside);
			TPZGeoElSide neigh = gelside.Neighbour();
			for( ; neigh != gelside ; neigh++){
				TPZGeoEl* neighel = neigh.Element();
				if(neighel->Dimension() != 3) continue;
				const int64_t neighmdomainindex = mSubdomainIndexGel[neighel->Index()];
				neighdomains.insert(neighmdomainindex);
			}
			if (neighdomains.size() > 1) {
				auto skel = gelside.HasNeighbour(mSimData.mTGeometry.m_skeletonMatId);
				if(!skel) DebugStop();
			}
		}
	}
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::VerifySubdomainIntegrity()
{
    TPZGeoMesh *gmesh = mGeometry;
    int64_t nel = gmesh->NElements();
    int interface_matid = mSimData.mTGeometry.mInterface_material_id;
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension()!=2) continue;
        int matid = gel->MaterialId();
        int domain = mSubdomainIndexGel[el];
        // the skeleton have no domain
        if((matid == mSimData.mTGeometry.m_skeletonMatId || matid == mSimData.mTGeometry.m_skeletonMatId-1) && domain != -1)
        {
            DebugStop();
        }
        if((matid == mSimData.mTGeometry.m_skeletonMatId || matid == mSimData.mTGeometry.m_skeletonMatId-1) )
        {
            continue;
        }
        TPZGeoElSide gelside(gel);
        std::set<int> neighdomain;
        auto neig = gelside;
        for(auto neig = gelside.Neighbour();neig != gelside; neig++)
        {
            if(neig.Element()->Dimension() == 3) neighdomain.insert(mSubdomainIndexGel[neig.Element()->Index()]);
        }
        if(neighdomain.find(domain) == neighdomain.end()) DebugStop();
        // skipping the boundary elements
        if(!IsFracMatId(matid)) continue;
        for (int side = gel->NCornerNodes(); side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            auto neigh = gelside.Neighbour();
            int neighdomain = mSubdomainIndexGel[neigh.Element()->Index()];
            int neighdim = neigh.Element()->Dimension();
            int neighmatid = neigh.Element()->MaterialId();
            if(neighdim == 1 && neighdomain != domain)
            {
                DebugStop();
            }
            if((neighdim == 1) && (neighdomain != domain)) DebugStop();
            // we have the side of intersecting fractures
            // the next neighbour needs to be an interface element
//            std::cout << "matid " << matid << " neighmatid " << neighmatid << " neighdim " << neighdim << std::endl;
            if(neighdim == 1 && neighmatid == matid+2)
            {
                auto neighinterface = neigh.Neighbour();
                int neighinterfacematid = neighinterface.Element()->MaterialId();
                int neighinterfacedim = neighinterface.Element()->Dimension();
                int64_t neighinterfaceindex = neighinterface.Element()->Index();
                int neighinterfacedomain = mSubdomainIndexGel[neighinterfaceindex];
                if(neighinterfacematid != interface_matid) DebugStop();
                if(neighinterfacedim != 1) DebugStop();
                if(neighinterfacedomain != domain) DebugStop();
            }
            std::set<int> neighdomains = {domain};
            neigh = gelside;
            for(neigh = gelside.Neighbour(); neigh!= gelside; neigh++)
            {
                int neighdim = neigh.Element()->Dimension();
                int neighmatid = neigh.Element()->MaterialId();
                int neighdomain = mSubdomainIndexGel[neigh.Element()->Index()];
                if(!IsFracMatId(neighmatid)) continue;
                if(neighdomain == -1) DebugStop();
                neighdomains.insert(neighdomain);
            }
            if(neighdomains.size() > 1) neighdomains.insert(-1);
            // lower dimensional elements have to be in the domain set
            neigh = gelside;
            for(neigh = gelside.Neighbour(); neigh!= gelside; neigh++)
            {
                int neighdim = neigh.Element()->Dimension();
                if(neighdim != 1) continue;
                int neighmatid = neigh.Element()->MaterialId();
                int neighdomain = mSubdomainIndexGel[neigh.Element()->Index()];
                if(neighdomains.find(neighdomain) == neighdomains.end()) DebugStop();
            }
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::IdentifySubdomainForLowdimensionElements(TPZCompMesh *fluxmesh)
{
    if(mGeometry->Reference() != fluxmesh)
    {
        mGeometry->ResetReference();
        fluxmesh->LoadReferences();
    }
    const int skeletonmatid = mSimData.mTGeometry.m_skeletonMatId;
    //TODO: LOOK AT THIS BEAUTIFUL PROGRAMMING!! NATHAN JOSE
    const int fineskeletonmatid = 18;
    mSubdomainIndexGel.Resize(mGeometry->NElements(), -1);
    int dim = mGeometry->Dimension();
    for(int lowdim = dim-1; lowdim >= 0; lowdim--)
    {
        for(auto gel : mGeometry->ElementVec())
        {
            if(!gel || gel->Dimension() != lowdim) continue;
            // need to modify a domain of a skeleton element
            if(gel->MaterialId() == skeletonmatid) continue;
            if(gel->MaterialId() == fineskeletonmatid) continue;
            auto index = gel->Index();
            // if the element has a subdomain index already, do nothing
            if(mSubdomainIndexGel[index] != -1) continue;
            TPZCompEl *cel = gel->Reference();
            if(!cel)
            {
                // we have a pressure lagrange element or an interface element
                // identify the domain indices of all geometric elements coupled to it
                // if there is only one domain index, then this is the domain index of the element
                // else the domain index is -1
                std::set<int> neighdomains;
                TPZGeoElSide gelside(gel);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    TPZGeoEl *neighgel = neighbour.Element();
                    if(neighgel->Dimension() == lowdim+1)
                    {
                        TPZCompEl *neighcel = neighgel->Reference();
                        if(neighcel && neighcel->NConnects() != 1) {
                            neighdomains.insert(mSubdomainIndexGel[neighgel->Index()]);
                        }
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(neighdomains.size() == 1)
                {
                    mSubdomainIndexGel[index] = *neighdomains.begin();
                }
            }
            else // cel != 0
            {
                // we have either an HDivBound, or a fracture element
                int nc = cel->NConnects();
                if(nc == 1) // this is an hdivbound element
                {
                    int64_t connectindex = cel->ConnectIndex(0);
                    // working an hdivbound element that is not a skeleton element
                    // look for a neighbour of higher dimension that shares this connect
                    // if found the domain index is equal to the domain index of the higher dimension element
                    TPZGeoElSide gelside(gel);
                    TPZGeoElSide neighbour(gelside.Neighbour());
                    while(neighbour != gelside)
                    {
                        TPZCompEl *neighcel = neighbour.Element()->Reference();
                        if(neighbour.Element()->Dimension() == lowdim+1 && neighcel)
                        {
                            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(neighcel);
                            if(!intel) DebugStop();
                            int neighside = neighbour.Side();
                            int nc = intel->NConnects();
                            if(nc != 1) // exclude hdiv bound neighbours
                            {
                                int64_t sideconnectindex = intel->SideConnectIndex(0, neighside);
                                if(sideconnectindex == connectindex) // we have a computational element that shares a connect
                                {
                                    int subdomain = mSubdomainIndexGel[neighbour.Element()->Index()];
                                    mSubdomainIndexGel[index] = subdomain;
                                }
                            }
                        }
                        neighbour = neighbour.Neighbour();
                    }
                } else // nc != 1 -> HDiv Collapsed element
                {
                    // we have a lower dimensional flux element, we assume of dimension dim-1
                    if(lowdim != dim-1) DebugStop();
                    // look for the flux elements connected to both connects and accumulate the domain ids
                    std::set<int> domidlower, domidupper;
                    int64_t connectlower = cel->ConnectIndex(nc-2);
                    int64_t connectupper = cel->ConnectIndex(nc-1);
                    TPZGeoElSide gelside(gel);
                    TPZGeoElSide neighbour(gelside.Neighbour());
                    while(neighbour != gelside) // loop over the neighbours
                    {
						TPZCompEl* neighcel = neighbour.Element()->Reference();
                        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(neighcel);
                        if(intel) {
							int neighside = neighbour.Side();
							int64_t neighindex = neighbour.Element()->Index();
							int64_t sideconnectindex = intel->SideConnectIndex(0, neighside);
							if(sideconnectindex == connectlower)
							{
								domidlower.insert(mSubdomainIndexGel[neighindex]);
							}
							if(sideconnectindex == connectupper)
							{
								domidupper.insert(mSubdomainIndexGel[neighindex]);
							}
						}
                        neighbour = neighbour.Neighbour();
                    }
					const int noDomain = -1;
                    // if the element is neighbour of a skeleton element then there will be two domains
                    if(domidlower.size() > 1 && domidlower.find(noDomain) == domidlower.end()) DebugStop();
                    if(domidupper.size() > 1 && domidupper.find(noDomain) == domidupper.end()) DebugStop();
                    int domidL = -1;
                    if(domidlower.size() == 1) domidL = *domidlower.begin();
                    int domidU = -1;
                    if(domidupper.size() == 1) domidU = *domidupper.begin();
                    if(domidL != -1 && domidU != -1 && domidL != domidU)
                    {
                        // there must be a skeleton element sharing a connect
                        for(neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++)
                        {
                            int matid = neighbour.Element()->MaterialId();
                            if(matid == skeletonmatid)
                            {
                                int skeldomain = mSubdomainIndexGel[neighbour.Element()->Index()];
                                if(skeldomain != -1) DebugStop();
                                TPZCompEl *cel = neighbour.Element()->Reference();
                                int64_t skelconindex = cel->ConnectIndex(0);
                                if(skelconindex == connectupper)
                                {
                                    domidU = -1;
                                }
                                else if(skelconindex == connectlower)
                                {
                                    domidL = -1;
                                }
                                else
                                {
                                    DebugStop();
                                }
                            }
                        }
                    }
                    if(domidL != -1)
                    {
                        mSubdomainIndexGel[index] = domidL;
                    }
                    if(domidU != -1)
                    {
                        mSubdomainIndexGel[index] = domidU;
                    }
                }
            }
        }
    }
	
    
#ifdef PZDEBUG
    // Function to check all the domains of all elements
//	const int64_t nsub = mSubdomainIndexGel.size();
//	for(int i = 0 ; i < nsub ; i++) {
//		const int macroindex = mSubdomainIndexGel[i];
//		TPZGeoEl* gel = mGeometry->Element(i);
//		const int eldim = gel->Dimension();
//		if(eldim == 3) continue;
//		const int matid = gel->MaterialId();
//		cout << "Element index " << i << " | dim " << eldim << " | matid " << matid << " | domain " << macroindex << endl;
//	}
#endif
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::SetInterfaceDomains(TPZStack<int64_t> &pressureindices,std::pair<int,int> &interfacematids)
{
    TPZCompMesh *cmesh = mMixedOperator;
    TPZGeoMesh *gmesh = cmesh->Reference();
    mSubdomainIndexGel.Resize(gmesh->NElements(), -1);
    for(auto index : pressureindices)
    {
        int dompressure = mSubdomainIndexGel[index];
        TPZGeoEl *gel = gmesh->Element(index);
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour(gelside.Neighbour());
        while(neighbour != gelside)
        {
            int neighmatid = neighbour.Element()->MaterialId();
            if(neighmatid == interfacematids.first || neighmatid == interfacematids.second)
            {
				TPZCompEl* celneigh = neighbour.Element()->Reference();
				const int64_t neighindex = neighbour.Element()->Index();
				TPZMultiphysicsInterfaceElement* mpinterface = dynamic_cast<TPZMultiphysicsInterfaceElement*>(celneigh);
				if(!mpinterface) DebugStop();
				TPZCompEl* celright = mpinterface->RightElement();
				TPZGeoEl* gelright = celright->Reference();
				
                int domneigh = mSubdomainIndexGel[gelright->Index()];
                if(domneigh == dompressure)
                {
                    mSubdomainIndexGel[neighindex] = domneigh;
                }
                else if(domneigh == -1)
                {
                    mSubdomainIndexGel[neighindex] = dompressure;
                }
                else if(dompressure == -1)
                {
                    mSubdomainIndexGel[neighindex] = domneigh;
                }
                else
                {
                    DebugStop();
                }
				cout << "Element index " << neighindex << " | dim " << gel->Dimension() << " | matid " << gel->MaterialId() << " | domain " << mSubdomainIndexGel[neighindex] << endl;
            }
            neighbour = neighbour.Neighbour();
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

bool TMRSApproxSpaceGenerator::IsFracMatId(int matiD){
    return (mSimData.mTFracProperties.m_fracprops.find(matiD) != mSimData.mTFracProperties.m_fracprops.end());
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::OrderFractures(TPZCompMesh *fluxmesh, TPZVec<TPZGeoElSide> &fracvec)
{
    // compute the normal direction of the 3D element with positive side direction
    TPZManVector<REAL,3> normal(3,0.);
    int fracdomain = mSubdomainIndexGel[fracvec[0].Element()->Index()];
    TPZGeoElSide first3D, last3D;
    {
        TPZGeoElSide gelside = fracvec[0];
        for (auto neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            if(neigh.Element()->Dimension() == 3)
            {
                int direction = neigh.Element()->NormalOrientation(neigh.Side());
                if(direction == 1)
                {
                    TPZManVector<REAL,2> centerksi(2,0.);
                    neigh.CenterPoint(centerksi);
                    neigh.Normal(centerksi, normal);
                    first3D = neigh;
                }
                else if(direction == -1)
                {
                    last3D = neigh;
                }
                else DebugStop();
            }
        }
        if(!first3D || !last3D) DebugStop();
    }
    
    // order the TPZGeoElSide as a function of their position in the normal direction
    std::multimap<REAL, TPZGeoElSide> ordered;
    {
        TPZManVector<REAL,3> xcenter(3);
        TPZGeoElSide gelside = fracvec[0];
        gelside.CenterX(xcenter);
        for (int i=0; i<fracvec.size(); i++) {
            int matid = fracvec[i].Element()->MaterialId();
            auto &dfn = mSimData.mTFracProperties.m_fracprops[matid].m_polydata;
            TPZManVector<REAL,3> projected(3);
            projected = dfn.GetProjectedX(xcenter);
            
            REAL dist = 0.;
            for(int c=0; c<3; c++) dist += (projected[c]-xcenter[c])*normal[c];
            ordered.insert({dist,fracvec[i]});
        }
    }
    last3D.RemoveConnectivity();
    for (int i = 0; i<fracvec.size(); i++) {
        fracvec[i].RemoveConnectivity();
    }
  
    
    int gluematid = mSimData.mTFracIntersectProperties.m_FractureGlueId;
    if(gluematid < 0) DebugStop();
    TPZGeoElSide prev = first3D;

    for (auto it = ordered.begin(); it != ordered.end(); it++) {
        if(prev != first3D)
        {
            std::cout<<"FindGlue:: "<<std::endl;
            TPZGeoElBC gbc(prev,gluematid);
            TPZGeoEl *glue = gbc.CreatedElement();
            auto glueindex = glue->Index();
            if(glueindex >= mSubdomainIndexGel.size())
            {
                mSubdomainIndexGel.Resize(glueindex+1, fracdomain);
            }
            mSubdomainIndexGel[glueindex] = fracdomain;
            {
                auto cel = fluxmesh->ApproxSpace().CreateCompEl(glue, *fluxmesh);
                auto cindex = cel->ConnectIndex(0);
                std::cout << "Glue connect index " << cindex << std::endl;
            }
            TPZGeoElSide gels(glue);
            gels.SetConnectivity(it->second);
            prev = it->second;
        }
        else {
            prev.SetConnectivity(it->second);
            prev = it->second;
        }
    }
    prev.SetConnectivity(last3D);
#ifdef PZDEBUG
//    {
//        int64_t index = first3D.Element()->Index();
//        std::cout << "el index " << first3D.Element()->Index() << " dim " << first3D.Element()->Dimension() <<
//        " matid " << first3D.Element()->MaterialId() <<
//        " domain " << this->mSubdomainIndexGel[index] << std::endl;
//        for(auto neigh = first3D.Neighbour(); neigh != first3D; neigh++)
//        {
//            int64_t index = neigh.Element()->Index();
//            std::cout << "el index " << neigh.Element()->Index() << " dim " << neigh.Element()->Dimension() <<
//            " matid " << neigh.Element()->MaterialId() <<
//            " domain " << this->mSubdomainIndexGel[index] << std::endl;
//        }
//    }
#endif
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::OrderOverlappingFractures()
{
    TPZCompMesh *cmesh = mGeometry->Reference();
    int64_t nel = mGeometry->NElements();
    int fracgluematid = mSimData.mTFracIntersectProperties.m_FractureGlueId;
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = mGeometry->Element(el);
        if(!gel) continue;
        int matid = gel->MaterialId();
        if(!IsFracMatId(matid)) continue;
        TPZGeoElSide gelside(gel);
        // if fracglue elements were inserted, the fracture elements are ordered already
        if(gelside.HasNeighbour(fracgluematid)) continue;
        
        
        TPZStack<TPZGeoElSide> allfracs;
        allfracs.Push(gelside);
        for (auto neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            if(IsFracMatId(neigh.Element()->MaterialId())) allfracs.Push(neigh);
        }
        if(allfracs.size() > 1)
        {
            // this will order the neighbour sequence starting with the 3D element, the fractures till the next 3D element
            OrderFractures(cmesh, allfracs);
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TMRSApproxSpaceGenerator::InitializeMemoryFractureGlue()
{
    int64_t nel = mMixedOperator->NElements();
    int matglueid = mSimData.mTFracIntersectProperties.m_FractureGlueId;
    TMRSDarcyFractureGlueFlowWithMem *matglue = dynamic_cast<TMRSDarcyFractureGlueFlowWithMem *>(mMixedOperator->FindMaterial(matglueid));
    if(!matglue) DebugStop();
    std::shared_ptr<TPZAdmChunkVector<TGlueMem>> memvec = matglue->GetMemory();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = mMixedOperator->Element(el);
        TPZMaterial *mat = cel->Material();
        
        if(mat == matglue)
        {
            TPZGeoEl *gel = cel->Reference();
            TPZGeoElSide gelside(gel);
            TPZGeoElSide prev(gelside);
            prev--;
            TPZGeoElSide next(gelside);
            next++;
            int prevmatid = prev.Element()->MaterialId();
            int nextmatid = next.Element()->MaterialId();
            auto &fracprev = mSimData.mTFracProperties.m_fracprops[prevmatid].m_polydata;
            auto &fracnext = mSimData.mTFracProperties.m_fracprops[nextmatid].m_polydata;
            
         
            TPZManVector<int64_t> memindices;
            cel->GetMemoryIndices(memindices);
            const TPZIntPoints &intpoints = cel->GetIntegrationRule();
            int npoints = intpoints.NPoints();
            if(memindices.size() != npoints) DebugStop();
            for (int ip = 0; ip < npoints; ip++) {
                REAL weight;
                TPZManVector<REAL,3> pt(gel->Dimension());
                intpoints.Point(ip, pt, weight);
                TPZManVector<REAL,3> xco(3), xprev(3), xnext(3);
                gel->X(pt, xco);
                xprev = fracprev.GetProjectedX(xco);
                xnext = fracnext.GetProjectedX(xco);
                REAL distance = dist(xprev,xnext);
                int64_t globindex = memindices[ip];
                (*memvec)[globindex].m_xco = xco;
                (*memvec)[globindex].m_dist = distance;
                (*memvec)[globindex].m_fracs = {prevmatid,nextmatid};
            }
            
        }
    }
}

void TMRSApproxSpaceGenerator::UpdateMemoryFractureGlue(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    int matglueid = mSimData.mTFracIntersectProperties.m_FractureGlueId;
    TMRSDarcyFractureGlueFlowWithMem *matglue = dynamic_cast<TMRSDarcyFractureGlueFlowWithMem *>(mMixedOperator->FindMaterial(matglueid));
    if(!matglue) DebugStop();
    std::shared_ptr<TPZAdmChunkVector<TGlueMem>> memvec = matglue->GetMemory();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel){continue;}
        TPZSubCompMesh *subcmesh=dynamic_cast<TPZSubCompMesh *>(cel);
        if (subcmesh) {
            UpdateMemoryFractureGlue(subcmesh);
        }
        TPZMaterial *mat = cel->Material();
        if(mat == matglue)
        {
            TPZGeoEl *gel = cel->Reference();
            TPZGeoElSide gelside(gel);
            TPZGeoElSide prev(gelside);
            prev--;
            TPZGeoElSide next(gelside);
            next++;
            int prevmatid = prev.Element()->MaterialId();
            int nextmatid = next.Element()->MaterialId();
            TPZManVector<int64_t> memindices;
            cel->GetMemoryIndices(memindices);
            const TPZIntPoints &intpoints = cel->GetIntegrationRule();
            int npoints = intpoints.NPoints();
            if(memindices.size() != npoints) DebugStop();
            for (int ip = 0; ip < npoints; ip++) {
                REAL weight;
                TPZManVector<REAL,3> pt(gel->Dimension());
                intpoints.Point(ip, pt, weight);
                TPZManVector<REAL,3> xco(3);
                gel->X(pt, xco);
                int64_t globindex = memindices[ip];
                std::pair<int, int> fracids = (*memvec)[globindex].m_fracs;
//              (*memvec)[globindex].m_xco = xco;
                REAL dist = (*memvec)[globindex].m_dist;
                TPZVec<STATE> solution1, solution2;
                TPZCompEl *cel2 =  next.Element()->Reference();
                next.Element()->Reference()->Solution(pt, 2, solution1); //pressure
                prev.Element()->Reference()->Solution(pt, 2, solution2); //pressure
                REAL dp =solution1[0] - solution2[0];
                (*memvec)[globindex].m_dp = dp;
                REAL flux = (*memvec)[globindex].m_flux;
                
            }
            
        }
    }
}
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
