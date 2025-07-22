
#include "TPZMHMixedMesh4SpacesControl.h"
#include "pzelementgroup.h"
#include "pzelmat.h"
#include "TPZFastCondensedElement.h"
#include "TMRSApproxSpaceGenerator.h"
#include "TPZReservoirTools.h"
#include "pzcompelwithmem.h"
#include "pzsmanal.h"
//#include "ConfigurateCase.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.MixedMeshSpaceControl"));
#endif

/**
 * @brief Build Computational Mesh
 * @param useSubstructure: wether use substructure
 */
void TPZMHMixedMesh4SpacesControl::BuildComputationalMesh(bool useSubstructure){
    //A check for the polynomial order
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        std::cout<<"Wrong polynomial order set"<<std::endl;
        DebugStop();
    }
    InsertPeriferalMaterialObjects();  // Insert BC objects that do not perform any actual computation
    CreateHDivMHMMesh();
    
    InsertPeriferalPressureMaterialObjects();
    
    if(fNState > 1){
        fRotationMesh = new TPZCompMesh(fGMesh);
        InsertPeriferalRotationMaterialObjects();
    }
#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    CreatePressureMHMMesh();
    CreateAverageFlux();
    CreateAveragePressure();
   
    if(fNState > 1)
    {
        CreateRotationMesh();
    }
    
    CreateHDivPressureMHMMesh();
    
    
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    std::cout<<"condensed: "<<fCMesh->NEquations()<<std::endl;
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
#ifdef PZDEBUG
    CheckMeshConsistency();
#endif
    
   
    if (useSubstructure) {
        HideTheElements();
    }
  
    
    fNumeq = fCMesh->NEquations();
}

/**
 * @brief Creates the HDiv Pressure 4 spaces MHM Mesh
 */
void TPZMHMixedMesh4SpacesControl::CreateHDivPressureMHMMesh()
{
    TPZManVector<TPZCompMesh *,4 > cmeshes(4);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fcmeshFluxAverg.operator->();
    cmeshes[3] = fcmeshPressureAverg.operator->();
    if (1){
        std::ofstream out("PressureMesh_MultiPhis.txt");
        cmeshes[1] ->Print(out);
        
        std::ofstream out2("FluxMesh_MultiPhis.txt");
        cmeshes[0] ->Print(out2);
        
        std::ofstream out3("FluxAverage_MultiPhis.txt");
        cmeshes[2] ->Print(out3);
        
        std::ofstream out4("PressureAverage_MultiPhis.txt");
        cmeshes[3] ->Print(out4);
    };
    
    
    if(fNState > 1) {
        cmeshes.Resize(3);
        cmeshes[2] = fRotationMesh.operator->();
    }
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    // Multiphysics mesh
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    gSinglePointMemory = true;

    BuildMultiPhysicsMesh();
    
    // neste ponto eh preciso inicializar as regras de integracao
    std::cout << "Initializing the material memory objects\n";
    
    {
        int dim = MixedFluxPressureCmesh->Dimension();
        int64_t nelem = MixedFluxPressureCmesh->NElements();
        for(int64_t el=0; el<nelem; el++)
        {
            TPZCompEl *cel = MixedFluxPressureCmesh->Element(el);
            if(!cel) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            if(gel->Dimension() == dim)
            {
                cel->PrepareIntPtIndices();
            }
        }
    }
    
    TPZManVector<TPZCompMesh * ,4> meshvector;
    
    meshvector = cmeshes;
  
    
    // Populate the connects to subdomain data structure for the multiphysics mesh
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
    //    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    //    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    // copy the solution of the atomic meshes to the multiphysics mesh
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
   
    
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
//    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
//    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
    MixedFluxPressureCmesh->ComputeNodElCon();
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    return;
}

/**
 * @brief Creates Distributed Flux mesh
 */
void TPZMHMixedMesh4SpacesControl::CreateAverageFlux()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    TPZCompMesh * cmeshtemp = new TPZCompMesh(gmesh);
    fcmeshFluxAverg = cmeshtemp;
    
    //The distributed flux mesh should be empty when calling this method
    int64_t nskeletonconnects = fcmeshFluxAverg->NConnects();
    if(nskeletonconnects != 0){     //Check that it is empty
        DebugStop();
    }

    // create and organize the distributed flux mesh
    TPZCompMesh * cmeshfluxavg = fcmeshFluxAverg.operator->();
    gmesh->ResetReference();
    cmeshfluxavg->SetName("DistributedFlux");
    cmeshfluxavg->SetDimModel(gmesh->Dimension());
    //
    cmeshfluxavg->SetAllCreateFunctionsDiscontinuous(); //AQUI
    //
    cmeshfluxavg->ApproxSpace().CreateDisconnectedElements(true);
    cmeshfluxavg->SetDefaultOrder(0);
    
//    generate elements for all material ids of mesh dim
    std::set<int> matids;

    TPZNullMaterial<STATE> * volume = new TPZNullMaterial(1);
    cmeshfluxavg->InsertMaterialObject(volume);
    TPZNullMaterial<STATE> * volume2 = new TPZNullMaterial(2);
    cmeshfluxavg->InsertMaterialObject(volume2);
    matids.insert(1);
    matids.insert(2);
    cmeshfluxavg->AutoBuild(matids);
    fcmeshFluxAverg->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("DistributedFluxMesh.txt");
        cmeshfluxavg->Print(out);
    }
    
    
#ifdef PZDEBUG
    // a very strange check!! Why does material id 1 need to be volumetric?
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshfluxavg->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshfluxavg->ConnectVec()[ic].SetLagrangeMultiplier(2); //Distributed flux
    }
    // associate the connects with the proper subdomain
    gmesh->ResetReference();
    int64_t nel = cmeshfluxavg->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshfluxavg->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        // if a computational element was created outside the range of material ids
        // something very strange happened...
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            DebugStop();
        }
        int domain = fGeoToMHMDomain[gel->Index()];
#ifdef PZDEBUG
        if (domain == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, domain);
    }
    
    return;
}

/**
 * @brief Creates Average Pressure mesh
 */
void TPZMHMixedMesh4SpacesControl::CreateAveragePressure()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();

    gmesh->ResetReference();
    TPZCompMesh * cmeshtemp = new TPZCompMesh(gmesh);
    fcmeshPressureAverg = cmeshtemp;
    
    // the pressure mesh should be empty when calling this method
    int64_t nskeletonconnects = fcmeshPressureAverg->NConnects();
    if(nskeletonconnects != 0){
        DebugStop();
    }
    
    // create and organize the pressure mesh
    // the pressure mesh is composed of discontinuous H1 elements
    TPZCompMesh * cmeshpressureavr = fcmeshPressureAverg.operator->();
    gmesh->ResetReference();
    cmeshpressureavr->SetName("PressureMeshAverage");
    cmeshpressureavr->SetDimModel(gmesh->Dimension());
    cmeshpressureavr->SetAllCreateFunctionsDiscontinuous(); //AQUI
    cmeshpressureavr->ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressureavr->SetDefaultOrder(0);
    
    // generate elements for all material ids of meshdim
    std::set<int> matids;
    //    for (auto it:fMaterialIds) {
    //        TPZMaterial *mat = fPressureFineMesh->FindMaterial(it);
    //        if (mat && mat->Dimension() == meshdim) {
    //            matids.insert(it);
    //            cmeshfluxavg->InsertMaterialObject(mat);
    //        }
    //    }
    TPZNullMaterial<STATE> * volume = new TPZNullMaterial(1);
    cmeshpressureavr->InsertMaterialObject(volume);
    TPZNullMaterial<STATE> * volume2 = new TPZNullMaterial(2);
    cmeshpressureavr->InsertMaterialObject(volume2);
    matids.insert(1);
    matids.insert(2);
    cmeshpressureavr->AutoBuild(matids);
    fcmeshPressureAverg->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("PressureAVERAGEMesh.txt");
        cmeshpressureavr->Print(out);
    }
    
    
#ifdef PZDEBUG
    // a very strange check!! Why does material id 1 needs to be volumetric?
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshpressureavr->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshpressureavr->ConnectVec()[ic].SetLagrangeMultiplier(3);
    }
    // associate the connects with the proper subdomain
    gmesh->ResetReference();
    int64_t nel = cmeshpressureavr->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshpressureavr->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        // if a computational element was created outside the range of material ids
        // something very strange happened...
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            DebugStop();
        }
        int domain = fGeoToMHMDomain[gel->Index()];
#ifdef PZDEBUG
        if (domain == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, domain);
    }

    return;
}

/**
 * @brief Builds the MultiPhysics Mesh
 */
void TPZMHMixedMesh4SpacesControl::BuildMultiPhysicsMesh()
{
    //Checks that the mesh is empty before creation
    if (fCMesh->NElements() != 0) {
        std::cout<<"The mesh has to be empty to build it at this stage"<<std::endl;
        DebugStop();
    }
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    fCMesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(fCMesh.operator->());
   
    int vecsize = 4;
    TPZManVector<TPZCompMesh *> meshvec(vecsize);
    meshvec[0] = fFluxMesh.operator->();
    meshvec[1] = fPressureFineMesh.operator->();
    meshvec[2] = fcmeshFluxAverg.operator->();
    meshvec[3] = fcmeshPressureAverg.operator->();
//    if(fNState > 1)
//    {
//        meshvec[2] = this->fRotationMesh.operator->();
//    }
    TPZManVector<int64_t> shouldcreate(fGMesh->NElements(),0);
    std::set<int> matids;
    for (auto it : fCMesh->MaterialVec()) {
        matids.insert(it.first);
    }
    int64_t nel = fFluxMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        // this means that all geometric elements associated with flux elements will generate a computational element
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    nel = fcmeshPressureAverg->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fcmeshPressureAverg->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    // define the intersection of the finest references
    nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (shouldcreate[el])
        {
            TPZGeoEl *fat = gel->Father();
            while(fat)
            {
                if(shouldcreate[fat->Index()] == 1)
                {
                    shouldcreate[fat->Index()] = 0;
                }
                fat = fat->Father();
            }
        }
    }
    TPZStack<int64_t> gelindexes;
    for (int64_t el=0; el<nel; el++) {
        if (shouldcreate[el])
        {
            gelindexes.Push(el);
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric indices for which we will create multiphysics elements" << std::endl;
        sout << gelindexes;
        //        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str());
    }
#endif

    mphysics->BuildMultiphysicsSpace(meshvec,gelindexes);
}

/**
 * @brief Hide The Elements, allocate the connects in the submeshes and condense
 */
void TPZMHMixedMesh4SpacesControl::HideTheElements()
{
   
    int KeepOneLagrangian = 3;
    if (fHybridize) {
        KeepOneLagrangian = false;
    }
    
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fCMesh->Reference();
    gmesh->ResetReference();
    fCMesh->LoadReferences();
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        int64_t domain = WhichSubdomain(cel);
        if (domain == -1) {
            continue;
        }
        ElementGroups[domain].insert(el);
    }
    if (ElementGroups.size() <= 5)
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<int64_t,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << fCMesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, KeepOneLagrangian);
    

    std::cout << "After putting in substructures\n";
    fMHMtoSubCMesh = submeshindices;
    fCMesh->ComputeNodElCon();
    fCMesh->CleanUpUnconnectedNodes();
   
    GroupandCondenseElements();
//    GroupandCondenseElementsEigen();
    
    std::cout << "Finished substructuring\n";
}

/**
 * @brief Group and Condense Elements
 */
void TPZMHMixedMesh4SpacesControl::GroupandCondenseElementsEigen()
{
//    int el = 1;
//    int nels = fCMesh->NElements();
    
    
    
    for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
//    for (int el=0; el<nels; el++) {
        TPZCompEl *cel = fCMesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
//            DebugStop();
            continue;
        }
        TPZElementMatrixT<STATE> ek;
        TPZElementMatrixT<STATE> ef;
       
//        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();

        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // Increment nelconnected of exterior connects
        
        int nel = subcmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = subcmesh->Element(el);
            if (!cel) {
                continue;
            }
            int nconnects = cel->NConnects();
            for (int icon=0; icon<nconnects; icon++) {
                TPZConnect &connect = cel->Connect(icon);

                int lagrangemult = connect.LagrangeMultiplier();
                //Increment the number of connected elements for the avg pressure in order to not condense them
                if (lagrangemult==3) {
                    connect.IncrementElConnected();
                }
            }
        }
      
        TPZReservoirTools::CreatedCondensedElements(subcmesh, false, true);
//        TPZCompMeshTools::CreatedCondensedElements(subcmesh, false);
        subcmesh->CleanUpUnconnectedNodes();
      
       
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
//        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
        
        subcmesh->SaddlePermute();
        subcmesh->PermuteExternalConnects();
        TPZSymetricSpStructMatrixEigen *strloc = new TPZSymetricSpStructMatrixEigen(subcmesh);
        TPZAutoPointer<TPZStructMatrix> str(strloc);
        str->SetNumThreads(0);
        int64_t numinternal = subcmesh->NumInternalEquations();
        str->EquationFilter().SetMinMaxEq(0, numinternal);
        TPZAutoPointer<TPZMatrix<STATE>> mat = dynamic_cast<TPZMatrix<STATE>*>(str->Create());
        str->EquationFilter().Reset();
        TPZAutoPointer<TPZLinearAnalysis> Analysis = new TPZSubMeshAnalysis(subcmesh);
        Analysis->SetStructuralMatrix(str);
        TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(mat);
        step->SetDirect(ELDLt);
        TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
        Analysis->SetSolver(autostep);
//        subcmesh->SetAnalysis(Analysis);
        
//        std::ofstream filehide2("subcmeshAfter.txt");
//        subcmesh->Print(filehide2);
    }
}

void TPZMHMixedMesh4SpacesControl::GroupandCondenseElements()
{
//    int el = 1;
//    int nels = fCMesh->NElements();
    
    
    
    for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
//    for (int el=0; el<nels; el++) {
        TPZCompEl *cel = fCMesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
//            DebugStop();
            continue;
        }
        TPZElementMatrixT<STATE> ek;
        TPZElementMatrixT<STATE> ef;
       
//        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();

        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // Increment nelconnected of exterior connects
        
        int nel = subcmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = subcmesh->Element(el);
            if (!cel) {
                continue;
            }
            int nconnects = cel->NConnects();
            for (int icon=0; icon<nconnects; icon++) {
                TPZConnect &connect = cel->Connect(icon);

                int lagrangemult = connect.LagrangeMultiplier();
                //Increment the number of connected elements for the avg pressure in order to not condense them
                if (lagrangemult==3) {
                    connect.IncrementElConnected();
                }
            }
        }
      
        TPZReservoirTools::CreatedCondensedElements(subcmesh, false, true);
//        TPZCompMeshTools::CreatedCondensedElements(subcmesh, false);
        subcmesh->CleanUpUnconnectedNodes();
      
        int numthreads = 0;
        int preconditioned = 0;
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        subcmesh->SetAnalysisSkyline(numthreads, preconditioned);
        
        
        
//        std::ofstream filehide2("subcmeshAfter.txt");
//        subcmesh->Print(filehide2);
    }
}

/**
 * @brief Puts the element set into a subcompmesh and make the connects internal
 * @param cmesh: Computational mesh
 * @param elindices: pointer with the Element indices
 * @param indices: pointer with the Element indices
 * @param KeepOneLagrangian: wether keep the one lagrangian or not
 */
void TPZMHMixedMesh4SpacesControl::PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, bool KeepOneLagrangian)

{
    for (std::map<int64_t,std::set<int64_t> >::iterator it = elindices.begin(); it != elindices.end(); it++) {
        TPZSubCompMesh *subcmesh = new TPZSubCompMesh(*cmesh);
        indices[it->first] = subcmesh->Index();
        for (std::set<int64_t>::iterator itloc = it->second.begin(); itloc != it->second.end(); itloc++) {
            subcmesh->TransferElement(cmesh, *itloc);
        }
    }
    cmesh->ComputeNodElCon();
    for (std::map<int64_t,int64_t>::iterator it = indices.begin(); it != indices.end(); it++) {
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh->Element(it->second));
        if (!subcmesh) {
            DebugStop();
        }
        int count = 0;
        if (KeepOneLagrangian)
        {
            int64_t nconnects = subcmesh->NConnects();
            for (int64_t ic=0; ic<nconnects; ic++) {
                TPZConnect &c = subcmesh->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    count++;
                    if(count == 1 && c.NState() == 1)
                    {
                        break;
                    }
                    else if(count == 2 && c.NState() == 2)
                    {
                        break;
                    }
                    else if(count == 3 && c.NState() == 3)
                    {
                        break;
                    }
                }
            }
        }
        subcmesh->MakeAllInternal();
    }
}
