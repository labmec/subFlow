#include "TPZReservoirTools.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("imrs.reservoirtools");
#endif

void TPZReservoirTools::CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix)
{
    //    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        if (KeepOneLagrangian) {
            int count = 0;
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    count++;
                    if(count == 1 && c.NState() == 1)
                    {
                        break;
                    } else if(count == 2 && c.NState() == 2)
                    {
                        break;
                    } else if(count == 3 && c.NState() == 3)
                    {
                        break;
                    }
                    
                }
            }
        }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
        {
            TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
            cond->SetLambda(1.0);
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}
/// create a condensed element and do not condense the connect with a given lagrange level
// the method does the same procedure as CreatedCondensedElements, but has different policy for
// keeping a connect out the condensation loop
void TPZReservoirTools::CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix)
{
    //    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        //        std::cout << "Element " << el << std::endl;
        TPZCompEl *cel = cmesh->Element(el);
     
        
        if (!cel) {
            continue;
        }
        
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcmesh)
        {
            CondenseElements(subcmesh, LagrangeLevelNotCondensed, keepmatrix);
            continue;
        }
      
        int nc = cel->NConnects();
        if(LagrangeLevelNotCondensed >=0)
        {
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                
//                char lag =c.LagrangeMultiplier();
//                if (lag==6) {
//                    int ncon =c.NElConnected();
//                    if (ncon==2) {
//                        c.DecrementElConnected();
//                    }
//                }
                
                if((c.LagrangeMultiplier() >= LagrangeLevelNotCondensed && c.NElConnected() == 1) )
                {
                    c.IncrementElConnected();
                }
            }
        }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) { // if the el has dependency or if it is a boundary
                continue;
            }
            break;
        }
        // we can condense if there is a connect without dependency or a connect with more than one element connected
        // maybe we should condense as a function of the material id?
        bool cancondense = (ic != (nc));
        if(cancondense)
        {
            TPZGeoEl *gel = cel->Reference();
            TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
            cond->SetLambda(1.0);
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}

//Condense by MATID

void TPZReservoirTools::CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix, std::set<int> matids)
{
    
    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        
        int nc = cel->NConnects();
        bool found = false;
        if(LagrangeLevelNotCondensed >=0)
          {
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                char lag =c.LagrangeMultiplier();
                if (lag==6) {
                    
                    int ncon =c.NElConnected();
                    if (ncon==2) {
                        c.DecrementElConnected();
                    }
                    
                }
                
                if((c.LagrangeMultiplier() >= LagrangeLevelNotCondensed && c.NElConnected() == 1))
                  {
                    c.IncrementElConnected();
                    found = true;
                  }
            }
          }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
          {
            int gel_matId = 0;
            if (gel) {
                gel_matId = gel->MaterialId();
            }
            
            //            if(LagrangeLevelNotCondensed >= 0 && !found) DebugStop();
            int verif = 0;
            
            for (auto matId:matids) {
                if (gel_matId==matId) {
                    verif=1;
                    break;
                }
            }
            if (verif==1) {
                TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
                cond->SetLambda(1.0);
            }
            else{
                TPZCondensedCompElT<REAL> *cond = new TPZCondensedCompElT<REAL>(cel, keepmatrix);
                
            }
          }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}

void TPZReservoirTools::FindCondensed(TPZCompEl *cel, TPZStack<TPZFastCondensedElement *> &condensedelements)
{
    TPZFastCondensedElement *cond = dynamic_cast<TPZFastCondensedElement *>(cel);
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    if(cond)
    {
        condensedelements.Push(cond);
        FindCondensed(cond->ReferenceCompEl(), condensedelements);
        return;
    }
    if(elgr)
    {
        const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
        int nel = elvec.size();
        for (int el = 0; el<nel; el++) {
            TPZCompEl *loccel = elvec[el];
            FindCondensed(loccel, condensedelements);
        }
    }
}

void TPZReservoirTools::AddDependency( std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons){
    
    for(auto pairfatherson:fatherAndSons){
        TPZCompEl *celfat = pairfatherson.first->Reference();
        TPZInterpolatedElement *intelFather = dynamic_cast<TPZInterpolatedElement*>(celfat);
        if(celfat->NConnects() !=1){
            DebugStop();
        }
        for(auto son:pairfatherson.second){
            TPZCompEl *celSon= son->Reference();
            if(celSon->NConnects() !=1){
                DebugStop();
            }
            TPZInterpolatedElement *intelSon = dynamic_cast<TPZInterpolatedElement*>(celSon);
            intelSon->RestrainSide(son->NSides()-1, intelFather, pairfatherson.first->NSides()-1);
        }
    }
}
void TPZReservoirTools::TakeFatherSonsCorrespondence(TPZCompMesh *fluxCmesh, TPZVec<int64_t> &subdomain, std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons){
    

    TPZGeoMesh *gmesh = fluxCmesh->Reference();
//    gmesh->ResetReference();
    fluxCmesh->LoadReferences();
    std::map<int, std::vector<TPZGeoEl* >> interfaces;
    std::vector<int> matIds;
    int skeletonCoarseId = 19;
    matIds.push_back(skeletonCoarseId);
    // the interfaces data structure will contain all geoelement pointers with skeletonCoarseId
    TakeElementsbyID(gmesh, interfaces, matIds);
    int count=-1;
#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        LOGPZ_DEBUG(logger, "number of skeleton elements " << interfaces.size())
    }
#endif
    for(auto gel:interfaces[skeletonCoarseId] ){
        TPZStack<TPZGeoEl *> Sons;
        // get the smallest partition of the skeleton element
        gel->YoungestChildren(Sons);
#ifdef PZ_LOG
        if(logger.isDebugEnabled()){
            std::stringstream sout;
            sout << "For skeleton element " << gel->Index() << " subdomain " << subdomain[gel->Index()] << " has comp element " << (void *)gel->Reference() << std::endl;
            for(auto el : Sons) sout << "subel index " << el->Index() << " subdomain " <<
                subdomain[el->Index()] <<std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
       
        std::vector<TPZGeoEl*> sons;
        TPZGeoElSide fatside(gel, gel->NSides()-1);
        TPZGeoElSide neigfatside = fatside.Neighbour();
        while(fatside!=neigfatside){
            TPZCompEl *cel = neigfatside.Element()->Reference();
            int matId =neigfatside.Element()->MaterialId();
            int nconnects =-1;
            if(cel){
                nconnects=cel->NConnects();
            }
            if(matId==19 && nconnects>0){
                std::pair fatson = std::make_pair(neigfatside.Element(),sons);
                fatherAndSons.push_back(fatson);
                count++;
            }
            neigfatside = neigfatside.Neighbour();
        }
        for(auto songel: Sons){
#ifdef PZ_LOG
        if(logger.isDebugEnabled()){
            std::stringstream sout;
            sout << "For son element " << songel->Index() << " subdomain " << subdomain[songel->Index()] << " has comp element " << (void *)songel->Reference() << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
            TPZGeoElSide gelside(songel, songel->NSides()-1);
            TPZGeoElSide neig = gelside.Neighbour();
            std::vector<TPZCompEl *> fluxEls;
            while(neig!=gelside){
                TPZGeoEl *neigGel = neig.Element();
                int matId = neigGel->MaterialId();
                TPZCompEl *celneih = neigGel->Reference();
                int nconnects = -1;
                if(celneih){
                    nconnects=celneih->NConnects();
                }
                
//                TPZCompEl *neigCel = neigGel->Reference();
                if( matId==40){
                    fatherAndSons[count].second.push_back(neigGel);
//                    TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *>(neigCel);
                }
                neig=neig.Neighbour();
            }
        }
    }
}
void TPZReservoirTools::TakeElementsbyID(TPZGeoMesh *mGeometry, std::map<int, std::vector<TPZGeoEl* >> &interfaces, std::vector<int> &matIds){
    int nels = mGeometry->NElements();
    for(int iel =0; iel <nels; iel++){
        TPZGeoEl *gel = mGeometry->Element(iel);
        int matId = gel->MaterialId();
        for(auto mat: matIds){
            if(matId==mat){
                interfaces[matId].push_back(gel);
                break;
            }
        }
    }
}
/// group elements that share a connect with the basis elements
void TPZReservoirTools::GroupNeighbourElements(TPZCompMesh *cmesh, const std::set<int64_t> &seed_elements, std::set<int64_t> &groupindexes)
{
    TPZVec<int64_t> connectgroup(cmesh->NConnects(),-1);
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> elhandled(nel,0);
    for(auto el : seed_elements)
    {
        elhandled[el] = 1;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
        const int64_t index = elgr->Index();
        if(index < nel) elhandled[index] = 1;
        groupindexes.insert(index);
        TPZCompEl *cel = cmesh->Element(el);
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++)
        {
            int64_t cindex = cel->ConnectIndex(ic);
#ifdef PZDEBUG
            // the elements in seed_elements should not share any connect
            if(connectgroup[cindex] != -1) DebugStop();
#endif
            connectgroup[cindex] = index;
        }
        elgr->AddElement(cel);
    }
    for(int64_t el = 0; el < nel; el++)
    {
        if(elhandled[el]) continue;
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        std::set<int64_t> groups;
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++)
        {
            int64_t cindex = cel->ConnectIndex(ic);
            int64_t group = connectgroup[cindex];
            if(group != -1) groups.insert(group);
        }
        if(groups.size()>1) DebugStop();
        if(groups.size())
        {
            int64_t group = *groups.begin();
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cmesh->Element(group));
            elgr->AddElement(cel);
        }
    }
}
void TPZReservoirTools::TakeSeedElements(TPZCompMesh *cmesh, std::set<int64_t> &seed_elements){
    
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel){
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() == dim) {
            std::set<int64_t> seed;
            seed_elements.insert(el);
        }
    }
}

// duplicate the flux elements and put each of them in a separate subdomain
void TPZReservoirTools::PutFluxElementsinSubdomain(TPZCompMesh *fluxCmesh, TPZVec<int64_t> &subdomain, std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons)
{
    fluxCmesh->SetAllCreateFunctionsHDiv();
    fluxCmesh->ApproxSpace().CreateDisconnectedElements(false);
    TPZGeoMesh *gmesh = fluxCmesh->Reference();
    if(subdomain.size() != gmesh->NElements()) DebugStop();
    for (auto skel : fatherAndSons) {
//      auto gelskel = skel.first;
        auto &gelsons = skel.second;
        for(auto gelson : gelsons)
        {
            //int64_t sondomain = subdomain[gelson->Index()];
            TPZGeoElSide gelside(gelson);
            TPZGeoElSide neighbour = gelside.Neighbour();
            std::set<int64_t> domains;
            while (neighbour != gelside) {
                auto neighdomain = subdomain[neighbour.Element()->Index()];
                if(neighdomain != -1) domains.insert(neighdomain);
                neighbour = neighbour.Neighbour();
            }
            if(domains.size() != 2) DebugStop();
            TPZGeoElBC gelbc(gelside,gelson->MaterialId());
            TPZCompEl *cel = gelson->Reference();
            if(!cel || cel->NConnects() != 1) DebugStop();
          
            TPZStack<std::pair<TPZGeoElSide, TPZCompElSide>> loadstruct;
            loadstruct.Push({gelside,gelside.Reference()});
            gelside.Element()->ResetReference();
            {
                neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    loadstruct.Push({neighbour,neighbour.Reference()});
                    neighbour.Element()->ResetReference();
                    neighbour = neighbour.Neighbour();
                }
            }
            TPZCompEl *cel2 = fluxCmesh->ApproxSpace().CreateCompEl(gelbc.CreatedElement(), *fluxCmesh);
            int64_t c2index = cel2->ConnectIndex(0);
            TPZConnect &c2 = cel2->Connect(0);
            c2.ResetElConnected();
            c2.RemoveDepend();
            c2.Reset();
            fluxCmesh->ConnectVec().SetFree(c2index);
            cel2->SetConnectIndex(0, cel->ConnectIndex(0));
            if(cel2->NConnects() != 1) DebugStop();
            for(int64_t i = 0; i< loadstruct.size(); i++)
            {
                if (loadstruct[i].second) {
                    loadstruct[i].second.Element()->LoadElementReference();
                }
            }
            int64_t gelindex1 = gelson->Index();
            int64_t gelindex2 = gelbc.CreatedElement()->Index();
            int64_t dom1 = *domains.begin();
            int64_t dom2 = *domains.rbegin();
            {
                auto nel = gmesh->NElements();
                subdomain.resize(nel);
            }

            subdomain[gelindex1] = dom1;
            subdomain[gelindex2] = dom2;
            
        }
    }
}

/// push the connect with LagrangeLevel with highest sequence number to LagrangeDest and apply saddle permute
void TPZReservoirTools::PushConnectBackward(TPZCompMesh *cmesh, char LagrangeStart, char LagrangeDest)
{
    // apply the method to all subcompmeshes
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcmesh)
        {
            PushConnectBackward(subcmesh, LagrangeStart, LagrangeDest);
        }
    }
    
    
    // renumber the connects
    {
        TPZLinearAnalysis an(cmesh);
    }
    // look for the connect with lagrange level LagrangeStart and highest sequance number
    int64_t nc = cmesh->NConnects();
    int64_t pushconnect = -1;
    int64_t pushconnectSeqnum = -1;
    for (int64_t ic = 0; ic<nc; ic++) {
        TPZConnect &c = cmesh->ConnectVec()[ic];
        if(c.NElConnected() && c.LagrangeMultiplier() == LagrangeStart)
        {
            int64_t seqnum = c.SequenceNumber();
            if(seqnum > pushconnectSeqnum)
            {
                pushconnect = ic;
                pushconnectSeqnum = seqnum;
            }
        }
    }
    if(pushconnect == -1) return;
    TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh);
    TPZConnect &c = cmesh->ConnectVec()[pushconnect];
    c.SetLagrangeMultiplier(LagrangeDest);
    
    // Renumber the connect such that the ordering is consistent
    cmesh->SaddlePermute();
    
    // If the mesh is a SubCompMesh put the external connects last
    if(subcmesh)
    {
        subcmesh->PermuteExternalConnects();
        auto subanalysis = subcmesh->Analysis();
        auto solver = subanalysis->Solver();
        solver->ResetMatrix();
    }
}

/// Split the H(div) connects and create an HDivBound next to the faces with
/// It is assumed that the H(div) mesh is the mesh referenced by the geometric mesh
/// @param : gelIntersect geometric element with intersection matid.
void TPZReservoirTools::SplitConnect(TPZGeoEl *gelIntersect)
{
    TPZGeoMesh *gmesh = gelIntersect->Mesh();
    TPZCompMesh *cmesh = gmesh->Reference();
    if(!cmesh) DebugStop();
    // identify the two neighbours
    auto connected = FindHdiv(gelIntersect);
    // verify is the side hasn't been hybridized yet
    TPZCompEl *celIntersect = gelIntersect->Reference();
    if(celIntersect) DebugStop();
    // reset the reference for the second neighbour
    if(connected.second)
    {
        TPZCompEl *cel = connected.second.Element();
        TPZConnect &c = cel->Connect(0);
        c.DecrementElConnected();
        int64_t cindex = cmesh->AllocateNewConnect(c.NShape(), c.NState(), c.Order());
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        int locindex = intel->SideConnectLocId(0, connected.second.Side());
        intel->SetConnectIndex(locindex, cindex);
        TPZConnect &c2 = cel->Connect(locindex);
        c2.IncrementElConnected();
        TPZGeoElSide gelside = connected.second.Reference();
        gelside.Element()->ResetReference();
    }
    // create the hdivbound for the first element/side
    {
        TPZGeoElSide gelside(gelIntersect);
        TPZGeoElSide volside = connected.first.Reference();
        // we assume the mesh is "loaded"
        if(!volside.Element()->Reference()) DebugStop();
        MakeFirstNeighbour(volside, gelside);
        cmesh->ApproxSpace().CreateCompEl(gelIntersect, *cmesh);
        celIntersect = gelIntersect->Reference();
        volside.Element()->ResetReference();
        gelside.Element()->ResetReference();
    }
    // creat hdivbound for the second element/side
    if(connected.second)
    {
        TPZGeoElSide volside = connected.second.Reference();
        connected.second.Element()->LoadElementReference();
        TPZGeoElBC gbc(volside,gelIntersect->MaterialId());
        TPZGeoElSide gelside(gbc.CreatedElement());
        cmesh->ApproxSpace().CreateCompEl(gelside.Element(), *cmesh);
    }
    // reload the references
    celIntersect->LoadElementReference();
    connected.first.Element()->LoadElementReference();
#ifdef PZDEBUG
    auto gelside1 = connected.first.Reference();
    auto neigh1 = gelside1.Neighbour();
    if(neigh1.Element()->Dimension() != 1) DebugStop();
    auto celside2 = connected.second;
    if(celside2)
    {
        auto gelside2 = connected.second.Reference();
        auto neigh2 = gelside2.Neighbour();
        if(neigh2.Element()->Dimension() != 1) DebugStop();
    }
#endif
}

/// Find one or two TPZCompElSides which are neighbour of the geometric element and have a dimension one higher
std::pair<TPZCompElSide,TPZCompElSide> TPZReservoirTools::FindHdiv(TPZGeoEl *gelintersect)
{
    int dim =gelintersect->Dimension();
    TPZGeoElSide gelside(gelintersect);
    TPZGeoElSide neighbour = gelside.Neighbour();
    std::pair<TPZCompElSide,TPZCompElSide> result;
    int count = 0;
    while(neighbour != gelside)
    {
        TPZGeoEl *neighgel = neighbour.Element();
        if(neighgel->Reference() && neighgel->Dimension() == dim+1)
        {
            if(count == 0) result.first = neighbour.Reference();
            else if(count == 1) result.second = neighbour.Reference();
            else DebugStop();
            count++;
        }
        neighbour = neighbour.Neighbour();
    }
    if(count == 0) DebugStop();
    return result;
}

/// Make the bound element first neighbour of the volume element
void TPZReservoirTools::MakeFirstNeighbour(TPZGeoElSide &volume, TPZGeoElSide &bound)
{
#ifdef PZDEBUG
    if(!bound.IsNeighbour(volume))
    {
        DebugStop();
    }
#endif
    bound.RemoveConnectivity();
    volume.SetConnectivity(bound);
}
