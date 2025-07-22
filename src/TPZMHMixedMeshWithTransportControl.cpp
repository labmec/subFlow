
#include "TPZMHMixedMeshWithTransportControl.h"
#include "TPZNullMaterial.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("iRMS.TPZMHMixedMeshWithTransportControl"));

#endif

void TPZMHMixedMeshWithTransportControl::BuildComputationalMesh(bool usersubstructure){
    
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    InsertPeriferalMaterialObjects();
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
    CreateTransport();
//    CreateAveragePressure();
    
    if(fNState > 1)
    {
        CreateRotationMesh();
    }
    
    CreateHDivPressureMHMMesh();
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG
    CheckMeshConsistency();
#endif
    //+
    
    TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(fCMesh.operator->());
    //Link Memory
    faproxspace->mMixedOperator= multcmesh;
    faproxspace->BuildTransportMultiPhysicsCompMesh();
    TPZCompMesh *cmeshtrans = faproxspace->GetTransportOperator();
    faproxspace->LinkMemory(multcmesh, cmeshtrans);
    //
    
    if (usersubstructure) {
        HideTheElements();
    }
    fNumeq = fCMesh->NEquations();
    
}
void TPZMHMixedMeshWithTransportControl::CreateHDivPressureMHMMesh()
{
    TPZManVector<TPZCompMesh *,2 > cmeshes(3);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fcmeshTransport.operator->();
   
    
    
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
    
    BuildMultiPhysicsMesh();
    TPZManVector<TPZCompMesh * ,3> meshvector;
    
    meshvector = cmeshes;
    
    // populate the connect to subdomain data structure for the multiphysics mesh
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
    //    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    //    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    // copy the solution of the atomic meshes to the multiphysics mesh
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
//    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
//    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);

    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
    return;
    
}


void TPZMHMixedMeshWithTransportControl::CreateTransport()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
//    gmesh->Print();
    gmesh->ResetReference();
    TPZCompMesh * cmeshtemp = new TPZCompMesh(gmesh);
    fcmeshTransport = cmeshtemp;
    // the pressure mesh should be empty when calling this method
    int64_t nskeletonconnects = fcmeshTransport->NConnects();
    if(nskeletonconnects != 0){
        DebugStop();
    }
    
    int porder = fpOrderInternal;
    // create and organize the pressure mesh
    // the pressure mesh is composed of discontinuous H1 elements
    TPZCompMesh * cmeshTransport = fcmeshTransport.operator->();
    gmesh->ResetReference();
    cmeshTransport->SetName("PressureMeshAverage");
    cmeshTransport->SetDimModel(gmesh->Dimension());
    cmeshTransport->SetAllCreateFunctionsDiscontinuous(); //AQUI
    cmeshTransport->ApproxSpace().CreateDisconnectedElements(true);
    cmeshTransport->SetDefaultOrder(0);
   
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
    cmeshTransport->InsertMaterialObject(volume);
    TPZNullMaterial<STATE> * volume2 = new TPZNullMaterial(2);
    cmeshTransport->InsertMaterialObject(volume2);
    matids.insert(1);
    matids.insert(2);
    cmeshTransport->AutoBuild(matids);
//    cmeshTransport->Print();
    cmeshTransport->ExpandSolution();
    
    if(1)
    {
        std::ofstream out("TransportCmesh.txt");
        cmeshTransport->Print(out);
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
    
    // the lagrange multiplier level is set to one
    int64_t nc = cmeshTransport->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshTransport->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    // associate the connects with the proper subdomain
    gmesh->ResetReference();
    int64_t nel = cmeshTransport->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshTransport->Element(el);
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
    
    {
        std::ofstream out("PressureFineMesh2AVERAGE.txt");
        cmeshTransport->Print(out);
    }
    
    return;
}
void TPZMHMixedMeshWithTransportControl::BuildMultiPhysicsMesh()
{
    if (fCMesh->NElements() != 0) {
        DebugStop();
    }
    fCMesh->SetAllCreateFunctionsMultiphysicElemWithMem();
//    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(fCMesh.operator->());
   
    int vecsize = 3;
    TPZManVector<TPZCompMesh *> meshvec(vecsize);
    meshvec[0] = fFluxMesh.operator->();
    meshvec[1] = fPressureFineMesh.operator->();
    meshvec[2] = fcmeshTransport.operator->();
//    meshvec[3] = fcmeshPressureAverg.operator->();
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
    nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
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
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZVec<int> active_approx_spaces(3,1);
    active_approx_spaces[2]=0;
    mphysics->BuildMultiphysicsSpaceWithMemory(active_approx_spaces, meshvec);
//    mphysics->BuildMultiphysicsSpace(active_approx_spaces, meshvec);
//    mphysics->BuildMultiphysicsSpace(meshvec,gelindexes);
    
}
