//
//  TMRSTransportAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSTransportAnalysis.h"
#include "pzfunction.h"
#include "TPZTracerFlow.h"
#include "pzvtkmesh.h"

// Uses the new vtk function developed by Fran
#define USENEWVTK

#ifdef USENEWVTK
#include "TPZVTKGenerator.h"
#endif
#include "TPZVTKGeoMesh.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif
#include "TPZSpStructMatrix_Eigen.h"
#include "TPZSpMatrixEigen.h"
TMRSTransportAnalysis::TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::~TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::TMRSTransportAnalysis(TPZCompMesh * cmesh_mult,
                                             const RenumType& renumtype) : TPZLinearAnalysis(cmesh_mult, renumtype){
    
    TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh_mult);
    if(cmesh){
        m_soltransportTransfer.BuildTransferData(cmesh_mult);
    }
}

void TMRSTransportAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
}

TMRSDataTransfer * TMRSTransportAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSTransportAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSTransportAnalysis::Configure(int n_threads, bool UsePardiso_Q){
    
    if (UsePardiso_Q) {
        TPZSpStructMatrix<> matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELU);
//        SetSolver(step);
        
//        TPZSpStructMatrixEigen matrix(Mesh());
//        matrix.SetNumThreads(n_threads);
//        SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELU);
//        SetSolver(step);
        
    }else{
        TPZSkylineNSymStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
}

void TMRSTransportAnalysis::Assemble(){
    if (m_parallel_execution_Q) {
        fTransportSpMatrix->Assemble();
//        Assemble_parallel();
    }else{
        Assemble_serial();
    }
}

void TMRSTransportAnalysis::Assemble_serial(){
//    TPZAnalysis::Assemble();
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
    if(!this->fSolver){
        DebugStop();
    }
    TPZMatrix<STATE> *mat = 0;
    TPZMatrixSolver<STATE> *matsol = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
    if(matsol && !matsol->Matrix())
    {
        TPZBaseMatrix *othermat = fStructMatrix->Create();
        mat = dynamic_cast<TPZMatrix<STATE> *>(othermat);
        matsol->SetMatrix(mat);
    }
    else if(matsol)
    {
        mat = matsol->Matrix().operator->();
    }
    else
    {
        DebugStop();
    }
//    mat->Redim(ncells, ncells);
//    Rhs().Redim(ncells,1);
    mat->Zero();
    Rhs().Zero();
    TPZFMatrix<STATE> &rhs = Rhs();
    //Volumetric Elements
    for (int ivol = 0; ivol<ncells; ivol++) {
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1,1);
        fAlgebraicTransport.Contribute(ivol, elmat, ef);
        mat->AddSub(eqindex, eqindex, elmat);
        rhs.AddSub(eqindex, 0,ef);
    }
   
    //Interface Elements
    int interID = m_sim_data->mTGeometry.mInterface_material_id;
    int ninterfaces = fAlgebraicTransport.fInterfaceData[interID].fFluxSign.size();
    if (ninterfaces<1) {
        DebugStop();
    }
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[interID].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];
        TPZVec<int64_t> destinationindex(2);
        destinationindex[0]=lefteq;
        destinationindex[1]=righteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterface(interf,elmat, ef);
        mat->AddKel(elmat, destinationindex);
        rhs.AddFel(ef, destinationindex);
    }
    
    int inlet_mat_id = -2;
    //INLET
    ninterfaces = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fFluxSign.size();
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_mat_id].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
//        int right = lrindex.second;
        TPZVec<int64_t> destinationindex(1);
        destinationindex[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(interf,ef);
        rhs.AddFel(ef, destinationindex);
    }
    
    int outlet_mat_id = -4;
    //outlet
    ninterfaces = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fFluxSign.size();
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[outlet_mat_id].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        TPZVec<int64_t> destinationindex(1);
        destinationindex[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterface(interf,elmat,ef);
        mat->AddKel(elmat, destinationindex);
        rhs.AddFel(ef, destinationindex);
    }
}

void TMRSTransportAnalysis::AssembleResidual(){
    if (m_parallel_execution_Q) {
        AssembleResidual_Eigen();
    }else{
        AssembleResidual_serial();
    }
}

void TMRSTransportAnalysis::RunTimeStep(){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    TPZCompMesh * cmesh2 = dynamic_cast<TPZCompMesh *>(Mesh());
    if (!cmesh &&!cmesh2) {
        DebugStop();
    }
    
    int n = m_sim_data->mTNumerics.m_max_iter_transport;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_transport;
    TPZFMatrix<STATE> dx(Solution()),x(Solution());
    TPZFMatrix<STATE> correction(Solution());
    correction.Zero();
    
    corr_norm = ComputeInitialGuess(x); // from the linear problem (tangent and residue)
//    bool QN_converge_Q = QuasiNewtonSteps(x,20); // assuming linear operator (tangent)
//    if(QN_converge_Q){
//        return;
//    }

    //Linear problem Benchmark
    res_norm = Norm(Rhs());
    if(res_norm < res_tol && corr_norm < corr_tol){
        std::cout << "Transport operator: Converged - (InitialGuess)" << std::endl;
        std::cout << "Number of iterations = " << 1 << std::endl;
        std::cout << "residue norm = " << res_norm << std::endl;
        return;
    }
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
       
        NewtonIteration();
        dx = Solution();
//        std::cout<<"Sol Correct: "<<std::endl;
        x += dx;
        

        LoadSolution(x);
//        cmesh->LoadSolutionFromMultiPhysics();
        // PostProcessTimeStep();
        fAlgebraicTransport.fCellsData.UpdateSaturations(x);
        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    
        AssembleResidual();
        corr_norm = Norm(dx);
        res_norm = Norm(Rhs());
        
#ifdef PZDEBUG
        {
            if(std::isnan(corr_norm) || std::isnan(res_norm))
            {
                DebugStop();
            }
        }
#endif
        std::cout << "res_norm " << res_norm << " corr_norm " << corr_norm << std::endl;
        stop_criterion_Q = (res_norm < res_tol);
        stop_criterion_corr_Q = (corr_norm < corr_tol);
        if (stop_criterion_Q && stop_criterion_corr_Q) {
            std::cout << "Transport operator: Converged" << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            std::cout << "residue norm = " << res_norm << std::endl;
            break;
        }

    }
    if (!stop_criterion_Q && !stop_criterion_corr_Q) DebugStop(); //failed to converge
}

REAL TMRSTransportAnalysis::ComputeInitialGuess(TPZFMatrix<STATE> &x){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    TPZCompMesh * cmesh2 = dynamic_cast<TPZCompMesh *>(Mesh());
    if (!cmesh &&!cmesh2) {
        DebugStop();
    }
    
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    
    LoadSolution(x);
    if(cmesh){
        cmesh->LoadSolutionFromMultiPhysics();
    }
    else{
        cmesh2->LoadSolution(x);
    }
 
    
    NewtonIteration();
    auto dx = Solution();
    REAL corr_norm = Norm(dx);
    x += dx;
    
//    std::cout<<"SOLUTION: "<<std::endl;
//    std::cout<<x<<std::endl;
    LoadSolution(x);
    if(cmesh){
        cmesh->LoadSolutionFromMultiPhysics();
    }
    else{
        cmesh2->LoadSolution(x);
    }
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    AssembleResidual();
    REAL res_norm = Norm(Rhs());
    std::cout << "Initial guess residue norm : " <<  res_norm << ", corr_norm : " << corr_norm << std::endl;

    return corr_norm;
}

bool TMRSTransportAnalysis::QuasiNewtonSteps(TPZFMatrix<STATE> &x, int n){
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    std::cout << "Quasi-Newton process : " <<  std::endl;
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        NewtonIteration();
        
        x += Solution();
        
#ifdef PZDEBUG
        {
            REAL norm = Norm(x);
            if(std::isnan(norm))
            {
                DebugStop();
            }
        }
#endif
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        
        
        if (m_sim_data->mTPetroPhysics.mKrModel == 0) {
            fAlgebraicTransport.fCellsData.UpdateSaturations(x);
            fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(0);
        }else{
            fAlgebraicTransport.fCellsData.UpdateSaturations(x);
            fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambdaQuasiNewton();
        }
        
        AssembleResidual();
        REAL res_norm = Norm(Rhs());
        std::cout << " Residue norm : " <<  res_norm << std::endl;
        
        res_norm = Norm(Rhs());
        
#ifdef PZDEBUG
        if(std::isnan(res_norm))
        {
            DebugStop();
        }
#endif
        
        bool stop_criterion_Q = (res_norm < res_tol);
        if (stop_criterion_Q) {
            std::cout << "Transport operator: Converged" << std::endl;
            std::cout << "Quasi-Newton iterations = " << m_k_iteration << std::endl;
            std::cout << "Residue norm = " << res_norm << std::endl;
            return true;
        }
    }
    
    return false;
}

void TMRSTransportAnalysis::NewtonIteration(){
    
    if (m_parallel_execution_Q) {
        NewtonIteration_Eigen();
    }else{
        NewtonIteration_serial();
    }
}

void TMRSTransportAnalysis::NewtonIteration_serial(){

    fTransportSpMatrix->Assemble();
    fTransportSpMatrix->Solve();

    auto ds = fTransportSpMatrix->Solution();
    TPZFMatrix<STATE> &sol = Solution();

    for (int i = 0; i < ds.rows(); i++)
    {
        sol(i, 0) = ds(i, 0);
    }
#ifdef PZDEBUG
    STATE norm = Norm(Solution());
    if (std::isnan(norm)) DebugStop();
#endif
}

void TMRSTransportAnalysis::AnalyzePattern(){
//    Assemble_mass_parallel();
//    Assemble_parallel();
//    m_transmissibility += m_mass;
//    m_analysis.analyzePattern(m_transmissibility);
    fTransportSpMatrix = new TPZAnalysisAuxEigen(&fAlgebraicTransport);
    fTransportSpMatrix->SetAlgebraicTransport(&fAlgebraicTransport);
    fTransportSpMatrix->AnalyzePattern();
    // Because in some parts this objects are needed.
    Solution().Resize(fTransportSpMatrix->NRhsRows(), 1);
    Rhs().Resize(fTransportSpMatrix->NRhsRows(), 1);
    

}

void TMRSTransportAnalysis::NewtonIteration_Eigen(){
    
//    Assemble_parallel();
    fTransportSpMatrix->Assemble();
    
#ifdef PZDEBUG
//    {
//        auto norm = fTransportSpMatrix->RhsNorm();
//        if(std::isnan(norm))
//        {
//            DebugStop();
//        }
//    }
#endif
    std::cout<<"\n ---------------------- Solve Transport----------------------"<<std::endl;
    fTransportSpMatrix->Solve();
    
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> ds = fTransportSpMatrix->Solution();
    assert(Solution().Rows() == ds.rows());
    TPZFMatrix<STATE> &sol = Solution();
#ifdef USING_TBB
    tbb::parallel_for(size_t(0), size_t(ds.rows()), size_t(1),
                      [&ds, &sol] (size_t & i){
        sol(i,0) = ds(i,0);
        }
    );
#else
    for (int i = 0; i < ds.rows(); i++) {
        sol(i,0) = ds(i,0);
    }
#endif
#ifdef PZDEBUG
    STATE norm = Norm(Solution());
    if(std::isnan(norm)) DebugStop();
#endif
}

void TMRSTransportAnalysis::AssembleResidual_serial()
{

    if (0)
    {
        int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
        if (!this->fSolver)
        {
            DebugStop();
        }
        Rhs().Resize(ncells, 1);
        Rhs().Zero();
        TPZFMatrix<STATE> &rhs = Rhs();

        // Volumetric Elements
        for (int ivol = 0; ivol < ncells; ivol++)
        {
            int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
            TPZFMatrix<double> ef;
            ef.Resize(1, 1);
            fAlgebraicTransport.ContributeResidual(ivol, ef);
            rhs.AddSub(eqindex, 0, ef);
        }

        // Interface Elements
        int interID = m_sim_data->mTGeometry.mInterface_material_id;
        int ninterfaces = fAlgebraicTransport.fInterfaceData[interID].fFluxSign.size();
        if (ninterfaces < 1)
        {
            DebugStop();
        }
        for (int interf = 0; interf < ninterfaces; interf++)
        {
            std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[interID].fLeftRightVolIndex[interf];
            int left = lrindex.first;
            int right = lrindex.second;
            int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
            int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];
            TPZVec<int64_t> destinationindex(2);
            destinationindex[0] = lefteq;
            destinationindex[1] = righteq;
            TPZFMatrix<double> ef;
            ef.Resize(2, 1);
            fAlgebraicTransport.ContributeInterfaceResidual(interf, ef);
            rhs.AddFel(ef, destinationindex);
        }

        int inlet_mat_id = -2;
        // INLET
        ninterfaces = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fFluxSign.size();
        for (int interf = 0; interf < ninterfaces; interf++)
        {
            std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fLeftRightVolIndex[interf];
            int left = lrindex.first;
            int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
            TPZVec<int64_t> destinationindex(1);
            destinationindex[0] = lefteq;
            TPZFMatrix<double> ef;
            ef.Resize(1, 1);
            fAlgebraicTransport.ContributeBCInletInterface(interf, ef);
            rhs.AddFel(ef, destinationindex);
        }

        int outlet_mat_id = -4;
        // outlet
        ninterfaces = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fFluxSign.size();
        for (int interf = 0; interf < ninterfaces; interf++)
        {
            std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fLeftRightVolIndex[interf];
            int left = lrindex.first;
            int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
            TPZVec<int64_t> destinationindex(1);
            destinationindex[0] = lefteq;
            TPZFMatrix<REAL> ef;
            ef.Resize(1, 1);
            fAlgebraicTransport.ContributeBCOutletInterfaceResidual(interf, ef);
            rhs.AddFel(ef, destinationindex);
        }
    }

    fTransportSpMatrix->AssembleResidual();
    TPZFMatrix<STATE> &rhs = fRhs;
    auto r = fTransportSpMatrix->Rhs().toDense();
    for (int i = 0; i < r.rows(); i++)
    {
        rhs(i, 0) = r(i, 0);
    }
}

void TMRSTransportAnalysis::AssembleResidual_Eigen(){
 
    fTransportSpMatrix->AssembleResidual();
    TPZFMatrix<STATE> &rhs = fRhs;
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> r = fTransportSpMatrix->Rhs().toDense();
    assert(Rhs().Rows() == r.rows());
    #ifdef USING_TBB2
        tbb::parallel_for(size_t(0), size_t(r.rows()), size_t(1),
            [this,&r] (size_t & i){
             rhs(i,0) = r(i,0);
            }
        );
    #else
        for (int i = 0; i < r.rows(); i++) {
            rhs(i,0) = r(i,0);
        }
//    std::cout<<"RHS= "<<Rhs()<<std::endl;
    #endif
    
}


void TMRSTransportAnalysis::PostProcessTimeStep(){
    
    TPZStack<std::string,10> scalnames, vecnames;
    scalnames = m_sim_data->mTPostProcess.m_scalnamesTransport;
    constexpr int vtkRes{0}; //resolucao do vtk
    int dim = Mesh()->Reference()->Dimension();
    std::string file = m_sim_data->mTPostProcess.m_file_name_transport;
//    std::ofstream file2("transport"+std::to_string(fpostprocessindex)+".vtk");
    
//    
//    std::set<int> matidsToPost;
//    std::map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
//    if (m_sim_data->mTGeometry.isThereFracture()) {
//        for (it = m_sim_data->mTFracProperties.m_fracprops.begin(); it != m_sim_data->mTFracProperties.m_fracprops.end(); it++)
//        {
//            int matfracid = it->first;
//            matidsToPost.insert(matfracid);
//        }
//        std::string file_frac("fracture_s.vtk");
//        
//#ifdef USENEWVTK
//        const std::string plotfile = file_frac.substr(0, file.find("."));
//        for (auto nm : vecnames) {
//            scalnames.Push(nm);
//        }
//        auto vtk = TPZVTKGenerator(fCompMesh, matidsToPost, scalnames, plotfile, vtkRes);
//        vtk.SetNThreads(8);
//        vtk.Do();
//#else
//        DefineGraphMesh(2,matidsToPost,scalnames,vecnames,file_frac);
//        PostProcess(vtkRes,2);
//#endif
//    }
//    
//    if(m_sim_data->mTFracIntersectProperties.isThereFractureIntersection()){
//        matidsToPost.clear();
//        matidsToPost.insert(m_sim_data->mTGeometry.m_pressureMatId);
//        std::string file_frac2("fracture_s1d.vtk");
//        
//#ifdef USENEWVTK
//        const std::string plotfile = file_frac2.substr(0, file.find("."));
//        for (auto nm : vecnames) {
//            scalnames.Push(nm);
//        }
//        auto vtk = TPZVTKGenerator(fCompMesh, matidsToPost, scalnames, plotfile, vtkRes);
//        vtk.SetNThreads(8);
//        vtk.Do();
//#else
//        DefineGraphMesh(1,matidsToPost,scalnames,vecnames,file_frac2);
//        PostProcess(vtkRes,1);
//#endif
//    }
//    
//#ifdef USENEWVTK
//    const std::string plotfile = file.substr(0, file.find(".")); //sem o .vtk no final
//    for (auto nm : vecnames) {
//        scalnames.Push(nm);
//    }
//    
//    auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dim);
//    vtk.SetNThreads(8);
//    vtk.Do();
//#else
//    DefineGraphMesh(dim,scalnames,vecnames,file);
//    PostProcess(vtkRes,dim);
//#endif
    int nels = fGeoMesh->NElements();
    fCompMesh->LoadReferences();
    TPZVec<TPZVec<REAL>> elData(nels);
    for (int i = 0; i < nels; i++) {
        elData[i].Resize(2, 0.0); // two scalars: saturation and density
    }
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
    for (int icell = 0; icell< ncells; icell++) {
        int indexGeo = fAlgebraicTransport.fCellsData.fGeoIndex[icell];
        REAL sat= fAlgebraicTransport.fCellsData.fSaturation[icell];
        REAL airdensity = fAlgebraicTransport.fCellsData.fDensityOil[icell];
        auto center = fAlgebraicTransport.fCellsData.fCenterCoordinate[icell];
        TPZGeoEl *gel = fCompMesh->Reference()->Element(indexGeo);
        if(!gel) DebugStop();
        
#ifdef PZDEBUG
        TPZManVector<REAL,3> xcenter(3,0.), xi(gel->Dimension(),0.);
        gel->CenterPoint(gel->NSides()-1, xi);
        gel->X(xi, xcenter);
        REAL diff = 0;
        std::vector<REAL> celcenter = fAlgebraicTransport.fCellsData.fCenterCoordinate[icell];
        for(int i=0; i<3; i++)
        {
            diff += fabs(xcenter[i]-celcenter[i]);
        }
        if(diff > 1.e-16)
        {
            DebugStop();
        }
# endif
        if(fabs(sat) < 1e-20) sat = 0.;
        elData[indexGeo][0]=sat;
        elData[indexGeo][1]=airdensity;
    }
    
    std::string fileAdjusted = file.substr(0,file.find("vtk")) + std::to_string(fpostprocessindex) + ".vtk";
    char* filename = const_cast<char *>(fileAdjusted.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(fCompMesh->Reference(), filename, elData);
    fpostprocessindex++;
}

void TMRSTransportAnalysis::UpdateInitialSolutionFromCellsData(){
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    TPZCompMesh * cmesh2 = dynamic_cast<TPZCompMesh *>(Mesh());
    if (!cmesh && !cmesh2) {
        DebugStop();
    }
    fAlgebraicTransport.fCellsData.UpdateSaturationsTo(Solution());
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    LoadSolution();
    //REVISAR JOSE EN EL CASO DE NO MULTIFISICO
//    cmesh->LoadSolutionFromMultiPhysics();
}
