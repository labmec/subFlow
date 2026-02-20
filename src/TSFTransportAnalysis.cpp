//
//  Created by Giovane Avancini on 04/11/25.
//

#include "TSFTransportAnalysis.h"
#include "TPZVTKGenerator.h"
#include <cmath>

TSFTransportAnalysis::TSFTransportAnalysis() : TPZLinearAnalysis() {}

TSFTransportAnalysis::TSFTransportAnalysis(TPZCompMesh *cmesh, const RenumType &renumtype) : TPZLinearAnalysis(cmesh, renumtype) {}

TSFTransportAnalysis::~TSFTransportAnalysis() {}

void TSFTransportAnalysis::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
  fAlgebraicTransport.fCellsData.SetProblemData(simData);
}

TSFProblemData *TSFTransportAnalysis::GetProblemData() {
  return fSimData;
}

int TSFTransportAnalysis::GetNumberOfIterations() {
  return fKiteration;
}

void TSFTransportAnalysis::Initialize() {
#ifdef PZ_USING_MKL
  TPZSpStructMatrix<STATE> matrix(fCompMesh);
#else
  TPZSkylineNSymStructMatrix<STATE> matrix(fCompMesh); // this will break. I'm waiting for MUMPS
#endif
  int n_threads = 0;
  matrix.SetNumThreads(n_threads);
  SetStructuralMatrix(matrix);
  TPZStepSolver<REAL> step;
  step.SetDirect(ELU);
  SetSolver(step);
  TPZLinearAnalysis::Assemble(); // to create the sparsity pattern of fStructMatrix

  // Initialize fTransmissibilityMatrix with sparsity pattern from fStructMatrix
  auto baseMatrix = fStructMatrix->Create();

#ifdef PZ_USING_MKL
  TPZFYsmpMatrixPardiso<REAL> *sparseMatrix = dynamic_cast<TPZFYsmpMatrixPardiso<REAL> *>(baseMatrix);
#else
  TPZFYsmpMatrix<REAL> *sparseMatrix = dynamic_cast<TPZFYsmpMatrix<REAL> *>(baseMatrix);
#endif

#if PZDEBUG
  if (!sparseMatrix) {
    PZError << "ERROR: Could not cast matrix to TPZFYsmpMatrix" << std::endl;
    DebugStop();
  }
#endif

  // Extract sparsity pattern from base matrix
  TPZVec<int64_t> IA, JA;
  TPZVec<REAL> A;
  sparseMatrix->GetData(IA, JA, A);

  // Create fTransmissibilityMatrix with same dimensions and sparsity pattern
  int64_t neq = fCompMesh->NEquations();
#ifdef PZ_USING_MKL
  fTransmissibilityMatrix = new TPZFYsmpMatrixPardiso<REAL>(neq, neq);
#else
  fTransmissibilityMatrix = new TPZFYsmpMatrix<REAL>(neq, neq);
#endif

  fTransmissibilityMatrix->SetData(std::move(IA), std::move(JA), std::move(A));
  delete baseMatrix;

  std::cout << "Number of Transport equations: " << fCompMesh->NEquations() << std::endl;
  std::cout << "Number of Transport elements: " << fCompMesh->NElements() << std::endl;
}

void TSFTransportAnalysis::RunTimeStep(std::ostream &out) {
  auto SetNearZeroEntriesToZero = [](TPZFMatrix<REAL> &matrix) {
    // this is to avoid round-off errors that can compromise the vtk export
    const int64_t nrows = matrix.Rows();
    for (int64_t i = 0; i < nrows; i++) {
      if (IsZero(matrix(i, 0))) {
        matrix(i, 0) = 0.0;
      }
    }
  };

  TPZCompMesh *cmesh = Mesh();
  int matIter = fSimData->fTNumerics.fMaxIterTransport;
  REAL res_norm = 1.0;
  REAL corr_norm = 1.0;
  REAL res_tol = fSimData->fTNumerics.fResTolTransport;
  REAL corr_tol = fSimData->fTNumerics.fCorrTolTransport;
  bool converged = false;

  TPZFMatrix<STATE> sol = Solution();
  for (fKiteration = 0; fKiteration < matIter; fKiteration++) {
    Assemble();
    // auto matrix = MatrixSolver<REAL>().Matrix();
    // auto rhs = Rhs();
    // {
    //   matrix->Print("Matrix: ", out, EMathematicaInput);
    //   rhs.Print("RHS: ", out, EMathematicaInput);
    // }
    // Check residual convergence
    if (fKiteration > 0) {
      TPZFMatrix<STATE> &rhs = Rhs();
      res_norm = Norm(rhs);
      out << "------Newton iteration: " << fKiteration << std::endl;
      out << "---------Residual norm: " << res_norm << std::endl;
      out << "---------Correction norm: " << corr_norm << std::endl;
      if (res_norm < res_tol || corr_norm < corr_tol) {
        out << "------Iterative method converged with res_norm: " << res_norm << std::endl;
        out << "------Number of iterations = " << fKiteration << std::endl;
        converged = true;
        fSolution = sol;
        break;
      }
    }
    Solve();
    TPZFMatrix<STATE>& dsol = Solution();
    corr_norm = Norm(dsol);
    sol += dsol;
    SetNearZeroEntriesToZero(sol);
    cmesh->LoadSolution(sol);
    fAlgebraicTransport.fCellsData.UpdateSaturations(sol);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(fSimData->fTPetroPhysics.fKrModel);
  }

  if (!converged) {
    std::cout << "------Iterative method did not converge. res_norm: " << res_norm << " corr_norm: " << corr_norm << std::endl;
    fSolution = sol;
  }
}

void TSFTransportAnalysis::PostProcessTimeStep(int dimToPost, int step) {
  auto start_time = std::chrono::steady_clock::now();

  TPZStack<std::string, 10> scalnames;

  scalnames = fSimData->fTPostProcess.fScalnamesTransport;

  if (dimToPost < 0) {
    dimToPost = Mesh()->Reference()->Dimension();
  }

  std::string file = fSimData->fTPostProcess.fFileNameTransport;
  const int vtkRes = fSimData->fTPostProcess.fvtkResolution;

  const std::string plotfile = file.substr(0, file.find(".")); // sem o .vtk no final

  auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dimToPost);
  vtk.SetStep(step);
  int nthreads = fSimData->fTPostProcess.fNThreads;
  vtk.SetNThreads(nthreads);
  vtk.Do();

  auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.;
  std::cout << "Total time post process = " << total_time << " seconds" << std::endl;
}

void TSFTransportAnalysis::Assemble() {
  auto start_time_ass = std::chrono::steady_clock::now();

  if (fIsFirstAssemble) {
    AssembleMass();
    fIsFirstAssemble = false;
  }

  fTransmissibilityMatrix->Zero();
  TPZFMatrix<STATE> &rhs = Rhs();
  rhs.Zero();

  TPZVec<int64_t> IA, JA;
  TPZVec<REAL> A;
  fTransmissibilityMatrix->GetData(IA, JA, A);

  auto AddMatrixEntry = [&IA, &JA, &A](int64_t row, int64_t col, REAL value) {
    // Find the position of (row, col) in the sparse matrix
    for (int64_t idx = IA[row]; idx < IA[row + 1]; idx++) {
      if (JA[idx] == col) {
        A[idx] += value;
        return;
      }
    }
  };

  // Volume contribution
  int64_t ncells = fAlgebraicTransport.fCellsData.fVolume.size();
  for (int64_t icell = 0; icell < ncells; icell++) {
    int64_t eqid = fAlgebraicTransport.fCellsData.fEqNumber[icell];
    TPZFNMatrix<1, REAL> ef(1, 1, 0.0);
    fAlgebraicTransport.ContributeResidual(icell, ef);
    rhs(eqid) += ef(0, 0);
    AddMatrixEntry(eqid, eqid, fMassMatrix(eqid, 0));
  }

  // Interface contribution
  int interface_matid = fAlgebraicTransport.fCellsData.fSimData->fTGeometry.fInterface_material_id;
  int64_t n_internal_interfaces = fAlgebraicTransport.fInterfaceData[interface_matid].fFluxSign.size();
  for (int64_t interface = 0; interface < n_internal_interfaces; interface++) {
    std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[interface_matid].fLeftRightVolIndex[interface];
    int64_t left = lrindex.first;
    int64_t right = lrindex.second;
    int64_t lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
    int64_t righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];

    TPZFNMatrix<4, REAL> ek(2, 2, 0.0);
    TPZFNMatrix<2, REAL> ef(2, 1, 0.0);

    fAlgebraicTransport.ContributeInterface(interface, interface_matid, ek, ef);

    rhs(lefteq) += ef(0, 0);
    rhs(righteq) += ef(1, 0);
    AddMatrixEntry(lefteq, lefteq, ek(0, 0));
    AddMatrixEntry(lefteq, righteq, ek(0, 1));
    AddMatrixEntry(righteq, lefteq, ek(1, 0));
    AddMatrixEntry(righteq, righteq, ek(1, 1));
  }

  // BC contribution
  std::set<int> bc_matids;
  for (auto it = fAlgebraicTransport.fboundaryCMatVal.begin(); it != fAlgebraicTransport.fboundaryCMatVal.end(); it++) {
    bc_matids.insert(it->first);
  }
  int64_t n_inout_faces = 0;
  for (auto it = bc_matids.begin(); it != bc_matids.end(); it++) {
    n_inout_faces += fAlgebraicTransport.fInterfaceData[*it].fFluxSign.size();
  }

  for (auto it = bc_matids.begin(); it != bc_matids.end(); it++) {
    int64_t nbcs = fAlgebraicTransport.fInterfaceData[*it].fFluxSign.size();
    for (int64_t ibc = 0; ibc < nbcs; ibc++) {
      std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[*it].fLeftRightVolIndex[ibc];
      int64_t left = lrindex.first;
      int64_t lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];

      TPZFNMatrix<4, REAL> ek(1, 1, 0.0);
      TPZFNMatrix<2, REAL> ef(1, 1, 0.0);

      fAlgebraicTransport.ContributeBC(ibc, *it, ek, ef);

      rhs(lefteq) += ef(0, 0);
      AddMatrixEntry(lefteq, lefteq, ek(0, 0));
    }
  }

  fTransmissibilityMatrix->SetData(std::move(IA), std::move(JA), std::move(A));

  // Set the IA, JA, A structure directly in the base matrix of fStructMatrix
  auto baseMatrix = MatrixSolver<REAL>().Matrix();

#ifdef PZ_USING_MKL
  TPZFYsmpMatrixPardiso<REAL> *sparseMatrix = dynamic_cast<TPZFYsmpMatrixPardiso<REAL> *>(baseMatrix.operator->());
#else
  TPZFYsmpMatrix<REAL> *sparseMatrix = dynamic_cast<TPZFYsmpMatrix<REAL> *>(baseMatrix.operator->());
#endif
#ifdef PZDEBUG
  if (!sparseMatrix) {
    PZError << "ERROR: Could not cast matrix to TPZFYsmpMatrix" << std::endl;
    DebugStop();
  }
#endif

  sparseMatrix->CopyFrom(fTransmissibilityMatrix);

  auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count() / 1000.;
  std::cout << "---------Time to assemble: " << total_time_ass << " seconds" << std::endl;
}

void TSFTransportAnalysis::AssembleMass() {
  // Assemble fMassMatrix with diagonal structure
  int64_t ncells = fAlgebraicTransport.fCellsData.fVolume.size();
  fMassMatrix.Resize(ncells, 1);

  for (int64_t icell = 0; icell < ncells; icell++) {
    int64_t eqid = fAlgebraicTransport.fCellsData.fEqNumber[icell];
    TPZFNMatrix<1, REAL> ek(1, 1, 0.0), ef(1, 1, 0.0);
    fAlgebraicTransport.Contribute(icell, ek, ef);
    fMassMatrix(icell, 0) = ek(0, 0);
  }
}

void TSFTransportAnalysis::Solve() {
  auto start_time_solve = std::chrono::steady_clock::now();

  TPZLinearAnalysis::Solve();

  auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count() / 1000.;
  std::cout << "---------Time to solve: " << total_time_solve << " seconds" << std::endl;
}

void TSFTransportAnalysis::SetInitialSaturation() {
  TPZFMatrix<STATE> &sol = Solution();
  fAlgebraicTransport.fCellsData.UpdateSaturationsTo(sol);
  fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(fSimData->fTPetroPhysics.fKrModel);
  LoadSolution();
}

void TSFTransportAnalysis::SetLastStateVariables() {
  fAlgebraicTransport.fCellsData.SetLastStateVariables();
}

void TSFTransportAnalysis::SetLastStateSaturation() {
  fAlgebraicTransport.fCellsData.SetLastStateSaturation();
}

void TSFTransportAnalysis::SetLastStateDensities() {
  fAlgebraicTransport.fCellsData.SetLastStateDensities();
}

void TSFTransportAnalysis::SetLastStateVolumeFactor() {
  fAlgebraicTransport.fCellsData.SetLastStateVolumeFactor();
}