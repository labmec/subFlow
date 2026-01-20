//
//  Created by Giovane Avancini on 04/11/25.
//

#include "TSFTransportAnalysis.h"
#include "TPZVTKGenerator.h"

using namespace std;

TSFTransportAnalysis::TSFTransportAnalysis() : TPZLinearAnalysis() {}

TSFTransportAnalysis::TSFTransportAnalysis(TPZCompMesh *cmesh, const RenumType &renumtype) : TPZLinearAnalysis(cmesh, renumtype) {}

TSFTransportAnalysis::~TSFTransportAnalysis() {}

void TSFTransportAnalysis::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
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
  TPZSkylineNSymStructMatrix<STATE> matrix(fCompMesh);
#endif
  int n_threads = fSimData->fTNumerics.fNThreadsDarcy;
  matrix.SetNumThreads(n_threads);
  SetStructuralMatrix(matrix);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  SetSolver(step);

  std::cout << "Number of equations: " << fCompMesh->NEquations() << std::endl;
  std::cout << "Number of elements: " << fCompMesh->NElements() << std::endl;
}

void TSFTransportAnalysis::RunTimeStep() {
  TPZCompMesh *cmesh = Mesh();

  int matIter = fSimData->fTNumerics.fMaxIterDarcy;
  REAL res_norm = 1.0;
  REAL corr_norm = 1.0;
  REAL res_tol = fSimData->fTNumerics.fResTolDarcy;
  REAL corr_tol = fSimData->fTNumerics.fCorrTolDarcy;
  bool converged = false;

  TPZFMatrix<STATE> sol = Solution();
  for (fKiteration = 0; fKiteration < matIter; fKiteration++) {
    Assemble();

    // Check residual convergence
    if (fKiteration > 0) {
      TPZFMatrix<STATE> rhs = Rhs();
      res_norm = Norm(rhs);
      std::cout << "------Newton iteration: " << fKiteration << std::endl;
      std::cout << "---------Residual norm: " << res_norm << std::endl;
      std::cout << "---------Correction norm: " << corr_norm << std::endl;
      if (res_norm < res_tol || corr_norm < corr_tol) {
        std::cout << "------Iterative method converged with res_norm: " << res_norm << std::endl;
        std::cout << "------Number of iterations = " << fKiteration << std::endl;
        converged = true;
        fSolution = sol;
        break;
      }
    }
    Solve();
    TPZFMatrix<STATE> dsol = Solution();
    corr_norm = Norm(dsol);
    sol += dsol;
    cmesh->LoadSolution(sol);
    cmesh->TransferMultiphysicsSolution();
  }

  if (!converged) {
    std::cout << "------Iterative method did not converge. res_norm: " << res_norm << " corr_norm: " << corr_norm << std::endl;
  }
}

void TSFTransportAnalysis::PostProcessTimeStep(int dimToPost, int step) {
  auto start_time_pp = std::chrono::steady_clock::now();

  TPZStack<std::string, 10> scalnames, vecnames;

  scalnames = fSimData->fTPostProcess.fScalnamesDarcy;
  vecnames = fSimData->fTPostProcess.fVecnamesDarcy;

  int div = 0;
  if (dimToPost < 0) {
    dimToPost = Mesh()->Reference()->Dimension();
  }

  std::string file = fSimData->fTPostProcess.fFileNameDarcy;
  const int vtkRes = fSimData->fTPostProcess.fvtkResolution;

  const std::string plotfile = file.substr(0, file.find(".")); // sem o .vtk no final
  for (auto nm : vecnames) {
    scalnames.Push(nm);
  }

  auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dimToPost);
  vtk.SetStep(step);
  int nthreads = fSimData->fTPostProcess.fNThreads;
  vtk.SetNThreads(nthreads);
  vtk.Do();

  auto total_time_pp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_pp).count() / 1000.;
  cout << "Total time post process = " << total_time_pp << " seconds" << endl;
}

void TSFTransportAnalysis::Assemble() {
  auto start_time_ass = std::chrono::steady_clock::now();

  TPZLinearAnalysis::Assemble();

  auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count() / 1000.;
  std::cout << "---------Time to assemble: " << total_time_ass << " seconds" << std::endl;
}

void TSFTransportAnalysis::Solve() {
  auto start_time_solve = std::chrono::steady_clock::now();

  TPZLinearAnalysis::Solve();

  auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count() / 1000.;
  std::cout << "---------Time to solve: " << total_time_solve << " seconds" << std::endl;
}

void TSFTransportAnalysis::FillNeumannBCMatids(std::set<int> &neumannMatids) {
  for (auto it = fSimData->fTBoundaryConditions.fBCDarcyMatIdToTypeValue.begin(); it != fSimData->fTBoundaryConditions.fBCDarcyMatIdToTypeValue.end(); ++it) {
    int matid = it->first;
    int bc_type = it->second.first;
    if (bc_type == 1)
      neumannMatids.insert(matid);
  }
}

void TSFTransportAnalysis::SetInitialSolution(std::set<int> &neumannMatids) {
  fCompMesh->LoadReferences();
  TPZGeoMesh *gmesh = fCompMesh->Reference();

  TPZFMatrix<STATE> &cmesh_sol = fCompMesh->Solution();
  for (auto el : gmesh->ElementVec()) {
    int elMatID = el->MaterialId();

    if (neumannMatids.find(elMatID) == neumannMatids.end())
      continue;

    TPZCompEl *compEl = el->Reference();

    int64_t nConnects = compEl->NConnects();

    if (nConnects != 1)
      DebugStop();

    int64_t seq = compEl->Connect(0).SequenceNumber();
    int ncorner = el->NCornerNodes();
    REAL volume = el->Volume();

    auto firstEq = fCompMesh->Block().Position(seq);

    int64_t blockSize = fCompMesh->Block().Size(seq);

    REAL val = fSimData->fTBoundaryConditions.fBCDarcyMatIdToTypeValue[elMatID].second;
    val *= volume / ncorner;
    for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++) {
      if (eq - firstEq < ncorner)
        cmesh_sol.PutVal(eq, 0, val);
    }
  }

  fCompMesh->TransferMultiphysicsSolution();

  // When the internal dofs are condensed, the analysis solution size is different from the cmesh solution size
  // Analysis only holds the independent equations, while cmesh holds all equations
  // The independent equations are stored first in the cmesh solution vector, so we just need to copy them to the analysis solution
  int cmesh_neq = fCompMesh->NEquations();
  TPZFMatrix<STATE> &sol = Solution();
  for (int i = 0; i < cmesh_neq; i++) {
    sol.PutVal(i, 0, cmesh_sol.GetVal(i, 0));
  }
}

void TSFTransportAnalysis::ApplyEquationFilter(std::set<int> &neumannMatids) {
  fCompMesh->LoadReferences();
  std::set<int64_t> removeEquations;
  TPZGeoMesh *gmesh = fCompMesh->Reference();
  TPZFMatrix<STATE> sol = fCompMesh->Solution();
  for (auto el : gmesh->ElementVec()) {
    int elMatID = el->MaterialId();

    if (neumannMatids.find(elMatID) == neumannMatids.end()) continue;

    TPZCompEl *compEl = el->Reference();

    int64_t nConnects = compEl->NConnects();

    if (nConnects != 1)
      DebugStop();

    int64_t seq = compEl->Connect(0).SequenceNumber();
    int ncorner = el->NCornerNodes();

    auto firstEq = fCompMesh->Block().Position(seq);

    int64_t blockSize = fCompMesh->Block().Size(seq);

    for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++) {
      removeEquations.insert(eq);
    }
  }

  TPZEquationFilter filter(fCompMesh->NEquations());
  filter.ExcludeEquations(removeEquations);
  fStructMatrix->EquationFilter() = filter;
}