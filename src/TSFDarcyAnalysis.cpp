//
//  Created by Giovane Avancini on 02/09/25.
//

#include "TSFDarcyAnalysis.h"
#include "TPZVTKGenerator.h"

using namespace std;

TSFDarcyAnalysis::TSFDarcyAnalysis() : TPZLinearAnalysis() {}

TSFDarcyAnalysis::TSFDarcyAnalysis(TPZMultiphysicsCompMesh *cmesh, const RenumType &renumtype) : TPZLinearAnalysis(cmesh, renumtype) {}

TSFDarcyAnalysis::~TSFDarcyAnalysis() {}

void TSFDarcyAnalysis::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
}

TSFProblemData *TSFDarcyAnalysis::GetProblemData() {
  return fSimData;
}

int TSFDarcyAnalysis::GetNumberOfIterations() {
  return fKiteration;
}

void TSFDarcyAnalysis::Initialize() {
#ifdef PZ_USING_MKL
  TPZSSpStructMatrix<STATE> matrix(fCompMesh);
#else
  TPZSkylineStructMatrix<STATE> matrix(fCompMesh);
#endif
  int n_threads = fSimData->fTNumerics.fNThreadsDarcy;
  matrix.SetNumThreads(n_threads);
  SetStructuralMatrix(matrix);
  std::set<int> neumannMatids;
  FillNeumannBCMatids(neumannMatids);
  SetInitialSolution(neumannMatids);
  ApplyEquationFilter(neumannMatids);
  int nreducedeq = fStructMatrix->NReducedEquations();
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  SetSolver(step);

  std::cout << "Number of equations: " << fCompMesh->NEquations() << std::endl;
  std::cout << "Number of elements: " << fCompMesh->NElements() << std::endl;
}

void TSFDarcyAnalysis::RunTimeStep() {
  TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
  if (!cmesh) DebugStop();

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

void TSFDarcyAnalysis::PostProcessTimeStep(int dimToPost, int step) {
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

void TSFDarcyAnalysis::Assemble() {
  auto start_time_ass = std::chrono::steady_clock::now();

  TPZLinearAnalysis::Assemble();

  auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count() / 1000.;
  std::cout << "---------Time to assemble: " << total_time_ass << " seconds" << std::endl;
}

void TSFDarcyAnalysis::Solve() {
  auto start_time_solve = std::chrono::steady_clock::now();

  TPZLinearAnalysis::Solve();

  auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count() / 1000.;
  std::cout << "---------Time to solve: " << total_time_solve << " seconds" << std::endl;
}

void TSFDarcyAnalysis::FillNeumannBCMatids(std::set<int> &neumannMatids) {
  for (auto it = fSimData->fTBoundaryConditions.fBCDarcyMatIdToTypeValue.begin(); it != fSimData->fTBoundaryConditions.fBCDarcyMatIdToTypeValue.end(); ++it) {
    int matid = it->first;
    int bc_type = it->second.first;
    if (bc_type == 1)
      neumannMatids.insert(matid);
  }
}

void TSFDarcyAnalysis::SetInitialSolution(std::set<int> &neumannMatids) {
  fCompMesh->LoadReferences();
  TPZGeoMesh *gmesh = fCompMesh->Reference();

  TPZFMatrix<STATE> sol = fCompMesh->Solution();
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
        sol.PutVal(eq, 0, val);
    }
  }

  fCompMesh->LoadSolution(sol);
  fCompMesh->TransferMultiphysicsSolution();
  Solution() = sol;
}

void TSFDarcyAnalysis::ApplyEquationFilter(std::set<int> &neumannMatids) {
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