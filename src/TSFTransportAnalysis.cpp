//
//  Created by Giovane Avancini on 04/11/25.
//

#include "TSFTransportAnalysis.h"
#include "TPZVTKGenerator.h"

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
  TPZSkylineNSymStructMatrix<STATE> matrix(fCompMesh);
#endif
  int n_threads = 0;
  matrix.SetNumThreads(n_threads);
  SetStructuralMatrix(matrix);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  SetSolver(step);

  std::cout << "Number of equations: " << fCompMesh->NEquations() << std::endl;
  std::cout << "Number of elements: " << fCompMesh->NElements() << std::endl;
}

void TSFTransportAnalysis::RunTimeStep(std::ostream &out) {
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
    TPZFMatrix<STATE> dsol = Solution();
    corr_norm = Norm(dsol);
    sol += dsol;
    cmesh->LoadSolution(sol);
    cmesh->TransferMultiphysicsSolution();
  }

  if (!converged) {
    out << "------Iterative method did not converge. res_norm: " << res_norm << " corr_norm: " << corr_norm << std::endl;
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

  if (fIsFirstAssemble) {//We call mother class Assemble() to compute the sparsity pattern, but the TSFTransportMaterial doesnt fill the matrix.
    AssembleMass();
    fIsFirstAssemble = false;
  }

  

  auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count() / 1000.;
  std::cout << "---------Time to assemble: " << total_time_ass << " seconds" << std::endl;
}

void TSFTransportAnalysis::AssembleMass() {

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
}

void TSFTransportAnalysis::SetLastStateSaturation() {
  // Implementation here
}

void TSFTransportAnalysis::SetLastStateDensities() {
  // Implementation here
}

void TSFTransportAnalysis::SetLastStateVolumeFactor() {
  // Implementation here
}