//
//  Created by Giovane Avancini on 23/12/25.
//

#include "TSFSFIAnalysis.h"

TSFSFIAnalysis::TSFSFIAnalysis() {
  fDarcyAnalysis = nullptr;
  fTransportAnalysis = nullptr;
  fKiteration = 0;
  fSimData = nullptr;
  fDarcySolution.Resize(0, 0);
  fTransportSolution.Resize(0, 0);
}

TSFSFIAnalysis::TSFSFIAnalysis(TPZMultiphysicsCompMesh *darcy_cmesh, TPZCompMesh *transport_cmesh, const RenumType &renumtype) {
  fDarcyAnalysis = TSFDarcyAnalysis(darcy_cmesh, renumtype);
  fTransportAnalysis = TSFTransportAnalysis(transport_cmesh, renumtype);
  fDataTransfer.SetMeshes(darcy_cmesh, transport_cmesh);
  fKiteration = 0;
  fSimData = nullptr;
  fDarcySolution.Resize(0, 0);
  fTransportSolution.Resize(0, 0);
}

TSFSFIAnalysis::~TSFSFIAnalysis() {}

void TSFSFIAnalysis::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
  fDarcyAnalysis.SetProblemData(simData);
  fTransportAnalysis.SetProblemData(simData);
}

void TSFSFIAnalysis::Initialize() {
  fDarcyAnalysis.Initialize();
  fTransportAnalysis.Initialize();
  fDataTransfer.Initialize();
  fDataTransfer.InitializeAlgebraicTransport(fTransportAnalysis.fAlgebraicTransport);
  fDarcyAnalysis.SetInitialSolution();
  fTransportAnalysis.SetInitialSaturation();
  TransferTransportToDarcy();
  TransferDarcyToTransport();
  UpdateLastStateVariables();
}

void TSFSFIAnalysis::PostProcessTimeStep(const int type, const int dim, int step) {
  auto start_time = std::chrono::steady_clock::now();

  if (type == 0) {
    fDarcyAnalysis.PostProcessTimeStep(dim, step);
    fTransportAnalysis.PostProcessTimeStep(dim, step);
  }
  if (type == 1) {
    fDarcyAnalysis.PostProcessTimeStep(dim, step);
  }
  if (type == 2) {
    fTransportAnalysis.PostProcessTimeStep(dim, step);
  }

  auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.;
  std::cout << "Total time post process = " << total_time << " seconds" << std::endl;
}

void TSFSFIAnalysis::Run(std::ostream &out) {
  const int dim = fSimData->fTGeometry.fDimension;
  const int maxIter = fSimData->fTNumerics.fMaxIterSFI;
  const REAL tol = fSimData->fTNumerics.fTolSFI;
  const REAL norm = 1.0;
  fKiteration = 0;
  const int nsteps = fSimData->fTNumerics.fNSteps;
  const REAL dt = fSimData->fTNumerics.fDt;

  // Instants where post-processing will be done
  TPZStack<REAL, 100> postProcessTimes = fSimData->fTPostProcess.fVecReportingTimes;
  REAL time = 0.0;
  REAL nextPostProcessTime = postProcessTimes[0];

  PostProcessTimeStep(fSimData->fTPostProcess.fProblemTypeInit, dim, 0);

  // Time stepping loop
  for (int tstep = 1; tstep <= nsteps; tstep++) {
    out << "\n=================================================================" << std::endl;
    out << "-------------------------- TIME Step " << tstep << " --------------------------" << std::endl;
    out << "=================================================================" << std::endl;
    time += dt;
    out << "Time = " << time << std::endl;
    fDarcyAnalysis.SetTime(time);
    fTransportAnalysis.SetTime(time);
    RunTimeStep(out);
    // Post-processing
    if (tstep == 1) {
      PostProcessTimeStep(fSimData->fTPostProcess.fProblemTypeInit, dim, tstep);
      nextPostProcessTime = postProcessTimes[tstep];
    } else if (time >= nextPostProcessTime - 1e-8) {
      PostProcessTimeStep(fSimData->fTPostProcess.fProblemType, dim, tstep);
      nextPostProcessTime = postProcessTimes[tstep];
    }
    // Check Mass Balance
    fTransportAnalysis.fAlgebraicTransport.CheckMassBalance(time, out);
  }
  out << "-------------------- Simulation Finished --------------------" << std::endl;
  std::ofstream darcy_sol("darcy_solution.txt");
  fDarcyAnalysis.Mesh()->Print(darcy_sol);
}

void TSFSFIAnalysis::RunTimeStep(std::ostream &out) {
  // Updating Darcy and Transport solutions
  fDarcySolution = fDarcyAnalysis.Solution();
  fTransportSolution = fTransportAnalysis.Solution();

  const int maxIter = fSimData->fTNumerics.fMaxIterSFI;
  const REAL tol = fSimData->fTNumerics.fTolSFI;
  REAL darcy_norm = 1.0;
  REAL transport_norm = 1.0;

  for (fKiteration = 1; fKiteration <= maxIter; fKiteration++) {
    out << "-------------------- SFI Iteration " << fKiteration << " --------------------" << std::endl;

    // Solve Darcy Problem
    if (fShouldSolveDarcy) {
      TransferTransportToDarcy(); // Update Darcy problem with new saturation
      out << "Solving Darcy Problem..." << std::endl;
      fDarcyAnalysis.RunTimeStep(out);
      TransferDarcyToTransport(); // Update transport problem with new densities and fluxes
      if (fSimData->fTNumerics.fIsLinearTrace) fShouldSolveDarcy = false;
    }

    // Solve Transport Problem
    out << "Solving Transport Problem..." << std::endl;
    fTransportAnalysis.RunTimeStep(out);

    // Computing norms
    darcy_norm = Norm(fDarcySolution - fDarcyAnalysis.Solution());
    transport_norm = Norm(fTransportSolution - fTransportAnalysis.Solution());
    out << "------Darcy norm: " << darcy_norm << " Transport norm: " << transport_norm << std::endl;

    // Updating Darcy and Transport solutions
    fDarcySolution = fDarcyAnalysis.Solution();
    fTransportSolution = fTransportAnalysis.Solution();

    // Check convergence
    if (darcy_norm < tol && transport_norm < tol) {
      out << "------SFI Converged in " << fKiteration << " iterations. Darcy norm: " << darcy_norm << " Transport norm: " << transport_norm << std::endl;
      break;
    } else if (fKiteration == maxIter) {
      out << "------SFI failed to converge. Darcy norm: " << darcy_norm << " Transport norm: " << transport_norm << std::endl;
    }
  }
  UpdateLastStateVariables();
}

void TSFSFIAnalysis::TransferTransportToDarcy() {
  fDataTransfer.TransferSaturation();
}

void TSFSFIAnalysis::TransferDarcyToTransport() {
  fDataTransfer.TransferPressures();
  fDataTransfer.TransferDarcyMeshMultiplyingCoefficients();
  fTransportAnalysis.fAlgebraicTransport.UpdateInterfacesIntegratedFlux(fSimData->fTGeometry.fInterface_material_id);

  std::set<int> bc_matids;
  for (auto it = fTransportAnalysis.fAlgebraicTransport.fboundaryCMatVal.begin(); it != fTransportAnalysis.fAlgebraicTransport.fboundaryCMatVal.end(); it++) {
    bc_matids.insert(it->first);
  }
  for (auto &bc_matid : bc_matids) {
    fTransportAnalysis.fAlgebraicTransport.UpdateInterfacesIntegratedFlux(bc_matid);
  }

  fTransportAnalysis.fAlgebraicTransport.fCellsData.UpdateDensitiesAndVolumeFactor();
}

void TSFSFIAnalysis::UpdateLastStateVariables() {
  fTransportAnalysis.SetLastStateVariables();
  fDarcyAnalysis.SetLastStateVariables();
}