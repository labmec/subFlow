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

  if (type == 0) { // Both
    fDarcyAnalysis.PostProcessTimeStep(dim, step);
    fTransportAnalysis.PostProcessTimeStep(dim, step);
  }
  else if (type == 1) { // Darcy
    fDarcyAnalysis.PostProcessTimeStep(dim, step);
  }
  else if (type == 2) { // Transport
    fTransportAnalysis.PostProcessTimeStep(dim, step);
  }
  else {
    std::cout << "ERROR: Unknown problem type for post-processing: " << type << std::endl;
    DebugStop();
  }

  auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.;
  std::cout << "Total time post process = " << total_time << " seconds" << std::endl;
}

void TSFSFIAnalysis::Run(std::ostream &out) {
  auto GetPlotBaseName = [](const std::string &fileName) {
    const size_t dotPos = fileName.find('.');
    if (dotPos == std::string::npos) {
      return fileName;
    }
    return fileName.substr(0, dotPos);
  };

  auto WritePVD = [](const std::string &pvdFileName, const std::vector<std::pair<REAL, std::string>> &entries) {
    std::ofstream outFile(pvdFileName);
    outFile << "<?xml version=\"1.0\"?>\n";
    outFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    outFile << "  <Collection>\n";
    for (const auto &entry : entries) {
      outFile << "    <DataSet timestep=\"" << entry.first << "\" group=\"\" part=\"0\" file=\"" << entry.second << "\"/>\n";
    }
    outFile << "  </Collection>\n";
    outFile << "</VTKFile>\n";
  };

  auto WriteVTKSeries = [](const std::string &seriesFileName, const std::vector<std::pair<REAL, std::string>> &entries) {
    std::ofstream outFile(seriesFileName);
    outFile << "{\n";
    outFile << "  \"file-series-version\": \"1.0\",\n";
    outFile << "  \"files\": [\n";
    for (size_t i = 0; i < entries.size(); i++) {
      outFile << "    { \"name\": \"" << entries[i].second << "\", \"time\": " << entries[i].first << " }";
      if (i + 1 < entries.size()) {
        outFile << ",";
      }
      outFile << "\n";
    }
    outFile << "  ]\n";
    outFile << "}\n";
  };

  const std::string darcyPlotBaseName = GetPlotBaseName(fSimData->fTPostProcess.fFileNameDarcy);
  const std::string transportPlotBaseName = GetPlotBaseName(fSimData->fTPostProcess.fFileNameTransport);
  std::vector<std::pair<REAL, std::string>> darcyPVDEntries;
  darcyPVDEntries.reserve(fSimData->fTPostProcess.fVecReportingTimes.size() + 2);
  std::vector<std::pair<REAL, std::string>> transportPVDEntries;
  transportPVDEntries.reserve(fSimData->fTPostProcess.fVecReportingTimes.size() + 2);

  auto RegisterPostProcess = [&](const int type, const REAL currentTime, const int step) {
    if (type == 0 || type == 1) {
      darcyPVDEntries.push_back({currentTime, darcyPlotBaseName + "." + std::to_string(step) + ".vtk"});
      // WritePVD(darcyPlotBaseName + ".pvd", darcyPVDEntries);
      WriteVTKSeries(darcyPlotBaseName + ".vtk.series", darcyPVDEntries);
    }
    if (type == 0 || type == 2) {
      transportPVDEntries.push_back({currentTime, transportPlotBaseName + "." + std::to_string(step) + ".vtk"});
      // WritePVD(transportPlotBaseName + ".pvd", transportPVDEntries);
      WriteVTKSeries(transportPlotBaseName + ".vtk.series", transportPVDEntries);
    }
  };

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
  RegisterPostProcess(fSimData->fTPostProcess.fProblemTypeInit, time, 0);

  // Time stepping loop
  int pos = 0;
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
      RegisterPostProcess(fSimData->fTPostProcess.fProblemTypeInit, time, tstep);
    } else if (time >= nextPostProcessTime - 1.0e-8) {
      PostProcessTimeStep(fSimData->fTPostProcess.fProblemType, dim, tstep);
      RegisterPostProcess(fSimData->fTPostProcess.fProblemType, time, tstep);
      nextPostProcessTime = postProcessTimes[++pos];
    }
    // Check Mass Balance
    fTransportAnalysis.fAlgebraicTransport.CheckMassBalance(time, out);
  }
  out << "-------------------- Simulation Finished --------------------" << std::endl;
  SaveSolution();
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
      if (fSimData->fTNumerics.fAnalysisType == 2) fShouldSolveDarcy = false;
    }

    // Solve Transport Problem
    out << "Solving Transport Problem..." << std::endl;
    if (fSimData->fTNumerics.fAnalysisType != 1) fTransportAnalysis.RunTimeStep(out);

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

void TSFSFIAnalysis::SaveSolution() {
  std::ofstream darcy_out("darcy-final-solution.txt");
  auto& darcy_sol = fDarcyAnalysis.Mesh()->Solution();
  darcy_sol.Print("Darcy final solution", darcy_out);

  std::ofstream transport_out("transport-final-solution.txt");
  auto& transport_sol = fTransportAnalysis.Solution();
  transport_sol.Print("Transport final solution", transport_out);
}