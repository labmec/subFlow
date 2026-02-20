//
//  Created by Giovane Avancini on 02/09/25.
//

#include "TSFDarcyAnalysis.h"
#include "TPZVTKGenerator.h"
#include "TPZInterfaceEl.h"

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
  SetInitialBCValue(neumannMatids);
  ApplyEquationFilter(neumannMatids);
  int nreducedeq = fStructMatrix->NReducedEquations();
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  SetSolver(step);

  std::cout << "Number of Darcy equations: " << fCompMesh->NEquations() << std::endl;
  std::cout << "Number of Darcy elements: " << fCompMesh->NElements() << std::endl;
}

void TSFDarcyAnalysis::RunTimeStep(std::ostream &out) {
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
    UpdateDensityAndCoefficients();
    Assemble();
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
    TPZFMatrix<STATE> &dsol = Solution();
    corr_norm = Norm(dsol);
    sol += dsol;
    cmesh->LoadSolution(sol);
    cmesh->TransferMultiphysicsSolution();
  }

  if (!converged) {
    std::cout << "------Iterative method did not converge. res_norm: " << res_norm << " corr_norm: " << corr_norm << std::endl;
    fSolution = sol;
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

void TSFDarcyAnalysis::SetInitialBCValue(std::set<int> &neumannMatids) {
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

void TSFDarcyAnalysis::SetInitialSolution() {
  if (!fSimData->fTReservoirProperties.fP0Func) return;
  fCompMesh->LoadReferences();
  TPZFMatrix<STATE> &cmesh_sol = fCompMesh->Solution();
  TPZMultiphysicsCompMesh *mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(fCompMesh);
  TPZGeoMesh *gmesh = fCompMesh->Reference();
  const int dim = gmesh->Dimension();

  TPZCompMesh *pressure_cmesh = mp_cmesh->MeshVector()[1];
  TPZFMatrix<STATE> &pressure_cmesh_sol = pressure_cmesh->Solution();
  for (TPZCompEl *cel : pressure_cmesh->ElementVec()) { // Pressure mesh
    TPZGeoEl *gel = cel->Reference();
    if (gel->Dimension() != dim) continue;
    int ncorner = gel->NCornerNodes();
    TPZMultiphysicsElement *mp_el = dynamic_cast<TPZMultiphysicsElement *>(gel->Reference());
#ifdef PZDEBUG
    if (!mp_el) {
      DebugStop();
    }
#endif
    int cindex = mp_el->ElementVec()[0].Element()->NConnects(); // pressure connects come after flux connects
    for (int i = 0; i < ncorner; i++) {
      TPZManVector<REAL, 3> coord(3, 0.0);
      gel->NodePtr(i)->GetCoordinates(coord);
      REAL val = fSimData->fTReservoirProperties.fP0Func(coord);
      TPZConnect &cloc = cel->Connect(i);
      TPZConnect &c = mp_el->Connect(cindex + i);
      int64_t seqloc = cloc.SequenceNumber(); //
      int64_t seq = c.SequenceNumber();
      int64_t firstEqLoc = pressure_cmesh->Block().Position(seqloc);
      int64_t firstEq = fCompMesh->Block().Position(seq);
      int blockSize = fCompMesh->Block().Size(seq);
#ifdef PZDEBUG
      int blockSizeloc = pressure_cmesh->Block().Size(seqloc);
      if (blockSize != blockSizeloc) {
        DebugStop();
      }
#endif
      for (int64_t eqloc = firstEqLoc; eqloc < firstEqLoc + blockSizeloc; eqloc++) { // atomic pressure solution
        pressure_cmesh_sol.PutVal(eqloc, 0, val);
      }
      for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++) { // multiphysics pressure solution
        cmesh_sol.PutVal(eq, 0, val);
      }
    }
  }

  TPZCompMesh *avgpressure_cmesh = mp_cmesh->MeshVector()[3];
  TPZFMatrix<STATE> &avgpressure_cmesh_sol = avgpressure_cmesh->Solution();
  for (TPZCompEl *cel : avgpressure_cmesh->ElementVec()) { // Average pressure elements
    TPZGeoEl *gel = cel->Reference();
    if (gel->Dimension() != dim) continue;
    int nsides = gel->NSides();
    TPZManVector<REAL, 3> center_local(dim, 0.0);
    TPZManVector<REAL, 3> center_global(3, 0.0);
    gel->CenterPoint(nsides - 1, center_local);
    gel->X(center_local, center_global);
    REAL val = fSimData->fTReservoirProperties.fP0Func(center_global);
    TPZMultiphysicsElement *mp_el = dynamic_cast<TPZMultiphysicsElement *>(gel->Reference());
#ifdef PZDEBUG
    if (!mp_el) {
      DebugStop();
    }
#endif
#ifdef PZDEBUG
    if (cel->NConnects() != 1)
      DebugStop();
#endif
    int cindex = mp_el->NConnects() - 1; // avg pressure connect is the last one
    TPZConnect &cloc = cel->Connect(0);
    TPZConnect &c = mp_el->Connect(cindex);
    int64_t seqloc = cloc.SequenceNumber(); //
    int64_t seq = c.SequenceNumber();
    int64_t firstEqLoc = avgpressure_cmesh->Block().Position(seqloc);
    int64_t firstEq = fCompMesh->Block().Position(seq);
    int blockSize = fCompMesh->Block().Size(seq);
#ifdef PZDEBUG
    int blockSizeloc = avgpressure_cmesh->Block().Size(seqloc);
    if (blockSize != blockSizeloc) {
      DebugStop();
    }
#endif
    for (int64_t eqloc = firstEqLoc; eqloc < firstEqLoc + blockSizeloc; eqloc++) { // atomic pressure solution
      avgpressure_cmesh_sol.PutVal(eqloc, 0, val);
    }
    for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++) { // multiphysics pressure solution
      cmesh_sol.PutVal(eq, 0, val);
    }
  }

  // Transfering multiphysics solution to the analysis solution
  int cmesh_neq = fCompMesh->NEquations();
  TPZFMatrix<STATE> &sol = Solution();
  for (int i = 0; i < cmesh_neq; i++) {
    sol.PutVal(i, 0, cmesh_sol.GetVal(i, 0));
  }
}

void TSFDarcyAnalysis::ApplyEquationFilter(std::set<int> &neumannMatids) {
  fCompMesh->LoadReferences();
  std::set<int64_t> removeEquations;
  TPZGeoMesh *gmesh = fCompMesh->Reference();
  TPZFMatrix<STATE> sol = fCompMesh->Solution();
  for (auto el : gmesh->ElementVec()) {
    int elMatID = el->MaterialId();

    if (neumannMatids.find(elMatID) == neumannMatids.end()) continue;

    TPZCompEl *compEl = el->Reference();

    TPZInterfaceElement *interfaceEl = dynamic_cast<TPZInterfaceElement *>(compEl);
    if (interfaceEl) continue; // skip interface elements, it only belongs to Transport mesh

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

void TSFDarcyAnalysis::SetLastStateVariables() {
  TPZCompMesh *cmesh = Mesh();
#ifdef PZDEBUG
  TPZMultiphysicsCompMesh *mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
  if (!cmesh)
    DebugStop();
#endif

  int nels = cmesh->NElements();
  for (int iel = 0; iel < nels; iel++) {
    TPZCompEl *cel = cmesh->Element(iel);
    TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
    if (!condensed) continue;

    REAL sw = condensed->GetSw();
    condensed->SetSwLast(sw);

    TPZCompEl *compel = condensed->ReferenceCompEl();
    int dim = compel->Dimension();
    TPZVec<REAL> qsi(dim, 0.0);
    TPZVec<STATE> sol(dim, 0.0);
    int presureindex = 2;
    compel->Solution(qsi, presureindex, sol);
    REAL pressure = sol[0];
    condensed->SetPressureLastState(pressure);
  }
}

void TSFDarcyAnalysis::SetLastStateSaturation() {
  TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
  if (!cmesh)
    DebugStop();

  int nels = cmesh->NElements();
  for (int iel = 0; iel < nels; iel++) {
    TPZCompEl *cel = cmesh->Element(iel);
    TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
    if (!condensed) continue;

    REAL sw = condensed->GetSw();
    condensed->SetSwLast(sw);
  }
}

void TSFDarcyAnalysis::SetLastStatePressure() {
  TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
  if (!cmesh)
    DebugStop();

  int nels = cmesh->NElements();
  for (int iel = 0; iel < nels; iel++) {
    TPZCompEl *cel = cmesh->Element(iel);
    TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
    if (!condensed) continue;

    TPZCompEl *compel = condensed->ReferenceCompEl();
    int dim = compel->Dimension();
    TPZVec<REAL> qsi(dim, 0.0);
    TPZVec<STATE> sol(dim, 0.0);
    int presureindex = 2;
    compel->Solution(qsi, presureindex, sol);
    REAL pressure = sol[0];
    condensed->SetPressureLastState(pressure);
  }
}

void TSFDarcyAnalysis::SetTime(REAL time) {
  this->fTime = time;
  TSFMixedDarcy::fTime = time;
}

void TSFDarcyAnalysis::UpdateDensityAndCoefficients() {
  TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
  if (!cmesh)
    DebugStop();

  auto fWaterDensityF = fSimData->fTFluidProperties.fWaterDensityFunc;
  auto fGasDensityF = fSimData->fTFluidProperties.fGasDensityFunc;
  int64_t nels = cmesh->NElements();
  for (int iel = 0; iel < nels; iel++) {
    TPZCompEl *cel = cmesh->Element(iel);
    TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
    if (!condensed) continue;

    // Getting the avg pressure current solution
    TPZCompEl *compel = condensed->ReferenceCompEl();
    TPZGeoEl *gel = compel->Reference();
    int dim = compel->Dimension();
    TPZVec<REAL> qsi(dim, 0.0);
    TPZVec<STATE> sol(dim, 0.0);
    int pressureindex = 2;
    compel->Solution(qsi, pressureindex, sol);
    REAL pressure = sol[0];

    REAL vol = gel->Volume();
    REAL rhowref = fSimData->fTFluidProperties.fWaterDensityRef;
    REAL rhogref = fSimData->fTFluidProperties.fGasDensityRef;

    // Update the density based on the pressure if functions are provided
    std::tuple<REAL, REAL> rhoWvalderiv, rhoGvalderiv;
    if (fWaterDensityF && fGasDensityF) {
      rhoWvalderiv = fWaterDensityF(pressure);
      rhoGvalderiv = fGasDensityF(pressure);
#ifdef PZDEBUG
      if (std::get<0>(rhoWvalderiv) < 0.0 || std::get<0>(rhoGvalderiv) < 0.0)
        std::cout << "Negative density value detected! " << std::endl;
#endif
    } else {
      rhoWvalderiv = std::make_tuple(rhowref, 0.0);
      rhoGvalderiv = std::make_tuple(rhogref, 0.0);
    }

    REAL rhow = std::get<0>(rhoWvalderiv);
    REAL drhow = std::get<1>(rhoWvalderiv);
    REAL rhog = std::get<0>(rhoGvalderiv);
    REAL drhog = std::get<1>(rhoGvalderiv);
    REAL bw = rhowref / rhow;
    REAL bg = rhogref / rhog;
    REAL dbwinv = drhow / rhowref;
    REAL dbginv = drhog / rhogref;

    // With the new density we update the coefficients
    int krModel = fSimData->fTPetroPhysics.fKrModel;
    auto lambdaWfunc = fSimData->fTPetroPhysics.fLambdaw[krModel];
    auto lambdaGfunc = fSimData->fTPetroPhysics.fLambdag[krModel];
    auto lambdaTotalfunc = fSimData->fTPetroPhysics.fLambdaTotal[krModel];
    auto fwfunc = fSimData->fTPetroPhysics.fFw[krModel];
    auto fgfunc = fSimData->fTPetroPhysics.fFg[krModel];

    REAL sw = condensed->GetSw();
    auto fwfvalderiv = fwfunc(sw, bw, bg);
    auto fgvalderiv = fgfunc(sw, bw, bg);
    auto lambdaWvalderiv = lambdaWfunc(sw, bw);
    auto lambdaGvalderiv = lambdaGfunc(sw, bg);
    auto lambdaTotalvalderiv = lambdaTotalfunc(sw, bw, bg);

    // Update the mixed density
    REAL fw = std::get<0>(fwfvalderiv);
    REAL fg = std::get<0>(fgvalderiv);

    REAL mixedDensity = rhow * fw + rhog * fg;
    condensed->SetMixedDensity(mixedDensity);

    // Update the coefficients
    REAL lambda = std::get<0>(lambdaTotalvalderiv);
    condensed->SetLambda(lambda);

    if (fWaterDensityF && fGasDensityF) {
      int matid = compel->Material()->Id();
      REAL porosity = 0.0;
      for (auto &domain : fSimData->fTReservoirProperties.fPorosityAndPermeability) {
        if (std::get<0>(domain) == matid) {
          porosity = std::get<1>(domain).first;
          break;
        }
      }
      REAL dt = fSimData->fTNumerics.fDt;
      REAL sg = 1.0 - sw;
      REAL swlast = condensed->GetSwLast();
      REAL sglast = 1.0 - swlast;
      REAL compterm = -vol * (porosity / dt) * (sw * dbwinv + sg * dbginv);

      REAL pressurelast = condensed->GetPressureLastState();
      REAL rhoWlast = std::get<0>(fWaterDensityF(pressurelast));
      REAL rhoGlast = std::get<0>(fGasDensityF(pressurelast));
      REAL bwlast = rhowref / rhoWlast;
      REAL bglast = rhogref / rhoGlast;
      REAL termrhscurrent = (sw / bw) + (sg / bg);
      REAL termrhslast = (swlast / bwlast) + (sglast / bglast);
      REAL comptermrhs = vol * (porosity / dt) * (termrhscurrent - termrhslast);

      condensed->SetCompressibiilityTerm(compterm, comptermrhs);
    }
  }
}