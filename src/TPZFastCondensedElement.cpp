//
//  TPZFastCondensedElement.cpp
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//

#include "TPZFastCondensedElement.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcondensedcompel"));
#endif

bool TPZFastCondensedElement::fSkipLoadSolution = false;

TPZFastCondensedElement::TPZFastCondensedElement(TPZCompEl *ref, bool keepmatrix) : TPZCondensedCompElT(ref, keepmatrix) {
  IdentifyConnectandEquations();
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZFastCondensedElement::CalcStiff(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef) {
  if (this->fMatrixComputed == false) {
    TPZCondensedCompElT::CalcStiff(fEK, fEF);
    ComputeBodyforceRefValues();
    ComputeConstantPressureValues();
    fMatrixComputed = true;
  }
  ek = fEK;
  ef = fEF;
  int nrows = ek.fMat.Rows();
  int ncols = ek.fMat.Cols();
  REAL Glambda = fMixedDensity;

  ek.fMat *= (1. / fLambda);
  for (int icol = 0; icol < ncols; icol++) {
    ek.fMat(nrows - 1, icol) *= fLambda;
  }
  for (int irow = 0; irow < nrows; irow++) {
    ek.fMat(irow, ncols - 1) *= fLambda;
  }

  TPZFNMatrix<30, STATE> solvec(fEK.fMat.Rows(), 1, 0.);
  GetSolutionVector(solvec);

  // When an initial solution is given, the residual contains not only the gravitational forces, but also K*sol0.
  // We cannot multiply it by Glambda.
  // The ideal fix would be adding a flag to the material refered to the use or not of FastCondensedCompel. If so, we only compute the body forces in Contribute(),
  // and the K*sol0 term is computed here during the MultAdd call. In this case, the MultAdd call should be called always, not only when fMatrixComputed is true.
  // PLEASE FIX ME
  ef.fMat *= 1.0 * Glambda;
  ef.fMat(nrows - 1) += fCompressibiilityRhsTerm; // should use the += operator since ef already has the flux divergence computed in the Contribute method

  /** @brief Computes z = alpha * opt(this)*x + beta * y */
  /** @note z and x cannot overlap in memory */
  //    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
  //                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
  STATE alpha = -1.;

  // If an initial solution is given, this should not be called at the first time step,
  // since K * solvec is already computed in ef during Contribute()
  ek.fMat.MultAdd(solvec, ef.fMat, ef.fMat, alpha, 1);
  ek.fMat(nrows - 1, ncols - 1) = fCompressibilityMatrixTerm; // this term is initially null as only the incompressible part is computed in the Contribute method
}

// extract the solution vector of the condensed element
void TPZFastCondensedElement::GetSolutionVector(TPZFMatrix<STATE> &solvec) {
  int nc = fEK.fConnect.size();
  TPZCompMesh *cmesh = Mesh();
  TPZFMatrix<STATE> &meshsol = cmesh->Solution();
  int64_t vecsize = fEK.fMat.Rows();
  int count = 0;
  for (int ic = 0; ic < nc; ic++) {
    int64_t cindex = fEK.fConnect[ic];
    TPZConnect &c = cmesh->ConnectVec()[cindex];
    int64_t seqnum = c.SequenceNumber();
    int blsize = c.NShape() * c.NState();
    for (int dof = 0; dof < blsize; dof++) {
      int ind = cmesh->Block().Index(seqnum, dof);
      solvec(count + dof, 0) = meshsol(ind, 0);
    }
    count += blsize;
  }
  if (count != vecsize) DebugStop();
}

/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZFastCondensedElement::CalcResidual(TPZElementMatrixT<STATE> &ef) {
  TPZElementMatrixT<STATE> ek;
  CalcStiff(ek, ef);
}

void TPZFastCondensedElement::SetLambda(REAL lambda) {
  fLambda = lambda;
}
REAL TPZFastCondensedElement::GetLambda() {
  return fLambda;
}
void TPZFastCondensedElement::SetSw(REAL sw) {
  fSw = sw;
}
REAL TPZFastCondensedElement::GetSw() {
  return fSw;
}
void TPZFastCondensedElement::SetSwLast(REAL swlast) {
  fSwLast = swlast;
}
REAL TPZFastCondensedElement::GetSwLast() {
  return fSwLast;
}
void TPZFastCondensedElement::SetPressureLastState(REAL pressureLastState) {
  fPressureLastState = pressureLastState;
}
REAL TPZFastCondensedElement::GetPressureLastState() {
  return fPressureLastState;
}
void TPZFastCondensedElement::SetMixedDensity(REAL mdensity) {
  fMixedDensity = mdensity;
}
void TPZFastCondensedElement::SetCompressibiilityTerm(REAL matrix, REAL rhs) {
  fCompressibilityMatrixTerm = matrix;
  fCompressibiilityRhsTerm = rhs;
}

REAL TPZFastCondensedElement::GetMixedDensity() {
  return fMixedDensity;
}
void TPZFastCondensedElement::SetPermTensorAndInv(TPZFNMatrix<9, REAL> &PermeabilityTensor, TPZFNMatrix<9, REAL> &InvPerm) {
  fPermeabilityTensor = PermeabilityTensor;
  fInvPerm = InvPerm;
}
TPZFMatrix<REAL> &TPZFastCondensedElement::GetPermTensor() {
  return fPermeabilityTensor;
}

/**
 * @brief Calculates the solution - sol - for the variable var
 * at point qsi, where qsi is expressed in terms of the
 * master element coordinates
 * @param qsi master element coordinate
 * @param var variable name
 * @param sol vetor for the solution
 */
void TPZFastCondensedElement::Solution(TPZVec<REAL> &qsi, int var, TPZVec<STATE> &sol) {
  switch (var) {
  case 7:
    sol[0] = fPermeabilityTensor(0, 0);
    break;
  case 8:
    sol[0] = fPermeabilityTensor(1, 1);
    break;
  case 9:
    sol[0] = fPermeabilityTensor(2, 2);
    break;
  case 10:
    sol[0] = fLambda;
    break;
  default:
    TPZCondensedCompElT::Solution(qsi, var, sol);
    break;
  }
}

/** @brief Loads the solution within the internal data structure of the element */
/**
 * Is used to initialize the solution of connect objects with dependency \n
 * Is also used to load the solution within SuperElements
 */
void TPZFastCondensedElement::LoadSolution() {
  if (fSkipLoadSolution) {
    TPZCompEl::LoadSolution();
  } else {
    ComputeInternalCoefficients();
  }
}

// global indices of the pressure equations
void TPZFastCondensedElement::PressureEquations(TPZVec<int64_t> &eqs) {
  int nconnects = fPressureConnects.size();
  int numpressure_equations = 0;
  TPZCompMesh *cmesh = Mesh();
  for (int ic = 0; ic < nconnects; ic++) {
    int64_t pressconnectindex = fPressureConnects[ic];
    TPZConnect &c = cmesh->ConnectVec()[pressconnectindex];
    int neq = c.NShape() * c.NState();
    numpressure_equations += neq;
  }
  TPZBlock &block = cmesh->Block();
  eqs.Resize(numpressure_equations, 0);
  int count = 0;
  for (int ic = 0; ic < nconnects; ic++) {
    int64_t pressconnectindex = fPressureConnects[ic];
    TPZConnect &c = cmesh->ConnectVec()[pressconnectindex];
    int neq = c.NShape() * c.NState();
    int64_t seqnum = c.SequenceNumber();
    int64_t firsteq = block.Position(seqnum);
    for (int i = 0; i < neq; i++) {
      eqs[count] = firsteq + i;
      count++;
    }
  }
}

// global indices of the pressure equations
void TPZFastCondensedElement::InternalFluxEquations(TPZVec<int64_t> &eqs) {
  int nconnects = fFluxConnects.size();
  int numflux_equations = 0;
  TPZCompMesh *cmesh = Mesh();
  for (int ic = 0; ic < nconnects; ic++) {
    int64_t fluxconnectindex = fFluxConnects[ic];
    TPZConnect &c = cmesh->ConnectVec()[fluxconnectindex];
    int neq = c.NShape() * c.NState();
    numflux_equations += neq;
  }
  TPZBlock &block = cmesh->Block();
  eqs.Resize(numflux_equations, 0);
  int count = 0;
  for (int ic = 0; ic < nconnects; ic++) {
    int64_t fluxconnectindex = fFluxConnects[ic];
    TPZConnect &c = cmesh->ConnectVec()[fluxconnectindex];
    int neq = c.NShape() * c.NState();
    int64_t seqnum = c.SequenceNumber();
    int64_t firsteq = block.Position(seqnum);
    for (int i = 0; i < neq; i++) {
      eqs[count] = firsteq + i;
      count++;
    }
  }
}

// global index of the average pressure equation
int64_t TPZFastCondensedElement::AveragePressureEquation() {
  if (fAveragePressureConnect == -1) DebugStop();

  TPZConnect &c = Mesh()->ConnectVec()[fAveragePressureConnect];
  int64_t seq_num = c.SequenceNumber();
#ifdef PZDEBUG
  if (Mesh()->Block().Size(seq_num) != 1) DebugStop();
#endif
  int64_t globeq = Mesh()->Block().Position(seq_num);
  return globeq;
}

// global indices of the boundary flux (external) equations
void TPZFastCondensedElement::BoundaryFluxEquations(TPZVec<int64_t> &eqs) {

  int nconnects = NConnects();
  int numflux_equations = 0;
  TPZCompMesh *cmesh = Mesh();
  // the last connect is the pressure connect
  for (int ic = 0; ic < nconnects; ic++) {
    int64_t cindex = ConnectIndex(ic);
    if (cindex == fAveragePressureConnect) {
      continue;
    }
    TPZConnect &c = Connect(ic);
    int neq = c.NShape() * c.NState();
    numflux_equations += neq;
    //        if(c.HasDependency())
    //        {
    //            int64_t cindex = c.FirstDepend()->fDepConnectIndex;
    //            TPZConnect &cdep = cmesh->ConnectVec()[cindex];
    //            neq = cdep.NShape()*cdep.NState();
    //            numflux_equations += neq;
    //        }
  }
  TPZBlock &block = cmesh->Block();
  eqs.Resize(numflux_equations, 0);
  int count = 0;
  for (int ic = 0; ic < nconnects; ic++) {
    int64_t cindex = ConnectIndex(ic);
    if (cindex == fAveragePressureConnect) {
      continue;
    }
    TPZConnect &c = Connect(ic);
    int neq = c.NShape() * c.NState();
    int64_t seqnum = c.SequenceNumber();
    int64_t firsteq = block.Position(seqnum);
    for (int i = 0; i < neq; i++) {
      eqs[count] = firsteq + i;
      count++;
    }
    //        if(c.HasDependency())
    //        {
    //            int64_t cindex = c.FirstDepend()->fDepConnectIndex;
    //            TPZConnect &cdep = cmesh->ConnectVec()[cindex];
    //            seqnum = cdep.SequenceNumber();
    //            firsteq = block.Position(seqnum);
    //            for(int i = 0; i<neq; i++)
    //            {
    //                eqs[count] = firsteq+i;
    //                count++;
    //            }
    //        }
  }
}

// adjust the multiplying coeficients of the pressure equations
void TPZFastCondensedElement::AdjustPressureCoefficients() {
  TPZManVector<int64_t, 20> fluxeqs, pressureqs;
  PressureEquations(pressureqs);
  BoundaryFluxEquations(fluxeqs);
  int64_t averagepressureq = AveragePressureEquation();
  TPZFMatrix<STATE> &solution = *(Mesh()->Block().Matrix<STATE>());
  int npres = pressureqs.size();
  int nflux = fluxeqs.size();
  TPZManVector<STATE> gravity_pressure(npres, 0.);
  TPZManVector<STATE> average_pressure(npres, 0.);
  TPZManVector<STATE> flux_pressure(npres, 0.);
  TPZManVector<STATE> boundary_fluxes(nflux, 0.);
  STATE average = solution(averagepressureq, 0);
  // store coeficients of the f
  for (int ifl = 0; ifl < nflux; ifl++) {
    boundary_fluxes[ifl] = solution(fluxeqs[ifl], 0);
  }
  // build the gravity pressure coefs
  solution(averagepressureq, 0) = 0.;
  // zero de boundary fluxes
  for (int ifl = 0; ifl < nflux; ifl++)
    solution(fluxeqs[ifl], 0) = 0.;
  TPZCondensedCompElT::LoadSolution();
  for (int ipr = 0; ipr < npres; ipr++)
    gravity_pressure[ipr] = solution(pressureqs[ipr], 0);

  // build the pressure for constant pressure
  // zero the boundary fluxes
  for (int ifl = 0; ifl < nflux; ifl++)
    solution(fluxeqs[ifl], 0) = 0.;
  solution(averagepressureq, 0) = average;
  TPZCondensedCompElT::LoadSolution();
  for (int ipr = 0; ipr < npres; ipr++)
    average_pressure[ipr] = solution(pressureqs[ipr], 0) - gravity_pressure[ipr];

  // build the pressure solution due to boundary fluxes
  // restore the boundary fluxes
  for (int ifl = 0; ifl < nflux; ifl++)
    solution(fluxeqs[ifl], 0) = boundary_fluxes[ifl];
  solution(averagepressureq, 0) = 0.;
  TPZCondensedCompElT::LoadSolution();
  for (int ipr = 0; ipr < npres; ipr++)
    flux_pressure[ipr] = solution(pressureqs[ipr], 0) - gravity_pressure[ipr];
  // compose the final pressure coeficients
  solution(averagepressureq, 0) = average;
  for (int ipr = 0; ipr < npres; ipr++) {
    solution(pressureqs[ipr], 0) = average_pressure[ipr] -
                                   fMixedDensity * gravity_pressure[ipr] +
                                   flux_pressure[ipr] / fLambda;
  }
}

/// compute internal coeficients as a function of the average pressure and boundary fluxes
void TPZFastCondensedElement::ComputeInternalCoefficients() {
  TPZManVector<int64_t, 20> fluxeqs, pressureqs;
  PressureEquations(pressureqs);
  InternalFluxEquations(fluxeqs);
  int64_t avpreseq = AveragePressureEquation();
  TPZCompMesh *cmesh = Mesh();
  TPZFMatrix<STATE> &sol = cmesh->Solution();
  STATE avpressure = sol(avpreseq);
  int npressure = pressureqs.size();
  int nfluxes = fluxeqs.size();
  // value of the pressure due to fluid flow
  TPZManVector<STATE, 100> PressureValues(npressure);
  // value of the fluxes due to fluid flow
  TPZManVector<STATE, 100> FluxValues(nfluxes);

  for (int i = 0; i < npressure; i++) {
    PressureValues[i] = fBodyForcePressureRef[i] * fMixedDensity + fConstantUnitPressure[i] * avpressure;
  }
  for (int i = 0; i < nfluxes; i++) {
    FluxValues[i] = fBodyForceFluxRef[i] * fMixedDensity / fLambda;
  }

  sol(avpreseq) = 0.;
  TPZManVector<STATE, 100> FlowFluxValues(nfluxes), FlowPressureValues(npressure);
  TPZCondensedCompElT::LoadSolution();
  for (int i = 0; i < npressure; i++) {
    FlowPressureValues[i] = sol(pressureqs[i]) / fLambda;
    sol(pressureqs[i]) = PressureValues[i] + FlowPressureValues[i];
  }
  for (int i = 0; i < nfluxes; i++) {
    FlowFluxValues[i] = sol(fluxeqs[i]);
    sol(fluxeqs[i]) = FluxValues[i] + FlowFluxValues[i];
  }
  sol(avpreseq) = avpressure;
}

static void GatherConnects(TPZCompEl *cel, std::set<std::pair<int64_t, int>> &connectset);

static void GatherConnects(TPZMultiphysicsElement *cel, std::set<std::pair<int64_t, int>> &connectset) {
  TPZManVector<std::pair<int64_t, int>> connectmesh;
  cel->GetConnectMeshPairs(connectmesh);
  int nc = connectmesh.size();
  connectset.insert(&connectmesh[0], &connectmesh[0] + nc);
}

static void GatherConnects(TPZElementGroup *elgr, std::set<std::pair<int64_t, int>> &connectset) {
  const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
  int nel = elvec.size();
  for (int el = 0; el < nel; el++) {
    TPZCompEl *cel = elvec[el];
    GatherConnects(cel, connectset);
  }
}

static void GatherConnects(TPZCompEl *cel, std::set<std::pair<int64_t, int>> &connectset) {
  if (!cel) DebugStop();
  TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
  TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
  TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
  if (cond) {
    GatherConnects(cond->ReferenceCompEl(), connectset);
    return;
  }
  if (elgr) {
    GatherConnects(elgr, connectset);
    return;
  }
  if (mphys) {
    GatherConnects(mphys, connectset);
    return;
  }
}

// Identify the connects and associated equations
void TPZFastCondensedElement::IdentifyConnectandEquations() {
  std::set<std::pair<int64_t, int>> connects;
  GatherConnects(fReferenceCompEl, connects);
  std::set<int64_t> externalconnects;
  int nconnects = NConnects();
  externalconnects.insert(&fActiveConnectIndexes[0], &fActiveConnectIndexes[0] + nconnects);
  std::map<int, int> nconnects_bymesh;

  int pressmesh = 1;
  int fluxmesh = 0;
  // int distfluxmesh = 2;
  //  int avpressmesh = 3;
  fAveragePressureConnect = -1;

  for (auto it : connects) {
    if (externalconnects.find(it.first) == externalconnects.end()) {
      nconnects_bymesh[it.second]++;
    }
    if (it.second == 3) {
      if (fAveragePressureConnect != -1) DebugStop();
      fAveragePressureConnect = it.first;
    }
  }

  if (fAveragePressureConnect == -1) DebugStop();

  fPressureConnects.resize(nconnects_bymesh[pressmesh]);
  fFluxConnects.resize(nconnects_bymesh[fluxmesh]);

  int iflux_connect = 0;
  int ipres_connect = 0;

  for (auto it : connects) {
    if (externalconnects.find(it.first) == externalconnects.end()) {
      int imesh = it.second;
      if (imesh == pressmesh) {
        fPressureConnects[ipres_connect++] = it.first;
      }
      if (imesh == fluxmesh) {
        fFluxConnects[iflux_connect++] = it.first;
      }
    }
  }

  if (ipres_connect != fPressureConnects.size()) DebugStop();
  if (iflux_connect != fFluxConnects.size()) DebugStop();
}

static void FindCondensed(TPZCompEl *cel, TPZStack<TPZCondensedCompElT<REAL> *> &condensedelements) {
  TPZCondensedCompElT<REAL> *cond = dynamic_cast<TPZCondensedCompElT<REAL> *>(cel);
  TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
  if (cond) {
    condensedelements.Push(cond);
    FindCondensed(cond->ReferenceCompEl(), condensedelements);
    return;
  }
  if (elgr) {
    const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
    int nel = elvec.size();
    for (int el = 0; el < nel; el++) {
      TPZCompEl *loccel = elvec[el];
      FindCondensed(loccel, condensedelements);
    }
  }
}

// Identify the condensed elements in this structure
void TPZFastCondensedElement::FindCondensed(TPZStack<TPZCondensedCompElT *> &condensedelements) {
  ::FindCondensed(this, condensedelements);
}

/// compute the body force reference values
void TPZFastCondensedElement::ComputeBodyforceRefValues() {
  int64_t avpress_eq = AveragePressureEquation();
  TPZManVector<int64_t, 25> bound_flux_eq;
  BoundaryFluxEquations(bound_flux_eq);
  TPZFMatrix<STATE> &sol = Mesh()->Solution();
  // copy the average pressure and set the solution value to zero
  STATE avpress = sol(avpress_eq, 0);
  sol(avpress_eq, 0) = 0.;
  TPZVec<STATE> bound_fluxes(bound_flux_eq.size(), 0.);
  // copy the values of the boundary fluxes and set the fluxes to zero
  for (int eq = 0; eq < bound_flux_eq.size(); eq++) {
    bound_fluxes[eq] = sol(bound_flux_eq[eq], 0);
    sol(bound_flux_eq[eq], 0) = 0.;
  }
  // after this call the internal solution will be due to the body forces
  TPZCondensedCompElT::LoadSolution();

  TPZVec<int64_t> pressure_eqs;
  TPZVec<int64_t> flux_eqs;
  PressureEquations(pressure_eqs);
  InternalFluxEquations(flux_eqs);

  fBodyForceFluxRef.resize(flux_eqs.size());
  fBodyForcePressureRef.resize(pressure_eqs.size());

  for (int eq = 0; eq < pressure_eqs.size(); eq++) {
    fBodyForcePressureRef[eq] = sol(pressure_eqs[eq], 0);
  }
  for (int eq = 0; eq < flux_eqs.size(); eq++) {
    fBodyForceFluxRef[eq] = sol(flux_eqs[eq], 0);
  }

  sol(avpress_eq, 0) = avpress;
  // copy the values of the boundary fluxes and set the fluxes to zero
  for (int eq = 0; eq < bound_flux_eq.size(); eq++) {
    sol(bound_flux_eq[eq], 0) = bound_fluxes[eq];
  }
}

/// compute pressure equation values with respect to a constant pressure
/// this will zero the body forces of the condensed elements
void TPZFastCondensedElement::ComputeConstantPressureValues() {
  TPZStack<TPZCondensedCompElT *> condensed;
  FindCondensed(condensed);
  for (int el = 0; el < condensed.size(); el++) {
    condensed[el]->Matrix().F0().Zero();
  }
  // avpress_eq is the equation of the constant pressure
  int64_t avpress_eq = AveragePressureEquation();
  TPZManVector<int64_t, 25> bound_flux_eq;
  BoundaryFluxEquations(bound_flux_eq);
  TPZFMatrix<STATE> &sol = Mesh()->Solution();
  // copy the average pressure and set the solution value to zero
  STATE avpress = sol(avpress_eq, 0);
  // average pressure is now equal one
  sol(avpress_eq, 0) = 1.;
  TPZManVector<STATE, 25> bound_fluxes(bound_flux_eq.size(), 0.);
  // copy the values of the boundary fluxes and set the fluxes to zero
  for (int eq = 0; eq < bound_flux_eq.size(); eq++) {
    bound_fluxes[eq] = sol(bound_flux_eq[eq], 0);
    sol(bound_flux_eq[eq], 0) = 0.;
  }
  // after this call the internal solution will be constant
  TPZCondensedCompElT::LoadSolution();

  TPZManVector<int64_t, 40> pressure_eqs;
  TPZManVector<int64_t, 40> flux_eqs;
  PressureEquations(pressure_eqs);
  InternalFluxEquations(flux_eqs);
  fConstantUnitPressure.resize(pressure_eqs.size());
  for (int eq = 0; eq < pressure_eqs.size(); eq++) {
    fConstantUnitPressure[eq] = sol(pressure_eqs[eq], 0);
  }
#ifdef PZDEBUG
  bool allok = true;
  TPZManVector<STATE, 40> fluxvals(flux_eqs.size());
  for (int eq = 0; eq < flux_eqs.size(); eq++) {
    fluxvals[eq] = sol(flux_eqs[eq]);
    if (abs(fluxvals[eq]) > 1.e-10) allok = false;
  }
  if (!allok) DebugStop();
  // NOTE: This can happen if BCs are now well set
#endif
  sol(avpress_eq, 0) = avpress;
  // copy the values of the boundary fluxes and set the fluxes to zero
  for (int eq = 0; eq < bound_flux_eq.size(); eq++) {
    sol(bound_flux_eq[eq], 0) = bound_fluxes[eq];
  }
}

void TPZFastCondensedElement::SetConnectIndex(int inode, int64_t index) {
  TPZCompEl *candidate = this->ReferenceCompEl();
  TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(candidate);
  this->SetMultiphysics(mphys);
  mphys->SetConnectIndex(inode, index);
}
