//
// Created by Giovane Avancini on 01/09/2025
//

#include "TSFMixedDarcy.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

TSFMixedDarcy::TSFMixedDarcy() : TPZRegisterClassId(&TSFMixedDarcy::ClassId), TBase(), fIsAxisymmetric(false), fFourSpaces(false), fGravity(3, 1, 0.0) {}

TSFMixedDarcy::TSFMixedDarcy(int id, int dim) : TPZRegisterClassId(&TSFMixedDarcy::ClassId), TBase(id, dim), fIsAxisymmetric(false), fFourSpaces(false), fGravity(3, 1, 0.0) {}

TSFMixedDarcy::TSFMixedDarcy(const TSFMixedDarcy &copy) : TBase(copy) {
  *this = copy;
}

TSFMixedDarcy &TSFMixedDarcy::operator=(const TSFMixedDarcy &copy) {
  TBase::operator=(copy);
  fIsAxisymmetric = copy.fIsAxisymmetric;
  fFourSpaces = copy.fFourSpaces;
  fGravity = copy.fGravity;
  return *this;
}

void TSFMixedDarcy::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

  if (datavec.size() < 2) DebugStop();

  TPZFNMatrix<110, REAL> &phiU = datavec[0].fDeformedDirections;
  TPZFNMatrix<20, REAL> &phiP = datavec[1].phi;
  TPZFNMatrix<20, REAL> &divU = datavec[0].divphi;
  TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

  int nphiU, nphiP;
  nphiU = datavec[0].fDeformedDirections.Cols();
  nphiP = phiP.Rows();

  TPZVec<STATE> Usol = datavec[0].sol[0];
  TPZFMatrix<STATE> UsolMat(Usol.size(), 1);
  for (int i = 0; i < Usol.size(); ++i) {
    UsolMat(i, 0) = Usol[i];
  }
  REAL psol = datavec[1].sol[0][0];
  REAL divUsol = datavec[0].divsol[0][0];

  REAL axiFactor = 1.0;
  if (fIsAxisymmetric) // Axisymmetric: assuming radius is aligned with the x axis
  {
    REAL r = datavec[0].x[0];
    axiFactor = 1.0 / (2.0 * M_PI * r);
  }

  // Tangent matrix
  REAL K = GetPermeability(datavec[0].x);
  REAL factor = weight * axiFactor / K;
  ek.AddContribution(0, 0, phiU, 1, phiU, 0, factor);      // A
  ek.AddContribution(nphiU, 0, phiP, 0, divU, 1, -weight); // B^T
  ek.AddContribution(0, nphiU, divU, 0, phiP, 1, -weight); // B

  // Residual vector constitutive equation (negative)
  ef.AddContribution(0, 0, phiU, 1, UsolMat, 0, -factor);
  factor = psol * weight;
  ef.AddContribution(0, 0, divU, 0, Aux, 1, factor);
  ef.AddContribution(0, 0, phiU, 1, fGravity, 0, weight);

  // Residual vector conservation equation (negative)
  factor = divUsol * weight;
  ef.AddContribution(nphiU, 0, phiP, 0, Aux, 0, factor);

  if (fFourSpaces && datavec.size() < 4) DebugStop();

  if (fFourSpaces) ContributeFourSpaces(datavec, weight, ek, ef);
}

void TSFMixedDarcy::ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

  int qb = 0;
  int pb = 1;
  int numactive = 0;
  for (auto &it : datavec)
    if (it.fActiveApproxSpace) numactive++;
  if (numactive % 2 != 0) DebugStop();

  int nAverage = (numactive - 2) / 2;
  TPZFNMatrix<20, REAL> &phiP = datavec[1].phi;
  int nPhiU = datavec[0].fVecShapeIndex.NElements();
  int nPhiP = phiP.Rows();
  if (nPhiU + nPhiP + nAverage * 2 != ek.Rows()) DebugStop();

  for (int iavg = 0; iavg < nAverage; iavg++) {
    int g_avgb = 2 + 2 * iavg;
    int p_avgb = 3 + 2 * iavg;

    int nphi_gb = datavec[g_avgb].phi.Rows();
    int nphi_pb = datavec[p_avgb].phi.Rows();

    if (nphi_gb != 1 || nphi_pb != 1) DebugStop();

    STATE p = datavec[pb].sol[0][0];
    STATE g_avg = datavec[g_avgb].sol[0][0];
    STATE p_avg = datavec[p_avgb].sol[0][0];

    for (int ip = 0; ip < nPhiP; ip++) {
      ef(nPhiU + ip, 0) += weight * g_avg * phiP(ip, 0);
      ek(nPhiU + ip, nPhiU + nPhiP + 2 * iavg) += weight * phiP(ip, 0);
      ek(nPhiU + nPhiP + 2 * iavg, nPhiU + ip) += weight * phiP(ip, 0);
    }

    ef(nPhiU + nPhiP + 1 + 2 * iavg, 0) += -weight * g_avg;
    ek(nPhiU + nPhiP + 1 + 2 * iavg, nPhiU + nPhiP + 2 * iavg) += -weight;

    ef(nPhiU + nPhiP + 2 * iavg, 0) += weight * (p - p_avg);
    ek(nPhiU + nPhiP + 2 * iavg, nPhiU + nPhiP + 1 + 2 * iavg) += -weight;
  }
}

void TSFMixedDarcy::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

  int dim = Dimension();

  REAL bigNumber = TPZMaterial::fBigNumber * 1.e-2;

  TPZFNMatrix<20, REAL> &phiU = datavec[0].phi;
  int nPhiU = phiU.Rows();
  TPZManVector<STATE, 3> Usol = datavec[0].sol[0];

  REAL axiFactor = 1.0;
  if (fIsAxisymmetric) // Axisymmetric: assuming radius is aligned with the x axis
  {
    REAL r = datavec[0].x[0];
    axiFactor = 1.0 / (2.0 * M_PI * r);
    Usol[0] *= axiFactor;
  }

  REAL v2 = bc.Val2()[0];
  REAL v1 = bc.Val1()(0, 0);
  REAL u_D = 0;
  REAL normflux = 0.;

  if (bc.HasForcingFunctionBC()) {
    TPZManVector<STATE> res(3);
    TPZFNMatrix<9, STATE> gradu(3, 1);
    bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

    const STATE perm = GetPermeability(datavec[0].x);

    for (int i = 0; i < 3; i++) {
      normflux += datavec[0].normal[i] * perm * gradu(i, 0);
    }

    if (bc.Type() == 0 || bc.Type() == 4) {
      v2 = res[0];
      u_D = res[0];
      normflux *= (-1.);
    } else if (bc.Type() == 1 || bc.Type() == 2) {
      v2 = -normflux;
      if (bc.Type() == 2) {
        v2 = -res[0] + v2 / v1;
      }
    } else if (bc.Type() == 5) {
      v2 = res[0];
    } else {
      DebugStop();
    }
  } else {
    v2 = bc.Val2()[0];
  }

  switch (bc.Type()) {
  case 0: // Dirichlet condition
    for (int i = 0; i < nPhiU; i++) {
      // the contribution of the Dirichlet boundary condition appears in the flow equation
      ef(i, 0) += (-1.) * v2 * phiU(i, 0) * weight;
    }
    break;

  case 1: // Neumann condition
    // for (int i = 0; i < nPhiU; i++) {
    //   ef(i, 0) += bigNumber * (v2 - Usol[0]) * phiU(i, 0) * weight;
    //   for (int j = 0; j < nPhiU; j++) {

    //     ek(i, j) += bigNumber * phiU(i, 0) * phiU(j, 0) * weight;
    //   }
    // }
    break;
  }
}

void TSFMixedDarcy::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {

  int nref = datavec.size();
  for (int i = 0; i < nref; i++) {
    datavec[i].SetAllRequirements(false);
    datavec[i].fNeedsSol = true;
  }
}

void TSFMixedDarcy::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {

  int nref = datavec.size();
  for (int iref = 0; iref < nref; iref++) {
    datavec[iref].SetAllRequirements(false);
    datavec[iref].fNeedsSol = true;
    datavec[iref].fNeedsNormal = true;
  }
}

void TSFMixedDarcy::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {

  Solout.Resize(this->NSolutionVariables(var));
  TPZManVector<STATE, 3> p, u;

  u = datavec[0].sol[0];
  p = datavec[1].sol[0];
  REAL div_u = datavec[0].divsol[0][0];

  REAL axiFactor = 1.0;
  if (fIsAxisymmetric) // Axisymmetric: assuming radius is aligned with the x axis
  {
    REAL r = datavec[0].x[0];
    axiFactor = 1.0 / (2.0 * M_PI * r);
    u[0] *= axiFactor;
    u[1] *= axiFactor;
  }

  if (var == 1) {
    for (int i = 0; i < 3; i++) {
      Solout[i] = u[i];
    }
    return;
  }

  if (var == 2) {
    Solout[0] = p[0];
    return;
  }

  if (var == 3) {
    Solout[0] = div_u;
    return;
  }

  if (fFourSpaces) {

    int g_avgb = 2;
    int p_avgb = 3;

    if (var == 5) {
      Solout[0] = datavec[g_avgb].sol[0][0];
      return;
    }
    if (var == 6) {
      Solout[0] = datavec[p_avgb].sol[0][0];
      return;
    }
  }

  DebugStop();
}