//
//  Created by Giovane Avancini on 10/01/26.
//  Copied from iMRS and is evolving to fit subFlow needs
//  PLEASE ADAPT THIS FILE
//

#include "TSFAlgebraicTransport.h"

/// Default constructor
TSFAlgebraicTransport::TSFAlgebraicTransport() {}

/// Copy constructor
TSFAlgebraicTransport::TSFAlgebraicTransport(const TSFAlgebraicTransport &other) {
  fNFluxCoefficients = other.fNFluxCoefficients;
  fNVolumesTransport = other.fNVolumesTransport;
  fCellsData = other.fCellsData;
  fInterfaceData = other.fInterfaceData;
  fNFluxCoefficients = other.fNFluxCoefficients;
  fboundaryCMatVal = other.fboundaryCMatVal;
}

/// Assignment operator
const TSFAlgebraicTransport &TSFAlgebraicTransport::operator=(const TSFAlgebraicTransport &other) {
  fNFluxCoefficients = other.fNFluxCoefficients;
  fNVolumesTransport = other.fNVolumesTransport;
  fCellsData = other.fCellsData;
  fInterfaceData = other.fInterfaceData;
  fboundaryCMatVal = other.fboundaryCMatVal;
  return *this;
}

/// Default destructor
TSFAlgebraicTransport::~TSFAlgebraicTransport() {}

void TSFAlgebraicTransport::TCellData::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
}

// TInterfaceData Methods
void TSFAlgebraicTransport::TInterfaceData::Print(std::ostream &out) {

  int ninterfaces = this->fFluxSign.size();
  out << "Material_ID: " << this->fMatid << std::endl;
  for (int iinter = 0; iinter < ninterfaces; iinter++) {
    out << "Left_Index: " << this->fLeftRightVolIndex[iinter].first << std::endl;
    out << "Right_Index: " << this->fLeftRightVolIndex[iinter].second << std::endl;
    out << "fCoefficientsFlux :";
    for (int icoe = 0; icoe < fCoefficientsFlux.size(); icoe++) {
      out << fCoefficientsFlux[icoe][iinter] << " ";
    }
    out << std::endl;
    out << "IntegralFluxFunctions: " << fIntegralFluxFunctions[iinter] << std::endl;
    out << "IntegralFlux: " << fIntegralFlux[iinter] << std::endl;
    out << "fFluxSign: " << fFluxSign[iinter] << std::endl;
    out << "fNormalDirection: ";
    out << std::get<0>(fNormalFaceDirection[iinter]) << " ";
    out << std::get<1>(fNormalFaceDirection[iinter]) << " ";
    out << std::get<2>(fNormalFaceDirection[iinter]) << " ";
    out << std::endl;
  }
}

// TCellData Methods
void TSFAlgebraicTransport::TCellData::UpdateSaturations(TPZFMatrix<STATE> &sol) {
  int ncells = fVolume.size();
  for (int icell = 0; icell < ncells; icell++) {
    int eq_number = fEqNumber[icell];
    fSaturation[icell] = sol(eq_number);
  }
}

void TSFAlgebraicTransport::TCellData::UpdateSaturationsLastState(TPZFMatrix<STATE> &sol) {
  int ncells = fVolume.size();
  for (int icell = 0; icell < ncells; icell++) {
    int eq_number = fEqNumber[icell];
    fSaturationLastState[icell] = sol(eq_number);
  }
}

void TSFAlgebraicTransport::TCellData::UpdateSaturationsTo(TPZFMatrix<STATE> &sol) {
  int ncells = fVolume.size();
  for (int icell = 0; icell < ncells; icell++) {
    int eq_number = fEqNumber[icell];
    sol(eq_number) = fSaturation[icell];
  }
}

void TSFAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambda(int krModel) {
  auto lambdaWfunc = fSimData->fTPetroPhysics.fLambdaw[krModel];
  auto lambdaGfunc = fSimData->fTPetroPhysics.fLambdag[krModel];
  auto lambdaTotalfunc = fSimData->fTPetroPhysics.fLambdaTotal[krModel];
  auto fwfunc = fSimData->fTPetroPhysics.fFw[krModel];
  auto fgfunc = fSimData->fTPetroPhysics.fFg[krModel];

  int nvols = this->fVolume.size();
  for (int ivol = 0; ivol < nvols; ivol++) {

    REAL sw = this->fSaturation[ivol];
    REAL rhow = this->fDensityWater[ivol];
    REAL rhoo = this->fDensityGas[ivol];
    auto fwfvalderiv = fwfunc(sw, rhow, rhoo);
    auto fgvalderiv = fgfunc(sw, rhow, rhoo);
    auto lambdaWvalderiv = lambdaWfunc(sw, rhow);
    auto lambdaGvalderiv = lambdaGfunc(sw, rhoo);
    auto lambdaTotalvalderiv = lambdaTotalfunc(sw, rhow, rhoo);
    this->fWaterfractionalflow[ivol] = std::get<0>(fwfvalderiv);
    this->fDerivativeWfractionalflow[ivol] = std::get<1>(fwfvalderiv);
    this->fGasfractionalflow[ivol] = std::get<0>(fgvalderiv);
    this->fDerivativeGfractionalflow[ivol] = std::get<1>(fgvalderiv);
    this->fLambda[ivol] = std::get<0>(lambdaTotalvalderiv);
    this->fDlambdaWaterdsw[ivol] = std::get<1>(lambdaWvalderiv);
    this->fDlambdaGasdsw[ivol] = std::get<1>(lambdaGvalderiv);
  }
}

void TSFAlgebraicTransport::TCellData::UpdateDensities() {
  int ncells = fVolume.size();
  auto WaterDensityFunc = fSimData->fTFluidProperties.fWaterDensityFunc;
  auto GasDensityFunc = fSimData->fTFluidProperties.fGasDensityFunc;

  if (WaterDensityFunc) {
    for (int icell = 0; icell < ncells; icell++) {
      REAL pressure = fPressure[icell];
      auto densityWvalderiv = WaterDensityFunc(pressure);
      auto densityGvalderiv = GasDensityFunc(pressure);
#ifdef PZDEBUG
      if (std::get<0>(densityWvalderiv) < 0.0 || std::get<0>(densityGvalderiv) < 0.0) {
        DebugStop();
      }
#endif
      fDensityWater[icell] = std::get<0>(densityWvalderiv);
      fDdensityWaterdp[icell] = std::get<1>(densityWvalderiv);
      fDensityGas[icell] = std::get<0>(densityGvalderiv);
      fDdensityGasdp[icell] = std::get<1>(densityGvalderiv);
      fVolumeFactorWater[icell] = fReferenceDensity[0] / fDensityWater[icell];
      fVolumeFactorGas[icell] = fReferenceDensity[1] / fDensityGas[icell];
    }
  }
}

void TSFAlgebraicTransport::TCellData::UpdateDensitiesLastState() {
  int ncells = fVolume.size();
  for (int icell = 0; icell < ncells; icell++) {
    fDensityWaterLastState[icell] = fDensityWater[icell];
    fDensityGasLastState[icell] = fDensityGas[icell];
  }
}

void TSFAlgebraicTransport::TCellData::UpdateMixedDensity() {
  int ncells = fVolume.size();
  for (int i = 0; i < ncells; i++) {
    REAL mixedDen = (fWaterfractionalflow[i] * fDensityWater[i]) + fGasfractionalflow[i] * fDensityGas[i];
    fMixedDensity[i] = mixedDen;
  }
}

void TSFAlgebraicTransport::TCellData::Print(std::ostream &out) {
  int nels = this->fVolume.size();
  for (int iel = 0; iel < nels; iel++) {
    out << "EqNumber: " << this->fEqNumber[iel] << std::endl;
    out << "Volume : " << this->fVolume[iel] << std::endl;
    out << "MatId: " << this->fMatId[iel] << std::endl;
    out << "GeoIndex: " << this->fGeoIndex[iel] << std::endl;
    out << "Pressure: " << this->fPressure[iel] << std::endl;
    out << "Saturation: " << this->fSaturation[iel] << std::endl;
    out << "Saturation Last State: " << this->fSaturationLastState[iel] << std::endl;
    out << "Water density: " << this->fDensityWater[iel] << std::endl;
    out << "Water density derivative: " << this->fDdensityWaterdp[iel] << std::endl;
    out << "Water density Last State: " << this->fDensityWaterLastState[iel] << std::endl;
    out << "Volume Factor Water: " << this->fVolumeFactorWater[iel] << std::endl;
    out << "Gas density: " << this->fDensityGas[iel] << std::endl;
    out << "Gas density derivative: " << this->fDdensityGasdp[iel] << std::endl;
    out << "Gas density Last State: " << this->fDensityGasLastState[iel] << std::endl;
    out << "Volume Factor Gas: " << this->fVolumeFactorGas[iel] << std::endl;
    out << "Mixed Density: " << this->fMixedDensity[iel] << std::endl;
    out << "Lambda: " << this->fLambda[iel] << std::endl;
    out << "Water Lambda derivative: " << this->fDlambdaWaterdsw[iel] << std::endl;
    out << "Gas Lambda derivative: " << this->fDlambdaGasdsw[iel] << std::endl;
    out << "Porosity: " << this->fPorosity[iel] << std::endl;
    out << "Permeability: " << this->fKappa[iel] << std::endl;
    out << "Water Fractional Flow: " << this->fWaterfractionalflow[iel] << std::endl;
    out << "Water Fractional Flow Derivative: " << this->fDerivativeWfractionalflow[iel] << std::endl;
    out << "Gas Fractional Flow: " << this->fGasfractionalflow[iel] << std::endl;
    out << "Gas Fractional Flow Derivative: " << this->fDerivativeGfractionalflow[iel] << std::endl;
    out << std::endl;
    out << "Center Cord: ";
    for (int ic = 0; ic < fCenterCoordinate[iel].size(); ic++) {
      out << fCenterCoordinate[iel][ic] << " ";
    }
    out << std::endl;
    out << "Water compressibility: " << fCompressibility[0] << std::endl;
    out << "Gas compressibility: " << fCompressibility[1] << std::endl;
    out << "Water viscosity: " << fViscosity[0] << std::endl;
    out << "Gas viscosity: " << fViscosity[1] << std::endl;
    out << "Water Reference Pressure: " << fReferencePressures[0] << std::endl;
    out << "Gas Reference Pressure: " << fReferencePressures[1] << std::endl;
    out << "Reference Water density: " << this->fReferenceDensity[0] << std::endl;
    out << "Reference Gas density: " << this->fReferenceDensity[1] << std::endl;
  }
}

void TSFAlgebraicTransport::CheckMassBalance(REAL time_step, std::ostream &out) {
  out << "Mass Balance Check:" << std::endl;
  out << "Total Water Inflow: " << fWaterMassIn << std::endl;
  out << "Total Water Outflow: " << fWaterMassOut << std::endl;
  out << "Internal Water Mass: " << fWaterMass << std::endl;
  out << "Total Gas Inflow: " << fGasMassIn << std::endl;
  out << "Total Gas Outflow: " << fGasMassOut << std::endl;
  out << "Internal Gas Mass: " << fGasMass << std::endl;
  REAL water_mass_error = fInitialWaterMass + fWaterMassIn - fWaterMassOut - fWaterMass;
  out << "Water Mass Balance Error: " << water_mass_error << std::endl;
  REAL gas_mass_error = fInitialGasMass + fGasMassIn - fGasMassOut - fGasMass;
  out << "Gas Mass Balance Error: " << gas_mass_error << std::endl;
}

void TSFAlgebraicTransport::UpdateInterfacesIntegratedFlux(int matid) {

  if (fInterfaceData.find(matid) == fInterfaceData.end()) return;
  int nels = fInterfaceData[matid].fCoefficientsFlux.size();
  if (nels == 0) return;
  std::vector<REAL> val = fInterfaceData[matid].fCoefficientsFlux[0];
  fInterfaceData[matid].fIntegralFlux = val;
  for (int i = 1; i < nels; i++) {
    int np = fInterfaceData[matid].fCoefficientsFlux[i].size();

    for (int index = 0; index < np; index++) {
      REAL val = fInterfaceData[matid].fCoefficientsFlux[i][index];
      fInterfaceData[matid].fIntegralFlux[index] += val;
    }
  }
}

void TSFAlgebraicTransport::Contribute(int cellId, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef){

    REAL sat = fCellsData.fSaturation[cellId];
    REAL satLast = fCellsData.fSaturationLastState[cellId];
    REAL densityWater = fCellsData.fDensityWater[cellId];
    REAL densityWaterLastState = fCellsData.fDensityWaterLastState[cellId];
    REAL phi = fCellsData.fPorosity[cellId];
#ifdef PZDEBUG
    if(std::abs(phi) < 1e-12) DebugStop();
#endif
    ef(0) = fCellsData.fVolume[cellId]*phi*(sat*densityWater-satLast*densityWaterLastState);
    ek(0,0) = fCellsData.fVolume[cellId]*phi*densityWater;
}

void TSFAlgebraicTransport::ContributeResidual(int cellId, TPZFMatrix<REAL> &ef){
  REAL sat = fCellsData.fSaturation[cellId];
  REAL satLast = fCellsData.fSaturationLastState[cellId];
  REAL densityWater = fCellsData.fDensityWater[cellId];
  REAL densityWaterLastState = fCellsData.fDensityWaterLastState[cellId];
  REAL phi = fCellsData.fPorosity[cellId];
  ef(0) = fCellsData.fVolume[cellId] * phi * (sat * densityWater - satLast * densityWaterLastState);
}

void TSFAlgebraicTransport::ContributeInterface(int interfaceId, int interfaceMatId, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef){
  std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceMatId].fLeftRightVolIndex[interfaceId];
    
  REAL fluxint  = 1.0*fInterfaceData[interfaceMatId].fIntegralFlux[interfaceId];
  REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
  REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];
  REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
  REAL dfwSw_R = fCellsData.fDerivativeWfractionalflow[lr_index.second];

  //upwind  matrix
  REAL beta = 0.0;
  if (fluxint > 0.0) {
    beta = 1.0;
  }

  REAL dt = fCellsData.fSimData->fTNumerics.fDt;

  ef(0) = +1.0 * (beta * fw_L + (1.0 - beta) * fw_R)*fluxint* dt;
  ef(1) = -1.0 * (beta * fw_L + (1.0 - beta) * fw_R)*fluxint* dt;
    
  ek(0,0) = +1.0 * dfwSw_L  * beta * fluxint * dt;
  ek(0,1) = +1.0 * dfwSw_R * (1.0 - beta) * fluxint * dt;
  ek(1,0) = -1.0 * dfwSw_L * beta * fluxint * dt;
  ek(1,1) = -1.0 * dfwSw_R * (1.0 - beta) * fluxint * dt;

  //IHU matrix
  REAL rhow_L = fCellsData.fDensityWater[lr_index.first];
  REAL rhow_R = fCellsData.fDensityWater[lr_index.second];
  REAL rhog_L = fCellsData.fDensityGas[lr_index.first];
  REAL rhog_R = fCellsData.fDensityGas[lr_index.second];
  REAL fg_L = fCellsData.fGasfractionalflow[lr_index.first];
  REAL fg_R = fCellsData.fGasfractionalflow[lr_index.second];
  REAL dfgSw_L = fCellsData.fDerivativeGfractionalflow[lr_index.first];
  REAL dfgSw_R = fCellsData.fDerivativeGfractionalflow[lr_index.second];
  REAL lambda_L = fCellsData.fLambda[lr_index.first];
  REAL lambda_R = fCellsData.fLambda[lr_index.second];
  REAL dlambda_L = fCellsData.fDlambdaWaterdsw[lr_index.first] + fCellsData.fDlambdaGasdsw[lr_index.first];
  REAL dlambda_R = fCellsData.fDlambdaWaterdsw[lr_index.second] + fCellsData.fDlambdaGasdsw[lr_index.second];

  std::tuple<REAL, REAL, REAL> normal = fInterfaceData[interfaceMatId].fNormalFaceDirection[interfaceId];

  TPZManVector<REAL, 3> n(3, 0.0);
  n[0] = std::get<0>(normal);
  n[1] = std::get<1>(normal);
  n[2] = std::get<2>(normal);

  TPZManVector<REAL, 3> gravity = fCellsData.fSimData->fTNumerics.fGravity;
  REAL g_dot_n = n[0] * gravity[0] + n[1] * gravity[1] + n[2] * gravity[2];
  REAL kappa =  fCellsData.fKappa[lr_index.first]; //for homogeneous permeability. PLESE IMPLEMENT HETEROGENEOUS CASE
  REAL kappa_dot_g_dot_n = kappa * g_dot_n;

  beta = 0.0;
  if (g_dot_n > 0.0) {
    beta = 1.0;
  }
  REAL fstar = beta * fw_L * fg_R + (1.0 - beta) * fw_R * fg_L;
  REAL dfstar_dswL = beta * dfwSw_L * fg_R + (1.0 - beta) * fw_R * dfgSw_L;
  REAL dfstar_dswR = beta * fw_L * dfgSw_R + (1.0 - beta) * dfwSw_R * fg_L;
  REAL lambda_star = beta * lambda_L + (1.0 - beta) * lambda_R;
  REAL dflambdaStar_dswL = beta * dlambda_L;
  REAL dflambdaStar_dswR = (1.0 - beta) * dlambda_R;

  ef(0) += fstar * lambda_star * kappa_dot_g_dot_n * (rhow_L - rhog_L) * dt;
  ef(1) += -fstar * lambda_star * kappa_dot_g_dot_n * (rhow_R - rhog_R) * dt;

  ek(0,0) += (dfstar_dswL * lambda_star + fstar * dflambdaStar_dswL) * kappa_dot_g_dot_n * (rhow_L - rhog_L) * dt;
  ek(0,1) += (dfstar_dswR * lambda_star + fstar * dflambdaStar_dswR) * kappa_dot_g_dot_n * (rhow_L - rhog_L) * dt;
  ek(1,0) += -(dfstar_dswL * lambda_star + fstar * dflambdaStar_dswL) * kappa_dot_g_dot_n * (rhow_R - rhog_R) * dt;
  ek(1,1) += -(dfstar_dswR * lambda_star + fstar * dflambdaStar_dswR) * kappa_dot_g_dot_n * (rhow_R - rhog_R) * dt;
}

void TSFAlgebraicTransport::ContributeInterfaceResidual(int interfaceId, int interfaceMatId, TPZFMatrix<REAL> &ef){
  std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceMatId].fLeftRightVolIndex[interfaceId];
  REAL fluxint  = fInterfaceData[interfaceMatId].fIntegralFlux[interfaceId];
  REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
  REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];

  REAL dt = fCellsData.fSimData->fTNumerics.fDt;

  //upwind
  REAL beta =0.0;
  if (fluxint>0.0) {
    beta = 1.0;
  }
    
  ef(0) = +1.0*(beta * fw_L + (1.0 - beta) * fw_R) * fluxint * dt;
  ef(1) = -1.0*(beta * fw_L  + (1.0 - beta) * fw_R) * fluxint * dt;

  //IHU
  REAL rhow_L = fCellsData.fDensityWater[lr_index.first];
  REAL rhow_R = fCellsData.fDensityWater[lr_index.second];
  REAL rhog_L = fCellsData.fDensityGas[lr_index.first];
  REAL rhog_R = fCellsData.fDensityGas[lr_index.second];
  REAL fg_L = fCellsData.fGasfractionalflow[lr_index.first];
  REAL fg_R = fCellsData.fGasfractionalflow[lr_index.second];
  REAL lambda_L = fCellsData.fLambda[lr_index.first];
  REAL lambda_R = fCellsData.fLambda[lr_index.second];

  std::tuple<REAL, REAL, REAL> normal = fInterfaceData[interfaceMatId].fNormalFaceDirection[interfaceId];

  TPZManVector<REAL, 3> n(3, 0.0);
  n[0] = std::get<0>(normal);
  n[1] = std::get<1>(normal);
  n[2] = std::get<2>(normal);

  TPZManVector<REAL, 3> gravity = fCellsData.fSimData->fTNumerics.fGravity;
  REAL g_dot_n = n[0] * gravity[0] + n[1] * gravity[1] + n[2] * gravity[2];
  REAL kappa =  fCellsData.fKappa[lr_index.first]; //for homogeneous permeability. PLESE IMPLEMENT HETEROGENEOUS CASE
  REAL kappa_dot_g_dot_n = kappa * g_dot_n;

  beta = 0.0;
  if (g_dot_n > 0.0) {
    beta = 1.0;
  }
  REAL fstar = beta * fw_L * fg_R + (1.0 - beta) * fw_R * fg_L;
  REAL lambda_star = beta * lambda_L + (1.0 - beta) * lambda_R;

  ef(0) += fstar * lambda_star * kappa_dot_g_dot_n * (rhow_L - rhog_L) * dt;
  ef(1) += -fstar * lambda_star * kappa_dot_g_dot_n * (rhow_R - rhog_R) * dt;
}