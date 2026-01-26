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