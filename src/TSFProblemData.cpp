//
//  Created by Giovane Avancini on 05/09/25.
//

#include "TSFProblemData.h"

TSFProblemData::TSFProblemData() : fTGeometry(), fTPetroPhysics(), fTFluidProperties(), fTReservoirProperties(), fTBoundaryConditions(), fTNumerics(), fTPostProcess() {}

TSFProblemData::TSFProblemData(const TSFProblemData &other) {
  fTGeometry = other.fTGeometry;
  fTPetroPhysics = other.fTPetroPhysics;
  fTFluidProperties = other.fTFluidProperties;
  fTReservoirProperties = other.fTReservoirProperties;
  fTBoundaryConditions = other.fTBoundaryConditions;
  fTNumerics = other.fTNumerics;
  fTPostProcess = other.fTPostProcess;
  fSimulationName = other.fSimulationName;
}

TSFProblemData &TSFProblemData::operator=(const TSFProblemData &other) {
  if (this != &other) // prevent self-assignment
  {
    fTGeometry = other.fTGeometry;
    fTPetroPhysics = other.fTPetroPhysics;
    fTFluidProperties = other.fTFluidProperties;
    fTReservoirProperties = other.fTReservoirProperties;
    fTBoundaryConditions = other.fTBoundaryConditions;
    fTNumerics = other.fTNumerics;
    fTPostProcess = other.fTPostProcess;
    fSimulationName = other.fSimulationName;
  }
  return *this;
}

TSFProblemData::~TSFProblemData() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void TSFProblemData::Write(TPZStream &buf, int withclassid) const {
  DebugStop();
}

void TSFProblemData::Read(TPZStream &buf, void *context) {
  DebugStop();
}

int TSFProblemData::ClassId() const {
  DebugStop();
  return 1;
}

void TSFProblemData::TFluidProperties::CreateLinearDensityFunction() {
  fWaterDensityFunc = [this](REAL &p) {
    REAL dp = p - fReferencePressure;
    REAL rho = fWaterDensityRef * (1. + fWaterCompressibility * dp);
    REAL drho_dp = fWaterCompressibility * fWaterDensityRef;
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };
  fGasDensityFunc = [this](REAL &p) {
    REAL dp = p - fReferencePressure;
    REAL rho = fGasDensityRef * (1. + fGasCompressibility * dp);
    REAL drho_dp = fGasCompressibility * fGasDensityRef;
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };
}

void TSFProblemData::TFluidProperties::CreateExponentialDensityFunction() {
  fWaterDensityFunc = [this](REAL &p) {
    REAL dp = p - fReferencePressure;
    REAL rho = fWaterDensityRef * std::exp(std::exp((fWaterCompressibility * (dp))));

    REAL drho_dp = fWaterCompressibility * fWaterDensityRef * std::exp(std::exp((fWaterCompressibility * (dp))));
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };

  fGasDensityFunc = [this](REAL &p) {
    REAL dp = p - fReferencePressure;
    REAL rho = fGasDensityRef * std::exp(((fGasCompressibility * (dp))));

    REAL drho_dp = fGasCompressibility * fGasDensityRef * std::exp(((fGasCompressibility * (dp))));
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };
}

void TSFProblemData::TPetroPhysics::CreateLinearKrModel() {
  REAL swr = fSwr;
  REAL sgr = fSgr;

  fKrw[0] = [swr](REAL &sw) {
    REAL krw = 0;
    REAL dkrw = 0;
    if (sw > swr) {
      krw = (sw - swr) / (1. - swr);
      dkrw = 1.0 / (1. - swr);
    }
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };

  fKrg[0] = [sgr](REAL &sw) {
    REAL sg = 1 - sw;
    REAL krg = 0;
    REAL dkrg = 0;
    if (sg > sgr) {
      krg = (sg - sgr) / (1. - sgr);
      dkrg = -1.0 / (1. - sgr);
    }
    std::tuple<REAL, REAL> valderiv(krg, dkrg);
    return valderiv;
  };
  UpdateLambdasAndFracFlows(0);
}

void TSFProblemData::TPetroPhysics::CreateQuadraticKrModel() {
  REAL swr = fSwr;
  REAL sgr = fSgr;

  fKrw[2] = [swr](REAL &sw) {
    REAL krw = 0;
    REAL dkrw = 0;
    if (sw > swr) {
      krw = (sw - swr) * (sw - swr) / ((1. - swr) * (1. - swr));
      dkrw = (2 * (sw - swr)) / ((1. - swr) * (1. - swr));
    }
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };

  fKrg[2] = [sgr](REAL &sw) {
    REAL sg = 1 - sw;
    REAL krg = 0;
    REAL dkrg = 0;
    if (sg > sgr) {
      krg = (sg - sgr) * (sg - sgr) / ((1. - sgr) * (1. - sgr));
      dkrg = -(2 * (sg - sgr)) / ((1. - sgr) * (1. - sgr));
    }
    std::tuple<REAL, REAL> valderiv(krg, dkrg);
    return valderiv;
  };
  UpdateLambdasAndFracFlows(1);
}

void TSFProblemData::TPetroPhysics::UpdateLambdasAndFracFlows(int krModel) {
  fLambdaw[krModel] = [this, krModel](REAL &sw, REAL &rhow) {
    std::tuple<REAL, REAL> krwvalderiv = fKrw[krModel](sw);
    auto krw = std::get<0>(krwvalderiv);
    REAL dkrw = std::get<1>(krwvalderiv);
    REAL lambdaw = krw / fWaterViscosity;
    REAL dlambdadsw = dkrw / fWaterViscosity;
    std::tuple<REAL, REAL> valderiv(lambdaw, dlambdadsw);
    return valderiv;
  };

  fLambdag[krModel] = [this, krModel](REAL &sw, REAL &rhog) {
    std::tuple<REAL, REAL> krwvalderiv = fKrg[krModel](sw);
    auto kro = std::get<0>(krwvalderiv);
    REAL dkro = std::get<1>(krwvalderiv);
    REAL lambdao = kro / fGasViscosity;
    REAL dlambdadso = dkro / fGasViscosity;
    std::tuple<REAL, REAL> valderiv(lambdao, dlambdadso);
    return valderiv;
  };

  fLambdaTotal[krModel] = [this, krModel](REAL &sw, REAL &rhow, REAL &rhog) {
    std::tuple<REAL, REAL> lambdaWvalderiv = fLambdaw[krModel](sw, rhow);
    std::tuple<REAL, REAL> lambdaOvalderiv = fLambdag[krModel](sw, rhog);
    REAL lw = std::get<0>(lambdaWvalderiv);
    REAL dlwdsw = std::get<1>(lambdaWvalderiv);
    REAL lg = std::get<0>(lambdaOvalderiv);
    REAL dlgdsw = std::get<1>(lambdaOvalderiv);
    REAL lambdaTotal = lw + lg;
    REAL dlambdaTotaldsw = dlwdsw + dlgdsw;
    std::tuple<REAL, REAL> valderiv(lambdaTotal, dlambdaTotaldsw);
    return valderiv;
  };

  fFw[krModel] = [this, krModel](REAL &sw, REAL &rhow, REAL &rhog) {
    std::tuple<REAL, REAL> lambdaWvalderiv = fLambdaw[krModel](sw, rhow);
    std::tuple<REAL, REAL> lambdaTotalvalderiv = fLambdaTotal[krModel](sw, rhow, rhog);
    REAL lw = std::get<0>(lambdaWvalderiv);
    REAL dlwdsw = std::get<1>(lambdaWvalderiv);
    REAL ltotal = std::get<0>(lambdaTotalvalderiv);
    REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
    REAL fracflow = lw / ltotal;
    REAL dfracflowdsw = (dlwdsw / ltotal) - ((lw * dltotaldsw) / (ltotal * ltotal));
    std::tuple<REAL, REAL> valderiv(fracflow, dfracflowdsw);
    return valderiv;
  };

  fFg[krModel] = [this, krModel](REAL &sw, REAL &rhow, REAL &rhog) {
    std::tuple<REAL, REAL> lambdaOvalderiv = fLambdag[krModel](sw, rhog);
    std::tuple<REAL, REAL> lambdaTotalvalderiv = fLambdaTotal[krModel](sw, rhow, rhog);
    REAL lg = std::get<0>(lambdaOvalderiv);
    REAL dlgdsw = std::get<1>(lambdaOvalderiv);
    REAL ltotal = std::get<0>(lambdaTotalvalderiv);
    REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
    REAL fracflow = lg / ltotal;
    REAL dfracflowdsw = (dlgdsw / ltotal) - ((lg * dltotaldsw) / (ltotal * ltotal));
    std::tuple<REAL, REAL> valderiv(fracflow, dfracflowdsw);
    return valderiv;
  };
}
