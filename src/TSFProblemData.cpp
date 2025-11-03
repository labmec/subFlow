//
//  Created by Giovane Avancini on 05/09/25.
//

#include "TSFProblemData.h"

TSFProblemData::TSFProblemData() : mTGeometry(), mTPetroPhysics(), mTFluidProperties(), mTReservoirProperties(), mTBoundaryConditions(), mTNumerics(), mTPostProcess() {}

TSFProblemData::TSFProblemData(const TSFProblemData &other) {
  mTGeometry = other.mTGeometry;
  mTPetroPhysics = other.mTPetroPhysics;
  mTFluidProperties = other.mTFluidProperties;
  mTReservoirProperties = other.mTReservoirProperties;
  mTBoundaryConditions = other.mTBoundaryConditions;
  mTNumerics = other.mTNumerics;
  mTPostProcess = other.mTPostProcess;
  mSimulationName = other.mSimulationName;
}

TSFProblemData &TSFProblemData::operator=(const TSFProblemData &other) {

  if (this != &other) // prevent self-assignment
  {
    mTGeometry = other.mTGeometry;
    mTPetroPhysics = other.mTPetroPhysics;
    mTFluidProperties = other.mTFluidProperties;
    mTReservoirProperties = other.mTReservoirProperties;
    mTBoundaryConditions = other.mTBoundaryConditions;
    mTNumerics = other.mTNumerics;
    mTPostProcess = other.mTPostProcess;
    mSimulationName = other.mSimulationName;
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
  mWaterDensityF = [this](REAL &p) {
    REAL dp = p - mReferencePressure;
    REAL rho = mWaterDensityRef * (1 + mWaterCompressibility * dp);
    REAL drho_dp = mWaterCompressibility * mWaterDensityRef;
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };
  mOilDensityF = [this](REAL &p) {
    REAL dp = p - mReferencePressure;
    REAL rho = mOilDensityRef * (1 + mOilCompressibility * dp);
    REAL drho_dp = mOilCompressibility * mOilDensityRef;
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };
}
void TSFProblemData::TFluidProperties::CreateExponentialDensityFunction() {
  mWaterDensityF = [this](REAL &p) {
    REAL dp = p - mReferencePressure;
    REAL rho = mWaterDensityRef * std::exp(std::exp((mWaterCompressibility * (dp))));

    REAL drho_dp = mWaterCompressibility * mWaterDensityRef * std::exp(std::exp((mWaterCompressibility * (dp))));
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };

  mOilDensityF = [this](REAL &p) {
    REAL dp = p - mReferencePressure;
    REAL rho = mOilDensityRef * std::exp(((mOilCompressibility * (dp))));

    REAL drho_dp = mOilCompressibility * mOilDensityRef * std::exp(((mOilCompressibility * (dp))));
    std::tuple<REAL, REAL> valderiv(rho, drho_dp);
    return valderiv;
  };
}

void TSFProblemData::TPetroPhysics::CreateLinearKrModel() {

  mKrw[0] = [](REAL &sw) {
    REAL krw = sw;
    REAL dkrw = 1;
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };
  mKro[0] = [](REAL &sw) {
    REAL krw = (1 - sw);
    REAL dkrw = -1;
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };
  UpdateLambdasAndFracFlows(0);
}
void TSFProblemData::TPetroPhysics::CreateQuadraticKrModel() {

  mKrw[1] = [](REAL &sw) {
    REAL krw = sw * sw;
    REAL dkrw = 2 * sw;
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };
  mKro[1] = [](REAL &sw) {
    REAL krw = (1 - sw) * (1 - sw);
    REAL dkrw = -2 * (1 - sw);
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };
  UpdateLambdasAndFracFlows(1);
}
void TSFProblemData::TPetroPhysics::CreateQuadraticResidualKrModel() {
  REAL swr = mSwr;
  REAL sor = mSor;
  mKrw[2] = [swr](REAL &sw) {
    REAL krw = 0;
    REAL dkrw = 0;
    if (sw > swr) {
      krw = (sw - swr) * (sw - swr) / ((1. - swr) * (1. - swr));
      dkrw = (2 * (sw - swr)) / ((1. - swr) * (1. - swr));
    }
    std::tuple<REAL, REAL> valderiv(krw, dkrw);
    return valderiv;
  };
  mKro[2] = [sor](REAL &sw) {
    REAL so = 1 - sw;
    REAL kro = 0;
    REAL dkro = 0;
    if (so > sor) {
      kro = (so - sor) * (so - sor) / ((1. - sor) * (1. - sor));
      dkro = -(2 * (so - sor)) / ((1. - sor) * (1. - sor));
    }
    std::tuple<REAL, REAL> valderiv(kro, dkro);
    return valderiv;
  };
  UpdateLambdasAndFracFlows(2);
}

void TSFProblemData::TPetroPhysics::UpdateLambdasAndFracFlows(int krModel) {
  mLambdaW[krModel] = [this, krModel](REAL &sw, REAL &rhow) {
    std::tuple<REAL, REAL> krwvalderiv = mKrw[krModel](sw);
    auto krw = std::get<0>(krwvalderiv);
    REAL dkrw = std::get<1>(krwvalderiv);
    REAL lambdaw = krw / mWaterViscosity;
    REAL dlambdadsw = dkrw / mWaterViscosity;
    std::tuple<REAL, REAL> valderiv(lambdaw, dlambdadsw);
    return valderiv;
  };
  mLambdaO[krModel] = [this, krModel](REAL &sw, REAL rhoo) {
    std::tuple<REAL, REAL> krwvalderiv = mKro[krModel](sw);
    auto kro = std::get<0>(krwvalderiv);
    REAL dkro = std::get<1>(krwvalderiv);
    REAL lambdao = kro / mOilViscosity;
    REAL dlambdadso = dkro / mOilViscosity;
    std::tuple<REAL, REAL> valderiv(lambdao, dlambdadso);
    return valderiv;
  };

  mLambdaTotal[krModel] = [this, krModel](REAL &sw, REAL &rhow, REAL &rhoo) {
    std::tuple<REAL, REAL> lambdaWvalderiv = mLambdaW[krModel](sw, rhow);
    std::tuple<REAL, REAL> lambdaOvalderiv = mLambdaO[krModel](sw, rhoo);
    REAL lw = std::get<0>(lambdaWvalderiv);
    REAL dlwdsw = std::get<1>(lambdaWvalderiv);
    REAL lo = std::get<0>(lambdaOvalderiv);
    REAL dlodsw = std::get<1>(lambdaOvalderiv);
    REAL lambdaTotal = lw + lo;
    REAL dlambdaTotaldsw = dlwdsw + dlodsw;
    std::tuple<REAL, REAL> valderiv(lambdaTotal, dlambdaTotaldsw);
    return valderiv;
  };
  mFw[krModel] = [this, krModel](REAL &sw, REAL &rhow, REAL &rhoo) {
    std::tuple<REAL, REAL> lambdaWvalderiv = mLambdaW[krModel](sw, rhow);
    std::tuple<REAL, REAL> lambdaTotalvalderiv = mLambdaTotal[krModel](sw, rhow, rhoo);
    REAL lw = std::get<0>(lambdaWvalderiv);
    REAL dlwdsw = std::get<1>(lambdaWvalderiv);
    REAL ltotal = std::get<0>(lambdaTotalvalderiv);
    REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
    REAL fracflow = lw / ltotal;
    REAL dfracflowdsw = (dlwdsw / ltotal) - ((lw * dltotaldsw) / (ltotal * ltotal));
    std::tuple<REAL, REAL> valderiv(fracflow, dfracflowdsw);
    return valderiv;
  };
  mFo[krModel] = [this, krModel](REAL &sw, REAL &rhow, REAL &rhoo) {
    std::tuple<REAL, REAL> lambdaOvalderiv = mLambdaO[krModel](sw, rhoo);
    std::tuple<REAL, REAL> lambdaTotalvalderiv = mLambdaTotal[krModel](sw, rhow, rhoo);
    REAL lo = std::get<0>(lambdaOvalderiv);
    REAL dlodso = std::get<1>(lambdaOvalderiv);
    REAL ltotal = std::get<0>(lambdaTotalvalderiv);
    REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
    REAL fracflow = lo / ltotal;
    REAL dfracflowdsw = (dlodso / ltotal) - ((lo * dltotaldsw) / (ltotal * ltotal));
    std::tuple<REAL, REAL> valderiv(fracflow, dfracflowdsw);
    return valderiv;
  };
}
