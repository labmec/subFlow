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

void TSFProblemData::ReadJSONFile(std::string filename) {
  using json = nlohmann::json;

  std::string jsonpath = std::string(INPUTDIR) + "/" + filename;
  std::ifstream filejson(jsonpath);
  json input = json::parse(filejson, nullptr, true, true); // to ignore comments in json file

  // ------------------------ Getting number of domains and fractures ------------------------
  if (input.find("UseGMsh") == input.end()) DebugStop();
  fTGeometry.fUseGMsh = input["UseGMsh"];
  if (input["UseGMsh"])
    fTGeometry.fGmshFile = input["MshFile"];
  if (input.find("Domains") == input.end()) DebugStop();
  fTGeometry.fDimension = input["Dimension"];

  // ------------------------ Reading 3D Domain matids ------------------------
  for (auto &domain : input["Domains"]) {
    if (domain.find("matid") == domain.end()) DebugStop();
    if (domain.find("name") == domain.end()) DebugStop();
    if (domain.find("K") == domain.end()) DebugStop();
    if (domain.find("phi") == domain.end()) DebugStop();
    const int matid = domain["matid"];
    const std::string name = domain["name"];
    const REAL permeability = domain["K"];
    const REAL phi = domain["phi"];
    fTGeometry.fDomainNameAndMatId[name] = matid;
    fTReservoirProperties.fPorosityAndPermeability[matid] = std::make_pair(phi, permeability);
  }

  // ------------------------ Reading 3D Domain BC matids ------------------------
  if (input.find("Boundary") == input.end()) DebugStop();
  std::map<int, std::pair<int, REAL>> &BCDarcyMatIdToTypeValue = fTBoundaryConditions.fBCDarcyMatIdToTypeValue;
  std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> &BCDarcyMatIdToFunctionId = fTBoundaryConditions.fBCDarcyMatIdToFunctionId;
  std::map<int, std::pair<int, REAL>> &BCTransportMatIdToTypeValue = fTBoundaryConditions.fBCTransportMatIdToTypeValue;
  std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> &BCTransportMatIdToFunctionId = fTBoundaryConditions.fBCTransportMatIdToFunctionId;
  for (auto &bc : input["Boundary"]) {
    if (bc.find("matid") == bc.end()) DebugStop();
    if (bc.find("type") == bc.end()) DebugStop();
    if (bc.find("value") == bc.end()) DebugStop();
    if (bc.find("name") == bc.end()) DebugStop();
    const int matid = bc["matid"];
    const int type = bc["type"];
    const REAL value = bc["value"];
    const std::string name = bc["name"];
    int functionID = 0;
    if (bc.find("functionID") != bc.end()) {
      functionID = bc["functionID"];
    }
    REAL external_saturation = 0.0;
    int saturation_functionID = 0;
    if (input["Numerics"]["AnalysisType"] != 0) { // If transport problem is being solved
      if (bc.find("ExternalSaturation") == bc.end()) DebugStop();
      external_saturation = bc["ExternalSaturation"];
      if (bc.find("SaturationFunctionID") != bc.end()) {
        saturation_functionID = bc["SaturationFunctionID"];
      }
    }
    fTBoundaryConditions.fDomainNameAndMatId[name] = matid;
    if (BCDarcyMatIdToTypeValue.find(matid) != BCDarcyMatIdToTypeValue.end()) DebugStop();
    BCDarcyMatIdToTypeValue[matid] = std::make_pair(type, value);
    TSFFunctionsGenerator functionGen;
    TSFFunctionsGenerator::EDarcyBCFunctionType darcy_functionType = bc["functionID"];
    functionGen.SetDarcyBCFuncType(darcy_functionType);
    auto darcyBCfunc = functionGen.CreateDarcyBC();
    BCDarcyMatIdToFunctionId[matid] = std::make_pair(functionID, darcyBCfunc);
    TSFFunctionsGenerator::ETransportBCFunctionType transport_functionType = bc["SaturationFunctionID"];
    functionGen.SetTransportBCFuncType(transport_functionType);
    auto transportBCfunc = functionGen.CreateTransportBC();
    BCTransportMatIdToTypeValue[matid] = std::make_pair(type, external_saturation);
    BCTransportMatIdToFunctionId[matid] = std::make_pair(saturation_functionID, transportBCfunc);
  }

  // ------------------------ Numerics Parameters ------------------------
  if (input.find("Numerics") != input.end()) {
    auto numerics = input["Numerics"];
    fTNumerics.fAnalysisType = numerics["AnalysisType"];
    if (numerics.find("DeltaT") == numerics.end()) DebugStop();
    fTNumerics.fDt = numerics["DeltaT"];
    if (numerics.find("NSteps") == numerics.end()) DebugStop();
    fTNumerics.fNSteps = numerics["NSteps"];
    if (numerics.find("Gravity") == numerics.end()) DebugStop();
    TPZManVector<REAL, 3> grav(3, 0.0);
    for (int i = 0; i < 3; i++) {
      grav[i] = numerics["Gravity"][i];
    }
    fTNumerics.fGravity = grav;

    fTNumerics.fResTolTransport = 1.e-6;
    fTNumerics.fCorrTolTransport = 1.e-6;
    fTNumerics.fResTolDarcy = 1.e-6;
    fTNumerics.fCorrTolDarcy = 1.e-6;
    fTNumerics.fFourApproxSpaces = true;
    fTNumerics.fMaxIterSFI = 10;
    fTNumerics.fTolSFI = 1.e-6;
    fTNumerics.fMaxIterDarcy = 10;
    fTNumerics.fMaxIterTransport = 10;
    if (numerics.find("IsAxisymmetric") != numerics.end()) {
      fTNumerics.fIsAxisymmetric = numerics["IsAxisymmetric"];
    }
    if (numerics.find("IsLinearTrace") != numerics.end()) {
      fTNumerics.fIsLinearTrace = numerics["IsLinearTrace"];
    }
    if (numerics.find("FourApproxSpaces") != numerics.end()) {
      fTNumerics.fFourApproxSpaces = numerics["FourApproxSpaces"];
    }
    if (numerics.find("NThreadsDarcy") != numerics.end()) {
      fTNumerics.fNThreadsDarcy = numerics["NThreadsDarcy"];
    }
    if (numerics.find("MaxIterSFI") != numerics.end()) {
      fTNumerics.fMaxIterSFI = numerics["MaxIterSFI"];
    }
    if (numerics.find("TolSFI") != numerics.end()) {
      fTNumerics.fTolSFI = numerics["TolSFI"];
    }
    if (numerics.find("MaxIterDarcy") != numerics.end()) {
      fTNumerics.fMaxIterDarcy = numerics["MaxIterDarcy"];
    }
    if (numerics.find("ResTolDarcy") != numerics.end()) {
      fTNumerics.fResTolDarcy = numerics["ResTolDarcy"];
    }
    if (numerics.find("CorrTolDarcy") != numerics.end()) {
      fTNumerics.fCorrTolDarcy = numerics["CorrTolDarcy"];
    }
    if (numerics.find("MaxIterTransport") != numerics.end()) {
      fTNumerics.fMaxIterTransport = numerics["MaxIterTransport"];
    }
    if (numerics.find("ResTolTransport") != numerics.end()) {
      fTNumerics.fResTolTransport = numerics["ResTolTransport"];
    }
    if (numerics.find("CorrTolTransport") != numerics.end()) {
      fTNumerics.fCorrTolTransport = numerics["CorrTolTransport"];
    }
  }

  // ------------------------ Fluids Properties ------------------------
  if (input.find("FluidProperties") != input.end()) {
    auto properties = input["FluidProperties"];
    if (properties.find("WaterDensity") == properties.end()) DebugStop();
    fTFluidProperties.fWaterDensityRef = properties["WaterDensity"];
    if (properties.find("WaterViscosity") == properties.end()) DebugStop();
    fTFluidProperties.fWaterViscosity = properties["WaterViscosity"];
    if (properties.find("WaterCompressibility") == properties.end()) DebugStop();
    fTFluidProperties.fWaterCompressibility = properties["WaterCompressibility"];
    if (properties.find("GasDensity") == properties.end()) DebugStop();
    fTFluidProperties.fGasDensityRef = properties["GasDensity"];
    if (properties.find("GasViscosity") == properties.end()) DebugStop();
    fTFluidProperties.fGasViscosity = properties["GasViscosity"];
    if (properties.find("GasCompressibility") == properties.end()) DebugStop();
    fTFluidProperties.fGasCompressibility = properties["GasCompressibility"];
    if (properties.find("DensityModel") == properties.end()) DebugStop();
    if (properties["DensityModel"] == 0) {
      fTFluidProperties.CreateLinearDensityFunction();
    } else {
      fTFluidProperties.CreateExponentialDensityFunction();
    }
    if (properties.find("ReferencePressure") != properties.end()) {
      fTFluidProperties.fReferencePressure = properties["ReferencePressure"];
    }
  }

  // ------------------------ Petro Physics ------------------------
  if (input.find("PetroPhysics") != input.end()) {
    auto petro = input["PetroPhysics"];
    if (petro.find("KrModel") == petro.end()) DebugStop();
    fTPetroPhysics.fKrModel = petro["KrModel"];
    fTPetroPhysics.fSwr = petro["Swr"];
    fTPetroPhysics.fSgr = petro["Sgr"];
    fTPetroPhysics.CreateLinearKrModel();
    fTPetroPhysics.CreateQuadraticKrModel(); // It is necessary to call these methods after the residual saturations are set
    fTPetroPhysics.fWaterViscosity = fTFluidProperties.fWaterViscosity;
    fTPetroPhysics.fGasViscosity = fTFluidProperties.fGasViscosity;
  }

  // ------------------------ Reservoir Properties ------------------------
  if (input.find("ReservoirProperties") != input.end()) {
    auto reservoir = input["ReservoirProperties"];
    if (reservoir.find("s0") != reservoir.end()) {
      auto s0 = reservoir["s0"];
      TSFFunctionsGenerator::ES0FunctionType s0_functionType = s0["functionType"];
      REAL constVal = s0["value"];
      TSFFunctionsGenerator functionGen;
      functionGen.SetS0FuncType(s0_functionType, constVal);
      auto s0func = functionGen.CreateS0();
      fTReservoirProperties.fS0 = s0func;
    }
  }

  // ------------------------ Post Processing ------------------------
  fTPostProcess.fFileNameDarcy = "postdarcy.vtk";
  fTPostProcess.fFileNameTransport = "posttransport.vtk";
  if (input.find("PostProcess") != input.end()) {
    auto postprocess = input["PostProcess"];
    if (postprocess.find("DarcyFileName") != postprocess.end()) {
      fTPostProcess.fFileNameDarcy = postprocess["DarcyFileName"];
    }
    if (postprocess.find("TransportFileName") != postprocess.end()) {
      fTPostProcess.fFileNameTransport = postprocess["TransportFileName"];
    }
    if (postprocess.find("PostProcessFrequency") != postprocess.end()) {
      fTPostProcess.fPostProcessFrequency = postprocess["PostProcessFrequency"];
    }
    if (postprocess.find("NThreads") != postprocess.end()) {
      fTPostProcess.fNThreads = postprocess["NThreads"];
    }
    if (postprocess.find("VTKResolution") != postprocess.end()) {
      fTPostProcess.fvtkResolution = postprocess["VTKResolution"];
    }
  }
  TPZStack<std::string, 10> scalnames, vecnames, scalnamesTransport;
  vecnames.Push("Flux");
  scalnames.Push("Pressure");
  scalnames.Push("DivFlux");
  if (fTNumerics.fFourApproxSpaces) {
    scalnames.Push("g_average");
    scalnames.Push("p_average");
  }
  scalnamesTransport.Push("Sw");
  scalnamesTransport.Push("Sg");

  fTPostProcess.fVecnamesDarcy = vecnames;
  fTPostProcess.fScalnamesDarcy = scalnames;
  fTPostProcess.fScalnamesTransport = scalnamesTransport;

  int n_steps = fTNumerics.fNSteps;
  fTPostProcess.fFileTimeStep = fTNumerics.fDt;
  REAL dt = fTNumerics.fDt;
  REAL time = fTPostProcess.fFileTimeStep;
  const int freq = fTPostProcess.fPostProcessFrequency;
  int n_reporting_times = (n_steps) / (time * freq / dt) + 1;
  TPZStack<REAL, 100> reporting_times;
  REAL r_time = 0.0;
  int j = freq;
  for (int i = 1; i <= n_reporting_times; i++) {

    r_time = j * dt * (time / dt);
    reporting_times.push_back(r_time);
    j += freq;
  }
  fTPostProcess.fVecReportingTimes = reporting_times;
}
