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

  std::ifstream filejson(filename);
  json input = json::parse(filejson, nullptr, true, true); // to ignore comments in json file

  // ------------------------ Getting number of domains and fractures ------------------------
  if (input.find("Domains") == input.end()) DebugStop();
  if (input.find("Mesh") == input.end()) DebugStop();
  const int ndom = input["Domains"].size();
  std::string mesh = input["Mesh"];
  fTGeometry.fGmeshFileName = mesh;

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
    if (input["Numerics"]["RunWithTransport"]) {
      if (bc.find("ExternalSaturation") == bc.end()) DebugStop();
      external_saturation = bc["ExternalSaturation"];
      if (bc.find("SaturationFunctionID") != bc.end()) {
        saturation_functionID = bc["SaturationFunctionID"];
      }
    }
    fTBoundaryConditions.fDomainNameAndMatId[name] = matid;
    if (BCDarcyMatIdToTypeValue.find(matid) != BCDarcyMatIdToTypeValue.end()) DebugStop();
    BCDarcyMatIdToTypeValue[matid] = std::make_pair(type, value);
    // BCDarcyMatIdToFunctionId[matid] = std::make_pair(functionID, forcingfunctionBC[functionID]); //REMEMBER TO FIX IT
    BCTransportMatIdToTypeValue[matid] = std::make_pair(type, external_saturation);
    // BCTransportMatIdToFunctionId[matid] = std::make_pair(saturation_functionID, forcingfunctionBC[saturation_functionID]); //REMEMBER TO FIX IT
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

    if (numerics.find("IsAxisymmetric") != numerics.end()) {
      fTNumerics.fIsAxisymmetric = numerics["IsAxisymmetric"];
    }
    if (numerics.find("IsLinearTrace") != numerics.end()) {
      fTNumerics.fIsLinearTrace = numerics["IsLinearTrace"];
    }
  }

  // ------------------------ Fluids Properties ------------------------
  if (input.find("FluidProperties") != input.end()) {
    auto properties = input["FluidProperties"];
    if (properties.find("WaterDensity") == properties.end()) DebugStop();
    sim_data.mTFluidProperties.mWaterDensityRef = properties["WaterDensity"];
    if (properties.find("WaterViscosity") == properties.end()) DebugStop();
    sim_data.mTFluidProperties.mWaterViscosity = properties["WaterViscosity"];
    if (properties.find("WaterCompressibility") == properties.end()) DebugStop();
    sim_data.mTFluidProperties.mWaterCompressibility = properties["WaterCompressibility"];
    if (properties.find("OilDensity") == properties.end()) DebugStop();
    sim_data.mTFluidProperties.mOilDensityRef = properties["OilDensity"];
    if (properties.find("OilViscosity") == properties.end()) DebugStop();
    sim_data.mTFluidProperties.mOilViscosity = properties["OilViscosity"];
    if (properties.find("OilCompressibility") == properties.end()) DebugStop();
    sim_data.mTFluidProperties.mOilCompressibility = properties["OilCompressibility"];
    if (properties.find("DensityModel") == properties.end()) DebugStop();
    if (properties["DensityModel"] == 0) {
      sim_data.mTFluidProperties.CreateLinearDensityFunction();
    } else {
      sim_data.mTFluidProperties.CreateExponentialDensityFunction();
    }
    if (properties.find("ReferencePressure") != properties.end()) {
      sim_data.mTFluidProperties.mReferencePressure = properties["ReferencePressure"];
    }
  }

  // ------------------------ Petro Physics ------------------------
  if (input.find("PetroPhysics") != input.end()) {
    auto petro = input["PetroPhysics"];
    if (petro.find("KrModel") == petro.end()) DebugStop();
    sim_data.mTPetroPhysics.mKrModel = petro["KrModel"];
    if (petro["KrModel"] == 2) {
      if (petro.find("Swr") == petro.end()) DebugStop();
      if (petro.find("Sor") == petro.end()) DebugStop();
      sim_data.mTPetroPhysics.mSwr = petro["Swr"];
      sim_data.mTPetroPhysics.mSor = petro["Sor"];
      sim_data.mTPetroPhysics.CreateQuadraticResidualKrModel(); // It is necessary to call this method after the residual saturations are set
    }
    sim_data.mTPetroPhysics.mWaterViscosity = sim_data.mTFluidProperties.mWaterViscosity;
    sim_data.mTPetroPhysics.mOilViscosity = sim_data.mTFluidProperties.mOilViscosity;
  }

  // ------------------------ Reservoir Properties ------------------------
  if (input.find("ReservoirProperties") != input.end()) {
    auto reservoir = input["ReservoirProperties"];
    if (reservoir.find("s0") != reservoir.end()) {
      auto s0 = reservoir["s0"];
      TMRSPropertiesFunctions::EFunctionType s0_functionType = s0["functionType"];
      REAL constant_val = s0["value"];
      TMRSPropertiesFunctions reservoir_properties;
      reservoir_properties.set_function_type_s0(s0_functionType, constant_val);
      auto s0_function = reservoir_properties.Create_s0();
      sim_data.mTReservoirProperties.s0 = s0_function;
    }
  }

  // ------------------------ Setting extra stuff that is still not in JSON ------------------------
  const int D_Type = 0, N_Type = 1, Mixed_Type = 2;
  // sim_data.mTGeometry.mInterface_material_id = 100;
  // sim_data.mTGeometry.mInterface_material_idFracInf = 102;
  // sim_data.mTGeometry.mInterface_material_idFracSup = 101;
  // sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
  // sim_data.mTGeometry.mInterface_material_idFracBound = 104;

  // sim_data.mTGeometry.mSkeletonDiv = 0;
  sim_data.mTNumerics.m_sfi_tol = 1.e-8;
  sim_data.mTNumerics.m_res_tol_transport = 1.e-8;
  sim_data.mTNumerics.m_corr_tol_transport = 1.e-8;
  sim_data.mTNumerics.m_res_tol_mixed = 1.e-8;
  sim_data.mTNumerics.m_corr_tol_mixed = 1.e-8;
  sim_data.mTNumerics.m_four_approx_spaces_Q = true;
  sim_data.mTNumerics.m_nThreadsMixedProblem = glob_n_threads;
  sim_data.mTNumerics.m_max_iter_sfi = 10;
  sim_data.mTNumerics.m_max_iter_mixed = 10;
  sim_data.mTNumerics.m_max_iter_transport = 10;

  sim_data.mTPostProcess.m_file_name_mixed = "postdarcy.vtk";
  sim_data.mTPostProcess.m_file_name_transport = "posttransport.vtk";
  TPZStack<std::string, 10> scalnames, vecnames, scalnamesTransport;
  vecnames.Push("Flux");
  scalnames.Push("Pressure");
  scalnames.Push("div_q");
  if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
    scalnames.Push("g_average");
    scalnames.Push("p_average");
  }
  scalnamesTransport.Push("Sw");
  scalnamesTransport.Push("So");

  sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
  sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
  sim_data.mTPostProcess.m_scalnamesTransport = scalnamesTransport;

  int n_steps = sim_data.mTNumerics.m_n_steps;
  sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
  REAL dt = sim_data.mTNumerics.m_dt;
  TPZStack<REAL, 100> reporting_times;
  REAL time = sim_data.mTPostProcess.m_file_time_step;
  int n_reporting_times = (n_steps) / (time * 1 / dt) + 1;
  REAL r_time = 0.0;
  int j = 1;
  for (int i = 1; i <= n_reporting_times; i++) {

    r_time = j * dt * (time / dt);
    reporting_times.push_back(r_time);
    j += 1;
  }
  sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
}
