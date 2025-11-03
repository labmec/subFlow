//
//  Created by Giovane Avancini on 05/09/25.
//

#pragma once

#include "Material/TPZMatTypes.h"
#include "TSFSavable.h"
#include "pzmanvector.h"
#include "pzstack.h"
#include <map>
#include <stdio.h>
#include <tuple>

/// Object that stores all data required to set up a reservoir problem
class TSFProblemData : public TSFSavable {
public:
  /// Default constructor
  TSFProblemData();

  std::string mSimulationName = "";

  /// Copy constructor
  TSFProblemData(const TSFProblemData &other);

  // Copy assignment operator
  TSFProblemData &operator=(const TSFProblemData &other);

  /// Destructor
  ~TSFProblemData();

  /// Write object state
  void Write(TPZStream &buf, int withclassid) const;

  /// Read object state
  void Read(TPZStream &buf, void *context);

  /// Read object state
  virtual int ClassId() const;

  /**
   * @brief Class that stores geometric information
   */
  class TGeometry : public TSFSavable {
  public:
    /** @brief
     Contains the  name and material for elements of a certain dimension.
     */
    std::map<std::string, int> mDomainNameAndMatId;

    /** @brief
     MaterialID of the interface element that will be inserted in the transport mesh
     */
    int mInterface_material_id = 100; // Interface material Volumetric-Volumetric

    // definition of the standard material ids
    /// material id for the Macro fluxes
    int m_skeletonMatId = 19;
    /// material id for H(div) wrap materials
    int m_HdivWrapMatId = 15;
    /// material id for the mortar pressure space
    int m_MortarMatId = 20;
    /// material id for a positive lagrange multiplier interface
    int m_posLagrangeMatId = 30;
    /// material id for a negative lagrange multiplier interface
    int m_negLagrangeMatId = 35;
    /// material id for a zero order H(div) boundary flux
    int m_zeroOrderHdivFluxMatId = 40;
    /// material id for pressure lagrange multiplier (same as mortar because they should not be in the same mesh)
    int m_pressureMatId = 20;

    std::string mGmeshFileName = "";

    /** @brief Default constructor */
    TGeometry() {}

    /** @brief Destructor */
    ~TGeometry() {}

    /** @brief Copy constructor */
    TGeometry(const TGeometry &other) {
      mDomainNameAndMatId = other.mDomainNameAndMatId;
      mInterface_material_id = other.mInterface_material_id;

      m_pressureMatId = other.m_pressureMatId;
      /// material id for the Macro fluxes
      m_skeletonMatId = other.m_skeletonMatId;
      /// material id for H(div) wrap materials
      m_HdivWrapMatId = other.m_HdivWrapMatId;
      /// material id for the mortar pressure space
      m_MortarMatId = other.m_MortarMatId;
      /// material id for a positive lagrange multiplier interface
      m_posLagrangeMatId = other.m_posLagrangeMatId;
      /// material id for a negative lagrange multiplier interface
      m_negLagrangeMatId = other.m_negLagrangeMatId;
      /// material id for a zero order H(div) boundary flux
      m_zeroOrderHdivFluxMatId = other.m_zeroOrderHdivFluxMatId;

      mGmeshFileName = other.mGmeshFileName;
    }
    /** @brief Copy assignment operator*/
    TGeometry &operator=(const TGeometry &other) {
      if (this != &other) // prevent self-assignment
      {
        mDomainNameAndMatId = other.mDomainNameAndMatId;
        mInterface_material_id = other.mInterface_material_id;
        m_pressureMatId = other.m_pressureMatId;
        /// material id for the Macro fluxes
        m_skeletonMatId = other.m_skeletonMatId;
        /// material id for H(div) wrap materials
        m_HdivWrapMatId = other.m_HdivWrapMatId;
        /// material id for the mortar pressure space
        m_MortarMatId = other.m_MortarMatId;
        /// material id for a positive lagrange multiplier interface
        m_posLagrangeMatId = other.m_posLagrangeMatId;
        /// material id for a negative lagrange multiplier interface
        m_negLagrangeMatId = other.m_negLagrangeMatId;
        /// material id for a zero order H(div) boundary flux
        m_zeroOrderHdivFluxMatId = other.m_zeroOrderHdivFluxMatId;

        mGmeshFileName = other.mGmeshFileName;
      }
      return *this;
    }
  };

  /**
   * @brief Class that stores PetroPhysics information, i.e. Interaction between fluid and rock.
   */

  class TPetroPhysics : public TSFSavable {
  public:
    REAL mOilViscosity;
    REAL mWaterViscosity;
    REAL mSwr;    // residual water saturation
    REAL mSor;    // residual oil saturation
    int mKrModel; // 0 - Linear, 1 - Quadratic, 2 - Quadratic with residual

    std::vector<std::function<std::tuple<REAL, REAL>(REAL &)>> mKro;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &)>> mKrw;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &, REAL &)>> mFo;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &, REAL &)>> mFw;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &)>> mLambdaW;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &)>> mLambdaO;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &, REAL &)>> mLambdaTotal;

    /** @brief Default constructor */
    TPetroPhysics() {
      mOilViscosity = 1.0;
      mWaterViscosity = 1.0;
      REAL mSwr = 0.0;
      REAL mSor = 0.0;
      int mKrModel = 0;
      mKro.resize(3);
      mKrw.resize(3);
      mFo.resize(3);
      mFw.resize(3);
      mLambdaW.resize(3);
      mLambdaO.resize(3);
      mLambdaTotal.resize(3);
      CreateLinearKrModel();
      CreateQuadraticKrModel();
      CreateQuadraticResidualKrModel();
    }

    /** @brief Destructor */
    ~TPetroPhysics() {
    }

    /** @brief Copy constructor */
    TPetroPhysics(const TPetroPhysics &other) {
      mOilViscosity = other.mOilViscosity;
      mWaterViscosity = other.mWaterViscosity;
      mSwr = other.mSwr;
      mSor = other.mSor;
      mKrModel = other.mKrModel;
      mKro = other.mKro;
      mKrw = other.mKrw;
      mFo = other.mFo;
      mFw = other.mFw;
      mLambdaW = other.mLambdaW;
      mLambdaO = other.mLambdaO;
      mLambdaTotal = other.mLambdaTotal;
    }

    /** @brief Copy assignment operator*/
    TPetroPhysics &operator=(const TPetroPhysics &other) {
      if (this != &other) // prevent self-assignment
      {
        mOilViscosity = other.mOilViscosity;
        mWaterViscosity = other.mWaterViscosity;
        mSwr = other.mSwr;
        mSor = other.mSor;
        mKrModel = other.mKrModel;
        mKro = other.mKro;
        mKrw = other.mKrw;
        mFo = other.mFo;
        mFw = other.mFw;
        mLambdaW = other.mLambdaW;
        mLambdaO = other.mLambdaO;
        mLambdaTotal = other.mLambdaTotal;
      }
      return *this;
    }
    void CreateLinearKrModel();
    void CreateQuadraticKrModel();
    void CreateQuadraticResidualKrModel();
    void UpdateLambdasAndFracFlows(int krModel);
  };

  /**
   * @brief Class that stores Fluid properties. For instance density of water and oil, compressibility of water and oil.
   */
  class TFluidProperties : public TSFSavable {
  public:
    REAL mOilViscosity;
    REAL mWaterViscosity;
    REAL mWaterCompressibility;
    REAL mOilCompressibility;
    REAL mOilDensityRef;
    REAL mWaterDensityRef;
    REAL mReferencePressure;
    std::function<std::tuple<REAL, REAL>(REAL &)> mOilDensityF;
    std::function<std::tuple<REAL, REAL>(REAL &)> mWaterDensityF;

    /** @brief Default constructor */
    TFluidProperties() {
      mOilViscosity = 1.0;
      mWaterViscosity = 1.0;
      mOilDensityRef = 1000;
      mWaterDensityRef = 1000;
      mWaterCompressibility = 1.0e-8;
      mOilCompressibility = 1.0e-7;
      mReferencePressure = 1.013e5;
      mOilDensityF = 0;
      mWaterDensityF = 0;
      CreateLinearDensityFunction();
    }

    /** @brief Destructor */
    ~TFluidProperties() {
    }

    /** @brief Copy constructor */
    TFluidProperties(const TFluidProperties &other) {
      mOilViscosity = other.mOilViscosity;
      mWaterViscosity = other.mWaterViscosity;
      mWaterCompressibility = other.mWaterCompressibility;
      mOilCompressibility = other.mOilCompressibility;

      mReferencePressure = other.mReferencePressure;
      mOilDensityRef = other.mOilDensityRef;
      mWaterDensityRef = other.mWaterDensityRef;
      mOilDensityF = other.mOilDensityF;
      mWaterDensityF = other.mWaterDensityF;
    }

    /** @brief Copy assignment operator*/
    TFluidProperties &operator=(const TFluidProperties &other) {
      mOilViscosity = other.mOilViscosity;
      mWaterViscosity = other.mWaterViscosity;
      mWaterCompressibility = other.mWaterCompressibility;
      mOilCompressibility = other.mOilCompressibility;
      mOilDensityRef = other.mOilDensityRef;
      mWaterDensityRef = other.mWaterDensityRef;
      mReferencePressure = other.mReferencePressure;
      mOilDensityF = other.mOilDensityF;
      mWaterDensityF = other.mWaterDensityF;
      return *this;
    }
    void CreateLinearDensityFunction();
    void CreateExponentialDensityFunction();
  };

  /**
   * @brief Class that stores Reservoir rock properties. For instance permeability and porosity of the rock
   */
  class TReservoirProperties : public TSFSavable {
  public:
    // associates with a material id the porosity and permeability
    std::map<int, std::pair<REAL, REAL>> mPorosityAndPermeability;

    // Function that given a point (x,y,z) returns s0 (initial saturation)
    std::function<REAL(const TPZVec<REAL> &)> s0;

    /** @brief Default constructor */
    TReservoirProperties() {
      mPorosityAndPermeability[0] = std::make_pair(1.0, 1.0);
    }

    /** @brief Destructor */
    ~TReservoirProperties() {
    }

    /** @brief Copy constructor */
    TReservoirProperties(const TReservoirProperties &other) {
      mPorosityAndPermeability = other.mPorosityAndPermeability;
      s0 = other.s0;
    }

    /** @brief Copy assignment operator*/
    TReservoirProperties &operator=(const TReservoirProperties &other) {
      mPorosityAndPermeability = other.mPorosityAndPermeability;
      s0 = other.s0;
      return *this;
    }
  };

  /**
   * @brief Class that stores the boundary conditions of the problem
   */
  class TBoundaryConditions : public TSFSavable {
  public:
    std::map<std::string, int> mDomainNameAndMatId;
    /**
     * @brief Contains the boundary conditions for Darcy problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> mBCDarcyMatIdToTypeValue;

    /**
     * @brief Contains the boundary conditions for Transport problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> mBCTransportMatIdToTypeValue;

    /**
     * @brief Contains the functions to be applied on the boundaries in a map. Key = matidOfBC, value = functionID (0 if no function is prescribed)
     */
    std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> mBCDarcyMatIdToFunctionId;

    /**
     * @brief Contains the saturation functions to be applied on the boundaries in a map. Key = matidOfBC, value = functionID (0 if no function is prescribed)
     */
    std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> mBCTransportMatIdToFunctionId;

    /** @brief Default constructor */
    TBoundaryConditions() {}

    /** @brief Destructor */
    ~TBoundaryConditions() {
    }

    /** @brief Copy constructor */
    TBoundaryConditions(const TBoundaryConditions &other) {
      mBCDarcyMatIdToTypeValue = other.mBCDarcyMatIdToTypeValue;
      mBCTransportMatIdToTypeValue = other.mBCTransportMatIdToTypeValue;
      mDomainNameAndMatId = other.mDomainNameAndMatId;
      mBCDarcyMatIdToFunctionId = other.mBCDarcyMatIdToFunctionId;
      mBCTransportMatIdToFunctionId = other.mBCTransportMatIdToFunctionId;
    }

    /** @brief Copy assignment operator*/
    TBoundaryConditions &operator=(const TBoundaryConditions &other) {
      if (this != &other) // prevent self-assignment
      {
        mBCDarcyMatIdToTypeValue = other.mBCDarcyMatIdToTypeValue;
        mBCTransportMatIdToTypeValue = other.mBCTransportMatIdToTypeValue;
        mDomainNameAndMatId = other.mDomainNameAndMatId;
        mBCDarcyMatIdToFunctionId = other.mBCDarcyMatIdToFunctionId;
        mBCTransportMatIdToFunctionId = other.mBCTransportMatIdToFunctionId;
      }
      return *this;
    }
  };

  /**
   * @brief Class that stores the numerical parameters of the simulation
   * stores the size of the timestep
   * stores the tolerances of the iterative solver
   * stores the maximum number of iterations
   * stores the type of approximation spaces to be used
   * stores the number of threads to be used in paralel computing
   * defines polynomial orders when applying hybridization
   * defines if the mesh is going to be substructured
   */
  class TNumerics : public TSFSavable {
  public:
    /**
     * @brief time step size
     */
    REAL m_dt;

    /**
     * @brief Flag to indicate the analysis type: 1 - only Darcy, 2 - only transport, 3 - coupled
     */
    int fAnalysisType;

    /**
     * @brief Residual tolerance for Darcy
     */
    REAL m_res_tol_darcy;

    /**
     * @brief Residual tolerance for transport
     */
    REAL m_res_tol_transport;

    /**
     * @brief Correction tolerance for darcy
     */
    REAL m_corr_tol_darcy;

    /**
     * @brief Correction tolerance for transport
     */
    REAL m_corr_tol_transport;

    REAL m_sfi_tol;
    /**
     * @brief Maximum number of iterations per time step for darcy
     */
    int m_max_iter_darcy;

    /**
     * @brief Maximum number of iterations per time step for transport
     */
    int m_max_iter_transport;

    /**
     * @brief Maximum number of Sequential Fully Implicit (SFI) iterations per time step
     */
    int m_max_iter_sfi;

    /**
     * @brief Number of time steps
     */
    int m_n_steps;

    // Order of approximation for the border of the Pressure element
    int m_MortarBorderElementPresOrder;
    // Order of approximation for the border of the Flux element
    int m_MortarBorderElementFluxOrder;
    /**
     * @brief Directive for the use of four spaces
     */
    bool m_four_approx_spaces;

    /**
     * @brief Axisymmetric flag
     */
    bool m_is_axisymmetric;

    /**
     * @brief Linear trace flag. If true, the darcy problem is ran only once
     */
    bool m_is_linearTrace;

    /**
     * @brief Approximation space "type"
     */
    enum MSpaceType { ENone,
                      E2Space,
                      E4Space,
                      E2SpaceMHM,
                      E4SpaceMHM,
                      E4SpaceMortar,
                      E4Space1Hybridization };

    MSpaceType m_SpaceType = ENone;

    std::vector<REAL> m_gravity;

    int m_nThreadsDarcyProblem = 0;

    bool m_UseSubstructures_Q;
    /** @brief Default constructor */
    TNumerics() {
      m_dt = 0.0;
      fAnalysisType = 3;
      m_res_tol_darcy = 1.0e-4;
      m_res_tol_transport = 1.0e-7;
      m_corr_tol_darcy = 1.0e-4;
      m_corr_tol_transport = 1.0e-7;
      m_max_iter_darcy = 0;
      m_max_iter_transport = 0;
      m_max_iter_sfi = 0;
      m_sfi_tol = 1.0e-3;
      m_n_steps = 0;
      m_four_approx_spaces = false;
      m_is_axisymmetric = false;
      m_is_linearTrace = true;
      m_SpaceType = ENone;
      m_gravity.resize(3, 0.0);
      m_nThreadsDarcyProblem = 0;
      m_MortarBorderElementPresOrder = 0;
      m_MortarBorderElementFluxOrder = 0;
      m_UseSubstructures_Q = false;
    }
    /** @brief Destructor */
    ~TNumerics() {
    }

    /** @brief Copy constructor */
    TNumerics(const TNumerics &other) {
      m_dt = other.m_dt;
      fAnalysisType = other.fAnalysisType;
      m_res_tol_darcy = other.m_res_tol_darcy;
      m_res_tol_transport = other.m_res_tol_transport;
      m_corr_tol_darcy = other.m_corr_tol_darcy;
      m_corr_tol_transport = other.m_corr_tol_transport;
      m_max_iter_darcy = other.m_max_iter_darcy;
      m_max_iter_transport = other.m_max_iter_transport;
      m_max_iter_sfi = other.m_max_iter_sfi;
      m_sfi_tol = other.m_sfi_tol;
      m_n_steps = other.m_n_steps;
      m_four_approx_spaces = other.m_four_approx_spaces;
      m_is_axisymmetric = other.m_is_axisymmetric;
      m_is_linearTrace = other.m_is_linearTrace;
      m_SpaceType = other.m_SpaceType;
      m_gravity = other.m_gravity;
      m_nThreadsDarcyProblem = other.m_nThreadsDarcyProblem;
      m_MortarBorderElementPresOrder = other.m_MortarBorderElementPresOrder;
      m_MortarBorderElementFluxOrder = other.m_MortarBorderElementFluxOrder;
      m_UseSubstructures_Q = other.m_UseSubstructures_Q;
    }

    /** @brief Copy assignment operator*/
    TNumerics &operator=(const TNumerics &other) {
      // check for self-assignment
      if (&other == this) {
        return *this;
      }

      m_dt = other.m_dt;
      fAnalysisType = other.fAnalysisType;
      m_res_tol_darcy = other.m_res_tol_darcy;
      m_res_tol_transport = other.m_res_tol_transport;
      m_corr_tol_darcy = other.m_corr_tol_darcy;
      m_corr_tol_transport = other.m_corr_tol_transport;
      m_max_iter_darcy = other.m_max_iter_darcy;
      m_max_iter_transport = other.m_max_iter_transport;
      m_max_iter_sfi = other.m_max_iter_sfi;
      m_sfi_tol = other.m_sfi_tol;
      m_n_steps = other.m_n_steps;
      m_four_approx_spaces = other.m_four_approx_spaces;
      m_is_axisymmetric = other.m_is_axisymmetric;
      m_is_linearTrace = other.m_is_linearTrace;
      m_SpaceType = other.m_SpaceType;
      m_gravity = other.m_gravity;
      m_nThreadsDarcyProblem = other.m_nThreadsDarcyProblem;
      m_MortarBorderElementPresOrder = other.m_MortarBorderElementPresOrder;
      m_MortarBorderElementFluxOrder = other.m_MortarBorderElementFluxOrder;
      m_UseSubstructures_Q = other.m_UseSubstructures_Q;
      return *this;
    }

    bool operator==(const TNumerics &other) {
      // check for self-assignment
      if (&other == this) {
        return true;
      }

      return m_dt == other.m_dt &&
             fAnalysisType == other.fAnalysisType &&
             m_res_tol_darcy == other.m_res_tol_darcy &&
             m_res_tol_transport == other.m_res_tol_transport &&
             m_corr_tol_darcy == other.m_corr_tol_darcy &&
             m_corr_tol_transport == other.m_corr_tol_transport &&
             m_max_iter_darcy == other.m_max_iter_darcy &&
             m_max_iter_transport == other.m_max_iter_transport &&
             m_max_iter_sfi == other.m_max_iter_sfi &&
             m_sfi_tol == other.m_sfi_tol &&
             m_n_steps == other.m_n_steps &&
             m_four_approx_spaces == other.m_four_approx_spaces &&
             m_is_axisymmetric == other.m_is_axisymmetric &&
             m_is_linearTrace == other.m_is_linearTrace &&
             m_SpaceType == other.m_SpaceType &&
             m_gravity == other.m_gravity &&
             m_nThreadsDarcyProblem == other.m_nThreadsDarcyProblem &&
             m_MortarBorderElementPresOrder == other.m_MortarBorderElementPresOrder &&
             m_MortarBorderElementFluxOrder == other.m_MortarBorderElementFluxOrder &&
             m_UseSubstructures_Q == other.m_UseSubstructures_Q;
    }

    void Write(TPZStream &buf, int withclassid) const { // ok
      buf.Write(&m_dt);
      buf.Write(&fAnalysisType);
      buf.Write(&m_res_tol_darcy);
      buf.Write(&m_res_tol_transport);
      buf.Write(&m_corr_tol_darcy);
      buf.Write(&m_max_iter_darcy);
      buf.Write(&m_max_iter_transport);
      buf.Write(&m_max_iter_sfi);
      buf.Write(&m_sfi_tol);
      buf.Write(&m_n_steps);
      int temp = m_four_approx_spaces;
      buf.Write(&temp);
      temp = m_is_axisymmetric;
      buf.Write(&temp);
      temp = m_is_linearTrace;
      buf.Write(&temp);
      temp = m_SpaceType;
      buf.Write(&temp);
      buf.Write(m_gravity);
      buf.Write(m_nThreadsDarcyProblem);
      buf.Write(m_UseSubstructures_Q);
    }

    void Read(TPZStream &buf, void *context) { // ok
      buf.Read(&m_dt);
      buf.Read(&fAnalysisType);
      buf.Read(&m_res_tol_darcy);
      buf.Read(&m_res_tol_transport);
      buf.Read(&m_corr_tol_darcy);
      buf.Read(&m_corr_tol_transport);
      buf.Read(&m_max_iter_darcy);
      buf.Read(&m_max_iter_transport);
      buf.Read(&m_max_iter_sfi);
      buf.Read(&m_sfi_tol);
      buf.Read(&m_n_steps);
      int temp;
      buf.Read(&temp);
      m_four_approx_spaces = temp;
      buf.Read(&temp);
      m_is_axisymmetric = temp;
      buf.Read(&temp);
      m_is_linearTrace = temp;
      buf.Read(&temp);
      m_SpaceType = (MSpaceType)temp;
      buf.Read(&m_nThreadsDarcyProblem);
      buf.Read(&m_UseSubstructures_Q);
    }

    virtual int ClassId() const {
      return Hash("TSFProblemData::TNumerics");
    }

    void Print() const {
      std::cout << m_dt << std::endl;
      std::cout << fAnalysisType << std::endl;
      std::cout << m_res_tol_darcy << std::endl;
      std::cout << m_res_tol_transport << std::endl;
      std::cout << m_corr_tol_darcy << std::endl;
      std::cout << m_corr_tol_transport << std::endl;
      std::cout << m_max_iter_darcy << std::endl;
      std::cout << m_max_iter_transport << std::endl;
      std::cout << m_max_iter_sfi << std::endl;
      std::cout << m_n_steps << std::endl;
      std::cout << m_four_approx_spaces << std::endl;
      std::cout << m_is_axisymmetric << std::endl;
      std::cout << m_is_linearTrace << std::endl;
      std::cout << m_SpaceType << std::endl;
      std::cout << m_nThreadsDarcyProblem << std::endl;
    }
  };

  /**
   * @brief Class that stores the PostProcess information
   */
  class TPostProcess : public TSFSavable {
  public:
    /**
     * @brief Mixed operator vtk file name
     */
    std::string m_file_name_darcy;

    /**
     * @brief Transpor operator vtk file name
     */
    std::string m_file_name_transport;

    /**
     * @brief Contains scalar variables that will be postprocessed
     */
    TPZStack<std::string, 10> m_scalnamesDarcy;

    /**
     * @brief Contains scalar variables that will be postprocessed
     */
    TPZStack<std::string, 10> m_scalnamesTransport;

    /**
     * @brief Contains vector variables that will be postprocessed
     */
    TPZStack<std::string, 10> m_vecnamesDarcy;

    /**
     * @brief Period of time post-processed data is printed
     */
    REAL m_file_time_step;

    /**
     * @brief Contains the times at which post-processed data is printed
     */
    TPZStack<REAL, 100> m_vec_reporting_times;

    /**
     * @brief Default constructor
     */
    TPostProcess() {
      m_file_name_darcy = "";
      m_file_name_transport = "";
      m_scalnamesDarcy.Resize(0);
      m_scalnamesDarcy.Resize(0);
      m_scalnamesTransport.Resize(2);
      m_scalnamesTransport.Push("Sw");
      m_scalnamesTransport.Push("So");
      m_file_time_step = 0.0;
      m_vec_reporting_times.Resize(0);
    }
    /**
     * @brief Destructor
     */
    ~TPostProcess() {
    }

    /**
     * @brief Copy constructor
     */
    TPostProcess(const TPostProcess &other) {
      m_file_name_darcy = other.m_file_name_darcy;
      m_file_name_transport = other.m_file_name_transport;
      m_scalnamesDarcy = other.m_scalnamesDarcy;
      m_scalnamesDarcy = other.m_scalnamesDarcy;
      m_scalnamesTransport = other.m_scalnamesTransport;
      m_file_time_step = other.m_file_time_step;
      m_vec_reporting_times = other.m_vec_reporting_times;
    }
    /**
     * @brief Copy assignment operator
     */
    TPostProcess &operator=(const TPostProcess &other) {
      // check for self-assignment
      if (&other == this) {
        return *this;
      }

      m_file_name_darcy = other.m_file_name_darcy;
      m_file_name_transport = other.m_file_name_transport;
      m_vecnamesDarcy = other.m_vecnamesDarcy;
      m_scalnamesDarcy = other.m_scalnamesDarcy;
      m_scalnamesTransport = other.m_scalnamesTransport;
      m_file_time_step = other.m_file_time_step;
      m_vec_reporting_times = other.m_vec_reporting_times;

      return *this;
    }

    bool operator==(const TPostProcess &other) {
      // check for self-assignment
      if (&other == this) {
        return true;
      }

      return m_file_name_darcy == other.m_file_name_darcy &&
             m_file_name_transport == other.m_file_name_transport &&
             m_vecnamesDarcy == other.m_vecnamesDarcy &&
             m_scalnamesDarcy == other.m_scalnamesDarcy &&
             m_scalnamesTransport == other.m_scalnamesTransport &&
             m_file_time_step == other.m_file_time_step &&
             m_vec_reporting_times == other.m_vec_reporting_times;
    }

    void Write(TPZStream &buf, int withclassid) const { // ok
      buf.Write(&m_file_name_darcy);
      buf.Write(&m_file_name_transport);
      buf.Write(m_vecnamesDarcy);
      buf.Write(m_scalnamesDarcy);
      buf.Write(m_scalnamesTransport);
      buf.Write(&m_file_time_step);
      buf.Write(m_vec_reporting_times);
    }

    void Read(TPZStream &buf, void *context) { // ok
      buf.Read(&m_file_name_darcy);
      buf.Read(&m_file_name_transport);
      buf.Read(m_vecnamesDarcy);
      buf.Read(m_scalnamesDarcy);
      buf.Read(m_scalnamesTransport);
      buf.Read(&m_file_time_step);
      buf.Read(m_vec_reporting_times);
    }

    virtual int ClassId() const {
      return Hash("TSFProblemData::TPostProcess");
    }

    void Print() const {
      std::cout << m_file_name_darcy << std::endl;
      std::cout << m_file_name_transport << std::endl;
      std::cout << m_file_time_step << std::endl;
      std::cout << m_vec_reporting_times << std::endl;
      // scalnames and vecnames
    }
  };

  // describes the material id of peripheral material objects (skeleton matid, hdivwrap matid, etc)
  // contains a data structure associating a matid with a string generated by gmesh
  TGeometry mTGeometry;
  // describes the properties of the fluid mixture (relative permeabilities)
  TPetroPhysics mTPetroPhysics;
  // describes the properties of the oil and water
  TFluidProperties mTFluidProperties;
  // describes the initial conditions of saturation
  // describes the values of permeability by material id
  // describes the porosity of the reservoir
  // describes the volume factor (?) by material id
  //
  TReservoirProperties mTReservoirProperties;

  // defines the matids and values associated with the boundary conditions
  TBoundaryConditions mTBoundaryConditions;
  /*
   * stores the size of the timestep
   * stores the tolerances of the iterative solver
   * stores the maximum number of iterations
   * stores the type of approximation spaces to be used
   * stores the number of threads to be used in paralel computing
   * defines polynomial orders when applying hybridization
   * defines if the mesh is going to be substructured
   */
  TNumerics mTNumerics;
  // define the post processing information
  TPostProcess mTPostProcess;
};
