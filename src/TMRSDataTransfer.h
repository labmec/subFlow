//
//  TMRSDataTransfer.hpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#ifndef TMRSDataTransfer_hpp
#define TMRSDataTransfer_hpp

#include <stdio.h>

#include <map>
#include <map>    // for map
#include <tuple>  // for tuple
#include <vector>

#include "DFNPolygon.h"
#include "TMRSSavable.h"
#include "TRSLinearInterpolator.h"
#include "pzmanvector.h"
#include "Material/TPZMatTypes.h"

/// Object that stores all data required to set up a reservoir problem
class TMRSDataTransfer : public TMRSSavable {
 public:
  /// Default constructor
  TMRSDataTransfer();

  std::string mSimulationName = "";

  /// Copy constructor
  TMRSDataTransfer(const TMRSDataTransfer &other);

  // Copy assignment operator
  TMRSDataTransfer &operator=(const TMRSDataTransfer &other);

  /// Destructor
  ~TMRSDataTransfer();

  /// Write object state
  void Write(TPZStream &buf, int withclassid) const;

  /// Read object state
  void Read(TPZStream &buf, void *context);

  /// Read object state
  virtual int ClassId() const;

  /**
   * @brief Class that stores geometric information
   */
  class TGeometry : public TMRSSavable {
   public:
    /** @brief
     Contains the  name and material for elements of a certain dimension.
     */
    std::map<std::string, int> mDomainNameAndMatId;

    /** @brief
     Contains the material id of the fractures by dimension.
     */
    std::map<std::string, int> mDomainFracNameAndMatId;
    std::map<std::string, int> mDomainFracIntersectionNameAndMatId;

    /** @brief
     MaterialID of the interface element that will be inserted in the transport mesh
     */
    int mInterface_material_id = 100;           // Interface material Volumetric-Volumetric
    int mInterface_material_idFracSup = 101;    // Interface material Fracture-Volumetric1
    int mInterface_material_idFracInf = 102;    // Interface material Fracture-Volumetric2
    int mInterface_material_idFracFrac = 103;   // Interface material Fracture-Fracture
    int mInterface_material_idFracBound = 104;  // Interface material Fracture-Boundary

    /// number of times the skeleton elements should be refined
    int mSkeletonDiv = 0;
    /// unused variable in this project space
    int mnLayers = 1;
    /// number of times the global mesh should be divided
    int mnref = 0;

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
    ~TGeometry() {
    }

    /** @brief Copy constructor */
    TGeometry(const TGeometry &other) {
      mDomainNameAndMatId = other.mDomainNameAndMatId;
      mDomainFracNameAndMatId = other.mDomainFracNameAndMatId;
      mDomainFracIntersectionNameAndMatId = other.mDomainFracIntersectionNameAndMatId;
      mInterface_material_id = other.mInterface_material_id;
      mInterface_material_idFracSup = other.mInterface_material_idFracSup;
      mInterface_material_idFracInf = other.mInterface_material_idFracInf;
      mInterface_material_idFracFrac = other.mInterface_material_idFracFrac;
      mInterface_material_idFracBound = other.mInterface_material_idFracBound;

      mSkeletonDiv = other.mSkeletonDiv;
      m_pressureMatId = other.m_pressureMatId;
      mnLayers = other.mnLayers;
      mnref = other.mnref;
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
      if (this != &other)  // prevent self-assignment
      {
        mDomainNameAndMatId = other.mDomainNameAndMatId;
        mDomainFracNameAndMatId = other.mDomainFracNameAndMatId;
        mDomainFracIntersectionNameAndMatId = other.mDomainFracIntersectionNameAndMatId;
        mInterface_material_id = other.mInterface_material_id;
        mInterface_material_idFracSup = other.mInterface_material_idFracSup;
        mInterface_material_idFracInf = other.mInterface_material_idFracInf;
        mInterface_material_idFracFrac = other.mInterface_material_idFracFrac;
        mInterface_material_idFracBound = other.mInterface_material_idFracBound;
        mSkeletonDiv = other.mSkeletonDiv;
        m_pressureMatId = other.m_pressureMatId;
        mnLayers = other.mnLayers;
        mnref = other.mnref;
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

    const bool isThereFracture() const { return mDomainFracNameAndMatId.size(); }
  };

  /**
   * @brief Class that stores PetroPhysics information, i.e. Interaction between fluid and rock.
   */

  class TPetroPhysics : public TMRSSavable {
   public:
    REAL mOilViscosity;
    REAL mWaterViscosity;
    REAL mSwr; // residual water saturation
    REAL mSor; // residual oil saturation
    int mKrModel; // 0 - Linear, 1 - Quadratic, 2 - Quadratic with residual

    /** @brief Contains the water relative permeability model for each layer */
    std::vector<TRSLinearInterpolator> mLayer_Krw_RelPerModel;
    /** @brief Contains the oil relative permeability model for each layer */
    std::vector<TRSLinearInterpolator> mLayer_Kro_RelPerModel;

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
      mLayer_Krw_RelPerModel.clear();
      mLayer_Kro_RelPerModel.clear();
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
      mLayer_Krw_RelPerModel = other.mLayer_Krw_RelPerModel;
      mLayer_Kro_RelPerModel = other.mLayer_Kro_RelPerModel;
    }

    /** @brief Copy assignment operator*/
    TPetroPhysics &operator=(const TPetroPhysics &other) {
      if (this != &other)  // prevent self-assignment
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
        mLayer_Krw_RelPerModel = other.mLayer_Krw_RelPerModel;
        mLayer_Kro_RelPerModel = other.mLayer_Kro_RelPerModel;
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
  class TFluidProperties : public TMRSSavable {
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
  class TReservoirProperties : public TMRSSavable {
   public:
    // associates with a material id the porosity and volume scale (volume scale is a constant, 1 for volume elements and lenght (maybe width?) for fractures)
    // the volume factor is the property of a fluid and the porosity is a property of the reservoir
    // why are these properties not defined by a map?
    std::vector<std::tuple<int, REAL, REAL>> mPorosityAndVolumeScale;
    bool fPropsFromPreProcess = false;

    // Function that given a point (x,y,z) returns kx, ky, kz and porosity
    // This is used in the transport problem for very heterogeneous cases.
    std::function<std::vector<REAL>(const TPZVec<REAL> &)> kappa_phi;

    std::map<int, REAL> m_permeabilitiesbyId;

    // Function that given a point (x,y,z) returns s0 (initial saturation)
    std::function<REAL(const TPZVec<REAL> &)> s0;
    std::string mPropsFileName = "";
    /** @brief Default constructor */
    TReservoirProperties() {
      mPorosityAndVolumeScale.resize(1);
      mPorosityAndVolumeScale[0] = std::make_tuple(1, 1.0, 0.1);
      fPropsFromPreProcess = false;
    }

    /** @brief Destructor */
    ~TReservoirProperties() {
    }

    /** @brief Copy constructor */
    TReservoirProperties(const TReservoirProperties &other) {
      mPorosityAndVolumeScale = other.mPorosityAndVolumeScale;
      fPropsFromPreProcess = other.fPropsFromPreProcess;
      kappa_phi = other.kappa_phi;
      s0 = other.s0;
      mPropsFileName = other.mPropsFileName;
      m_permeabilitiesbyId = other.m_permeabilitiesbyId;
    }

    /** @brief Copy assignment operator*/
    TReservoirProperties &operator=(const TReservoirProperties &other) {
      mPorosityAndVolumeScale = other.mPorosityAndVolumeScale;
      fPropsFromPreProcess = other.fPropsFromPreProcess;
      kappa_phi = other.kappa_phi;
      s0 = other.s0;
      mPropsFileName = other.mPropsFileName;
      m_permeabilitiesbyId = other.m_permeabilitiesbyId;
      return *this;
    }
  };

  /**
   * @brief Class that stores the boundary conditions of the problem
   */
  class TBoundaryConditions : public TMRSSavable {
   public:


    std::map<std::string, int> mDomainNameAndMatId;
    /**
     * @brief Contains the boundary conditions for Flow problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> mBCFlowMatIdToTypeValue;

    /**
     * @brief Contains the fracture boundary conditions for Flow problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> mBCFlowFracMatIdToTypeValue;

    /**
     * @brief Contains the boundary conditions for Transport problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> mBCTransportMatIdToTypeValue;

    /**
     * @brief Contains the functions to be applied on the boundaries in a map. Key = matidOfBC, value = functionID (0 if no function is prescribed)
     */
    std::map<int, std::pair<int,ForcingFunctionBCType<REAL>>> mBCFlowMatIdToFunctionId;

    /**
     * @brief Contains the saturation functions to be applied on the boundaries in a map. Key = matidOfBC, value = functionID (0 if no function is prescribed)
     */
    std::map<int, std::pair<int,ForcingFunctionBCType<REAL>>> mBCTransportMatIdToFunctionId;

    /** @brief Default constructor */
    TBoundaryConditions() {}

    /** @brief Destructor */
    ~TBoundaryConditions() {
    }

    /** @brief Copy constructor */
    TBoundaryConditions(const TBoundaryConditions &other) {
      mBCFlowMatIdToTypeValue = other.mBCFlowMatIdToTypeValue;
      mBCFlowFracMatIdToTypeValue = other.mBCFlowFracMatIdToTypeValue;
      mBCTransportMatIdToTypeValue = other.mBCTransportMatIdToTypeValue;
      mDomainNameAndMatId = other.mDomainNameAndMatId;
      mBCFlowMatIdToFunctionId = other.mBCFlowMatIdToFunctionId;
      mBCTransportMatIdToFunctionId = other.mBCTransportMatIdToFunctionId;
    }

    /** @brief Copy assignment operator*/
    TBoundaryConditions &operator=(const TBoundaryConditions &other) {
      if (this != &other)  // prevent self-assignment
      {
        mBCFlowMatIdToTypeValue = other.mBCFlowMatIdToTypeValue;
        mBCFlowFracMatIdToTypeValue = other.mBCFlowFracMatIdToTypeValue;
        mBCTransportMatIdToTypeValue = other.mBCTransportMatIdToTypeValue;
        mDomainNameAndMatId = other.mDomainNameAndMatId;
        mBCFlowMatIdToFunctionId = other.mBCFlowMatIdToFunctionId;
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
  class TNumerics : public TMRSSavable {
   public:
    /**
     * @brief time step size
     */
    REAL m_dt;

    /**
     * @brief if true, runs the transport problems. If not, just runs the flow/pressure problem
     */
    bool m_run_with_transport;

    /**
     * @brief Residual tolerance for mixed operator
     */
    REAL m_res_tol_mixed;

    /**
     * @brief Residual tolerance for transport operator
     */
    REAL m_res_tol_transport;

    /**
     * @brief Correction tolerance for mixed operator
     */
    REAL m_corr_tol_mixed;

    /**
     * @brief Correction tolerance for transport operator
     */
    REAL m_corr_tol_transport;

    REAL m_sfi_tol;
    /**
     * @brief Maximum number of iterations per time step for mixed operator
     */
    int m_max_iter_mixed;

    /**
     * @brief Maximum number of iterations per time step for transport operator
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
    bool m_four_approx_spaces_Q;

    /**
     * @brief Directive MHM mixed approximation
     */
    bool m_mhm_mixed_Q;

    /**
     * @brief Directive MHM mixed approximation
     */
    bool m_need_merge_meshes_Q;

    /**
     * @brief Axisymmetric flag
     */
    bool m_is_axisymmetric;

    /**
     * @brief Linear trace flag
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

    int m_nThreadsMixedProblem = 0;

    bool m_UseSubstructures_Q;
    /** @brief Default constructor */
    TNumerics() {
      m_dt = 0.0;
      m_run_with_transport = false;
      m_res_tol_mixed = 1.0e-4;
      m_res_tol_transport = 1.0e-7;
      m_corr_tol_mixed = 1.0e-4;
      m_corr_tol_transport = 1.0e-7;
      m_max_iter_mixed = 0;
      m_max_iter_transport = 0;
      m_max_iter_sfi = 0;
      m_sfi_tol = 1.0e-3;
      m_n_steps = 0;
      m_four_approx_spaces_Q = false;
      m_mhm_mixed_Q = false;
      m_is_axisymmetric = false;
      m_is_linearTrace = true;
      m_need_merge_meshes_Q = true;
      m_SpaceType = ENone;
      m_gravity.resize(3, 0.0);
      m_nThreadsMixedProblem = 0;
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
      m_run_with_transport = other.m_run_with_transport;
      m_res_tol_mixed = other.m_res_tol_mixed;
      m_res_tol_transport = other.m_res_tol_transport;
      m_corr_tol_mixed = other.m_corr_tol_mixed;
      m_corr_tol_transport = other.m_corr_tol_transport;
      m_max_iter_mixed = other.m_max_iter_mixed;
      m_max_iter_transport = other.m_max_iter_transport;
      m_max_iter_sfi = other.m_max_iter_sfi;
      m_sfi_tol = other.m_sfi_tol;
      m_n_steps = other.m_n_steps;
      m_four_approx_spaces_Q = other.m_four_approx_spaces_Q;
      m_is_axisymmetric = other.m_is_axisymmetric;
      m_is_linearTrace = other.m_is_linearTrace;
      m_mhm_mixed_Q = other.m_mhm_mixed_Q;
      m_need_merge_meshes_Q = other.m_need_merge_meshes_Q;
      m_SpaceType = other.m_SpaceType;
      m_gravity = other.m_gravity;
      m_nThreadsMixedProblem = other.m_nThreadsMixedProblem;
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
      m_run_with_transport = other.m_run_with_transport;
      m_res_tol_mixed = other.m_res_tol_mixed;
      m_res_tol_transport = other.m_res_tol_transport;
      m_corr_tol_mixed = other.m_corr_tol_mixed;
      m_corr_tol_transport = other.m_corr_tol_transport;
      m_max_iter_mixed = other.m_max_iter_mixed;
      m_max_iter_transport = other.m_max_iter_transport;
      m_max_iter_sfi = other.m_max_iter_sfi;
      m_sfi_tol = other.m_sfi_tol;
      m_n_steps = other.m_n_steps;
      m_four_approx_spaces_Q = other.m_four_approx_spaces_Q;
      m_is_axisymmetric = other.m_is_axisymmetric;
      m_is_linearTrace = other.m_is_linearTrace;
      m_mhm_mixed_Q = other.m_mhm_mixed_Q;
      m_need_merge_meshes_Q = other.m_need_merge_meshes_Q;
      m_SpaceType = other.m_SpaceType;
      m_gravity = other.m_gravity;
      m_nThreadsMixedProblem = other.m_nThreadsMixedProblem;
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
             m_run_with_transport == other.m_run_with_transport &&
             m_res_tol_mixed == other.m_res_tol_mixed &&
             m_res_tol_transport == other.m_res_tol_transport &&
             m_corr_tol_mixed == other.m_corr_tol_mixed &&
             m_corr_tol_transport == other.m_corr_tol_transport &&
             m_max_iter_mixed == other.m_max_iter_mixed &&
             m_max_iter_transport == other.m_max_iter_transport &&
             m_max_iter_sfi == other.m_max_iter_sfi &&
             m_sfi_tol == other.m_sfi_tol &&
             m_n_steps == other.m_n_steps &&
             m_four_approx_spaces_Q == other.m_four_approx_spaces_Q &&
             m_is_axisymmetric == other.m_is_axisymmetric &&
             m_is_linearTrace == other.m_is_linearTrace &&
             m_mhm_mixed_Q == other.m_mhm_mixed_Q &&
             m_need_merge_meshes_Q == other.m_need_merge_meshes_Q &&
             m_SpaceType == other.m_SpaceType &&
             m_gravity == other.m_gravity &&
             m_nThreadsMixedProblem == other.m_nThreadsMixedProblem &&
             m_MortarBorderElementPresOrder == other.m_MortarBorderElementPresOrder &&
             m_MortarBorderElementFluxOrder == other.m_MortarBorderElementFluxOrder &&
             m_UseSubstructures_Q == other.m_UseSubstructures_Q;
    }

    void Write(TPZStream &buf, int withclassid) const {  // ok
      buf.Write(&m_dt);
      buf.Write(m_run_with_transport);
      buf.Write(&m_res_tol_mixed);
      buf.Write(&m_res_tol_transport);
      buf.Write(&m_corr_tol_mixed);
      buf.Write(&m_max_iter_mixed);
      buf.Write(&m_max_iter_transport);
      buf.Write(&m_max_iter_sfi);
      buf.Write(&m_sfi_tol);
      buf.Write(&m_n_steps);
      int temp = m_four_approx_spaces_Q;
      buf.Write(&temp);
      temp = m_is_axisymmetric;
      buf.Write(&temp);
      temp = m_is_linearTrace;
      buf.Write(&temp);
      temp = m_mhm_mixed_Q;
      buf.Write(&temp);
      temp = m_need_merge_meshes_Q;
      buf.Write(&temp);
      temp = m_SpaceType;
      buf.Write(&temp);
      buf.Write(m_gravity);
      buf.Write(m_nThreadsMixedProblem);
      buf.Write(m_UseSubstructures_Q);
    }

    void Read(TPZStream &buf, void *context) {  // ok
      buf.Read(&m_dt);
      buf.Read(m_run_with_transport);
      buf.Read(&m_res_tol_mixed);
      buf.Read(&m_res_tol_transport);
      buf.Read(&m_corr_tol_mixed);
      buf.Read(&m_corr_tol_transport);
      buf.Read(&m_max_iter_mixed);
      buf.Read(&m_max_iter_transport);
      buf.Read(&m_max_iter_sfi);
      buf.Read(&m_sfi_tol);
      buf.Read(&m_n_steps);
      int temp;
      buf.Read(&temp);
      m_four_approx_spaces_Q = temp;
      buf.Read(&temp);
      m_is_axisymmetric = temp;
      buf.Read(&temp);
      m_is_linearTrace = temp;
      buf.Read(&temp);
      m_mhm_mixed_Q = temp;
      buf.Read(&temp);
      m_need_merge_meshes_Q = temp;
      buf.Read(&temp);
      m_SpaceType = (MSpaceType)temp;
      buf.Read(&m_nThreadsMixedProblem);
      buf.Read(m_UseSubstructures_Q);
    }

    virtual int ClassId() const {
      return Hash("TMRSDataTransfer::TNumerics");
    }

    void Print() const {
      std::cout << m_dt << std::endl;
      std::cout << m_run_with_transport << std::endl;
      std::cout << m_res_tol_mixed << std::endl;
      std::cout << m_res_tol_transport << std::endl;
      std::cout << m_corr_tol_mixed << std::endl;
      std::cout << m_corr_tol_transport << std::endl;
      std::cout << m_max_iter_mixed << std::endl;
      std::cout << m_max_iter_transport << std::endl;
      std::cout << m_max_iter_sfi << std::endl;
      std::cout << m_n_steps << std::endl;
      std::cout << m_four_approx_spaces_Q << std::endl;
      std::cout << m_is_axisymmetric << std::endl;
      std::cout << m_is_linearTrace << std::endl;
      std::cout << m_mhm_mixed_Q << std::endl;
      std::cout << m_need_merge_meshes_Q << std::endl;
      std::cout << m_SpaceType << std::endl;
      std::cout << m_nThreadsMixedProblem << std::endl;
    }
  };

  /**
   * @brief Class that stores the PostProcess information
   */
  class TPostProcess : public TMRSSavable {
   public:
    /**
     * @brief Mixed operator vtk file name
     */
    std::string m_file_name_mixed;

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
      m_file_name_mixed = "";
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
      m_file_name_mixed = other.m_file_name_mixed;
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

      m_file_name_mixed = other.m_file_name_mixed;
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

      return m_file_name_mixed == other.m_file_name_mixed &&
             m_file_name_transport == other.m_file_name_transport &&
             m_vecnamesDarcy == other.m_vecnamesDarcy &&
             m_scalnamesDarcy == other.m_scalnamesDarcy &&
             m_scalnamesTransport == other.m_scalnamesTransport &&
             m_file_time_step == other.m_file_time_step &&
             m_vec_reporting_times == other.m_vec_reporting_times;
    }

    void Write(TPZStream &buf, int withclassid) const {  // ok
      buf.Write(&m_file_name_mixed);
      buf.Write(&m_file_name_transport);
      buf.Write(m_vecnamesDarcy);
      buf.Write(m_scalnamesDarcy);
      buf.Write(m_scalnamesTransport);
      buf.Write(&m_file_time_step);
      buf.Write(m_vec_reporting_times);
    }

    void Read(TPZStream &buf, void *context) {  // ok
      buf.Read(&m_file_name_mixed);
      buf.Read(&m_file_name_transport);
      buf.Read(m_vecnamesDarcy);
      buf.Read(m_scalnamesDarcy);
      buf.Read(m_scalnamesTransport);
      buf.Read(&m_file_time_step);
      buf.Read(m_vec_reporting_times);
    }

    virtual int ClassId() const {
      return Hash("TMRSDataTransfer::TPostProcess");
    }

    void Print() const {
      std::cout << m_file_name_mixed << std::endl;
      std::cout << m_file_name_transport << std::endl;
      std::cout << m_file_time_step << std::endl;
      std::cout << m_vec_reporting_times << std::endl;
      // scalnames and vecnames
    }
  };

  class TFracProperties : public TMRSSavable {
   public:
    /// Structure with all fractures properties
    struct FracProp {
      //			int m_fracbc; // fracture boundary condition matid
      std::set<int> m_fracbc;    // fracture boundary condition matid
      int m_fracIntersectMatID;  // material id for intersections
      REAL m_perm;               // permeability of the fracture
      REAL m_width;              // fracture opening
      REAL m_porosity;           // porosity of the fracture
      DFNPolygon m_polydata;     // geometry of the fracture

      FracProp() : m_fracIntersectMatID(0), m_perm(0.), m_width(0.), m_porosity(0.) {
        m_fracbc.clear();
      }

      FracProp(const FracProp &copy) = default;

      FracProp &operator=(const FracProp &copy) = default;
    };

    /// Map <matid, FractureProperties>
    /// index is the fracture material id
    std::map<int, FracProp> m_fracprops;

    /**
     * @brief Default constructor
     */
    TFracProperties() {
    }
    /**
     * @brief Destructor
     */
    ~TFracProperties() {}

    /**
     * @brief Copy constructor
     */
    TFracProperties(const TFracProperties &other) {
      m_fracprops = other.m_fracprops;
    }
    /**
     * @brief Copy assignment operator
     */
    TFracProperties &operator=(const TFracProperties &other) {
      // check for self-assignment
      if (&other == this) {
        return *this;
      }

      m_fracprops = other.m_fracprops;
      return *this;
    }

    void Write(TPZStream &buf, int withclassid) const {
      buf.Write(&m_fracprops);
    }

    void Read(TPZStream &buf, void *context) {  // ok
      //            buf.Read(&m_fracprops);
    }

    virtual int ClassId() const {
      return Hash("TMRSDataTransfer::TFracProperties");
    }

    void Print() {
      std::cout << "NFractures= " << m_fracprops.size() << std::endl;
      std::cout << "Fractures Information: " << std::endl;
      std::map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
      for (it = m_fracprops.begin(); it != m_fracprops.end(); it++) {
        int matfracid = it->first;
        TFracProperties::FracProp prop = it->second;
        std::cout << "Material: " << matfracid << std::endl;
        std::cout << "bcMaterial: ";
        for (auto bcid : prop.m_fracbc) {
          std::cout << bcid << "\t";
        }
        std::cout << std::endl;
        std::cout << "Inersection material id" << prop.m_fracIntersectMatID << std::endl;
        std::cout << "Permeability: " << prop.m_perm << std::endl;
        std::cout << "Width: " << prop.m_width << std::endl;
        std::cout << "Porosity: " << prop.m_porosity << std::endl;
        std::cout << "***************************************" << std::endl;
      }
    }
    const bool isFracMatId(const int matid) {
      return m_fracprops.find(matid) != m_fracprops.end();
    }
    const bool isFracBCMatId(const int matid) {
      std::map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
      for (it = m_fracprops.begin(); it != m_fracprops.end(); it++) {
        int matfracid = it->first;
        TFracProperties::FracProp prop = it->second;
        for (auto bcid : prop.m_fracbc) {
          if (bcid == matid) {
            return true;
          }
        }
      }
      return false;
    }
  };

  class TFracIntersectProperties : public TMRSSavable {
   public:
    // material id of the intersection as given by the mesh generator
    int m_IntersectionId;
    // material id for representing permeability between overlapping fractures
    int m_FractureGlueId = -10000;
    // permeability of the porous media between the fractures
    // maybe this value should come from the original problem description
    REAL m_FractureGluePerm = 0.;
    // fracture ids corresponding to this intersection
    std::pair<int, int> fractureids;
    // fracture material ids corresponding to this intersection
    std::pair<int, int> fracturematids;
    // the material id of the HDiv wrappers
    int m_IntersectionPressureLossId;

    /**
     * @brief Default constructor
     */
    TFracIntersectProperties() {
      m_IntersectionId = -10000;
      m_IntersectionPressureLossId = -10000;
    }
    /**
     * @brief Destructor
     */
    ~TFracIntersectProperties() {
    }

    /**
     * @brief Copy constructor
     */
    TFracIntersectProperties(const TFracIntersectProperties &other) = default;
    /**
     * @brief Copy assignment operator
     */
    TFracIntersectProperties &operator=(const TFracIntersectProperties &other) = default;

    bool operator==(const TFracIntersectProperties &other) {
      // check for self-assignment
      if (&other == this) {
        return true;
      }

      return m_IntersectionId == other.m_IntersectionId &&
             m_IntersectionPressureLossId == other.m_IntersectionPressureLossId &&
             m_FractureGlueId == other.m_FractureGlueId &&
             m_FractureGluePerm == other.m_FractureGluePerm;
    }

    void Write(TPZStream &buf, int withclassid) const {  // ok
      buf.Write(&m_IntersectionId);
      buf.Write(&m_IntersectionPressureLossId);
      DebugStop();
    }

    void Read(TPZStream &buf, void *context) {  // ok
      buf.Read(&m_IntersectionId);
      buf.Read(&m_IntersectionPressureLossId);
      DebugStop();
    }

    virtual int ClassId() const {
      return Hash("TMRSDataTransfer::TFracIntersectProperties");
    }

    void Print() const {
      std::cout << "m_IntersectionId: " << m_IntersectionId << std::endl;
      std::cout << "m_IntersectionPressureLossId: " << m_IntersectionPressureLossId << std::endl;
    }

    const bool isThereFractureIntersection() {
      return m_IntersectionId != -10000;
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
  // @TODO defines exclusively the permeability of the fractures
  // what about fracture width, fracture porosity? How do you compute the fracture volume
  TFracProperties mTFracProperties;
  // @TODO contains a unique material id of the pressure lagrange multiplier and a material id of an H(div) wrapper?
  TFracIntersectProperties mTFracIntersectProperties;
};

#endif /* TMRSDataTransfer_h */
