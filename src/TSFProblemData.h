//
//  Created by Giovane Avancini on 05/09/25.
//

#pragma once

#include "Material/TPZMatTypes.h"
#include "TPZStream.h"
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

  std::string fSimulationName = "";

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
    std::map<std::string, int> fDomainNameAndMatId;

    /** @brief
     MaterialID of the interface element that will be inserted in the transport mesh
     */
    int fInterface_material_id = 100; // Interface material Volumetric-Volumetric

    // definition of the standard material ids
    /// material id for the Macro fluxes
    int fSkeletonMatId = 19;
    /// material id for H(div) wrap materials
    int fHdivWrapMatId = 15;
    /// material id for the mortar pressure space
    int fMortarMatId = 20;
    /// material id for a positive lagrange multiplier interface
    int fPosLagrangeMatId = 30;
    /// material id for a negative lagrange multiplier interface
    int fNegLagrangeMatId = 35;
    /// material id for a zero order H(div) boundary flux
    int fZeroOrderHdivFluxMatId = 40;
    /// material id for pressure lagrange multiplier (same as mortar because they should not be in the same mesh)
    int fPressureMatId = 20;

    std::string fGmeshFileName = "";

    /** @brief Default constructor */
    TGeometry() {}

    /** @brief Destructor */
    ~TGeometry() {}

    /** @brief Copy constructor */
    TGeometry(const TGeometry &other) {
      fDomainNameAndMatId = other.fDomainNameAndMatId;
      fInterface_material_id = other.fInterface_material_id;

      fPressureMatId = other.fPressureMatId;
      /// material id for the Macro fluxes
      fSkeletonMatId = other.fSkeletonMatId;
      /// material id for H(div) wrap materials
      fHdivWrapMatId = other.fHdivWrapMatId;
      /// material id for the mortar pressure space
      fMortarMatId = other.fMortarMatId;
      /// material id for a positive lagrange multiplier interface
      fPosLagrangeMatId = other.fPosLagrangeMatId;
      /// material id for a negative lagrange multiplier interface
      fNegLagrangeMatId = other.fNegLagrangeMatId;
      /// material id for a zero order H(div) boundary flux
      fZeroOrderHdivFluxMatId = other.fZeroOrderHdivFluxMatId;

      fGmeshFileName = other.fGmeshFileName;
    }
    /** @brief Copy assignment operator*/
    TGeometry &operator=(const TGeometry &other) {
      if (this != &other) // prevent self-assignment
      {
        fDomainNameAndMatId = other.fDomainNameAndMatId;
        fInterface_material_id = other.fInterface_material_id;
        fPressureMatId = other.fPressureMatId;
        /// material id for the Macro fluxes
        fSkeletonMatId = other.fSkeletonMatId;
        /// material id for H(div) wrap materials
        fHdivWrapMatId = other.fHdivWrapMatId;
        /// material id for the mortar pressure space
        fMortarMatId = other.fMortarMatId;
        /// material id for a positive lagrange multiplier interface
        fPosLagrangeMatId = other.fPosLagrangeMatId;
        /// material id for a negative lagrange multiplier interface
        fNegLagrangeMatId = other.fNegLagrangeMatId;
        /// material id for a zero order H(div) boundary flux
        fZeroOrderHdivFluxMatId = other.fZeroOrderHdivFluxMatId;

        fGmeshFileName = other.fGmeshFileName;
      }
      return *this;
    }
  };

  /**
   * @brief Class that stores PetroPhysics information, i.e. Interaction between fluid and rock.
   */

  class TPetroPhysics : public TSFSavable {
  public:
    REAL fGasViscosity;
    REAL fWaterViscosity;
    REAL fSwr;    // residual water saturation
    REAL fSgr;    // residual gas saturation
    int fKrModel; // 0 - Linear, 1 - Quadratic, 2 - Quadratic with residual

    std::vector<std::function<std::tuple<REAL, REAL>(REAL &)>> fKrg;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &)>> fKrw;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &, REAL &)>> fFg;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &, REAL &)>> fFw;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &)>> fLambdaw;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &)>> fLambdag;
    std::vector<std::function<std::tuple<REAL, REAL>(REAL &, REAL &, REAL &)>> fLambdaTotal;

    /** @brief Default constructor */
    TPetroPhysics() {
      fGasViscosity = 1.0;
      fWaterViscosity = 1.0;
      fSwr = 0.0;
      fSgr = 0.0;
      fKrModel = 0;
      fKrg.resize(3);
      fKrw.resize(3);
      fFg.resize(3);
      fFw.resize(3);
      fLambdaw.resize(3);
      fLambdag.resize(3);
      fLambdaTotal.resize(3);
      CreateLinearKrModel();
      CreateQuadraticKrModel();
    }

    /** @brief Destructor */
    ~TPetroPhysics() {
    }

    /** @brief Copy constructor */
    TPetroPhysics(const TPetroPhysics &other) {
      fGasViscosity = other.fGasViscosity;
      fWaterViscosity = other.fWaterViscosity;
      fSwr = other.fSwr;
      fSgr = other.fSgr;
      fKrModel = other.fKrModel;
      fKrg = other.fKrg;
      fKrw = other.fKrw;
      fFg = other.fFg;
      fFw = other.fFw;
      fLambdaw = other.fLambdaw;
      fLambdag = other.fLambdag;
      fLambdaTotal = other.fLambdaTotal;
    }

    /** @brief Copy assignment operator*/
    TPetroPhysics &operator=(const TPetroPhysics &other) {
      if (this != &other) // prevent self-assignment
      {
        fGasViscosity = other.fGasViscosity;
        fWaterViscosity = other.fWaterViscosity;
        fSwr = other.fSwr;
        fSgr = other.fSgr;
        fKrModel = other.fKrModel;
        fKrg = other.fKrg;
        fKrw = other.fKrw;
        fFg = other.fFg;
        fFw = other.fFw;
        fLambdaw = other.fLambdaw;
        fLambdag = other.fLambdag;
        fLambdaTotal = other.fLambdaTotal;
      }
      return *this;
    }
    void CreateLinearKrModel();
    void CreateQuadraticKrModel();
    void UpdateLambdasAndFracFlows(int krModel);
  };

  /**
   * @brief Class that stores Fluid properties. For instance density of water and Gas, compressibility of water and Gas.
   */
  class TFluidProperties : public TSFSavable {
  public:
    REAL fGasViscosity;
    REAL fWaterViscosity;
    REAL fWaterCompressibility;
    REAL fGasCompressibility;
    REAL fGasDensityRef;
    REAL fWaterDensityRef;
    REAL fReferencePressure;
    std::function<std::tuple<REAL, REAL>(REAL &)> fGasDensityFunc;
    std::function<std::tuple<REAL, REAL>(REAL &)> fWaterDensityFunc;

    /** @brief Default constructor */
    TFluidProperties() {
      fGasViscosity = 1.0;
      fWaterViscosity = 1.0;
      fGasDensityRef = 1000;
      fWaterDensityRef = 1000;
      fWaterCompressibility = 1.0e-8;
      fGasCompressibility = 1.0e-7;
      fReferencePressure = 1.013e5;
      fGasDensityFunc = 0;
      fWaterDensityFunc = 0;
      CreateLinearDensityFunction();
    }

    /** @brief Destructor */
    ~TFluidProperties() {
    }

    /** @brief Copy constructor */
    TFluidProperties(const TFluidProperties &other) {
      fGasViscosity = other.fGasViscosity;
      fWaterViscosity = other.fWaterViscosity;
      fWaterCompressibility = other.fWaterCompressibility;
      fGasCompressibility = other.fGasCompressibility;

      fReferencePressure = other.fReferencePressure;
      fGasDensityRef = other.fGasDensityRef;
      fWaterDensityRef = other.fWaterDensityRef;
      fGasDensityFunc = other.fGasDensityFunc;
      fWaterDensityFunc = other.fWaterDensityFunc;
    }

    /** @brief Copy assignment operator*/
    TFluidProperties &operator=(const TFluidProperties &other) {
      fGasViscosity = other.fGasViscosity;
      fWaterViscosity = other.fWaterViscosity;
      fWaterCompressibility = other.fWaterCompressibility;
      fGasCompressibility = other.fGasCompressibility;
      fGasDensityRef = other.fGasDensityRef;
      fWaterDensityRef = other.fWaterDensityRef;
      fReferencePressure = other.fReferencePressure;
      fGasDensityFunc = other.fGasDensityFunc;
      fWaterDensityFunc = other.fWaterDensityFunc;
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
    std::map<int, std::pair<REAL, REAL>> fPorosityAndPermeability;

    // Function that given a point (x,y,z) returns s0 (initial saturation)
    std::function<REAL(const TPZVec<REAL> &)> fS0;

    /** @brief Default constructor */
    TReservoirProperties() {
      fPorosityAndPermeability[0] = std::make_pair(1.0, 1.0);
    }

    /** @brief Destructor */
    ~TReservoirProperties() {
    }

    /** @brief Copy constructor */
    TReservoirProperties(const TReservoirProperties &other) {
      fPorosityAndPermeability = other.fPorosityAndPermeability;
      fS0 = other.fS0;
    }

    /** @brief Copy assignment operator*/
    TReservoirProperties &operator=(const TReservoirProperties &other) {
      fPorosityAndPermeability = other.fPorosityAndPermeability;
      fS0 = other.fS0;
      return *this;
    }
  };

  /**
   * @brief Class that stores the boundary conditions of the problem
   */
  class TBoundaryConditions : public TSFSavable {
  public:
    std::map<std::string, int> fDomainNameAndMatId;
    /**
     * @brief Contains the boundary conditions for Darcy problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> fBCDarcyMatIdToTypeValue;

    /**
     * @brief Contains the boundary conditions for Transport problem in a map. Key = matidOfBC, value = pair<typeOfBC,valueOfBC>
     */
    std::map<int, std::pair<int, REAL>> fBCTransportMatIdToTypeValue;

    /**
     * @brief Contains the functions to be applied on the boundaries in a map. Key = matidOfBC, value = functionID (0 if no function is prescribed)
     */
    std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> fBCDarcyMatIdToFunctionId;

    /**
     * @brief Contains the saturation functions to be applied on the boundaries in a map. Key = matidOfBC, value = functionID (0 if no function is prescribed)
     */
    std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> fBCTransportMatIdToFunctionId;

    /** @brief Default constructor */
    TBoundaryConditions() {}

    /** @brief Destructor */
    ~TBoundaryConditions() {
    }

    /** @brief Copy constructor */
    TBoundaryConditions(const TBoundaryConditions &other) {
      fBCDarcyMatIdToTypeValue = other.fBCDarcyMatIdToTypeValue;
      fBCTransportMatIdToTypeValue = other.fBCTransportMatIdToTypeValue;
      fDomainNameAndMatId = other.fDomainNameAndMatId;
      fBCDarcyMatIdToFunctionId = other.fBCDarcyMatIdToFunctionId;
      fBCTransportMatIdToFunctionId = other.fBCTransportMatIdToFunctionId;
    }

    /** @brief Copy assignment operator*/
    TBoundaryConditions &operator=(const TBoundaryConditions &other) {
      if (this != &other) // prevent self-assignment
      {
        fBCDarcyMatIdToTypeValue = other.fBCDarcyMatIdToTypeValue;
        fBCTransportMatIdToTypeValue = other.fBCTransportMatIdToTypeValue;
        fDomainNameAndMatId = other.fDomainNameAndMatId;
        fBCDarcyMatIdToFunctionId = other.fBCDarcyMatIdToFunctionId;
        fBCTransportMatIdToFunctionId = other.fBCTransportMatIdToFunctionId;
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
    REAL fDt;

    /**
     * @brief Flag to indicate the analysis type: 1 - only Darcy, 2 - only transport, 3 - coupled
     */
    int fAnalysisType;

    /**
     * @brief Residual tolerance for Darcy
     */
    REAL fResTolDarcy;

    /**
     * @brief Residual tolerance for transport
     */
    REAL fResTolTransport;

    /**
     * @brief Correction tolerance for darcy
     */
    REAL fCorrTolDarcy;

    /**
     * @brief Correction tolerance for transport
     */
    REAL fCorrTolTransport;

    REAL fSfiTol;
    /**
     * @brief Maximum number of iterations per time step for darcy
     */
    int fMaxIterDarcy;

    /**
     * @brief Maximum number of iterations per time step for transport
     */
    int fMaxIterTransport;

    /**
     * @brief Maximum number of Sequential Fully Implicit (SFI) iterations per time step
     */
    int fMaxIterSfi;

    /**
     * @brief Number of time steps
     */
    int fNSteps;

    // Order of approximation for the border of the Pressure element
    int fMortarBorderElementPresOrder;
    // Order of approximation for the border of the Flux element
    int fMortarBorderElementFluxOrder;
    /**
     * @brief Directive for the use of four spaces
     */
    bool fFourApproxSpaces;

    /**
     * @brief Axisymmetric flag
     */
    bool fIsAxisymmetric;

    /**
     * @brief Linear trace flag. If true, the darcy problem is ran only once
     */
    bool fIsLinearTrace;

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

    MSpaceType fSpaceType = ENone;

    std::vector<REAL> fGravity;

    int fNThreadsDarcyProblem = 0;

    /** @brief Default constructor */
    TNumerics() {
      fDt = 0.0;
      fAnalysisType = 3;
      fResTolDarcy = 1.0e-4;
      fResTolTransport = 1.0e-7;
      fCorrTolDarcy = 1.0e-4;
      fCorrTolTransport = 1.0e-7;
      fMaxIterDarcy = 0;
      fMaxIterTransport = 0;
      fMaxIterSfi = 0;
      fSfiTol = 1.0e-3;
      fNSteps = 0;
      fFourApproxSpaces = false;
      fIsAxisymmetric = false;
      fIsLinearTrace = true;
      fSpaceType = ENone;
      fGravity.resize(3, 0.0);
      fNThreadsDarcyProblem = 0;
      fMortarBorderElementPresOrder = 0;
      fMortarBorderElementFluxOrder = 0;
    }
    /** @brief Destructor */
    ~TNumerics() {
    }

    /** @brief Copy constructor */
    TNumerics(const TNumerics &other) {
      fDt = other.fDt;
      fAnalysisType = other.fAnalysisType;
      fResTolDarcy = other.fResTolDarcy;
      fResTolTransport = other.fResTolTransport;
      fCorrTolDarcy = other.fCorrTolDarcy;
      fCorrTolTransport = other.fCorrTolTransport;
      fMaxIterDarcy = other.fMaxIterDarcy;
      fMaxIterTransport = other.fMaxIterTransport;
      fMaxIterSfi = other.fMaxIterSfi;
      fSfiTol = other.fSfiTol;
      fNSteps = other.fNSteps;
      fFourApproxSpaces = other.fFourApproxSpaces;
      fIsAxisymmetric = other.fIsAxisymmetric;
      fIsLinearTrace = other.fIsLinearTrace;
      fSpaceType = other.fSpaceType;
      fGravity = other.fGravity;
      fNThreadsDarcyProblem = other.fNThreadsDarcyProblem;
      fMortarBorderElementPresOrder = other.fMortarBorderElementPresOrder;
      fMortarBorderElementFluxOrder = other.fMortarBorderElementFluxOrder;
    }

    /** @brief Copy assignment operator*/
    TNumerics &operator=(const TNumerics &other) {
      // check for self-assignment
      if (&other == this) {
        return *this;
      }

      fDt = other.fDt;
      fAnalysisType = other.fAnalysisType;
      fResTolDarcy = other.fResTolDarcy;
      fResTolTransport = other.fResTolTransport;
      fCorrTolDarcy = other.fCorrTolDarcy;
      fCorrTolTransport = other.fCorrTolTransport;
      fMaxIterDarcy = other.fMaxIterDarcy;
      fMaxIterTransport = other.fMaxIterTransport;
      fMaxIterSfi = other.fMaxIterSfi;
      fSfiTol = other.fSfiTol;
      fNSteps = other.fNSteps;
      fFourApproxSpaces = other.fFourApproxSpaces;
      fIsAxisymmetric = other.fIsAxisymmetric;
      fIsLinearTrace = other.fIsLinearTrace;
      fSpaceType = other.fSpaceType;
      fGravity = other.fGravity;
      fNThreadsDarcyProblem = other.fNThreadsDarcyProblem;
      fMortarBorderElementPresOrder = other.fMortarBorderElementPresOrder;
      fMortarBorderElementFluxOrder = other.fMortarBorderElementFluxOrder;
      return *this;
    }

    bool operator==(const TNumerics &other) {
      // check for self-assignment
      if (&other == this) {
        return true;
      }

      return fDt == other.fDt &&
             fAnalysisType == other.fAnalysisType &&
             fResTolDarcy == other.fResTolDarcy &&
             fResTolTransport == other.fResTolTransport &&
             fCorrTolDarcy == other.fCorrTolDarcy &&
             fCorrTolTransport == other.fCorrTolTransport &&
             fMaxIterDarcy == other.fMaxIterDarcy &&
             fMaxIterTransport == other.fMaxIterTransport &&
             fMaxIterSfi == other.fMaxIterSfi &&
             fSfiTol == other.fSfiTol &&
             fNSteps == other.fNSteps &&
             fFourApproxSpaces == other.fFourApproxSpaces &&
             fIsAxisymmetric == other.fIsAxisymmetric &&
             fIsLinearTrace == other.fIsLinearTrace &&
             fSpaceType == other.fSpaceType &&
             fGravity == other.fGravity &&
             fNThreadsDarcyProblem == other.fNThreadsDarcyProblem &&
             fMortarBorderElementPresOrder == other.fMortarBorderElementPresOrder &&
             fMortarBorderElementFluxOrder == other.fMortarBorderElementFluxOrder;
    }

    void Write(TPZStream &buf, int withclassid) const { // ok
      buf.Write(&fDt);
      buf.Write(&fAnalysisType);
      buf.Write(&fResTolDarcy);
      buf.Write(&fResTolTransport);
      buf.Write(&fCorrTolDarcy);
      buf.Write(&fMaxIterDarcy);
      buf.Write(&fMaxIterTransport);
      buf.Write(&fMaxIterSfi);
      buf.Write(&fSfiTol);
      buf.Write(&fNSteps);
      int temp = fFourApproxSpaces;
      buf.Write(&temp);
      temp = fIsAxisymmetric;
      buf.Write(&temp);
      temp = fIsLinearTrace;
      buf.Write(&temp);
      temp = fSpaceType;
      buf.Write(&temp);
      buf.Write(fGravity);
      buf.Write(fNThreadsDarcyProblem);
    }

    void Read(TPZStream &buf, void *context) { // ok
      buf.Read(&fDt);
      buf.Read(&fAnalysisType);
      buf.Read(&fResTolDarcy);
      buf.Read(&fResTolTransport);
      buf.Read(&fCorrTolDarcy);
      buf.Read(&fCorrTolTransport);
      buf.Read(&fMaxIterDarcy);
      buf.Read(&fMaxIterTransport);
      buf.Read(&fMaxIterSfi);
      buf.Read(&fSfiTol);
      buf.Read(&fNSteps);
      int temp;
      buf.Read(&temp);
      fFourApproxSpaces = temp;
      buf.Read(&temp);
      fIsAxisymmetric = temp;
      buf.Read(&temp);
      fIsLinearTrace = temp;
      buf.Read(&temp);
      fSpaceType = (MSpaceType)temp;
      buf.Read(&fNThreadsDarcyProblem);
    }

    virtual int ClassId() const {
      return Hash("TSFProblemData::TNumerics");
    }

    void Print() const {
      std::cout << fDt << std::endl;
      std::cout << fAnalysisType << std::endl;
      std::cout << fResTolDarcy << std::endl;
      std::cout << fResTolTransport << std::endl;
      std::cout << fCorrTolDarcy << std::endl;
      std::cout << fCorrTolTransport << std::endl;
      std::cout << fMaxIterDarcy << std::endl;
      std::cout << fMaxIterTransport << std::endl;
      std::cout << fMaxIterSfi << std::endl;
      std::cout << fNSteps << std::endl;
      std::cout << fFourApproxSpaces << std::endl;
      std::cout << fIsAxisymmetric << std::endl;
      std::cout << fIsLinearTrace << std::endl;
      std::cout << fSpaceType << std::endl;
      std::cout << fNThreadsDarcyProblem << std::endl;
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
    std::string fFileNameDarcy;

    /**
     * @brief Transpor operator vtk file name
     */
    std::string fFileNameTransport;

    /**
     * @brief Contains scalar variables that will be postprocessed
     */
    TPZStack<std::string, 10> fScalnamesDarcy;

    /**
     * @brief Contains scalar variables that will be postprocessed
     */
    TPZStack<std::string, 10> fScalnamesTransport;

    /**
     * @brief Contains vector variables that will be postprocessed
     */
    TPZStack<std::string, 10> fVecnamesDarcy;

    /**
     * @brief Period of time post-processed data is printed
     */
    REAL fFileTimeStep;

    /**
     * @brief Contains the times at which post-processed data is printed
     */
    TPZStack<REAL, 100> fVecReportingTimes;

    /**
     * @brief Default constructor
     */
    TPostProcess() {
      fFileNameDarcy = "";
      fFileNameTransport = "";
      fScalnamesDarcy.Resize(0);
      fScalnamesDarcy.Resize(0);
      fScalnamesTransport.Resize(2);
      fScalnamesTransport.Push("Sw");
      fScalnamesTransport.Push("Sg");
      fFileTimeStep = 0.0;
      fVecReportingTimes.Resize(0);
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
      fFileNameDarcy = other.fFileNameDarcy;
      fFileNameTransport = other.fFileNameTransport;
      fScalnamesDarcy = other.fScalnamesDarcy;
      fScalnamesDarcy = other.fScalnamesDarcy;
      fScalnamesTransport = other.fScalnamesTransport;
      fFileTimeStep = other.fFileTimeStep;
      fVecReportingTimes = other.fVecReportingTimes;
    }
    /**
     * @brief Copy assignment operator
     */
    TPostProcess &operator=(const TPostProcess &other) {
      // check for self-assignment
      if (&other == this) {
        return *this;
      }

      fFileNameDarcy = other.fFileNameDarcy;
      fFileNameTransport = other.fFileNameTransport;
      fVecnamesDarcy = other.fVecnamesDarcy;
      fScalnamesDarcy = other.fScalnamesDarcy;
      fScalnamesTransport = other.fScalnamesTransport;
      fFileTimeStep = other.fFileTimeStep;
      fVecReportingTimes = other.fVecReportingTimes;

      return *this;
    }

    bool operator==(const TPostProcess &other) {
      // check for self-assignment
      if (&other == this) {
        return true;
      }

      return fFileNameDarcy == other.fFileNameDarcy &&
             fFileNameTransport == other.fFileNameTransport &&
             fVecnamesDarcy == other.fVecnamesDarcy &&
             fScalnamesDarcy == other.fScalnamesDarcy &&
             fScalnamesTransport == other.fScalnamesTransport &&
             fFileTimeStep == other.fFileTimeStep &&
             fVecReportingTimes == other.fVecReportingTimes;
    }

    void Write(TPZStream &buf, int withclassid) const { // ok
      buf.Write(&fFileNameDarcy);
      buf.Write(&fFileNameTransport);
      buf.Write(fVecnamesDarcy);
      buf.Write(fScalnamesDarcy);
      buf.Write(fScalnamesTransport);
      buf.Write(&fFileTimeStep);
      buf.Write(fVecReportingTimes);
    }

    void Read(TPZStream &buf, void *context) { // ok
      buf.Read(&fFileNameDarcy);
      buf.Read(&fFileNameTransport);
      buf.Read(fVecnamesDarcy);
      buf.Read(fScalnamesDarcy);
      buf.Read(fScalnamesTransport);
      buf.Read(&fFileTimeStep);
      buf.Read(fVecReportingTimes);
    }

    virtual int ClassId() const {
      return Hash("TSFProblemData::TPostProcess");
    }

    void Print() const {
      std::cout << fFileNameDarcy << std::endl;
      std::cout << fFileNameTransport << std::endl;
      std::cout << fFileTimeStep << std::endl;
      std::cout << fVecReportingTimes << std::endl;
      // scalnames and vecnames
    }
  };

  // describes the material id of peripheral material objects (skeleton matid, hdivwrap matid, etc)
  // contains a data structure associating a matid with a string generated by gmesh
  TGeometry fTGeometry;
  // describes the properties of the fluid mixture (relative permeabilities)
  TPetroPhysics fTPetroPhysics;
  // describes the properties of the Gas and water
  TFluidProperties fTFluidProperties;
  // describes the initial conditions of saturation
  // describes the values of permeability by material id
  // describes the porosity of the reservoir
  // describes the volume factor (?) by material id
  //
  TReservoirProperties fTReservoirProperties;

  // defines the matids and values associated with the boundary conditions
  TBoundaryConditions fTBoundaryConditions;
  /*
   * stores the size of the timestep
   * stores the tolerances of the iterative solver
   * stores the maximum number of iterations
   * stores the type of approximation spaces to be used
   * stores the number of threads to be used in paralel computing
   * defines polynomial orders when applying hybridization
   * defines if the mesh is going to be substructured
   */
  TNumerics fTNumerics;
  // define the post processing information
  TPostProcess fTPostProcess;
};
