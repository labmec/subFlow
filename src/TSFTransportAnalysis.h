//
//  Created by Giovane Avancini on 04/11/25.
//

#pragma once

#include "TPZLinearAnalysis.h"
#include "TPZYSMPMatrix.h"
#ifdef PZ_USING_MKL
#include "TPZSpStructMatrix.h"
#include "TPZYSMPPardiso.h"
#endif
#include "TSFAlgebraicTransport.h"
#include "TSFProblemData.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include <stdio.h>

class TSFTransportAnalysis : public TPZLinearAnalysis {

public:
  /// Data transfer object
  TSFProblemData *fSimData;

  /// Algebraic transport object
  TSFAlgebraicTransport fAlgebraicTransport;

  /// Number of iterations
  int fKiteration = 0;

  bool fIsFirstAssemble = true;
#ifdef PZ_USING_MKL
  TPZFYsmpMatrixPardiso<REAL> *fTransmissibilityMatrix;
#else
  TPZFYsmpMatrix<REAL> *fTransmissibilityMatrix;
#endif
  TPZFMatrix<REAL> fMassMatrix;

  /// Default constructor
  TSFTransportAnalysis();

  /// Constructor based on a cmesh and optimization band directive
  TSFTransportAnalysis(TPZCompMesh *cmesh, const RenumType &renumtype = RenumType::EDefault);

  /// Default destructor
  ~TSFTransportAnalysis();

  /// Configurates iternal members
  void Initialize();

  /// Set data transfer object
  void SetProblemData(TSFProblemData *simData);

  /// Get data transfer object
  TSFProblemData *GetProblemData();

  /// Get the number of iterations
  int GetNumberOfIterations();

  /// Run a time step
  void RunTimeStep(std::ostream &out = std::cout);

  /// Render a vtk file with requested variables for a time step
  void PostProcessTimeStep(int dimToPost = -1, int step = -1);

  /// override solve to have timers
  void Solve() override;

  /// @brief Assemble the mass matrix. This is done only once at the first Assemble() call.
  void AssembleMass();

  /// @brief Assemble the transmissibility matrix and the rhs.
  void Assemble() override;

  /// Verifies if the sum of the fluxes over all faces of an element is zero
  void VerifyElementFluxes();

  /// Update the density and coefficients
  void UpdateDensityAndCoefficients();

  /// SetInitialSaturation
  void SetInitialSaturation();

  /// Update the last state variables (saturation and pressure)
  void SetLastStateVariables();
  void SetLastStateSaturation();
  void SetLastStateDensities();
  void SetLastStateVolumeFactor();
};
