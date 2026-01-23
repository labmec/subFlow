//
//  Created by Giovane Avancini on 04/11/25.
//

#pragma once

#include "TPZLinearAnalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZYSMPMatrix.h"
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

  TPZFYsmpMatrix<STATE> fMassMatrix;

  TPZFYsmpMatrix<STATE> fTransmissibilityMatrix;

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

  /// Perform a Newton iteration
  void NewtonIteration();

  /// override assemble to have timers
  void Assemble() override;

  /// override solve to have timers
  void Solve() override;

  void AssembleMass();

  void AssembleTransmissibility();

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
