//
//  Created by Giovane Avancini on 02/09/25.
//

#pragma once

#include "TPZFastCondensedElement.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "TSFMixedDarcy.h"
#include "TSFProblemData.h"
#include "pzmultiphysicselement.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <stdio.h>

class TSFDarcyAnalysis : public TPZLinearAnalysis {

public:
  /// Data transfer object
  TSFProblemData *fSimData;

  /// Number of iterations
  int fKiteration = 0;

  bool fIsFirstAssemble = true;

  /// Default constructor
  TSFDarcyAnalysis();

  /// Constructor based on a cmesh and optimization band directive
  TSFDarcyAnalysis(TPZMultiphysicsCompMesh *cmesh, const RenumType &renumtype = RenumType::EDefault);

  /// Default destructor
  ~TSFDarcyAnalysis();

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

  /// Fill the set of material ids related to Neumann BCs
  /// This method is used to apply Neumann BCs without using BigNumbers
  void FillNeumannBCMatids(std::set<int> &neumannMatids);

  /// Set the initial value of the boundary conditions
  /// This method is used to apply Neumann BCs without using BigNumbers
  void SetInitialBCValue(std::set<int> &neumannMatids);

  /// Set the initial solution
  void SetInitialSolution();

  /// Remove equations related to Neumann BCs from the system
  /// This method is used to apply Neumann BCs without using BigNumbers
  void ApplyEquationFilter(std::set<int> &neumannMatids);

  /// Verifies if the sum of the fluxes over all faces of an element is zero
  void VerifyElementFluxes();

  /// Update the density and coefficients
  void UpdateDensityAndCoefficients();

  /// Update the last state variables (saturation and pressure)
  void SetLastStateVariables();
  void SetLastStateSaturation();
  void SetLastStatePressure();

  /// Sets the current time
  void SetTime(REAL time);
};
