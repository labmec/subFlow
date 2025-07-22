//
//  TMRSMixedAnalysis.hpp
//
//  Created by Omar Durán on 10/15/19.
//

#ifndef TMRSMixedAnalysis_h
#define TMRSMixedAnalysis_h

#include <stdio.h>
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TMRSDataTransfer.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZSymetricSpStructMatrixEigen.h"
#include "TPZSSpMatrixEigen.h"

class TMRSMixedAnalysis : public TPZLinearAnalysis {
    
    

    
    
public:
    
    /// Data transfer object
    TMRSDataTransfer * m_sim_data;
    TPZMFSolutionTransfer fsoltransfer;
    
    /// Number of iterations
    int m_k_iteration;
    
    
    bool  mIsFirstAssembleQ=true;
    
    
    
    /// Default constructor
    TMRSMixedAnalysis();
    
    /// Default destructor
    ~TMRSMixedAnalysis();
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSMixedAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,
                      const RenumType& renumtype = RenumType::EDefault);
    
    /// Configurates iternal members
    void Configure(int n_threads, bool UsePardiso_Q, bool UsePZ=false);
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer * sim_data);
    
    /// Get data transfer object
    TMRSDataTransfer * GetDataTransfer();
    
    /// Get the number of iterations
    int GetNumberOfIterations();
    
    /// Run a time step
    void RunTimeStep();
    
    /// Render a vtk file with requested variables for a time step
    void PostProcessTimeStep(int dimToPost = -1, int step = -1);
    
    /// Perform a Newton iteration
    void NewtonIteration();
   
    /// override assemble to have timers
    void Assemble() override;
    
    /// override solve to have timers
    void Solve() override;
	
	/// Verifies if the sum of the fluxes over all faces of an element is zero
	void VerifyElementFluxes();
    
    /// Update the density and coefficients
    void UpdateDensityAndCoefficients();

    /// Update the last state variables (saturation and pressure)
    void SetLastStateVariables();
    void SetLastStateSaturation();
    void SetLastStatePressure();
};

#endif /* TMRSMixedAnalysis_h */
