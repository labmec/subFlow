//
//  Created by Giovane Avancini on 23/12/25.
//

#pragma once

#include "TPZLinearAnalysis.h"
#include "TSFDarcyAnalysis.h"
#include "TSFDataTransfer.h"
#include "TSFProblemData.h"
#include "TSFTransportAnalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include <stdio.h>

class TSFSFIAnalysis : public TPZLinearAnalysis {

public:
  /// Data transfer object
  TSFProblemData *fSimData;

  /// Number of iterations
  int fKiteration = 0;

  /// Flag used to solve Darcy problem only once if Linear Tracer is considered
  bool fShouldSolveDarcy = true;

  TSFDarcyAnalysis fDarcyAnalysis;

  TSFTransportAnalysis fTransportAnalysis;

  TSFDataTransfer fDataTransfer;

  /// Hold the pressure/flux solution
  TPZFMatrix<STATE> fDarcySolution;

  /// Hold the saturation solution
  TPZFMatrix<STATE> fTransportSolution;

  /// Constructor based on a cmesh and optimization band directive
  TSFSFIAnalysis(TPZMultiphysicsCompMesh *darcy_cmesh, TPZCompMesh *transport_cmesh, const RenumType &renumtype = RenumType::EDefault);

  ~TSFSFIAnalysis();

  /// Set data transfer object
  void SetProblemData(TSFProblemData *simData);
};