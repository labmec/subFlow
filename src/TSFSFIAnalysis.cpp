//
//  Created by Giovane Avancini on 23/12/25.
//

#include "TSFSFIAnalysis.h"

TSFSFIAnalysis::TSFSFIAnalysis() {
  fDarcyAnalysis = nullptr;
  fTransportAnalysis = nullptr;
  fKiteration = 0;
  fSimData = nullptr;
  fDarcySolution.Resize(0, 0);
  fTransportSolution.Resize(0, 0);
}

TSFSFIAnalysis::TSFSFIAnalysis(TPZMultiphysicsCompMesh *darcy_cmesh, TPZCompMesh *transport_cmesh, const RenumType &renumtype) {
  fDarcyAnalysis = TSFDarcyAnalysis(darcy_cmesh, renumtype);
  fTransportAnalysis = TSFTransportAnalysis(transport_cmesh, renumtype);
  fDataTransfer.SetMeshes(darcy_cmesh, transport_cmesh);
  fKiteration = 0;
  fSimData = nullptr;
  fDarcySolution.Resize(0, 0);
  fTransportSolution.Resize(0, 0);
}

TSFSFIAnalysis::~TSFSFIAnalysis() {}

void TSFSFIAnalysis::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
  fDarcyAnalysis.SetProblemData(simData);
  fTransportAnalysis.SetProblemData(simData);
}