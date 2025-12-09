#include "TPZHDivApproxCreator.h"
#include "TSFMixedDarcy.h"
#include "TSFProblemData.h"

class TSFApproxCreator : public TPZHDivApproxCreator {
public:
  TSFApproxCreator();

  TSFApproxCreator(TPZGeoMesh *gmesh);

  ~TSFApproxCreator();

  void SetProblemData(TSFProblemData *simData);

  TSFProblemData *GetProblemData();

  /// Create a discontinuous mesh for transport
  TPZCompMesh *TransportCmesh();

  /// Creates auxiliary cmesh for transport. Mostly used to identify interfaces
  void BuildAuxTransportCmesh();

  /// insert the necessary interface elements
  void InsertInterfaceElements();

protected:
  TSFProblemData *fSimData;

  TPZCompMesh *fTransportMesh;
};