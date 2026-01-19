#include "TPZFastCondensedElement.h"
#include "TPZHDivApproxCreator.h"
#include "TPZInterfaceEl.h"
#include "TPZMatInterfaceSingleSpace.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZVTKGeoMesh.h"
#include "TSFMixedDarcy.h"
#include "TSFProblemData.h"
#include "TSFTransportMaterial.h"

class TSFApproxCreator : public TPZHDivApproxCreator {
public:
  TSFApproxCreator();

  TSFApproxCreator(TPZGeoMesh *gmesh);

  ~TSFApproxCreator();

  void SetProblemData(TSFProblemData *simData);

  TSFProblemData *GetProblemData();

  void ConfigureDarcySpace();

  void AddDarcyMaterials();

  /// Driver function. Will create the Darcy atomic meshes (HDiv, L2, and rigid bodies) and an associate multiphysics mesh.
  /// Also, it creates the auxiliary transport mesh to post-process the results.
  TPZMultiphysicsCompMesh *CreateApproximationSpace() override;

  void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix = true);

  /// Creates auxiliary cmesh for transport. Mostly used to identify interfaces
  void BuildTransportCmesh();

  /// insert the necessary interface elements
  void CreateInterfaceElements();

  TPZCompMesh *GetTransportCmesh();

protected:
  TSFProblemData *fSimData;

  TPZCompMesh *fTransportMesh;
};