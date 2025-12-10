#include "TSFApproxCreator.h"

TSFApproxCreator::TSFApproxCreator() : TPZHDivApproxCreator(), fSimData(nullptr), fTransportMesh(nullptr) {}

TSFApproxCreator::TSFApproxCreator(TPZGeoMesh *gmesh) : TPZHDivApproxCreator(gmesh), fSimData(nullptr), fTransportMesh(nullptr) {}

TSFApproxCreator::~TSFApproxCreator() {}

void TSFApproxCreator::SetProblemData(TSFProblemData *simData) {
  fSimData = simData;
}

TSFProblemData *TSFApproxCreator::GetProblemData() {
  return fSimData;
}

void TSFApproxCreator::ConfigureDarcySpace() {
  TPZHDivApproxCreator::ProbType() = ProblemType::EDarcy;
  TPZHDivApproxCreator::SetDefaultOrder(fSimData->fTNumerics.fFluxOrder);
  if (fSimData->fTNumerics.fFourApproxSpaces) {
    TPZHDivApproxCreator::IsRigidBodySpaces() = true;
  } else {
    TPZHDivApproxCreator::SetShouldCondense(false);
  }
}

void TSFApproxCreator::AddDarcyMaterials() {
  // Add materials to the HDivApproxCreator
  int dim = fSimData->fTGeometry.fDimension;
  TSFMixedDarcy *matDarcy = nullptr;

  std::map<std::string, int> &domainNameAndMatId = fSimData->fTGeometry.fDomainNameAndMatId;
  for (const auto &domainPair : domainNameAndMatId) {
    std::string name = domainPair.first;
    int matId = domainPair.second;
    matDarcy = new TSFMixedDarcy(matId, dim);
    matDarcy->SetAxisymmetry(fSimData->fTNumerics.fIsAxisymmetric);
    REAL permeability = fSimData->fTReservoirProperties.fPorosityAndPermeability[matId].second;
    matDarcy->SetConstantPermeability(permeability);
    matDarcy->SetFourSpaces(fSimData->fTNumerics.fFourApproxSpaces);
    TPZHDivApproxCreator::InsertMaterialObject(matDarcy);
  }

  std::map<int, std::pair<int, REAL>> &BCDarcyMatIdToTypeValue = fSimData->fTBoundaryConditions.fBCDarcyMatIdToTypeValue;
  std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>> &BCDarcyMatIdToFunctionId = fSimData->fTBoundaryConditions.fBCDarcyMatIdToFunctionId;
  for (const auto &bcPair : BCDarcyMatIdToTypeValue) {
    int bcMatId = bcPair.first;
    int bcType = bcPair.second.first;
    REAL bcValue = bcPair.second.second;
    TPZManVector<REAL, 1> val2(1, bcValue); // Part that goes to the RHS vector
    TPZFMatrix<REAL> val1(1, 1, 0.);        // Part that goes to the Stiffnes matrix
    TPZBndCondT<REAL> *bcond = matDarcy->CreateBC(matDarcy, bcMatId, bcType, val1, val2);
    int functionID = BCDarcyMatIdToFunctionId[bcMatId].first;
    if (functionID != 0) {
      auto functionBC = BCDarcyMatIdToFunctionId[bcMatId].second;
      bcond->SetForcingFunctionBC(functionBC, 2);
    }
    TPZHDivApproxCreator::InsertMaterialObject(bcond);
  }
}

TPZMultiphysicsCompMesh *TSFApproxCreator::CreateApproximationSpace() {
  int lagmultilevel = 1;
  int numDarcyMeshes = TPZHDivApproxCreator::NumMeshes();
  TPZManVector<TPZCompMesh *, 7> meshvec(numDarcyMeshes);
  TPZHDivApproxCreator::CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
  TPZMultiphysicsCompMesh *cmesh = nullptr;
  TPZHDivApproxCreator::CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);
  CondenseElements(cmesh, lagmultilevel, true);

  return cmesh;
}

void TSFApproxCreator::CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix) {
  int64_t nel = cmesh->NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    if (!cel) continue;

    int nc = cel->NConnects();
    if (LagrangeLevelNotCondensed >= 0) {
      for (int ic = 0; ic < nc; ic++) {
        TPZConnect &c = cel->Connect(ic);
        if ((c.LagrangeMultiplier() >= LagrangeLevelNotCondensed && c.NElConnected() == 1)) {
          c.IncrementElConnected();
        }
      }
    }
    int ic;
    for (ic = 0; ic < nc; ic++) {
      TPZConnect &c = cel->Connect(ic);
      if (c.HasDependency() || c.NElConnected() > 1) { // if the el has dependency or if it is a boundary
        continue;
      }
      break;
    }

    bool cancondense = (ic != (nc));
    if (cancondense) {
      TPZGeoEl *gel = cel->Reference();
      TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
      cond->SetLambda(1.0);
    }
  }

  cmesh->CleanUpUnconnectedNodes();
}

void TSFApproxCreator::BuildAuxTransportCmesh() {
  int dimension = fSimData->fTGeometry.fDimension;

  fTransportMesh = new TPZCompMesh(TPZHDivApproxCreator::fGeoMesh); // This mesh is only used to set the interface elements
}