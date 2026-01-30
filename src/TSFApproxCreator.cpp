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
  }
  TPZHDivApproxCreator::SetShouldCondense(false); // the static condensation will be done inside this class, using TPZFastCondensedElements
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
  CondenseElements(cmesh, lagmultilevel - 1, true);

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
        auto laglevel = c.LagrangeMultiplier();
        auto nelcon = c.NElConnected();
        if ((laglevel >= LagrangeLevelNotCondensed && nelcon == 1)) {
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

void TSFApproxCreator::BuildTransportCmesh() {
  int dimension = fSimData->fTGeometry.fDimension;

  fTransportMesh = new TPZCompMesh(TPZHDivApproxCreator::fGeoMesh);
  fTransportMesh->SetDimModel(dimension);
  fTransportMesh->SetDefaultOrder(0);
  fTransportMesh->SetAllCreateFunctionsDiscontinuous();

  // Insert materials
  TSFTransportMaterial *matTransport = nullptr;
  std::map<std::string, int> &domainNameAndMatId = fSimData->fTGeometry.fDomainNameAndMatId;
  for (const auto &domainPair : domainNameAndMatId) {
    std::string name = domainPair.first;
    int matId = domainPair.second;
    matTransport = new TSFTransportMaterial(matId, dimension);
    fTransportMesh->InsertMaterialObject(matTransport);
  }

  // Insert boundary conditions
  std::map<int, std::pair<int, REAL>> &BCTransportMatIdToTypeValue = fSimData->fTBoundaryConditions.fBCTransportMatIdToTypeValue;
  for (const auto &bcPair : BCTransportMatIdToTypeValue) {
    int bcMatId = bcPair.first;
    int bcType = bcPair.second.first;
    REAL bcValue = bcPair.second.second;
    TPZManVector<REAL, 1> val2(1, bcValue); // Part that goes to the RHS vector
    TPZFMatrix<REAL> val1(1, 1, 0.);        // Part that goes to the Stiffnes matrix
    TPZBndCondT<REAL> *bcond = matTransport->CreateBC(matTransport, bcMatId, bcType, val1, val2);
    fTransportMesh->InsertMaterialObject(bcond);
  }

  fTransportMesh->AutoBuild();
  fTransportMesh->LoadReferences();

  // Insert interface
  int interfaceMatId = fSimData->fTGeometry.fInterface_material_id;
  TSFTransportMaterial *matInterface = new TSFTransportMaterial(interfaceMatId, dimension - 1);
  fTransportMesh->InsertMaterialObject(matInterface);
  CreateInterfaceElements();
}

void TSFApproxCreator::CreateInterfaceElements() {
  TPZGeoMesh *gmesh = TPZHDivApproxCreator::fGeoMesh;
  int dim = gmesh->Dimension();
  int64_t nel = gmesh->NElements();
  int interfaceMatId = fSimData->fTGeometry.fInterface_material_id;

  std::map<int, std::pair<int, REAL>> &BCTransportMatIdToTypeValue = fSimData->fTBoundaryConditions.fBCTransportMatIdToTypeValue;
  std::set<int> bcMatIds;
  for (const auto &bcPair : BCTransportMatIdToTypeValue) {
    int bcMatId = bcPair.first;
    bcMatIds.insert(bcMatId);
  }

  // Interface Element between boundary and domain elements
  for (auto const &BcMatID : bcMatIds) {
    for (int64_t el = 0; el < nel; el++) {
      int meshDim = gmesh->Dimension();

      TPZGeoEl *geoEl = gmesh->Element(el);
      int matID = geoEl->MaterialId();

      if (matID != BcMatID) continue;

      int nSides = geoEl->NSides();
      TPZGeoElSide geoElSide(geoEl, nSides - 1);
      TPZCompElSide compElSide = geoElSide.Reference();

      TPZStack<TPZGeoElSide> neighbourSet;
      geoElSide.AllNeighbours(neighbourSet);

      int64_t nneighs = neighbourSet.size();

      for (int stack_i = 0; stack_i < nneighs; stack_i++) {
        TPZGeoElSide neighbour = neighbourSet[stack_i];
        int neighMatID = neighbour.Element()->MaterialId();
        TPZCompElSide compElNeigh = neighbour.Reference();

        int64_t neighIndex = neighbour.Element()->Index();

        if (neighbour.Element()->Dimension() != meshDim) continue;

        if (neighbour.Element()->HasSubElement())
          DebugStop();

        TPZGeoElBC gbc(neighbour, interfaceMatId);
        TPZInterfaceElement *mp_interface_el = new TPZInterfaceElement(*fTransportMesh, gbc.CreatedElement(), compElSide, compElNeigh);
      }
    }
  }

  // Interface Element between domain neighbour elements
  for (int64_t el = 0; el < nel; el++) {
    TPZGeoEl *geoEl = gmesh->Element(el);

    if (!geoEl) continue;
    if (geoEl->HasSubElement()) continue;
    if (geoEl->Dimension() != dim) continue;

    int nside = geoEl->NSides();

    for (int side = 0; side < nside - 1; side++) {
      if (geoEl->SideDimension(side) != dim - 1) continue;

      TPZGeoElSide geoElSide(geoEl, side);
      TPZCompElSide compElSide = geoElSide.Reference();
      TPZGeoElSide neighbour = geoElSide.Neighbour();

      TPZStack<TPZGeoElSide> neighbourSet;
      neighbour.AllNeighbours(neighbourSet);
      bool hasInterface = false;
      for (auto const &neigh : neighbourSet) {
        if (neigh.Element()->MaterialId() == interfaceMatId) {
          hasInterface = true;
          break;
        }
      }
      if (hasInterface) continue;
      if (neighbour == geoElSide) continue;
      if (neighbour.Element()->HasSubElement()) continue;

      while (neighbour != geoElSide) {
        if (neighbour.Element()->MaterialId() == interfaceMatId) {
          break;
        } else if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
          TPZGeoElBC gbc(geoElSide, interfaceMatId);
          TPZCompElSide compElNeigh = neighbour.Reference();
          TPZInterfaceElement *mp_interface_el = new TPZInterfaceElement(*fTransportMesh, gbc.CreatedElement(), compElSide, compElNeigh);
          break;
        }
        neighbour = neighbour.Neighbour();
      }
    }
  }

  gmesh->BuildConnectivity();
}

TPZCompMesh *TSFApproxCreator::GetTransportCmesh() {
  return fTransportMesh;
}