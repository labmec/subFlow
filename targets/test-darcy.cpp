
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZHDivApproxCreator.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZRefPatternDataBase.h"
#include "TPZSSpStructMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "TSFMixedDarcy.h"
#include "pzfstrmatrix.h"
#include "pzintel.h"
#include "pzlog.h"
#include "pzmultiphysicselement.h"
#include "pzstepsolver.h"
#include "pzvec_extras.h"
#include <iostream>
#include "TPZEquationFilter.h"

// ================
// Global variables
// ================

// Exact solution
TLaplaceExample1 gexact;

// Material IDs for domain and boundaries
enum EnumMatIds {
  EMatId = 1,
  EBottom = 2,
  ERight = 3,
  ETop = 4,
  ELeft = 5,
  EFront = 6,
  EBack = 7
};

int nthreads = 0;

auto SetBoundaryCondition = [](TPZVec<REAL> x, TPZVec<REAL> u, TPZFMatrix<REAL> &du) {
  u[0] = x[1];
};

// Creates a geometric mesh using TPZGenGrid3D
TPZGeoMesh *createGeoMesh3D(const TPZManVector<int, 3> &nelDiv, const TPZManVector<REAL, 3> &minX, const TPZManVector<REAL, 3> &maxX);

// Creates a 2D geometric mesh using TPZGenGrid2D
TPZGeoMesh *createGeoMesh2D(const TPZManVector<int, 2> &nelDiv, const TPZManVector<REAL, 2> &minX, const TPZManVector<REAL, 2> &maxX);

// Creates a computational mesh for mixed approximation
TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, int order = 1, bool isNL = false);

void SolveLinear(int order, TPZGeoMesh *gmesh);

void SolveNonLinear(int order, TPZGeoMesh *gmesh);

void ApplyEquationFilter(TPZGeoMesh *gmesh, TPZCompMesh *cmesh_m, TPZStructMatrixT<STATE> &stmat);

int main(int argc, char *const argv[]) {

// Initialize logger
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Initialize uniform refinements for 1D and 2D elements
  gRefDBase.InitializeUniformRefPattern(EOned);
  gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
  gRefDBase.InitializeUniformRefPattern(ETriangle);

  // --- Solve darcy problem ---

  int order = 1; // Polynomial order

  // Initial geometric mesh
  bool isNL = true; // Non-linear flag
  TPZGeoMesh *gmesh = createGeoMesh3D({3, 3, 3}, {0., 0., 0.}, {1., 1., 1.});
  // TPZGeoMesh* gmesh = createGeoMesh2D({3, 3}, {0., 0.}, {1., 1.});

  if (isNL) {
    std::cout << "Solving non-linear Darcy problem..." << std::endl;
    SolveNonLinear(order, gmesh);
  } else {
    std::cout << "Solving linear Darcy problem..." << std::endl;
    SolveLinear(order, gmesh);
  }
}

// =========
// Functions
// =========

TPZGeoMesh *createGeoMesh3D(const TPZManVector<int, 3> &nelDiv, const TPZManVector<REAL, 3> &minX, const TPZManVector<REAL, 3> &maxX) {

  TPZGenGrid3D generator(minX, maxX, nelDiv, MMeshType::EHexahedral);

  generator.BuildVolumetricElements(EMatId);
  TPZGeoMesh *gmesh = generator.BuildBoundaryElements(EBottom, ELeft, EFront, ERight, EBack, ETop);

  return gmesh;
}

TPZGeoMesh *createGeoMesh2D(const TPZManVector<int, 2> &nelDiv, const TPZManVector<REAL, 2> &minX, const TPZManVector<REAL, 2> &maxX) {

  TPZGeoMesh *gmesh = new TPZGeoMesh;
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  generator.SetElementType(MMeshType::EQuadrilateral);
  generator.Read(gmesh, EMatId);
  generator.SetBC(gmesh, 4, EFront);
  generator.SetBC(gmesh, 5, ERight);
  generator.SetBC(gmesh, 6, EBack);
  generator.SetBC(gmesh, 7, ELeft);

  return gmesh;
}

TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, int order, bool isNL) {

  const int dim = gmesh->Dimension();

  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.SetDefaultOrder(order);
  hdivCreator.SetShouldCondense(true);

  // Add materials (weak formulation)
  TPZMixedDarcyFlow *matDarcy = nullptr;
  if (isNL) {
    TSFMixedDarcy *temp_mat = new TSFMixedDarcy(EMatId, gmesh->Dimension());
    temp_mat->SetConstantPermeability(1.0);
    matDarcy = temp_mat;
  } else {
    matDarcy->SetConstantPermeability(1.0);
  }
  hdivCreator.InsertMaterialObject(matDarcy);

  // Create, set and add boundary conditions
  TPZManVector<REAL, 1> val2(1, 0.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Part that goes to the Stiffnes matrix
  val2[0] = 0.;
  auto bcond = matDarcy->CreateBC(matDarcy, ERight, 1, val1, val2);
  hdivCreator.InsertMaterialObject(bcond);

  bcond = matDarcy->CreateBC(matDarcy, ELeft, 1, val1, val2);
  hdivCreator.InsertMaterialObject(bcond);

  val2[0] = 0.;
  bcond = matDarcy->CreateBC(matDarcy, EFront, 0, val1, val2);
  hdivCreator.InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = matDarcy->CreateBC(matDarcy, EBack, 0, val1, val2);
  hdivCreator.InsertMaterialObject(bcond);

  if (gmesh->Dimension() == 3) {
    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, EBottom, 1, val1, val2);
    hdivCreator.InsertMaterialObject(bcond);

    bcond = matDarcy->CreateBC(matDarcy, ETop, 1, val1, val2);
    hdivCreator.InsertMaterialObject(bcond);
  }

  int lagmultilevel = 1;
  TPZManVector<TPZCompMesh *, 7> meshvec(hdivCreator.NumMeshes());
  hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
  TPZMultiphysicsCompMesh *cmesh = nullptr;
  hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);

  {
    std::ofstream out("cmeshMP.txt");
    cmesh->Print(out);
  }

  return cmesh;
}

void SolveNonLinear(int order, TPZGeoMesh *gmesh) {

  TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, order, true);

  // Mixed solver
  TPZLinearAnalysis anMixed(cmeshMixed, RenumType::EMetis);
#ifdef PZ_USING_MKL
  TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
  TPZFStructMatrix<STATE> matMixed(cmeshMixed);
#endif
  matMixed.SetNumThreads(nthreads);
  ApplyEquationFilter(gmesh, cmeshMixed, matMixed);
  anMixed.SetStructuralMatrix(matMixed);
  TPZStepSolver<STATE> stepMixed;
  stepMixed.SetDirect(ELDLt);
  anMixed.SetSolver(stepMixed);

  const int maxIter = 10;
  const REAL tol = 1.e-6;
  TPZFMatrix<STATE> Sol(cmeshMixed->NEquations(), 1, 0.);
  for (int iteration = 0; iteration < maxIter; iteration++) {
    std::cout << "Non-linear iteration " << iteration << std::endl;
    anMixed.Assemble();
    // check convergence
    if (iteration > 0) {
      TPZMatrix<STATE> &rhs = anMixed.Rhs();
      REAL norm = 0.;
      for (int i = 0; i < rhs.Rows(); i++) {
        norm += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
      }
      norm = sqrt(norm);
      std::cout << "norm = " << norm << std::endl;
      if (norm < tol) {
        std::cout << "Converged!" << std::endl;
        break;
      }
    }
    anMixed.Solve();
    TPZMatrix<STATE> &dsol = anMixed.Solution();
    Sol += dsol;
    // anMixed.LoadSolution(Sol);
    cmeshMixed->LoadSolution(Sol);
    cmeshMixed->TransferMultiphysicsSolution();
  }

  // ---- Plotting ---

  {
    const std::string plotfile = "darcy_mixed";
    constexpr int vtkRes{0};
    TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
    auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
    vtk.Do();
  }

  // --- Clean up ---

  delete cmeshMixed;
}

void SolveLinear(int order, TPZGeoMesh *gmesh) {
  TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, order, false);

  // Mixed solver
  TPZLinearAnalysis anMixed(cmeshMixed, RenumType::EMetis);
#ifdef PZ_USING_MKL
  TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
  TPZFStructMatrix<STATE> matMixed(cmeshMixed);
#endif
  matMixed.SetNumThreads(nthreads);
  anMixed.SetStructuralMatrix(matMixed);
  TPZStepSolver<STATE> stepMixed;
  stepMixed.SetDirect(ELDLt);
  anMixed.SetSolver(stepMixed);
  anMixed.Run();

  // ---- Plotting ---

  {
    const std::string plotfile = "darcy_mixed";
    constexpr int vtkRes{0};
    TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
    auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
    vtk.Do();
  }

  // --- Clean up ---

  delete cmeshMixed;
}

void ApplyEquationFilter(TPZGeoMesh *gmesh, TPZCompMesh *cmesh_m, TPZStructMatrixT<STATE> &strmat)
{
    std::set<int64_t> removeConnectSeq;
    std::set<int64_t> removeEquations;

    cmesh_m->LoadReferences();
    
    for (auto el : gmesh->ElementVec())
    {        
        int elMatID = el->MaterialId();
        
        if (elMatID != ERight && elMatID != ELeft && elMatID != EBottom && elMatID != ETop)
            continue;

        TPZCompEl *compEl = el->Reference();
        
        int64_t nConnects = compEl->NConnects();
        
        if (nConnects != 1)
            DebugStop();
        
        removeConnectSeq.insert(compEl->Connect(0).SequenceNumber());
    }
    
    for (auto blockNumber : removeConnectSeq)
    {
        auto firstEq = cmesh_m->Block().Position(blockNumber);
        
        int64_t blockSize = cmesh_m->Block().Size(blockNumber);
        
        for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++)
            removeEquations.insert(eq);
    }

    TPZEquationFilter filter(cmesh_m->NEquations());
    filter.ExcludeEquations(removeEquations);
    strmat.EquationFilter() = filter;
}