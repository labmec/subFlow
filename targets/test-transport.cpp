
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "TPZEquationFilter.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZGmshReader.h"
#include "TPZHDivApproxCreator.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "TPZRefPatternDataBase.h"
#include "TPZSSpStructMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "TSFApproxCreator.h"
#include "TSFDarcyAnalysis.h"
#include "TSFMixedDarcy.h"
#include "TSFProblemData.h"
#include "pzfstrmatrix.h"
#include "pzintel.h"
#include "pzlog.h"
#include "pzmultiphysicselement.h"
#include "pzstepsolver.h"
#include "pzvec_extras.h"
#include <iostream>

TPZGeoMesh *ReadMeshFromGmsh(TSFProblemData &simData);

int main(int argc, char *const argv[]) {

// Initialize logger
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG("subFlow-log.cfg");
#endif

  // Initial geometric mesh
  bool isNL = true; // Non-linear flag
  TSFProblemData simData;
  simData.ReadJSONFile("test-transport.json");
  TPZGeoMesh *gmesh = ReadMeshFromGmsh(simData);
  {
    std::ofstream out("gmesh-before-interfaces.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  TSFApproxCreator approxCreator(gmesh);
  approxCreator.SetProblemData(&simData);
  approxCreator.ConfigureDarcySpace();
  approxCreator.AddDarcyMaterials();
  TPZMultiphysicsCompMesh *darcy_cmesh = approxCreator.CreateApproximationSpace();
  {
    std::ofstream out("darcy-cmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(darcy_cmesh, out);
  }

  if (simData.fTNumerics.fAnalysisType == 0) { // Darcy
    TSFDarcyAnalysis darcyAnalysis(darcy_cmesh);
    darcyAnalysis.SetProblemData(&simData);
    darcyAnalysis.Initialize();
    darcyAnalysis.RunTimeStep();
    darcyAnalysis.PostProcessTimeStep(gmesh->Dimension(), 0);
  } else {
    approxCreator.BuildTransportCmesh();
    TPZCompMesh *transport_cmesh = approxCreator.GetTransportCmesh();
    {
      std::ofstream out("transport-cmesh.vtk");
      TPZVTKGeoMesh::PrintCMeshVTK(transport_cmesh, out);
    }

    delete transport_cmesh;
  }

  delete darcy_cmesh;
  delete gmesh;

  return 0;
}

TPZGeoMesh *ReadMeshFromGmsh(TSFProblemData &simData) {
  // read mesh from gmsh
  TPZGeoMesh *gmesh;
  gmesh = new TPZGeoMesh();
  std::string file_name = std::string(INPUTDIR) + "/" + simData.fTGeometry.fGmshFile;
  {
    TPZGmshReader reader;
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}