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
#include "TSFSFIAnalysis.h"
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
  TSFProblemData simData;
  simData.ReadJSONFile("module-radial-filtration.json");
  TPZGeoMesh *gmesh = ReadMeshFromGmsh(simData);

  // Create computational meshes
  TSFApproxCreator approxCreator(gmesh);
  approxCreator.SetProblemData(&simData);
  approxCreator.ConfigureDarcySpace();
  approxCreator.AddDarcyMaterials();
  TPZMultiphysicsCompMesh *darcy_cmesh = approxCreator.CreateApproximationSpace();
  approxCreator.BuildTransportCmesh();
  TPZCompMesh *transport_cmesh = approxCreator.GetTransportCmesh();

  // Run SFI analysis
  TSFSFIAnalysis SFIAnalysis(darcy_cmesh, transport_cmesh);
  SFIAnalysis.SetProblemData(&simData);
  SFIAnalysis.Initialize();
  SFIAnalysis.Run();

  delete transport_cmesh;
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