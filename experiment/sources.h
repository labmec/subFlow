#include <iostream>
#include <fstream>
#include <filesystem>
#include <json.hpp>
#include <math.h>

using json = nlohmann::json;

#include "pzcmesh.h"
#include "TPZElementMatrixT.h"
#include "pzmanvector.h"
#include "TPZGeoMeshTools.h"
#include "TPZCompMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzvec.h"
#include "TPZVTKGenerator.h"
#include <fstream>
#include "TPZMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZNullMaterialCS.h"
#include "TPZNullMaterial.h"
#include "TPZAnalysis.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHDivApproxCreator.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzintel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

enum MatID{
    EMatId = 1,
    EVugId = 600, 
    EVugBcId = 100, 
    EFracId = 700, 
    EFracBcId = 200,
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcTop = 4,
    EbcBottom = 5,
    ELagrange = 6
};

struct DomData {
    std::string name = "none"; // name of the gmesh tag
    int matId = 0;
    double permeability = 0.0;
};

struct BcData {
    std::string name = "none"; // name of the bc
    int matId = 0;
    int type = 0; 
    TPZManVector<double, 3>  value = {0.0, 0.0, 0.0}; 
};

extern std::set<int> vugBcIds;
extern std::set<int> vugIds;

extern std::set<int> fracBcIds;
extern std::set<int> fracIds;

class ReadJson
{
public:

    //Constructor
    ReadJson(std::string fileName);
    
    //Methods
    std::string MeshName();

    std::string MeshFile();

    std::vector<DomData> DomainData();

    int approxType();

    int problemType();

    int pressOrder();

    int dim();

    int resolution();

    // double visc();

    std::vector<BcData> BCInput();

    std::vector<BcData> FracBCInput();


private:

    json fInputFile;
    
    std::string fMeshName;
    
    std::string fMeshDirectory;

    std::vector<DomData> fDomainDataVec;
    
    int fApproxType;

    int fProblemType;
    
    int fPresspOrder;

    int fDim;
    
    int fResolution;
    
    double fPerm;

    double fVugPerm;

    double fFracPerm;

    double fVisc;

    std::vector<BcData> fBcDataVec;

    std::vector<BcData> fFracBcVec;

};

TPZGeoMesh* generateGMeshWithPhysTagVec(ReadJson inputData, std::string filename, std::string meshName);

TPZCompMesh *CreateMesh(TPZGeoMesh* gmesh, ReadJson inputData);

void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);

void RefineElement(TPZGeoMesh *gmesh, int matIdtoFind, int maxlevel);

void GetCompEls(TPZGeoMesh* gmesh, TPZCompMesh *cmeshH1, TPZCompMesh *cmeshHdiv, TPZVec<TPZCompEl*> &celVecH1, TPZVec<TPZCompEl*> &celVecHdiv);

REAL ComputeErrorH1Hdiv(TPZVec<TPZCompEl*> &celVecH1, TPZVec<TPZCompEl*> &celVecHdiv, std::set<int> matId, TPZFMatrix<STATE> &elSolMat);

REAL CalcElementError(TPZCompEl* celH1, TPZCompEl* celHdiv);

REAL CalcEnergy(TPZCompEl* celH1, TPZCompEl* celHdiv);

void SetAnalysis(TPZLinearAnalysis* an, TPZCompMesh* cmesh, ReadJson inputData);

void PrintCompMesh(TPZCompMesh *cmesh);

void PrintGeoMesh(TPZGeoMesh *gmesh);