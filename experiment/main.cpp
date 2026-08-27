#include "sources.h"
#include "pzlog.h"
#include "pzcheckgeom.h"

void Hdiv_MixedCT(int refLevel, std::ofstream &outfile){

    ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/Obstruction.json");
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    std::string meshName = inputData.MeshName();

    std::string filename = inputData.MeshFile();

    int problemDim = inputData.dim();

    int approxType = inputData.approxType();

    int problemType = inputData.problemType();

    int pressOrder = inputData.pressOrder();

    int refUniform = 0;

    int refDir = refLevel;


    std::string approxName;

    // if(approxType == 0){
    //     approxName = "_H1";
    // }
    // else if(approxType == 1){
    //     approxName = "_Mixed";
    // }
    if(problemType == 0){
        approxName = "_0";
    }
    else if(problemType == 1){
        approxName = "_1";
    }

    gmesh = generateGMeshWithPhysTagVec(inputData, filename, meshName);
    {
        TPZCheckGeom check(gmesh);
        check.UniformRefine(refUniform);
    }

    PrintGeoMesh(gmesh);

    if(!inputData.FracBCInput().empty()){
        std::set<int> fracMatId = {inputData.FracBCInput()[0].matId};
        RefineElement(gmesh, inputData.FracBCInput()[0].matId, refDir); // Vary refinament
    }

    PrintGeoMesh(gmesh);

    //! MeshWithSegment(inputData, gmesh);

    PrintGeoMesh(gmesh);

    std::set<int> elsId, bcId;
    GetAtomicIds(gmesh, elsId, bcId);


    if(approxType == 1){
        return;
    }
    else{

        TPZCompMesh *cmesh = CreateMesh(gmesh, inputData);
        PrintCompMesh(cmesh);
    
        TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);

        SetAnalysis(Analisys, cmesh, inputData);

        {
            TPZFMatrix<STATE> &rhs = Analisys->Rhs();
            rhs.Print(std::cout);
        }
        Analisys->Assemble();
        
        
        // BCInitialSolution(Analisys, cmesh, bcId, inputData, 0);
        // ApplyEquationFilter(Analisys, cmesh, bcId, inputData, 0);
        // NewtonMethod(cmesh, 5, 1.0e-6, 1.0e-6, Analisys);

        TPZFMatrix<STATE> &rhs = Analisys->Rhs();
        Analisys->PrintVectorByElement(std::cout, rhs, 1.e-6);

        PrintCompMesh(cmesh);
        {
            const std::string plotfile = meshName + approxName;
            constexpr int vtkRes{0};
            TPZManVector<std::string, 3> fields = {"Flux", "Pressure", "GradU"};
            auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
            vtk.Do();
        }
    }
}


int main (){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    gRefDBase.InitializeRefPatterns(2);
    
    std::ofstream fileErrors("errorsH1Hdiv.txt", std::ios::app);
    fileErrors << "\nnEq " << " Error"<< std::endl;
    int refMax = 0;
    for(int ref = 0; ref <= refMax; ref++)
        Hdiv_MixedCT(ref, fileErrors);
    return 0;
}
