#ifndef SOURCES_H
#define SOURCES_H

#include "sources.h"

std::set<int> vugBcIds;
std::set<int> vugIds;

std::set<int> fracBcIds;
std::set<int> fracIds;

TPZGeoMesh* generateGMeshWithPhysTagVec(ReadJson inputData, std::string filename, std::string meshName){

    TPZGmshReader GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4);

    std::vector<DomData> domainData = inputData.DomainData();
    std::vector<BcData> bcData = inputData.BCInput();
    std::vector<BcData> FracbcData = inputData.FracBCInput();
    int Dim = inputData.dim();

    for (auto bc: bcData){
        int elDim = Dim-1;
        dim_name_and_physical_tag[elDim][bc.name] = bc.matId;
    }

    for (auto bc: FracbcData){
        int elDim = 0;
        dim_name_and_physical_tag[elDim][bc.name] = bc.matId;
    }

    for (auto info: domainData){
        int elDim = Dim;
        //if (info.matId == 5) elDim = 1;
        dim_name_and_physical_tag[elDim][info.name] = info.matId;
    }

    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tag);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);
    gmeshFine->SetName(meshName);
    return gmeshFine;
}

void RefineElement(TPZGeoMesh *gmesh, int matIdtoFind, int maxlevel){

    int nels = gmesh->NElements();
    std::set<int> matids2 = {matIdtoFind};
  
    for(int el = 0; el < nels; el++){
        TPZGeoEl *gel = gmesh->Element(el);

        if(gel->MaterialId() == matIdtoFind){
            for(int il = 1; il <= maxlevel; il++){
                TPZGeoElSide gelside(gel);
                TPZGeoElSide neigh = gelside.Neighbour();
                while(neigh != gelside){
                    TPZGeoEl *neighgel = neigh.Element();
                    if(neighgel->Level() < il){
                        TPZRefPatternTools::RefineDirectional(neigh.Element(), matids2);
                    }
                    neigh = neigh.Neighbour();
                }
            }
        }
    }
}

void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &elsId, std::set<int> &bcId){
    int dim = geomesh->Dimension();
    for (auto gel: geomesh->ElementVec()) {
        if (! gel) {
            continue;
        }
        int matId = gel->MaterialId();
        int geldim = gel->Dimension();
        if (geldim == dim || (fracIds.find(matId) != fracIds.end())) {
            elsId.insert(matId);
        }
        if ((geldim == dim-1) && (fracIds.find(matId) == fracIds.end())) {
            bcId.insert(matId);
        }
    }
}

TPZCompMesh *CreateMesh(TPZGeoMesh* gmesh, ReadJson inputData){

    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
    cmesh->SetName("CompMesh");
    cmesh->SetDefaultOrder(2);

    std::vector<DomData> domainData = inputData.DomainData();
    std::map<int,double> matId_perm;
    for(auto info: domainData){
        matId_perm[info.matId] = info.permeability;
    }

    std::vector<BcData> bcData = inputData.BCInput();
    int bcType = 0;
    int bcId = 0;
    TPZManVector<double, 3> bcValue;
    
    // Add materials (weak formulation)
    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(EMatId, inputData.dim());
    matDarcy->SetNonLinearContext(false);
    matDarcy->SetConstantPermeability(matId_perm.at(EMatId));
    cmesh->InsertMaterialObject(matDarcy);

    TPZDarcyFlow *matDarcyFrac = nullptr;
    for(auto Id: fracIds) {
        matDarcyFrac = new TPZDarcyFlow(Id, 1);
        matDarcyFrac->SetNonLinearContext(true);
        if(inputData.problemType() == 0){
            matDarcyFrac->SetConstantPermeability(matId_perm.at(5));
        }
        cmesh->InsertMaterialObject(matDarcyFrac);
    }

    TPZDarcyFlow *matDarcyVug = nullptr;
    for(auto Id: vugIds) {
        matDarcyVug = new TPZDarcyFlow(Id, inputData.dim());
        matDarcyVug->SetNonLinearContext(true);
        if(inputData.problemType() == 0) {
            matDarcyVug->SetConstantPermeability(matId_perm.at(6));
        }
        cmesh->InsertMaterialObject(matDarcyVug);
    }

    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);

    for (auto bc: bcData){
        bcId = bc.matId;
        bcType = bc.type;
        bcValue = bc.value;

        val2[0] = bcValue[0];
        TPZBndCond * face = matDarcy->CreateBC(matDarcy, bcId, bcType, val1, val2);
        cmesh->InsertMaterialObject(face);
    }

    val2[0] = 0;
    // for(auto bcId: fracBcIds) { // Frac boundary elements
    //     TPZBndCond *faceFrac = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
    //     cmesh->InsertMaterialObject(faceFrac);
    // }

    // for(auto bcId: vugBcIds) { // Vug boundary elements 
    //     TPZBndCond *faceVug = matDarcyVug->CreateBC(matDarcyVug, bcId, 0, val1, val2);
    //     cmesh->InsertMaterialObject(faceVug);
    // }

    cmesh->AutoBuild();

    if(inputData.problemType() == 1) SetUniqueVugConnect(gmesh, cmesh);

    cmesh->CleanUpUnconnectedNodes();

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    cmesh->ExpandSolution();

    return cmesh;
}

void PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting comp mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}


void PrintGeoMesh(TPZGeoMesh *gmesh)
{
    std::cout << "\nPrinting geo mesh in .txt and .vtk formats...\n";

    std::ofstream VTKGeoMeshFile(gmesh->Name() + ".vtk");
    std::ofstream TextGeoMeshFile(gmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, VTKGeoMeshFile);
    gmesh->Print(TextGeoMeshFile);
}

void PrintElement(TPZCompMesh *cmesh, int matId){

    int nel = cmesh->NElements();
    std::ofstream fileElements("FractureEls.txt", std::ios::app);

    for(int el = 0; el < nel; el++){
        TPZCompEl* cel = cmesh->Element(el);
        TPZGeoEl* gel = cel->Reference();
        if(gel->MaterialId() == matId){
            gel->Print(fileElements);
            fileElements << "\n\n";
            cel->Print(fileElements);
            fileElements << "\n\n";
        }
    }
    fileElements << "------------------------------------------------------------------------------------------------\n\n";
}

void GetCompEls(TPZGeoMesh* gmesh, TPZCompMesh *cmeshH1, TPZCompMesh *cmeshHdiv, TPZVec<TPZCompEl*> &celVecH1, TPZVec<TPZCompEl*> &celVecHdiv){

    int64_t nels = gmesh->NElements();
    int64_t nels1 = cmeshH1->NElements();
    int64_t nels2 = cmeshHdiv->NElements();
    // celVecH1.Resize(nels);
    // celVecHdiv.Resize(nels);
    

    for(int64_t el = 0; el < nels1; el++){

        TPZCompEl* celH1 = cmeshH1->Element(el);
        TPZGeoEl *gel = celH1->Reference();
        if(!gel) DebugStop();
        int64_t gelId = gel->Index();

        celVecH1[gelId] = celH1;
    }

    for(int64_t el = 0; el < nels2; el++){

        TPZCompEl* celHdiv = cmeshHdiv->Element(el);
        TPZGeoEl *gel = celHdiv->Reference();
        if(!gel) DebugStop();
        int64_t gelId = gel->Index();

        celVecHdiv[gelId] = celHdiv;
    }

}

REAL ComputeErrorH1Hdiv(TPZVec<TPZCompEl*> &celVecH1, TPZVec<TPZCompEl*> &celVecHdiv, std::set<int> matId, TPZFMatrix<STATE> &elSolMat){

    REAL error = 0;
    REAL elerror = 0;
    int nels = celVecH1.size();

    for(int el = 0; el < nels; el++){

        if(!celVecH1[el] || !celVecHdiv[el]) continue;

        int elMatId = celVecH1[el]->Reference()->MaterialId();

        if (matId.find(elMatId) == matId.end()) continue;

        int64_t H1index = celVecH1[el]->Index();

        elerror = CalcElementError(celVecH1[el], celVecHdiv[el]);
        //elerror = CalcEnergy(celVecH1[el], celVecHdiv[el]);

        elSolMat(H1index,0) = std::sqrt(elerror);

        error += elerror;
    }
    return std::sqrt(error);
}

REAL CalcElementError(TPZCompEl* celH1, TPZCompEl* celHdiv){

    REAL result = 0.;
    if(!celH1 || !celHdiv) {
        std::cout << "No computational element found!\n";
        return result;
    }

    int dim = celH1->Dimension();
    int matId = celH1->Reference()->MaterialId();

 
    TPZMaterial *matH1 = celH1->Material();
    TPZMaterial *matHdiv = celHdiv->Material();

    //TPZMaterial *mat = cmeshH1->FindMaterial(EMatId);
    TPZDarcyFlow *matDarcy = dynamic_cast<TPZDarcyFlow*>(matH1); 
    //downcasting—converting a base class pointer to a derived class pointer

    //const TPZIntPoints intrule = celH1->GetIntegrationRule(); //! Ask
    TPZGeoEl *gel = celH1->Reference();
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides()-1,4); //! not sure order
    int npoints = intrule->NPoints();

    TPZVec<REAL> xi(dim, 0.0);
    TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3); //jacobian
    REAL detjac;

    TPZVec<STATE> solH1(dim, 0.0);
    TPZVec<STATE> solH1x(1, 0.0);
    TPZVec<STATE> solH1y(1, 0.0);
    TPZVec<STATE> solHdiv(dim, 0.0);
    TPZFMatrix<STATE> Perm(dim, dim, 0.0), InvPerm(dim, dim, 0.0); 

    // Performing numerical integration
    for (int point = 0; point < npoints; point++) {
        REAL weight;
        intrule->Point(point, xi, weight);
        gel->Jacobian(xi, jac, axes, detjac, jacinv);

        celH1->Solution(xi, 2, solH1); //GradP
        celH1->Solution(xi, 3, solH1x); //kGradP
        celH1->Solution(xi, 4, solH1y); //kGradP
        celHdiv->Solution(xi, 1, solHdiv); //Flux

        const STATE perm = matDarcy->GetPermeability(xi);
        const STATE inv_perm = 1 / perm;
        for (int i = 0; i < dim; i++) {
            Perm(i, i) = perm;
            InvPerm(i, i) = inv_perm; //! Check
        }

        TPZFMatrix<STATE> fluxH1(dim, 1, 0.0);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                fluxH1(i,0) += Perm(i,j)*solH1[j];
            }
        }

        // fluxH1(0,0) = solH1x[0];
        // if (dim > 1) fluxH1(1,0) = solH1y[0];

        REAL aux = 0.;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                aux += (fluxH1(i,0)+solHdiv[i])*InvPerm(i,j)*(fluxH1(j,0)+solHdiv[j]);
            }
        }

        result += aux*weight*fabs(detjac);
    }
    delete intrule;
    return result;
}

REAL CalcEnergy(TPZCompEl* celH1, TPZCompEl* celHdiv){

    REAL result = 0.;
    if(!celH1 || !celHdiv) {
        std::cout << "No computational element found!\n";
        return result;
    }

    int dim = celH1->Dimension();
    int matId = celH1->Reference()->MaterialId();

 
    TPZMaterial *matH1 = celH1->Material();
    TPZMaterial *matHdiv = celHdiv->Material();

    //TPZMaterial *mat = cmeshH1->FindMaterial(EMatId);
    TPZDarcyFlow *matDarcy = dynamic_cast<TPZDarcyFlow*>(matH1); 
    //downcasting—converting a base class pointer to a derived class pointer

    //const TPZIntPoints intrule = celH1->GetIntegrationRule(); //! Ask
    TPZGeoEl *gel = celH1->Reference();
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides()-1,4); //! not sure order
    int npoints = intrule->NPoints();

    TPZVec<REAL> xi(dim, 0.0);
    TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3); //jacobian
    REAL detjac;

    TPZVec<STATE> solH1(dim, 0.0);
    TPZVec<STATE> solH1x(1, 0.0);
    TPZVec<STATE> solH1y(1, 0.0);
    TPZVec<STATE> solHdiv(dim, 0.0);
    TPZFMatrix<STATE> Perm(dim, dim, 0.0), InvPerm(dim, dim, 0.0); 

    // Performing numerical integration
    for (int point = 0; point < npoints; point++) {
        REAL weight;
        intrule->Point(point, xi, weight);
        gel->Jacobian(xi, jac, axes, detjac, jacinv);

        celH1->Solution(xi, 2, solH1); //GradP
        celH1->Solution(xi, 3, solH1x); //kGradP
        celH1->Solution(xi, 4, solH1y); //kGradP
        celHdiv->Solution(xi, 1, solHdiv); //Flux

        const STATE perm = matDarcy->GetPermeability(xi);
        const STATE inv_perm = 1 / perm;
        for (int i = 0; i < dim; i++) {
            Perm(i, i) = perm;
            InvPerm(i, i) = inv_perm; //! Check
        }

        TPZFMatrix<STATE> fluxH1(dim, 1, 0.0);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                fluxH1(i,0) += Perm(i,j)*solH1[j];
            }
        }

        // fluxH1(0,0) = solH1x[0];
        // if (dim > 1) fluxH1(1,0) = solH1y[0];

        REAL aux = 0.;
        // for (int i = 0; i < dim; i++) {
        //     for (int j = 0; j < dim; j++) {
        //         aux += (solHdiv[i])*InvPerm(i,j)*(solHdiv[j]);
        //     }
        // }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                aux += solH1[i]*Perm(i,j)*solH1[j];
            }
        }

        result += aux*weight*fabs(detjac);
    }
    delete intrule;
    return result;
}

void SetAnalysis(TPZLinearAnalysis* an, TPZCompMesh* cmesh, ReadJson inputData){

    // an->LoadSolution();


    #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmesh);
    #else
    TPZFStructMatrix<STATE> matMixed(cmesh);
    #endif
    matMixed.SetNumThreads(4);
    an->SetStructuralMatrix(matMixed);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    an->Assemble();
    an->Solve();
}

ReadJson::ReadJson(std::string fileName){
    
    std::ifstream filejson(fileName);

    fInputFile = json::parse(filejson,nullptr,true,true,true); 

    fMeshName = fInputFile["MeshName"];
    
    fMeshDirectory = fInputFile["MeshDirectory"];
    
    fApproxType = fInputFile["Simulation"]["ApproximationType"];

    fProblemType = fInputFile["Simulation"]["ProblemType"];
    
    fPresspOrder = fInputFile["PressurepOrder"];

    fDim = fInputFile["Dimension"];
    
    fResolution = fInputFile["Resolution"];

    // fVisc = fInputFile["ProblemData"]["FluidViscosity"];

    DomData domInput;

    for(auto& domain : fInputFile["DomainData"]){
        if(domain.find("Name") == domain.end()) DebugStop(); // check if the information exists

        domInput.name = domain["Name"];
        domInput.matId = domain["matId"];
        domInput.permeability = domain["Permeability"];

        fDomainDataVec.push_back(domInput);
    }

    BcData bcInput;

    for(auto& bcjson : fInputFile["BCs"]){
        if(bcjson.find("Name") == bcjson.end()) DebugStop(); // check if the information exists

        bcInput.name = bcjson["Name"];
        bcInput.matId = bcjson["matId"];
        bcInput.type = bcjson["Type"];
        bcInput.value[0] = bcjson["Value"];
        
        fBcDataVec.push_back(bcInput);
    }

    BcData fracBcInput;

    for(auto& bcjson : fInputFile["FracBCs"]){
        if(bcjson.find("Name") == bcjson.end()) DebugStop(); // check if the information exists

        fracBcInput.name = bcjson["Name"];
        fracBcInput.matId = bcjson["matId"];
        fracBcInput.type = bcjson["Type"];
        fracBcInput.value[0] = bcjson["Value"];
        
        fFracBcVec.push_back(fracBcInput);
    }
}
    
std::string ReadJson::MeshName(){
    return fMeshName;
}

std::string ReadJson::MeshFile(){

    std::string meshPath = fMeshDirectory + fMeshName + ".msh";
    return meshPath;
}

std::vector<DomData> ReadJson::DomainData(){
    return fDomainDataVec;
}

int ReadJson::approxType(){
    return fApproxType;
}

int ReadJson::problemType(){
    return fProblemType;
}

int ReadJson::pressOrder(){
    return fPresspOrder;
}

int ReadJson::dim(){
    return fDim;
}

int ReadJson::resolution(){
    return fResolution;
}
    
// double ReadJson::visc(){
//     return fVisc;
// }

std::vector<BcData> ReadJson::BCInput(){
    return fBcDataVec;
}

std::vector<BcData> ReadJson::FracBCInput(){
    return fFracBcVec;
}

#endif
