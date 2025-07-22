/**
 * @file TestOneFrac.cpp
 * @brief Define a Unit Test for a one fracture problem
 *
 */

//#define USEMAIN

#ifndef USEMAIN
#include <catch2/catch.hpp>
#endif

#include <TMRSApproxSpaceGenerator.h>
#include <tpzgeoelrefpattern.h>
#include <TPZGenGrid3D.h>
#include "imrs_config.h"
#include "pzlog.h"

// ---- Enum for materials ----
enum EMatid {ENone, EDomain, EInlet, EOutlet, ENoflux, EPressure, EIntersection, EIntersectionEnd, EVolume, EFaceBCPressure};
int globFracID = 10;

// ---- Functions ----
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase1(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase2(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase3(std::string &filename);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void ChangeBCsToNoFlux(TPZGeoMesh* gmesh);
void ExtractAllCompEls(TPZElementGroup* elgr, std::list<TPZCompEl*>& ellist);
void GetListOfCompEls(TPZCompEl* cel, std::list<TPZCompEl*>& ellist);
const STATE GetPressureAtCenter(TPZCompEl* cel, TPZVec<REAL>& xcent);

// ---- Driver Function ----
namespace onefractest{
    void TestOneFrac(const int& caseToSim);
};

using namespace std;

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

#ifndef USEMAIN
// ---- Test 0 ----
TEST_CASE("constant_pressure","[onefrac_test]"){
    onefractest::TestOneFrac(0);
}

// ---- Test 1 ----
TEST_CASE("linear_pressure","[onefrac_test]"){
    onefractest::TestOneFrac(1);
}

// ---- Test 2 ----
TEST_CASE("constant_pressure_frac_at_dom_bound","[onefrac_test]"){
    onefractest::TestOneFrac(2);
}

// ---- Test 3 ----
TEST_CASE("linear_pressure_transport","[onefrac_test]"){
    onefractest::TestOneFrac(3);
}

#else
int main(){
    onefractest::TestOneFrac(3);
    return 0;
}
#endif

// ---- Driver Function Implementation ----
void onefractest::TestOneFrac(const int& caseToSim) {
    
    string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
    string logpath = basemeshpath + "/../DFNIMRS/log4cxx.cfg";
    TPZLogger::InitializePZLOG(logpath);
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif

    
    cout << "\n\n\t\t---------------------- Start of Simulation " << caseToSim <<  " ------------------------\n" << endl;
    
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmesh = nullptr;
    TMRSDataTransfer sim_data;
    CreateGMeshAndDataTransfer(gmesh,sim_data,caseToSim);
    const bool printgmesh = true;
    if (printgmesh) {
        std::ofstream name("GeoMeshFractures.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name);
    }
    
    // ----- Approximation space -----
    TMRSApproxSpaceGenerator aspace;
	//    sim_data.mTFracProperties.m_matid = globFracID; // This is now set for each fracture
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    
    // ----- Setting gmesh -----
    aspace.InitMatIdForMergeMeshes() = 8;
    aspace.SetGeometry(gmesh,gmesh);
    
    // ----- Setting the global data transfer -----
//    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
    aspace.SetDataTransfer(sim_data);    
    
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    REAL mass =0.0;
    if(caseToSim == 3){
        aspace.BuildAuxTransportCmesh();
        TPZCompMesh * transport_operator = aspace.GetTransportOperator();
        TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
        sfi_analysis->SetDataTransferAndBuildAlgDatStruct(&sim_data);
        sfi_analysis->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
        int n_steps = sim_data.mTNumerics.m_n_steps;
        REAL dt = sim_data.mTNumerics.m_dt;
        
        TPZStack<REAL,100> reporting_times;
        reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
        REAL sim_time = 0.0;
        int pos =0;
        REAL current_report_time = reporting_times[pos];
        int npos = reporting_times.size();
        
        sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
        REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        std::cout << "Mass report at time : " << 0.0 << std::endl;
        std::cout << "Mass integral :  " << initial_mass << std::endl;
        std::ofstream fileCilamce("IntegratedSat.txt");
        TPZFastCondensedElement::fSkipLoadSolution = false;
        bool first=true;
        for (int it = 1; it <= n_steps; it++) {
            sim_time = it*dt;
            sfi_analysis->m_transport_module->SetCurrentTime(dt);
            const int typeToPostProc = 2; // only saturation
            sfi_analysis->PostProcessTimeStep(2);
            sfi_analysis->RunTimeStep();
            mixed_operator->LoadSolution(mixed_operator->Solution());
            if (sim_time >=  current_report_time) {
                std::cout << "Time step number:  " << it << std::endl;
                std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
                mixed_operator->UpdatePreviousState(-1.);
                const int typeToPostProcEnd = 0; // p/flux and transport
                sfi_analysis->PostProcessTimeStep(typeToPostProcEnd);
                pos++;
                current_report_time =reporting_times[pos];
                REAL InntMassFrac=sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(10);
                fileCilamce<<current_report_time/(86400*365)<<", "<<InntMassFrac<<std::endl;
                mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
                std::cout << "Mass report at time : " << sim_time << std::endl;
                std::cout << "Mass integral :  " << mass << std::endl;
            }
        }
    }
    else{
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // ----- Post processing -----
    // Just doing this part in case the code crashes during it
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    const int dimToPost = 3;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    }
    // ----- Checking if results match -----
#ifndef USEMAIN
    gmesh->ResetReference();
    mixed_operator->LoadReferences();
    std::list<TPZCompEl *> ellist;
    ellist.clear();
    for (auto cel : mixed_operator->ElementVec()) {
        ellist.clear();
        if (!cel) continue;
        GetListOfCompEls(cel,ellist); // in case it is condensedcompel
        for (auto celgr : ellist) {
            TPZGeoEl* gel = celgr->Reference();
            if (!gel) continue;
            if (gel->Dimension() != 3) continue;
            TPZManVector<REAL,3> xcent(3,0.);
            const STATE pInCentOfEl = GetPressureAtCenter(celgr,xcent);
            if (caseToSim == 0 || caseToSim == 2) {
                REQUIRE( pInCentOfEl == Approx( 1.0 ) );
            }
            else if (caseToSim == 1){
                const STATE pexact = 2. - (xcent[2]+1);
                Approx targetCase1 = Approx(pexact).epsilon(1.e-3); //from catch2 lib
                REQUIRE( pInCentOfEl == targetCase1 );
            }
            else if (caseToSim == 3){
                const STATE pexact = 2. - (xcent[2]+1);
                Approx targetCase1 = Approx(pexact).epsilon(1.e-3); //from catch2 lib
                REQUIRE( pInCentOfEl == targetCase1 );
                const STATE massExact = 1.67997;
                Approx targetCaseSat = Approx(massExact).epsilon(1.e-3); //from catch2 lib
                REQUIRE( mass == targetCaseSat );
            }
            else{
            }
        } // ellist
    } // elementvec
#endif
   
    // ----- Cleaning up -----
    delete gmesh;
    
    cout << "\n\n\t\t****************** End of Simulation " << caseToSim <<  " ******************\n" << endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

const STATE GetPressureAtCenter(TPZCompEl* cel, TPZVec<REAL>& xcent) {
    TPZGeoEl *gel = cel->Reference();
    TPZManVector<REAL,3> qsicent(3,0.);
    TPZManVector<REAL,1> sol(1);
    gel->CenterPoint(gel->NSides()-1, qsicent);
    gel->X(qsicent, xcent);
    const int varp = cel->Material()->VariableIndex("Pressure");
    cel->Solution(qsicent, varp, sol);
    const STATE pInCentOfEl = sol[0];
    return pInCentOfEl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void GetListOfCompEls(TPZCompEl* cel, std::list<TPZCompEl*>& ellist) {
    TPZCondensedCompEl* condcel = dynamic_cast<TPZCondensedCompEl*>(cel);
    if (condcel) {
        cel = condcel->ReferenceCompEl();
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (!elgr){
            ellist.push_back(cel);
        }
        else{
            ExtractAllCompEls(elgr,ellist);
        }
    }
    else{
        ellist.push_back(cel);
    }

}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void ExtractAllCompEls(TPZElementGroup* elgr, std::list<TPZCompEl*>& ellist) {
    const TPZVec<TPZCompEl *>& elvec = elgr->GetElGroup();
    for (auto grcel : elvec) {
        TPZCondensedCompEl* condcel = dynamic_cast<TPZCondensedCompEl*>(grcel);
        if (condcel) {
            grcel = condcel->ReferenceCompEl();
            TPZElementGroup *elgrloc = dynamic_cast<TPZElementGroup *> (grcel);
            if (!elgrloc){
                ellist.push_back(grcel);
            }
            else{
                ExtractAllCompEls(elgrloc,ellist);
            }
        }
        else{
            ellist.push_back(grcel);
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim) {
    string basemeshpath(FRACMESHES);
    std::string filename = basemeshpath + "/2DMeshes/1fracNoBnd.msh";
    switch (caseToSim) {
        case 0: {
            gmesh = ReadFractureMeshCase0(filename);
            
        }
            break;
        case 1: {
            gmesh = ReadFractureMeshCase1(filename);
        }
            break;
        case 2: {
            std::string filename2 = basemeshpath + "/verifications/1frac2el.msh";
            gmesh = ReadFractureMeshCase2(filename2);
        }
            break;
        case 3: {
            std::string filename2 = basemeshpath + "/verifications/1frac2el.msh";
            gmesh = ReadFractureMeshCase3(filename2);
        }
            break;
        default:
            DebugStop();
    }
            
    sim_data  = SettingFracturesSimple(caseToSim);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // domain
    std::string volbase = "c";

    // Volumes
    for (int ivol = 1; ivol < 28; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
            
    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer SettingFracturesSimple(const int caseToSim){
    
    // Fracture material
    TMRSDataTransfer sim_data;
    
    int D_Type = 0;
    int N_Type = 1;
    REAL pressure_in = 1.0 ;
    REAL zero_flux = 0.0;
    
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fractures"] = globFracID;

    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
    
    // Boundary conditions
    if (caseToSim < 2) {
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,1.);

		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EPressure] = std::make_pair(N_Type, zero_flux);
    
    }
    else if (caseToSim < 3){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,1.);
        sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,1.);
        sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,1.);
        sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(D_Type,1.);
        
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EPressure] = std::make_pair(D_Type, 1.);
    }
    else if (caseToSim < 4){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
		        
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EPressure] = std::make_pair(N_Type, zero_flux);
                
		sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EInlet] = std::make_pair(D_Type, 1.);
		sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EOutlet] = std::make_pair(D_Type, 0.);
		sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[ENoflux] = std::make_pair(5, zero_flux);
//		sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EFaceBCPressure] = std::make_pair(5, 1.);
        
    }
    else {
        DebugStop();
    }
    
    // Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.1;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 1;
    sim_data.mTNumerics.m_max_iter_sfi = 1;
    
    // BorderElementOrder
    sim_data.mTNumerics.m_MortarBorderElementPresOrder=1;
    sim_data.mTNumerics.m_MortarBorderElementFluxOrder=1;
    
    // Other properties
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    sim_data.mTNumerics.m_n_steps = 10;
    sim_data.mTNumerics.m_dt      = 0.5;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
    // Fracture permeability
//    sim_data.mTFracProperties.m_Permeability = 1.e4;
	TMRSDataTransfer::TFracProperties::FracProp fracprop;
	fracprop.m_perm = 1.e4;
	fracprop.m_width = 1.;
	fracprop.m_fracbc.insert(EPressure);
	fracprop.m_fracIntersectMatID = EIntersection;
	sim_data.mTFracProperties.m_fracprops[globFracID] = fracprop;

    
    //FracAndReservoirProperties
    REAL kappa=1.0;
    int  id1=EVolume;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa;
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    REAL resPorosity = 0.2;
    REAL fracPorosity = 0.1;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(2);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[0] = std::make_tuple(EVolume, resPorosity,1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[1] = std::make_tuple(globFracID, fracPorosity, 0.2);
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    
    return sim_data;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase1(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // domain
    std::string volbase = "c";

    // Volumes
    for (int ivol = 1; ivol < 28; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
        
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    ChangeBCsToNoFlux(gmesh);
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase2(std::string &filename){
    
    const bool fromgmesh = 0;
    TPZGeoMesh* gmesh = nullptr;
    if (fromgmesh) {
        
        TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
        
        // domain
        std::string volbase = "c";
        
        // Volumes
        for (int ivol = 1; ivol < 3; ivol++) {
            std::string ivolstring = volbase + to_string(ivol);
            dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
        }
        dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
        
        // Fractures
        dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;
        
        // Fractures BCs
        dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
        
        gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    }
    else{
        // ----- Create Geo Mesh -----
        const TPZVec<REAL> minX = {-1.,-1.,-1.};
        const TPZVec<REAL> maxX = {1.,1.,1.};
        const TPZVec<int> nelDiv = {1,1,2};
        const MMeshType elType = MMeshType::EHexahedral;

        TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
        gmesh = gen3d.BuildVolumetricElements(EVolume);

        // ----- Fracture element -----
        int64_t index;
        TPZManVector<int64_t,2> nodesId = {4,5};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        nodesId = {5,7};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        nodesId = {7,6};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        nodesId = {6,4};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        TPZManVector<int64_t,4> nodesIdVec = {4,6,7,5};
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);

        // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
        gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
        
        gmesh->BuildConnectivity();
        std::ofstream out("meshbad.txt");
        gmesh->Print(out);
    }
    
    return gmesh;
}
TPZGeoMesh *ReadFractureMeshCase3(std::string &filename){
    
    const bool fromgmesh = 0;
    TPZGeoMesh* gmesh = nullptr;
   
        // ----- Create Geo Mesh -----
        const TPZVec<REAL> minX = {-1.,-1.,-1.};
        const TPZVec<REAL> maxX = {1.,1.,1.};
        const TPZVec<int> nelDiv = {1,1,2};
        const MMeshType elType = MMeshType::EHexahedral;

        TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
        gmesh = gen3d.BuildVolumetricElements(EVolume);

        // ----- Fracture element -----
        int64_t index;
        TPZManVector<int64_t,2> nodesId = {4,5};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        nodesId = {5,7};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        nodesId = {7,6};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        nodesId = {6,4};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
        TPZManVector<int64_t,4> nodesIdVec = {4,6,7,5};
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);

        // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
        gmesh = gen3d.BuildBoundaryElements(EInlet, ENoflux,ENoflux , ENoflux,ENoflux , EOutlet);
        
        gmesh->BuildConnectivity();
        std::ofstream out("meshbad.txt");
        gmesh->Print(out);
    std::ofstream file("FineMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    
    return gmesh;
}
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine){
            
    // Creating gmsh reader
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);

    return gmeshFine;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ChangeBCsToNoFlux(TPZGeoMesh* gmesh) {
    for(auto gel : gmesh->ElementVec()){
        if(gel->Dimension() != 2) continue;
        if(gel->MaterialId() != EFaceBCPressure) continue;

        const int pos = 2;
        TPZManVector<REAL,2> qsicent(2,-1.);
        TPZManVector<REAL,3> cent(3,-1.);
        gel->CenterPoint(gel->NSides()-1, qsicent);
        gel->X(qsicent, cent);
        const REAL distM1 = fabs(cent[pos] + 1.);
        const REAL distP1 = fabs(cent[pos] - 1.);
        const REAL tol = 1.e-8;
        if (distM1 < tol) {
            gel->SetMaterialId(EInlet);
            continue;
        }
        if (distP1 < tol) {
            gel->SetMaterialId(EOutlet);
            continue;
        }
        gel->SetMaterialId(ENoflux);
    }
}
