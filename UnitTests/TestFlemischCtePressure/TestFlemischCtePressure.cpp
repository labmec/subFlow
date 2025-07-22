// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>
#include "TPZTimer.h"

#define USEMAIN

#ifndef USEMAIN
// Unit test includes
#include <catch2/catch.hpp>
#endif


// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(const int simcase);
TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void fixPossibleMissingIntersections(TPZGeoMesh* gmesh);
const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);

void ReadMeshesFlemischCase1(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase3(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase4LF(string& filenameFine, string& filenameCoarse,
                               TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);


enum EMatid {/*0*/ENone, EVolume, EVolume2, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EFracture, EIntersection, EIntersectionEnd, EPLossAtIntersect,
    /*14 HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
};

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

#ifndef USEMAIN
// 1: Flemisch case 1
// 2: Flemisch case 2
// 3: Flemisch case 3
// 4: Flemisch case 4 with less fractures (and no overlap).
TEST_CASE("case_1","[flemisch_cte_pressure]"){
    RunProblem(1);
}
TEST_CASE("case_2","[flemisch_cte_pressure]"){
    RunProblem(2);
}
TEST_CASE("case_3","[flemisch_cte_pressure]"){
    RunProblem(3);
}
// Test 4 is not run by default because it takes too long!
//TEST_CASE("case_4","[flemisch_cte_pressure]"){
//    RunProblem(4);
//}
#else
int main(){

#ifdef PZ_LOG
    string basemeshpath(FRACMESHES);
    string logpath = basemeshpath + "/../DFNIMRS/log4cxx.cfg";
    TPZLogger::InitializePZLOG(logpath);
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif

    int simcase = 1;
    RunProblem(simcase);
    return 0;
}
// ----------------- End of Main -----------------
#endif

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
void RunProblem(const int simcase)
{
    auto start_time = std::chrono::steady_clock::now();
    
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    string filenameCoarse, filenameFine;
    string basemeshpath(FRACMESHES);
    if (simcase == 1){
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case1_coarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case1_fine_NotSoFine.msh";
        ReadMeshesFlemischCase1(filenameFine,filenameCoarse,gmeshfine,gmeshcoarse);
    }
    else if (simcase == 2){
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case2_coarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case2_fine.msh";
        ReadMeshesFlemischCase2(filenameFine,filenameCoarse,gmeshfine,gmeshcoarse);
    }
    else if (simcase == 3){
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case3_coarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case3_fine.msh";
        ReadMeshesFlemischCase3(filenameFine,filenameCoarse,gmeshfine,gmeshcoarse);
    }
    else if (simcase == 4){
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse_lf.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine_lf.msh";
        ReadMeshesFlemischCase4LF(filenameFine,filenameCoarse,gmeshfine,gmeshcoarse);
    }
    else
        DebugStop();
        
    fixPossibleMissingIntersections(gmeshfine); // read about this in the function
    
    
    TMRSDataTransfer sim_data;
	// ----- Approximation space -----
	sim_data.mTNumerics.m_four_approx_spaces_Q = true;
#ifndef USEMAIN
	sim_data.mTNumerics.m_mhm_mixed_Q = GENERATE(false,true);
#else
    sim_data.mTNumerics.m_mhm_mixed_Q = true;
#endif
	sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    
	FillDataTransfer(sim_data);
	
    // ----- Printing gmesh -----
#ifdef PZDEBUG
    const bool printgmesh = true;
    if (printgmesh) {
        if(gmeshfine){
            std::ofstream name("GeoMesh_Fine_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
        }
        if(gmeshcoarse){
            std::ofstream name("GeoMesh_Coarse_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, name);
        }
    }
#endif
    

    
    // ----- Setting gmesh -----
    // Code takes a fine and a coarse mesh to generate MHM data structure
    TMRSApproxSpaceGenerator aspace;
    aspace.InitMatIdForMergeMeshes() = EInitVolumeMatForMHM;
  //  sim_data.mTFracProperties.m_matid = EFracture;
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
	
	// ----- Setting the global data transfer -----
	aspace.SetDataTransfer(sim_data);
	
    aspace.SetGeometry(gmeshfine,gmeshcoarse);
		
//    aspace.SetGeometry(gmeshfine);
    

    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true; // This makes solving very fast!
    int n_threads = 8;
    bool UsingPzSparse = true; // Necessary to use multithread for now...
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    cout << "\n--------------------- Assembling ---------------------\n" << endl;
    cout << "Number of equations: " << mixed_operator->NEquations() << endl;
    mixAnalisys->Assemble();

    cout << "\n--------------------- Solving ---------------------\n" << endl;
    mixAnalisys->Solve();
    
    TPZFastCondensedElement::fSkipLoadSolution = false; // So we can postprocess variables correctly
    mixed_operator->LoadSolution(mixed_operator->Solution());
	
	// The system is solve as non linear, so have to multiply by -1
	mixed_operator->UpdatePreviousState(-1.);
    
    // ----- Post processing -----
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    if(simcase == 1){
        // Post processing takes to long for unit testing purposes,
        // so only doing it in the most simple case
        int dimToPost = 2;
        mixAnalisys->PostProcessTimeStep(dimToPost);
        dimToPost = 3;
        mixAnalisys->PostProcessTimeStep(dimToPost);
    }

    // ----- Compute integral of pressure and flux over domain and compare with analytical solution -----
    const std::string pvarname = "Pressure";
    const STATE integratedpressure = ComputeIntegralOverDomain(mixed_operator,pvarname);
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(mixed_operator,qvarname);
    if (fabs(integratedflux) < 1.e-11 ) integratedflux = 0.; // to make Approx(0.) work
    if(simcase == 4){ // case has a HUGE domain
        if (fabs(integratedflux) < 1.e-7 ) integratedflux = 0.; // to make Approx(0.) work
    }
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;
    
#ifndef USEMAIN
    // ----- Comparing with analytical solution -----
    // Imposing unit pressure in all domains. Then, the integrated pressure should be equal
    // to the volume of each domain
    if(simcase == 1) // domain is 100x100x100
        REQUIRE( integratedpressure == Approx( 1.e6 ) );
    else if(simcase == 2) // domain is 1x1x1
        REQUIRE( integratedpressure == Approx( 1.0 ) );
    else if(simcase == 3) // domain is 1x2.25x1
        REQUIRE( integratedpressure == Approx( 2.25 ) );
    else if(simcase == 4) // domain is 850x1400x600
        REQUIRE( integratedpressure == Approx( 7.14e8 ) );

    // Flux should always be zero since all cases have cte pressure
    REQUIRE( integratedflux == Approx( 0.) );
#endif
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()/1000.;
    cout << "\n\n\t--------- Total time of simulation = " << total_time << " seconds -------\n" << endl;

    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, unit_pressure = 1.0, zero_pressure = 0.;
    
    // Domain material
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    
    // Domain boundary conditions
    
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,unit_pressure);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,zero_pressure);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,unit_pressure);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fractures"] = EFracture;
    
    // Fracture boundary conditions
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracNoFlux] = std::make_pair(N_Type, zero_flux);

	sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
        
    //FracAndReservoirProperties
    REAL kappa=1.0;
    int  id1=EVolume;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa;
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
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
    return sim_data;
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

void ReadMeshesFlemischCase1(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    for (int ivol = 33; ivol <= 44; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc6"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 33; ivol <= 44; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-33);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc6"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux;
    
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
            
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    for (int ivol = 1; ivol <= 512; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
        
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 1; ivol <= 512; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 8; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_2"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_4"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_8"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_2"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_8"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_4"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_4"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_8"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_8"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_8"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_7_8"] = EIntersection;
            
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesFlemischCase3(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    const int initc = 229, endc = 444;
    for (int ivol = initc; ivol <= endc; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc6"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
        
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = initc; ivol <= endc; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-initc);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc6"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture12"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture13"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture14"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture15"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture16"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture17"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 7; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;

    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
    
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesFlemischCase4LF(string& filenameFine, string& filenameCoarse,
                               TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    const int initc = 1, endc = 1000;
    for (int ivol = initc; ivol <= endc; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
        
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = initc; ivol <= endc; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-initc);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 36; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_1_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_31"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_5" ] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_19"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_29"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_7_11"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_8_34"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_17_22"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_25_36"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_26_27"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_29_32"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_30_33"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_30_36"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_31_33"] = -121321313;
    dim_name_and_physical_tagFine[1]["fracIntersection_31_36"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_34_35"] = EIntersection;

    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void fixPossibleMissingIntersections(TPZGeoMesh* gmesh){
    // it may happen that the supplied mesh from DFN has missing intersections where 2 fractures intersect.
    // We treat this here. But it would be ideal to make DFN robust enough so we could remove this function whatsoever
    // The main idea is to check if a fracture element has more than 1 fracture neightbor. If so, an intersection element is
    // needed to hybridzie the region.
    cout << "\n\t==================== Searching for problematic intersection elements ====================\n" << endl;
    
    int nInterCreated = 0;
    for (auto gel : gmesh->ElementVec()) {
        if (!gel || gel->MaterialId() != EFracture) {
            continue;
        }
        if (gel->Dimension() != 2) {
            DebugStop(); // Should it work with 2D problems? (1d fractures)
        }
        
        const int firstedge = gel->FirstSide(1);
        const int lastedge = gel->FirstSide(2);
        for (int iside = firstedge; iside < lastedge; iside++) {
            TPZStack<TPZGeoEl*> interEls;
            int nFracElsForSide = 0;
            TPZGeoElSide gside(gel,iside);
            TPZGeoElSide neig = gside.Neighbour();
            for (; gside != neig; neig++) {
                TPZGeoEl* neigel = neig.Element();
                if (neigel->MaterialId() == EIntersection) {
                    interEls.push_back(neigel);
                }
                if (neigel->MaterialId() == EFracture) {
                    nFracElsForSide++;
                }
            }
            if (interEls.size() > 1) {
                if (interEls.size() > 2) {
                    DebugStop(); // there are 3 intersection elements in the same place! Please check why...
                }
                cout << "Found two intersection elements at the same place!" << endl;
                cout << "Manually deleting intersection geoel..." << endl;
                TPZGeoEl* gelduplicate = interEls[1];
                const int64_t duplicateIndex = gelduplicate->Index();
                gelduplicate->RemoveConnectivities();
                delete gelduplicate;
                gmesh->ElementVec()[duplicateIndex] = nullptr;
            }
            if (nFracElsForSide > 1 && !interEls.size()) {
                cout << "nfracs for this side: " << nFracElsForSide << endl;
                cout << "Manually creating intersection geoel..." << endl;
                TPZGeoElBC(gel, iside, EIntersection);
                nInterCreated++;
            }
            
        }
    }
    if (nInterCreated) {
        cout << "\n\t\t==================== MESSAGE ====================" << endl;
        cout << "\t\t=================================================" << endl;
        cout << "\n\t- Imported mesh has " << nInterCreated << " missing intersection elements" << endl;
        cout << "\t- These were created here..\n" << endl;
    }
        
    gmesh->BuildConnectivity();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname) {
    std::set<int> matids;
    matids.insert(EVolume);
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences(); // compute integral in the multiphysics mesh
    TPZVec<STATE> vecint = cmesh->Integrate(varname, matids);
    if ((varname == "Pressure" && vecint.size() != 1) ||
        (varname == "Flux" && vecint.size() != 3)){
        DebugStop();
    }
    if (varname == "Pressure")
        return vecint[0];
    else if (varname == "Flux")
        return vecint[1];
    else
        DebugStop();
    
    return -100000; // default value so compiler does not complain
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
