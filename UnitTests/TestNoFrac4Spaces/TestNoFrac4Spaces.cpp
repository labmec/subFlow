// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>

// Unit test includes
#include <catch2/catch.hpp>

// ----- Functions -----
void RunTest(const int caseToSim);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TPZGeoMesh*& gmeshcoarse,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);
void ModifyBcsForTransport(TPZGeoMesh* gmesh, int inletBc, int outletBc);

TPZGeoMesh *ReadFractureMeshCase0();
void ReadFractureMeshCase1(TPZGeoMesh*& gmesh,TPZGeoMesh*& gmeshcoarse);


enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EDONTUSEGLOBALFRACID, EIntersection, EIntersectionEnd, EPLossAtIntersect,
    /*14 HAS TO BE LAST*/EInitVolumeMatForMHM
};

auto exactSol = [](const TPZVec<REAL>& loc,
                   TPZVec<STATE>& u,
                   TPZFMatrix<STATE>& gradU){
    const auto& x = loc[0];
    const auto& y = loc[1];
    const auto& z = loc[2];
    u[0] = 1. - x;
    gradU(0,0) = 0.; // not used
};

// ----- End of Functions -----

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Log -----
#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------

// ----- Test cases -----
// ---- Test 1 ----
// 1: 2x1x1 domain with a coarse and a fine mesh to generate MHM structure. Cte pressure
TEST_CASE("constant_pressure","[test_nofrac_3D]"){
    RunTest(1);
}
// ---- Test 2 ----
// 2x1x1 domain with a coarse and a fine mesh to generate MHM structure. Linear pressure
TEST_CASE("linear_pressure","[test_nofrac_3D]"){
    RunTest(2);
}

// ---- Test 3 ----
// 2x1x1 domain with a coarse and a fine mesh to generate MHM structure. Linear pressure and transport
TEST_CASE("transport_linear_pressure","[test_nofrac_3D]"){
    RunTest(3);
}

//int main(){
//    RunTest(3);
//    return 0;
//}

void RunTest(const int caseToSim)
{
    cout << "\n\n\t\t---------------------- Start of Simulation " << caseToSim <<  " ------------------------\n" << endl;
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    TMRSDataTransfer sim_data;
    CreateGMeshAndDataTransfer(gmeshfine,gmeshcoarse,sim_data,caseToSim);
    
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
    
    // ----- Approximation space -----
    TMRSApproxSpaceGenerator aspace;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    
    // ----- Setting gmesh -----
//    aspace.SetGeometry(gmeshfine,gmeshcoarse);
    aspace.SetGeometry(gmeshfine); // Not testing MHM yet
    
    // ----- Setting the global data transfer -----
    aspace.SetDataTransfer(sim_data);
    if(caseToSim == 2 || caseToSim == 3){
        aspace.SetForcingFunctionBC(exactSol);
    }
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    REAL mass=0.0;
    if(caseToSim == 3){
        ModifyBcsForTransport(aspace.mGeometry, EInlet,  EOutlet);
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
        TPZFastCondensedElement::fSkipLoadSolution = false;
        bool first=true;
        const int typeOfPostProc = 2; // 0: both, 1: p/sigma, 2: transport
        for (int it = 1; it <= n_steps; it++) {
            sim_time = it*dt;
            sfi_analysis->m_transport_module->SetCurrentTime(dt);
            sfi_analysis->PostProcessTimeStep(typeOfPostProc);
            sfi_analysis->RunTimeStep();
            mixed_operator->LoadSolution(mixed_operator->Solution());
            if (sim_time >=  current_report_time) {
                std::cout << "Time step number:  " << it << std::endl;
                std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
                mixed_operator->UpdatePreviousState(-1.);
                const int typeOfPostProcEnd = 0; // both p/flux and transport
                sfi_analysis->PostProcessTimeStep(typeOfPostProcEnd);
                pos++;
                current_report_time = reporting_times[pos];
                REAL InntMassFrac = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(10);
               
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
        mixAnalisys->fsoltransfer.TransferFromMultiphysics();
        const int dimToPost = 3;
        mixAnalisys->PostProcessTimeStep(dimToPost);
    }
    
    // ----- Compute integral of pressure and flux over domain and compare with analytical solution -----
    const std::string pvarname = "Pressure";
    STATE integratedpressure = ComputeIntegralOverDomain(mixed_operator,pvarname);
    if (fabs(integratedpressure) < 1.e-14 ) integratedpressure = 0.; // to make Approx(0.) work
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(mixed_operator,qvarname);
    if (fabs(integratedflux) < 1.e-14 ) integratedflux = 0.; // to make Approx(0.) work
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;
    
    // ----- Comparing with analytical solution -----
    // Results are intuitive by looking at paraview plots of the pressure and flux
    if (caseToSim == 2){ // linear pressure variation
        REQUIRE( integratedpressure == Approx( 0. ) ); // Approx is from catch2 lib
        REQUIRE( integratedflux == Approx( 8.0 ) ); // Approx is from catch2 lib
    }
    else if (caseToSim == 1){
        REQUIRE( integratedpressure == Approx( 8.0 ) ); // Approx is from catch2 lib
        REQUIRE( integratedflux == Approx( 0.) ); // Approx is from catch2 lib
    }
    else if(caseToSim == 3){
        REQUIRE( mass == Approx( 0.8) ); // Approx is from catch2 lib
    }
    else{
        DebugStop();
    }
    
    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
    
    cout << "\n\n\t\t****************** End of Simulation " << caseToSim <<  " ******************\n" << endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TPZGeoMesh*& gmeshcoarse, TMRSDataTransfer &sim_data, const int caseToSim) {
    switch (caseToSim) {
        case 0: {
            gmesh = ReadFractureMeshCase0();
        }
            break;
        case 1: {
            ReadFractureMeshCase1(gmesh,gmeshcoarse);
        }
        case 2: {
            ReadFractureMeshCase1(gmesh, gmeshcoarse);
        }
            break;
        case 3: {
            ReadFractureMeshCase1(gmesh, gmeshcoarse);
        }
            break;
        default:
            DebugStop();
    }
            
    sim_data  = SettingFracturesSimple(caseToSim);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer SettingFracturesSimple(const int caseToSim){
    
    // Fracture material
    TMRSDataTransfer sim_data;
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0;
    
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    
    // Boundary conditions Mixed Problem
    if (caseToSim < 4) {
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,1.);

    }
    else {
        DebugStop();
    }
    
    // Boundary conditions Transport Problemm
    if (caseToSim < 4) {
		sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EInlet] = std::make_pair(D_Type, 1.);
		sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EOutlet] = std::make_pair(N_Type, 1.);
    }
    else {
        DebugStop();
    }
    
    // Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.1;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
            
    // Other properties?
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_n_steps = 10;
    sim_data.mTNumerics.m_dt      = 0.5;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
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

TPZGeoMesh *ReadFractureMeshCase0(){
            
    TPZGeoMesh* gmesh = nullptr;

    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,1,2};
    const MMeshType elType = MMeshType::EHexahedral;
    
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);
    gmesh = gen3d.BuildBoundaryElements(EInlet,EFaceBCPressure,EFaceBCPressure,EOutlet,EFaceBCPressure,EFaceBCPressure);
    
    gmesh->BuildConnectivity();
        
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void ReadFractureMeshCase1(TPZGeoMesh*& gmeshfine,TPZGeoMesh*& gmeshcoarse){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    string basemeshpath(FRACMESHES);
    
    // ===================> Coarse mesh
    
    // Gmsh file
    std::string filenamecoarse = basemeshpath + "/verificationMHMNoHybrid/twoElCoarse.msh";
    
    // Domain
    std::string volbase = "c";
    for (int ivol = 1; ivol <= 2; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
    gmeshcoarse = generateGMeshWithPhysTagVec(filenamecoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
    
    // Gmsh file
    std::string filenamefine = basemeshpath + "/verificationMHMNoHybrid/twoElFine.msh";
        
    // Domain
    for (int ivol = 1; ivol <= 2; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
//        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc2"] = EInlet;
            
    gmeshfine = generateGMeshWithPhysTagVec(filenamefine,dim_name_and_physical_tagFine);
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
        return vecint[0];
    else
        DebugStop();
    
    return -100000; // default value so compiler does not complain
}
void ModifyBcsForTransport(TPZGeoMesh* gmesh, int inletBc, int outletBc) {
    int nels = gmesh->NElements();
    int meshdim = gmesh->Dimension();
    for (int iel=0; iel<nels; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() == meshdim) {
            continue;
        }
        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        if(xcenter[0]!=0){
            gel->SetMaterialId(outletBc);
        }
        else if (xcenter[0]==0) {
            gel->SetMaterialId(inletBc);
        }
        else{
            DebugStop();
        }
    }
}
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
