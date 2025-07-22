// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>

//#define USEMAIN 

#ifndef USEMAIN
// Unit test includes
#include <catch2/catch.hpp>
#endif

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(const int& caseToSim);
TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data, int const simcase);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void modifyBCsForInletOutlet(TPZGeoMesh* gmesh);
void modifyBCsForInletOutletFracTransport(TPZGeoMesh* gmesh);
void MapFractureIntersection(const std::string &filenameBase, std::map<std::string,int> &matmap, int firstfracintersect, std::map<int,std::pair<int,int>> &matidtoFractures);
void CreateIntersectionElementForEachFrac(TPZGeoMesh* gmeshfine,
                                          std::map<int,std::pair<int,int>>& matidtoFractures,
                                          const int fracInitMatId, const int fracinc, const int FractureHybridPressureMatId);


void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse, int const  caseToSim);
const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);

enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFractureHybridPressure, EFracPressure,
    /*10*/EFracture, EFracNoFlux, EIntersection, EIntersectionEnd, EPLossAtIntersect, EFracture2, EFracNoFlux2, EIntersection2,
    /*HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
};

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

#ifdef USEMAIN
int main(){
    int CaseToSim = 7;
    RunProblem(CaseToSim);
    return 0;
}
#else

// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("Two_El_One_Frac_cte_pressure","[test_frac_4spaces_3D]"){
    RunProblem(0);
}
// ---- Test 1 ----
TEST_CASE("Two_El_One_Frac_lin_pressure","[test_frac_4spaces_3D]"){
    RunProblem(1);
}
// ---- Test 2 ----
TEST_CASE("Four_El_Two_Frac_cte_pressure","[test_frac_4spaces_3D]"){
    RunProblem(2);
}
// ---- Test 3 ----
TEST_CASE("Four_El_Two_Frac_lin_pressure","[test_frac_4spaces_3D]"){
    RunProblem(3);
}
// ---- Test 4 ----
TEST_CASE("Two_El_One_Frac_cte_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(4);
}
// ---- Test 5 ----
TEST_CASE("Two_El_One_Frac_lin_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(5);
}
// ---- Test 6 ----
TEST_CASE("Four_El_Two_Frac_cte_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(6);
}
// ---- Test 7 ----
TEST_CASE("Four_El_Two_Frac_lin_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(7);
}

#endif

// Case0: TwoElements, one fracture constant pressure
// Case1: TwoElements, one fracture linear pressure
// Case2: FourElements, two intersecting fractures constant pressure
// Case3: FourElements, two intersecting fractures linear pressure

// Case4: TwoElements, one fracture constant pressure with tranport
// Case5: TwoElements, one fracture linear pressure with tranport
// Case6: FourElements, two intersecting fractures constant pressure with tranport
// Case7: FourElements, two intersecting fractures linear pressure with tranport

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------

void RunProblem( const int &caseToSim){
    
    cout << "\n\n\t\t---------------------- Start of Simulation " << caseToSim <<  " ------------------------\n" << endl;
    
    string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
    string logpath = basemeshpath + "/../DFNIMRS/log4cxx.cfg";
    TPZLogger::InitializePZLOG(logpath);
//    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    string filenameCoarse,filenameFine;
    if (caseToSim==0 || caseToSim==1 || caseToSim==4 || caseToSim==5) {
        //OneFractures
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/TwoElsUnitTestCoarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/TwoElsUnitTestFine.msh";
    }else if (caseToSim==2 || caseToSim==3 || caseToSim==6 || caseToSim==7){
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/intersectCoarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/intersectFine.msh";
    }
    else{
        DebugStop();
    }
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    ReadMeshes(filenameFine,filenameCoarse,gmeshfine,gmeshcoarse,caseToSim);
    TMRSDataTransfer sim_data;
    FillDataTransfer(sim_data, caseToSim);
    
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
    
    // ----- Approximation space -----
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    
    // ----- Setting gmesh -----
    // Code takes a fine and a coarse mesh to generate MHM data structure
    TMRSApproxSpaceGenerator aspace;
    aspace.InitMatIdForMergeMeshes() = EInitVolumeMatForMHM;
//    sim_data.mTFracProperties.m_matid = EFracture;
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
    aspace.SetGeometry(gmeshfine,gmeshcoarse);
//  aspace.SetGeometry(gmeshfine);
    aspace.PrintGeometry("geomesh_aftermerge.vtk");
    
    // ----- Setting the global data transfer -----
    aspace.SetDataTransfer(sim_data);
    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = false;
    int n_threads = 0;
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    REAL mass=0.0;
    if(caseToSim == 4 || caseToSim == 5 || caseToSim == 6 || caseToSim == 7){
        modifyBCsForInletOutletFracTransport(aspace.GetGeometry());
        aspace.BuildAuxTransportCmesh();
        TPZCompMesh * transport_operator = aspace.GetTransportOperator();
        
        std::ofstream name("TransportOperator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(transport_operator, name);
        
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
        REAL InntMassFrac=0.0;
        for (int it = 1; it <= n_steps; it++) {
            sim_time = it*dt;
            sfi_analysis->m_transport_module->SetCurrentTime(dt);
            const int typeToPostProc = 2; // only transport
            sfi_analysis->PostProcessTimeStep(typeToPostProc);
            sfi_analysis->RunTimeStep(); // runs mixed and transport problems
            mixed_operator->LoadSolution(mixed_operator->Solution());
            if (sim_time >=  current_report_time) {
                std::cout << "Time step number:  " << it << std::endl;
                std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
                mixed_operator->UpdatePreviousState(-1.);
                const int typeToPostProcEnd = 0; // p/flux and transport
                sfi_analysis->PostProcessTimeStep(typeToPostProcEnd);
                pos++;
                current_report_time =reporting_times[pos];
                
                // computes integral of saturation for elements of certain matid
//                InntMassFrac=sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(aspace.FractureMatId());
                fileCilamce<<current_report_time/(86400*365)<<", "<<InntMassFrac<<std::endl;
               
                mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
                std::cout << "Mass report at time : " << sim_time << std::endl;
                std::cout << "Mass integral :  " << mass << std::endl;
            }
        }
    }
    else{
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
    const STATE integratedpressure = ComputeIntegralOverDomain(mixed_operator,pvarname);
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;

    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(mixed_operator,qvarname);
    if (fabs(integratedflux) < 1.e-14 ) integratedflux = 0.; // to make Approx(0.) work
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;

    if(caseToSim > 3){
        std::cout << "\nMass  = " << mass << std::endl;
    }

#ifndef USEMAIN
    // ----- Comparing with analytical solution -----
    REQUIRE( integratedpressure == Approx( 8.0 ) ); // Approx is from catch2 lib
    if (caseToSim == 1 || caseToSim == 3 || caseToSim == 5 || caseToSim == 7) // linear pressure variation
        // DeltaP = 2, DeltaS = 2. + fracwidth = 3 | q = - DeltaP/DeltaS = 2/3. Integrated = 2./3. * 8/ = 16./3.
        REQUIRE( integratedflux == Approx( 16./3. ) ); // Approx is from catch2 lib
    if (caseToSim == 0 || caseToSim == 2 || caseToSim == 4 || caseToSim == 6)  // cte pressure
        REQUIRE( integratedflux == Approx( 0.) ); // Approx is from catch2 lib
    if(caseToSim==5)
        REQUIRE( mass == Approx(2.38975).epsilon(1.e-3) ); // Approx is from catch2 lib
    if(caseToSim==7)
        REQUIRE( mass == Approx(3.58463).epsilon(1.e-3) ); // Approx is from catch2 lib
    if(caseToSim == 4 || caseToSim == 6)
        REQUIRE( mass == Approx( 0.) ); // Approx is from catch2 lib
#endif
    
    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
    
    cout << "\n\n\t\t****************** End of Simulation " << caseToSim <<  " ******************\n" << endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data, int const simcase){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, unit_pressure = 1.0, zero_pressure = 0., pressure_two = 2.;
    
    // Domain material
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fracture"] = EFracture;
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fracture2"] = EFracture2;
//    sim_data.mTGeometry.mDomainFracIntersectionNameAndMatId["EIntersection"] = EIntersection;
    sim_data.mTGeometry.m_pressureMatId = EFractureHybridPressure;
    
    // Domain boundary conditions
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,pressure_two);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,zero_pressure);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,unit_pressure);
	
    // Fracture boundary conditions
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracNoFlux] = std::make_pair(N_Type, zero_flux);
    sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracNoFlux2] = std::make_pair(N_Type, zero_flux);    
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracInlet] = std::make_pair(D_Type, pressure_two);
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracOutlet] = std::make_pair(D_Type, zero_pressure);
	
    if(simcase == 2 || simcase == 6){
        sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EIntersection] = std::make_pair(D_Type, zero_pressure);
        sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EIntersection2] = std::make_pair(D_Type, zero_pressure);
    }
    else if(simcase == 3 || simcase == 7){
        sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EIntersection] = std::make_pair(Mixed_Type, 2.);
        sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EIntersection2] = std::make_pair(Mixed_Type, 2.);
    }

    
//    sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EPLossAtIntersect] = std::make_pair(Mixed_Type, 0.5);
    
	
	
    sim_data.mTFracIntersectProperties.m_IntersectionPressureLossId = EPLossAtIntersect;
    
    if (simcase==0 || simcase==2 || simcase==4 || simcase==6 ) {
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,pressure_two);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracOutlet] = std::make_pair(D_Type, pressure_two);
    }
    

    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
//    sim_data.mTFracIntersectProperties.m_IntersectionId = EPLossAtIntersect;
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EOutlet] = std::make_pair(N_Type, 1.);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EInlet] = std::make_pair(D_Type, 1.);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[ENoflux] = std::make_pair(5, 1.);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EFaceBCPressure] = std::make_pair(5, 1.);
//	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EFracNoFlux] = std::make_pair(5, 1.);
//    sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EFracNoFlux2] = std::make_pair(5, 1.);
//	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EFracInlet] = std::make_pair(5, 1.);
//	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[EFracOutlet] = std::make_pair(5, 1.);
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
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
    sim_data.mTNumerics.m_n_steps = 5;
    sim_data.mTNumerics.m_dt      = 1.;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
    // Fracture permeability
    const REAL fracwidth = 1., fracperm = 1.;
    if(simcase == 0 || simcase == 1 || simcase == 4 || simcase == 5){
        TMRSDataTransfer::TFracProperties::FracProp fracprop;
        fracprop.m_perm = fracperm;
        fracprop.m_width = fracwidth;
        fracprop.m_fracbc.insert(EFracNoFlux);
        fracprop.m_fracIntersectMatID = EIntersection;
        sim_data.mTFracProperties.m_fracprops[EFracture] = fracprop;
    }else{
        {
            TMRSDataTransfer::TFracProperties::FracProp fracprop;
            fracprop.m_perm = fracperm;
            fracprop.m_width = fracwidth;
            fracprop.m_fracbc.insert(EFracNoFlux);
            fracprop.m_fracbc.insert(EFracInlet);
            fracprop.m_fracbc.insert(EFracOutlet);
            fracprop.m_fracIntersectMatID = EIntersection;
            sim_data.mTFracProperties.m_fracprops[EFracture] = fracprop;
        }
        {
            TMRSDataTransfer::TFracProperties::FracProp fracprop;
            fracprop.m_perm = fracperm;
            fracprop.m_width = fracwidth;
            fracprop.m_fracbc.insert(EFracNoFlux2);
            fracprop.m_fracIntersectMatID = EIntersection2;
            sim_data.mTFracProperties.m_fracprops[EFracture2] = fracprop;
        }
    }
    
    //FracAndReservoirProperties
    REAL kappa = 1.0;
    int  id1 = EVolume;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa;
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    REAL resPorosity = 0.2;
    REAL fracPorosity = 0.2;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(4);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[0] = std::make_tuple(EVolume, resPorosity,1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[1] = std::make_tuple(EFracture, fracPorosity, fracwidth);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[2] = std::make_tuple(EFractureHybridPressure, fracPorosity, fracwidth*fracwidth);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[3] = std::make_tuple(EFracture2, fracPorosity, fracwidth);

    
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

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int const caseToSim){
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    for (int ivol = 1; ivol <= 4; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 1; ivol <= 4; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = EFracture;
    
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture2;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EFracture;

    // Fractures BCs
    if(caseToSim == 0 || caseToSim == 1 || caseToSim == 4 || caseToSim == 5){
        dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux;
    }
    else{
        dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux2;
        dim_name_and_physical_tagFine[1]["BCfrac1"] = EFracNoFlux;
    }
    
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection2;
    
    std::map<int,std::pair<int,int>> matidtoFractures;
    if(caseToSim == 2 || caseToSim == 3 || caseToSim == 6 || caseToSim == 7){
        int firstfracintersect = EInitVolumeMatForMHM + 4;
        firstfracintersect = (1+firstfracintersect/100)*100;
        MapFractureIntersection(filenameFine, dim_name_and_physical_tagFine[1], firstfracintersect, matidtoFractures);
    }
            
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);

    
    if(caseToSim == 1 || caseToSim == 3 || caseToSim == 5 || caseToSim == 7){
        modifyBCsForInletOutlet(gmeshfine);
    }
    
    if(caseToSim == 2 || caseToSim == 3 || caseToSim == 6 || caseToSim == 7){
        const int fracinc = 5;
        CreateIntersectionElementForEachFrac(gmeshfine,matidtoFractures,EFracture,fracinc,EFractureHybridPressure);
    }
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateIntersectionElementForEachFrac(TPZGeoMesh* gmeshfine,
                                          std::map<int,std::pair<int,int>>& matidtoFractures,
                                          const int fracInitMatId, const int fracinc, const int FractureHybridPressureMatId) {
    int64_t nelem = gmeshfine->NElements();
    for(int64_t el = 0; el<nelem; el++)
    {
        TPZGeoEl *gel = gmeshfine->Element(el);
        if(!gel) continue;
        int matid = gel->MaterialId();
        auto it = matidtoFractures.find(matid);
        if(it != matidtoFractures.end())
        {
            if(gel->Dimension() == 3) DebugStop();
            TPZGeoElSide gelside(gel);
            int frac1 = it->second.first;
            int matid1 = EIntersection;
//            int frac1bcid = fracInitMatId + frac1*fracinc + 1;
            std::set<int> frac1bcids = {EFracNoFlux,EFracInlet,EFracOutlet};
            auto frac1bc = gelside.HasNeighbour(frac1bcids);
            if(frac1bc)
            {
                frac1bc.Element()->SetMaterialId(matid1);
            }
            else
            {
                TPZGeoElBC(gelside,matid1);
            }
            int frac2 = it->second.second;
            int matid2 = EIntersection2;
//            int frac2bcid = fracInitMatId + frac2*fracinc + 1;
            auto frac2bc = gelside.HasNeighbour(EFracNoFlux2);
            if(frac2bc)
            {
                frac2bc.Element()->SetMaterialId(matid2);
            }
            else
            {
                TPZGeoElBC(gelside,matid2);
            }
            // create the geometric element for receiving the hybridized pressure element
            auto fracintersect = gelside.HasNeighbour(FractureHybridPressureMatId);
            if(!fracintersect)
            {
                TPZGeoElBC(gelside,FractureHybridPressureMatId);
            }
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void MapFractureIntersection(const std::string &filenameBase, std::map<std::string,int> &matmap, int firstfracintersect, std::map<int,std::pair<int,int>> &matidtoFractures)
{
    // Creating gmsh reader
    TPZGmshReader  Geometry;
    std::string filename = filenameBase;
    std::ifstream input(filename);
    Geometry.ReadPhysicalProperties4(input); // There are some automatic associations that are not used. So dont worry if the code alerts about it
    auto physicaltags = Geometry.GetDimPhysicalTagName();
    // loop over all physical tags of dimension 1
    const std::string begins("fracIntersection");
    const int beglength = begins.size();
    for(auto iter : physicaltags[1])
    {
        auto name = iter.second;
        if(name.compare(0,beglength,begins) == 0)
        {
            auto pos1 = name.find_first_of("_",0);
            auto pos2 = name.find_last_of("_");
            std::string strfac1 = name.substr(pos1+1,pos2-pos1-1);
            std::string strfac2 = name.substr(pos2+1);
            int frac1 = std::stoi(strfac1);
            int frac2 = std::stoi(strfac2);
            matmap[name] = firstfracintersect;
            matidtoFractures[firstfracintersect] = std::pair<int,int>(frac1,frac2);
            firstfracintersect++;
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void modifyBCsForInletOutlet(TPZGeoMesh* gmesh) {
    
    const REAL cmin = -1, cmax = 1.;
    const int pos = 1;
    for (auto gel : gmesh->ElementVec()) {
        const int gelmatid = gel->MaterialId();
        if (gelmatid != EFaceBCPressure && gelmatid != EFracNoFlux)
            continue;
        
        TPZManVector<REAL,3> centqsi(2,0.), cent(3,0.);
        if (gelmatid == EFracNoFlux)
            centqsi.Resize(1, 0.);
        
        gel->CenterPoint(gel->NSides()-1, centqsi);
        gel->X(centqsi, cent);
        const REAL minsub = fabs(cent[pos] - cmin);
        const REAL maxsub = fabs(cent[pos] - cmax);
        if (minsub < ZeroTolerance()) {
            if(gelmatid == EFracNoFlux){
                gel->SetMaterialId(EFracInlet);
            }
            else{
                gel->SetMaterialId(EInlet);
            }
            
        }
        else if (maxsub < ZeroTolerance()) {
            if(gelmatid == EFracNoFlux){
                gel->SetMaterialId(EFracOutlet);
            }
            else {
                gel->SetMaterialId(EOutlet);
            }
            
        }
        else {
            if(gelmatid == EFracNoFlux){
                gel->SetMaterialId(EFracNoFlux);
            }
            else{
                gel->SetMaterialId(ENoflux);
            }
            
        }
    }
}
void modifyBCsForInletOutletFracTransport(TPZGeoMesh* gmesh) {
    
    const REAL cmin = -1, cmax = 1.;
    const int pos = 1;
    for (auto gel : gmesh->ElementVec()) {
        const int gelmatid = gel->MaterialId();
        if(gelmatid == EFracInlet){
            gel->SetMaterialId(EInlet);
        }
        if(gelmatid == EFracOutlet){
            gel->SetMaterialId(EOutlet);
        }
    }
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
        return vecint[1];
    else
        DebugStop();
    
    return -100000; // default value so compiler does not complain
}
