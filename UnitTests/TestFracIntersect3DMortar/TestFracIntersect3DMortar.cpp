#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

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
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
TPZGeoMesh *ReadFractureMeshCase0();
TPZGeoMesh *ReadFractureMeshCase1();
const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);

enum EMatid {ENone, EVolume, EInlet, EOutlet, ENoflux, EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure, EDONTUSEGLOBALFRACID, EIntersection, EIntersectionEnd, EPLossAtIntersect};
int globFracID = 10;
// ----- End of Functions -----

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("constant_pressure","[test_intersection_3D_Mortar]"){
    RunTest(0);
}
// ---- Test 0 ----
TEST_CASE("linear_pressure","[test_intersection_3D_Mortar]"){
    RunTest(1);
}

// ----- Main Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("testintersect3d");
#endif
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
void RunTest(const int caseToSim)
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for TestIntersect3D target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
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
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
    
    // ----- Setting gmesh -----
    aspace.SetGeometry(gmesh);
    
    // ----- Setting the global data transfer -----
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
	//  sim_data.mTFracProperties.m_matid = globFracID; // This is now set for each fracture
    aspace.SetDataTransfer(sim_data);

    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    mixAnalisys->Assemble();
//    const char* name = "RHS";
//    TPZMatrixSolver<STATE>* matsol = dynamic_cast<TPZMatrixSolver<STATE>*>(mixAnalisys->Solver());
//    std::ofstream oo("matbad.txt");
//    matsol->Matrix()->Print(oo);
//    mixAnalisys->Rhs().Print(name);
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // ----- Post processing -----
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    const int dimToPost = 3;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    
    // ----- Compute integral of pressure and flux over domain and compare with analytical solution -----
    const std::string pvarname = "Pressure";
    const STATE integratedpressure = ComputeIntegralOverDomain(mixed_operator,pvarname);
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(mixed_operator,qvarname);
    if (fabs(integratedflux) < 1.e-14 ) integratedflux = 0.; // to make Approx(0.) work
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;
    
    // ----- Comparing with analytical solution -----
    // Domain volume is 2*2*2=8. If p cte: 1*8 = 8. If p varies linearly from 2 to 0: ((2-0)/2) * 8 = 8
    REQUIRE( integratedpressure == Approx( 8.0 ) ); // Approx is from catch2 lib
    if (caseToSim == 1) // linear pressure variation
        REQUIRE( integratedflux == Approx( 8./3. ) ); // Approx is from catch2 lib
    else // cte ppressyre
        REQUIRE( integratedflux == Approx( 0.) ); // Approx is from catch2 lib
    
    // ----- Cleaning up -----
    delete gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim) {
    string basemeshpath(FRACMESHES);
    std::string filename = basemeshpath + "/2DMeshes/1fracNoBnd.msh";
    switch (caseToSim) {
        case 0: {
            gmesh = ReadFractureMeshCase0();
        }
            break;
        case 1: {
            gmesh = ReadFractureMeshCase1();
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
    REAL pressure_in = 1.0 ;
    REAL zero_flux = 0.0;
    
    
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fractures"] = globFracID;

    // NS: What are these?
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    

    
    // Boundary conditions
    if (caseToSim == 0){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,1.);
        
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracPressure] = std::make_pair(D_Type, 1.);
    }
    else if (caseToSim == 1){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
		
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracInlet] = std::make_pair(D_Type, 2.);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracOutlet] = std::make_pair(D_Type, 0.);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracNoFlux] = std::make_pair(N_Type, zero_flux);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EPLossAtIntersect] = std::make_pair(Mixed_Type, 0.5);
        
        sim_data.mTFracIntersectProperties.m_IntersectionPressureLossId = EPLossAtIntersect;
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
    
    // Other properties?
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0; //*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    
    //FracAndReservoirProperties
//    sim_data.mTFracProperties.m_Permeability = 1.;
	TMRSDataTransfer::TFracProperties::FracProp fracprop;
	fracprop.m_perm = 1.;
	fracprop.m_width = 1.;
	fracprop.m_fracbc.insert(100000); // Does not matter for mortar spaces
	fracprop.m_fracIntersectMatID = EIntersection;
	sim_data.mTFracProperties.m_fracprops[globFracID] = fracprop;
	
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
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    return sim_data;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase0(){
    
    TPZGeoMesh* gmesh = nullptr;
    
    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,2,2};
    const MMeshType elType = MMeshType::EHexahedral;

    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);

    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {6,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {7,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {9,11};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {10,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {8,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);

    nodesId = {9,15};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {14,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {8,2};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {2,3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {3,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);

    nodesId = {8,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EIntersection,*gmesh,index);
    
    TPZManVector<int64_t,4> nodesIdVec = {6,7,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {2,3,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
    gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
    
    gmesh->BuildConnectivity();

    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase1(){
    
    TPZGeoMesh* gmesh = nullptr;
    
    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,2,2};
    const MMeshType elType = MMeshType::EHexahedral;

    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);

    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {6,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracInlet,*gmesh,index);
    nodesId = {7,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {9,11};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracOutlet,*gmesh,index);
    nodesId = {10,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {8,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);

    nodesId = {9,15};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {14,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {8,2};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {2,3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {3,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);

    nodesId = {8,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EIntersection,*gmesh,index);
    
    TPZManVector<int64_t,4> nodesIdVec = {6,7,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {2,3,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
    gmesh = gen3d.BuildBoundaryElements(ENoflux, ENoflux, EInlet, ENoflux, EOutlet, ENoflux);
    
    gmesh->BuildConnectivity();

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
