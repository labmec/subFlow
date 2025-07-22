#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// C++ includes

// PZ includes
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"

// Unit test includes
#include <catch2/catch.hpp>

// ----- Functions -----

void CaseOnlyFractures(const int caseToSim);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);

// ---- Driver Function ----
void Test2frac(const int& caseToSim);

int globFracID = 10;
enum EMatid {ENone, EDomain, EInlet, EOutlet, ENoflux, EPressure, EIntersection, EIntersectionEnd};

// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("2frac","[onlyfracintersect_test]"){
    Test2frac(0);
}

// ----- Namespaces -----
using namespace std;

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("onlyfractures");
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
void Test2frac(const int& caseToSim)
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for OnlyFractures target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    // -------------- Reading ONLY the fractures of mesh from DFN --------------
    TPZGeoMesh *gmesh = nullptr;
    TMRSDataTransfer sim_data;
    CreateGMeshAndDataTransfer(gmesh,sim_data,caseToSim);
    
    const bool printgmesh = true;
    if (printgmesh) {
        std::ofstream name("GeoMeshFractures.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name);
    }
        
    // -------------- Approximation space --------------
    TMRSApproxSpaceGenerator aspace;
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E2Space;
    
    // -------------- Setting gmesh --------------
    aspace.SetGeometry(gmesh);
    
    // -------------- Setting the global data transfer --------------
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
//    sim_data.mTFracProperties.m_matid = globFracID; // This is now set for each fracture
	
    aspace.SetDataTransfer(sim_data);
    
    // -------------- Creates de multiphysics compmesh --------------
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
        
    // -------------- Analysis parameters --------------
    bool must_opt_band_width_Q = false;
    int n_threads = 0;
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    
    // -------------- Setting analysis --------------
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    // -------------- End of running problem --------------
    
    // -------------- Post processing --------------
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    mixAnalisys->PostProcessTimeStep(2);
    
    // -------------- Integral of pressure over domain --------------
    const std::string varname = "Pressure";
    std::set<int> matids;
    matids.insert(globFracID);
    mixed_operator->Reference()->ResetReference();
    mixed_operator->LoadReferences();
    TPZVec<STATE> vecint = mixed_operator->Integrate(varname, matids);
    if (vecint.size() != 1){
        DebugStop();
    }
    const STATE integratedpressure = vecint[0];
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    // Analytic solution is unit constant pressure (ctepressuresol = 1)
    // There are two 2x2 fracture. Therefore total area is eight (area = 8)
    // Thus, the integral of pressure should be area*ctepressuresol = 8
    REQUIRE( integratedpressure == Approx( 8.0 ) ); // Approx is from catch2 lib
    
    // Cleaning up
    delete gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim) {
    string basemeshpath(FRACMESHES);
    switch (caseToSim) {
        case 0: {
            std::string filename = basemeshpath + "/2DMeshes/2fracFromDfn.msh";
            gmesh = ReadFractureMeshCase0(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        default:
            DebugStop();
    }
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer SettingFracturesSimple(const int caseToSim){
    
    // Fracture material
    TMRSDataTransfer sim_data;
    sim_data.mTGeometry.mDomainNameAndMatId["Fractures"] = globFracID;

    // NS: What are these?
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
    int D_Type = 0;
    int N_Type = 1;
    REAL pressure_in = 1.0 ;
    
    // Boundary conditions
    if (caseToSim == 0) {
        int bcfracid = EPressure;
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[bcfracid] = std::make_pair(D_Type,pressure_in);
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
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    
    //FracAndReservoirProperties
//    sim_data.mTFracProperties.m_Permeability = 0.00001;
	TMRSDataTransfer::TFracProperties::FracProp fracprop;
	fracprop.m_perm = 0.00001;
	fracprop.m_width = 1.;
	fracprop.m_fracbc.insert(EPressure);
	fracprop.m_fracIntersectMatID = EIntersection;
	sim_data.mTFracProperties.m_fracprops[globFracID] = fracprop;
	
    REAL kappa=1.0;
    int  id1=1;
    int  id2=2;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa;
    idPerm[id2]= kappa;
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

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
        
    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
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

