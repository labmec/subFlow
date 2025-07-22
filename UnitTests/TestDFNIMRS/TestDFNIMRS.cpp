// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include "TPZSimpleTimer.h"
#include <tpzgeoelrefpattern.h>
#include "json.hpp"
#include <filesystem>
#include <libInterpolate/Interpolate.hpp>
#include <libInterpolate/AnyInterpolator.hpp>
// include dfn filereader
#include "filereader.h"

#include "DFNMesh.h"

//#define USEMAIN

#ifndef USEMAIN
// Unit test includes
#include <catch2/catch.hpp>
#endif


// ----- Namespaces -----
using namespace std;
namespace fs = std::filesystem;

// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(string& filenameBase, const bool isLinPressure);
void ReadMeshesDFN(string& filenameBase, TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int& initVolForMergeMeshes, bool& isMHM, bool& needsMergeMeshes);
void CreateIntersectionElementForEachFrac(TPZGeoMesh* gmeshfine,
										  std::map<int,std::pair<int,int>>& matidtoFractures,
										  const int fracInitMatId, const int fracinc, const int FractureHybridPressureMatId);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void fixPossibleMissingIntersections(TMRSDataTransfer& sim_data, TPZGeoMesh* gmesh);

void FillDataTransferDFN(string& filenameBase, string& outputFolder, TMRSDataTransfer& sim_data);
void FilterZeroNeumann(std::string& outputFolder, TMRSDataTransfer& sim_data, TPZAutoPointer<TPZStructMatrix> strmat, TPZCompMesh* cmesh);

void FillPCteSol(TPZMultiphysicsCompMesh* mpcmesh, const REAL pcte);

bool fileExists(const fs::path& p, fs::file_status s = fs::file_status{});
void CreateOutputFolders(std::string& outputFolder);
void CopyInputFilesToOutputFolderAndFixFilename(std::string& filenameBase, std::string& outputFolder);

const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);

std::map<int,TPZVec<STATE>> computeIntegralOfNormalFlux(const std::set<int> &bcMatId, TPZMultiphysicsCompMesh *cmesh);

// TODO: Delete this enum? (May 2022)
enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EFracture, EIntersection, EIntersectionEnd, EPLossAtIntersect, EVolume2,
    /*15 HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
};

// ----- Logger -----
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
#ifdef USEMAIN
int main(int argc, char* argv[]){
    string basemeshpath(FRACMESHES);
    string filenameBase = basemeshpath + "/dfnimrs/intersectSnap";
    const bool isLinPressure = false;
    RunProblem(filenameBase,isLinPressure);
    return 0;
}
#else

// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("Four_El_Two_Frac_Snap_cte_pressure","[test_dfnimrs]"){
    string basemeshpath(FRACMESHES);
    string filenameBase = basemeshpath + "/dfnimrs/intersectSnap";
    const bool isLinPressure = false;
    RunProblem(filenameBase,isLinPressure);
}

// ---- Test 1 ----
TEST_CASE("Four_El_Two_Frac_cte_pressure","[test_dfnimrs]"){
    string basemeshpath(FRACMESHES);
    string filenameBase = basemeshpath + "/dfnimrs/twoElCoarse";
    const bool isLinPressure = false;
    RunProblem(filenameBase,isLinPressure);
}

// ---- Test 2 ----
TEST_CASE("Six_El_One_Frac_lin_pressure","[test_dfnimrs]"){
    string basemeshpath(FRACMESHES);
    string filenameBase = basemeshpath + "/dfnimrs/boxPerpFlux6el";
    const bool isLinPressure = true;
    RunProblem(filenameBase,isLinPressure);
}

#endif

// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void RunProblem(string& filenameBase, const bool isLinPressure)
{
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
    
    auto start_time = std::chrono::steady_clock::now();
    
	// ----- Simulation and printing parameters -----
    const bool isRefineMesh = false;
    const bool isPostProc = true;
	const bool isRunWithTranport = false;
	bool isMHM = true;
    bool needsMergeMeshes = true;
    const bool isFilterZeroNeumann = GENERATE(true,false);
	const int n_threads = 8;
    
    // ----- output folder stuff -----
    if(filenameBase.back() != '/') filenameBase = filenameBase + "/";
    std::string outputFolder = "Output/"+filenameBase.substr(filenameBase.find("dfnimrs/") + 8);
    CreateOutputFolders(outputFolder);
    CopyInputFilesToOutputFolderAndFixFilename(filenameBase,outputFolder);
    outputFolder = outputFolder.substr(0,outputFolder.find_last_of("/"));
    outputFolder = outputFolder + "/";

    
    // ----- Creating gmesh and data transfer -----
	int initVolForMergeMeshes = -1000000;
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
	ReadMeshesDFN(filenameBase, gmeshfine, gmeshcoarse, initVolForMergeMeshes,isMHM,needsMergeMeshes);
    
    // ----- Printing gmesh -----
#ifdef PZDEBUG
    if (1) {
        if(gmeshfine){
            gmeshfine->SetDimension(3);
            std::ofstream name(outputFolder + "GeoMesh_Fine_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
            std::ofstream name2(outputFolder + "GeoMesh_Fine_Initial.txt");
            gmeshfine->Print(name2);
        }
        if(gmeshcoarse){
            std::ofstream name(outputFolder + "GeoMesh_Coarse_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, name);
        }
    }
#endif
	
	// ----- Approximation space -----
    TMRSDataTransfer sim_data;
	sim_data.mTNumerics.m_four_approx_spaces_Q = true;
	sim_data.mTNumerics.m_mhm_mixed_Q = isMHM;
	sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;

	// ----- Filling data transfer from json input file -----
	FillDataTransferDFN(filenameBase, outputFolder, sim_data);
	
	// ----- Creating intersection GeoElBCs -----
	fixPossibleMissingIntersections(sim_data,gmeshfine); // read about this in the function itself
    
	// ----- Creating approximation space creator and manager -----
	// Code takes a fine and a coarse mesh to generate MHM data structure
    TMRSApproxSpaceGenerator aspace;
	aspace.InitMatIdForMergeMeshes() = initVolForMergeMeshes;
	
	// ----- Setting the global data transfer in approximation space -----
	aspace.SetDataTransfer(sim_data);

	// ----- Refining mesh (always done in fine mesh) -----
	if(isRefineMesh){
		cout << "\n---------------------- Uniformly refining geomesh ----------------------" << endl;
		gRefDBase.InitializeUniformRefPattern(ECube);
		gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
		gRefDBase.InitializeUniformRefPattern(EOned);
		for (auto gel : gmeshfine->ElementVec()){
			if(!gel || gel->HasSubElement()) continue;
			TPZManVector<TPZGeoEl*,10> children;
			gel->Divide(children);
		}
	}
	
	// ----- Setting gmesh -----
	// Code takes a fine and a coarse mesh to generate MHM data structure
	aspace.SetGeometry(gmeshfine,gmeshcoarse);

	{
        {
            int dim = gmeshfine->Dimension();
            int64_t nel = gmeshfine->NElements();
            int error = 0;
            for (int64_t el = 0; el<nel; el++) {
                TPZGeoEl *gel = gmeshfine->Element(el);
                if(gel->Dimension() != dim) continue;
                int firstside = gel->FirstSide(2);
                for(int is = firstside; is< gel->NSides()-1; is++)
                {
                    TPZGeoElSide gelside(gel,is);
                    auto neigh = gelside.Neighbour();
                    if(neigh == gelside) {
                        std::cout << gelside << "Has no neighbour\n";
                        error = 1;
                        TPZGeoElBC gelbc(gelside,-1000);
                    }
                }
            }
            if(error == 1){
                // NOTE: This may happen if setting the boundary conditions wrongly. I suggest checking them
                std::cout << "3D element without neighbour \n";
                std::ofstream name(outputFolder + "GeoMesh_Fine_AfterMergeMeshes.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
                DebugStop();
            }
        }
		std::ofstream name(outputFolder + "GeoMesh_Fine_AfterMergeMeshes.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
//        std::ofstream name2(outputFolder + "GeoMesh_MHM_domain.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name2,aspace.mSubdomainIndexGel);
	}
        
    // ----- Creates the multiphysics compmesh -----
	const int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true; // This makes solving very fast!
    bool UsingPzSparse = true; // Necessary to use multithread for now...
    bool UsePardiso_Q = true; // lighting fast!
    
    cout << "\n---------------------- Creating Analysis (Might optimize bandwidth) ----------------------" << endl;
    // -------------- Running problem --------------
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    if(isFilterZeroNeumann) FilterZeroNeumann(outputFolder,sim_data,mixAnalisys->StructMatrix(),mixed_operator);
    //        {
    //            std::ofstream out(outputFolder + "mixedCMesh.txt");
    //            mixed_operator->Print(out);
    //        }
    mixAnalisys->Assemble();
    
    // Testing if constant pressure leads to zero residual in cte pressure problem
    const bool testCtePressure = false;
    if(testCtePressure){
        const int neq = mixAnalisys->Mesh()->NEquations();
        TPZFMatrix<STATE> res(neq,1,0.);
        FillPCteSol(mixed_operator,1.);
        mixAnalisys->fsoltransfer.TransferFromMultiphysics();
        //			mixAnalisys->PostProcessTimeStep(2);
        //			mixAnalisys->PostProcessTimeStep(3);
        TPZMatrix<STATE>* mat = mixAnalisys->MatrixSolver<STATE>().Matrix().operator->();
        mat->Multiply(mixed_operator->Solution(), res);
        res = res + mixAnalisys->Rhs();
        std::ofstream out(outputFolder + "problematicEls.txt");
        mixAnalisys->PrintVectorByElement(out, res, 1.e-6);
    }
    
    // Solving problem
    mixAnalisys->Solve();
    mixAnalisys->VerifyElementFluxes();
    
    
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // The problem is linear, and therefore, we can just call assemble and solve once.
    // However, the system is assembled in a "nonlinear" fashion, thus, the solution represents
    // -DeltaU. So, to obtain the correct solution, we multiply it by -1.
    // Note: This has to be done after LoadSolution()!
    mixed_operator->UpdatePreviousState(-1.);
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    
    {
        // print the flux mesh with the solutions "loaded"
        TPZCompMesh *fluxmesh = mixed_operator->MeshVector()[0];
        std::ofstream flux("fluxmesh.txt");
        fluxmesh->Print(flux);
    }
    
    //        TPZCompMesh *pressure = mixed_operator->MeshVector()[1];
    //        pressure->Solution().Print("pressure multipliers");
    // Computes the integral of the normal flux on the boundaries.
    // To use, change the inletMatId and outletMatId according to problem
    
    
    // ----- Post processing -----
    if (isPostProc) {
        mixAnalisys->fsoltransfer.TransferFromMultiphysics();
        int dimToPost = 3;
        mixAnalisys->PostProcessTimeStep(dimToPost);
        dimToPost = 2;
        mixAnalisys->PostProcessTimeStep(dimToPost);
    }
    
    // ----- Compute integral of pressure and flux over domain and compare with analytical solution -----
    const std::string pvarname = "Pressure";
    const STATE integratedpressure = ComputeIntegralOverDomain(mixed_operator,pvarname);
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;

    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(mixed_operator,qvarname);
    if (fabs(integratedflux) < 1.e-12 ) integratedflux = 0.; // to make Approx(0.) work
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;
    
#ifndef USEMAIN
    // ----- Comparing with analytical solution -----
    if(isLinPressure){
        REQUIRE( integratedpressure == Approx( 4.0 ) ); // Approx is from catch2 lib
        REQUIRE( integratedflux == Approx( 8.*1./2.1) ); // Domain is 2x2x2 and flux is dp/dx = (1-0)/(1+0.1+1)
    }
    else{
        REQUIRE( integratedpressure == Approx( 8.0 ) );
        REQUIRE( integratedflux == Approx( 0.) );
    }
#endif
    
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()/1000.;
    cout << "\n\n\t--------- Total time of simulation = " << total_time << " seconds -------\n" << endl;

    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
}



// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransferDFN(string& filenameBase, string& outputFolder, TMRSDataTransfer& sim_data) {
	
	using json = nlohmann::json;
	std::string filenamejson = filenameBase + ".json";
	
	std::ifstream filejson(filenamejson);
	json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file

	
	// ------------------------ Getting number of domains and fractures ------------------------
	if(input.find("Domains") == input.end()) DebugStop();
	const int ndom = input["Domains"].size();
	if(input.find("Fractures") == input.end());
	const int nfrac = input["Fractures"].size();
	sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(ndom+nfrac+1);
	int countPhi = 0;
	
	// ------------------------ Reading 3D Domain matids ------------------------
	std::map<int, REAL> idPerm;
	for(auto& domain : input["Domains"]){
		if(domain.find("matid") == domain.end()) DebugStop();
		if(domain.find("name") == domain.end()) DebugStop();
		if(domain.find("K") == domain.end()) DebugStop();
		if(domain.find("phi") == domain.end()) DebugStop();
		const int matid = domain["matid"];
		const string name = domain["name"];
		const REAL permeability = domain["K"];
		const REAL phi = domain["phi"];
		sim_data.mTGeometry.mDomainNameAndMatId[name] = matid;
		idPerm[matid]= permeability;
		sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(matid,phi, 1.0);
	}
    
	sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
	
	// ------------------------ Reading 3D Domain BC matids ------------------------
	if(input.find("Boundary") == input.end()) DebugStop();
	std::map<int,std::pair<int,REAL>>& BCFlowMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue;
    std::map<int,std::pair<int,REAL>>& BCTransportMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue;
	for(auto& bc : input["Boundary"]){
		if(bc.find("matid") == bc.end()) DebugStop();
		if(bc.find("type") == bc.end()) DebugStop();
		if(bc.find("value") == bc.end()) DebugStop();
		const int matid = bc["matid"];
		const int type = bc["type"];
		const REAL value = bc["value"];
		
		if(BCFlowMatIdToTypeValue.find(matid) != BCFlowMatIdToTypeValue.end()) DebugStop();
		BCFlowMatIdToTypeValue[matid] = std::make_pair(type, value);
        BCTransportMatIdToTypeValue[matid] = std::make_pair(type, value);
        
	}
	
	// ------------------------ Reading fractures and fracture bcs matids ------------------------
    int initfracmatid =input["FractureInitMatId"];
    int actualfracid = initfracmatid;
    REAL phiintersec=0;
    REAL inersecLenght =0;
    for(auto& fracture : input["Fractures"]){
        const int i = fracture["Index"];
        std::string name = "Fracture" + std::to_string(i);
        const int matid = actualfracid;
        const REAL permerm = fracture["K"];
		const REAL fracWidth = fracture["width"];
//		const REAL fracWidth = 0.1;
		TMRSDataTransfer::TFracProperties::FracProp fracprop;
        sim_data.mTGeometry.mDomainFracNameAndMatId[name] = matid;
        actualfracid +=5;
		
		// Fracture properties
		fracprop.m_perm = permerm;
		fracprop.m_width = fracWidth;
		fracprop.m_fracbc.insert(matid + 1);
        fracprop.m_fracIntersectMatID = matid + 2;
        
		// Fracture polygon
        int npoints = fracture["Nodes"].size();
        TPZFMatrix<REAL> polygonmatrix(3,npoints);
        for(int j=0; j<npoints; j++){
            for(int k=0; k<3; k++){
                polygonmatrix(k,j) = (REAL)fracture["Nodes"][j][k];
            }
        }
        fracprop.m_polydata.SetCornersX(polygonmatrix);
        
		// Setting data structure
        sim_data.mTFracProperties.m_fracprops[matid] = fracprop;
        
		// Fracture boundary conditions
        const REAL zero_flux = 0., zero_pressure = 0.;
        const int bcFracType_Dirichlet = 0, bcFracType_Neumann = 1;
        if(fracprop.m_fracbc.size() == 1){
            sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[*fracprop.m_fracbc.begin()] = std::make_pair(bcFracType_Neumann, zero_flux);
        }
        else{
            DebugStop(); // for now, dfnimrs only allows one bc which is zero flux
        }
        sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[fracprop.m_fracIntersectMatID] = std::make_pair(bcFracType_Dirichlet, zero_pressure);
		
		// Setting fracture porosity in reservoir properties data structure
		// NOTE: Width is stored in two place for now. In mPorosityAndVolumeScale and in fracprop.m_width. Dangerous...
		if(fracture.find("phi") == fracture.end()) DebugStop();
		const REAL phifrac = fracture["phi"];
		sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(matid, phifrac, fracWidth);
        
        REAL phiintersec=phifrac;
        REAL inersecLenght =fracWidth;

    }
//    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(299,0.2, 0.01);
	

    int FractureHybridPressureMatId = input["FractureHybridPressureMatId"];
    // this is the material id of the pressure between hybridized fluxes of intersecting fractures
    sim_data.mTGeometry.m_pressureMatId = FractureHybridPressureMatId;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(FractureHybridPressureMatId, phiintersec, inersecLenght*inersecLenght);
    
    
    // the material id of the conductivity between fractures
    if(input.find("FractureGlueMatId") != input.end())
    {
        int fractureGlueMatId = input["FractureGlueMatId"];
        sim_data.mTFracIntersectProperties.m_FractureGlueId = fractureGlueMatId;
        if(input.find("FractureGluePerm") == input.end()) DebugStop();
        REAL fractureGluePerm = input["FractureGluePerm"];
        sim_data.mTFracIntersectProperties.m_FractureGluePerm = fractureGluePerm;
    }
	
	// ------------------------ Setting extra stuff that is still not in JSON ------------------------
	const int D_Type = 0, N_Type = 1, Mixed_Type = 2;
	sim_data.mTGeometry.mInterface_material_id = 100;
	sim_data.mTGeometry.mInterface_material_idFracInf = 102;
	sim_data.mTGeometry.mInterface_material_idFracSup = 101;
	sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
	sim_data.mTGeometry.mInterface_material_idFracBound = 104;
	
	// Other properties
	sim_data.mTGeometry.mSkeletonDiv = 0;
	sim_data.mTNumerics.m_sfi_tol = 0.0001;
	sim_data.mTNumerics.m_res_tol_transport = 0.0001;
	sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
	sim_data.mTNumerics.m_n_steps = 1 ;
	sim_data.mTNumerics.m_dt      = 1.0; //*day;
	sim_data.mTNumerics.m_four_approx_spaces_Q = true;
	std::vector<REAL> grav(3,0.0);
	grav[1] = 0.0;//-9.8*(1.0e-6); // hor
	sim_data.mTNumerics.m_gravity = grav;
	sim_data.mTNumerics.m_ISLinearKrModelQ = true;
	sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
	
    
    sim_data.mTNumerics.m_n_steps = 100;
	sim_data.mTNumerics.m_dt      = 0.01;//*day;
	sim_data.mTNumerics.m_max_iter_sfi=1;
	sim_data.mTNumerics.m_max_iter_mixed=1;
	sim_data.mTNumerics.m_max_iter_transport=1;
	
	// PostProcess controls
//	std::string vtkfilename = filenameBase.substr(filenameBase.find("dfnimrs/") + 8);
//	std::replace( vtkfilename.begin(), vtkfilename.end(), '/', '_');
    std::string vtkfilename = filenameBase.substr(filenameBase.find_last_of("/") + 1);
	std::string pressurevtk = outputFolder + vtkfilename + ".vtk";
	std::string transportvtk = outputFolder + vtkfilename + "_transport.vtk";
	std::cout << "\n===> PostProcess file name: " << vtkfilename << std::endl;
	
	
	
	sim_data.mTPostProcess.m_file_name_mixed = pressurevtk;
	sim_data.mTPostProcess.m_file_name_transport = transportvtk;
	TPZStack<std::string,10> scalnames, vecnames, scalnamesTransport;
	vecnames.Push("Flux");
	scalnames.Push("Pressure");
    scalnames.Push("div_q");
	if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
		scalnames.Push("g_average");
		scalnames.Push("p_average");
	}
	scalnamesTransport.Push("Sw");
	scalnamesTransport.Push("So");
	
	sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
	sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
	sim_data.mTPostProcess.m_scalnamesTransport = scalnamesTransport;
	
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

void MapFractureIntersection(const std::string &filenameBase, std::map<std::string,int> &matmap, int firstfracintersect, std::map<int,std::pair<int,int>> &matidtoFractures)
{
    // Creating gmsh reader
    TPZGmshReader  Geometry;
    std::string filename = filenameBase + "_fine.msh";
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


void ReadMeshesDFN(string& filenameBase, TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int& initVolForMergeMeshes, bool& isMHM, bool& needsMergeMeshes) {
	using json = nlohmann::json;
	std::string filenamejson = filenameBase + ".json";
	
	std::ifstream filejson(filenamejson);
    if(!filejson) DebugStop();
	json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file

	// ===================> Coarse mesh <=======================
	TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
	TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
	
	std::set<int> allmatids; // used to check for repeated matids and highest matid
	
    // ------------------------ Check if isMHM is set on file ------------------------
    isMHM = true;
    if(input.find("useMHM") != input.end()) {
        isMHM = input["useMHM"];
    }
    
    if(input.find("needsMerge") != input.end()) {
        needsMergeMeshes = input["needsMerge"];
    }
    
    if(isMHM && !needsMergeMeshes) DebugStop();
    
    int SimulationDim= 3;
    if(input.find("SimulationDim") != input.end()) {
        SimulationDim = input["SimulationDim"];
    }
    
	// ------------------------ Get matids of 3D domain ------------------------
	if(input.find("Domains") == input.end()) DebugStop();
	for(auto& domain : input["Domains"]){
		if(domain.find("matid") == domain.end()) DebugStop();
		if(domain.find("name") == domain.end()) DebugStop();
		const int matid = domain["matid"];
		const string name = domain["name"];
		dim_name_and_physical_tagCoarse[SimulationDim][name] = matid;
		const bool is_in = allmatids.find(matid) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(matid);
	}
	
	// ------------------------ Get matids of BCs of 3D domain ------------------------
	// NOTE: This might not be needed for the coarse mesh. Please check...
	if(input.find("Boundary") == input.end()) DebugStop();
	for(auto& bc : input["Boundary"]){
		if(bc.find("name") == bc.end()) DebugStop();
		if(bc.find("matid") == bc.end()) DebugStop();
		const string name = bc["name"];
		const int matid = bc["matid"];
		dim_name_and_physical_tagCoarse[SimulationDim-1][name] = matid;
		dim_name_and_physical_tagFine[SimulationDim-1][name] = matid;
		const bool is_in = allmatids.find(matid) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(matid);
	}

	// ------------------------ Generate gmesh coarse ------------------------
    if(input.find("Mesh") == input.end()) DebugStop();
    std::string meshfile = input["Mesh"];
    {
        auto lastslash = filenamejson.find_last_of("/");
        auto meshdirname = filenamejson.substr(0,lastslash+1);
//        meshFile = meshFile.substr(meshFile.find("examples/") + 9,meshFile.length());
        meshfile = meshdirname + meshfile;
    }
    int ncoarse_vol = 0;
    if(needsMergeMeshes){
        gmeshcoarse = generateGMeshWithPhysTagVec(meshfile,dim_name_and_physical_tagCoarse);
        
        int64_t nelcoarse = gmeshcoarse->NElements();
        for(int64_t el = 0; el<nelcoarse; el++)
        {
            TPZGeoEl *gel = gmeshcoarse->Element(el);
            if(gel && gel->Dimension()==3) ncoarse_vol++;
        }
    }
	
//    string filenameCoarse = filenameBase + "_coarse.msh";
//    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);

   
	// ===================> Fine mesh <=======================
	// Note that boundary elements have been added previously to the set dim_name_and_physical_tagFine

	if(input.find("FractureInitMatId") == input.end()) DebugStop();
	const int fracInitMatId = input["FractureInitMatId"];
    const int fracinc = 5;
	
	// ------------------------ Loop over fractures in Json ------------------------
	int fracCounter = 0;
    bool isFracSim=false;
    if(input.find("Fractures") != input.end()){
	for(auto& frac : input["Fractures"]){
		const int matid = fracInitMatId + fracCounter*fracinc;
        const int bcmatid = fracInitMatId + fracCounter*fracinc+1;
		string fracname = "Fracture" + to_string(fracCounter);
		string bcfracname = "BCfrac" + to_string(fracCounter);
		const int currentFracId = fracInitMatId + fracinc * fracCounter;
        dim_name_and_physical_tagFine[2][fracname] = currentFracId;
        dim_name_and_physical_tagFine[1][bcfracname] = (currentFracId)+1;
        bool is_in = allmatids.find(currentFracId) != allmatids.end();
        if(is_in) DebugStop();
        allmatids.insert((currentFracId));
        is_in = allmatids.find((currentFracId)+1) != allmatids.end();
        if(is_in) DebugStop();
        allmatids.insert((currentFracId)+1);
		fracCounter++;
	}
        isFracSim=true;
    }
	
	// ------------------------ Adding volume physical tags------------------------
	const int maxMatId = *allmatids.rbegin();
	initVolForMergeMeshes = (1+maxMatId/100)*100;
	if(input.find("NCoarseGroups") == input.end()) DebugStop();
	const int nCoarseGroups = input["NCoarseGroups"];
    
    std::string volbase = "c";
    if(needsMergeMeshes){
        if(nCoarseGroups != ncoarse_vol) DebugStop();
        for (int ivol = 0; ivol < nCoarseGroups; ivol++) {
            std::string ivolstring = volbase + to_string(ivol);
            dim_name_and_physical_tagFine[SimulationDim][ivolstring] = initVolForMergeMeshes + ivol;
        }
    }
	// ------------------------ Generate gmesh fine ------------------------
	string filenameFine = filenameBase + "_fine.msh";
    
    /// add the intersection material ids to be read
    /// read only the header of the msh file
    /// identify the intersection groups by identify the substring
    if(input.find("FractureHybridPressureMatId") == input.end()) DebugStop();
    int FractureHybridPressureMatId = input["FractureHybridPressureMatId"];
    
    int firstfracintersect = initVolForMergeMeshes + ncoarse_vol;
    firstfracintersect = (1+firstfracintersect/100)*100;
    std::map<int,std::pair<int,int>> matidtoFractures;
    
    if(isFracSim){
        MapFractureIntersection(filenameBase, dim_name_and_physical_tagFine[1], firstfracintersect, matidtoFractures);
    }
    
    if (!isMHM && !needsMergeMeshes) {
        gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagCoarse);
    }
    else{
        gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    }

    // for each intersection element create an intersection element for each fracture specifically
    if (isFracSim) {
        CreateIntersectionElementForEachFrac(gmeshfine,matidtoFractures,fracInitMatId,fracinc,FractureHybridPressureMatId);
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
			int matid1 = fracInitMatId + frac1*fracinc + 2;
			int frac1bcid = fracInitMatId + frac1*fracinc + 1;
			auto frac1bc = gelside.HasNeighbour(frac1bcid);
			if(frac1bc)
			{
				frac1bc.Element()->SetMaterialId(matid1);
			}
			else
			{
				TPZGeoElBC(gelside,matid1);
			}
			int frac2 = it->second.second;
			int matid2 = fracInitMatId + frac2*fracinc + 2;
			int frac2bcid = fracInitMatId + frac2*fracinc + 1;
			auto frac2bc = gelside.HasNeighbour(frac2bcid);
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

void fixPossibleMissingIntersections(TMRSDataTransfer& sim_data, TPZGeoMesh* gmesh){
    // it may happen that the supplied mesh from DFN has missing intersections where 2 fractures intersect.
    // We treat this here. But it would be ideal to make DFN robust enough so we could remove this function whatsoever
    // The main idea is to check if a fracture element has more than 1 fracture neightbor. If so, an intersection element is
    // needed to hybridzie the region.
    cout << "\n---------------------- Searching for problematic intersection elements ----------------------" << endl;
    
    int nInterCreated = 0;
    for (auto gel : gmesh->ElementVec()) {
        if (!gel || !sim_data.mTFracProperties.isFracMatId(gel->MaterialId())) {
            continue;
        }
        int gelmatid = gel->MaterialId();
        if (gel->Dimension() != 2) {
            DebugStop(); // Should it work with 2D problems (1d fractures)?
        }
		
        auto fracprop = sim_data.mTFracProperties.m_fracprops[gelmatid];
//        const int matidintersect = sim_data.mTFracIntersectProperties.m_IntersectionId;
        const int matidintersect = fracprop.m_fracIntersectMatID;

        const int firstedge = gel->FirstSide(1);
        const int lastedge = gel->FirstSide(2);
        // loop over the sides of dimension 1
        for (int iside = firstedge; iside < lastedge; iside++) {
            TPZStack<TPZGeoEl*> interEls;
            TPZStack<TPZGeoEl*> bcEls;
            int nFracElsForSide = 0;
            // loop over the neighbours
            TPZGeoElSide gside(gel,iside);
            TPZGeoElSide neig = gside.Neighbour();
            for (; gside != neig; neig++) {
                TPZGeoEl* neigel = neig.Element();
				const int neighelmatid = neigel->MaterialId();
                if (neighelmatid == matidintersect) {
                    interEls.push_back(neigel);
                }
                // compute the number of fracture elements linked to this side
                if (sim_data.mTFracProperties.isFracMatId(neighelmatid)) {
                    nFracElsForSide++;
                }
                // compute the number of boundary condition elements linked to this side
                if (sim_data.mTFracProperties.isFracBCMatId(neighelmatid)) {
                    bcEls.push_back(neigel);
                }
            }
//            if(bcEls.size() > 1) DebugStop();
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
            if (nFracElsForSide > 1 && !interEls.size() && !bcEls.size()) {
                cout << "nfracs for this side: " << nFracElsForSide << endl;
                cout << "This is an error" << endl;
//                DebugStop();
            }
            if (0 && (interEls.size() || nFracElsForSide > 1) && bcEls.size()) {
                cout << "PLEASE CHECK CAREFULLY! An element may have a boundary even if it intersects\n";
                cout << "Domains intersecting and boundary through snap\n Intersection matids ";
                for(auto it : interEls) cout << it->MaterialId() << " ";
                cout << "\n Boundary condition material ids ";
                for(auto it : bcEls) cout << it->MaterialId() << " ";
                cout << std::endl;
                cout << "matid of all neighbours " << gel->MaterialId();
                neig = gside.Neighbour();
                for (; gside != neig; neig++) {
                    TPZGeoEl* neigel = neig.Element();
                    const int neighelmatid = neigel->MaterialId();
                    cout << " " << neighelmatid;
                }
                cout << std::endl;
//                if(bcEls.size() > 1)
//                    DebugStop(); // This is rather odd. Please check why there are two bcs at the same boundary
//                for (auto bcgeoel : bcEls) {
//                    bcgeoel->RemoveConnectivities();
//                    delete bcgeoel;
//                }
            }
            
            
        }
    }
    if (nInterCreated) {
        cout << "\n\t\t==================== MESSAGE =================" << endl;
        cout << "\n\t- Imported mesh has " << nInterCreated << " missing intersection elements" << endl;
        cout << "\t- These were created here and this shouldn't be a problem...\n" << endl;
    }
    else {
        cout << "\n====> Ok! No problematic intersections found" << endl;
    }

	// No need to buildConnectivity. It is already correct
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillPCteSol(TPZMultiphysicsCompMesh* mpcmesh, const REAL pcte) {
	
	int64_t nc = mpcmesh->NConnects();
	TPZFMatrix<STATE> &sol = mpcmesh->Solution();
	for (int64_t ic = 0; ic<nc; ic++) {
		TPZConnect &c = mpcmesh->ConnectVec()[ic];
		int64_t seqnum = c.SequenceNumber();
		if(seqnum < 0) continue;
		unsigned char lagrange = c.LagrangeMultiplier();
		STATE fill = 0.;
		if(lagrange == 1 || lagrange == 3 || lagrange == 5)
		{
			fill = pcte;
		}
		int ndof = c.NShape();
		for (int idf = 0; idf < ndof ; idf++) {
			int64_t index = mpcmesh->Block().Index(seqnum, idf);
			sol(index,0) = fill;
		}

	}
		
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

std::map<int,TPZVec<STATE>> computeIntegralOfNormalFlux(const std::set<int> &bcMatId, TPZMultiphysicsCompMesh *cmesh) {
	
	// modifying bc of top outlet of case 3 to compute integral separately
	TPZGeoMesh* gmesh = cmesh->Reference();
	const REAL zerotol = ZeroTolerance();
    std::map<int,TPZVec<STATE>> result;
	std::set<int> matidsInlet, matidsOutlet;
	std::string varname="BCNormalFlux";
	cmesh->Reference()->ResetReference();
	cmesh->LoadReferences(); // compute integral in the multiphysics mesh
	int nels = cmesh->NElements();
	for(auto it: bcMatId)
    {
        std::set<int> matid = {it};
        TPZVec<STATE> vecint = cmesh->Integrate(varname, matid);
        result[it] = vecint;
        std::cout << "Integral of normal flux for matid " << it << " = " << vecint << std::endl;
        
    }
    return result;
	
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

bool fileExists(const fs::path& p, fs::file_status s) {
    if(fs::status_known(s) ? fs::exists(s) : fs::exists(p))
        return true;
    else
        return false;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateOutputFolders(std::string& outputFolder) {

    std::string folders = outputFolder;
    char c = '/';
    std::string folderToCreate = "";
    int nfolders = 0;
    while (folders.find("/") != std::string::npos) {
        if(nfolders == 0) folderToCreate = folders.substr(0,folders.find("/"));
        else folderToCreate = folderToCreate + "/" + folders.substr(0,folders.find("/"));
        folders = folders.substr(folders.find("/")+1);
        if(!fileExists(folderToCreate)){
            if (!fs::create_directory(folderToCreate))
                DebugStop();
            else
                cout << "Directory created with name " << folderToCreate << endl;
        }
        nfolders++;
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CopyInputFilesToOutputFolderAndFixFilename(std::string& filenameBase, std::string& outputFolder){
    
    // Copying all files from input folder to output folder
    std::string onlyFolder = outputFolder.substr(0,outputFolder.find_last_of("/"));
    fs::copy(filenameBase , onlyFolder, fs::copy_options::update_existing | fs::copy_options::recursive);
    
    // Adding the json name to filenameBase
    int njson = 0;
    for (auto const& dir_entry : std::filesystem::directory_iterator{filenameBase}){
        std::string filename = dir_entry.path().string();
//        std::cout << dir_entry.path().string() << endl;
        
        if (filename.find(".") != std::string::npos){
            if(filename.substr(filename.find(".")) == ".json"){
                if(njson == 0){
                    filenameBase = filename.substr(0,filename.find("."));
                    njson++;
                }
                else {
                    cout << "\n\n=====> ERROR! There are two json files in the provided input folder" << endl;
                    DebugStop();
                }
            }
        }
    }

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

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FilterZeroNeumann(std::string& outputFolder, TMRSDataTransfer& sim_data, TPZAutoPointer<TPZStructMatrix> strmat, TPZCompMesh* cmesh) {

    std::cout << "\n---------------------- Filtering zero neumann equations ----------------------" << std::endl;
    TPZSimpleTimer timer_filter("Timer Filter Equations");
    
   
    std::set<int64_t> matidset;
    
    // First find all the zero neumann in in the 3d domain
    std::cout << "Domain BC matids: ";
    for(auto &chunk : sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
        const int bc_id = chunk.first;
        const std::pair<int,REAL>& typeAndVal = chunk.second;
        const int bc_type = typeAndVal.first;
        const REAL val = typeAndVal.second;
        if(bc_type == 1 && fabs(val) < ZeroTolerance()){
            std::cout << bc_id << " ";
            matidset.insert(bc_id);
        }
    }
    
    // Then all the zero neumann in the fractures
    std::cout << "\nFracture BC matids: ";
    for (auto& chunk : sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
        const int bc_id   = chunk.first;
        const std::pair<int,REAL>& typeAndVal = chunk.second;
        const int bc_type = typeAndVal.first;
        const REAL val = typeAndVal.second;
        if(bc_type == 1 && fabs(val) < ZeroTolerance()){
            std::cout << bc_id << " ";
            matidset.insert(bc_id);
        }
    }
    
    // Set equations to be filted
    std::set<int64_t> eqset;
    cmesh->GetEquationSetByMat(matidset, eqset);
    if (eqset.size()) {
        strmat->EquationFilter().ExcludeEquations(eqset);
    }
    
    int count = 0;
    for(auto cel : cmesh->ElementVec()){
        TPZSubCompMesh* subcmesh = dynamic_cast<TPZSubCompMesh*>(cel);
        if (subcmesh) {
            std::cout << "\n\t------- Submesh " << count << " -------" << std::endl;
//            if(count == 1 || count == 3){
//                count++;
//                std::cout << "===> Skipping filter" << std::endl;
//                continue;
//            }
            count++;
#ifdef PZDEBUG
            {
                std::string filename = outputFolder + "submesh_" + to_string(count) + ".vtk";
                std::ofstream out(filename);
                TPZVTKGeoMesh::PrintCMeshVTK(subcmesh, out);
            }
#endif
            std::set<int64_t> eqsetsub;
            subcmesh->GetEquationSetByMat(matidset, eqsetsub);
            const int64_t ninteq = subcmesh->NumInternalEquations();
            const int64_t neq = subcmesh->TPZCompMesh::NEquations();
            if (eqsetsub.size()) {
                subcmesh->Analysis()->StructMatrix()->EquationFilter().SetNumEq(neq); // Setting again bcz els could have been condensed
                subcmesh->Analysis()->StructMatrix()->EquationFilter().ExcludeEquations(eqsetsub); // will set as active all eqs (internal and external) that are not zero neumann
                auto& activeeq = subcmesh->Analysis()->StructMatrix()->EquationFilter().GetActiveEquations();
                
                std::cout << "size eqsetsub = " << eqsetsub.size() << " | eqsetsub = " << eqsetsub;
                std::cout << "size activeeq = " << activeeq.size() << " | active eq = " << activeeq << std::endl;
                auto& excludedeq = subcmesh->Analysis()->StructMatrix()->EquationFilter().GetExcludedEquations();
                std::cout << "size excludedeq = " << excludedeq.size() << " | excluded eq = " << excludedeq << std::endl;
                const int64_t nexteq = neq - ninteq;
#ifdef PZDEBUG
                // All the external equations must be present in the active equations
                int64_t count = 0;
                for(auto& eq : activeeq){
                    if(eq > ninteq-1) count++;
                }
                if(count < nexteq) DebugStop();
#endif
                // Now we want to set as active only the ones that are internal and not zero Neumann
                // Since activeeq is ordered, we simply pick equations until we reach the first non internal equation
                TPZVec<int64_t> activeinternal(ninteq);
                int64_t i = 0;
                for(auto& eq : activeeq){
                    if(eq > ninteq-1) break;
                    activeinternal[i++] = eq;
                }
                activeinternal.Resize(i);
                subcmesh->Analysis()->StructMatrix()->EquationFilter().Reset();
                subcmesh->Analysis()->StructMatrix()->EquationFilter().SetActiveEquations(activeinternal); // sets new filter
                // Note that, in subcmesh, after we create the data structure of the matred matrices, we need to the
                // set the external equations as active again
            }
        }
    }
    
    std::cout << "\n==> Total Filter time: " << timer_filter.ReturnTimeDouble()/1000. << " seconds" << std::endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
