//
//  TMRSSFIAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSSFIAnalysis.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZDarcyMemory.h"
#include "TPZFastCondensedElement.h"
#include "TPZDarcyFlowWithMem.h"
#include "TPZCompElHDivCollapsed.h"
#include <pzshapequad.h>
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "TPZSimpleTimer.h"

TMRSSFIAnalysis::TMRSSFIAnalysis(){
    m_sim_data = nullptr;
    m_k_iteration = 0;
    m_mixed_module = nullptr;
    m_transport_module = nullptr;
    m_x_mixed.Resize(0, 0);
    m_x_transport.Resize(0, 0);
}

TMRSSFIAnalysis::~TMRSSFIAnalysis(){
    
}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZCompMesh * cmesh_transport,
                                 const RenumType& renumtype){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,renumtype);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,renumtype);
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    
    
}

void TMRSSFIAnalysis::BuildAlgebraicDataStructure(){
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    FillProperties();
    fAlgebraicDataTransfer.TransferSaturation();
    m_mixed_module->SetLastStateVariables();
}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed,
                                 TPZMultiphysicsCompMesh * cmesh_transport,
                                 std::function<REAL(const TPZVec<REAL> & )> & kx,
                                 std::function<REAL(const TPZVec<REAL> & )> & ky,
                                 std::function<REAL(const TPZVec<REAL> & )> & kz,
                                 std::function<REAL(const TPZVec<REAL> & )> & phi,
                                 std::function<REAL(const TPZVec<REAL> & )> & s0,
                                 const RenumType& renumtype){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,renumtype);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,renumtype);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    bool propsfromPre = m_sim_data->mTReservoirProperties.fPropsFromPreProcess;
    if (propsfromPre==false) {
            fAlgebraicDataTransfer.fkx = kx;
            fAlgebraicDataTransfer.fky = ky;
            fAlgebraicDataTransfer.fkz = kz;
            fAlgebraicDataTransfer.fphi = phi;
            fAlgebraicDataTransfer.fs0 = s0;
    }

    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
//    FillMaterialMemoryDarcy(1);
//    std::string fileprops("Props.txt");
//    std::ofstream file(fileprops);
//    if(file && !propsfromPre){
//        FillProperties(fileprops, &m_transport_module->fAlgebraicTransport);
//    }

}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed,
                                 TPZMultiphysicsCompMesh * cmesh_transport,
                                 std::function<std::vector<REAL>(const TPZVec<REAL> & )> & kappa_phi,
                                 std::function<REAL(const TPZVec<REAL> & )> & s0,
                                 const RenumType& renumtype){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,renumtype);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,renumtype);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    
    fAlgebraicDataTransfer.fkappa_phi = kappa_phi;
    fAlgebraicDataTransfer.fs0 = s0;
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    
    FillMaterialMemoryDarcy(1);
    std::string fileprops("Props.txt");
    FillProperties(fileprops, &m_transport_module->fAlgebraicTransport);
}

void TMRSSFIAnalysis::Configure(int n_threads, bool UsePardiso_Q, bool usepz){
    m_mixed_module->Configure(n_threads, UsePardiso_Q, usepz);
    m_transport_module->Configure(n_threads, UsePardiso_Q);
    
}
void TMRSSFIAnalysis::FillMaterialMemoryDarcy(int material_id){
  
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }
    // Initialize integration points memory
    
    TPZCompMesh * cmesh = m_mixed_module->Mesh();
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    TPZDarcyFlowWithMem *mat = dynamic_cast<TPZDarcyFlowWithMem *>(material);

    if (!material)
        DebugStop();

    TPZMatWithMem<TPZDarcyMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPZDarcyMemory> * >(material);
    if (!mat_with_memory)
        DebugStop();

    std::shared_ptr<TPZAdmChunkVector<TPZDarcyMemory>> & memory_vector = mat_with_memory->GetMemory();

//    for(auto &meshit : fAlgebraicDataTransfer.fTransportMixedCorrespondence) {
//        int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
//#ifdef PZDEBUG
//        if(meshit.fTransport == 0)
//            DebugStop();
//
//#endif
//        for (int icell = 0; icell < ncells; icell++) {
//            int64_t algbindex = meshit.fAlgebraicTransportCellIndex[icell];
//            TPZFastCondensedElement * fastCond = meshit.fMixedCell[icell];
//            if(fastCond->Reference()->MaterialId()!=material_id){
//                continue;
//            }
//            TPZCompEl *celcomp = fastCond->ReferenceCompEl();
//
//
//            TPZManVector<int64_t> indices;
//            celcomp->GetMemoryIndices(indices);
//            indices.Print();
//            for (int index = 0; index<indices.size(); index++) {
//                int valIndex = indices[index];
//                TPZDarcyMemory &mem = memory_vector.get()->operator [](valIndex);
//                mem.fTransportCellIndex = algbindex;
//                m_transport_module->fAlgebraicTransport.fCellsData.fKx[algbindex] = 1.0;
//                m_transport_module->fAlgebraicTransport.fCellsData.fKy[algbindex] = 1.0;
//                m_transport_module->fAlgebraicTransport.fCellsData.fKz[algbindex] = 1.0;
//                m_transport_module->fAlgebraicTransport.fCellsData.fporosity[algbindex] = 0.1;
//            }
//        }
//    }
    
}
void TMRSSFIAnalysis::FillProperties(){
    
    bool propsfromPre = m_sim_data->mTReservoirProperties.fPropsFromPreProcess;
    if (propsfromPre==false) {
        if (m_sim_data->mTReservoirProperties.kappa_phi) {
            fAlgebraicDataTransfer.fkappa_phi = m_sim_data->mTReservoirProperties.kappa_phi;
            fAlgebraicDataTransfer.fs0=m_sim_data->mTReservoirProperties.s0;
            FillProperties(&m_transport_module->fAlgebraicTransport);
        }
        else{
            std::vector<REAL> kappa_phi(4,1.0);

            if (m_sim_data->mTReservoirProperties.kappa_phi){
                fAlgebraicDataTransfer.fkappa_phi = m_sim_data->mTReservoirProperties.kappa_phi;
            }
            if (m_sim_data->mTReservoirProperties.s0){
                fAlgebraicDataTransfer.fs0 = m_sim_data->mTReservoirProperties.s0;
            }
            
            m_transport_module->fAlgebraicTransport.interfaceid = m_sim_data->mTGeometry.mInterface_material_id;
            m_transport_module->fAlgebraicTransport.fgravity = m_sim_data->mTNumerics.m_gravity;
            
            //Set initial properties
            m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[0] = m_sim_data->mTFluidProperties.mWaterViscosity;
            m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[1] = m_sim_data->mTFluidProperties.mOilViscosity;
            int ncells = m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil.size();
            REAL rhow = m_sim_data->mTFluidProperties.mWaterDensityRef;
            REAL rhoo = m_sim_data->mTFluidProperties.mOilDensityRef;
            m_transport_module->fAlgebraicTransport.initialMass = 0.0;
          
            for (int icell =0; icell<ncells; icell++) {
                m_transport_module->fAlgebraicTransport.fCellsData.fDensityWater[icell]= rhow; m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil[icell]= rhoo;
                int matid = m_transport_module->fAlgebraicTransport.fCellsData.fMatId[icell];
                bool fountmat = false;
                for (auto i:m_sim_data->mTReservoirProperties.mPorosityAndVolumeScale) {
                    int id = std::get<0>(i);
                    REAL porosity =std::get<1>(i);
                    REAL volfactor =std::get<2>(i);
                    if(id == matid){
                        m_transport_module->fAlgebraicTransport.fCellsData.fporosity[icell] = porosity;
                        m_transport_module->fAlgebraicTransport.fCellsData.fVolumefactor[icell] = volfactor;
                        m_transport_module->fAlgebraicTransport.fCellsData.fVolume[icell] *= volfactor;
                        fountmat =true;
                        break;
                    }
                }
                REAL satW = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation[icell];
                REAL satO = 1.0 - satW;
                REAL rhoW = m_transport_module->fAlgebraicTransport.fCellsData.fDensityWater[icell];
                REAL rhoO = m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil[icell];
                REAL phi = m_transport_module->fAlgebraicTransport.fCellsData.fporosity[icell];
                REAL vol = m_transport_module->fAlgebraicTransport.fCellsData.fVolume[icell];
                m_transport_module->fAlgebraicTransport.initialMass += (satW*rhoW+satO*rhoO) * phi * vol;

                REAL kappa = m_sim_data->mTReservoirProperties.m_permeabilitiesbyId[matid]; //if no permeability function is set, all cells must have the same value (for IHU)
                m_transport_module->fAlgebraicTransport.fCellsData.fKx[icell] = kappa;
                m_transport_module->fAlgebraicTransport.fCellsData.fKy[icell] = kappa;
                m_transport_module->fAlgebraicTransport.fCellsData.fKz[icell] = kappa;
                
                if (!fountmat){
                    DebugStop();
                }
//
                
            }
            m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensitiesLastState();
            m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensities();
            
            m_transport_module->fAlgebraicTransport.fdt = m_sim_data->mTNumerics.m_dt;
            bool foundinlet = false;
			for (auto& chunk : m_sim_data->mTBoundaryConditions.mBCTransportMatIdToTypeValue) {
				const int idVal   = chunk.first;
				std::pair<int,REAL>& typeAndVal = chunk.second;
				const int idType = typeAndVal.first;
				const REAL idValue   = typeAndVal.second;
				std::pair<int, REAL> bccond = std::make_pair(idType, idValue);
				m_transport_module->fAlgebraicTransport.fboundaryCMatVal[idVal] =bccond;
			}
        }
    }
    else{
        if (m_sim_data->mTReservoirProperties.mPropsFileName=="") {
            std::string propsname ="PreProcess/props/"+ m_sim_data->mSimulationName + "Props_nL_"+ std::to_string(m_sim_data->mTGeometry.mnLayers)  +"_nRef_"+std::to_string(m_sim_data->mTGeometry.mnref)+".txt" ;
            m_sim_data->mTReservoirProperties.mPropsFileName = propsname;
        }
        
        std::ifstream file(m_sim_data->mTReservoirProperties.mPropsFileName);
        
        if(file){
            FillProperties(m_sim_data->mTReservoirProperties.mPropsFileName, &m_transport_module->fAlgebraicTransport);
        }
        else{
            
            std::cout<<"The properties of the reservoir have not been loaded. Please enter the properties in a text file or set the functions in the simulationdata object"<<std::endl;
            std::vector<REAL> kappa_phi(4,1.0e-7);
            kappa_phi[3]=0.1;
            FillProperties(&m_transport_module->fAlgebraicTransport, kappa_phi);
        }
        
    }
}
void TMRSSFIAnalysis::FillProperties(std::string fileprops, TPZAlgebraicTransport *algebraicTransport){
    
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }
    std::vector<REAL> Kx, Ky, Kz, Phi;
    ReadProperties(fileprops, false, Kx, Ky, Kz, Phi);
    
    TPZCompMesh * cmesh = m_mixed_module->Mesh();
    
//    for(auto &meshit : fAlgebraicDataTransfer.fTransportMixedCorrespondence)
//    {
//        int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
//#ifdef PZDEBUG
//        if(meshit.fTransport == 0)
//        {
//            DebugStop();
//        }
//#endif
//        for (int icell = 0; icell < ncells; icell++)
//        {
//            int64_t algbindex = meshit.fAlgebraicTransportCellIndex[icell];
//            TPZFastCondensedElement * fastCond = meshit.fMixedCell[icell];
//            TPZCompEl *celcomp = fastCond->ReferenceCompEl();
//            int geoIndexMixed = celcomp->Reference()->Index();
//            if (Kx[geoIndexMixed] <0) {
//                DebugStop();
//            }
//            algebraicTransport->fCellsData.fKx[algbindex] = Kx[geoIndexMixed];
//            algebraicTransport->fCellsData.fKy[algbindex] = Ky[geoIndexMixed];
//            algebraicTransport->fCellsData.fKz[algbindex] = Kz[geoIndexMixed];
//            algebraicTransport->fCellsData.fporosity[algbindex] = Phi[geoIndexMixed];
//        }
//    }
    m_transport_module->fAlgebraicTransport.fHasPropQ=true;
    
}
void TMRSSFIAnalysis::FillProperties(TPZAlgebraicTransport *algebraicTransport){
    
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }

    
    TPZCompMesh * cmesh = m_mixed_module->Mesh();
    
//    for(auto &meshit : fAlgebraicDataTransfer.fTransportMixedCorrespondence)
//    {
//        int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
//#ifdef PZDEBUG
//        if(meshit.fTransport == 0)
//        {
//            DebugStop();
//        }
//#endif
//        for (int icell = 0; icell < ncells; icell++)
//        {
//            int64_t algbindex = meshit.fAlgebraicTransportCellIndex[icell];
//            TPZFastCondensedElement * fastCond = meshit.fMixedCell[icell];
//            TPZCompEl *celcomp = fastCond->ReferenceCompEl();
//            TPZGeoEl *gel = celcomp->Reference();
//            int dim= gel->Dimension();
//            TPZVec<REAL> ximasscent(dim);
//            gel->CenterPoint(gel->NSides()-1, ximasscent);
//            std::vector<REAL> center(3,0.0);
//            TPZManVector<REAL,3> coord(3,0.0);
//            gel->X(ximasscent, coord);
//            
//            std::vector<REAL> kappa_phi = m_sim_data->mTReservoirProperties.kappa_phi(coord);
//                algebraicTransport->fCellsData.fKx[algbindex]  = kappa_phi[0];
//                algebraicTransport->fCellsData.fKx[algbindex]  = kappa_phi[1];
//                algebraicTransport->fCellsData.fKx[algbindex]  = kappa_phi[2];
//                algebraicTransport->fCellsData.fKx[algbindex] = kappa_phi[3] + 0.01;
//        }
//    }
    m_transport_module->fAlgebraicTransport.fHasPropQ=true;
}

void TMRSSFIAnalysis::FillProperties(TPZAlgebraicTransport *algebraicTransport, std::vector<REAL> kappa_phi){
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }
    
    int ncells = algebraicTransport->fCellsData.fVolume.size();
    for (int icell =0; icell<ncells; icell++) {
        algebraicTransport->fCellsData.fKx[icell]  = kappa_phi[0];
        algebraicTransport->fCellsData.fKy[icell]  = kappa_phi[1];
        algebraicTransport->fCellsData.fKz[icell]  = kappa_phi[2];
//        algebraicTransport->fCellsData.fporosity[icell] = kappa_phi[3] ;
    }
    
   
    m_transport_module->fAlgebraicTransport.fHasPropQ=true;
    
}


void TMRSSFIAnalysis::SetDataTransferAndBuildAlgDatStruct(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    m_mixed_module->SetDataTransfer(sim_data);
    m_transport_module->SetDataTransfer(sim_data);
    m_transport_module->fAlgebraicTransport.fCellsData.SetDataTransfer(sim_data);
    
    // Creates the interface data structure which holds data such as integrated fluxes at
    // interfaces and information of neighborhood
    // Also creats the TCellData which holds information such as the saturation, density,
    // porosity, etc. in each element
    BuildAlgebraicDataStructure();
    
    
    m_transport_module->AnalyzePattern();
}

TMRSDataTransfer * TMRSSFIAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSSFIAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSSFIAnalysis::RunTimeStep(){
    
    m_x_mixed = m_mixed_module->Solution();
    m_x_transport = m_transport_module->Solution();
    
    int n_iterations = m_sim_data->mTNumerics.m_max_iter_sfi;
    REAL eps_tol = m_sim_data->mTNumerics.m_sfi_tol;
    bool stop_criterion_Q = false;
    REAL error_rel_mixed = 1.0;
    REAL error_rel_transport = 1.0;
    
    for (m_k_iteration = 1; m_k_iteration <= 1; m_k_iteration++) {
        std::cout << "\nSFI iteration: " << m_k_iteration << std::endl;
        
        SFIIteration();        
        error_rel_mixed = Norm(m_x_mixed - m_mixed_module->Solution());
        
        if(iszero(Norm(m_transport_module->Solution()))){
            error_rel_transport = Norm(m_x_transport - m_transport_module->Solution());
        }else{
            error_rel_transport = Norm(m_x_transport - m_transport_module->Solution());
        }

        stop_criterion_Q = m_sim_data->mTNumerics.m_is_linearTrace? true : (error_rel_transport < eps_tol && error_rel_mixed < eps_tol); // Stop by saturation variation
        if (stop_criterion_Q && m_k_iteration >= 1) {
            std::cout << "SFI converged with error_rel_transport: " << error_rel_transport << " and error_rel_mixed: " << error_rel_mixed << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
            m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensitiesLastState(); //this should be called only once per time step
            m_mixed_module->SetLastStateVariables();
            break;
        }
     
        m_x_mixed = m_mixed_module->Solution();
        m_x_transport = m_transport_module->Solution();
    }
    
    if (!stop_criterion_Q) {
        std::cout << "SFI fail to converge " << std::endl;
        std::cout << "Number of iterations = " << n_iterations << std::endl;
        std::cout << "SFI will continue with : " << std::endl;
        std::cout << "Mixed problem variation = " << error_rel_mixed << std::endl;
        std::cout << "Transport problem variation = " << error_rel_transport << std::endl;
        m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
        m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensitiesLastState(); //this should be called only once per time step
        m_mixed_module->SetLastStateVariables();
        return;
    }
    
    
}

void TMRSSFIAnalysis::PostProcessTimeStep(const int type, const int dim, int step){

    std::cout << "\nTMRSSFIAnalysis Post Process" << std::endl;
    TPZSimpleTimer timer_pp("Timer SFIAnalysis Post Process");
    if (type == 0) {
        m_mixed_module->PostProcessTimeStep(dim, step);
        m_transport_module->PostProcessTimeStep();
    }
    if (type == 1) {
        m_mixed_module->PostProcessTimeStep(dim, step);
    }
    if (type == 2) {
        m_transport_module->PostProcessTimeStep();
    }

    std::cout << "TMRSSFIAnalysis Post Process total time: " << timer_pp.ReturnTimeDouble()/1000. << " seconds" << std::endl;    
}

void TMRSSFIAnalysis::SFIIteration(){
    

    TPZSimpleTimer timer_sfi("Timer SFI Iteration");
    
    m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    m_transport_module->fAlgebraicTransport.fCellsData.UpdateMixedDensity();
    // fAlgebraicDataTransfer.TransferLambdaCoefficients();
    fAlgebraicDataTransfer.TransferSaturation();

    if(shouldSolveDarcy){
        std::cout << "---Running Darcy problem" << std::endl;
        m_mixed_module->RunTimeStep(); // Newton iterations for mixed problem are done here till convergence
        VerifyElementFluxes();
        UpdateAllFluxInterfaces();
        if (m_sim_data->mTNumerics.m_is_linearTrace) {
            shouldSolveDarcy = false;
        }
    }
    fAlgebraicDataTransfer.TransferPressures();
    m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensities();
    
    std::cout << "---Running Transport problem" << std::endl;
    // Solves the transport problem
    m_transport_module->RunTimeStep();
    // m_transport_module->PostProcessTimeStep();
    // fAlgebraicDataTransfer.TransferSaturation();
    // m_mixed_module->SetLastStateSaturation();
    
    std::cout << "SFIIteration time: " << timer_sfi.ReturnTimeDouble()/1000 << " seconds" << std::endl;
}

void TMRSSFIAnalysis::UpdateAllFluxInterfaces(){
    // to compute the transport of saturations in an algebraic way
    fAlgebraicDataTransfer.TransferMixedMeshMultiplyingCoefficients();
    auto& idata = m_transport_module->fAlgebraicTransport.fInterfaceData;
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_id);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracInf);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracSup);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracFrac);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracBound);

    std::set<int> bc_matids;
    for (auto it = m_transport_module->fAlgebraicTransport.fboundaryCMatVal.begin(); it != m_transport_module->fAlgebraicTransport.fboundaryCMatVal.end(); it++) {
        bc_matids.insert(it->first);
    }
    for (auto& bc : bc_matids) {
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(bc);
    }
    
    // m_transport_module->fAlgebraicTransport.VerifyElementFLuxes();
}

void TMRSSFIAnalysis::VerifyElementFluxes(){
    const REAL tol = 1.e-12;
    TPZMultiphysicsCompMesh *mixedmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh()) ;
    TPZCompMesh *cmesh =mixedmesh->MeshVector()[0];
//    std::ofstream file("fuxmesh.txt");
    const TPZFMatrix<STATE> &meshSol = cmesh->Solution();
//    mixedmesh->MeshVector()[0]->Print(file);
    int nels =mixedmesh->MeshVector()[0]->NElements();
    for (int iel =0; iel<nels-1; iel++) {
        TPZCompEl *cel = mixedmesh->MeshVector()[0]->Element(iel);
//            if(cel->Dimension() != cmesh->Dimension()) continue;
        TPZCompElHDiv<pzshape::TPZShapeCube> *hdivel = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(cel);
        TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivelq = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(cel);
        TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollaps = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivelq);
        if (!hdivel && !hdivCollaps) {
            continue;
        }
        if (!cel) continue;
        int ncon = cel->NConnects();
        int nCorners = cel->Reference()->NCornerNodes();
        int nsides1 = cel->Reference()->NSides(1);
        int nsides = cel->Reference()->NSides();
        REAL sumel=0.0;
        for (int icon=0; icon<ncon-1; icon++) {
            TPZConnect &con = cel->Connect(icon);
            int sideOrient =0;
            if(!hdivCollaps){
                 sideOrient = hdivel->GetSideOrient(nCorners+nsides1+icon);
            }
            else{
                if (hdivCollaps && icon < 5) {
                    sideOrient = hdivCollaps->GetSideOrient(nCorners+icon);
                }
                 
            }
            if (hdivCollaps && icon == 4) {
                ncon++;
                continue;
            }
            if (hdivCollaps && icon == 5) {
                sideOrient = hdivCollaps->GetSideOrient(8);
            }
            if (hdivCollaps && icon == 6) {
                sideOrient = hdivCollaps->GetSideOrient(9);
            }
            
           
            int sequence =con.fSequenceNumber;
            int64_t pos = cmesh->Block().Position(sequence);
            if(sequence==-1) continue;
            int blocksize = cmesh->Block().Size(sequence);
            for(int ieq=0; ieq< blocksize; ieq++)
            {
            sumel += sideOrient*meshSol.GetVal(cmesh->Block().Index(sequence,ieq),0);
            }
        }
        if(std::abs(sumel)> tol ){
            std::cout << "------ERROR! Conservation of element index " << cel->Reference()->Index() << " is " << sumel << std::endl;
            DebugStop();
        }
    }
    std::cout << "------All flux elements satisfy conservation up to tolerance " << tol << std::endl;
}

void TMRSSFIAnalysis::TransferToTransportModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }
    
    //     m_mixed_module->m_soltransportTransfer.TransferFromMultiphysics();
    mixed_cmesh->LoadSolutionFromMultiPhysics();
    
    // flux and pressure are transferred to transport module
    int q_b = 0;
    int p_b = 1;
    TPZFMatrix<STATE> & q_dof = mixed_cmesh->MeshVector()[q_b]->Solution();
    TPZFMatrix<STATE> & p_dof = mixed_cmesh->MeshVector()[p_b]->Solution();
    transport_cmesh->MeshVector()[q_b]->LoadSolution(q_dof);
    transport_cmesh->MeshVector()[p_b]->LoadSolution(p_dof);
    
    if (m_sim_data->mTNumerics.m_four_approx_spaces_Q) {
        // average flux and pressure are transferred to transport module
        int qavg_b = 2;
        int pavg_b = 3;
        TPZFMatrix<STATE> & q_dof = mixed_cmesh->MeshVector()[qavg_b]->Solution();
        TPZFMatrix<STATE> & p_dof = mixed_cmesh->MeshVector()[pavg_b]->Solution();
        transport_cmesh->MeshVector()[qavg_b]->LoadSolution(q_dof);
        transport_cmesh->MeshVector()[pavg_b]->LoadSolution(p_dof);
    }
    
    //    m_transport_module->m_soltransportTransfer.TransferToMultiphysics();
    transport_cmesh->LoadSolutionFromMeshes();
    
}

void TMRSSFIAnalysis::TransferToMixedModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }
    
    //    m_transport_module->m_soltransportTransfer.TransferFromMultiphysics();
    transport_cmesh->LoadSolutionFromMultiPhysics();
    
    // Saturations are transferred to mixed module
    int s_b = 2;
    TPZFMatrix<STATE> & s_dof = transport_cmesh->MeshVector()[s_b]->Solution();
    mixed_cmesh->MeshVector()[s_b]->LoadSolution(s_dof);
    
    
    //     m_mixed_module->m_soltransportTransfer.TransferToMultiphysics();
    mixed_cmesh->LoadSolutionFromMeshes();
    
}

void TMRSSFIAnalysis::UpdateMemoryMixedModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    
    if (!mixed_cmesh) {
        DebugStop();
    }
    
    // Updating memory
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, mixed_cmesh, true);
    m_mixed_module->AssembleResidual();
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, mixed_cmesh, false);
    
}

void TMRSSFIAnalysis::UpdateMemoryTransportModule(){
    
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    if (!transport_cmesh) {
        DebugStop();
    }
    
    // Updating memory
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, transport_cmesh, true);
    m_transport_module->AssembleResidual();
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, transport_cmesh, false);
}

void TMRSSFIAnalysis::UpdateMemoryInModules(){
    //    UpdateMemoryMixedModule();
    UpdateMemoryTransportModule();
}

#include "TPZFastCondensedElement.h"
// transfer the permeability and lambda to the element solution for post processing
void TMRSSFIAnalysis::SetMixedMeshElementSolution(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    cmesh->ElementSolution().Redim(nel, 4);
    TPZFMatrix<STATE> &elsol = cmesh->ElementSolution();
    int64_t count = 0;
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(submesh)
        {
            SetMixedMeshElementSolution(submesh);
            continue;
        }
        TPZFastCondensedElement *fast = dynamic_cast<TPZFastCondensedElement *>(cel);
        if(!fast) continue;
        count++;
        elsol(el,0) = fast->GetPermTensor()(0,0);
        elsol(el,1) = fast->GetPermTensor()(1,1);
        elsol(el,2) = fast->GetPermTensor()(2,2);
        elsol(el,3) = fast->GetLambda();
    }
}
void TMRSSFIAnalysis::ReadProperties(std::string name, bool print_table_Q, std::vector<REAL> &Kx, std::vector<REAL> &Ky, std::vector<REAL> &Kz, std::vector<REAL> &Phi){
    
    
    std::ifstream file;
    file.open(name);
    int i=1;
    
    
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        if(l != '/'){
            i=i+1;
//            int val = i%15;
//            if(val ==0){
                double a, b, c, d;
                if(iss >> a >> b >> c >> d) ;
                Kx.push_back(a);
                Ky.push_back(b);
                Kz.push_back(c);
                Phi.push_back(d);
//            };
        };
    };
    
    if(Kx.size() == 0){
        std::cout<<"No data read."<<std::endl;
        DebugStop();
    }
    if(print_table_Q){
        std::cout<<"*************************"<<std::endl;
        std::cout<<"Reading file... ok!"<<std::endl;
        std::cout<<"*************************"<<std::endl;
    }
    file.close();
}
