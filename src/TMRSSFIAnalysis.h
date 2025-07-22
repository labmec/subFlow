//
//  TMRSSFIAnalysis.h
//
//  Created by Omar Durán on 10/15/19.
//

#ifndef TMRSSFIAnalysis_h
#define TMRSSFIAnalysis_h

#include <stdio.h>
#include "TPZMultiphysicsCompMesh.h"
#include "TMRSDataTransfer.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TMRSApproxSpaceGenerator.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZAlgebraicDataTransfer.h"


/**
 * @class TMRSSFIAnalysis Manager for solving pressure/flux problem and tranport problem. Also takes care of exchanging data between them
 * @brief SFI stands for Sequential Full Implicit
 */
class TMRSSFIAnalysis {
    
private:
    
    /// Data transfer object. Holds all the problem parameters
    TMRSDataTransfer * m_sim_data;

    /// Hold the pressure/flux solution at every step of the non linear iterations. If it is linear, holds the final solution
    TPZFMatrix<STATE> m_x_mixed;
    
    /// Hold the saturation solution at every step of the non linear iterations. If it is linear, holds the final solution
    TPZFMatrix<STATE> m_x_transport;
    
    /// Boolean that indicates if the problem is linear. If so, does not need to run newton iterations to converge
    bool shouldSolveDarcy = true;
    
public:
    
    /// Number of iterations
    int m_k_iteration = 0;
    
    /// Mixed analysis for managing pressure/flux problem
    TMRSMixedAnalysis * m_mixed_module;
    
    /// Transport analysis for managing the transport problem
    TMRSTransportAnalysis * m_transport_module;
    
    /// Object that takes care of solving the transport problem in a purely algebraic way and transfering information between problem also algebraically
    TPZAlgebraicDataTransfer fAlgebraicDataTransfer;
    
    /// Default constructor
    TMRSSFIAnalysis();
    
    /// Default destructor
    ~TMRSSFIAnalysis();
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed,
                    TPZCompMesh * cmesh_transport,
                    const RenumType& renumtype = RenumType::EDefault);
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed,
                    TPZMultiphysicsCompMesh * cmesh_transport,
                    std::function<REAL(const TPZVec<REAL> & )> & kx,
                    std::function<REAL(const TPZVec<REAL> & )> & ky,
                    std::function<REAL(const TPZVec<REAL> & )> & kz,
                    std::function<REAL(const TPZVec<REAL> & )> & phi,
                    std::function<REAL(const TPZVec<REAL> & )> & s0,
                    const RenumType& renumtype = RenumType::EDefault);
    
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed,
                    TPZMultiphysicsCompMesh * cmesh_transport,
                    std::function<std::vector<REAL>(const TPZVec<REAL> & )> & kappa_phi,
                    std::function<REAL(const TPZVec<REAL> & )> & s0,
                    const RenumType& renumtype = RenumType::EDefault);
    void BuildAlgebraicDataStructure();
    
    // Deprecated methods since Apr 2022
    void FillMaterialMemoryDarcy(int material_id);
    void FillProperties();
    void FillProperties(std::string fileprops, TPZAlgebraicTransport *algebraicTransport);
    void FillProperties(TPZAlgebraicTransport *algebraicTransport);
    void FillProperties(TPZAlgebraicTransport *algebraicTransport, std::vector<REAL> kappa_phi);
    
    static  void ReadProperties(std::string name, bool print_table_Q, std::vector<REAL> &Kx, std::vector<REAL> &Ky, std::vector<REAL> &Kz, std::vector<REAL> &Phi);
    /// Configurates iternal members
    void Configure(int n_threads, bool UsePardiso_Q, bool UsePZ=false);
    
    /// Set data transfer object
    void SetDataTransferAndBuildAlgDatStruct(TMRSDataTransfer * sim_data);
    
    /// Get data transfer object
    TMRSDataTransfer * GetDataTransfer();
    
    /// Get the number of iterations
    int GetNumberOfIterations();
    
    /// Run a time step
    void RunTimeStep();
    
    /// Render a vtk file with requested variables for a time step
    void PostProcessTimeStep(const int type, const int dim = 3, int step = -1);
    
    /// Perform a SFI iteration
    void SFIIteration();
    void UpdateAllFluxInterfaces();
    void VerifyElementFluxes();
    
    void TransferToTransportModule();
    
    void TransferToMixedModule();
    
    void UpdateMemoryMixedModule();
    
    void UpdateMemoryTransportModule();
    
    void UpdateMemoryInModules();
    void PostProcessResevoirProp();
    
    // transfer the permeability and lambda to the element solution for post processing
    void SetMixedMeshElementSolution(TPZCompMesh *cmesh);
    
};

#endif /* TMRSSFIAnalysis_h */
