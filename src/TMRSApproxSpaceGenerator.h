//
//  TMRSApproxSpaceGenerator.h
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#ifndef TMRSApproxSpaceGenerator_h
#define TMRSApproxSpaceGenerator_h

#include "Projection/TPZL2Projection.h"
#include <stdio.h>
#include "TMRSSavable.h"
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMixedDarcyFlow.h"
#include "TMRSDarcyFlowWithMem_impl.h"
#include "TMRSMultiphaseFlow_impl.h"
#include "TMRSMemory.h"
#include "TMRSDataTransfer.h"
#include "TPZTracerFlow.h"
#include "TPZCompMeshTools.h"
#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "TMRSSFIAnalysis.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZHybridizeHDiv.h"

class TPZAlgebraicDataTransfer;

/// Approximation space generator and manager class for iMRS pressure/flow approximation
class TMRSApproxSpaceGenerator : public TMRSSavable {
    
private:
    
	
	/// Adds materials to multiphysics compmesh and sets which ones need memory
	/// @param order order of approximation
	/// @param MatsWithmem set of materials with memory
	/// @param MatsWitOuthmem set of materials without memory
    void AddMultiphysicsMaterialsToCompMesh(const int order, std::set<int> &MatsWithmem, std::set<int> &MatsWitOuthmem);
	
	
	/// Returns the transport materials
	/// @param MatsWithmem set of materials with memory
	/// @param MatsWitOuthmem set of materials without memory
    void GetTransportMaterials(std::set<int> &MatsWithmem, std::set<int> &MatsWitOuthmem);
    
	/// Associates a lagrange multiplier level to the connects of the meshes
    void SetLagrangeMultiplier4Spaces(TPZVec<TPZCompMesh *>& mesh_vec);
	
	
	/// Adds needed materials to the atomic meshes. Mostly creates TPZNullMaterials
	/// @param dim dimension
	/// @param cmesh compmesh to receive materials
	/// @param matids material ids to be added in the atomic cmesh
	/// @param bcmatids boundary material ids to be added to the cmesh
    void AddAtomicMaterials(const int dim, TPZCompMesh* cmesh,
                            std::set<int>& matids,
                            std::set<int>& bcmatids,
                            const bool isInsertBCs = true);
	
	
	/// Method to generate the Hdiv compmesh for problems that have fractures
	/// @param cmesh Hdiv compmesh
	/// @param matids 3D material ids
	/// @param bcids 3D boundary material ids
	/// @param matids_dim2 2D materials ids (fractures)
	/// @param bcids_dim2 2D boundary material ids (fracture bcs)
    void CreateFractureHDivCompMesh(TPZCompMesh* cmesh,
                                    std::set<int>& matids, std::set<int>& bcids,
                                    std::set<int>& matids_dim2, std::set<int>& bcids_dim2);
    
    /// Order overlapping fracture elements such that their normal follows their position in the original fractures
    void OrderOverlappingFractures();
    
    /// Order the overlapping fracture elements such that they correspond to the order of the fracture planes
	/// This will update the neighbour information between 3D elements and between fracture elements
    /// create HDivBound glue elements between the fractures
    void OrderFractures(TPZCompMesh *cmesh, TPZVec<TPZGeoElSide> &fracvec);
    
    
    /// Creates the H(div) spaces of the fracture elements
    void CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh);	
	
	/// Splits the connect in a certain element interface. Mostly used for when there is a 2D fracture in between two 3D elements. Then, the connects are
	/// split with this function, and later (in another methods) an hdivcollapsed is created at that location and the split connects are set as top and bottom of the hdivcollapsed
	/// @param compside compelside that will have its connect split in two
    void SplitConnectsAtInterface(TPZCompElSide& compside);
    
    /// Initialize the integration point information for fracture glue
    void InitializeMemoryFractureGlue();
    
	/// Makes sure all boundary elements of fractures have positive orientation
    void AdjustOrientBoundaryEls(TPZCompMesh* cmesh, std::set<int>& buildmatids);
    
	/// Initial matid for merge meshes. This is used in MergeMeshes() to set which fine elements belong to a certain coarse element.
	/// After MergeMeshes(), this is not needed anymore
    int fInitMatIdForMergeMeshes = -1000;
    
public:
    
	/// Geometric mesh
    TPZGeoMesh * mGeometry;
    
	/// Simulation data: stores all the parameters for the simulation (material, numerical, etc.)
    TMRSDataTransfer mSimData;
    
	/// This is the multiphysics compmesh for the problem
    TPZMultiphysicsCompMesh * mMixedOperator;
    
	/// Compmesh for the transport problem. Not really needed for computations but helpful for creating and storing data
    TPZCompMesh * mTransportOperator;
    
    /// This vector contains the MHM domain index for each geometric element
    TPZVec<int64_t> mSubdomainIndexGel;
    
	/// Class that manages hybridization. Used for hybridizing intersections
	/// NOT USED anymore for BuildMixed4SpacesMultiPhysicsCompMesh only for Mortar mesh because it is overly complicated for our needs
    TPZHybridizeHDiv* mHybridizer;
 
	/// Forcing function to be applied in boundary condition
	/// Jun 2022: Check if this is applied everywhere it should be
    ForcingFunctionBCType<STATE> mForcingFunctionBC;
    
public:
    
	/// Constructor
    TMRSApproxSpaceGenerator();
    
    /// Copy assignment operator
    TMRSApproxSpaceGenerator &operator=(const TMRSApproxSpaceGenerator &other);
    
    /// Copy constructor
	///	Should not be called
    TMRSApproxSpaceGenerator(const TMRSApproxSpaceGenerator &copy) { DebugStop(); }
    
    /// Destructor
    ~TMRSApproxSpaceGenerator();
    
    /// Write object state
    void Write(TPZStream &buf, int withclassid) const;
    
    /// Read object state
    void Read(TPZStream &buf, void *context);
    
    /// Read object state
    virtual int ClassId() const;
    
    void SetGeometry(TPZGeoMesh * geometry);
    
    /// Get method for fInitMatIdForMergeMeshes
    const int& InitMatIdForMergeMeshes() const {return fInitMatIdForMergeMeshes;}
	
	/// Set method for fInitMatIdForMergeMeshes
    int& InitMatIdForMergeMeshes() {return fInitMatIdForMergeMeshes;}
		    
    /// For MHM
    /// Sets the geometry based on a fine and a coarse mesh. It creates a list of subdomains based on that
    /// and fills mSubdomainIndexGel vector that will be used to set the macro domains
    void SetGeometry(TPZGeoMesh * gmeshfine, TPZGeoMesh * gmeshcoarse);
    
	/// Checks if there is a skeleton element between volume elements of different domains
	void CheckMeshIntegrity(TPZGeoMesh* gmesh);
    
    /// Checks if the sideorient between fracture elements is correct
    void CheckSideOrientOfFractureEls();
	
	/// Set method for mSubdomainIndexGel
    void SetSubdomainIndexes(TPZVec<int64_t> &subIndexes){
        mSubdomainIndexGel = subIndexes;
    }
	
	/// Get method for mSubdomainIndexGel (creates copy on return)
    TPZVec<int64_t> GetSubdomainIndexes(){
        return mSubdomainIndexGel;
    }
    
	/// Sets the forcing function for boundary condtions
    void SetForcingFunctionBC(ForcingFunctionBCType<STATE> f){
        mForcingFunctionBC = f;
    }
	
	/// Returns true if there is a forcing function for bcs
    bool HasForcingFunctionBC() const {
        return (bool)mForcingFunctionBC;
    }
	
	/// Get method for forcing function bc
    const ForcingFunctionBCType<STATE> &ForcingFunctionBC() const {
        return mForcingFunctionBC;
    }
	
	/// Get/set method for forcing function bc
    ForcingFunctionBCType<STATE> &ForcingFunctionBC() {
        return mForcingFunctionBC;
    }
    
	/// Returns true if there are fracture intersection in the mesh
    const bool isThereFracIntersection() const;
	
    /// Split the connects of the fluxmesh, create HDivBound elements and pressure elements
	/// NOT USED in 4Space mesh since Junm 2022. Still used in mortar mesh though
    void HybridizeIntersections(TPZVec<TPZCompMesh *>& mesh_vec);
    
	/// Assign a subdomain to the lower level elements
    void IdentifySubdomainForLowdimensionElements(TPZCompMesh *fluxmesh);
    
    /// Adjust the neighbouring information such that the boundary of the fracture elements is the first boundary
    /// Verify if the assigned subdomains are consistent
    void VerifySubdomainIntegrity();
    
    /// Identify the domain indices of the interface elements
    void SetInterfaceDomains(TPZStack<int64_t> &pressureindices,std::pair<int,int> &interfacematids);

	/// Creates interface multiphysics element between L2 pressure 1D element at intersection and flux connect on the border of a fracture
    void CreateIntersectionInterfaceElements(TPZVec<TPZCompMesh *>& meshvec_Hybrid);

	/// Creates interface multiphysics element between L2 pressure 1D element at intersection and flux connect on the border of a fracture
	/// This function does not use mHybridizer deprecated class for generating hibridization
    void CreateIntersectionInterfaceElements();
	
	/// Apply uniform refinements with quantity nelref
    void ApplyUniformRefinement(int nelref);
	
	/// Apply uniform refinement with quantity mSimData.mTGeometry.mnref
    void ApplyUniformRefinement();
    
	/// Helper function to print geomesh to vtk
    void PrintGeometry(std::string name,  bool vtkFile=true, bool textfile=false);
    
	/// Return mGeometry GeoMesh
	TPZGeoMesh * GetGeometry() { return mGeometry; }
    
	/// Sets the data transfer
    void SetDataTransfer(TMRSDataTransfer & DataTransfer);
    
	/// Returns data transfer
    TMRSDataTransfer & GetDataTransfer();
    
	/// Return the hybridizer class. Only used in Mortar simulations
    const TPZHybridizeHDiv* Hybridizer() const {return mHybridizer;}
	
	/// Get/Set the hybridizer class. Only used in Mortar simulations
    TPZHybridizeHDiv* Hybridizer() {return mHybridizer;}
 
    /// Create the HDiv computational mesh
    TPZCompMesh * HdivFluxCmesh(int order);
    
    /// Create hybridized hdiv computational mesh
    TPZCompMesh * HdivFluxMeshHybridized(int order);
    
    /// Create a discontinuous mesh
    TPZCompMesh * DiscontinuousCmesh(int order, char lagrange);
    
    /// Group the connects of the discontinuous mesh such that all connects of a subdomain are identical
    void GroupConnectsBySubdomain(TPZCompMesh *cmesh);
    
    /// Create a discontinuous mesh for transport
    TPZCompMesh * TransportCmesh();
    
	/// Creates auxiliary cmesh for transport. Mostly used to identify interfaces
    void BuildAuxTransportCmesh();
    
    /// create an HDiv mesh for mortar approximation
    TPZCompMesh *HDivMortarFluxCmesh(char mortarlagrange);
    
    /// create a pressure with mortar elements
    TPZCompMesh *PressureMortarCmesh(char firstlagrangepressure, char lagrangepressure, char lagrangemortar);
    
    /// insert the necessary interface elements
    void InsertInterfaceElements();
    
    // build a discontinuous mesh with the order of the elements according to the
    // algebraic transport mesh. The order is always 0
    TPZCompMesh * DiscontinuousCmesh(TPZAlgebraicDataTransfer &Atransfer);
	
	/// Builds the multiphysics cmesh and atomics cmeshes based on booleans previously set
	/// @param order approximation order
    void BuildMixedMultiPhysicsCompMesh(int order);
    
	/// Generates a pressure/flow problem with H(div) space for flow and L2 space for pressure.
	/// May 2022: Only works for problems without fractures and is only tested for 2D domains living in 3D
	/// @param order approximation order
    void BuildMixed2SpacesMultiPhysicsCompMesh(int order);
    
    void BuildMixed4SpacesMultiPhysicsCompMesh(int order);
        
    /// Creates the mesh with 4 spaces and 1 hybridization between the flux elements
    /// First only the 3D elements are hybridized. TODO: hybridize the fracture elements.
    void BuildMixed4SpacesHybridized(int order);
    
    /// build a multiphysics mesh corresponding to a zero order mortar approximation
    void BuildMixed4SpacesMortarMesh();
    
    /// insert wrapper elements necessary for creating the (hybridized) mortar spaces. Only used for building Mortar spaces for now
    void InsertGeoWrappersForMortar();
    
    /// insert wrapper elements necessary for creating the hybridized spaces. Not used for now
    void InsertGeoWrappers();
    
	/// Create GeoEl wrappers for a side when building Mortar spaces
    void GeoWrappersForMortarGelSide(TPZGeoElSide &gelside, std::set<int> bcids);
	
	/// Returns the subdomain of the neighbour element through side gelside
    int FindNeighSubDomain(TPZGeoElSide &gelside);

	/// Used in TPZMHMixedMeshWithTransportControl::BuildComputationalMesh(). Can we delete it? (Jun 2022)
    void BuildTransportMultiPhysicsCompMesh();
    
	/// Builds transport computational mesh with 2 spaces. DEPRECATED Jun 2022
    void BuildTransport2SpacesMultiPhysicsCompMesh();
	
	/// Builds transport computational mesh with 4 spaces. DEPRECATED Jun 2022
    void BuildTransport4SpacesMultiPhysicsCompMesh();
    
    /// Return the material ids and boundary condition material ids
    void GetMaterialIds(int dim, std::set<int> &matids, std::set<int> &bcmatids);
    
	/// Returns the multiphysics compmesh
    TPZMultiphysicsCompMesh * GetMixedOperator() { return mMixedOperator; }
    
	/// Returns the transport compmesh
    TPZCompMesh * GetTransportOperator() { return mTransportOperator; }
    
    /// Linking the memory between the operators. DEPRECATED Jun 2022
    void LinkMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZCompMesh * TransportOperator);
    
	/// Adjust integration rules. DEPRECATED Jun 2022
    static void AdjustMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
	/// Set same memory for both meshes. DEPRECATED Jun 2022
    static void UnifyMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
	/// Fills memory with data such as permeability, porosity, etc.
    static void FillMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
	/// Fills memory with data such as permeability, porosity, etc.
    static void FillMaterialMemoryDarcy(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZAlgebraicTransport *algebraicTranspor);
    
	/// Sets to update (true or false depending on update_memory_Q) the memory of material with matid = material_id
    static void SetUpdateMaterialMemory(int material_id, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q = true);
    
	/// Sets to update (true or false depending on update_memory_Q) the memory of all 3D materials on the mesh
    static void SetUpdateMemory(int dimension, TMRSDataTransfer & sim_data, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q = true);

	/// Initialize all fracture properties (TMRSMemory)
    void InitializeFracProperties(TPZMultiphysicsCompMesh * MixedOperator);
   
    /// Returns in vector neihside all the neighbour elements of gelside that have matid equal to any matid in VolMatIds
    void findNeighElementbyMatId(TPZGeoElSide &gelside, std::vector<TPZGeoElSide > &neihside, std::set<int> VolMatIds);
		
	/// Creates interface elements between elements. Used by transport model to calculate flux between elements
    void CreateInterfaces(TPZCompMesh *cmesh);

	/// Creates interface elements between 3D elements and 3D elements and boundaries. Used by transport model to calculate flux between elements
	void CreateElementInterfaces(TPZGeoEl *gel);

	/// Creates interface elements between 2D elements and 2D elements and boundaries. Used by transport model to calculate flux between elements
    void CreateFracInterfaces(TPZGeoEl *gel);
		
	/// Wrapper function to create a interface element based on two geoelsides
    void CreateInterfaceElements(TPZGeoElSide &gelside, TPZGeoElSide &gelneig, int matid);
	
	/// Creates a compel for the transport model
    void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);
	
	/// Groups the elements belonging to a same subdomain and put them in subcompmeshes
    void HideTheElements(TPZCompMesh *cmesh);
    
	/// Takes a fine and coarse mesh and
	/// 1) Creates the parenthood data structure between macro and micro elements -- stores in mSubdomainIndexGel
	/// 2) Creates skeleton elements between macro domains
    void MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh);
	
	/// Returns true if there are any fractures in the model
    const bool isFracSim() const {return mSimData.mTGeometry.mDomainFracNameAndMatId.size();}
    
	/// Returns true if the matid is a fracture matid
    bool IsFracMatId(int matiD);
    
    // update dp in MemoryFractureGlue
    void UpdateMemoryFractureGlue(TPZCompMesh *cmesh);
    
	/// Returns true if the matid is a fracture bc matid
    bool IsFracBCMatId(int matiD) { return IsFracMatId(matiD-1); }
	
	/// Returns true if the matid is a fracture intersection matid
    bool IsFracIntersectMatId(int matiD) { return IsFracMatId(matiD-2); }
 
};

#endif /* TMRSApproxSpaceGenerator_h */
