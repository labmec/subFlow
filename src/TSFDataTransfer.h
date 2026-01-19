//
//  Created by Giovane Avancini on 08/01/26.
//  Just copied from iMRS and needs to be adapted to subFlow
//  I'm keeping the original structure for now just to have a working version
//  PLEASE ADAPT THIS FILE
//

#pragma once

#include "TPZAlgebraicTransport.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
class TPZSubCompMesh;
class TPZFastCondensedElement;
class TPZAlgebraicTransport;

class TSFDataTransfer {

public:
  // Spatial reservoir properties and initial saturation
  std::function<REAL(const TPZVec<REAL> &)> fkx;
  std::function<REAL(const TPZVec<REAL> &)> fky;
  std::function<REAL(const TPZVec<REAL> &)> fkz;
  std::function<REAL(const TPZVec<REAL> &)> fphi;
  std::function<REAL(const TPZVec<REAL> &)> fs0;

  std::function<std::vector<REAL>(const TPZVec<REAL> &)> fkappa_phi;

  TPZMultiphysicsCompMesh *fFluxMesh;

  TPZCompMesh *fTransportMesh;

  struct TInterfaceWithVolume {
    // geometric element index of the interface element
    int64_t fInterface_gelindex;
    // computational element index of the interface element
    int64_t fInterface_celindex;
    // left right geometric element index of the associated volume elements
    std::pair<TPZGeoElSideIndex, TPZGeoElSideIndex> fLeftRightGelSideIndex;
    // left right volume index in the AlgebraicTransport data structure
    std::pair<int64_t, int64_t> fLeftRightVolIndex;

    std::pair<int64_t, int64_t> fLeftRightGelIndex;

    TInterfaceWithVolume() : fInterface_gelindex(-1), fInterface_celindex(-1), fLeftRightGelSideIndex(), fLeftRightVolIndex(-1, -1), fLeftRightGelIndex(-1, -1) {
    }
  };

  struct TFromMixedToTransport {
    /// material id associated with this transfer structure
    int fMatid;
    /// integer of the flux data structure
    int flux_sequence;
    /// gather vector with respect to the solution of the mixed mesh
    std::vector<int64_t> fGather;
    /// scatter vector with respect to the transport interface
    std::vector<int64_t> fScatter;
    /// pointer to the matrix data structure
    TPZFMatrix<REAL> *fFrom;
    /// pointer to the vector data structure
    std::vector<REAL> *fTarget;
    /// pointer to the computational mesh
    TPZCompMesh *fMixedMesh;

    TFromMixedToTransport();

    TFromMixedToTransport(const TFromMixedToTransport &copy);

    TFromMixedToTransport &operator=(const TFromMixedToTransport &copy);

    ~TFromMixedToTransport() {
    }

    void Print(std::ostream &out);

    void TransferSolution();
  };

  struct TransportToMixedCorrespondence {

    TPZCompMesh *fMixedMesh;

    // Algebraic cell index
    std::vector<int64_t> fAlgebraicTransportCellIndex;

    // The equation number in the original transport mesh
    std::vector<int64_t> fEqNum;

    TPZAlgebraicTransport *fTransport;
    std::vector<TPZCondensedCompEl *> fMixedCell;

    TransportToMixedCorrespondence() : fMixedMesh(0), fAlgebraicTransportCellIndex(0), fEqNum(0), fTransport(0) {}

    TransportToMixedCorrespondence(const TransportToMixedCorrespondence &cp) {
      fMixedMesh = cp.fMixedMesh;
      fAlgebraicTransportCellIndex = cp.fAlgebraicTransportCellIndex;
      fTransport = cp.fTransport;
      fEqNum = cp.fEqNum;

      fMixedCell = cp.fMixedCell;
    }

    TransportToMixedCorrespondence &operator=(const TransportToMixedCorrespondence &cp) {
      fMixedMesh = cp.fMixedMesh;
      fAlgebraicTransportCellIndex = cp.fAlgebraicTransportCellIndex;
      fTransport = cp.fTransport;
      fMixedCell = cp.fMixedCell;
      fEqNum = cp.fEqNum;
      return *this;
    }

    void Print(std::ostream &out = std::cout);
  };

  // Interface data structure, one material at a time
  std::map<int, TPZVec<TInterfaceWithVolume>> fInterfaceGelIndexes;

  // The index of the computational volume elements in the transport mesh identified by material id
  // for each algebraic cell
  TPZVec<int64_t> fVolumeElements;

  // The index of the computational volume elements in the transport mesh identified by material id
  TPZVec<int64_t> fConnectsByInterfaceMatID;

  /// Interface element associated with each volumetric geometric element/side when applicable
  // the size of this data structure is Number Geometric Elements x 6
  // the interface index is related to the algebraic data structure
  TPZFMatrix<int64_t> fInterfaceByGeom;

  /// List of transfer information from the mixed mesh to the transport mesh
  std::map<int, std::list<TFromMixedToTransport>> fTransferMixedToTransport;

  /// List of correspondence of transport volumes and mixed elements for each mesh
  std::list<TransportToMixedCorrespondence> fTransportMixedCorrespondence;

public:
  /// Default constructor
  TSFDataTransfer();

  /// Copy constructor
  TSFDataTransfer(const TSFDataTransfer &other);

  /// Assignement constructor
  const TSFDataTransfer &operator=(const TSFDataTransfer &other);

  /// Default desconstructor
  ~TSFDataTransfer();

  void SetMeshes(TPZMultiphysicsCompMesh &fluxmesh, TPZCompMesh &transportmesh) {
    fFluxMesh = &fluxmesh;
    fTransportMesh = &transportmesh;
  }

  // compute the data transfer data structures between the fluxmesh and transport class
  void BuildTransportDataStructure(TPZAlgebraicTransport &transport);

  // Initialize the datastructures of the transport object
  void InitializeTransportDataStructure(TPZAlgebraicTransport &transport);

  // Build the data structure which defines the correspondence between
  // algebraic transport cells and indexes of mixed fast condensed elements
  // @param Volume_indexes : for each geometric element the algebraic volume index, if applicable, otherwise -10
  void BuildTransportToMixedCorrespondenceDatastructure(TPZCompMesh *fluxmesh, TPZVec<int64_t> &Volume_indexes);

  // Identify the geometric elements corresponding to interface elements. Order them as
  // a function of the number of corner nodes
  void IdentifyInterfaceGeometricElements();

  // Identify volume information to the interface data structure (TInterfaceWithVolume)
  // will loop over the interface elements, identify left and right elements and
  // initialize the fInterfaceByGeom data structure
  void IdentifyVolumeGeometricElements();
  void IdentifyVolumeGeometricElements2();

  // print the datastructure
  void Print(std::ostream &out = std::cout);

  // identify material of a face which is connected to a given geometric element
  std::vector<int> IdentifyMaterial(TPZGeoElSideIndex gelsideindex, int64_t faceindex);

  // identify material of a face which is connected to a given geometric element
  std::vector<int> IdentifyMaterial(const TPZGeoElSide &gelside, int delta = 0);

  // find the neighbouring interface element
  TPZGeoElSide IdentifyInterfaceElement(const TPZGeoElSide &gelside);

  /// extract the list of computational element from a substructured computational mesh
  // this method also searches for elements in element groups and condensed elements
  // each computational element has an associated geometric element
  static void GetElementAndSubmeshPointers(TPZCompMesh &mixedmesh, std::list<TPZCompEl *> &elpointers, std::list<TPZSubCompMesh *> &submeshes);

  /// build the data structure from mixed to transport
  // this method finds the correspondence between the connects of the volume
  // elements and the interfaces of the transport mesh
  // this method fills in the fTransferMixedToTransport variable
  void BuildMixedToTransportDataStructures(TPZCompMesh *fluxmesh);

  // Initialize the pointers to the transport data structure in the
  // list fTransferMixedToTransport
  void InitializeVectorPointersMixedToTransport(TPZAlgebraicTransport &transport);

  // Initialize the pointers from the transport data structure in the list TransportToMixedCorrespondence
  void InitializeVectorPointersTranportToMixed(TPZAlgebraicTransport &transport);

  // transfer the solution from the mixed mesh fluxes to the interfaces
  void TransferMixedMeshMultiplyingCoefficients();

  void TransferPressures();

  // transfer the permeability multiplier from the transport mesh to the mixed mesh elements
  void TransferLambdaCoefficients();
  void TransferPermeabiliyTensor();
  void TransferSaturation();

  // verify the correspondence of the mixed elements and the algebraic cells
  void CheckDataTransferTransportToMixed();
  void TakeOrientationAndLowerIndex(TPZCompElSide &celSide, int &orientation, int &lowerIndex, int matId);
  void TakeOrientationAndLowerIndexDimVolDimFrac(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matid);
  void TakeOrientationAndLowerIndexDimDim(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matid);
  void TakeOrientationAndLowerIndexFracFrac(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matid);
  TPZMultiphysicsElement *findMultiphysics(TPZElementGroup *group);

  std::pair<int, int> FindMortar(TPZGeoElSide &gelside);
  std::pair<int, int> FindMortar(TPZGeoElSide &gelside, int targetId);
  void TestSideOrient(TPZCompMesh *MultFlux);
};