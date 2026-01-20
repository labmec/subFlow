//
//  Created by Giovane Avancini on 08/01/26.
//  Just copied from iMRS and needs to be adapted to subFlow
//  I'm keeping the original structure for now just to have a working version,
//  with minor modifications
//  PLEASE ADAPT THIS FILE
//

#pragma once

#include "TPZMultiphysicsCompMesh.h"
#include "TSFAlgebraicTransport.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
class TPZSubCompMesh;
class TPZFastCondensedElement;
class TSFAlgebraicTransport;

class TSFDataTransfer {

public:
  // Spatial reservoir properties and initial saturation
  std::function<REAL(const TPZVec<REAL> &)> fkx;
  std::function<REAL(const TPZVec<REAL> &)> fky;
  std::function<REAL(const TPZVec<REAL> &)> fkz;
  std::function<REAL(const TPZVec<REAL> &)> fphi;
  std::function<REAL(const TPZVec<REAL> &)> fs0;

  std::function<std::vector<REAL>(const TPZVec<REAL> &)> fkappa_phi;

  TPZMultiphysicsCompMesh *fDarcyCmesh;

  TPZCompMesh *fTransportCmesh;

  struct TInterfaceWithVolume {
    // geometric element index of the interface element
    int64_t fGelIndex;
    // computational element index of the interface element
    int64_t fCelIndex;
    // left right geometric element index of the associated volume elements
    std::pair<TPZGeoElSideIndex, TPZGeoElSideIndex> fLeftRightGelSideIndex;
    // left right volume index in the AlgebraicTransport data structure
    std::pair<int64_t, int64_t> fLeftRightVolIndex;
    // left right geometric element index of the associated volume elements
    std::pair<int64_t, int64_t> fLeftRightGelIndex;

    TInterfaceWithVolume() : fGelIndex(-1), fCelIndex(-1), fLeftRightGelSideIndex(), fLeftRightVolIndex(-1, -1), fLeftRightGelIndex(-1, -1) {}
  };

  struct TFromDarcyToTransport {
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
    TPZCompMesh *fDarcyCmesh;

    TFromDarcyToTransport();

    TFromDarcyToTransport(const TFromDarcyToTransport &copy);

    TFromDarcyToTransport &operator=(const TFromDarcyToTransport &copy);

    ~TFromDarcyToTransport() {
    }

    void Print(std::ostream &out);

    void TransferSolution();
  };

  struct TTransferCorrespondence {

    TPZCompMesh *fDarcyCmesh;

    // Algebraic cell index
    std::vector<int64_t> fAlgebraicTransportCellIndex;

    // The equation number in the original transport mesh
    std::vector<int64_t> fEqNum;

    // The Algebraic transport pointer
    TSFAlgebraicTransport *fTransport;

    // Darcy condensed computational elements associated with each algebraic transport cell
    std::vector<TPZCondensedCompEl *> fDarcyCels;

    TTransferCorrespondence() : fDarcyCmesh(0), fAlgebraicTransportCellIndex(0), fEqNum(0), fTransport(0) {}

    TTransferCorrespondence(const TTransferCorrespondence &cp) {
      fDarcyCmesh = cp.fDarcyCmesh;
      fAlgebraicTransportCellIndex = cp.fAlgebraicTransportCellIndex;
      fTransport = cp.fTransport;
      fEqNum = cp.fEqNum;

      fDarcyCels = cp.fDarcyCels;
    }

    TTransferCorrespondence &operator=(const TTransferCorrespondence &cp) {
      fDarcyCmesh = cp.fDarcyCmesh;
      fAlgebraicTransportCellIndex = cp.fAlgebraicTransportCellIndex;
      fTransport = cp.fTransport;
      fDarcyCels = cp.fDarcyCels;
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

  /// List of transfer information from the Darcy mesh to the transport mesh
  std::map<int, std::list<TFromDarcyToTransport>> fTransferDarcyToTransport;

  /// List of correspondence of transport cells and darcy elements for each mesh
  std::list<TTransferCorrespondence> fCorrespondence;

public:
  /// Default constructor
  TSFDataTransfer();

  /// Copy constructor
  TSFDataTransfer(const TSFDataTransfer &other);

  /// Assignement constructor
  const TSFDataTransfer &operator=(const TSFDataTransfer &other);

  /// Default desconstructor
  ~TSFDataTransfer();

  void SetMeshes(TPZMultiphysicsCompMesh *darcy_cmesh, TPZCompMesh *transport_cmesh) {
    fDarcyCmesh = darcy_cmesh;
    fTransportCmesh = transport_cmesh;
  }

  // compute the data transfer data structures between the darcy and transport problems
  void Initialize();

  // Initialize the datastructures of the transport object
  void InitializeAlgebraicTransport(TSFAlgebraicTransport &transport);

  // Build the data structure which defines the correspondence between
  // algebraic transport cells and indexes of darcy fast condensed elements
  // @param transport_cell_ids : for each geometric element the algebraic volume index, if applicable, otherwise -10
  void BuildTransferCorrespondenceDatastructure(TPZCompMesh *darcy_cmesh, TPZVec<int64_t> &transport_cell_ids);

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

  /// extract the list of condensed compels from a computational mesh
  // this method searches only for condensed elements
  // each computational element has an associated geometric element
  static void GetCondensedElements(TPZCompMesh *cmesh, TPZStack<TPZCompEl *> &cels);

  /// build the data structure from darcy to transport
  // this method finds the correspondence between the connects of the volume
  // elements and the interfaces of the transport mesh
  // this method fills in the fTransferDarcyToTransport variable
  void BuildDarcyToTransportDataStructures(TPZCompMesh *darcy_cmesh);
  // Initialize the pointers to the transport data structure in the
  // list fTransferDarcyToTransport
  void InitializeVectorPointersDarcyToTransport(TSFAlgebraicTransport &transport);

  // Initialize the pointers from the transport data structure in the list TTransferCorrespondence
  void InitializeVectorPointersTransportToDarcy(TSFAlgebraicTransport &transport);

  // transfer the solution from the darcy mesh fluxes to the interfaces
  void TransferDarcyMeshMultiplyingCoefficients();

  void TransferPressures();

  // transfer the permeability multiplier from the transport mesh to the darcy mesh elements
  void TransferLambdaCoefficients();
  void TransferPermeabiliyTensor();
  void TransferSaturation();

  // verify the correspondence of the darcy elements and the algebraic cells
  void CheckDataTransferTransportToDarcy();
  void TakeOrientationAndLowerIndex(TPZCompElSide &celSide, int &orientation, int &lowerIndex, int matId);
  void TakeOrientationAndLowerIndexDimVolDimFrac(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matid);
  void TakeOrientationAndLowerIndexDimDim(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matid);
  void TakeOrientationAndLowerIndexFracFrac(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matid);
  TPZMultiphysicsElement *findMultiphysics(TPZElementGroup *group);
  void TestSideOrient(TPZCompMesh *MultFlux);
};