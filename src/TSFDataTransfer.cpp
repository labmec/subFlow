//
//  Created by Giovane Avancini on 08/01/26.
//  Just copied from iMRS and needs to be adapted to subFlow
//  I'm keeping the original structure for now just to have a working version
//  PLEASE ADAPT THIS FILE
//

#include "TSFDataTransfer.h"
#include "TPZCompElHDivCollapsed.h"
#include "TPZFastCondensedElement.h"
#include "TPZInterfaceEl.h"
#include "TPZLagrangeMultiplier.h" // for TPZLagrangeMultiplier
#include "TPZVTKGeoMesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzmultiphysicselement.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapepiram.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"
#include "pzsubcmesh.h"
/// Default constructor
TSFDataTransfer::TSFDataTransfer() : fFluxMesh(0), fTransportMesh(0) {
}

/// Copy constructor
TSFDataTransfer::TSFDataTransfer(const TSFDataTransfer &other) : fFluxMesh(other.fFluxMesh), fTransportMesh(other.fTransportMesh), fInterfaceGelIndexes(other.fInterfaceGelIndexes), fVolumeElements(other.fVolumeElements), fConnectsByInterfaceMatID(other.fVolumeElements) {
}

/// Assignement constructor
const TSFDataTransfer &TSFDataTransfer::operator=(const TSFDataTransfer &other) {
  fFluxMesh = other.fFluxMesh;
  fTransportMesh = other.fTransportMesh;
  fInterfaceGelIndexes = other.fInterfaceGelIndexes;
  fVolumeElements = other.fVolumeElements;
  fConnectsByInterfaceMatID = other.fConnectsByInterfaceMatID;
  return *this;
}

/// Default desconstructor
TSFDataTransfer::~TSFDataTransfer() {
}

// compute the data transfer data structures between the fluxmesh and transport class
void TSFDataTransfer::BuildTransportDataStructure(TPZAlgebraicTransport &transport) {
  IdentifyInterfaceGeometricElements();
  IdentifyVolumeGeometricElements2();
  BuildMixedToTransportDataStructures(fFluxMesh);
  TPZVec<int64_t> Volume_Index;
  BuildTransportToMixedCorrespondenceDatastructure(fFluxMesh, Volume_Index);
  InitializeTransportDataStructure(transport);
  InitializeVectorPointersMixedToTransport(transport);
  InitializeVectorPointersTranportToMixed(transport);
  CheckDataTransferTransportToMixed();
}

// Identify the geometric elements corresponding to interface elements. Order them as
// a function of the number of corner nodes
void TSFDataTransfer::IdentifyInterfaceGeometricElements() {
  // look for the geometric elements corresponding to interface elements
  // order them as a function of the number of corner nodes
  TPZGeoMesh *gmesh = fTransportMesh->Reference();
  std::pair<int, int64_t> defpair(100, -1);
  fTransportMesh->LoadReferences();
  int64_t neltr = fTransportMesh->NElements();
  TPZVec<std::pair<int, int64_t>> interfaces(neltr, defpair);
  int64_t num_interfaces = 0;
  for (int64_t el = 0; el < neltr; el++) {
    TPZCompEl *cel = fTransportMesh->Element(el);
    TPZMultiphysicsInterfaceElement *interf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
    TPZInterfaceElement *interf2 = dynamic_cast<TPZInterfaceElement *>(cel);
    if (!interf && !interf2)
      continue;

    TPZGeoEl *gel = NULL;
    if (interf)
      gel = interf->Reference();
    else
      gel = interf2->Reference();

    int ncorner = gel->NCornerNodes();
    if (ncorner > 4)
      DebugStop(); // an interface element should not have more than 4 corner nodes

    interfaces[num_interfaces++] = std::pair<int, int64_t>(ncorner, el);
  }
  // this vector has all interface elements
  interfaces.Resize(num_interfaces);
  std::sort(&interfaces[0], (&interfaces[0] + num_interfaces));

  // number of interface elements by matid and number of corner nodes
  std::map<int, std::map<int, int64_t>> numinterfaces_matid;
  // number of interface elements by matid
  std::map<int, int64_t> total_internal_interfaces_matid;
  for (int64_t el = 0; el < num_interfaces; el++) {
    if (interfaces[el].first > 4)
      DebugStop(); // an interface element should not have more than 4 corner nodes
    int64_t celindex = interfaces[el].second;
    TPZCompEl *cel = fTransportMesh->Element(celindex);
    int64_t gelindex = cel->Reference()->Index();
    int matid = gmesh->Element(gelindex)->MaterialId();
    total_internal_interfaces_matid[matid]++;
    numinterfaces_matid[matid][interfaces[el].first]++;
  }

  // then the number of interfaces elements will be the size of this data structure
  for (auto it = total_internal_interfaces_matid.begin(); it != total_internal_interfaces_matid.end(); it++) {
    fInterfaceGelIndexes[it->first].Resize(it->second);
  }
  // compute an index within the vector of interfaces
  std::map<int, std::map<int, int64_t>> count_interfaces;
  for (auto matit = numinterfaces_matid.begin(); matit != numinterfaces_matid.end(); matit++) {
    int64_t prev_count = 0;
    int matid = matit->first;
    for (auto it = matit->second.rbegin(); it != matit->second.rend(); it++) {
      int ncorner_nodes = it->first;
      count_interfaces[matid][ncorner_nodes] = prev_count;
      int64_t numfaces = it->second;
      prev_count += numfaces;
    }
  }

  for (int64_t el = 0; el < num_interfaces; el++) {
    if (interfaces[el].first > 4) DebugStop();
    int64_t celindex = interfaces[el].second;
    TPZCompEl *cel = fTransportMesh->Element(celindex);
    TPZGeoEl *gel = cel->Reference();
    int64_t gelindex = gel->Index();
    int matid = gmesh->Element(gelindex)->MaterialId();
    int ncornernodes = gel->NCornerNodes();
    TInterfaceWithVolume &intface = fInterfaceGelIndexes[matid][count_interfaces[matid][ncornernodes]];
    intface.fInterface_gelindex = gelindex;
    intface.fInterface_celindex = celindex;
    count_interfaces[matid][ncornernodes]++;
  }

#ifdef PZDEBUG
  {
    for (auto &it : fInterfaceGelIndexes) {
      TPZVec<TInterfaceWithVolume> &vec = it.second;
      int64_t nel = vec.size();
      for (int64_t el = 0; el < nel; el++) {
        if (vec[el].fInterface_gelindex == -1)
          DebugStop();
      }
    }
  }
#endif
}

int SideLowerIndex(TPZGeoEl *gel, int side) {
  int dim = gel->Dimension();
  if (dim == 1) return side;
  if (dim == 2) return side - gel->NCornerNodes();
  return side - gel->NCornerNodes() - gel->NSides(1);
}
int SideOriginalIndex(TPZGeoEl *gel, int side) {
  int dim = gel->Dimension();
  if (dim == 1) return side;
  if (dim == 2) return side + gel->NCornerNodes();
  return side + gel->NCornerNodes() + gel->NSides(1);
}

// Identify volume information to the interface data structure (TInterfaceWithVolume)
void TSFDataTransfer::IdentifyVolumeGeometricElements2() {
  TPZGeoMesh *gmesh = fTransportMesh->Reference();

  int64_t nel = gmesh->NElements();
  fInterfaceByGeom.Redim(nel, 6);
  for (int64_t el = 0; el < nel; el++) {
    for (int f = 0; f < 6; f++) {
      fInterfaceByGeom(el, f) = -1;
    }
  }
  int globcount = 0;
  TPZVec<int64_t> geometricvolume(nel, 0);
  int64_t volumecount = 0;
  TPZVec<int64_t> VolumeElementIndex(fTransportMesh->NElements(), -1);
  fFluxMesh->LoadReferences();
  fConnectsByInterfaceMatID.resize(fFluxMesh->NConnects());
  for (int ipos = 0; ipos < fFluxMesh->NConnects(); ipos++) {
    fConnectsByInterfaceMatID[ipos] = -10000;
  }

  for (auto it = fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++) {
    TPZVec<TInterfaceWithVolume> &facevec = it->second;
    int64_t nfaces = facevec.size();
    for (int64_t iface = 0; iface < nfaces; iface++) {
      int64_t celindex = facevec[iface].fInterface_celindex;
      TPZCompEl *cel = fTransportMesh->Element(celindex);
      int meshdim = fTransportMesh->Dimension();
      TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
      if (!intface) {
        DebugStop();
      }

      TPZCompElSide leftside = intface->LeftElementSide();
      TPZCompElSide rightside = intface->RightElementSide();
      TPZCompEl *left = leftside.Element();
      TPZGeoEl *leftgel = left->Reference();
      TPZGeoElSideIndex leftgelside(leftgel->Index(), leftside.Side());
      int dimL = leftgel->Dimension();

      TPZCompEl *right = rightside.Element();
      if (right->NConnects() > 1) DebugStop();
      TPZGeoEl *rightgel = right->Reference();
      TPZGeoElSideIndex rightgelside(rightgel->Index(), rightside.Side());
      int64_t rightindex = right->Index();
      int dimR = rightgel->Dimension();
      if (it->first == 101) {
        int ok = 0;
      }
      int orientationL = 0, lowerIndexL = 0, orientationR = 0, lowerIndexR = 0;
      if (dimL > dimR) {
        TakeOrientationAndLowerIndexDimVolDimFrac(leftside, rightside, orientationL, lowerIndexL, orientationR, lowerIndexR, it->first);
      } else if (dimL < dimR) {
        TakeOrientationAndLowerIndexDimVolDimFrac(rightside, leftside, orientationR, lowerIndexR, orientationL, lowerIndexL, it->first);
      } else if ((dimL == dimR) && (dimR == meshdim)) {
        TakeOrientationAndLowerIndexDimDim(rightside, leftside, orientationR, lowerIndexR, orientationL, lowerIndexL, it->first);
      } else if ((dimL == dimR) && dimR == meshdim - 1) {
        TakeOrientationAndLowerIndexFracFrac(rightside, leftside, orientationR, lowerIndexR, orientationL, lowerIndexL, it->first);
        //                DebugStop();
      } else {
        DebugStop();
      }
      fInterfaceByGeom(leftgel->Index(), lowerIndexL) = iface;
      int64_t leftindex = left->Index();
      if (VolumeElementIndex[leftindex] == -1 && (left->NConnects() == 1)) {
        VolumeElementIndex[leftindex] = volumecount;
        volumecount++;
      }
      fInterfaceByGeom(rightgel->Index(), lowerIndexR) = iface;
      if (VolumeElementIndex[rightindex] == -1 && right->NConnects() == 1) {
        VolumeElementIndex[rightindex] = volumecount;
        volumecount++;
      }
      if (orientationL == 1) {
        if (orientationR == -1 || orientationR == 0) {
          facevec[iface].fLeftRightGelSideIndex = {leftgelside, rightgelside};
          facevec[iface].fLeftRightVolIndex = {VolumeElementIndex[leftindex], VolumeElementIndex[rightindex]};
          facevec[iface].fLeftRightGelIndex = {leftgel->Index(), rightgel->Index()};
        } else {
          std::cout << "Warning" << std::endl;
          facevec[iface].fLeftRightGelSideIndex = {rightgelside, leftgelside};
          facevec[iface].fLeftRightVolIndex = {VolumeElementIndex[rightindex], VolumeElementIndex[leftindex]};
          facevec[iface].fLeftRightGelIndex = {rightgel->Index(), leftgel->Index()};
          //                    DebugStop();
        }
      } else {
        facevec[iface].fLeftRightGelSideIndex = {rightgelside, leftgelside};
        facevec[iface].fLeftRightVolIndex = {VolumeElementIndex[rightindex], VolumeElementIndex[leftindex]};
        facevec[iface].fLeftRightGelIndex = {rightgel->Index(), leftgel->Index()};
      }
      globcount++;
    }
  }
  fVolumeElements.Resize(volumecount);

  int64_t nelcomp = fTransportMesh->NElements();
  for (int64_t el = 0; el < nelcomp; el++) {
    if (VolumeElementIndex[el] == -1) continue;
    int64_t cellindex = VolumeElementIndex[el];
    fVolumeElements[cellindex] = el;
  }

  int ncon = fFluxMesh->NConnects();
  int count = 0;
  for (int icon = 0; icon < ncon; icon++) {
    if (fConnectsByInterfaceMatID[icon] != -10000) {
      count++;
    }
    //        std::cout<<"conect: "<<icon<<" material: "<<fConnectsByInterfaceMatID[icon]<<std::endl;
  }
  if (count != globcount) {
    DebugStop();
  }
  int nnterfaces = fInterfaceByGeom.Rows();
}
void TSFDataTransfer::TakeOrientationAndLowerIndexDimVolDimFrac(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matID) {

  TPZCompEl *celL = celSideL.Element();
  if (celL->NConnects() > 1) DebugStop();
  TPZGeoEl *gelL = celL->Reference();
  int SideL = celSideL.Side();
  int nsidesm = gelL->NSides();
  orientationL = 0;

  TPZCompEl *CompelMixedL = gelL->Reference();
  TPZMultiphysicsElement *celmultL = dynamic_cast<TPZMultiphysicsElement *>(CompelMixedL);
  TPZCompEl *hdivBoundL = celmultL->Element(0);

  TPZCompElHDiv<pzshape::TPZShapeCube> *hdivboundC = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(hdivBoundL);
  TPZCompElHDiv<pzshape::TPZShapeTetra> *hdivboundT = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(hdivBoundL);
  TPZCompElHDiv<pzshape::TPZShapePrism> *hdivboundP = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapePrism> *>(hdivBoundL);
  TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivboundQ = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(hdivBoundL);
  TPZCompElHDiv<pzshape::TPZShapeTriang> *hdivboundTrian = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTriang> *>(hdivBoundL);
  int connect_index = 1;
  int firstsideL = gelL->NSides() - gelL->NSides(2) - 1;
  lowerIndexL = SideLowerIndex(gelL, SideL);
  if (hdivboundC) {
    connect_index = hdivboundC->ConnectIndex(SideL - firstsideL);
    orientationL = hdivboundC->SideOrient(SideL - firstsideL);
  } else if (hdivboundT) {
    connect_index = hdivboundT->ConnectIndex(SideL - firstsideL);
    orientationL = hdivboundT->SideOrient(SideL - firstsideL);
  } else if (hdivboundP) {
    connect_index = hdivboundP->ConnectIndex(SideL - firstsideL);
    orientationL = hdivboundP->SideOrient(SideL - firstsideL);
  } else if (hdivboundQ) {
    firstsideL = gelL->NSides() - gelL->NSides(1) - 1;
    connect_index = hdivboundQ->ConnectIndex(SideL - firstsideL);
    orientationL = hdivboundQ->SideOrient(SideL - firstsideL);
  } else if (hdivboundTrian) {
    firstsideL = gelL->NSides() - gelL->NSides(1) - 1;
    connect_index = hdivboundTrian->ConnectIndex(SideL - firstsideL);
    orientationL = hdivboundTrian->SideOrient(SideL - firstsideL);
  }

  TPZCompEl *celR = celSideR.Element();
  if (celR->NConnects() > 1) DebugStop();
  TPZGeoEl *gelR = celR->Reference();
  int SideR = celSideR.Side();
  orientationR = 0;
  TPZCompEl *CompelMixedR = gelR->Reference();
  TPZMultiphysicsElement *celmultR = dynamic_cast<TPZMultiphysicsElement *>(CompelMixedR);
  TPZCompEl *hdivBoundR = nullptr;
  if (celmultR) {
    hdivBoundR = celmultR->Element(0);
  }
  if (!celmultR || !hdivBoundR) {
    std::cout << "TransportElement Without fluxelement: INTERESECTION?" << std::endl;
    orientationR = gelR->NormalOrientation(SideR);
    orientationR = 0.0; //-orientationL;
    std::cout << "orientation" << orientationR << std::endl;
    lowerIndexR = SideLowerIndex(gelR, SideR);
    if (fConnectsByInterfaceMatID[connect_index] > 0) {
      DebugStop();
    } else {
      fConnectsByInterfaceMatID[connect_index] = matID;
    }
    return;
  }
  //    TPZCompEl *hdivBoundR = celmultR->Element(0);
  TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollapsedCR = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivBoundR);
  TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *hdivCollapsedTR = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *>(hdivBoundR);

  TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivBoundCR = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(hdivBoundR);
  TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivBoundTR = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(hdivBoundR);
  TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivBoundLinR = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(hdivBoundR);

  lowerIndexR = SideLowerIndex(gelR, SideR);
  bool found = false;
  int connect_indexR = -1;
  if (hdivCollapsedCR) {
    int index2 = hdivCollapsedCR->ConnectIndex(5);
    if (index2 == connect_index) {
      orientationR = hdivCollapsedCR->GetSideOrient(gelR->NSides() - 1);
      connect_indexR = hdivCollapsedCR->ConnectIndex(5);
      lowerIndexR = 4;
      found = true;
    }
    int index3 = hdivCollapsedCR->ConnectIndex(6);
    if (index3 == connect_index) {
      orientationR = hdivCollapsedCR->GetSideOrient(gelR->NSides());
      connect_indexR = hdivCollapsedCR->ConnectIndex(6);
      lowerIndexR = 5;
      found = true;
    }
  } else if (hdivCollapsedTR) {
    int index2 = hdivCollapsedTR->ConnectIndex(4);
    if (index2 == connect_index) {
      orientationR = hdivCollapsedTR->GetSideOrient(gelR->NSides() - 1);
      connect_indexR = hdivCollapsedTR->ConnectIndex(4);
      lowerIndexR = 3;
      found = true;
    }
    int index3 = hdivCollapsedTR->ConnectIndex(5);
    if (index3 == connect_index) {
      orientationR = hdivCollapsedTR->GetSideOrient(gelR->NSides());
      connect_indexR = hdivCollapsedTR->ConnectIndex(5);
      lowerIndexR = 4;
      found = true;
    }
  } else if (hdivBoundCR) {
    orientationR = 0;
    connect_indexR = hdivBoundCR->ConnectIndex(0);
    int ncon = hdivBoundCR->NConnects();
    lowerIndexR = SideLowerIndex(gelR, SideR);
    found = true;
  } else if (hdivBoundTR) {
    orientationR = 0;
    connect_indexR = hdivBoundTR->ConnectIndex(0);
    lowerIndexR = SideLowerIndex(gelR, SideR);
    found = true;

  } else if (hdivBoundLinR) {
    orientationR = 0;
    connect_indexR = hdivBoundLinR->ConnectIndex(0);
    lowerIndexR = SideLowerIndex(gelR, SideR);
    found = true;

  } else {
    DebugStop();
  }
  if (!found) {
    DebugStop();
  }
  if (fInterfaceByGeom(gelL->Index(), lowerIndexL) != -1) {
    DebugStop();
  }
  if (fInterfaceByGeom(gelR->Index(), lowerIndexR) != -1) {
    DebugStop();
  }
  if (connect_index != connect_indexR) {
    DebugStop();
  }
  if (fConnectsByInterfaceMatID[connect_index] > 0) {
    DebugStop();
  } else {
    fConnectsByInterfaceMatID[connect_index] = matID;
  }
}
void TSFDataTransfer::TakeOrientationAndLowerIndexDimDim(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matID) {

  TPZCompEl *celL = celSideL.Element();
  if (celL->NConnects() > 1) DebugStop();
  TPZGeoEl *gelL = celL->Reference();
  const int domaindim = gelL->Mesh()->Dimension();
  int SideL = celSideL.Side();
  int nsidesm = gelL->NSides();
  orientationL = 0;

  TPZCompEl *CompelMixedL = gelL->Reference();
  TPZMultiphysicsElement *celmultL = dynamic_cast<TPZMultiphysicsElement *>(CompelMixedL);
  TPZCompEl *hdivL = celmultL->Element(0);

  TPZCompElHDiv<pzshape::TPZShapeCube> *hdivCL = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(hdivL);
  TPZCompElHDiv<pzshape::TPZShapeTetra> *hdivTL = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(hdivL);
  TPZCompElHDiv<pzshape::TPZShapePrism> *hdivPL = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapePrism> *>(hdivL);
  TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivQuadL = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(hdivL);
  lowerIndexL = SideLowerIndex(gelL, SideL);
  int firstsideL = gelL->NSides() - gelL->NSides(domaindim - 1) - 1;
  int connect_indexL = -1;
  if (hdivCL) {
    orientationL = hdivCL->SideOrient(SideL - firstsideL);
    connect_indexL = hdivCL->ConnectIndex(SideL - firstsideL);
  } else if (hdivTL) {
    orientationL = hdivTL->SideOrient(SideL - firstsideL);
    connect_indexL = hdivTL->ConnectIndex(SideL - firstsideL);
  } else if (hdivPL) {
    orientationL = hdivPL->SideOrient(SideL - firstsideL);
    connect_indexL = hdivPL->ConnectIndex(SideL - firstsideL);
  } else if (hdivQuadL) {
    orientationL = hdivQuadL->SideOrient(SideL - firstsideL);
    connect_indexL = hdivQuadL->ConnectIndex(SideL - firstsideL);
  } else {
    DebugStop();
  }

  int connect_indexR = -2;
  TPZCompEl *celR = celSideR.Element();
  if (celR->NConnects() > 1) DebugStop();
  TPZGeoEl *gelR = celR->Reference();
  int SideR = celSideR.Side();
  int nsidesmR = gelR->NSides();
  orientationR = 0;
  lowerIndexR = SideLowerIndex(gelR, SideR);
  TPZCompEl *CompelMixedR = gelR->Reference();
  TPZMultiphysicsElement *celmultR = dynamic_cast<TPZMultiphysicsElement *>(CompelMixedR);
  TPZCompEl *hdivR = celmultR->Element(0);

  TPZCompElHDiv<pzshape::TPZShapeCube> *hdivCR = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(hdivR);
  TPZCompElHDiv<pzshape::TPZShapeTetra> *hdivTR = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(hdivR);
  TPZCompElHDiv<pzshape::TPZShapePrism> *hdivPR = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapePrism> *>(hdivR);
  TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivQuadR = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(hdivR);
  int firstsideR = gelR->NSides() - gelR->NSides(domaindim - 1) - 1;
  if (hdivCR) {
    orientationR = hdivCR->SideOrient(SideR - firstsideR);
    connect_indexR = hdivCR->ConnectIndex(SideR - firstsideR);
  } else if (hdivTR) {
    orientationR = hdivTR->SideOrient(SideR - firstsideR);
    connect_indexR = hdivTR->ConnectIndex(SideR - firstsideR);
  } else if (hdivPR) {
    orientationR = hdivPR->SideOrient(SideR - firstsideR);
    connect_indexR = hdivPR->ConnectIndex(SideR - firstsideR);
  } else if (hdivQuadR) {
    orientationR = hdivQuadR->SideOrient(SideR - firstsideR);
    connect_indexR = hdivQuadR->ConnectIndex(SideR - firstsideR);
  } else {
    DebugStop();
  }
  if (fInterfaceByGeom(gelL->Index(), lowerIndexL) != -1) {
    DebugStop();
  }
  if (fInterfaceByGeom(gelR->Index(), lowerIndexR) != -1) {
    DebugStop();
  }

  if (connect_indexR != connect_indexL) {
    DebugStop();
  }
  if (fConnectsByInterfaceMatID[connect_indexL] > 0) {
    DebugStop();
  } else {
    fConnectsByInterfaceMatID[connect_indexL] = matID;
  }
}

void TSFDataTransfer::TakeOrientationAndLowerIndexFracFrac(TPZCompElSide &celSideL, TPZCompElSide &celSideR, int &orientationL, int &lowerIndexL, int &orientationR, int &lowerIndexR, int matID) {

  TPZCompEl *celL = celSideL.Element();
  if (celL->NConnects() > 1) DebugStop();
  TPZGeoEl *gelL = celL->Reference();
  int SideL = celSideL.Side();
  int nsidesm = gelL->NSides();
  orientationL = 0;

  TPZCompEl *CompelMixedL = gelL->Reference();
  TPZMultiphysicsElement *celmultL = dynamic_cast<TPZMultiphysicsElement *>(CompelMixedL);
  TPZCompEl *hdivBoundL = celmultL->Element(0);

  TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollapsedCL = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivBoundL);
  TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *hdivCollapsedTL = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *>(hdivBoundL);

  int connect_index = 1;
  int firstsideL = gelL->NSides() - gelL->NSides(1) - 1;
  int connect_indexL = -1;
  if (hdivCollapsedCL) {
    orientationL = hdivCollapsedCL->SideOrient(SideL - firstsideL);
    connect_indexL = hdivCollapsedCL->ConnectIndex(SideL - firstsideL);
  } else if (hdivCollapsedTL) {
    orientationL = hdivCollapsedTL->SideOrient(SideL - firstsideL);
    connect_indexL = hdivCollapsedTL->ConnectIndex(SideL - firstsideL);
  } else {
    DebugStop();
  }
  lowerIndexL = SideLowerIndex(gelL, SideL);

  TPZCompEl *celR = celSideR.Element();
  if (celR->NConnects() > 1) DebugStop();
  TPZGeoEl *gelR = celR->Reference();
  int SideR = celSideR.Side();
  orientationR = 0;
  int connect_indexR = -2;
  TPZCompEl *CompelMixedR = gelR->Reference();
  TPZMultiphysicsElement *celmultR = dynamic_cast<TPZMultiphysicsElement *>(CompelMixedR);
  TPZCompEl *hdivBoundR = celmultR->Element(0);
  TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollapsedCR = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivBoundR);
  TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *hdivCollapsedTR = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *>(hdivBoundR);

  TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivBoundCR = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(hdivBoundR);
  TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivBoundTR = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(hdivBoundR);
  TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivBoundLinR = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(hdivBoundR);

  int firstsideR = gelR->NSides() - gelL->NSides(1) - 1;
  if (hdivCollapsedCR) {
    orientationR = hdivCollapsedCR->SideOrient(SideR - firstsideR);
    connect_indexR = hdivCollapsedCR->ConnectIndex(SideR - firstsideR);
  } else if (hdivCollapsedTR) {
    orientationR = hdivCollapsedTR->SideOrient(SideR - firstsideR);
    connect_indexR = hdivCollapsedTR->ConnectIndex(SideR - firstsideR);
  } else {
    DebugStop();
  }

  lowerIndexR = SideLowerIndex(gelR, SideR);
  if (fInterfaceByGeom(gelL->Index(), lowerIndexL) != -1) {
    DebugStop();
  }
  if (fInterfaceByGeom(gelR->Index(), lowerIndexR) != -1) {
    DebugStop();
  }
  if (connect_indexR != connect_indexL) {
    DebugStop();
  }
  if (fConnectsByInterfaceMatID[connect_indexL] > 0) {
    DebugStop();
  } else {
    fConnectsByInterfaceMatID[connect_indexL] = matID;
  }
}
// identify material of a face which is connected to a given geometric element
std::vector<int> TSFDataTransfer::IdentifyMaterial(TPZGeoElSideIndex gelindex, int64_t faceindex) {
  std::vector<int> vector;
  bool check = false;
  for (auto it = fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++) {
    TPZVec<TInterfaceWithVolume> &faces = it->second;
    if (faceindex >= faces.size()) continue;
    TInterfaceWithVolume &intface = faces[faceindex];
    if (intface.fLeftRightGelSideIndex.first == gelindex || intface.fLeftRightGelSideIndex.second == gelindex) {
      vector.push_back(it->first);
      check = true;
      //            return it->first;
    }
  }
  if (check) {
    return vector;
  }
  std::cout << __PRETTY_FUNCTION__ << "couldnt find a material corresponding to gelindex " << gelindex << " and faceindex " << faceindex << std::endl;
  DebugStop();
  return vector;
}

// identify material of a face which is connected to a given geometric element
std::vector<int> TSFDataTransfer::IdentifyMaterial(const TPZGeoElSide &gelside, int delta) {
  int64_t el = gelside.Element()->Index();
  int i = SideLowerIndex(gelside.Element(), gelside.Side());
  int64_t interface = fInterfaceByGeom(el, i);
  int nsides = gelside.Element()->NSides();
  TPZGeoElSide refside;
  if (gelside.Side() >= nsides - 1) {
    TPZGeoElSide gelsideaux(gelside.Element(), nsides - 1);
    refside = gelsideaux;
  } else {
    TPZGeoElSide gelsideaux(gelside);
    refside = gelsideaux;
  }

  if (interface < 0) {
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
      i = SideLowerIndex(neighbour.Element(), neighbour.Side()) + delta;
      interface = fInterfaceByGeom(neighbour.Element()->Index(), i);
      if (interface >= 0) {
        refside = neighbour;
        el = refside.Element()->Index();
        break;
      }
      neighbour = neighbour.Neighbour();
    }
    if (interface < 0) DebugStop();
  }
  TPZGeoElSideIndex gelsideindex(el, refside.Side());
  return IdentifyMaterial(gelsideindex, interface);
}

// find the neighbouring interface element
TPZGeoElSide TSFDataTransfer::IdentifyInterfaceElement(const TPZGeoElSide &gelside) {

  int64_t el = gelside.Element()->Index();
  int i = SideLowerIndex(gelside.Element(), gelside.Side());
  int64_t interface = fInterfaceByGeom(el, i);
  TPZGeoElSide refside(gelside);
  if (interface < 0) {
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
      i = SideLowerIndex(neighbour.Element(), neighbour.Side());
      interface = fInterfaceByGeom(neighbour.Element()->Index(), i);
      if (interface >= 0) {
        refside = neighbour;
        return refside;
      }
      neighbour = neighbour.Neighbour();
    }
    if (interface < 0) DebugStop();
  }
  return gelside;
}

static void ExtractElement(TPZCompEl *cel, std::list<TPZCompEl *> &ellist) {
  TPZFastCondensedElement *condF = dynamic_cast<TPZFastCondensedElement *>(cel);
  TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
  TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
  if (cond || condF) {
    ellist.push_back(cel);
  }
  if (elgr) {
    std::cout << "Error" << std::endl;
  }
  //    else if(gel)
  //    {
  //            ellist.push_back(cel);
  //
  //    }
  //
  //    else if(elgr)
  //    {
  //        const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
  //        int64_t nel = elvec.size();
  //        for(int64_t el=0; el<nel; el++)
  //        {
  //            ExtractElement(elvec[el],ellist);
  //        }
  //    }
  //    else
  //    {
  //        DebugStop();
  //    }
}
/// extract the list of computational element from a substructured computational mesh
// this method also searches for elements in element groups and condensed elements
// each computational element has an associated geometric element
void TSFDataTransfer::GetElementAndSubmeshPointers(TPZCompMesh &mixedmesh, std::list<TPZCompEl *> &elpointers, std::list<TPZSubCompMesh *> &submeshes) {
  int64_t nel = mixedmesh.NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mixedmesh.Element(el);
    if (!cel) continue;

    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
    if (submesh) {
      submeshes.push_back(submesh);
    } else {
      ExtractElement(cel, elpointers);
    }
  }
}

/// build the data structure from mixed to transport
void TSFDataTransfer::BuildMixedToTransportDataStructures(TPZCompMesh *fluxmesh) {
  if (!fluxmesh) DebugStop();

  std::list<TPZCompEl *> cellist;
  std::list<TPZSubCompMesh *> submeshes;
  GetElementAndSubmeshPointers(*fluxmesh, cellist, submeshes);
  TPZVec<int64_t> shouldtransfer(fluxmesh->NConnects(), 0);
  TPZVec<int64_t> targetindex(fluxmesh->NConnects(), 0);
  fluxmesh->LoadReferences();
  if (cellist.size()) {
    // compute the number of connects that will transfer information
    std::map<int, TPZManVector<int64_t, 4>> ncontransfer;
    std::map<int, TPZStack<int64_t>> connectindexes;
    // build the connect list to be transferred
    // identify the face index in the algebraic data structure
    for (auto it = cellist.begin(); it != cellist.end(); it++) {
      TPZCompEl *cel = *it;
      TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
      TPZCondensedCompEl *cnd = dynamic_cast<TPZCondensedCompEl *>(cel);
      TPZCompEl *candidate = NULL;
      TPZElementGroup *group = NULL;
      TPZMultiphysicsElement *mphys = NULL;
      if (condensed) {
        candidate = condensed->ReferenceCompEl();
        group = dynamic_cast<TPZElementGroup *>(candidate);
      }
      if (cnd) {
        mphys = dynamic_cast<TPZMultiphysicsElement *>(cnd->ReferenceCompEl());
        cel = mphys;
      }

      if (group) {
        mphys = findMultiphysics(group);
        condensed->SetIsGroup(true);
        condensed->SetMultiphysics(mphys);
        cel = mphys;
      }
      if (condensed && !group) {
        mphys = dynamic_cast<TPZMultiphysicsElement *>(candidate);
        condensed->SetMultiphysics(mphys);
        cel = mphys;
      }
      TPZGeoEl *gel = cel->Reference();
      TPZCompEl *hdiv = mphys->ReferredElement(0);
      if (!hdiv) {
        DebugStop();
      }
      TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *collapsedQuad = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdiv);
      TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *collapsedTriangle = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *>(hdiv);
      int nc = hdiv->NConnects();
      if (nc == 1) {
        // a boundary condition element
        // boundary elements are not considered
        continue;
      } else {
        int nsides = gel->NSides();
        int dim = gel->Dimension();
        int numfaces = gel->NSides(dim - 1);
        int firstside = nsides - numfaces - 1;
        if (collapsedQuad) {
          numfaces = 6;
        }
        if (collapsedTriangle) {
          numfaces = 5;
        }
        if (nc != numfaces + 1) DebugStop();
        int warning = 0;
        for (int ic = 0; ic < nc - 1; ic++) {

          if (collapsedQuad && ic == 4) {
            warning = 1;
            nc++;
            continue;
          }

          if (collapsedTriangle && ic == 3) {
            warning = 1;
            nc++;
            continue;
          }

          int side = firstside + ic;
          if (warning) {
            side = side - 1;
          }
          TPZGeoElSide gelside(gel, side);
          int nshape = cel->Connect(ic).NShape();
          int64_t cindex = cel->ConnectIndex(ic);

          if (shouldtransfer[cindex] != 0) continue;

          TPZSubCompMesh *scmesh = dynamic_cast<TPZSubCompMesh *>(fluxmesh);
          int matid = -10000;
          if (scmesh) {
            std::map<int64_t, int64_t>::iterator it = scmesh->LocalToFather().find(cindex);
            if (it != scmesh->LocalToFather().end()) {
              const int fatherConnect = it->second;
              matid = fConnectsByInterfaceMatID[fatherConnect];
            } else {
              DebugStop();
            }

          } else
            matid = fConnectsByInterfaceMatID[cindex];

          if (matid == -10000) {
            DebugStop();
          }

          if (ncontransfer.find(matid) == ncontransfer.end()) {
            ncontransfer[matid].Resize(4, 0);
          }

          connectindexes[matid].Push(cindex);

          for (int is = 0; is < nshape; is++) {
            ncontransfer[matid][is]++;
          }
          shouldtransfer[cindex] = 1;

          TPZGeoElSide refside = IdentifyInterfaceElement(gelside);
          int64_t elindex = refside.Element()->Index();
          int sidelower = SideLowerIndex(refside.Element(), refside.Side());
          int targetind = fInterfaceByGeom(elindex, sidelower);
          targetindex[cindex] = targetind;
          if (targetindex[cindex] < 0) DebugStop();
        }
      }
    }
    for (auto &it : ncontransfer) {
      int matid = it.first;
      int nc = connectindexes[matid].size();
      for (int is = 0; is < 4; is++) {
        if (ncontransfer[matid][is] == 0) continue;
        TFromMixedToTransport transport;
        transport.fMatid = matid;
        transport.fGather.resize(ncontransfer[matid][is]);
        transport.fScatter.resize(ncontransfer[matid][is]);
        transport.fMixedMesh = fluxmesh;
        transport.flux_sequence = is;
        int64_t count = 0;
        for (int64_t ic = 0; ic < nc; ic++) {
          int64_t cindex = connectindexes[matid][ic];
          int64_t seqnum = fluxmesh->ConnectVec()[cindex].SequenceNumber();
          int64_t position = fluxmesh->Block().Position(seqnum);
          int64_t blsize = fluxmesh->Block().Size(seqnum);
          if (is >= blsize) continue;
          if (count >= ncontransfer[matid][is]) DebugStop();
          int gather = position + is;
          transport.fGather[count] = gather;
          int scatter = targetindex[cindex];
          transport.fScatter[count] = scatter;
          count++;
        }
        if (count != ncontransfer[matid][is]) DebugStop();
        fTransferMixedToTransport[matid].push_back(transport);
      }
    }
  }
  if (submeshes.size()) {
    for (auto it = submeshes.begin(); it != submeshes.end(); it++) {
      BuildMixedToTransportDataStructures(*it);
    }
  }
}

// Initialize the pointers to the transport data structure in the
// list fTransferMixedToTransport
void TSFDataTransfer::InitializeVectorPointersMixedToTransport(TPZAlgebraicTransport &transport) {
  for (auto &mat_iter : fTransferMixedToTransport) {
    int matid = mat_iter.first;
    if (transport.fInterfaceData.find(matid) == transport.fInterfaceData.end()) DebugStop();
    for (auto &list_iter : mat_iter.second) {
      int flux_index = list_iter.flux_sequence;
      if (transport.fInterfaceData[matid].fCoefficientsFlux.size() <= flux_index) DebugStop();
      auto vecptr = &transport.fInterfaceData[matid].fCoefficientsFlux[flux_index];
      list_iter.fTarget = vecptr;
      TPZFMatrix<STATE> &matptr = (list_iter.fMixedMesh->Solution());
      list_iter.fFrom = &matptr;
    }
  }
}

// Initialize the pointers from the transport data structure in the list TransportToMixedCorrespondence
void TSFDataTransfer::InitializeVectorPointersTranportToMixed(TPZAlgebraicTransport &transport) {
  for (auto &mesh_iter : fTransportMixedCorrespondence) {
    mesh_iter.fTransport = &transport;
  }
}

// print the datastructure
void TSFDataTransfer::Print(std::ostream &out) {
  TPZGeoMesh *gmesh = fTransportMesh->Reference();
  out << "Number of interface materials " << fInterfaceGelIndexes.size() << std::endl;
  for (auto it = fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++) {
    out << "Element indexes for material id " << it->first << std::endl;
    TPZVec<TInterfaceWithVolume> &gelindex = it->second;
    int64_t nel = gelindex.NElements();
    for (int64_t el = 0; el < nel; el++) {
      TInterfaceWithVolume &intface = gelindex[el];
      int64_t gindex = intface.fInterface_gelindex;
#ifdef PZDEBUG
      if (gindex < 0) DebugStop();
#endif
      TPZGeoEl *gel = gmesh->Element(gindex);
      int matid = gel->MaterialId();
      int ncorner = gel->NCornerNodes();
      out << "el = " << el << " gel index " << intface.fInterface_gelindex << " ncorner " << ncorner << " matid " << matid << std::endl;
      out << "         cel index " << intface.fInterface_celindex << " leftright gelindex " << intface.fLeftRightGelSideIndex.first << " " << intface.fLeftRightGelSideIndex.second << std::endl;
      out << "         leftright alg vol index " << intface.fLeftRightVolIndex << std::endl;
    }
  }
  out << "fVolumeElements Number of volume materials " << fVolumeElements.size() << std::endl;

  out << "Volume elements" << std::endl;
  TPZVec<int64_t> &volel_indexes = fVolumeElements;
  int64_t nel = volel_indexes.size();
  for (int64_t el = 0; el < nel; el++) {
    int64_t celindex = volel_indexes[el];
    TPZCompEl *cel = fTransportMesh->Element(celindex);
    TPZGeoEl *gel = cel->Reference();
    int64_t gelindex = gel->Index();
    out << "el = " << el << " gel index " << gelindex << " matid " << gel->MaterialId() << " dim " << gel->Dimension() << std::endl;
  }

  if (fInterfaceByGeom.Rows() == 0) {
    out << "fInterfaceByGeom not initialized\n";
  } else {
    out << "fInterfaceByGeom For each geometric element, which are the algebraic faces connected to it\n";
    int64_t nel_geo = gmesh->NElements();
    for (int64_t el = 0; el < nel_geo; el++) {
      TPZGeoEl *gel = gmesh->Element(el);
      bool hasface = false;
      for (int i = 0; i < 6; i++)
        if (fInterfaceByGeom(el, i) != -1) {
          hasface = true;
        }
      if (hasface) {
        out << "gel index = " << el << std::endl;
        for (int i = 0; i < 6; i++)
          if (fInterfaceByGeom(el, i) != -1) {
            int side = SideOriginalIndex(gel, i);
            TPZGeoElSideIndex gelside(el, side);
            //                    out << "i = " << i <<  " side " << side << " face index " << fInterfaceByGeom(el,i) << " matid "
            //                    << IdentifyMaterial(gelside,fInterfaceByGeom(el,i)) << std::endl;
          }
      }
    }
  }
  if (fTransferMixedToTransport.size()) {
    out << "Gather scatter to bring the flux data to the transport mesh\n";
    for (auto it : fTransferMixedToTransport) {
      int matid = it.first;
      out << "Flux data transfer for material id : " << matid << std::endl;
      std::list<TFromMixedToTransport> &list = it.second;
      for (auto itlist : list) {
        TFromMixedToTransport &transport = itlist;
      }
    }
  }
}

// Build the data structure which defines the correspondence between
// algebraic transport cells and indexes of mixed fast condensed elements
void TSFDataTransfer::BuildTransportToMixedCorrespondenceDatastructure(TPZCompMesh *fluxmesh, TPZVec<int64_t> &Alg_Cell_Index) {
  // build a vector of the size of the geometric elements
  // where applicable the value of the vector is the index of the
  // algebraic transport volume
  TPZGeoMesh *gmesh = fTransportMesh->Reference();
  int64_t nel_geo = gmesh->NElements();
  //    TPZVec<int64_t> Volume_Index(nel_geo,-10);
  if (Alg_Cell_Index.size() != nel_geo) {
    Alg_Cell_Index.Resize(nel_geo);
    Alg_Cell_Index.Fill(-10);
    for (int64_t vec_it = 0; vec_it < fVolumeElements.size(); vec_it++) {
      int64_t celindex = fVolumeElements[vec_it];
      TPZCompEl *cel = fTransportMesh->Element(celindex);
      if (cel->NConnects() != 0) {
        TPZGeoEl *gel = cel->Reference();
        int64_t geo_index = gel->Index();
        Alg_Cell_Index[geo_index] = vec_it;
      }
    }
  }
  if (!fluxmesh) DebugStop();
  std::list<TPZCompEl *> cellist;
  std::list<TPZSubCompMesh *> submeshes;
  GetElementAndSubmeshPointers(*fluxmesh, cellist, submeshes);
  if (cellist.size()) {
    int64_t num_elements = 0;
    // first count the number of volume cells in the mesh
    for (auto cel : cellist) {
      // skip boundary elements
      if (cel->NConnects() == 1) continue;
      TPZFastCondensedElement *fastel = dynamic_cast<TPZFastCondensedElement *>(cel);
      TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
      if (!fastel && !cond) DebugStop();

      TPZCompEl *celaux = NULL;
      if (fastel) {
        celaux = fastel->GetMultiphysics();
      } else {
        celaux = dynamic_cast<TPZMultiphysicsElement *>(cond->ReferenceCompEl());
      }

      int64_t celindex2 = cel->Index();
      int64_t gelindex = celaux->Reference()->Index();
      if (Alg_Cell_Index[gelindex] < 0) continue;
      TPZCompEl *celcondensed = fluxmesh->Element(celindex2);
      TPZFastCondensedElement *fast = dynamic_cast<TPZFastCondensedElement *>(celcondensed);
      TPZCondensedCompEl *cnd = dynamic_cast<TPZCondensedCompEl *>(celcondensed);
      if (!fast && !cnd) DebugStop();
      num_elements++;
    }
    if (num_elements > 0) {
      TransportToMixedCorrespondence transport;
      transport.fMixedMesh = fluxmesh;
      transport.fMixedCell.resize(num_elements);
      transport.fAlgebraicTransportCellIndex.resize(num_elements);
      transport.fEqNum.resize(num_elements);

      int64_t count = 0;
      for (auto cel : cellist) {
        // skip boundary elements
        if (cel->NConnects() == 1) continue;
        TPZFastCondensedElement *fastel = dynamic_cast<TPZFastCondensedElement *>(cel);
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (!fastel && !cond) DebugStop();
        TPZCompEl *celaux = NULL;
        if (fastel) {
          celaux = fastel->GetMultiphysics();
        } else {
          celaux = dynamic_cast<TPZMultiphysicsElement *>(cond->ReferenceCompEl());
        }

        int64_t celindex = cel->Index();
        int64_t gelindex = celaux->Reference()->Index();
        if (Alg_Cell_Index[gelindex] < 0) continue;
        TPZCompEl *celcondensed = fluxmesh->Element(celindex);
        TPZFastCondensedElement *fast = dynamic_cast<TPZFastCondensedElement *>(celcondensed);
        TPZCondensedCompEl *condensed = dynamic_cast<TPZCondensedCompEl *>(celcondensed);
        if (!fast && !condensed) continue;
        if (fast) {
          transport.fMixedCell[count] = fast;
        } else {
          transport.fMixedCell[count] = condensed;
        }

        int64_t AlgCelIndex = Alg_Cell_Index[gelindex];
        int64_t TransportCelIndex = fVolumeElements[AlgCelIndex];
        TPZGeoEl *fastgel = celaux->Reference();
        int dim = fastgel->Dimension();
        int side = fastgel->NSides() - 1;
        TPZVec<REAL> coordFast(3), xifast(dim, 0);
        fastgel->CenterPoint(side, xifast);
        fastgel->X(xifast, coordFast);

        TPZCompEl *transComp = fTransportMesh->Element(TransportCelIndex);
        TPZGeoEl *transportGel = transComp->Reference();
        TPZVec<REAL> coordTransp(3), xiTransp(dim, 0);
        transportGel->CenterPoint(side, xiTransp);
        transportGel->X(xiTransp, coordTransp);
        REAL diff = 0.0;
        for (int i = 0; i < 3; i++) {
          diff += fabs(coordTransp[i] - coordFast[i]);
        }
        if (diff > 1.e-16) {
          DebugStop();
        }

        transport.fAlgebraicTransportCellIndex[count] = Alg_Cell_Index[gelindex];
        count++;
      }
      if (count != num_elements) DebugStop();
      fTransportMixedCorrespondence.push_back(transport);
    }
  }
  if (submeshes.size() != 0) {
    for (auto sub : submeshes) {
      BuildTransportToMixedCorrespondenceDatastructure(sub, Alg_Cell_Index);
    }
  }
}

void TSFDataTransfer::InitializeTransportDataStructure(TPZAlgebraicTransport &transport) {

  TPZGeoMesh *gmesh = fTransportMesh->Reference();
  for (auto mat_iter : fInterfaceGelIndexes) {
    TPZManVector<int, 6> numfaces(4, 0);
    TPZAlgebraicTransport::TInterfaceDataTransport &InterfaceVec = transport.fInterfaceData[mat_iter.first];
    int ncormax = 0;
    for (auto face_it : mat_iter.second) {

      TPZGeoEl *gel = gmesh->Element(face_it.fInterface_gelindex);
      TPZCompEl *cel = gel->Reference();
      TPZMultiphysicsInterfaceElement *intel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
      TPZInterfaceElement *intel2 = dynamic_cast<TPZInterfaceElement *>(cel);
      if (!intel && !intel2) {
        DebugStop();
      }
      TPZVec<REAL> normal;
      if (intel) {
        intel->ComputeCenterNormal(normal);
      } else {
        intel2->CenterNormal(normal);
        //                intel2->ComputeCenterNormal(normal);
      }

      std::tuple<REAL, REAL, REAL> norm = std::make_tuple(normal[0], normal[1], normal[2]);
      InterfaceVec.fNormalFaceDirection.push_back(norm);
      int ncorner = gel->NCornerNodes();
      std::pair<int, int> lr = face_it.fLeftRightVolIndex;
      InterfaceVec.fLeftRightVolIndex.push_back(lr);
      InterfaceVec.fcelindex.push_back(face_it.fInterface_celindex);

      std::pair<int, int> lrgel = face_it.fLeftRightGelIndex;
      InterfaceVec.fLeftRightGelIndex.push_back(lrgel);

      if (ncorner > ncormax) ncormax = ncorner;
      for (int i = 0; i < ncorner; i++)
        numfaces[i]++;
    }

    InterfaceVec.fMatid = mat_iter.first;
    InterfaceVec.fCoefficientsFlux.resize(ncormax);
    for (int i = 0; i < ncormax; i++)
      InterfaceVec.fCoefficientsFlux[i].resize(numfaces[i]);
    InterfaceVec.fIntegralFlux.resize(numfaces[0], 0.0);
    InterfaceVec.fFluxSign.resize(numfaces[0]);
    InterfaceVec.fNormalFaceDirection.resize(numfaces[0]);
    InterfaceVec.fcelindex.resize(numfaces[0]);
    InterfaceVec.fIntegralFluxFunctions.resize(numfaces[0]);
  }

  auto &volData = fVolumeElements;

  int64_t nvols = volData.size();

  transport.fCellsData.SetNumCells(nvols);
  transport.fCellsData.fViscosity.resize(2);
  transport.fCellsData.fViscosity[0] = 1.0;
  transport.fCellsData.fViscosity[1] = 1.0;
  std::cout << __PRETTY_FUNCTION__ << "Filling in property map\n";

  for (int64_t i = 0; i < nvols; i++) {
    int64_t celindex = volData[i];

    TPZCompEl *cel = fTransportMesh->Element(celindex);
    TPZGeoEl *gel = cel->Reference();
    int indexgeo = gel->Index();
    int geldim = gel->Dimension();
    int matId = gel->MaterialId();

    REAL volume = gel->Volume();
    if (transport.fCellsData.fsim_data->mTNumerics.m_is_axisymmetric) {
      int ncorner = gel->NCornerNodes();
      REAL rmin = std::numeric_limits<REAL>::max();
      REAL rmax = std::numeric_limits<REAL>::min();
      for (int ic = 0; ic < ncorner; ic++) {
        REAL x = gel->NodePtr(ic)->Coord(0);
        if (fabs(x) < rmin) rmin = x;
        if (fabs(x) > rmax) rmax = x;
      }
      if (rmin * rmax < 0.0) DebugStop();
      REAL h = volume / (rmax - rmin);
      volume = h * M_PI * (rmax * rmax - rmin * rmin);
    }
    int side = gel->NSides() - 1;
    transport.fCellsData.fVolume[i] = volume;
    transport.fCellsData.fMatId[i] = matId;
    if (matId == 299) {
      std::cout << "IntersecElement: " << std::endl;
    }
    if (cel->NConnects() != 1) {
      DebugStop();
    }
    TPZConnect &con = cel->Connect(0);
    int block_num = con.SequenceNumber();
    int eq_number = fTransportMesh->Block().Position(block_num);

    transport.fCellsData.fEqNumber[i] = eq_number;
    transport.fCellsData.fDensityOil[i] = 1000.00;
    transport.fCellsData.fDensityWater[i] = 1000.00;

    int dim = gel->Dimension();
    transport.fCellsData.fCenterCoordinate[i].resize(3);
    TPZVec<REAL> ximasscent(dim);
    gel->CenterPoint(side, ximasscent);
    std::vector<REAL> center(3, 0.0);
    TPZManVector<REAL, 3> coord(3, 0.0);
    gel->X(ximasscent, coord);

    REAL s0_v = 0.0;
    auto s0_func = transport.fCellsData.fsim_data->mTReservoirProperties.s0;
    if (s0_func) {
      fs0 = s0_func;
      s0_v = fs0(coord);
    }

    transport.fCellsData.fSaturation[i] = s0_v;
    transport.fCellsData.fSaturationLastState[i] = s0_v;

    for (int ic = 0; ic < 3; ic++) {
      center[ic] = coord[ic];
    };
    transport.fCellsData.fCenterCoordinate[i] = center;
    transport.fCellsData.fGeoIndex[i] = indexgeo;
  }
  //    transport.fCellsData.fMatId = 1;

  // transport.fCellsData.UpdateFractionalFlowsAndLambda(true);
  this->InitializeVectorPointersTranportToMixed(transport);
}

TSFDataTransfer::TFromMixedToTransport::TFromMixedToTransport() : fMatid(-1),
                                                                  flux_sequence(-1), fFrom(0), fTarget(0), fMixedMesh(0) {
}

TSFDataTransfer::TFromMixedToTransport::TFromMixedToTransport(const TFromMixedToTransport &copy) : fMatid(copy.fMatid), flux_sequence(copy.flux_sequence),
                                                                                                   fGather(copy.fGather), fScatter(copy.fScatter), fFrom(copy.fFrom),
                                                                                                   fTarget(copy.fTarget), fMixedMesh(copy.fMixedMesh) {
  fMixedMesh = copy.fMixedMesh;
}

TSFDataTransfer::TFromMixedToTransport &TSFDataTransfer::TFromMixedToTransport::operator=(const TFromMixedToTransport &copy) {
  fMatid = copy.fMatid;
  flux_sequence = copy.flux_sequence;
  fFrom = copy.fFrom;
  fTarget = copy.fTarget;
  fMixedMesh = copy.fMixedMesh;
  return *this;
}

void TSFDataTransfer::TFromMixedToTransport::Print(std::ostream &out) {
  out << "FromMixedToTransport mesh = " << (void *)fMixedMesh << " flux index " << flux_sequence << " matid " << fMatid << std::endl;
  out << "Gather vector ";
  for (auto const &value : fGather)
    out << value << ' ';
  out << std::endl;
  out << "Scatter vector ";
  for (auto const &value : fScatter)
    out << value << ' ';
  out << std::endl;
}

void TSFDataTransfer::TransportToMixedCorrespondence::Print(std::ostream &out) {
  out << "TransportToMixedCorrespondence mesh = " << (void *)fMixedMesh << std::endl;
  out << "fTransportCell vector ";
  for (auto value : fAlgebraicTransportCellIndex)
    out << value << ' ';
  out << std::endl;
  out << "fMixedCell compel indices ";
  for (auto value : fMixedCell)
    out << value->Index() << ' ';
  out << std::endl;
}

// transfer the solution from the mixed mesh fluxes to the interfaces
void TSFDataTransfer::TransferMixedMeshMultiplyingCoefficients() {
  for (auto &matit : fTransferMixedToTransport) {
    for (auto &mixed : matit.second) {
      int64_t nel = mixed.fGather.size();
      for (int64_t el = 0; el < nel; el++) {

        (*mixed.fTarget)[mixed.fScatter[el]] = (*mixed.fFrom)(mixed.fGather[el], 0);
      }
    }
  }
}
void TSFDataTransfer::TransferPressures() {
  for (auto &meshit : fTransportMixedCorrespondence) {
    int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
#ifdef PZDEBUG
    if (meshit.fTransport == 0) {
      DebugStop();
    }
#endif
    for (int icell = 0; icell < ncells; icell++) {
      TPZCompEl *cel = meshit.fMixedCell[icell];
      TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
      if (!condensed) {
        TPZCondensedCompEl *condcompel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condcompel) {
          std::cout << "Element " << cel->Index() << " is not FastCondensed" << std::endl;
          continue;
        } else {
          DebugStop();
        }
      }
      int64_t cellindex = meshit.fAlgebraicTransportCellIndex[icell];
      TPZCompEl *compel = condensed->ReferenceCompEl();
      int dim = compel->Dimension();
      TPZVec<REAL> qsi(dim, 0.0);
      TPZVec<STATE> sol(dim, 0.0);
      int presureindex = 2;
      compel->Solution(qsi, presureindex, sol);
      meshit.fTransport->fCellsData.fPressure[cellindex] = sol[0];
    }
  }
}
// transfer the permeability multiplier from the transport mesh to the mixed mesh elements
void TSFDataTransfer::TransferLambdaCoefficients() {

  for (auto &meshit : fTransportMixedCorrespondence) {
    int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
#ifdef PZDEBUG
    if (meshit.fTransport == 0) {
      DebugStop();
    }
#endif
    for (int icell = 0; icell < ncells; icell++) {
      TPZCompEl *cel = meshit.fMixedCell[icell];
      TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
      if (!condensed) {
        TPZCondensedCompEl *condcompel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condcompel) {
          std::cout << "Element " << cel->Index() << " is not FastCondensed" << std::endl;
          continue;
        } else {
          DebugStop();
        }
      }
      int64_t cellindex = meshit.fAlgebraicTransportCellIndex[icell];
      REAL lambda = meshit.fTransport->fCellsData.flambda[cellindex];
      condensed->SetLambda(lambda);

      REAL mixedDensity = meshit.fTransport->fCellsData.fMixedDensity[cellindex];
      condensed->SetMixedDensity(mixedDensity);

      REAL porosity = meshit.fTransport->fCellsData.fporosity[cellindex];
      REAL dt = meshit.fTransport->fdt;
      REAL sw = meshit.fTransport->fCellsData.fSaturation[cellindex];
      REAL so = 1 - sw;
      REAL drhoWdp = meshit.fTransport->fCellsData.fdDensityWaterdp[cellindex];
      REAL drhoOdp = meshit.fTransport->fCellsData.fdDensityOildp[cellindex];
      REAL compterm = -(porosity / dt) * ((sw * drhoWdp) + (so * drhoOdp)); // negative sign in accordance with the lyx

      REAL swlast = meshit.fTransport->fCellsData.fSaturationLastState[cellindex];
      REAL solast = 1.0 - swlast;
      REAL plast = meshit.fTransport->fCellsData.fPressure[cellindex];
      REAL rhoW = meshit.fTransport->fCellsData.fDensityWater[cellindex];
      REAL rhoO = meshit.fTransport->fCellsData.fDensityOil[cellindex];
      REAL rhoWlast = meshit.fTransport->fCellsData.fDensityWaterLastState[cellindex];
      REAL rhoOlast = meshit.fTransport->fCellsData.fDensityOilLastState[cellindex];
      REAL termrhscurrent = (sw * rhoW) + (so * rhoO);
      REAL termrhslast = (swlast * rhoWlast) + (solast * rhoOlast);
      REAL comptermrhs = (porosity / dt) * (termrhscurrent - termrhslast) + compterm * plast;

      condensed->SetCompressibiilityTerm(compterm, comptermrhs);
    }
  }
}

void TSFDataTransfer::TransferSaturation() {

  for (auto &meshit : fTransportMixedCorrespondence) {
    int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
#ifdef PZDEBUG
    if (meshit.fTransport == 0) {
      DebugStop();
    }
#endif
    for (int icell = 0; icell < ncells; icell++) {
      TPZCompEl *cel = meshit.fMixedCell[icell];
      TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
      if (!condensed) {
        TPZCondensedCompEl *condcompel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condcompel) {
          std::cout << "Element " << cel->Index() << " is not FastCondensed" << std::endl;
          continue;
        } else {
          DebugStop();
        }
      }
      int64_t cellindex = meshit.fAlgebraicTransportCellIndex[icell];
      REAL sw = meshit.fTransport->fCellsData.fSaturation[cellindex];
      condensed->SetSw(sw);
    }
  }
}

void TSFDataTransfer::TransferPermeabiliyTensor() {
  CheckDataTransferTransportToMixed();
  for (auto &meshit : fTransportMixedCorrespondence) {
    TPZAlgebraicTransport::TCellData &celldata = meshit.fTransport->fCellsData;
    int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
    for (int icell = 0; icell < ncells; icell++) {
#ifdef PZDEBUG
      if (meshit.fTransport == 0 || meshit.fEqNum[icell] >= ncells) {
        DebugStop();
      }
#endif
      TPZCompEl *cel = meshit.fMixedCell[icell];
      TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
      if (!condensed) {
        TPZCondensedCompEl *condcompel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condcompel) {
          std::cout << "Element " << cel->Index() << " is not FastCondensed" << std::endl;
          continue;
        } else {
          DebugStop();
        }
      }

      TPZFNMatrix<9, REAL> PermeabilityT, InvPerm;
      PermeabilityT.Redim(3, 3);
      InvPerm.Redim(3, 3);
      int64_t transportcell = meshit.fAlgebraicTransportCellIndex[icell];
      PermeabilityT(0, 0) = (celldata.fKx)[transportcell];
      PermeabilityT(1, 1) = (celldata.fKy)[transportcell];
      PermeabilityT(2, 2) = (celldata.fKz)[transportcell];
      InvPerm(0, 0) = 1.0 / (PermeabilityT(0, 0));
      InvPerm(1, 1) = 1.0 / (PermeabilityT(1, 1));
      InvPerm(2, 2) = 1.0 / (PermeabilityT(2, 2));
      TPZGeoEl *gel = cel->Reference();
      //       meshit.fMixedCell[icell]->SetPermTensorAndInv(PermeabilityT,InvPerm) ;
      condensed->SetPermTensorAndInv(PermeabilityT, InvPerm);
      condensed->SetMixedDensity((celldata.fMixedDensity)[transportcell]);
    }
  }
}

// verify the correspondence of the mixed elements and the algebraic cells
void TSFDataTransfer::CheckDataTransferTransportToMixed() {
  for (auto &meshit : fTransportMixedCorrespondence) {
    TPZAlgebraicTransport::TCellData &celldata = meshit.fTransport->fCellsData;
    int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
    for (int icell = 0; icell < ncells; icell++) {
#ifdef PZDEBUG
      if (meshit.fTransport == 0 || meshit.fEqNum[icell] >= ncells) {
        DebugStop();
      }
#endif
      int64_t transportcell = meshit.fAlgebraicTransportCellIndex[icell];
      TPZCompEl *cel = meshit.fMixedCell[icell];

#ifdef PZDEBUG
      TPZFastCondensedElement *condensed = dynamic_cast<TPZFastCondensedElement *>(cel);
      if (!condensed) {
        TPZCondensedCompEl *condcompel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (condcompel) {
          std::cout << "Element " << cel->Index() << " is not FastCondensed" << std::endl;
          continue;
        } else {
          DebugStop();
        }
      }

      TPZCompEl *celaux = condensed->GetMultiphysics();
      TPZGeoEl *gel = celaux->Reference();
      if (!gel) DebugStop();
      TPZManVector<REAL, 3> xcenter(3, 0.), xi(gel->Dimension(), 0.);
      gel->CenterPoint(gel->NSides() - 1, xi);
      gel->X(xi, xcenter);
      REAL diff = 0;
      std::vector<REAL> celcenter = celldata.fCenterCoordinate[transportcell];
      for (int i = 0; i < 3; i++) {
        diff += fabs(xcenter[i] - celcenter[i]);
      }
      if (diff > 1.e-16) {
        DebugStop();
      }
#endif
    }
  }
}

TPZMultiphysicsElement *TSFDataTransfer::findMultiphysics(TPZElementGroup *group) {

  TPZVec<TPZCompEl *> groupels = group->GetElGroup();
  int dim_mesh = group->Mesh()->Dimension();
  TPZMultiphysicsElement *mult;
  for (auto els : groupels) {
    int dim_el = els->Dimension();
    if (dim_el == dim_mesh) {
      TPZCondensedCompEl *condensed = dynamic_cast<TPZCondensedCompEl *>(els);
      if (condensed) {
        TPZCompEl *cel = condensed->ReferenceCompEl();
        if (cel) {
          TPZElementGroup *elgroup2 = dynamic_cast<TPZElementGroup *>(cel);
          if (elgroup2) {
            mult = findMultiphysics(elgroup2);
          }
        }
      } else {
        mult = dynamic_cast<TPZMultiphysicsElement *>(els);
        if (!mult) {
          DebugStop();
        }

        return mult;
      }
    }
  }
  return mult;
}
std::pair<int, int> TSFDataTransfer::FindMortar(TPZGeoElSide &gelside) {

  int mortarId = 40;
  std::set<int> idstargets;
  idstargets.insert(mortarId);
  idstargets.insert(-1);
  idstargets.insert(-2);
  idstargets.insert(-3);
  idstargets.insert(-4);

  std::pair<int, int> indextargetPair(-1, -1);
  for (auto idtar : idstargets) {
    std::pair<int, int> candidatetarPair = FindMortar(gelside, idtar);

    int candidatetar = std::get<0>(candidatetarPair);
    if (candidatetar != -1) {
      indextargetPair = candidatetarPair;
      break;
    }
  }
  return indextargetPair;
}
std::pair<int, int> TSFDataTransfer::FindMortar(TPZGeoElSide &gelside, int targetId) {

  int indexcon = -1;
  int nshape = -1;
  std::pair<int, int> pairret(-1, -1);
  int nsides = gelside.Element()->NSides();
  gelside.Element()->Type();
  int nsidesgelside = gelside.Side();
  auto elType = gelside.Element()->Type();
  bool first = true;
  if (nsidesgelside == 9 && elType == EQuadrilateral) {
    first = false;
  }
  if (nsidesgelside == 7 && elType == ETriangle) {
    first = false;
  }
  TPZGeoElSide gelsideAux;
  if (nsidesgelside < 8 && elType == EQuadrilateral) {
    std::pair<int, int> pairt(-1, -1);
    return pairt;
  }
  if (nsidesgelside < 6 && elType == ETriangle) {
    std::pair<int, int> pairt(-1, -1);
    return pairt;
  }
  if (nsidesgelside >= nsides && gelside.Element()->Dimension() == 2) {
    if (elType == EQuadrilateral) {
      gelsideAux = TPZGeoElSide(gelside.Element(), 8);
    }
    if (elType == ETriangle) {
      gelsideAux = TPZGeoElSide(gelside.Element(), 6);
    }

  } else {
    gelsideAux = gelside;
  }
  TPZGeoElSide neig = gelsideAux.Neighbour();
  int count = 0;
  while (gelsideAux != neig) {
    int elId = neig.Element()->MaterialId();
    //
    if (elId == targetId) {
      TPZCompEl *comp = neig.Element()->Reference();
      TPZMultiphysicsInterfaceElement *mult = dynamic_cast<TPZMultiphysicsInterfaceElement *>(comp);
      if (comp && !mult) {
        count++;
        if (!first) {
          neig = neig.Neighbour();
          first = true;
          continue;
        }
        int ncon = comp->NConnects();
        if (ncon != 1) {
          DebugStop();
        }
        indexcon = comp->ConnectIndex(0);
        nshape = comp->Connect(0).NShape();
        std::pair<int, int> pairp(indexcon, nshape);
        pairret = pairp;
        break;
      }
    }
    neig = neig.Neighbour();
  }

  return pairret;
  //    compelside
}

void TSFDataTransfer::TestSideOrient(TPZCompMesh *MultFlux) {
  int nels = MultFlux->NElements();
  //    for(int iel =0; iel < nels; iel++){
  //        TPZCompEl * cel = MultFlux->Element(iel);
  std::list<TPZCompEl *> elpointers;
  std::list<TPZSubCompMesh *> submeshes;
  GetElementAndSubmeshPointers(*MultFlux, elpointers, submeshes);
  MultFlux->LoadReferences();
  int count = 0;
  for (auto cel : elpointers) {

    //        TPZCompEl *cel =neigel->Reference();
    //        TPZMultiphysicsElement *cmulint = dynamic_cast<TPZMultiphysicsElement *>(cel);
    //        TPZCompEl *hdivBound = cmulint->Element(0);
    //        TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeQuad>*>(cel);
    //        int sideOrient =hdivbound->GetSideOrient(8);
    //        std::cout<<"SideOrientSide1: "<<sideOrient <<std::endl;

    if (count == 0) {
      TPZGeoEl *gel = cel->Reference();
      int nsides = gel->NSides();
      TPZGeoElSide gelside(gel, nsides - 1);
      TPZGeoElSide neig = gelside.Neighbour();

      while (gelside != neig) {
        TPZGeoEl *neigel = neig.Element();
        int GelIndex = neigel->Index();
        if (neigel->Dimension() == 3) {
          std::cout << "GelIndex: " << GelIndex << std::endl;
          TPZVec<REAL> coordFast(3), xifast(3, 0);
          neigel->CenterPoint(neigel->NSides() - 1, xifast);
          neigel->X(xifast, coordFast);
          std::cout << "CoordX: " << coordFast[0] << std::endl;
          std::cout << "CoordY: " << coordFast[1] << std::endl;
          std::cout << "CoordZ: " << coordFast[2] << std::endl;
          int sideOrient = neigel->NormalOrientation(neig.Side());
          std::cout << "SideOrientVol: " << sideOrient << std::endl;
        }
        int mat_id = neigel->MaterialId();
        int sideOrient = neigel->NormalOrientation(neig.Side());
        if (mat_id == 40 || mat_id == 15) {
          TPZCompEl *cel = neigel->Reference();
          TPZMultiphysicsElement *cmulint = dynamic_cast<TPZMultiphysicsElement *>(cel);
          TPZCompEl *hdivBound = cmulint->Element(0);
          TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(hdivBound);
          int sideOrient = hdivbound->GetSideOrient(neig.Side());
          std::cout << "SideOrientSide1: " << sideOrient << std::endl;
        }
        std::cout << "SideOrientSide: " << sideOrient << std::endl;
        std::cout << "PrevId: " << mat_id << std::endl;
        if (mat_id == 30 || mat_id == 35) {
          TPZCompEl *cel = neigel->Reference();
          TPZMultiphysicsInterfaceElement *cmulint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);

          TPZMaterial *mat = cmulint->Material();
          //                TPZBndCondT<STATE> *bc = dynamic_cast<TPZBndCondT<STATE> *> (mat);
          TPZLagrangeMultiplier<STATE> *lag = dynamic_cast<TPZLagrangeMultiplier<STATE> *>(mat);
          std::cout << "Mult= " << lag->Multiplier() << std::endl;
          int nconects = cmulint->NConnects();
          for (int icon = 0; icon < nconects; icon++) {
            int lagrange = cmulint->Connect(icon).LagrangeMultiplier();
            std::cout << "LagrangeMult= " << lagrange << std::endl;
          }

          int ok = 10;
        }
        neig = neig.Neighbour();
      }
      count++;
    }
  }
  //    }
}
