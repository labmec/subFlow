//
//  Created by Giovane Avancini on 10/01/26.
//  Copied from iMRS and is evolving to fit subFlow needs
//  PLEASE ADAPT THIS FILE
//

#pragma once

#include "TPZMultiphysicsCompMesh.h"
#include "TSFProblemData.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

class TSFAlgebraicTransport {

public:
  struct TInterfaceData {
    int64_t fMatid;
    std::vector<int64_t> fcelindex;
    // vector for each multiplying coefficient of fluxes
    std::vector<std::vector<REAL>> fCoefficientsFlux;
    // left right volume index in the AlgebraicTransport data structure
    std::vector<std::pair<int64_t, int64_t>> fLeftRightVolIndex;
    // left right gel index in the AlgebraicTransport data structure
    std::vector<std::pair<int64_t, int64_t>> fLeftRightGelIndex;
    // Integral of the flux functions associated with faces
    std::vector<REAL> fIntegralFluxFunctions;
    // Integral of the flux
    std::vector<REAL> fIntegralFlux;
    // Sign with respect to the natural orientation of the flux (+1 or -1)
    std::vector<REAL> fFluxSign;
    // Direction of the normal to the face
    std::vector<std::tuple<REAL, REAL, REAL>> fNormalFaceDirection;

    TInterfaceData() : fMatid(0), fcelindex(0), fCoefficientsFlux(0), fIntegralFluxFunctions(0), fLeftRightVolIndex(0), fIntegralFlux(0), fFluxSign(0), fNormalFaceDirection(0), fLeftRightGelIndex(0) {
    }

    TInterfaceData(const TInterfaceData &copy) = default;
    TInterfaceData &operator=(const TInterfaceData &copy) = default;

    void Print(std::ostream &out);
  };

  struct TCellData {

    TSFProblemData *fSimData;

    std::vector<REAL> fEqNumber;
    std::vector<REAL> fVolume;
    std::vector<REAL> fVolumefactor;
    std::vector<REAL> fMatId;
    std::vector<REAL> fGeoIndex;
    std::vector<REAL> fPressure;
    std::vector<REAL> fSaturation;
    std::vector<REAL> fSaturationLastState;
    std::vector<REAL> fDensityOil;
    std::vector<REAL> fdDensityOildp;
    std::vector<REAL> fDensityOilLastState;
    std::vector<REAL> fDensityWater;
    std::vector<REAL> fdDensityWaterdp;
    std::vector<REAL> fDensityWaterLastState;
    std::vector<REAL> fMixedDensity;
    std::vector<REAL> flambda;
    std::vector<REAL> fdlambdawdsw;
    std::vector<REAL> fdlambdaodsw;
    std::vector<REAL> fporosity;
    std::vector<REAL> fKx;
    std::vector<REAL> fKy;
    std::vector<REAL> fKz;

    // Fractional flow and its derivative
    // un vetor para fracional y outro para a derivada.
    //  center o cord x, y, z tres vetores.
    std::vector<REAL> fWaterfractionalflow;
    std::vector<REAL> fDerivativeWfractionalflow;
    std::vector<REAL> fOilfractionalflow;
    std::vector<REAL> fDerivativeOfractionalflow;

    //        std::vector<std::vector<REAL>> fOilfractionalflow;
    std::vector<std::vector<REAL>> fCenterCoordinate;

    // fCompressibility[0] = water, fCompressibility[1]=oil, fCompressibility[2]=gas etc.
    std::vector<REAL> fCompressibility;
    std::vector<REAL> fViscosity;
    std::vector<REAL> fReferencePressures;
    std::vector<REAL> fReferenceDensity;

    TCellData() : fSimData(0), fEqNumber(0), fVolume(0), fVolumefactor(0), fMatId(0), fGeoIndex(0), fSaturation(0), fPressure(0), fSaturationLastState(0), fDensityOil(0), fdDensityOildp(0), fDensityOilLastState(0), fDensityWater(0), fdDensityWaterdp(0), fDensityWaterLastState(0), fMixedDensity(0), flambda(0), fdlambdawdsw(0), fdlambdaodsw(0), fporosity(0), fKx(0), fKy(0), fKz(0), fWaterfractionalflow(0), fDerivativeWfractionalflow(0), fOilfractionalflow(0), fDerivativeOfractionalflow(0), fCenterCoordinate(0),
                  fCompressibility(0), fViscosity(0), fReferencePressures(0),
                  fReferenceDensity(0) {
    }
    TCellData(const TCellData &copy) = default;
    TCellData &operator=(const TCellData &copy) = default;

    void SetNumCells(int64_t ncells) {
      fVolume.resize(ncells);
      fVolumefactor.resize(ncells);
      fMatId.resize(ncells);
      fGeoIndex.resize(ncells);
      fEqNumber.resize(ncells);
      fSaturation.resize(ncells);
      fPressure.resize(ncells);
      fSaturationLastState.resize(ncells);
      fporosity.resize(ncells);
      fKx.resize(ncells);
      fKy.resize(ncells);
      fKz.resize(ncells);
      fPressure.resize(ncells);
      fDensityOil.resize(ncells);
      fdDensityOildp.resize(ncells);
      fDensityOilLastState.resize(ncells);
      fDensityWater.resize(ncells);
      fdDensityWaterdp.resize(ncells);
      fDensityWaterLastState.resize(ncells);
      fMixedDensity.resize(ncells);
      flambda.resize(ncells);
      fdlambdawdsw.resize(ncells);
      fdlambdaodsw.resize(ncells);
      fWaterfractionalflow.resize(ncells);
      fDerivativeWfractionalflow.resize(ncells);
      fOilfractionalflow.resize(ncells);
      fDerivativeOfractionalflow.resize(ncells);
      fCenterCoordinate.resize(ncells);
    }
    void UpdateSaturations(TPZFMatrix<STATE> &dsx);
    void UpdateSaturationsLastState(TPZFMatrix<STATE> &sw);
    void UpdateSaturationsTo(TPZFMatrix<STATE> &sw);
    void UpdateFractionalFlowsAndLambda(int krModel);
    void UpdateFractionalFlowsAndLambdaQuasiNewton();
    void UpdateDensities();
    void UpdateDensitiesLastState();
    void UpdateMixedDensity();
    void SetDataTransfer(TSFProblemData *simdata);
    void Print(std::ostream &out);
  };

  std::map<int, std::pair<int, REAL>> fboundaryCMatVal;
  int fNFluxCoefficients;
  REAL massOut = 0.0;
  REAL initialMass = 0.0;

  // number of volumetric elements in the transport mesh
  int fNVolumesTransport;
  // Cells data structure, one material at a time
  TCellData fCellsData;
  // Interface data structure, by material, element and side
  std::map<int, TInterfaceData> fInterfaceData;

public:
  /// Default constructor
  TSFAlgebraicTransport();

  /// Copy constructor
  TSFAlgebraicTransport(const TSFAlgebraicTransport &other);

  /// Assignement constructor
  const TSFAlgebraicTransport &operator=(const TSFAlgebraicTransport &other);

  /// Default desconstructor
  ~TSFAlgebraicTransport();

  void BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh);

  // PLEASE IMPLEMENT THE CONTRIBUTE METHODS
};