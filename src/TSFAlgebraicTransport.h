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
    std::vector<REAL> fMatId;
    std::vector<REAL> fGeoIndex;
    std::vector<REAL> fPressure;
    std::vector<REAL> fSaturation;
    std::vector<REAL> fSaturationLastState;
    std::vector<REAL> fDensityWater;
    std::vector<REAL> fDdensityWaterdp;
    std::vector<REAL> fDensityWaterLastState;
    std::vector<REAL> fVolumeFactorWater;
    std::vector<REAL> fDensityGas;
    std::vector<REAL> fDdensityGasdp;
    std::vector<REAL> fDensityGasLastState;
    std::vector<REAL> fVolumeFactorGas;
    std::vector<REAL> fMixedDensity;
    std::vector<REAL> fLambda;
    std::vector<REAL> fDlambdaWaterdsw;
    std::vector<REAL> fDlambdaGasdsw;
    std::vector<REAL> fPorosity;
    std::vector<REAL> fKappa;

    // Fractional flow and its derivative
    // un vetor para fracional y outro para a derivada.
    //  center o cord x, y, z tres vetores.
    std::vector<REAL> fWaterfractionalflow;
    std::vector<REAL> fDerivativeWfractionalflow;
    std::vector<REAL> fGasfractionalflow;
    std::vector<REAL> fDerivativeGfractionalflow;

    std::vector<std::vector<REAL>> fCenterCoordinate;

    // [0] = water, [1] = gas
    std::vector<REAL> fCompressibility;
    std::vector<REAL> fViscosity;
    std::vector<REAL> fReferencePressures;
    std::vector<REAL> fReferenceDensity;

    TCellData() : fSimData(0), fEqNumber(0), fVolume(0), fVolumeFactorWater(0), fMatId(0), fGeoIndex(0), fSaturation(0), fPressure(0), fSaturationLastState(0), fDensityGas(0), fDdensityGasdp(0), fDensityGasLastState(0), fDensityWater(0), fDdensityWaterdp(0), fDensityWaterLastState(0), fMixedDensity(0), fLambda(0), fDlambdaWaterdsw(0), fDlambdaGasdsw(0), fPorosity(0), fKappa(0), fWaterfractionalflow(0), fDerivativeWfractionalflow(0), fGasfractionalflow(0), fDerivativeGfractionalflow(0), fCenterCoordinate(0),
                  fCompressibility(0), fViscosity(0), fReferencePressures(0),
                  fReferenceDensity(0) {
    }
    TCellData(const TCellData &copy) = default;
    TCellData &operator=(const TCellData &copy) = default;

    void SetNumCells(int64_t ncells) {
      fVolume.resize(ncells);
      fVolumeFactorWater.resize(ncells);
      fVolumeFactorGas.resize(ncells);
      fMatId.resize(ncells);
      fGeoIndex.resize(ncells);
      fEqNumber.resize(ncells);
      fSaturation.resize(ncells);
      fPressure.resize(ncells);
      fSaturationLastState.resize(ncells);
      fPorosity.resize(ncells);
      fKappa.resize(ncells);
      fPressure.resize(ncells);
      fDensityGas.resize(ncells);
      fDdensityGasdp.resize(ncells);
      fDensityGasLastState.resize(ncells);
      fDensityWater.resize(ncells);
      fDdensityWaterdp.resize(ncells);
      fDensityWaterLastState.resize(ncells);
      fMixedDensity.resize(ncells);
      fLambda.resize(ncells);
      fDlambdaWaterdsw.resize(ncells);
      fDlambdaGasdsw.resize(ncells);
      fWaterfractionalflow.resize(ncells);
      fDerivativeWfractionalflow.resize(ncells);
      fGasfractionalflow.resize(ncells);
      fDerivativeGfractionalflow.resize(ncells);
      fCenterCoordinate.resize(ncells);
    }
    void UpdateSaturations(TPZFMatrix<STATE> &dsx);
    void UpdateSaturationsLastState(TPZFMatrix<STATE> &sw);
    void UpdateSaturationsTo(TPZFMatrix<STATE> &sw);
    void UpdateFractionalFlowsAndLambda(int krModel);
    void UpdateDensities();
    void UpdateDensitiesLastState();
    void UpdateMixedDensity();
    void SetProblemData(TSFProblemData *simdata);
    void Print(std::ostream &out);
  };

  std::map<int, std::pair<int, REAL>> fboundaryCMatVal;
  int fNFluxCoefficients;
  REAL fWaterMassOut = 0.0;
  REAL fWaterMassIn = 0.0;
  REAL fWaterMass = 0.0;
  REAL fGasMassOut = 0.0;
  REAL fGasMassIn = 0.0;
  REAL fGasMass = 0.0;
  REAL fInitialWaterMass = 0.0;
  REAL fInitialGasMass = 0.0;

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

  /// Check Mass Balance
  void CheckMassBalance(REAL time_step, std::ostream &out = std::cout);

  /// Update the integrated fluxes at the interfaces for a given material id
  void UpdateInterfacesIntegratedFlux(int matid);

  // PLEASE IMPLEMENT THE CONTRIBUTE METHODS
  void Contribute(int cellId, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);

  void ContributeResidual(int cellId, TPZFMatrix<REAL> &ef);
  
  void ContributeInterface(int interfaceId, int interfaceMatId, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);

  void ContributeInterfaceResidual(int interfaceId, int interfaceMatId, TPZFMatrix<REAL> &ef);
};