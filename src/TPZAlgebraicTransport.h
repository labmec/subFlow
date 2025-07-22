//
//  AlgebraicTransport.hpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#ifndef AlgebraicTransport_h
#define AlgebraicTransport_h
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "TMRSDataTransfer.h"

#include "TPZMultiphysicsCompMesh.h"

class TPZAlgebraicTransport {
    
public:
    
    struct TInterfaceDataTransport
    {
        int64_t fMatid;
        std::vector<int64_t> fcelindex;
        //vector for each multiplying coefficient of fluxes
        std::vector<std::vector<REAL> > fCoefficientsFlux;
        // left right volume index in the AlgebraicTransport data structure
        std::vector<std::pair<int64_t, int64_t>> fLeftRightVolIndex;
        
        // left right volume index in the AlgebraicTransport data structure
        std::vector<std::pair<int64_t, int64_t>> fLeftRightGelIndex;
        //Integral of the flux functions associated with faces
        std::vector<REAL> fIntegralFluxFunctions;
        //  Integral of the flux
        std::vector<REAL>  fIntegralFlux;
        //Sign with respect to the natural orientation of the flux (+1 or -1)
        std::vector<REAL> fFluxSign;
        //Direction of the normal to the face
        // prencher do elemento interpol
        std::vector<std::tuple<REAL,REAL,REAL>> fNormalFaceDirection;
        
        //dado x calcula a permeabilidade... calcular a permeabilidade
        // transferir las densidades al problema mixto en transport to mixed
        
        TInterfaceDataTransport() : fMatid(0), fcelindex(0), fCoefficientsFlux(0), fIntegralFluxFunctions(0), fLeftRightVolIndex(0),fIntegralFlux(0),fFluxSign(0), fNormalFaceDirection(0), fLeftRightGelIndex(0) {
           
        }

        TInterfaceDataTransport(const TInterfaceDataTransport &copy) = default;
        TInterfaceDataTransport &operator=(const TInterfaceDataTransport &copy) = default;
        
        void Print(std::ostream &out);
    };
    
    // CELL DATA
    struct TCellData{
       
        TMRSDataTransfer *fsim_data;
        
        
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
      
        //Fractional flow and its derivative
        //un vetor para fracional y outro para a derivada.
        // center o cord x, y, z tres vetores.
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
        
        TCellData() : fsim_data(0), fEqNumber(0),fVolume(0), fVolumefactor(0),fMatId(0),fGeoIndex(0),fSaturation(0), fPressure(0), fSaturationLastState(0),  fDensityOil(0),fdDensityOildp(0), fDensityOilLastState(0),fDensityWater(0), fdDensityWaterdp(0),fDensityWaterLastState(0),fMixedDensity(0), flambda(0), fdlambdawdsw(0),fdlambdaodsw(0),fporosity(0),fKx(0),fKy(0),fKz(0), fWaterfractionalflow(0),fDerivativeWfractionalflow(0),fOilfractionalflow(0), fDerivativeOfractionalflow(0),fCenterCoordinate(0),
        fCompressibility(0),fViscosity(0),fReferencePressures(0),
        fReferenceDensity(0)
        {
            
        }
        TCellData(const TCellData &copy) = default;
        TCellData &operator=(const TCellData &copy) = default;

        void SetNumCells(int64_t ncells)
        {
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
        void SetDataTransfer(TMRSDataTransfer *simdata);
        void Print(std::ostream &out);
    };
    
    
    REAL fdt = 0.1;
    std::map<int, std::pair<int ,REAL>> fboundaryCMatVal;
    std::vector<REAL> fgravity;
    int fNFluxCoefficients;
    int inletmatid;
    int outletmatid;
    int interfaceid;
    bool fHasPropQ=false;
    REAL massOut = 0.0;
    REAL initialMass = 0.;
    
    //number of volumetric elements in the transport mesh
    int fNVolumesTransport = 0;
    // Cells data structure, one material at a time
    TCellData fCellsData;

    
    // Interface data structure, by material, element and side
    std::map<int, TInterfaceDataTransport> fInterfaceData;
public:
    
     /// Default constructor
    TPZAlgebraicTransport();
    
    /// Copy constructor
    TPZAlgebraicTransport(const TPZAlgebraicTransport & other);
    
    /// Assignement constructor
    const TPZAlgebraicTransport & operator=(const TPZAlgebraicTransport & other);
    

    /// Default desconstructor
    ~TPZAlgebraicTransport();
    
    void BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh);
    
    void Contribute(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef);
    void ContributeResidual(int index, TPZFMatrix<double> &ef);
    
    void ContributeInterface(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef, int interfaceId=100);
    void ContributeInterfaceIHU(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef);
    
    void ContributeInterfaceResidual(int index, TPZFMatrix<double> &ef, int interfaceID=100);
    void ContributeInterfaceIHUResidual(int index, TPZFMatrix<double> &ef);
    
    void ContributeBCInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef, int matid);
    void ContributeBCInletInterface(int index,TPZFMatrix<double> &ef, int InletID = -2);
    void ContributeBCOutletInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef, int OutId=-4);
    void ContributeBCOutletInterfaceResidual(int index, TPZFMatrix<double> &ef, int matId=4);

    
    // IHU auxiliary functions
    std::pair<REAL, std::pair<REAL, REAL>> f_star(std::pair<REAL, REAL> foL, std::pair<REAL, REAL> foR, std::pair<REAL, REAL> fwL, std::pair<REAL, REAL> fwR, REAL g_dot_n);
    std::pair<REAL, std::pair<REAL, REAL>> lambda_o_star(std::pair<REAL, REAL> lambda_L, std::pair<REAL, REAL> lambda_R, REAL g_dot_n, REAL rho_ratio);
    std::pair<REAL, std::pair<REAL, REAL>> lambda_w_star(std::pair<REAL, REAL> lambda_L, std::pair<REAL, REAL> lambda_R, REAL g_dot_n, REAL rho_ratio);
    
    void UpdateIntegralFlux(int matid);
    REAL FLuxIntegralbyID(int mat_id);
    REAL CalculateMass();
    REAL CalculateMassById(int matId);
    REAL CalculateMassById2(int matId);
    REAL CalculateMassByCoord();
    REAL CalculateMassByCoord2();
    std::pair<REAL, REAL> FLuxWaterOilIntegralbyID(int mat_id);
    void VerifyElementFLuxes();
    void VerifyConservation(int itime);
    void ZeroFluxes();
    void PrintFluxes();
};

#endif /* AlgebraicTransport_h */
