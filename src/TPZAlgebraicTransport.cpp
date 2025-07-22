//
//  AlgebraicTransport.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "TPZAlgebraicTransport.h"
#include "pzelmat.h"
//#include <Eigen/PardisoSupport>

/// Default constructor
TPZAlgebraicTransport::TPZAlgebraicTransport(){
    
    
}

/// Copy constructor
TPZAlgebraicTransport::TPZAlgebraicTransport(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
    fNFluxCoefficients =other.fNFluxCoefficients;
    inletmatid=other.inletmatid;
    outletmatid=other.outletmatid;
    interfaceid=other.interfaceid;
    fgravity=other.fgravity;
    fHasPropQ = other.fHasPropQ;
    fboundaryCMatVal =other.fboundaryCMatVal;
}

/// Assignement constructor
const TPZAlgebraicTransport & TPZAlgebraicTransport::operator=(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
    inletmatid=other.inletmatid;
    outletmatid=other.outletmatid;
    interfaceid=other.interfaceid;
    fgravity=other.fgravity;
    fHasPropQ = other.fHasPropQ;
    fboundaryCMatVal =other.fboundaryCMatVal;
    return *this;
}
TPZAlgebraicTransport::~TPZAlgebraicTransport(){
    
}

void TPZAlgebraicTransport::BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh)
{

}

void TPZAlgebraicTransport::Contribute(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef){
    
    REAL sat =fCellsData.fSaturation[index];
    REAL satLast = fCellsData.fSaturationLastState[index];
    REAL densityWater = fCellsData.fDensityWater[index];
    REAL densityWaterLastState = fCellsData.fDensityWaterLastState[index];
    REAL phi = fCellsData.fporosity[index];
#ifdef PZDEBUG
    if(std::abs(phi) < 1e-12) DebugStop();
#endif
    ef(0) = fCellsData.fVolume[index]*phi*(sat*densityWater-satLast*densityWaterLastState);
    ek(0,0) = fCellsData.fVolume[index]*phi*densityWater;
}

void TPZAlgebraicTransport::ContributeResidual(int index, TPZFMatrix<double> &ef){
    
    REAL sat =fCellsData.fSaturation[index];
    REAL satLast = fCellsData.fSaturationLastState[index];
    REAL densityWater = fCellsData.fDensityWater[index];
    REAL densityWaterLastState = fCellsData.fDensityWaterLastState[index];
    REAL phi = fCellsData.fporosity[index];
    ef(0) = fCellsData.fVolume[index]*phi*(sat*densityWater-satLast*densityWaterLastState);
}

void TPZAlgebraicTransport::ContributeInterface(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef, int interfaceId){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceId].fLeftRightVolIndex[index];
#ifdef PZDEBUG
    const REAL phiL = fCellsData.fporosity[lr_index.first];
    const REAL phiR = fCellsData.fporosity[lr_index.second];
    const REAL vfL = fCellsData.fVolumefactor[lr_index.first];
    const REAL vfR = fCellsData.fVolumefactor[lr_index.second];
    int matIdL = fCellsData.fporosity[lr_index.first];
    if(std::abs(phiL) < 1e-12) DebugStop();
    if(std::abs(phiR) < 1e-12) DebugStop();
    if(std::abs(vfL) < 1e-12) DebugStop();
    if(std::abs(vfR) < 1e-12) DebugStop();
#endif
    
    REAL fluxint  = 1.0*fInterfaceData[interfaceId].fIntegralFlux[index];
    REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
    REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];
    REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    REAL dfwSw_R = fCellsData.fDerivativeWfractionalflow[lr_index.second];
    REAL beta =0.0;
 //upwind
    if (fluxint>0.0) {
        beta = 1.0;
    }
    
    ef(0) = +1.0*(beta*fw_L + (1-beta)*fw_R)*fluxint* fdt;
    ef(1) = -1.0*(beta*fw_L  + (1-beta)*fw_R)*fluxint* fdt;
    
    ek(0,0) = +1.0*dfwSw_L  * beta * fluxint*fdt;
    ek(0,1) = +1.0*dfwSw_R * (1-beta) * fluxint*fdt;
    ek(1,0) = -1.0*dfwSw_L* beta * fluxint*fdt;
    ek(1,1) = -1.0*dfwSw_R * (1-beta)*fluxint* fdt;
   
//    Gravity fluxes contribution
   ContributeInterfaceIHU(index, ek, ef);
}

void TPZAlgebraicTransport::ContributeInterfaceResidual(int index, TPZFMatrix<double> &ef, int interfaceID){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceID].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[interfaceID].fIntegralFlux[index];
   
    REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
    REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];

    REAL beta =0.0;
//upwind
    if (fluxint>0.0) {
        beta = 1.0;
    }
    
    ef(0) = +1.0*(beta*fw_L + (1-beta)*fw_R)*fluxint* fdt;
    ef(1) = -1.0*(beta*fw_L  + (1-beta)*fw_R)*fluxint* fdt;
    
// Gravity fluxes contribution
   ContributeInterfaceIHUResidual(index, ef);
    
#ifdef PZDEBUG
    if(std::isnan(Norm(ef)))
    {
        std::cout << __PRETTY_FUNCTION__ << " nan" << std::endl;
    }
#endif
}

void TPZAlgebraicTransport::ContributeInterfaceIHU(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceid].fLeftRightVolIndex[index];
    std::tuple<REAL, REAL, REAL> normal = fInterfaceData[interfaceid].fNormalFaceDirection[index];
    
    std::vector<REAL> n(3,0.0);
    n[0] = std::get<0>(normal);
    n[1] = std::get<1>(normal);
    n[2] = std::get<2>(normal);
    
    REAL g_dot_n = n[0]*fgravity[0]+n[1]*fgravity[1]+n[2]*fgravity[2];
    
    REAL lambdaL = fCellsData.flambda[lr_index.first];
    REAL lambdaR = fCellsData.flambda[lr_index.second];
    
    std::pair<REAL, REAL> fwL = {fCellsData.fWaterfractionalflow[lr_index.first], fCellsData.fDerivativeWfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> fwR = {fCellsData.fWaterfractionalflow[lr_index.second],
        fCellsData.fDerivativeWfractionalflow[lr_index.second]};
    std::pair<REAL, REAL> foL = {fCellsData.fOilfractionalflow[lr_index.first],
            fCellsData.fDerivativeOfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> foR = {fCellsData.fOilfractionalflow[lr_index.second],
                fCellsData.fDerivativeOfractionalflow[lr_index.second]};
    
    std::pair<REAL, REAL> lambda_wL = {fwL.first * lambdaL, fCellsData.fdlambdawdsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_wR = {fwR.first * lambdaR, fCellsData.fdlambdawdsw[lr_index.second]};
    std::pair<REAL, REAL> lambda_oL = {foL.first * lambdaL, fCellsData.fdlambdaodsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_oR = {foR.first * lambdaR, fCellsData.fdlambdaodsw[lr_index.second]};
    
    REAL rho_wL = fCellsData.fDensityWater[lr_index.first];
    REAL rho_wR = fCellsData.fDensityWater[lr_index.second];
    REAL rho_oL = fCellsData.fDensityOil[lr_index.first];
    REAL rho_oR = fCellsData.fDensityOil[lr_index.second];
    
    // The upwinding logic should be the same for each function
    std::pair<REAL, std::pair<REAL, REAL>> fstarL = f_star(foL, foR, fwL, fwR, g_dot_n);
    std::pair<REAL, std::pair<REAL, REAL>> fstarR = f_star(foR, foL, fwR, fwL, -g_dot_n);
    
//    REAL rho_ratio_wL = ((rho_wL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_wR = ((rho_wR - rho_oR)/(rho_wR - rho_oR));
//    REAL rho_ratio_oL = ((rho_oL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_oR = ((rho_oR - rho_oR)/(rho_wR - rho_oR));
    REAL rho_ratio_wL = 1.;
    REAL rho_ratio_wR = 1.;
    REAL rho_ratio_oL = 0.;
    REAL rho_ratio_oR = 0.;

    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starL = lambda_w_star(lambda_wL, lambda_wR, g_dot_n, rho_ratio_wL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starR = lambda_w_star(lambda_wR, lambda_wL, -g_dot_n, rho_ratio_wR);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starL = lambda_o_star(lambda_oL, lambda_oR, g_dot_n, rho_ratio_oL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starR = lambda_o_star(lambda_oR, lambda_oL, -g_dot_n, rho_ratio_oR);
    
    // Harmonic permeability mean
    REAL Kx_L =  fCellsData.fKx[lr_index.first];
    REAL Ky_L =  fCellsData.fKy[lr_index.first];
    REAL Kz_L =  fCellsData.fKz[lr_index.first];
    
    REAL Kx_R =  fCellsData.fKx[lr_index.first];
    REAL Ky_R =  fCellsData.fKy[lr_index.first];
    REAL Kz_R =  fCellsData.fKz[lr_index.first];
    
    REAL K_x = 2.0*(Kx_L * Kx_R)/(Kx_L + Kx_R);
    REAL K_y = 2.0*(Ky_L * Ky_R)/(Ky_L + Ky_R);
    REAL K_z = 2.0*(Kz_L * Kz_R)/(Kz_L + Kz_R);
    
    // Beacuse we assume diagonal abs. perm tensor
    REAL K_times_g_dot_n = (K_x * n[0] * fgravity[0] + K_y * n[1] * fgravity[1] + K_z * n[2] * fgravity[2]);
    
    REAL res1 = fstarL.first * (lamba_w_starL.first + lamba_o_starL.first) * K_times_g_dot_n * (rho_wL - rho_oL) * fdt;
    REAL res2 = fstarR.first * (lamba_w_starR.first + lamba_o_starR.first) * K_times_g_dot_n * (rho_wR - rho_oR) * fdt;
    ef(0) += res1;
    ef(1) -= res2;
    
    REAL dGLdSL = fstarL.second.first * (lamba_w_starL.first + lamba_o_starL.first) + fstarL.first * (lamba_w_starL.second.first + lamba_o_starL.second.first);
    REAL dGLdSR = fstarL.second.second * (lamba_w_starL.first + lamba_o_starL.first) + fstarL.first * (lamba_w_starL.second.second + lamba_o_starL.second.second);
    REAL dGRdSL = fstarR.second.first * (lamba_w_starR.first + lamba_o_starR.first) + fstarR.first * (lamba_w_starR.second.first + lamba_o_starR.second.first);
    REAL dGRdSR = fstarR.second.second * (lamba_w_starR.first + lamba_o_starR.first) + fstarR.first * (lamba_w_starR.second.second + lamba_o_starR.second.second);
    
    ek(0,0) += dGLdSL * K_times_g_dot_n * (rho_wL - rho_oL) * fdt;
    ek(0,1) += dGLdSR * K_times_g_dot_n * (rho_wL - rho_oL) * fdt;
    
    ek(1,1) -= dGRdSL * K_times_g_dot_n * (rho_wR - rho_oR) * fdt;
    ek(1,0) -= dGRdSR * K_times_g_dot_n * (rho_wR - rho_oR) * fdt;
}

void TPZAlgebraicTransport::ContributeInterfaceIHUResidual(int index, TPZFMatrix<double> &ef){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceid].fLeftRightVolIndex[index];
    std::tuple<REAL, REAL, REAL> normal = fInterfaceData[interfaceid].fNormalFaceDirection[index];
    
    std::vector<REAL> n(3,0.0);
    n[0] = std::get<0>(normal);
    n[1] = std::get<1>(normal);
    n[2] = std::get<2>(normal);
    
    REAL g_dot_n = n[0]*fgravity[0]+n[1]*fgravity[1]+n[2]*fgravity[2];
    
    REAL lambdaL = fCellsData.flambda[lr_index.first];
    REAL lambdaR = fCellsData.flambda[lr_index.second];
    
    std::pair<REAL, REAL> fwL = {fCellsData.fWaterfractionalflow[lr_index.first], fCellsData.fDerivativeWfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> fwR = {fCellsData.fWaterfractionalflow[lr_index.second],
        fCellsData.fDerivativeWfractionalflow[lr_index.second]};
    std::pair<REAL, REAL> foL = {fCellsData.fOilfractionalflow[lr_index.first],
            fCellsData.fDerivativeOfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> foR = {fCellsData.fOilfractionalflow[lr_index.second],
                fCellsData.fDerivativeOfractionalflow[lr_index.second]};
    
    std::pair<REAL, REAL> lambda_wL = {fwL.first * lambdaL, fCellsData.fdlambdawdsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_wR = {fwR.first * lambdaR, fCellsData.fdlambdawdsw[lr_index.second]};
    std::pair<REAL, REAL> lambda_oL = {foL.first * lambdaL, fCellsData.fdlambdaodsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_oR = {foR.first * lambdaR, fCellsData.fdlambdaodsw[lr_index.second]};
    
    REAL rho_wL = fCellsData.fDensityWater[lr_index.first];
    REAL rho_wR = fCellsData.fDensityWater[lr_index.second];
    REAL rho_oL = fCellsData.fDensityOil[lr_index.first];
    REAL rho_oR = fCellsData.fDensityOil[lr_index.second];
    
    // The upwinding logic should be the same for each function
    std::pair<REAL, std::pair<REAL, REAL>> fstarL = f_star(foL, foR, fwL, fwR, g_dot_n);
    std::pair<REAL, std::pair<REAL, REAL>> fstarR = f_star(foR, foL, fwR, fwL, -g_dot_n);
    
//    REAL rho_ratio_wL = ((rho_wL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_wR = ((rho_wR - rho_oR)/(rho_wR - rho_oR));
//    REAL rho_ratio_oL = ((rho_oL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_oR = ((rho_oR - rho_oR)/(rho_wR - rho_oR));
    REAL rho_ratio_wL = 1.;
    REAL rho_ratio_wR = 1.;
    REAL rho_ratio_oL = 0.;
    REAL rho_ratio_oR = 0.;

    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starL = lambda_w_star(lambda_wL, lambda_wR, g_dot_n, rho_ratio_wL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starR = lambda_w_star(lambda_wR, lambda_wL, -g_dot_n, rho_ratio_wR);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starL = lambda_o_star(lambda_oL, lambda_oR, g_dot_n, rho_ratio_oL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starR = lambda_o_star(lambda_oR, lambda_oL, -g_dot_n, rho_ratio_oR);
    
    // Harmonic permeability mean
    REAL Kx_L =  fCellsData.fKx[lr_index.first];
    REAL Ky_L =  fCellsData.fKy[lr_index.first];
    REAL Kz_L =  fCellsData.fKz[lr_index.first];
    
    REAL Kx_R =  fCellsData.fKx[lr_index.first];
    REAL Ky_R =  fCellsData.fKy[lr_index.first];
    REAL Kz_R =  fCellsData.fKz[lr_index.first];
    
    REAL K_x = 2.0*(Kx_L * Kx_R)/(Kx_L + Kx_R);
    REAL K_y = 2.0*(Ky_L * Ky_R)/(Ky_L + Ky_R);
    REAL K_z = 2.0*(Kz_L * Kz_R)/(Kz_L + Kz_R);
    
    // Beacuse we assume diagonal abs. perm tensor
    REAL K_times_g_dot_n = (K_x*n[0]*fgravity[0]+K_y*n[1]*fgravity[1]+K_z*n[2]*fgravity[2]);
    
    REAL res1 = fstarL.first * (lamba_w_starL.first + lamba_o_starL.first) * K_times_g_dot_n * (rho_wL - rho_oL) * fdt;
    REAL res2 = fstarR.first * (lamba_w_starR.first + lamba_o_starR.first) * K_times_g_dot_n * (rho_wR - rho_oR) * fdt;
    ef(0) += res1;
    ef(1) -= res2;

}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::f_star(std::pair<REAL, REAL> foL, std::pair<REAL, REAL> foR, std::pair<REAL, REAL> fwL, std::pair<REAL, REAL> fwR, REAL g_dot_n){
    REAL fstar, dfstardsL, dfstardsR;
    if( g_dot_n < 0.0){
        fstar = foL.first * fwR.first;
        dfstardsL = foL.second * fwR.first;
        dfstardsR = foL.first * fwR.second;
    }else{
        fstar = foR.first * fwL.first;
        dfstardsL = foR.first * fwL.second;
        dfstardsR = foR.second * fwL.first;
    }
    return std::make_pair(fstar, std::make_pair(dfstardsL, dfstardsR));
}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::lambda_o_star(std::pair<REAL, REAL> lambda_L, std::pair<REAL, REAL> lambda_R, REAL g_dot_n, REAL rho_ratio){
    REAL lambda_star, dlambda_starL, dlambda_starR;
    if( g_dot_n < 0.0){
        lambda_star = rho_ratio * lambda_L.first  + (1-rho_ratio) * lambda_R.first;
        dlambda_starL = rho_ratio * lambda_L.second;
        dlambda_starR = (1-rho_ratio) * lambda_R.second;
    }else{
        lambda_star = rho_ratio * lambda_R.first  + (1-rho_ratio) * lambda_L.first;
        dlambda_starL = (1-rho_ratio) * lambda_L.second;
        dlambda_starR = rho_ratio * lambda_R.second;
    }
    return std::make_pair(lambda_star, std::make_pair(dlambda_starL, dlambda_starR));
}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::lambda_w_star(std::pair<REAL, REAL> lambda_L, std::pair<REAL, REAL> lambda_R, REAL g_dot_n, REAL rho_ratio){
    REAL lambda_star, dlambda_starL, dlambda_starR;
    if( g_dot_n > 0.0){
        lambda_star = rho_ratio * lambda_L.first  + (1-rho_ratio) * lambda_R.first;
        dlambda_starL = rho_ratio * lambda_L.second;
        dlambda_starR = (1-rho_ratio) * lambda_R.second;
    }else{
        lambda_star = rho_ratio * lambda_R.first  + (1-rho_ratio) * lambda_L.first;
        dlambda_starL = (1-rho_ratio) * lambda_L.second;
        dlambda_starR = rho_ratio * lambda_R.second;
    }
    return std::make_pair(lambda_star, std::make_pair(dlambda_starL, dlambda_starR));
}

void TPZAlgebraicTransport::ContributeBCInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef, int matid){

    REAL fluxint = fInterfaceData[matid].fIntegralFlux[index];
    REAL zerotol = 1.e-10;
    REAL rhoWRef = fCellsData.fsim_data->mTFluidProperties.mWaterDensityRef;
    REAL rhoORef = fCellsData.fsim_data->mTFluidProperties.mOilDensityRef;
    int type = fboundaryCMatVal[matid].first;
    const bool noflux = (fabs(fluxint) < zerotol && type == 1) ? true : false;
    if (fluxint < 0.0 && !noflux){ //inlet
        REAL s_inlet = fboundaryCMatVal[matid].second; //external saturation
        //What is correct: to use the inlet saturation or the fractional flow at the inlet?
        int krmodel = fCellsData.fsim_data->mTPetroPhysics.mKrModel;
        // krmodel = 0; //using the inlet saturation direcly
        auto fwf = fCellsData.fsim_data->mTPetroPhysics.mFw[krmodel];
        REAL fw_inlet = std::get<0>(fwf(s_inlet, rhoWRef, rhoORef));
        ef(0,0) = fw_inlet*fluxint*fdt;
    }
    else if (!noflux){ //outlet

        std::pair<int64_t, int64_t> lr_index = fInterfaceData[matid].fLeftRightVolIndex[index];
        REAL fw_L= fCellsData.fWaterfractionalflow[lr_index.first];
        REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
        ef(0,0) = fw_L*fluxint* fdt;
        ek(0,0) = dfwSw_L*fluxint* fdt;
    }
}

void TPZAlgebraicTransport::ContributeBCInletInterface(int index, TPZFMatrix<double> &ef, int inId){
   
    REAL s_inlet = 1.0;// fboundaryCMatVal[inId].second;;
    REAL fluxint  = 1.0*fInterfaceData[inId].fIntegralFlux[index];
    ef(0,0) = 1.0*s_inlet*fluxint* fdt;
}
void TPZAlgebraicTransport::ContributeBCOutletInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef, int outID){
  
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[outID].fLeftRightVolIndex[index];
    REAL fluxint  = 1.0*fInterfaceData[outID].fIntegralFlux[index];
    REAL fw_L= fCellsData.fWaterfractionalflow[lr_index.first];
    REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    ef(0,0) = fw_L*fluxint* fdt;
    ek(0,0) = dfwSw_L*fluxint* fdt;
}

void TPZAlgebraicTransport::ContributeBCOutletInterfaceResidual(int index, TPZFMatrix<double> &ef, int outId){
  
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[outId].fLeftRightVolIndex[index];
    REAL fluxint  = 1.0*fInterfaceData[outId].fIntegralFlux[index];
    REAL fw_L= fCellsData.fWaterfractionalflow[lr_index.first];
    ef(0,0) = fw_L*fluxint* fdt;
}

void TPZAlgebraicTransport::TCellData::Print(std::ostream &out){
    
    int nels = this->fVolume.size();
    for (int iel =0; iel< nels; iel++) {
//        out<<"Material_Id: "<<this->fMatId<<std::endl;
        out<<"Center Cord: ";
        for (int ic=0 ; ic<fCenterCoordinate[iel].size(); ic++) {
            out<<fCenterCoordinate[iel][ic]<<" ";
        }
        out<<std::endl;
        out<<"Volume : "<<this->fVolume[iel]<<std::endl;
        out<<"Oil density: "<<this->fDensityOil[iel]<<std::endl;
        out<<"Water density: "<<this->fDensityWater[iel]<<std::endl;
        out<<"Pressure: "<<this->fPressure[iel]<<std::endl;
        out<<"Saturation: "<<this->fSaturation[iel]<<std::endl;
        out<<"WaterFracFlow: "<<this->fWaterfractionalflow[iel]<<std::endl;
        out<<"dWaterFracFlowdSw: "<<this->fDerivativeWfractionalflow[iel]<<std::endl;
        out<<"OilFracFlow: "<<this->fOilfractionalflow[iel]<<std::endl;
        out<<"dOilFracFlowdSw: "<<this->fDerivativeOfractionalflow[iel]<<std::endl;
        out<<"Lambda: "<<this->flambda[iel]<<std::endl;
        out<<std::endl;
    }

}

void TPZAlgebraicTransport::UpdateIntegralFlux(int matid){
    
    if(fInterfaceData.find(matid) == fInterfaceData.end()) return;
    int nels = fInterfaceData[matid].fCoefficientsFlux.size();
    if(nels == 0) return;
    std::vector<REAL> val =fInterfaceData[matid].fCoefficientsFlux[0];
    fInterfaceData[matid].fIntegralFlux=val;
    for (int i = 1; i<nels; i++) {
        int np =fInterfaceData[matid].fCoefficientsFlux[i].size();
        
        for (int index = 0; index<np; index++) {
            REAL val = fInterfaceData[matid].fCoefficientsFlux[i][index];
            fInterfaceData[matid].fIntegralFlux[index] +=val;
//            std::cout<<"FluxInte: "<< val<<std::endl;
        }
      
    }
}
void TPZAlgebraicTransport::TInterfaceDataTransport::Print(std::ostream &out){
    
    int ninterfaces = this->fFluxSign.size();
    out<<"Material_ID: "<<this->fMatid<<std::endl;
    for (int iinter = 0 ; iinter <ninterfaces; iinter++) {
        out<<"Left_Index: "<<this->fLeftRightVolIndex[iinter].first<<std::endl;
        out<<"Right_Index: "<<this->fLeftRightVolIndex[iinter].second<<std::endl;
        out << "fCoefficientsFlux :";
        for (int icoe=0; icoe<fCoefficientsFlux.size(); icoe++) {
            out<<fCoefficientsFlux[icoe][iinter]<<" ";
        }
        out<<std::endl;
        out<<"IntegralFluxFunctions: "<<fIntegralFluxFunctions[iinter]<<std::endl;
        out<<"IntegralFlux: "<<fIntegralFlux[iinter]<<std::endl;
        out<<"fFluxSign: "<<fFluxSign[iinter]<<std::endl;
        out<<"fNormalDirection: ";
        out<<std::get<0>(fNormalFaceDirection[iinter])<<" ";
        out<<std::get<1>(fNormalFaceDirection[iinter])<<" ";
        out<<std::get<2>(fNormalFaceDirection[iinter])<<" ";
        out<<std::endl;
    }
   

}
void TPZAlgebraicTransport::TCellData::SetDataTransfer(TMRSDataTransfer *simdata){
    fsim_data = simdata;
}
void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambda(int krModel){
    
    auto labdaWf = fsim_data->mTPetroPhysics.mLambdaW[krModel];
    auto labdaOf = fsim_data->mTPetroPhysics.mLambdaO[krModel];
    auto lambdaTotalf = fsim_data->mTPetroPhysics.mLambdaTotal[krModel];
    auto fwf = fsim_data->mTPetroPhysics.mFw[krModel];
    auto fof = fsim_data->mTPetroPhysics.mFo[krModel];
    
    int nvols = this->fVolume.size();
    for (int ivol =0 ; ivol< nvols; ivol++) {

        REAL sw = this->fSaturation[ivol];
        REAL rhow = this->fDensityWater[ivol];
        REAL rhoo = this->fDensityOil[ivol];
        auto fwfvalderiv = fwf(sw, rhow, rhoo);
        auto fovalderiv  = fof(sw, rhow, rhoo);
        auto lambdaWvalderiv = labdaWf(sw, rhow);
        auto lambdaOvalderiv = labdaOf(sw, rhoo);
        auto lambdaTotalvalderiv = lambdaTotalf(sw, rhow, rhoo);
        
        this->fWaterfractionalflow[ivol] = std::get<0>(fwfvalderiv);
        this->fDerivativeWfractionalflow[ivol] = std::get<1>(fwfvalderiv);
        this->fOilfractionalflow[ivol] = std::get<0>(fovalderiv);
        this->fDerivativeOfractionalflow[ivol] = std::get<1>(fovalderiv);
        this->flambda[ivol] = std::get<0>(lambdaTotalvalderiv);
        this->fdlambdawdsw[ivol] =std::get<1>(lambdaWvalderiv);
        this->fdlambdaodsw[ivol] = std::get<1>(lambdaOvalderiv);
    
    }
    
}

void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambdaQuasiNewton(){
    
    //Essa funcao tá estranha, lembrar de checar
        int krModel = fsim_data->mTPetroPhysics.mKrModel;
        int nvols = this->fVolume.size();
        auto labdaWf = fsim_data->mTPetroPhysics.mLambdaW[krModel];
        auto labdaOf = fsim_data->mTPetroPhysics.mLambdaO[krModel];
        auto lambdaTotalf = fsim_data->mTPetroPhysics.mLambdaTotal[krModel];
        auto fwf = fsim_data->mTPetroPhysics.mFw[krModel];
        auto fof = fsim_data->mTPetroPhysics.mFo[krModel];
    
        auto fwLinearf =fsim_data->mTPetroPhysics.mFo[0]; //linear krmodel
    
        for (int ivol =0 ; ivol< nvols; ivol++) {
            REAL sw = this->fSaturation[ivol];
            REAL rhow = this->fDensityWater[ivol];
            REAL rhoo = this->fDensityOil[ivol];
            auto fwfvalderiv = fwf(sw, rhow, rhoo);
            auto fwflinearvalderiv = fwLinearf(sw, rhow, rhoo);
            auto fovalderiv  = fof(sw, rhow, rhoo);
            auto lambdaWvalderiv = labdaWf(sw, rhow);
            auto lambdaOvalderiv = labdaOf(sw, rhoo);
            auto lambdaTotalvalderiv = lambdaTotalf(sw, rhow, rhoo);
            this->fWaterfractionalflow[ivol] = std::get<0>(fwfvalderiv);
            this->fOilfractionalflow[ivol] = std::get<0>(fovalderiv);
            this->fDerivativeOfractionalflow[ivol] = std::get<1>(fovalderiv);
            this->flambda[ivol] = std::get<0>(lambdaTotalvalderiv);
            this->fdlambdawdsw[ivol] =std::get<1>(lambdaWvalderiv);
            this->fdlambdaodsw[ivol] = std::get<1>(lambdaOvalderiv);
            
            fsim_data->mTPetroPhysics.CreateLinearKrModel();
            
            this->fDerivativeWfractionalflow[ivol] =std::get<1>(fwflinearvalderiv);
        }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturations(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        fSaturation[icell] = sw(eq_number);
    }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturationsTo(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        sw(eq_number) = fSaturation[icell];
    }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturationsLastState(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        fSaturationLastState[icell] = sw(eq_number);
    }
}
void TPZAlgebraicTransport::TCellData::UpdateMixedDensity(){
    int ncells = fVolume.size();
    for (int i =0; i< ncells; i++) {
        REAL mixedDen = (fWaterfractionalflow[i]*fDensityWater[i])+fOilfractionalflow[i]*fDensityOil[i];
        fMixedDensity[i] = mixedDen;
    }
}
void TPZAlgebraicTransport::TCellData::UpdateDensities(){
    int ncells = fVolume.size();
    auto fWaterDensityF = fsim_data->mTFluidProperties.mWaterDensityF;
    auto fOilDensityF = fsim_data->mTFluidProperties.mOilDensityF;
    if(fWaterDensityF){
    for (int icell = 0; icell<ncells; icell++) {
        REAL pressure = fPressure[icell];
        auto densityWvalderiv = fWaterDensityF(pressure);
        auto densityOvalderiv = fOilDensityF(pressure);
        #ifdef PZDEBUG
        if (std::get<0>(densityWvalderiv) < 0.0 || std::get<0>(densityOvalderiv) < 0.0) {
            // DebugStop();
        }
        #endif
        fDensityWater[icell] = std::get<0>(densityWvalderiv);
        fdDensityWaterdp[icell]= std::get<1>(densityWvalderiv);
        fDensityOil[icell] = std::get<0>(densityOvalderiv);
        fdDensityOildp[icell] = std::get<1>(densityOvalderiv);
    }
        
    }
    else{
        for (int icell = 0; icell<ncells; icell++) {
            fDensityWater[icell] = fReferenceDensity[0];
            fDensityOil[icell] = fReferenceDensity[1];
        }
    }
    
}
void TPZAlgebraicTransport::TCellData::UpdateDensitiesLastState(){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        fDensityWaterLastState[icell] = fDensityWater[icell];
        fDensityOilLastState[icell] = fDensityOil[icell];
    }
    
}



REAL TPZAlgebraicTransport::CalculateMass(){
    int ncells = fCellsData.fVolume.size();
    REAL intMass = 0.0;
    for (int icel = 0; icel < ncells; icel++) {
        REAL sat = fCellsData.fSaturation[icel];
        REAL phi = fCellsData.fporosity[icel];
        REAL vol = fCellsData.fVolume[icel];
        intMass += sat*phi*vol;
    }
    return intMass;
}

REAL TPZAlgebraicTransport::CalculateMassById(int matId){
    int ncells = fCellsData.fVolume.size();
    REAL intMass = 0.0;
    REAL volfrac= 0.0;
    
    for (int icel = 0; icel < ncells; icel++) {
        int celId = fCellsData.fMatId[icel];
        if(celId==matId){
            REAL sat = fCellsData.fSaturation[icel];
            REAL phi = fCellsData.fporosity[icel];
            REAL vol = fCellsData.fVolume[icel];
            intMass += sat*phi*vol;
            volfrac += vol*phi;
        }
    }
    if (volfrac==0.0) {
        return  0.0;
    }
    return intMass/volfrac;
}

void TPZAlgebraicTransport::VerifyConservation(int itime){
    int ncells = fCellsData.fVolume.size();
    REAL intMass = 0.0;
    REAL intWaterMass = 0.0;
    REAL intOilMass = 0.0;
    
    for (int icel = 0; icel < ncells; icel++) {
        REAL satW = fCellsData.fSaturation[icel];
        REAL satO = 1.0 - satW; //assuming two phase flow
        REAL rhoW = fCellsData.fDensityWater[icel];
        REAL rhoO = fCellsData.fDensityOil[icel];
        REAL phi = fCellsData.fporosity[icel];
        REAL vol = fCellsData.fVolume[icel];
        intWaterMass += satW*phi*vol*rhoW;
        intOilMass += satO*phi*vol*rhoO;
    }

    intMass = intWaterMass + intOilMass;

    //Separating the boundary conditions in outlet and inlet here instead of using the inlet and outlet matid defined in TMRSSFIAnalysis::FillProperties()
    std::set<int> bc_matids;
    for (auto it = fboundaryCMatVal.begin(); it != fboundaryCMatVal.end(); it++) {
        bc_matids.insert(it->first);
    }
    
    REAL OutletWaterMass=0.0;
    REAL OutletOilMass=0.0;
    REAL InletWaterMass=0.0;
    REAL InletOilMass=0.0;
    REAL fluxIntegratedNoFlux=0.0;
    REAL zerotol = 1.e-10;
    for (auto it = bc_matids.begin(); it != bc_matids.end(); it++) {
        
        int matid = *it;
        int type = fboundaryCMatVal[matid].first;
        int nInterfaces = fInterfaceData[matid].fIntegralFlux.size();
        for (int i = 0; i < nInterfaces; i++) {
            REAL fluxint = fInterfaceData[matid].fIntegralFlux[i];
            if (fabs(fluxint) < zerotol && type == 1) { //noflux
                std::pair<int64_t, int64_t> left_right = fInterfaceData[matid].fLeftRightVolIndex[i];
                int64_t cell_id = left_right.first;
                REAL int_saturation = fCellsData.fSaturation[cell_id];
                fluxIntegratedNoFlux += int_saturation*fluxint*fdt;
            }
            else if (fluxint < 0.0) { //inlet
                REAL s_inlet = fboundaryCMatVal[matid].second;
                REAL rhoWRef = fCellsData.fsim_data->mTFluidProperties.mWaterDensityRef;
                REAL rhoORef = fCellsData.fsim_data->mTFluidProperties.mOilDensityRef;
                int krmodel = fCellsData.fsim_data->mTPetroPhysics.mKrModel;
                // krmodel = 0; //using the inlet saturation direcly
                auto fwf = fCellsData.fsim_data->mTPetroPhysics.mFw[krmodel];
                auto fof = fCellsData.fsim_data->mTPetroPhysics.mFo[krmodel];
                REAL fw_inlet = std::get<0>(fwf(s_inlet, rhoWRef, rhoORef));
                REAL fo_inlet = std::get<0>(fof(s_inlet, rhoWRef, rhoORef));
                InletWaterMass += fw_inlet*fluxint*fdt; //before, this was multiplied by itime, which makes no sense
                InletOilMass += fo_inlet*fluxint*fdt;
            }
            else { //outlet
                std::pair<int64_t, int64_t> left_right = fInterfaceData[matid].fLeftRightVolIndex[i];
                int64_t cell_id = left_right.first;
                REAL int_saturation = fCellsData.fSaturation[cell_id];
                REAL fwL = fCellsData.fWaterfractionalflow[cell_id];
                REAL foL = fCellsData.fOilfractionalflow[cell_id];
                REAL sum = fwL + foL;
                OutletWaterMass += fwL*fluxint*fdt;
                OutletOilMass += foL*fluxint*fdt;
            }
        }
    }

    REAL InletMass = InletWaterMass + InletOilMass;
    REAL OutletMass = OutletWaterMass + OutletOilMass;

    REAL massConservation = intMass + InletMass + OutletMass - initialMass;
    std::cout << "\nGlobal Conservation Diagnostics" << std::endl;
    std::cout << "---Inlet water mass: " << std::setprecision(14) << InletWaterMass << std::endl;
    std::cout << "---Inlet oil mass: " << std::setprecision(14) << InletOilMass << std::endl;
    std::cout << "---Inlet mass: " << InletMass << std::endl;
    std::cout << "---Outlet water mass: " << std::setprecision(14) << OutletWaterMass << std::endl;
    std::cout << "---Outlet oil mass: " << std::setprecision(14) << OutletOilMass << std::endl;
    std::cout << "---Outlet mass: " << OutletMass << std::endl;
    std::cout << "---Inlet - Outlet: " << InletMass + OutletMass << std::endl;
    if(fabs(fluxIntegratedNoFlux) > 1.e-10 ){
        std::cout << "---WARNING! Flux through no flux bc is significant: " << fluxIntegratedNoFlux << std::endl;
    }
    else{
        std::cout << "---NoFlux mass: " << fluxIntegratedNoFlux << std::endl;
    }    
    std::cout << "---System water mass: " << intWaterMass << std::endl;
    std::cout << "---System oil mass: " << intOilMass << std::endl;
    std::cout << "---System mass: " << intMass << std::endl;
    std::cout << "---Initial mass: " << initialMass << std::endl;
    std::cout << "---System mass - Initial mass: " << intMass - initialMass << std::endl;
    std::cout << "---Accumulated outlet mass: " << massOut << std::endl;

    if(std::abs(massConservation) < 1.0e-8 ){
        std::cout << "---Global mass conservation is ok! Total mass loss: " << std::setprecision(14) << massConservation << std::endl;
    }
    else{
        std::cout << "---WARNING! Global mass conservation NOT ok! Total mass loss: " << std::setprecision(14) << massConservation << std::endl;
    }
    massOut += OutletMass;
    initialMass = intMass; //initialMass now stands for the mass at the end of the previous time step
}





REAL TPZAlgebraicTransport::CalculateMassByCoord(){
    int ncells = fCellsData.fVolume.size();
    REAL intMass=0.0;
    int ncomputed = 0;
    std::vector<std::vector<REAL>>& centerCordvvec = fCellsData.fCenterCoordinate;
    for (int icel =0; icel<ncells; icel++) {
        std::vector<REAL> & centerCoord= centerCordvvec[icel];
        const REAL x =centerCoord[0], y=centerCoord[1], z=centerCoord[2];
//      if (x>0.5 && x<1.0 && y>0.0 && y<0.5 && z>0.0 && z<0.5) {
//      if (x>0.5 && x<0.75 && y>0.5 && y<0.75 && z>0.75 && z<=1.0) {
        if (x>0.75 && x<1.0 && y>0.75 && y<1.0 && z>0.5 && z<0.75) {
            REAL sat = fCellsData.fSaturation[icel];
            intMass +=sat;
            ncomputed++;
        }
    }
    return intMass/ncomputed;
}



REAL TPZAlgebraicTransport::CalculateMassById2(int matId){
    int ncells = fCellsData.fVolume.size();
    REAL intMass = 0.0;
    int nfracels=0;
    
    for (int icel = 0; icel < ncells; icel++) {
        int celId = fCellsData.fMatId[icel];
        if(celId==matId){
            REAL sat = fCellsData.fSaturation[icel];
            REAL phi = fCellsData.fporosity[icel];
            REAL vol = fCellsData.fVolume[icel];
            intMass += sat;
            nfracels++;
        }
    }
    if (nfracels==0) {
        return  0.0;
    }
    return intMass/nfracels;
}
std::pair<REAL, REAL> TPZAlgebraicTransport::FLuxWaterOilIntegralbyID(int mat_id){
    
    REAL WaterIntegral =0.0;
    REAL OilIntegral = 0.0;
    int ninter = fInterfaceData[mat_id].fIntegralFlux.size();
    for (int iface =0 ; iface < ninter; iface++) {
        int LeftElIndex = fInterfaceData[mat_id].fLeftRightVolIndex[iface].first;
        REAL fracFluxWater = fCellsData.fWaterfractionalflow[LeftElIndex];
        REAL fracFluxOil = fCellsData.fOilfractionalflow[LeftElIndex];
        REAL FluxInttegral = fInterfaceData[mat_id].fIntegralFlux[iface];
        WaterIntegral += FluxInttegral*fracFluxWater;
        OilIntegral += FluxInttegral*fracFluxOil;
        
    }
    return std::make_pair(WaterIntegral, OilIntegral);
}
REAL TPZAlgebraicTransport::FLuxIntegralbyID(int mat_id){
    
    REAL FluxInt =0.0;
    
    int ninter = fInterfaceData[mat_id].fIntegralFlux.size();
    for (int iface =0 ; iface < ninter; iface++) {
        REAL FluxInttegral = fInterfaceData[mat_id].fIntegralFlux[iface];
        FluxInt += FluxInttegral;
    }
    return FluxInt;
}
void TPZAlgebraicTransport::VerifyElementFLuxes(){
    
    int nels = fCellsData.fVolume.size();
    std::vector<int> nInterfacesByElement(nels);
    std::vector<int> numZeroFluxByElement(nels);
    std::vector<int> IndexGels(nels,-1);
    std::vector<REAL> SumFluxByElement(nels);
    for(auto& interdata: fInterfaceData){
     
        TInterfaceDataTransport &transport = interdata.second;
        int neles = transport.fLeftRightVolIndex.size();
        for(int iel = 0; iel<neles; iel++){
            int leftIndex = transport.fLeftRightVolIndex[iel].first;
            int rightIndex = transport.fLeftRightVolIndex[iel].second;
            if(IndexGels[leftIndex] == -1){
                IndexGels[leftIndex] = 1;
//                std::cout<<"IndexVec: "<<leftIndex<<" GelIndex: "<<transport.fLeftRightGelIndex[iel].first<<std::endl;
            }
            if(IndexGels[rightIndex] == -1){
                IndexGels[rightIndex] = 1;
//                std::cout<<"IndexVec: "<<rightIndex<<" GelIndex: "<<transport.fLeftRightGelIndex[iel].second<<std::endl;
            }
                
            REAL fluxInt = transport.fIntegralFlux[iel];
            if(IsZero(fluxInt)){
                if(leftIndex >= 0)
                    numZeroFluxByElement[leftIndex]++;
                if(rightIndex >= 0)
                    numZeroFluxByElement[rightIndex]++;
            }
            
            SumFluxByElement[leftIndex] += fluxInt;
            nInterfacesByElement[leftIndex]++;
            if(rightIndex < 0){
                continue;
            }
            nInterfacesByElement[rightIndex]++;
            SumFluxByElement[rightIndex] += -fluxInt;
        }
    }
    for(int iel = 0; iel<nInterfacesByElement.size(); iel++){
        if(abs(SumFluxByElement[iel])>1.0e-8){
            std::cout << "The sum of the flows on each element must be zero. Element: " << iel << " has a value of " << SumFluxByElement[iel] << std::endl;
            std::cout << "The problematic element in the flux mesh is " << fCellsData.fGeoIndex[iel] << std::endl;
            DebugStop();
        }
    }
    std::cout << "The sum of the flows over the elements is zero. This is correct!" << std::endl;
}
void TPZAlgebraicTransport::PrintFluxes(){
    for (auto interfaID: fInterfaceData) {
        std::cout<<std::endl;
        std::cout<<"Id: "<<interfaID.first<<std::endl;
        std::cout<<std::endl;
        int nFLuxes = interfaID.second.fIntegralFlux.size();
        for (int iflux =0; iflux<nFLuxes; iflux++) {
            std::pair<int, int> lefrig= interfaID.second.fLeftRightVolIndex[iflux];
            std::cout<<interfaID.second.fIntegralFlux[iflux]<<std::endl;
            std::cout<<fCellsData.fEqNumber[lefrig.first]<<" ,"<<fCellsData.fEqNumber[lefrig.second]<<std::endl;
        }
    }
}
void TPZAlgebraicTransport::ZeroFluxes(){
    for( auto interdata: fInterfaceData){
        TInterfaceDataTransport &transport  =interdata.second;
        int neles = transport.fLeftRightVolIndex.size();
        for(int iel = 0; iel<neles; iel++){
            transport.fIntegralFlux[iel] = 0.0;
        }
    }
}
