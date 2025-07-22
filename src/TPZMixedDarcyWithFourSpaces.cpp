//
//  TPZMixedDarcyWithFourSpaces.cpp
//  reservoirlib
//
//  Created by Omar Durán on 7/12/19.
//

#include "TPZMixedDarcyWithFourSpaces.h"


TPZMixedDarcyWithFourSpaces::TPZMixedDarcyWithFourSpaces() : TPZMixedDarcyFlow() {
    
}

TPZMixedDarcyWithFourSpaces::~TPZMixedDarcyWithFourSpaces(){
    
}

TPZMixedDarcyWithFourSpaces::TPZMixedDarcyWithFourSpaces(int mat_id, int dim) : TPZMixedDarcyFlow(mat_id,dim){
    
}

TPZMixedDarcyWithFourSpaces::TPZMixedDarcyWithFourSpaces(const TPZMixedDarcyWithFourSpaces &other) : TPZMixedDarcyFlow(other){
    
}

void TPZMixedDarcyWithFourSpaces::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    TPZMixedDarcyFlow::Contribute(datavec, weight, ek, ef);
    
    int qb = 0;
    int pb = 1;
    int g_avgb = 2;
    int p_avgb = 3;
    
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    
    int nphi_gb = datavec[g_avgb].phi.Rows();
    int nphi_ub = datavec[p_avgb].phi.Rows();
    if(nphi_q+nphi_p+nphi_gb+nphi_ub != ek.Rows())
    {
        DebugStop();
    }
    TPZManVector<STATE> force(1,0.);
    if(fForcingFunction) {
        fForcingFunction(datavec[qb].x,force);
    }
    
    for(int ip=0; ip<nphi_p; ip++)
    {
        ef(nphi_q+ip,0) += weight * force[0]*phi_ps(ip,0);
        ek(nphi_q+ip,nphi_q+nphi_p) += phi_ps(ip,0)*weight;
        ek(nphi_q+nphi_p,nphi_q+ip) += phi_ps(ip,0)*weight;
    }
    ek(nphi_q+nphi_p+1,nphi_q+nphi_p) += -weight;
    ek(nphi_q+nphi_p,nphi_q+nphi_p+1) += -weight;
    
}

void TPZMixedDarcyWithFourSpaces::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ekfake, ef);
}


int TPZMixedDarcyWithFourSpaces::VariableIndex(const std::string &name) const{
    if(!strcmp("g_average",name.c_str()))        return  5;
    if(!strcmp("p_average",name.c_str()))        return  6;
    if(!strcmp("kxx",name.c_str()))        return  100;
    if(!strcmp("kyy",name.c_str()))        return  101;
    if(!strcmp("kzz",name.c_str()))        return  102;
    if(!strcmp("lambda",name.c_str()))     return  103;
    return TPZMixedDarcyFlow::VariableIndex(name);
}

int TPZMixedDarcyWithFourSpaces::NSolutionVariables(int var)const {
    if(var == 5) return 1;
    if(var == 6) return 1;
    if(var == 7) return 100;
    if(var == 8) return 101;
    if(var == 9) return 102;
    if(var == 10) return 103;
    return TPZMixedDarcyFlow::NSolutionVariables(var);
}

void TPZMixedDarcyWithFourSpaces::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout){

    if(var >= 100 && var <= 103) DebugStop();
    int g_avgb = 2;
    int p_avgb = 3;
    
    if(var ==5)
    {
        Solout[0] = datavec[g_avgb].sol[0][0];
        return;
    }
    if(var ==6)
    {
        Solout[0] = datavec[p_avgb].sol[0][0];
        return;
    }
    
    TPZMixedDarcyFlow::Solution(datavec, var, Solout);
}
