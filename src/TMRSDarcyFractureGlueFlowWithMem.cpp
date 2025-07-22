//
//  TMRSDarcyFractureGlueFlowWithMem.cpp
//
//  Created by José Villegas on 30/05/22.

#include "TMRSMemory.h"

#include "TMRSDarcyFractureGlueFlowWithMem.h"


TMRSDarcyFractureGlueFlowWithMem::TMRSDarcyFractureGlueFlowWithMem() : TBase() {
}


TMRSDarcyFractureGlueFlowWithMem::TMRSDarcyFractureGlueFlowWithMem(int mat_id, REAL perm) : TBase(mat_id), m_permeability(perm) {

}


TMRSDarcyFractureGlueFlowWithMem::~TMRSDarcyFractureGlueFlowWithMem(){
    
}

void TMRSDarcyFractureGlueFlowWithMem::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        //datavec[idata].fNormalVec = true;
    }
    datavec[0].fNeedsSol = true;
}

void TMRSDarcyFractureGlueFlowWithMem::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
    }
}

void TMRSDarcyFractureGlueFlowWithMem::SetDataTransfer(TMRSDataTransfer & SimData){
    //this->mSimData = SimData;
}

void TMRSDarcyFractureGlueFlowWithMem::Print(std::ostream &out) const{
    TPZMaterial::Print(out);
}

int TMRSDarcyFractureGlueFlowWithMem::VariableIndex(const std::string &name) const{
    if(!strcmp("Flux",name.c_str()))            return  1;
    return TPZMaterial::VariableIndex(name);
}

int TMRSDarcyFractureGlueFlowWithMem::NSolutionVariables(int var) const{
    if(var == 1) return 1;
    return TBase::NSolutionVariables(var);
}

void TMRSDarcyFractureGlueFlowWithMem::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    DebugStop();
}

void TMRSDarcyFractureGlueFlowWithMem::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int qb = 0;
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    // Get the data at integrations points
    long gp_index = datavec[qb].intGlobPtIndex;
    if(gp_index < 0) DebugStop();
    auto & memory = this->GetMemory().get()->operator[](gp_index);
    REAL kappaNormal = m_permeability;
    // Caso normal kappaNormal = kappaNormal
    // Caso comunicante kappaNormal = 1.e12*kappaNormal
    // Caso nao comunicante kappaNormal = 1.e-8*kappaNormal
    REAL dist = memory.m_dist;
    memory.m_flux = q[0];

    REAL fact  = weight*dist/kappaNormal;
    auto &phi = datavec[0].phi;
    int nphi = datavec[qb].phi.Rows();
    for (int i=0; i<nphi; i++) {
        ef(i,0) += -1.0*phi(i)*fact*q[0];
        for (int j=0; j<nphi; j++) {
            ek(i,j) += fact*phi(i)*phi(j);
        }
    }
}


void TMRSDarcyFractureGlueFlowWithMem::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    int qb = 0;
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    // Get the data at integrations points
    long gp_index = datavec[qb].intGlobPtIndex;
    if(gp_index < 0) DebugStop();
    auto & memory = this->GetMemory().get()->operator[](gp_index);
    REAL kappaNormal = m_permeability;
    REAL dist = memory.m_dist;
    memory.m_flux = q[0];
    REAL fact  = weight*dist/kappaNormal;
    auto &phi = datavec[0].phi;
    int nphi = datavec[qb].phi.Rows();
    for (int i=0; i<nphi; i++) {
        ef(i,0) += -1.0*phi(i)*fact*q[0];
    }
}

void TMRSDarcyFractureGlueFlowWithMem::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    DebugStop();
}

void TMRSDarcyFractureGlueFlowWithMem::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    DebugStop();
}
