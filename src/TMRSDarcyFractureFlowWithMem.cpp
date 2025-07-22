//
//  TMRSDarcyFlowWithMem.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.

#include "TMRSMemory.h"

#include "TMRSDarcyFractureFlowWithMem.h"

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::TMRSDarcyFractureFlowWithMem() : TBase() {
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::TMRSDarcyFractureFlowWithMem(int mat_id, int dimension) : TBase(mat_id, dimension) {
    
    
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::TMRSDarcyFractureFlowWithMem(const TMRSDarcyFractureFlowWithMem & other) : TMRSDarcyFlowWithMem<TMEM>(other){
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM> & TMRSDarcyFractureFlowWithMem<TMEM>::operator=(const TMRSDarcyFractureFlowWithMem & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    return *this;
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::~TMRSDarcyFractureFlowWithMem(){
    
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fDeformedDirections = true;
    }
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::SetDataTransfer(TMRSDataTransfer & SimData){
    this->mSimData = SimData;
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Print(std::ostream &out) const{
    TPZMaterial::Print(out);
}

template <class TMEM>
int TMRSDarcyFractureFlowWithMem<TMEM>::VariableIndex(const std::string &name) const{
    return TMRSDarcyFlowWithMem<TMEM>::VariableIndex(name);
}

template <class TMEM>
int TMRSDarcyFractureFlowWithMem<TMEM>::NSolutionVariables(int var) const{
    return TMRSDarcyFlowWithMem<TMEM>::NSolutionVariables(var);
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    if(datavec[0].divsol.size() > 0) {
#ifdef PZDEBUG
        if(datavec[0].divsol[0].size() != 2) DebugStop();
#endif
        // overwrite the total divergence for the inplane divergence
        datavec[0].divsol[0][0] = datavec[0].divsol[0][1];
    }
    TMRSDarcyFlowWithMem<TMEM>::Solution(datavec,var,Solout);
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int qb = 0;
    int pb = 1;
    
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].fDeformedDirections; // fDeformedDirections stores the shape function for flux
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    
    TPZFNMatrix<40, REAL> div_phi = datavec[qb].divphi;
    REAL div_q = datavec[qb].divsol[0][0];
    
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    int nvecs       = phi_qs.Cols();
    
    // first index in fVecShapeIndex corresponding to the first transverse flux
    int first_transverse_q = 0;
    // first index in fVecShapeIndex corresponding to the second transverse flux
    int second_transverse_q = 0;
    int nconnects = datavec[qb].fHDiv.fNumConnectShape.size();    
    second_transverse_q = nvecs-datavec[qb].fHDiv.fNumConnectShape[nconnects-1];
    first_transverse_q = second_transverse_q-datavec[qb].fHDiv.fNumConnectShape[nconnects-2];

    if(first_transverse_q == 0 || second_transverse_q == 0 || first_transverse_q == second_transverse_q)
    {
        DebugStop();
    }
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
  
    // Get the data at integrations points
    long gp_index = datavec[qb].intGlobPtIndex;
    if(gp_index < 0) DebugStop();
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    TPZFNMatrix<3,STATE> phi_q_i(3,1,0.0), kappa_inv_phi_q_j(3,1,0.0), kappa_inv_q(3,1,0.0);
    
    REAL ad = memory.m_ad;
    for (int i = 0; i < 3; i++) {
        if(memory.m_kappa_inv(i,i) != memory.m_kappa_inv((i+1)%3,(i+1)%3)) DebugStop();
        for (int j = 0; j < 3; j++) {
            if(i!=j && memory.m_kappa_inv(i,j) != 0.) DebugStop();
            kappa_inv_q(i,0) += memory.m_kappa_inv(i,j)*q[j];
        }
    }
    
    for (int iq = 0; iq < first_transverse_q; iq++)
    {

        STATE kappa_inv_q_dot_phi_q_i = 0.0;
        for (int i = 0; i < 3; i++) {
            kappa_inv_q_dot_phi_q_i += kappa_inv_q(i,0)*phi_qs(i,iq);
        }
       // @TODO: Verify fracture ad here!!!
//        ef(iq + first_q) += weight*(-kappa_inv_q_dot_phi_q_i /ad + p * div_phi(iq,0));  // terms j and k in docs/Formulation.lyx
        ef(iq + first_q) += weight*(-kappa_inv_q_dot_phi_q_i  + p * div_phi(iq,0));  // terms j and k in docs/Formulation.lyx
        for (int jq = 0; jq < first_transverse_q; jq++) {

                kappa_inv_phi_q_j.Zero();
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        kappa_inv_phi_q_j(i,0) += (memory.m_kappa_inv(i,j) / ad) * phi_qs(j,jq);
                    }
                }
                
                STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
                for (int j = 0; j < 3; j++) {
                    kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_qs(j,iq);
                }
                
                ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i; // term s in docs/Formulation.lyx
            }
            
            for (int jp = 0; jp < nphi_p; jp++) {
                ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term t in docs/Formulation.lyx
                ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term x in docs/Formulation.lyx
            }
    }
    for (int ip = 0; ip < nphi_p; ip++) {
//        ef(ip + first_p) = weight * phi_ps(ip,0) * div_q; // term p in docs/Formulation.lyx
    }
    
//    
    // compute the contribution of the hdivbound
    for (int iq = first_transverse_q; iq < second_transverse_q; iq++)
    {
        //datavec[qb].Print(std::cout);
        STATE kappa_inv_q_dot_phi_q_i = 0.0;

        for (int i = 0; i < 3; i++) {
            kappa_inv_q(i,0) += memory.m_kappa_inv(i,i)*q[i];
        }
        
        for (int i = 0; i < 3; i++) {
            kappa_inv_q_dot_phi_q_i += kappa_inv_q(i,0)*phi_qs(i,iq);
        }

        ef(iq + first_q) += weight*(-0.5*ad*kappa_inv_q_dot_phi_q_i + p * div_phi(iq,0)); // terms c2 and c1 in docs/Formulation.lyx
//        ef(iq + first_q) += weight*(-kappa_inv_q_dot_phi_q_i + p * div_phi(iq,0)); // terms c2 and c1 in docs/Formulation.lyx

        for (int jq = first_transverse_q; jq < second_transverse_q; jq++) {

            kappa_inv_phi_q_j.Zero();

            for (int j = 0; j < 3; j++) {
                REAL KappaInvVal = memory.m_kappa_inv(j,j);
                kappa_inv_phi_q_j(j,0) = 0.5 * ad * KappaInvVal * phi_qs(j,jq);
            }

            // kappa_inv_phi_q_j_dot_phi_q_i is the dot product of both vectors
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_qs(j,iq);
            }
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i; // term h2 in docs/Formulation.lyx
        }
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term y in docs/Formulation.lyx (they really contribute as a divergent!)
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term h1 in docs/Formulation.lyx
        }
    }
//    // compute the contribution of the hdivbound
    for (int iq = second_transverse_q; iq < nphi_q; iq++) {

        STATE kappa_inv_q_dot_phi_q_i = 0.0;

        for (int i = 0; i < 3; i++) {
            kappa_inv_q(i,0) += memory.m_kappa_inv(i,i)*q[i];
        }
        
        for (int i = 0; i < 3; i++) {
            kappa_inv_q_dot_phi_q_i += kappa_inv_q(i,0)*phi_qs(i,iq);
        }
        
        ef(iq + first_q) += weight *((-0.5*ad)*kappa_inv_q_dot_phi_q_i + p * div_phi(iq,0)); // part of terms c2 and c1 in docs/Formulation.lyx
//        ef(iq + first_q) += weight *(-kappa_inv_q_dot_phi_q_i + p * div_phi(iq,0)); // part of terms c2 and c1 in docs/Formulation.lyx

        for (int jq = second_transverse_q; jq < nphi_q; jq++)
        {
            kappa_inv_phi_q_j.Zero();

            for (int j = 0; j < 3; j++) {
                REAL KappaInvVal = memory.m_kappa_inv(j,j);
                kappa_inv_phi_q_j(j,0) += 0.5 * ad * KappaInvVal * phi_qs(j,jq);
            }
            // kappa_inv_phi_q_j_dot_phi_q_i is the dot product of both vectors
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_qs(j,iq);
            }
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i; // term h2 in docs/Formulation.lyx
        }
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term y in docs/Formulation.lyx
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term h1 in docs/Formulation.lyx
        }
    }
    
    if(this->mSimData.mTNumerics.m_four_approx_spaces_Q){
        ContributeFourSpaces(datavec,weight,ek,ef);
    }
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    int qb = 0;
    int pb = 1;
    int numactive = 0;
    for(auto &it : datavec) if(it.fActiveApproxSpace) numactive++;
    if(numactive%2 != 0) DebugStop();
    int numaverage = (numactive-2)/2;
    TPZFMatrix<REAL> &phi_ps = datavec[pb].phi;
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    if(nphi_q+nphi_p+numaverage*2 != ek.Rows())
    {
        DebugStop();
    }
    for(int iavg = 0; iavg<numaverage; iavg++)
    {
        int g_avgb = 2+2*iavg;
        int p_avgb = 3+2*iavg;
    
    
        int nphi_gb = datavec[g_avgb].phi.Rows();
        int nphi_pb = datavec[p_avgb].phi.Rows();
    
        if(nphi_gb != 1 || nphi_pb != 1) DebugStop();
    
    
        STATE p     = datavec[pb].sol[0][0];
        STATE g_avg = datavec[g_avgb].sol[0][0];
        STATE p_avg = datavec[p_avgb].sol[0][0];
    
        for(int ip=0; ip<nphi_p; ip++)
        {
            ef(nphi_q+ip,0) += weight * g_avg * phi_ps(ip,0);
            ek(nphi_q+ip,nphi_q+nphi_p+2*iavg) += weight * phi_ps(ip,0);
            ek(nphi_q+nphi_p+2*iavg,nphi_q+ip) += weight * phi_ps(ip,0);
        }
        
        ef(nphi_q+nphi_p+1+2*iavg,0) += -weight * g_avg;
        ek(nphi_q+nphi_p+1+2*iavg,nphi_q+nphi_p+2*iavg) += -weight;
        
        ef(nphi_q+nphi_p+2*iavg,0) += weight * (p - p_avg);
        ek(nphi_q+nphi_p+2*iavg,nphi_q+nphi_p+1+2*iavg) += -weight;
    }

    
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){

    this->Contribute(datavec, weight, ef);
    
    if(TMRSDarcyFractureFlowWithMem<TMEM>::fUpdateMem){
        int qb = 0;
        int pb = 1;
        long gp_index = datavec[pb].intGlobPtIndex;
        TMEM & memory = this->GetMemory().get()->operator[](gp_index);
        TPZVec<REAL> q_n = {datavec[qb].sol[0][0]};
        memory.m_flux = q_n;
        
        REAL p_n = datavec[pb].sol[0][0];
        memory.m_p = p_n;
    }
    
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    
    REAL gBigNumber = 1.0e12; //TPZMaterial::gBigNumber;
    
    int qb = 0;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()[0];
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += weight * (-1.0)*p_D * phi_qs(iq,0); // term o in docs/Formulation.lyx
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                REAL qn_N = bc_data[0];
                REAL qn = 0.0;
                qn = q[0];
                
                ef(iq + first_q) += weight * gBigNumber * (-1.0)*(qn - qn_N) * phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    ek(iq + first_q,jq + first_q) += weight * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
                }
                
            }
            
        }
            break;
            
        case 2 :    // Mixed BC
        {
            const REAL v1 = bc.Val1()(0,0);
            for (int iq = 0; iq < nphi_q; iq++){
                const REAL qn_N = bc_data[0];
                REAL qn = 0.0;
                qn = q[0];
                ef(iq + first_q) += weight *(-1.0)* (qn - qn_N) * phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++){
                    ek(iq + first_q,jq + first_q) += weight/v1 * phi_qs(jq,0) * phi_qs(iq,0);
                }
            }
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

template class TMRSDarcyFractureFlowWithMem<TMRSMemory>;

























