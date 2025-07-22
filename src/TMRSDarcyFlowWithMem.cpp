//
//  TMRSDarcyFlowWithMem.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.

#include "TMRSDarcyFlowWithMem.h"
#include "TRMSpatialPropertiesMap.h"

#include "TMRSMemory.h"
#include "pzvec_extras.h"

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::TMRSDarcyFlowWithMem() : TBase(), mSimData() {
    m_dimension = 0;
    m_is_four_spaces_Q = false;
    m_scale_pressure = 1;
    m_scale_flux = 1;
    m_is_axisymmetric = false;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::TMRSDarcyFlowWithMem(int mat_id, int dimension) : TBase(mat_id), mSimData(){
//    DebugStop();
    m_dimension = dimension;
    m_is_four_spaces_Q = false;
    m_scale_pressure = 1;
    m_scale_flux = 1;
    m_is_axisymmetric = false;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::TMRSDarcyFlowWithMem(const TMRSDarcyFlowWithMem & other) : TBase(other){
    m_dimension         = other.m_dimension;
    m_scale_pressure    = other.m_scale_pressure;
    m_scale_flux        = other.m_scale_flux;
    m_is_four_spaces_Q  = other.m_is_four_spaces_Q;
    m_is_axisymmetric   = other.m_is_axisymmetric;
    mSimData  = other.mSimData;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM> & TMRSDarcyFlowWithMem<TMEM>::operator=(const TMRSDarcyFlowWithMem & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    m_dimension         = other.m_dimension;
    m_scale_pressure    = other.m_scale_pressure;
    m_scale_flux        = other.m_scale_flux;
    m_is_four_spaces_Q  = other.m_is_four_spaces_Q;
    m_is_axisymmetric   = other.m_is_axisymmetric;
    mSimData  = other.mSimData;
    return *this;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::~TMRSDarcyFlowWithMem(){
    
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec)const {
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
//        datavec[idata].fDeformedDirections = true;
    }
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::SetDataTransfer(TMRSDataTransfer & SimData){
    mSimData = SimData;
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Print(std::ostream &out) const{
//    TPZMatWithMem<TMEM>::Print(out);
    out << m_dimension << std::endl;
    out << m_scale_pressure << std::endl;
    out << m_scale_flux << std::endl;
    out << m_is_four_spaces_Q << std::endl;
    out << m_is_axisymmetric << std::endl;
}

template <class TMEM>
int TMRSDarcyFlowWithMem<TMEM>::VariableIndex(const std::string &name) const {
    if(!strcmp("Flux",name.c_str()))            return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    if(!strcmp("div_q",name.c_str()))           return  3;
    if(!strcmp("kappa",name.c_str()))           return  4;
    if(!strcmp("g_average",name.c_str()))       return  5;
    if(!strcmp("p_average",name.c_str()))       return  6;
    if(!strcmp("GradPressure", name.c_str()))   return 7;
    if(!strcmp("BCNormalFlux", name.c_str()))   return 8;
    if(!strcmp("FluxNorm", name.c_str()))   return 9;
    if(!strcmp("div_positive", name.c_str()))   return 10;
    if(!strcmp("div_negative", name.c_str()))   return 11;
    return TPZMaterial::VariableIndex(name);
}

template <class TMEM>
int TMRSDarcyFlowWithMem<TMEM>::NSolutionVariables(int var) const{
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    if(var == 7) return 3;
    if(var == 8) return 1;
    if(var == 9) return 1;
    if(var == 10) return 1;
    if(var == 11) return 1;
    return TPZMaterial::NSolutionVariables(var);
}

template <class TMEM>
int TMRSDarcyFlowWithMem<TMEM>::ClassId() const {
    return Hash("TMRSDarcyFlowWithMem") ^ TBase::ClassId() << 1;
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int qb = 0;
    int pb = 1;
    Solout.Resize( this->NSolutionVariables(var));
    TPZManVector<STATE,3> p, q;
    
    q = datavec[qb].sol[0];
    p = datavec[pb].sol[0];

    // Assuming radius is aligned with the x axis
    REAL axi_weight = 1.0;
    if (m_is_axisymmetric)
    {
        REAL r = datavec[qb].x[0];
        axi_weight = 1.0 / (2.0 * M_PI * r);
        q[0] *= axi_weight;
        q[1] *= axi_weight;
    }

    if(var == 8){ // normal flux. Only works for bc materials
        Solout[0] = q[0];
        return;
    }
    
    if(var == 1){
        for (int i=0; i < 3; i++)
        {
            Solout[i] = q[i];
        }
        return;
    }
    
    // fluxnorm
    if(var == 9) {
        Solout[0] = 0.;
        for (int i=0; i < 3; i++)
        {
            Solout[0] += q[i]*q[i];
        }
        Solout[0] = sqrt(Solout[0]);
        return;

    }
    
    if(var == 2){
        Solout[0] = p[0];
        return;
    }
    
    if(var == 3){
        REAL div_q = datavec[qb].divsol[0][0];
        Solout[0] = div_q;
        return;
    }
    // div positive
    if(var == 10) {
        REAL div_q = datavec[qb].divsol[0][0];
        Solout[0] = div_q > 0.? div_q: 0.;
        return;
    }
    // div negative
    if(var == 11) {
        REAL div_q = datavec[qb].divsol[0][0];
        Solout[0] = div_q < 0.? div_q: 0.;
        return;
    }

    if(var == 4){
     
            TPZManVector<double, 3> point;
        
            point = datavec[qb].XCenter;
        
            int val = rand() % 100;
            
            REAL kappa = 1.0;
            if (val<75) {
                kappa =  10000;
                
            }
            else{
                kappa =  10;
            }
        
            Solout[0] = kappa;
        
        return;
    }
    
    if (mSimData.mTNumerics.m_four_approx_spaces_Q) {
    
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
    }
    
    if (var == 9) {

        if(datavec[1].fShapeType == TPZMaterialData::EEmpty) return;
        TPZFNMatrix<9, REAL> dsoldx(3, 1.,0.);
        TPZFNMatrix<9, REAL> dsoldaxes(m_dimension, 1,0.);

        dsoldaxes = datavec[1].dsol[0];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        for (int i = 0; i < m_dimension; i++) {
            Solout[i] = dsoldx(i, 0);
        }

        return;
    }
    
//    DebugStop();
    
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // Note: All the matrix and residue terms can be compared with the equations in docs/Formulation.lyx.
    //       These terms are accompanied by letter a,b,c,... which related the equatiosn with the c++ code.
    // Note2: The saturation dependent term are implemente in TPZFastCondensedCompEl
    
    // Flux is 0 and pressure is 1
    int qb = 0;
    int pb = 1;
    
    // Getting the shape functions for flux and pressure
//    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<110,REAL> phi_qs       = datavec[qb].fDeformedDirections; // fDeformedDirections stores the shape function for flux
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    // Getting divergent of phi and of flux (in case it is non linear)
    TPZFNMatrix<40, REAL> div_phi = datavec[qb].divphi;
    REAL div_q = datavec[qb].divsol[0][0];
    
    // Number of shape function for flux and pressure
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    
    // Flux and pressure solution (for non linear analyzes)
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];

    //Assuming radius is aligned with the x axis
    REAL axi_weight = 1.0;
    if (m_is_axisymmetric)
    {
        REAL r = datavec[qb].x[0];
        axi_weight = 1.0 / (2.0*M_PI*r);
    }
    
    // Get the memory at the integrations point
    long gp_index = datavec[qb].intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    // Helper variables
    TPZFNMatrix<3,STATE> kappa_inv_phi_q_j(3,1,0.0), kappa_inv_q(3,1,0.0);
            
    // kappa_inv_q = K^-1 * q
    // Used for non linear problems
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            kappa_inv_q(i,0) += memory.m_kappa_inv(i,j)*q[j];
        }
    }
    
    std::vector<REAL>& gravity = mSimData.mTNumerics.m_gravity;
    for (int iq = 0; iq < nphi_q; iq++) {
        
        STATE kappa_inv_q_dot_phi_q_i = 0.0;
        STATE g_dot_phi_q_i = 0.0;
        for (int i = 0; i < 3; i++) {
            kappa_inv_q_dot_phi_q_i += kappa_inv_q(i,0)*phi_qs(i,iq)*axi_weight;
            g_dot_phi_q_i += gravity[i]*phi_qs(i,iq);
        }
        
        ef(iq + first_q) += weight * ( - kappa_inv_q_dot_phi_q_i + p * div_phi(iq,0)); // terms a and b in docs/Formulation.lyx
        ef(iq + first_q) += weight * (  g_dot_phi_q_i );

        for (int jq = 0; jq < nphi_q; jq++) {
            kappa_inv_phi_q_j.Zero();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    kappa_inv_phi_q_j(i,0) += memory.m_kappa_inv(i,j) * phi_qs(j,jq);
                }
            }
            
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_qs(j,iq) * axi_weight;
            }
            
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;  // term f in docs/Formulation.lyx
        }
        
        for (int jp = 0; jp < nphi_p; jp++) {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0); // term g in docs/Formulation.lyx
        }
        
    }
    
    for (int ip = 0; ip < nphi_p; ip++) {
        
        ef(ip + first_p) += weight * (div_q) * phi_ps(ip,0); // term e in docs/Formulation.lyx
        
        for (int jq = 0; jq < nphi_q; jq++) {
            ek(ip + first_p, jq + first_q) += -1.0 * weight * div_phi(jq,0) * phi_ps(ip,0); // term i in docs/Formulation.lyx
        }
        
    }
    if(mSimData.mTNumerics.m_four_approx_spaces_Q && datavec.size() < 4) DebugStop();
    
    if(mSimData.mTNumerics.m_four_approx_spaces_Q){
        ContributeFourSpaces(datavec,weight,ek,ef);
    }
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
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
void TMRSDarcyFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
//    this->Contribute(datavec, weight, ekfake, ef);
    // circular call?
    this->Contribute(datavec, weight, ef);
    
    if(TPZMatWithMem<TMEM>::fUpdateMem){
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
void TMRSDarcyFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    

    REAL gBigNumber = 1.0e12; //TPZMaterial::gBigNumber;

    int qb = 0, pb = 1;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    int bctype = bc.Type();
    if (phi_qs.Rows() == 0 && phi_ps.Rows() != 0) {
        bctype = 10;
    }
    
    int nphi_q       = phi_qs.Rows();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];

    //Assuming radius is aligned with the x axis
    REAL axi_weight = 1.0;
    if (m_is_axisymmetric)
    {
        REAL r = datavec[qb].x[0];
        axi_weight = 1.0 / (2.0*M_PI*r);
        q[0] *= axi_weight;
    }
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()[0];
    const REAL v1 = bc.Val1()(0,0);
    
    if (bc.HasForcingFunctionBC()) {
        TPZFNMatrix<1, STATE> mat_val(1, 1);
        if (bctype == 10) {
            bc.ForcingFunctionBC()(datavec[1].x, bc_data, mat_val);
        }
        else{
            bc.ForcingFunctionBC()(datavec[0].x, bc_data, mat_val);
        }
        
    }
    
    switch (bctype) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += weight *(-1.0)* p_D * phi_qs(iq,0); // term d in docs/Formulation.lyx
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

                ef(iq + first_q) += weight * gBigNumber *(-1.0)* (qn - qn_N) * phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    ek(iq + first_q,jq + first_q) += weight *gBigNumber * phi_qs(jq,0) * phi_qs(iq,0) * axi_weight;
                }
                
            }
            
        }
            break;
        case 2 :    // Mixed BC
        {
            for (int iq = 0; iq < nphi_q; iq++){
                const REAL qn_N = bc_data[0];
                REAL qn = 0.0;
                qn = q[0];
                ef(iq + first_q) += weight * (qn - qn_N) *(-1.0)* phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++){
                    ek(iq + first_q,jq + first_q) += weight/v1 * phi_qs(jq,0) * phi_qs(iq,0) * axi_weight;
                }
            }
        }
            break;
        case 10:
        {
            for (int ip = 0; ip < nphi_p; ip++)
            {
              STATE p_D = bc_data[0];

              ef(ip + first_p) += - weight * gBigNumber * p_D *(-1.0)* phi_ps(ip,0);
              for (int jp = 0; jp < nphi_p; jp++)
              {
                  ek(ip + first_p,jp + first_p) += weight * gBigNumber * phi_ps(jp,0) * phi_ps(ip,0);
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

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {

    /**
     * datavec[0]= Flux
     * datavec[1]= Pressure
     *
     * Errors:
     * [0] L2 for pressure
     * [1] L2 for flux
     * [2] L2 for div(flux)
     * [3] Grad pressure (Semi H1)
     * [4] Hdiv norm
    **/
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1,0);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    auto dsol = data[1].dsol;

    TPZVec<STATE> divsigma(1);

    TPZVec<STATE> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    if(data[1].sol[0].size())
        pressurefem[0] = data[1].sol[0][0];

    long gp_index = data[0].intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    const STATE perm = memory.m_kappa(0,0); // assuming isotropic!
    const STATE inv_perm = 1. / perm;

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexactneg(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexactneg[i] = -perm * gradpressure[i];
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < 3; i++) {
        L2flux += (fluxfem[i] + fluxexactneg[i]) * inv_perm * (fluxfem[i] + fluxexactneg[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
}

template class TMRSDarcyFlowWithMem<TMRSMemory>;
