
#include "TPZPostProcessResProp.h"


TPZPostProcessResProp::TPZPostProcessResProp(){

}

/** @brief Constructor based on a material id */
TPZPostProcessResProp::TPZPostProcessResProp(int matid, int dimension) : TBase(matid) {

    m_mat_id = matid;
    m_dimension = dimension;
    m_mass_matrix_Q = false;
    m_dt = 0.1;
    m_phi = 1.0;
    m_fracture_epsilon = 0.0;
    
}

/** @brief Constructor based on a TRMMultiphase object */
TPZPostProcessResProp::TPZPostProcessResProp(const TPZPostProcessResProp &other) : TBase(other) {
    m_mat_id = other.m_mat_id;
    m_dimension = other.m_dimension;
    m_mass_matrix_Q = other.m_mass_matrix_Q;
    m_dt = other.m_dt;
    m_phi = other.m_phi;
    m_fracture_epsilon = other.m_fracture_epsilon;
}

TPZPostProcessResProp & TPZPostProcessResProp::operator=(const TPZPostProcessResProp &other){
    if (this != & other) // prevent self-assignment
    {
        TPZMaterial::operator=(other);
        m_mat_id = other.m_mat_id;
        m_dimension = other.m_dimension;
        m_mass_matrix_Q = other.m_mass_matrix_Q;
        m_dt = other.m_dt;
        m_phi = other.m_phi;
    }
    return *this;
}

/** @brief Default destructor */
TPZPostProcessResProp::~TPZPostProcessResProp(){

}

/** @brief Set the required data at each integration point */
void TPZPostProcessResProp::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

/** @brief Set the required data at each integration point */
void TPZPostProcessResProp::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec)const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TPZPostProcessResProp::FillDataRequirementsInterface(TPZMaterialData &data)const{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
//    if(fLinearContext == false){
//        data.fNeedsNeighborSol = true;
//    }
}

void TPZPostProcessResProp::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    int nref_left = datavec_left.size();
    for(int iref = 0; iref<nref_left; iref++){
        datavec_left[iref].SetAllRequirements(false);
        datavec_left[iref].fNeedsSol = true;
        datavec_left[iref].fNeedsNormal = true;
    }
    int nref_right = datavec_right.size();
    for(int iref = 0; iref<nref_right; iref++){
        datavec_right[iref].SetAllRequirements(false);
        datavec_right[iref].fNeedsSol = true;
    }
    
}


/** @brief Print out the data associated with the material */
void TPZPostProcessResProp::Print(std::ostream &out) const{
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

/** @brief Returns the variable index associated with the name */
int TPZPostProcessResProp::VariableIndex(const std::string &name)const{

    if (!strcmp("Porosity", name.c_str())) return 0;
    if (!strcmp("Permeability_x", name.c_str())) return 1;
    if (!strcmp("Permeability_y", name.c_str())) return 2;
    if (!strcmp("Permeability_z", name.c_str())) return 3;

    return TPZMaterial::VariableIndex(name);
}

/** @brief Returns the number of variables associated with varindex */
int TPZPostProcessResProp::NSolutionVariables(int var)const{
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Scalar
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
 

    }
    return TPZMaterial::NSolutionVariables(var);
}



// Contribute Methods being used
/** @brief Returns the solution associated with the var index */
void TPZPostProcessResProp::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout){

    int phiBlock    = 0;
    int KxBlock    = 1;
    int KyBlock    = 2;
    int KzBlock    = 3;
   
    
    REAL phi = datavec[phiBlock].sol[0][0];
    REAL Kx = datavec[KxBlock].sol[0][0];
    REAL Ky = datavec[KyBlock].sol[0][0];
    REAL Kz = datavec[KzBlock].sol[0][0];
    Solout.Resize(this->NSolutionVariables(var));

    switch(var) {
        case 0:
        {
            Solout[0] = phi;
        }
            break;
        case 1:
        {
            Solout[0] = Kx;
        }
            break;
        case 2:
        {
            Solout[0] = Ky;
        }
            break;
        case 3:
        {
            Solout[0] = Kz;
        }
            break;
    }


}

REAL TPZPostProcessResProp::FractureFactor(TPZMaterialData & data){
    return 0.;
}

void TPZPostProcessResProp::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    return;
}

void TPZPostProcessResProp::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->Contribute(datavec, weight, ek_fake, ef);
    
}


/**
 * Unique identifier for serialization purposes
 */
int TPZPostProcessResProp::ClassId() const{
    DebugStop();
    return  -1942;
}

/**
 * Save the element data to a stream
 */
void TPZPostProcessResProp::Write(TPZStream &buf, int withclassid){
    DebugStop();
}

/**
 * Read the element data from a stream
 */
void TPZPostProcessResProp::Read(TPZStream &buf, void *context){
    DebugStop();
}

