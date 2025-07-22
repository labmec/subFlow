//
//  TPZMixedDarcyFlow.h
//  HDiv
//
//  Created by Omar Durán on 3/29/19.
//

#ifndef TPZMixedDarcyFlow_h
#define TPZMixedDarcyFlow_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#include "TPZMaterialDataT.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
// #include "TPZMatErrorCombinedSpaces.h"
// #include "TPZDarcyFlowInterface.h"

class TPZMixedDarcyFlow : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>> {
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>>;

private:
    
    /** @brief Absolute permeability */
    TPZFNMatrix<9,REAL> m_kappa;
    
    /** @brief Inver of absolute permeability */
    TPZFNMatrix<9,REAL> m_kappa_inv;
    
    /** @brief Dimensional factor */
    REAL m_d;
    
    /** @brief Material dimension */
    int m_dim;
    
public:
    
    /** @brief default constructor  */
    TPZMixedDarcyFlow();
    
    /** @brief default desconstructor  */
    ~TPZMixedDarcyFlow();
    
    /** @brief constructor based on material id and dimension */
    TPZMixedDarcyFlow(int mat_id, int dim);
    
    /** @brief constructor based on object copy */
    TPZMixedDarcyFlow(const TPZMixedDarcyFlow &other);
    
    /** @brief return a material object from a this object */
    TPZMaterial * NewMaterial() const override;
    
    /** @brief assignment operator */
    TPZMixedDarcyFlow &operator=(const TPZMixedDarcyFlow &other);
    
    /** @brief return the euclidean dimension of the weak statement */
    int Dimension() const override;
    
    /** @brief Sets material dimension */
    void SetDimension(int dim) { m_dim = dim; }
    
    /** @brief return the number of state variables associated with each trial function */
    virtual int NStateVariables() const override;
    
    /** @brief print all material information */
    void Print(std::ostream & out) const override;
    
    /** @brief print material name */
    std::string Name() const override;
    
    /** @brief fill requirements for volumetric contribute methods multiphsycis mesh */

    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /** @brief fill requirements for boundary contribute methods multiphsycis mesh */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    
    /** @brief unique class identifier */
    virtual int ClassId() const override;
    
    /** @brief write class in disk */
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief write class from disk */
    void Read(TPZStream &buf, void *context) override;
    

   
    
    /** @brief Volumetric contribute with jacobian matrix */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;
    
    /** @brief Boundary contribute */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc) override;

    
    /** @brief Variable index based on variable naming */
    int VariableIndex(const std::string &name) const override;
    
    /** @brief size of the current variable (1 -> scalar, 3-> vector, 9 ->  Tensor ) */
    int NSolutionVariables(int var) const override;
    
    /** @brief Postprocess required variables multiphysics */

    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;
    /** @brief Set Permeability */
    void SetPermeability(REAL kappa){
        m_kappa.Zero();
        m_kappa_inv.Zero();
        for (int i = 0; i < 3; i++) {
            m_kappa(i,i) = kappa;
            m_kappa_inv(i,i) = 1.0/kappa;
        }
    }
    
    /** @brief Set Permeability */
    void SetPermeability(TPZFNMatrix<9,REAL> & kappa){
        
        if (kappa.Rows() != 3 && kappa.Cols() != 3) {
            DebugStop();
        }
        
        m_kappa = kappa;
        m_kappa.Inverse(m_kappa_inv, ELU);
    }
    /** @brief Gravity field */
    std::vector<REAL> m_gravity;
    /** @brief Sets Gravity field */
    void SetGravity(std::vector<REAL> & gravity){
        m_gravity = gravity;
    }
    
    /** @brief Set fracture dim */
    void SetDimensionalFactor(REAL d){
        m_d = d;
    }
    
};

#endif /* TPZMixedDarcyFlow_h */
