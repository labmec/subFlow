//
//  AlgebraicDataTransfer.hpp
//  ALL_BUILD
//
//  Created by Jose on 6/17/20.
//



#ifndef TPZPostProcessResProp_h
#define TPZPostProcessResProp_h

#include <stdio.h>
// #include "pzdiscgal.h"
#include "TPZBndCondT.h"

#include "TPZMatBase.h"
#include "TPZMaterialDataT.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatInterfaceSingleSpace.h"
class TPZPostProcessResProp : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatInterfaceSingleSpace<STATE>> {
    
    using TBase =TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatInterfaceSingleSpace<STATE>>;
private:
    
    /** @brief material dimension */ // TODO:: Candidate for deletion
    int m_dimension;
    
    /** @brief material dimension */
    int m_mat_id;
    
    /** @brief Directive that stands for Mass matrix assembly  */
    bool m_mass_matrix_Q;
    
    /** @brief Regular time step size  */
    REAL m_dt;
    
    /** @brief Porosity  */
    REAL m_phi;
    
    REAL m_fracture_epsilon;
    
public:
    
    /** @brief Default constructor */
    TPZPostProcessResProp();
    
    /** @brief Constructor based on a material id */
    TPZPostProcessResProp(int matid, int dimension);
    
    /** @brief Constructor based on a TPZTracerFlow object */
    TPZPostProcessResProp(const TPZPostProcessResProp &other);
    
    /** @brief Assignment operator */
    TPZPostProcessResProp &operator=(const TPZPostProcessResProp &other);
    
    /** @brief Default destructor */
    ~TPZPostProcessResProp();
    
    /** @brief Set the required data at each integration point */
    virtual void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const override ;
    
    /** @brief Set the required data at each integration point */
    virtual void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data) const override;
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right) ;
    
    /** @brief Returns the name of the material */
    virtual std::string Name() const override{
        return "TPZPostProcessResProp";
    }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const {return m_dimension;}
    
    /** @brief Sets material dimension */
    void SetDimension(int dim) { m_dimension = dim; }
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() const {return 1;} // Deprecated, must to be removed
    
    /** @brief Returns material copied form this object */
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TPZPostProcessResProp(*this);
    }
    
    /** @brief Print out the data associated with the material */
    void Print(std::ostream &out = std::cout) const override;
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
    
    /** @brief Returns the number of variables associated with varindex */
    int NSolutionVariables(int var) const override;
    
    /** @brief Not used contribute methods */
    
    
    // Contribute Methods being used
    
    /** @brief Returns the solution associated with the var index */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    
    /**
     * Unique identifier for serialization purposes
     */
    int ClassId() const;
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context) override;
    
    /** @brief Set directive that stands for Mass matrix assembly  */
    void SetMassMatrixAssembly(bool mass_matrix_Q){
        m_mass_matrix_Q = mass_matrix_Q;
    }
    
    /** @brief Set directive that stands for Mass matrix assembly  */
    bool GetMassMatrixAssembly(){
        return m_mass_matrix_Q;
    }
    
    /** @brief Set regular time step size  */
    void SetTimeStep(REAL dt){
        m_dt = dt;
    }
    
    /** @brief Get regular time step size  */
    REAL GetTimeStep(){
        return m_dt;
    }
    
    /** @brief Set porosity  */
    void SetPorosity(REAL phi){
        m_phi = phi;
    }
    
    /** @brief Get porosity  */
    REAL GetPorosity(){
        return m_phi;
    }
    
    /** @brief Set fracture cross length  */
    void SetFractureCrossLength(REAL fracture_epsilon){
        m_fracture_epsilon = fracture_epsilon;
    }
    
    /** @brief Get fracture cross length  */
    REAL GetFractureCrossLength(){
        return m_fracture_epsilon;
    }
    
    REAL FractureFactor(TPZMaterialData & data);
    
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc)override{
        DebugStop();
        
    }
    virtual void
    SolutionInterface(const TPZMaterialDataT<STATE> &data,
                      const TPZMaterialDataT<STATE> &dataleft,
                      const TPZMaterialDataT<STATE> &dataright,
                      const int var,
                      TPZVec<STATE> &Solout)override{
        DebugStop();
    }
    virtual void
    ContributeInterface(const TPZMaterialDataT<STATE> &data,
                        const TPZMaterialDataT<STATE> &dataleft,
                        const TPZMaterialDataT<STATE> &dataright,
                        REAL weight, TPZFMatrix<STATE> &ek,
                        TPZFMatrix<STATE> &ef)override{
        DebugStop();
    }
    virtual void
    ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                          const TPZMaterialDataT<STATE> &dataleft, REAL weight,
                          TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                          TPZBndCondT<STATE> &bc) override{
        DebugStop();
    }

};

#endif /* defined(TPZTracerFlow) */
