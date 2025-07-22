//
//  TMRSDarcyFlowWithMem.h
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSDarcyFlowWithMem_h
#define TMRSDarcyFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "TPZMatWithMem.h"

template <class TMEM>
class TMRSDarcyFlowWithMem : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>, TPZMatErrorCombinedSpaces<STATE> >{
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>, TPZMatErrorCombinedSpaces<STATE> >;
    
protected:
    /// Dimension
    int m_dimension;
    
    /// Scale factor for pressure
    STATE m_scale_pressure = 1.0;//1.e-6;

    /// Scale factor for flux variable
    STATE m_scale_flux = 1.0;
    
    /// Directive that stands for the use of four approximations spaces (iterative method)
    bool m_is_four_spaces_Q;
    
    TMRSDataTransfer mSimData;

    bool m_is_axisymmetric;
    
public:
    
    /// Default constructor
    TMRSDarcyFlowWithMem();
    
    /// Constructor based on a material id
    TMRSDarcyFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFlowWithMem(const TMRSDarcyFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFlowWithMem &operator=(const TMRSDarcyFlowWithMem & other);
    
    /// Default destructor
    ~TMRSDarcyFlowWithMem();
    
    /// Set the required data at each integration point
    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;
    
    /**
     * @brief Returns the number of errors to be evaluated
     *
     * Returns the number of errors to be evaluated, that is, the number of error norms associated
     * with the problem.
     */
    int NEvalErrors() const override { return 5; }
    
    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    
    /// Returns the name of the material
    std::string Name() const override {
        return "TMRSDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial()  const override
    {
        return new TMRSDarcyFlowWithMem(*this);
    }
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) const override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) const override;
    
    /// returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var) const override;
    
    /// Set the axisymmetry flag
    void SetAxisymmetry(bool IsAxisymmetric) {m_is_axisymmetric = IsAxisymmetric;}

    /// Returns the axisymmetry flag
    bool IsAxisymmetric() const {return m_is_axisymmetric;}

    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout)  override;
    
   
    
    // Contribute Methods being used
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override ;
    
    void ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) ;
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
};

#endif /* TMRSDarcyFlowWithMem_h */
