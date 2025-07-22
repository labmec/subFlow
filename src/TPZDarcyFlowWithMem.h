//
//  TPZDarcyFlowWithMem.h
//
//  Created by José Villegas on 07/07/20.
//

#ifndef TPZDarcyFlowWithMem_h
#define TPZDarcyFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"
#include "TPZAlgebraicTransport.h"
#include "TPZDarcyMemory.h"

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "TPZMatWithMem.h"


class TPZDarcyFlowWithMem : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TPZDarcyMemory> > {
    
    using TBase=TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TPZDarcyMemory> >;
    
    std::vector<REAL> m_gravity;
    /// Dimension
    int m_dimension;
    
    /// Scale factor for pressure
    STATE m_scale_pressure = 1.0;//1.e-6;

    /// Scale factor for flux variable
    STATE m_scale_flux = 1.0;
    
    /// Directive that stands for the use of four approximations spaces (iterative method)
    bool m_is_four_spaces_Q;
    
    TPZAlgebraicTransport *fAlgebraicTransport;
    
public:
    TMRSDataTransfer mSimData;
    /// Default constructor
    TPZDarcyFlowWithMem();
    
    /// Constructor based on a material id
    TPZDarcyFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TPZDarcyFlowWithMem(const TPZDarcyFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TPZDarcyFlowWithMem &operator=(const TPZDarcyFlowWithMem & other);
    
    /// Default destructor
    ~TPZDarcyFlowWithMem();
    
    void SetAlgebraicTransport(TPZAlgebraicTransport *algebraicTransport){
        fAlgebraicTransport=algebraicTransport;
    }
    TPZAlgebraicTransport *GetAlgebraicTransport(){
        return fAlgebraicTransport;
    }
    
    void SetGravity(std::vector<REAL> gravity){
        m_gravity =gravity;
    }
    std::vector<REAL> GetGravity( ){
        return m_gravity;
    }
    /// Set the required data at each integration point
    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Returns the name of the material
    std::string Name() const override {
        return "TPZDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TPZDarcyFlowWithMem(*this);
    }
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) const override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) const override;
    
    /// returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var) const override;
    
   
    
    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    
    // Contribute Methods being used
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)override;
    
    void ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef)override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)override;
    
};

#endif /* TPZDarcyFlowWithMem_h */
