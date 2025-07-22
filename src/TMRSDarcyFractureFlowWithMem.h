//
//  TMRSDarcyFlowWithMem.h
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSDarcyFractureFlowWithMem_h
#define TMRSDarcyFractureFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"
#include "TMRSDarcyFlowWithMem.h"
#include "TPZMaterialDataT.h"

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "TPZMatWithMem.h"
template <class TMEM>
class TMRSDarcyFractureFlowWithMem : public TMRSDarcyFlowWithMem<TMEM> {
    using TBase =TMRSDarcyFlowWithMem<TMEM>;
    
public:
    
    /// Default constructor
    TMRSDarcyFractureFlowWithMem();
    
    /// Constructor based on a material id
    TMRSDarcyFractureFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureFlowWithMem(const TMRSDarcyFractureFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureFlowWithMem &operator=(const TMRSDarcyFractureFlowWithMem & other);
    
    /// Default destructor
    ~TMRSDarcyFractureFlowWithMem();
    
    /// Set the required data at each integration point
    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec)const override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Returns the name of the material
    std::string Name() const override {
        return "TMRSDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return this->m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TMRSDarcyFractureFlowWithMem<TMEM>(*this);
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
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    void ContributeFourSpaces( const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)override;
    
};

#endif /* TMRSDarcyFlowWithMem_h */
