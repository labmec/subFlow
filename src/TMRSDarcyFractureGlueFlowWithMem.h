//
//  RSDarcyGlueFractureFlowWithMem.h
//
//  Created by José Villegas on 30/05/22.
//

#ifndef TMRSDarcyFractureGlueFlowWithMem_h
#define TMRSDarcyFractureGlueFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZMaterialDataT.h"

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "TPZMatWithMem.h"

struct TGlueMem
{
    TPZManVector<REAL,3> m_xco;
    std::pair<int,int> m_fracs;
    REAL m_dist;
    REAL m_flux;
    REAL m_dp;
    
    /// Class name
    const std::string Name() const
    {
        return "TGlueMem";
    }
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const
    {
        DebugStop();
    }
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context)
    {
        DebugStop();
    }
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const
    {
        out << "xco = " << m_xco << " dist " << m_dist << " frac ids " << m_fracs.first << " " << m_fracs.second;
    }
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TGlueMem & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const
    {
        return Hash("TGlueMem");

    }


};


class TMRSDarcyFractureGlueFlowWithMem : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TGlueMem> >{
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TGlueMem> >;

protected:
    
    REAL m_permeability;
public:
    
    /// Default constructor
    TMRSDarcyFractureGlueFlowWithMem();
    
    /// Constructor based on a material id
    TMRSDarcyFractureGlueFlowWithMem(int mat_id, REAL permeability);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureGlueFlowWithMem(const TMRSDarcyFractureGlueFlowWithMem & other) = default;
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureGlueFlowWithMem &operator=(const TMRSDarcyFractureGlueFlowWithMem & other) = default;
    
    /// Default destructor
    ~TMRSDarcyFractureGlueFlowWithMem();
    
    /// Set the required data at each integration point
    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec)const override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Returns the name of the material
    std::string Name() const override {
        return "TMRSDarcyFractureGlueFlowWithMem";
    }
    
    int Dimension() const override
    {
        return 2;
    }
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TMRSDarcyFractureGlueFlowWithMem(*this);
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

#endif /* TMRSDarcyFractureGlueFlowWithMem_h */
