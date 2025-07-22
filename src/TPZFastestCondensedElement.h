//
//  TPZFastCondensedElement.h
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//

#ifndef TPZFASTESTCONDENSEDELEMENT_h
#define TPZFASTESTCONDENSEDELEMENT_h

#include "TPZFastCondensedElement.h"
#include "pzelmat.h"

/// This class implements a nonconforming H(div) element based on a constant
/// pressure mortar element
class TPZFastestCondensedElement : public TPZFastCondensedElement
{
    
protected:
    
    // this will be the multiplying factor for the condensed stiffness matrix K11
    
public:
    
    // extract the solution vector of the condensed element
    void GetSolutionVector(TPZFMatrix<STATE> &solvec);
    
    // global indices of the pressure equations
    void PressureEquations(TPZVec<int64_t> &eqs);
    
    // global index of the average pressure equation
    int64_t AveragePressureEquation();
    
    // global indices of the boundary flux equations
    void BoundaryFluxEquations(TPZVec<int64_t> &eqs);
    
    // adjust the multiplying coeficients of the pressure equations
    void AdjustPressureCoefficients();

public:
    
    TPZFastestCondensedElement(TPZCompEl *ref, bool keepmatrix = false) :
        TPZFastCondensedElement(ref,keepmatrix)
    {
       
        
    }
    
    
    /// Assignement constructor
    const TPZFastestCondensedElement & operator=(const TPZFastestCondensedElement & other){
        DebugStop();
        fLambda = other.fLambda;
        return *this;
    }
    
    /** @brief create a copy of the condensed computational element in the other mesh */
    TPZFastestCondensedElement(const TPZFastestCondensedElement &copy, TPZCompMesh &mesh) :
        TPZFastCondensedElement(copy, mesh)
    {
        fEK = copy.fEK;
        fEF = copy.fEF;
        fLambda = copy.fLambda;
        fMatrixComputed = copy.fMatrixComputed;
        DebugStop();
    }
    
    
    virtual ~TPZFastestCondensedElement()
    {
        
    }

    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override;
    
    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrixT<STATE> &ef) override;
    
    
};
#endif /* TPZFASTESTCONDENSEDELEMENT_h */
