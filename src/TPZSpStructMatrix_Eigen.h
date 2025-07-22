//
//  TPZSpStructMatrix_Eeigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//

#ifndef TPZSpStructMatrix_Eeigen_h
#define TPZSpStructMatrix_Eeigen_h

#include <stdio.h>
//#include "pzysmp.h"
#include "TPZYSMPPardiso.h"

#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"
//#include "pzstrmatrix.h"


/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSpStructMatrixEigen : public TPZStructMatrixT<STATE>,
        public TPZStructMatrixOR<STATE>
{
public:
    using TPZStructMatrixT<STATE>::TPZStructMatrixT;
    virtual TPZStructMatrix * Clone() override;
    virtual TPZMatrix<STATE> * Create() override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

//using TPZStructMatrix::CreateAssemble;
//    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
//    virtual TPZStructMatrix * Clone() override;
    
public:
    int ClassId() const override;
    
    
    virtual TPZSpStructMatrixEigen *NewMatrix() const{
        return new TPZSpStructMatrixEigen();
    }
    private :
    TPZSpStructMatrixEigen();
    
    friend TPZPersistenceManager;
};
#endif /* TPZSpStructMatrix_Eeigen_hpp */
