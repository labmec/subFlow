//
//  TPZSymetricSpStructMatrixEigenEigen.h
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#ifndef TPZSymetricSpStructMatrixEigenEigen_h
#define TPZSymetricSpStructMatrixEigenEigen_h

#include <stdio.h>
//#include "pzysmp.h"
#include "TPZYSMPPardiso.h"

#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpMatrixEigen.h"
#include "TPZFastCondensedElement.h"

#include "boost/date_time/posix_time/posix_time.hpp"

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSymetricSpStructMatrixEigen : public TPZStructMatrixT<STATE>,
       public TPZStructMatrixOR<STATE> {
    
           
    using TPZStructMatrixT<STATE>::TPZStructMatrixT;
    std::vector<Eigen::Triplet<REAL> > m_triplets;
    
public:

    TPZSymetricSpStructMatrixEigen(TPZCompMesh *cmesh) : TPZStructMatrixT<STATE>(cmesh)
       {
           
       }
    ~TPZSymetricSpStructMatrixEigen()
    {
        std::cout << "\nvirtual TPZSymetricSpStructMatrixEigen::~TPZSymetricSpStructMatrixEigen()" << std::endl;
    }
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration fAsTotalCalcStifSub = tsim2-tsim2;
    boost::posix_time::time_duration fAsTotalAdkelsSub = tsim2-tsim2;
    boost::posix_time::time_duration fAsTotaAssembleSub = tsim2-tsim2;
    boost::posix_time::time_duration fAsTotaCondensedSub = tsim2-tsim2;
    
    virtual TPZMatrix<STATE> * Create() override;
    
    virtual TPZMatrix<STATE> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
        
    virtual TPZStructMatrix * Clone() override;
    void Serial_Assemble(TPZBaseMatrix & stiffness, TPZBaseMatrix & rhs) override;
    void Serial_AssembleSub(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs);
    void Serial_AssembleGlob(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs);
    
    /** Used only for testing */
    static int main();
    
   //@{
   //!Read and Write methods
   int ClassId() const override;

   void Read(TPZStream& buf, void* context) override;

   void Write(TPZStream& buf, int withclassid) const override;
   //@}

    private :
    
 
    void Swap(int64_t *a, int64_t *b)
    {
        int64_t aux = *a;
        *a = *b;
        *b = aux;
    }
    
    friend TPZPersistenceManager;
};

#endif /* TPZSymetricSpStructMatrixEigenEigen_h */
