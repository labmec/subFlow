//
//  TPZSSpMatrixEigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#ifndef TPZSSpMatrixEigen_hpp
#define TPZSSpMatrixEigen_hpp
#define EIGEN_SUPERLU_SUPPORT

#include <stdio.h>
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZSSpMatrixEigen.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include<Eigen/SparseCholesky>
//#include <Eigen/SparseLDLt>
//#include <Eigen/PardisoSupport>

//#include <Eigen/SuperLUSupport>

//#include <Eigen/>
#include "TPZAnalysisAuxEigen.h"
#include "TPZSpMatrixEigen.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include<Eigen/SparseCholesky>

#ifdef PZ_USING_MKL
#include "TPZPardisoSolver.h"
#endif
/**
 * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
 * @ingroup matrix
 */


template<class TVar>
class TPZSYsmpMatrixEigen : public TPZMatrix<TVar>{
    
#ifdef PZ_USING_MKL
    friend class TPZPardisoSolver<TVar>;
#endif
    
    public :
    /** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrixEigen();
    /** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrixEigen(const int64_t rows, const int64_t cols );
    /** @brief Copy constructor */
    TPZSYsmpMatrixEigen(const TPZSYsmpMatrixEigen &cp){
 
    }
//    TPZSYsmpMatrixEigen(const TPZSYsmpMatrixEigen<TVar> &cp) :
//    TPZRegisterClassId(&TPZSYsmpMatrixEigen::ClassId),
//    TPZMatrix<TVar>(cp), fIA(cp.fIA), fJA(cp.fJA), fA(cp.fA), fDiag(cp.fDiag)
//    {
//        
//       
//    }
    
    TPZSYsmpMatrixEigen &operator=(const TPZSYsmpMatrixEigen<TVar> &copy);
    
    CLONEDEF(TPZSYsmpMatrixEigen)
    /** @brief Destructor */
    virtual ~TPZSYsmpMatrixEigen();
    
    /** @brief Checks if the current matrix is symmetric */
//    virtual int IsSimetric() const  override { return 1; }
    /** @brief Checks if current matrix is square */
    inline int IsSquare() const { return 1;}
    
    /** @brief To create a matrix of the same type */
    inline TPZSYsmpMatrixEigen<TVar>* NewMatrix() const override{
        return new TPZSYsmpMatrixEigen<TVar>{};
    }
    
    /** @brief Creates a copy from a given matrix of arbitrary storage format.
     Every implementation should check for type compatibility */
    virtual void CopyFrom(const TPZMatrix<TVar> *mat){
        DebugStop();
    }
    /** @brief Number of entries storaged in the Matrix*/
    virtual int64_t Size() const{
        DebugStop();
        return 1;
    }
    
    virtual TVar* &Elem(){
        DebugStop();
    }
    virtual const TVar* Elem() const{
        DebugStop();
        return nullptr;
    }
    
    /** @brief Zeroes the matrix */
    virtual int Zero() override {
        fsparse_eigen=0.0*fsparse_eigen;
        for (int ival=0; ival<fA.size(); ival++) {
           fsparse_eigen.valuePtr()[ival]=0;
        }
                TPZMatrix<TVar>::fDecomposed = ENoDecompose;
        
//        fA.Fill(0.);
//        fDiag.Fill(0.);
#ifndef PZ_USING_MKL
//        TPZMatrix<TVar>::fDecomposed = ENoDecompose;
#endif
        return 0;
        
    }
    
    /** @brief Zeroes the matrix */
    virtual int Redim(int64_t rows, int64_t cols) override
    {
        if(rows == this->fRow && cols == this->fCol)
        {
            Zero();
        }
        else
        {
            DebugStop();
        }
        return 0;
    }
    
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
    
    /** @brief Get the matrix entry at (row,col) without bound checking */
    virtual const TVar GetVal(const int64_t row, const int64_t col ) const override;
    
    /** @brief Put values without bounds checking \n
     *  This method is faster than "Put" if DEBUG is defined.
     */
    virtual int PutVal(const int64_t /*row*/,const int64_t /*col*/,const TVar & val ) override;
    int PutVal(const int64_t /*row*/,const int64_t /*col*/,const TVar & val, int &posfa );
    int PutVal(int posfa, const TVar & val );
    /** @brief Computes z = beta * y + alpha * opt(this)*x */
    /** @note z and x cannot overlap in memory */
    virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                         const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
    
    /** @brief Sets data to the class */
    virtual void SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A );
    
    /// Access function for the coefficients
    TPZVec<TVar> &A()
    {
        return fA;
    }
    
    TPZVec<int64_t> &IA()
    {
        return fIA;
    }
    
    TPZVec<int64_t> &JA()
    {
        return fJA;
    }
    
    /** @brief Print the matrix along with a identification title */
    virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const override;
    
#ifdef PZ_USING_MKL
    /**
     * @name Factorization
     * @brief Those member functions perform the matrix factorization
     * @{
     */
    
    /** @brief decompose the system of equations acording to the decomposition
      * scheme */
    virtual int Decompose(const DecomposeType dt) override {
        // Only allowing for LU
        switch (dt) {
            case ELU:
                return Decompose_LU();
                break;
            default:
                DebugStop();
                break;
        }
        return -1;
    }
    
    /**
     * @brief Solves the linear system using Direct methods
     * @param F The right hand side of the system and where the solution is stored.
     * @param dt Indicates type of decomposition
     */
    int SolveDirect( TPZFMatrix<TVar> &B , const DecomposeType dt) override {
        
        switch ( dt ) {
            case ELU:
                return( Solve_LU( &B)  );
            case ECholesky:
                return( Solve_Cholesky( &B )  );
            case ELDLt:
                return( Solve_LDLt( &B )  );
            default:
                DebugStop();
                break;
        }
        return ( 0 );
    }
    int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override
    {
        if(this->fDecomposed != dt) DebugStop();
        switch ( dt ) {
            case ELU:
                return( Substitution( &F)  );
            case ECholesky:
                return ( Subst_Forward(&F) && Subst_Backward(&F) );
            case ELDLt:
                return( Subst_LForward( &F ) && Subst_Diag( &F ) && Subst_LBackward( &F ) );
            default:
                DebugStop();
                break;
        }
        return ( 0 );
    }
    
    int Solve_LU( TPZFMatrix<TVar>*B ) {
        return ( ( !Decompose_LU() )?  0 : Substitution( B )  );
    }
    
    /**********************/
    /*** Solve Cholesky ***/
    //
    //  Se nao conseguir resolver por Cholesky retorna 0 e a matriz
    //   sera' modificada (seu valor perdera' o sentido).
    //
    int Solve_Cholesky( TPZFMatrix<TVar>* B )
    {
        return(
               ( !Decompose_Cholesky() )?  0 :( Subst_Forward( B ) && Subst_Backward( B ) )
               );
    }

    int Solve_Cholesky( TPZFMatrix<TVar>* B, std::list<int64_t> &singular ) {
        return(
               ( !Decompose_Cholesky(singular) )?  0 :( Subst_Forward( B ) && Subst_Backward( B ) )
               );
    }

    /******************/
    /*** Solve LDLt ***/

    int Solve_LDLt( TPZFMatrix<TVar>* B ) {
        
        return(
               ( !Decompose_LDLt() )? 0 :
               ( Subst_LForward( B ) && Subst_Diag( B ) && Subst_LBackward( B ) )
               );
    }

    /**
     * @brief Decomposes the current matrix using LDLt. \n
     * The current matrix has to be symmetric.
     * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
     */
    virtual int Decompose_LDLt(std::list<int64_t> &singular);
    /** @brief Decomposes the current matrix using LDLt. */
    virtual int Decompose_LDLt();
    
    /** @brief Decomposes the current matrix using Cholesky method. The current matrix has to be symmetric. */
    virtual int Decompose_Cholesky();
    /**
     * @brief Decomposes the current matrix using Cholesky method.
     * @param singular
     */
    virtual int Decompose_Cholesky(std::list<int64_t> &singular);
    
    virtual int Decompose_LU();
    virtual int Substitution( TPZFMatrix<TVar> * B ) const;
    
    /** @} */
    
    /**
     * @name Substitutions
     * @brief Substitutions forward and backward
     * @{
     */
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LForward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LBackward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
     * @param b right hand side and result after all
     */
    virtual int Subst_Diag( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Forward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Backward( TPZFMatrix<TVar>* b ) const;
    
    
    /** @} */
    
    
#endif
public:
    int ClassId() const override;
    
    void ComputeDiagonal();
    
    void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex) override;
    bool isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col);
    
   
    int hastriplets=0;
//    void SetFromTriplets(int ok) override;
    
    mutable std::vector<Eigen::Triplet<REAL> > m_triplets;
    mutable Eigen::SparseMatrix<REAL> fsparse_eigen;

    

//    mutable Eigen::PardisoLDLT<Eigen::SparseMatrix<REAL>, Eigen::Lower> fanalysis;
    mutable Eigen::SparseLU<Eigen::SparseMatrix<REAL>> fanalysis;
    //     mutable Eigen::SimplicialLDLT<Eigen::SparseMatrix<REAL>> fanalysis;


    void SetFromTriplets();
    
//    mutable Eigen::SuperLU<Eigen::SparseMatrix<REAL>> fanalysis;
//    Eigen::SuperLU<Eigen::SparseMatrix<double> > slu;
//    slu.compute(A);
//    x = slu.solve(b);
    void FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat) const;
    
    void FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const;
   
public:
    
    
    TPZVec<int64_t>  fIA;
    TPZVec<int64_t>  fJA;
    TPZVec<TVar> fA;
    
#ifdef PZ_USING_MKL
//    TPZPardisoControl<TVar> fPardisoControl;
#endif
    
    TPZVec<TVar> fDiag;
};

template<class TVar>
inline void TPZSYsmpMatrixEigen<TVar>::SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A )
{
    //
//    std::cout<<IA<<std::endl;
//        std::cout<<"*************"<<std::endl;
//    std::cout<<JA<<std::endl;
    fsparse_eigen.setZero();
    int ncols = IA.size()-1;
    fsparse_eigen.resize(ncols,ncols);
    std::vector<Eigen::Triplet<REAL> > triplets(2*(A.size()-ncols)+ncols);
    int count = 0;
    for (int icol = 0; icol < ncols; icol++) {
        
        for(int ivalk = IA[icol];ivalk<IA[icol+1]; ivalk++){
            REAL val = A[ivalk];
            int  c = icol;
            int i = JA[ivalk];
            Eigen::Triplet<REAL> trip(i,c,val);
            triplets[count]=trip;
             count++;           
        }
    }

    fsparse_eigen.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();
    fanalysis.analyzePattern(fsparse_eigen);

    fIA = IA;
    fJA = JA;
    fA  =  A;
    
    ComputeDiagonal();
    
   
}

#endif /* TPZSSpMatrixEigen_hpp */
