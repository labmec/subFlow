//
//  TPZSpMatrixEigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//

#ifndef TPZSpMatrixEigen_h
#define TPZSpMatrixEigen_h

#include <stdio.h>
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TMRSDataTransfer.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZAlgebraicTransport.h"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
//#include <Eigen/PardisoSupport>
//#include <Eigen/SuperLUSupport>
#include "TPZAnalysisAuxEigen.h"
#include "TPZSpMatrixEigen.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

/*
template<typename StorageIndex = typename Eigen::SparseMatrix<REAL>::StorageIndex >
class Eigen::Triplet
{
public:
    
    Eigen::Triplet() : m_row(0), m_col(0), m_value(0) {}
    
    Eigen::Triplet(const StorageIndex& i, const StorageIndex& j, const REAL& v = REAL(0))
    : m_row(i), m_col(j),m_value(v)
    {}
    
    const StorageIndex& row() const { return m_row; }
    const StorageIndex& col() const { return m_col; }
    const REAL & value() const { return m_value; }
    
protected:
    
    StorageIndex m_row, m_col;
    REAL m_value;
};
*/
//
    template<class TVar>
    class TPZSpMatrixEigen : public TPZMatrix<TVar> {
    
        
    public:
        
        TPZSpMatrixEigen();
        TPZSpMatrixEigen(const int64_t rows,const int64_t cols );
        TPZSpMatrixEigen(const TPZSpMatrixEigen &cp);
        TPZSpMatrixEigen &operator=(const TPZSpMatrixEigen<TVar> &copy);
        CLONEDEF(TPZSpMatrixEigen)
        
        virtual ~TPZSpMatrixEigen();
        
       
        
        /** @brief Fill matrix storage with randomic values */
        /** This method use GetVal and PutVal which are implemented by each type matrices */
        void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
        
        inline TPZSpMatrixEigen<TVar>* NewMatrix() const override{
            return new TPZSpMatrixEigen<TVar>();
        }
        
        /** @brief Get the matrix entry at (row,col) without bound checking */
        virtual const TVar GetVal(const int64_t row,const int64_t col ) const override;
        
        
        bool isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col);

        
        int PutVal(const int64_t row, const int64_t col, const TVar &Value) override;
        
        virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                             const TVar alpha=1.,const TVar beta = 0., const int opt = 0) const override;
        

        
        virtual int GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
                           const int64_t colSize, TPZFMatrix<TVar> & A ) const override;
        
        
        /** @brief Pass the data to the class. */
        virtual void SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A );
        
        /** @brief Print the matrix along with a identification title */
        virtual void Print(const char *title, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
        
        
        virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex) override;
        
        virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex) override;
        
        
        virtual int Zero() override;
        
        /**
         * @name Factorization
         * @brief Those member functions are related to matrices factorization
         */
        //@{
        /**
         * @brief Decomposes the current matrix using LU decomposition.
         */
        
        virtual int Decompose_LU(std::list<int64_t> &singular);
        virtual int Decompose_LU();

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
        
        int Solve_LU( TPZFMatrix<TVar>*B, std::list<int64_t> &singular) {
            return ( ( !Decompose_LU(singular) )?  0 : Substitution( B )  );
        }

        int Solve_LU( TPZFMatrix<TVar>*B ) {
            return ( ( !Decompose_LU() )?  0 : Substitution( B )  );
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
                default:
                    DebugStop();
                    break;
            }
            return 0;
        }
        int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override
        {
            if(this->fDecomposed != dt) DebugStop();
            switch ( dt ) {
                case ELU:
                    return( Substitution( &F)  );
                default:
                    DebugStop();
                    break;
            }
            return 0;
        }

        
        //@}
        
        /**
         * @name Substitutions
         * @brief Substitutions forward and backward
         */
        //@{
        /**
         * @brief Computes Forward and Backward substitution for a "LU" decomposed matrix.
         * @param B right hand side and result after all
         */
        virtual int Substitution( TPZFMatrix<TVar> * B ) const;
        
        void FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat) const;
       
        void FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const;
        //@}
        
    public:
        int ClassId() const override;
        Eigen::SparseMatrix<REAL> fsparse_eigen;
        mutable std::vector<Eigen::Triplet<REAL> > ftriplets;
//        mutable Eigen::PardisoLU<Eigen::SparseMatrix<REAL>> m_analysis;
//        mutable Eigen::SparseLU<Eigen::SparseMatrix<REAL>> m_analysis;
        mutable Eigen::SparseLU<Eigen::SparseMatrix<REAL>> m_analysis;
        protected:
        int   fSymmetric;
    
        void InitializeData();
        
        /** @brief Number of entries storaged in the Matrix*/
        virtual int64_t Size() const override{
            DebugStop();
            return 1;
        }
        /** @{ */
        /** @brief Pointer to the beginning of the storage of the matrix*/
        virtual TVar* &Elem() override{
            DebugStop();
        }
        virtual const TVar* Elem() const override{
            DebugStop();
            return nullptr;
        }
        virtual void CopyFrom(const TPZMatrix<TVar> *mat) override{
            DebugStop();
        }
    };
    

    
#endif


