//
//  TPZSpMatrixEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//
// #include "pz_config.h"
#include "TPZSpMatrixEigen.h"
#include "TPZSSpMatrixEigen.h"
//#include "pzysmp.h"
#include "TPZYSMPPardiso.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <memory.h>
#include <string>
#include <map>
#include <pthread.h>

#include "tpzverysparsematrix.h"
// #include "pz_pthread.h"
#include "pzstack.h"

using namespace std;

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************
template<class TVar>
void TPZSpMatrixEigen<TVar>::InitializeData(){}



// ****************************************************************************
//
// Constructor
//
// ****************************************************************************

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen(const TPZSpMatrixEigen &cp){
    
}

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen() : TPZRegisterClassId(&TPZSpMatrixEigen::ClassId),
TPZMatrix<TVar>(),  fsparse_eigen(0,0)
{
     
}

template<class TVar>
TPZSpMatrixEigen<TVar> &TPZSpMatrixEigen<TVar>::operator=(const TPZSpMatrixEigen<TVar> &cp) {
    TPZMatrix<TVar>::operator=(cp);
    fsparse_eigen = cp.fsparse_eigen;
    return *this;
}


template<class TVar>
int TPZSpMatrixEigen<TVar>::PutVal(const int64_t row, const int64_t col, const TVar &Value){
//    std::cout<<fsparse_eigen.toDense()<<std::endl;
    if (!isNull(fsparse_eigen, row, col)) {
        fsparse_eigen.coeffRef(row,col)=Value;
    }
    else{
        std::cout<<"Non existing position on sparse matrix: line =" << row << " column =" << col << std::endl;
    }
    
    return 1;
}
template<class TVar>
void TPZSpMatrixEigen<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex){
    
    int64_t i,j;
    TVar value=0.;
    int64_t ipos,jpos;
    for(i=0;i<elmat.Rows();i++){
        for(j=0;j<elmat.Rows();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(i,j);
            if (!isNull(fsparse_eigen, ipos, jpos)) {
                Eigen::Triplet<REAL> triplet(ipos, jpos, value);
                ftriplets.push_back(triplet);
//                fsparse_eigen.coeffRef(ipos, jpos) += value;
            }
            else{
                std::cout<<"Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;
            }
        }
    }
}

template<class TVar>
void TPZSpMatrixEigen<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex){
    int64_t i,j = 0;
    TVar value=0.;
    int64_t ipos,jpos;
    for(i=0;i<sourceindex.NElements();i++){
        for(j=0;j<sourceindex.NElements();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(sourceindex[i],sourceindex[j]);
            if (!isNull(fsparse_eigen, ipos, jpos)) {
                fsparse_eigen.coeffRef(ipos, jpos) += value;
            }
            else{
                std::cout<<"Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;
            }
        }
    }
}

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen(const int64_t rows,const int64_t cols ) :
TPZRegisterClassId(&TPZSpMatrixEigen::ClassId),TPZMatrix<TVar>(rows,cols) {
  
    Eigen::SparseMatrix<REAL> sparse_eigen(rows, cols);
    fsparse_eigen = sparse_eigen;
    fsparse_eigen.setZero();
    fSymmetric = 0;
   
#ifdef CONSTRUCTOR
    cerr << "TPZSpMatrixEigen(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZSpMatrixEigen<TVar>::~TPZSpMatrixEigen() {
    // Deletes everything associated with a TPZSpMatrixEigen
#ifdef DESTRUCTOR
    cerr << "~TPZSpMatrixEigen()\n";
#endif
}


template<class TVar>
const TVar  TPZSpMatrixEigen<TVar>::GetVal(const int64_t row,const int64_t col ) const {
    TVar val = fsparse_eigen.coeff(row, col);
    return val;

}

template<class TVar>
void TPZSpMatrixEigen<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                   TPZFMatrix<TVar> &z,
                                   const TVar alpha,const TVar beta,const int opt) const {
    
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> xeigen,yeigen,zeigen;
    this->FromPZtoEigen(x, xeigen);
    this->FromPZtoEigen(y, yeigen);
    zeigen = beta*yeigen + alpha*xeigen;
    this->FromEigentoPZ(z, zeigen);

    
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZSpMatrixEigen<TVar>::Print(const char *title, ostream &out ,const MatrixOutputFormat form) const {
    // Print the matrix along with a identification title
    if(form != EInputFormat) {
        out << "\nSparse Eigen Matrix Print: " << title << '\n';
        out<<fsparse_eigen<<'\n';
    }
}




template<class TVar>
int TPZSpMatrixEigen<TVar>::Zero()
{
    fsparse_eigen = 0.0*fsparse_eigen;

    return 1;
}


template<class TVar>
int TPZSpMatrixEigen<TVar>::GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
                                 const int64_t colSize, TPZFMatrix<TVar> & A ) const {
    
    Eigen::SparseMatrix<REAL> sub = fsparse_eigen.block(sRow, sCol, rowSize, colSize);
    A.Resize(sub.innerSize(), sub.outerSize());
    A.Zero();
    for (int col=0; col<sub.outerSize(); col++) {
        for (Eigen::SparseMatrix<REAL>::InnerIterator it(sub, col); it; ++it) {
            A(it.row(),col)=sub.coeffRef(it.row(), col);
        }
    }

    return 0;
}


/*
 * Perform row update of the sparse matrix
 */



/**
 * Decomposes the current matrix using LU decomposition.
 */
template<class TVar>
int TPZSpMatrixEigen<TVar>::Decompose_LU(std::list<int64_t> &singular)
{
    return Decompose_LU();
}
template<class TVar>
int TPZSpMatrixEigen<TVar>::Decompose_LU()
{
    if(this->IsDecomposed() == ELU) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    
    m_analysis.analyzePattern(fsparse_eigen);
    m_analysis.factorize(fsparse_eigen);
    
    this->SetIsDecomposed(ELU);
    return 1;
}

template<class TVar>
int TPZSpMatrixEigen<TVar>::Substitution( TPZFMatrix<TVar> *B ) const
{
    TPZFMatrix<TVar> x(*B);
    int nrows = x.Rows();
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> rhs(nrows,1);
    this->FromPZtoEigen(x, rhs);
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> dsol = m_analysis.solve(rhs);
    this->FromEigentoPZ(x, dsol);
    *B = x;
    return 1;
}
template<class TVar>
bool TPZSpMatrixEigen<TVar>::isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col)
{
    for (Eigen::SparseMatrix<REAL>::InnerIterator it(mat, col); it; ++it) {
        if (it.row() == row) return false;
    }
    return true;
}

template<class TVar>
int TPZSpMatrixEigen<TVar>::ClassId() const{
    return Hash("TPZSpMatrixEigen") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template<class TVar>
void TPZSpMatrixEigen<TVar>::FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
    int nrows = pzmat.Rows();
    int ncols = pzmat.Cols();
    if (nrows<0 || ncols<0) {
        DebugStop();
    }
    eigenmat.resize(nrows, ncols);
    eigenmat.setZero();
    for (int i=0; i< nrows; i++) {
        for (int j=0; j< ncols; j++) {
            eigenmat(i, j) = pzmat.Get(i, j);
        }
    }
}
template<class TVar>
void TPZSpMatrixEigen<TVar>::FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
    int nrows = pzmat.Rows();
    int ncols = pzmat.Cols();
    if (nrows<0 || ncols<0) {
        DebugStop();
    }
    eigenmat.resize(nrows, ncols);
    for (int i=0; i< nrows; i++) {
        for (int j=0; j< ncols; j++) {
         pzmat(i, j)= eigenmat(i, j) ;
        };
    }
}
/** @brief Pass the data to the class. */
template<class TVar>
inline void TPZSpMatrixEigen<TVar>::SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A ){
    fsparse_eigen.setZero();
//    auto valores =fsparse_eigen.innerIndexPtr();
  
    std::vector<Eigen::Triplet<REAL> > triplets(A.size());
    int nrows = IA.size()-1;
    int count =0;
    for (int irow = 0; irow < nrows; irow++) {
        for(int k=IA[irow]; k<IA[irow+1]; k++){
            int row= irow;
            int col = JA[k];
            REAL val = A[k];
            Eigen::Triplet<REAL> trip(row, col, val);
            triplets[count] = trip;
            count++;
        }
    }
    fsparse_eigen.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();
    m_analysis.analyzePattern(fsparse_eigen);
    
    
    if (IA.size() != this->Rows() + 1 ) {
        DebugStop();
    }
    
    if (JA.size() != IA[this->Rows()]) {
        DebugStop();
    }
    
}
;
template class TPZSpMatrixEigen<long double>;
template class TPZSpMatrixEigen<double>;
template class TPZSpMatrixEigen<float>;
