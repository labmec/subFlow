//
//  TPZSSpMatrixEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#include "TPZSSpMatrixEigen.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include<Eigen/SparseCholesky>

#include "boost/date_time/posix_time/posix_time.hpp"
/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrixEigen methods.
 */

#include <memory.h>

#include "TPZSYSMPMatrix.h"
//#include "pzsysmp.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.eigen"));
#endif

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrixEigen<TVar>::TPZSYsmpMatrixEigen() : TPZRegisterClassId(&TPZSYsmpMatrixEigen::ClassId),
TPZMatrix<TVar>() {
//fanalysis.pardisoParameterArray()[59] = 1;
}

template<class TVar>
TPZSYsmpMatrixEigen<TVar>::TPZSYsmpMatrixEigen(const int64_t rows,const int64_t cols ) : TPZRegisterClassId(&TPZSYsmpMatrixEigen::ClassId),
TPZMatrix<TVar>(rows,cols) {
    Eigen::SparseMatrix<REAL> sparse_eigen(rows, cols);
    fsparse_eigen = sparse_eigen;
    fsparse_eigen.setZero();
//    fanalysis.pardisoParameterArray()[59] = 1;
}

template<class TVar>
TPZSYsmpMatrixEigen<TVar>::~TPZSYsmpMatrixEigen() {
    // Deletes everything associated with a TPZSYsmpMatrixEigen
#ifdef DESTRUCTOR
    cerr << "~TPZSYsmpMatrixEigen()\n";
#endif
}

template<class TVar>
TPZSYsmpMatrixEigen<TVar> &TPZSYsmpMatrixEigen<TVar>::operator=(const TPZSYsmpMatrixEigen<TVar> &copy)
{
    TPZMatrix<TVar>::operator=(copy);
    fsparse_eigen = copy.fsparse_eigen;
    fIA =copy.fIA;
    fJA = copy.fJA;
    fA = copy.fA;
    fDiag = copy.fDiag;

    return *this;
}


template<class TVar>
const TVar TPZSYsmpMatrixEigen<TVar>::GetVal(const int64_t r,const int64_t c ) const {
    // Get the matrix entry at (row,col) without bound checking
    
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col ) return fsparse_eigen.valuePtr()[ic];
    }
    return this->gZero;
  
   
}

/** @brief Put values without bounds checking \n
 *  This method is faster than "Put" if DEBUG is defined.
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & val )
{
   
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }

    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col )
        {
            fsparse_eigen.valuePtr()[ic] = val;
            return 0;
        }
    }
    if (val != (TVar(0.))) {
        DebugStop();
    }
    return 0;
 
}

template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & val, int &posfa )
{
    
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
    
    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col )
        {
            fsparse_eigen.valuePtr()[ic] = val;
            posfa = ic;
            return 0;
        }
    }
    if (val != (TVar(0.))) {
        DebugStop();
    }
    return 0;
    
}

template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::PutVal(int posfa, const TVar & val )
{
    fsparse_eigen.valuePtr()[posfa] = val;
    return 0;
    
    
}

template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::Print(const char *title, std::ostream &out ,const MatrixOutputFormat form) const {
    // Print the matrix along with a identification title
    // Print the matrix along with a identification title
    if(form != EInputFormat) {
        out << "\nSparse Eigen Matrix Print: " << title << '\n';
        out<<fsparse_eigen<<'\n';
    }
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::ComputeDiagonal() {
    if(!fDiag.size()) fDiag.resize(this->Rows());
    for(int ir=0; ir<this->Rows(); ir++) {
        fDiag[ir] = fsparse_eigen.coeff(ir,ir);
    }
}
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex){
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
/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
{
    if (!symmetric || nrow != ncol) {
        DebugStop();
    }
    TPZFMatrix<TVar> orig;
    SymProp sprop = SymProp::NonSym;
    if(symmetric) sprop = SymProp::Sym;
    orig.AutoFill(nrow,ncol,sprop);
    
    TPZVec<int64_t> IA(nrow+1);
    TPZStack<int64_t> JA;
    TPZStack<TVar> A;
    IA[0] = 0;
    TPZVec<std::set<int64_t> > eqs(nrow);
    for (int64_t row=0; row<nrow; row++) {
        eqs[row].insert(row);
        for (int64_t col = 0; col<ncol; col++) {
            REAL test = rand()*1./RAND_MAX;
            if (test > 0.5) {
                eqs[row].insert(col);
                if (symmetric) {
                    eqs[col].insert(row);
                }
            }
        }
    }
   
    for (int64_t row=0; row< nrow; row++) {
        for (std::set<int64_t>::iterator col = eqs[row].begin(); col != eqs[row].end(); col++) {
            if(*col >= row)
            {
                JA.Push(*col);
                A.Push(orig(row,*col));
            }
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}


/**
 * @name Factorization
 * @brief Those member functions perform the matrix factorization
 * @{
 */


/**
 * @brief Decomposes the current matrix using LDLt. \n
 * The current matrix has to be symmetric.
 * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    Decompose_LDLt();
    return 1;
}


/** @brief Decomposes the current matrix using Cholesky method. The current matrix has to be symmetric. */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_Cholesky()
{
    if(this->IsDecomposed() == ECholesky) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    this->SetIsDecomposed(ECholesky);
    return 1;
}
/**
 * @brief Decomposes the current matrix using Cholesky method.
 * @param singular
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    return Decompose_Cholesky();
}



/** @} */
/** @brief Decomposes the current matrix using LDLt. */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_LDLt()
{
    
    if(this->IsDecomposed() == ELDLt) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout << "Eigen control parameters\n";
//        for(int i=0; i<64; i++) sout << fanalysis.pardisoParameterArray()[i] << ' ';
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
    fanalysis.factorize(fsparse_eigen);
//    if(!fanalysis.m_factorizationIsOk){
//        std::cout << "\n ==> ERROR! Could not factorize matrix!" << std::endl;
//    }
//#ifdef LOG4CXX
//    if(logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout << "Eigen control parameters\n";
//        for(int i=0; i<64; i++) sout << fanalysis.pardisoParameterArray()[i] << ' ';
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
    this->SetIsDecomposed(ELDLt);
    return 1;
    
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    //    fanalysis.compute(fsparse_eigen);
    
    //    fanalysis.pardisoParameterArray()[59] = 1;
    
    TPZFMatrix<TVar> x(*b);
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> rhs;

    FromPZtoEigen(x, rhs);
    
    Eigen::SparseMatrix<REAL> sparse_eigen_c = fsparse_eigen + Eigen::SparseMatrix<REAL>(fsparse_eigen.transpose());
     for(int i=0; i<fsparse_eigen.innerSize(); i++){
         sparse_eigen_c.coeffRef(i, i) *= 0.5;
    }

//    Eigen::PardisoLDLT<Eigen::SparseMatrix<REAL>> solverpar;
//    Eigen::SparseLU<Eigen::SparseMatrix<REAL>> solverpar;
//    std::cout<<"NEQSUB= "<<fsparse_eigen.innerSize()<<std::endl;
//    std::cout<<"******"<<std::endl;
//    std::cout<<fsparse_eigen.toDense()<<std::endl;
//    std::cout<<"******"<<std::endl;
    fanalysis.compute(sparse_eigen_c);
//    std::cout<< sparse_eigen_c.toDense()<< std::endl;
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> dsolS= fanalysis.solve(rhs);
    
    this->FromEigentoPZ(x, dsolS);
//    x.Print("SolEeigen=",std::cout, EMathematicaInput);
    *b = x;
    return 1;
}

//template<class TVar>
//void TPZSYsmpMatrixEigen<TVar>::SetFromTriplets(){
//    fsparse_eigen.setFromTriplets(m_triplets.begin(), m_triplets.end());
//}
//
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_LU()
{
    if(this->IsDecomposed() == ELU) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    //
    fanalysis.analyzePattern(fsparse_eigen);
    fanalysis.factorize(fsparse_eigen);
    //
    this->SetIsDecomposed(ELU);
    return 1;
}
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::SetFromTriplets(){
    if(hastriplets==1){
        fsparse_eigen.setZero();
        fsparse_eigen.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
    }
}

template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Substitution( TPZFMatrix<TVar> *B ) const
{
    TPZFMatrix<TVar> x(*B);
    int nrows = x.Rows();
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> rhs(nrows,1);
    this->FromPZtoEigen(x, rhs);
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> dsol = fanalysis.solve(rhs);
    this->FromEigentoPZ(x, dsol);
    *B = x;
    return 1;
}
//
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                        TPZFMatrix<TVar> &z,
                                        const TVar alpha,const TVar beta,const int opt) const {
    
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> xeigen,yeigen,zeigen;
    this->FromPZtoEigen(x, xeigen);
    this->FromPZtoEigen(y, yeigen);
    zeigen = beta*yeigen + alpha*xeigen;
    this->FromEigentoPZ(z, zeigen);
    
    
}
/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
    //    fPardisoControl.Solve(*b,x);
    DebugStop();
    *b = x;
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}



template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrixEigen") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template<class TVar>
bool TPZSYsmpMatrixEigen<TVar>::isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col)
{
    for (Eigen::SparseMatrix<REAL>::InnerIterator it(mat, col); it; ++it) {
        if (it.row() == row) return false;
    }
    return true;
}
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
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
void TPZSYsmpMatrixEigen<TVar>::FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
    int nrows = pzmat.Rows();
    int ncols = pzmat.Cols();
    if (nrows<0 || ncols<0) {
        DebugStop();
    }
    pzmat.Zero();
    //    eigenmat.resize(nrows, ncols);
    for (int i=0; i< nrows; i++) {
        for (int j=0; j< ncols; j++) {
            pzmat(i, j)= eigenmat(i, j) ;
        };
    }
}
template class TPZSYsmpMatrixEigen<double>;
template class TPZSYsmpMatrixEigen<float>;
template class TPZSYsmpMatrixEigen<long double>;
