//
//  TPZSparceMatrixEigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/17/20.
//

#ifndef TPZSparceMatrixEigen_h
#define TPZSparceMatrixEigen_h

#include <stdio.h>
#include "TPZSpStructMatrix.h"

#include "pzgeoelbc.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

#include "TPZLinearAnalysis.h"

////#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzdxmesh.h"
#include <fstream>

#include "pzelmat.h"

//#include "pzysmp.h"
#include "TPZYSMPPardiso.h"
#include "TPZBndCondT.h"
#include "TPZAlgebraicTransport.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "Eigen/IterativeLinearSolvers"

//#include <Eigen/PardisoSupport>

//typedef Triplet2 Eigen::Triplet;
/*
template<typename StorageIndex = typename Eigen::SparseMatrix<REAL>::StorageIndex >
class Triplet2
{
public:
    
   Triplet2() : m_row(0), m_col(0), m_value(0) {}
    
   Triplet2(const StorageIndex& i, const StorageIndex& j, const REAL& v = REAL(0))
    : m_row(i), m_col(j), m_value(v)
    {}
    
    const StorageIndex& row() const { return m_row; }
    const StorageIndex& col() const { return m_col; }
    const REAL & value() const { return m_value; }
    
protected:
    
    StorageIndex m_row, m_col;
    
    REAL m_value;
};
*/
//using namespace Eigen;
class TPZAnalysisAuxEigen {
private:
    
    bool isFirst = false;
    TPZAlgebraicTransport *fAlgebraicTransport;
    Eigen::SparseMatrix<REAL> m_mass;
    Eigen::SparseMatrix<REAL> m_transmissibility;
    Eigen::SparseMatrix<REAL> m_rhs;
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> m_ds;
    
    std::vector<Eigen::Triplet<REAL> >           m_trans_triplets;
    std::vector<Eigen::Triplet<REAL> >           m_rhs_triplets;
    std::vector<Eigen::Triplet<REAL> >           m_mass_triplets;
//    Eigen::PardisoLU<Eigen::SparseMatrix<REAL>>  m_analysis;
    Eigen::SparseLU<Eigen::SparseMatrix<REAL>>  m_analysis;
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<REAL, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> m_analysis2;

    
public:
    TPZAnalysisAuxEigen(){
        
    }
    TPZAnalysisAuxEigen(TPZAlgebraicTransport *algebraicTransport){
        fAlgebraicTransport = algebraicTransport;

    }
    void SetAlgebraicTransport(TPZAlgebraicTransport *algebraicTransport){
        fAlgebraicTransport = algebraicTransport;
    }
    void AssembleMass();
    void Assemble();
    void AssembleResidual();
    int  NRhsRows();
    Eigen::SparseMatrix<REAL> Rhs();
    REAL RhsNorm(){
        return m_rhs.norm();
    }
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> Solution(){
        return m_ds;
    }
    
    void AnalyzePattern();
    void Solve();
    
};


#endif /* TPZSparceMatrixEigen_hpp */
