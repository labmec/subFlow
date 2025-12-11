//
//  TPZFastCondensedElement.h
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//

#pragma once

#include "TPZElementMatrixT.h"
#include "pzcmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzelmat.h"
#include "pzmultiphysicscompel.h"

class TPZFastCondensedElement : public TPZCondensedCompElT<REAL> {

protected:
  // this will be the multiplying factor for the condensed stiffness matrix K11
  TPZCompEl *fcelmultiphysics;
  bool fIsGroup = false;
  REAL fLambda = 1.0;
  REAL fSw = 1.0;
  REAL fSwLast = 1.0;
  REAL fPressureLastState = 0.0;
  TPZFNMatrix<9, REAL> fPermeabilityTensor;
  TPZFNMatrix<9, REAL> fInvPerm;

  // mixture density = rhow*fw + rhoo*fo
  REAL fMixedDensity = 1.;
  REAL fCompressibilityMatrixTerm = 1.0;
  REAL fCompressibiilityRhsTerm = 1.0;
  // this constant contains the source term
  REAL fSource = 0.;
  // this flag indicates whether the matrix of the father element has been computed
  bool fMatrixComputed = false;

  // reference stifness matrix and rhs
  TPZElementMatrixT<STATE> fEK, fEF;

  // vector of pressure connects
  TPZManVector<int64_t, 85> fPressureConnects;
  // vector of flux connects
  TPZManVector<int64_t> fFluxConnects;
  // connect index of the average pressure
  int64_t fAveragePressureConnect;

  // value of pressure values for unit body force
  TPZManVector<STATE, 100> fBodyForcePressureRef;
  // value of flux values for unit body force
  TPZManVector<STATE, 100> fBodyForceFluxRef;

  // value of pressures for constant unit pressure
  TPZManVector<STATE, 100> fConstantUnitPressure;

public:
  static bool fSkipLoadSolution;

  bool hasIndexes = false;

  std::vector<int> faIndexK00;
  std::vector<int> faIndexK01;
  std::vector<int> faIndexK11;

  std::vector<int> faValK00;
  std::vector<int> faValK01;
  std::vector<int> faValK11;

protected:
  // extract the solution vector of the condensed element
  void GetSolutionVector(TPZFMatrix<STATE> &solvec);

  // global indices of the pressure equations
  void PressureEquations(TPZVec<int64_t> &eqs);

  // global indices of the pressure equations
  void InternalFluxEquations(TPZVec<int64_t> &eqs);

  // global index of the average pressure equation
  int64_t AveragePressureEquation();

  // global indices of the boundary flux equations
  void BoundaryFluxEquations(TPZVec<int64_t> &eqs);

  // adjust the multiplying coeficients of the pressure equations
  void AdjustPressureCoefficients();

  /// compute internal coeficients as a function of the average pressure and boundary fluxes
  void ComputeInternalCoefficients();

  // Identify the connects and associated equations
  void IdentifyConnectandEquations();

  // Identify the condensed elements in this structure
  void FindCondensed(TPZStack<TPZCondensedCompElT *> &condensedelements);

public:
  TPZFastCondensedElement(TPZCompEl *ref, bool keepmatrix = false);

  /// Assignement constructor
  const TPZFastCondensedElement &operator=(const TPZFastCondensedElement &other) {
    fLambda = other.fLambda;
    return *this;
  }

  /** @brief create a copy of the condensed computational element in the other mesh */
  TPZFastCondensedElement(const TPZFastCondensedElement &copy, TPZCompMesh &mesh) : TPZCondensedCompElT(copy, mesh) {
    fEK = copy.fEK;
    fEF = copy.fEF;
    fLambda = copy.fLambda;
    fMatrixComputed = copy.fMatrixComputed;
  }

  void SetMultiphysics(TPZCompEl *celmult) {
    fcelmultiphysics = celmult;
  }
  TPZCompEl *GetMultiphysics() {
    return fcelmultiphysics;
  }
  void SetIsGroup(bool isgroup) {
    fIsGroup = isgroup;
  }
  int GetIsGroup() {
    return fIsGroup;
  }

  virtual ~TPZFastCondensedElement() {
  }

  /// compute the body force reference values
  void ComputeBodyforceRefValues();

  /// compute pressure equation values with respect to a constant pressure
  /// this will zero the body forces of the condensed elements
  void ComputeConstantPressureValues();
  /**
   * @brief Computes the element stifness matrix and right hand side
   * @param ek element stiffness matrix
   * @param ef element load vector
   */
  virtual void CalcStiff(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef) override;

  /**
   * @brief Computes the element right hand side
   * @param ef element load vector(s)
   */
  virtual void CalcResidual(TPZElementMatrixT<STATE> &ef) override;

  /** @brief Loads the solution within the internal data structure of the element */
  /**
   * Is used to initialize the solution of connect objects with dependency \n
   * Is also used to load the solution within SuperElements
   */
  virtual void LoadSolution() override;

  void SetLambda(REAL lambda);
  REAL GetLambda();

  void SetSw(REAL sw);
  REAL GetSw();

  void SetSwLast(REAL swlast);
  REAL GetSwLast();

  void SetPressureLastState(REAL pressurelast);
  REAL GetPressureLastState();

  void SetMixedDensity(REAL density);
  void SetCompressibiilityTerm(REAL matrix, REAL rhs);
  REAL GetMixedDensity();

  /**
   * @brief Calculates the solution - sol - for the variable var
   * at point qsi, where qsi is expressed in terms of the
   * master element coordinates
   * @param qsi master element coordinate
   * @param var variable name
   * @param sol vetor for the solution
   */
  virtual void Solution(TPZVec<REAL> &qsi, int var, TPZVec<STATE> &sol) override;

  void SetPermTensorAndInv(TPZFNMatrix<9, REAL> &PermeabilityTensor, TPZFNMatrix<9, REAL> &InvPerm);
  TPZFMatrix<REAL> &GetPermTensor();
  void SetConnectIndex(int inode, int64_t index);
};