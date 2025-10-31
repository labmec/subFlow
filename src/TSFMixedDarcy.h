//
// Created by Giovane Avancini on 01/09/2025
//

#pragma once

#include "Material/DarcyFlow/TPZMixedDarcyFlow.h"

class TSFMixedDarcy : public TPZMixedDarcyFlow {

  using TBase = TPZMixedDarcyFlow;

public:
  /**
   * @brief Default constructor
   */
  TSFMixedDarcy();

  /**
   * @brief Class constructor
   * @param [in] id material id
   * @param [in] dim problem dimension
   */
  [[maybe_unused]] TSFMixedDarcy(int id, int dim);

  /**
           copy constructor
   */
  TSFMixedDarcy(const TSFMixedDarcy &copy);
  /**
           assignment operator
   */
  TSFMixedDarcy &operator=(const TSFMixedDarcy &copy);
  /**
   * @brief Returns a 'std::string' with the name of the material
   */
  [[nodiscard]] std::string Name() const override { return "TSFMixedDarcy"; }

  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
   * @param[in] datavec stores all input data
   * @param[in] weight is the weight of the integration rule
   * @param[out] ek is the element matrix
   * @param[out] ef is the rhs vector
   */
  void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

  /**
   * @brief It computes a contribution of additional spaces to the stiffness matrix and load vector at one integration point to allow static condensation
   * @param[in] datavec stores all input data
   * @param[in] weight is the weight of the integration rule
   * @param[out] ek is the element matrix
   * @param[out] ef is the rhs vector
   */
  void ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
   * @param[in] datavec stores all input data
   * @param[in] weight is the weight of the integration rule
   * @param[out] ek is the element matrix
   * @param[out] ef is the rhs vector
   * @param[in] bc is the boundary condition material
   */
  void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

  /*
   * @brief Fill requirements for volumetric contribute
   */
  void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

  /*
   * @brief Fill requirements for boundary contribute
   */
  void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

  /// Set the axisymmetry flag
  void SetAxisymmetry(bool IsAxisymmetric) { fIsAxisymmetric = IsAxisymmetric; }

  /// Returns the axisymmetry flag
  bool IsAxisymmetric() const { return fIsAxisymmetric; }

  /// Set gravity vector
  void SetGravity(const TPZFNMatrix<3, REAL> &gravity) { fGravity = gravity; }

  /// Get gravity vector
  const TPZFNMatrix<3, REAL> &GetGravity() const { return fGravity; }

  /// Enable/Disable four spaces formulation
  void SetFourSpaces(bool fourSpaces) { fFourSpaces = fourSpaces; }

  /// Returns the four spaces flag
  bool IsFourSpaces() const { return fFourSpaces; }

  /**
   * @brief Returns the solution associated with the variable var
   * @param[in] datavec stores all input data
   * @param[in] var index of the post-processing variable, according to VariableIndex method.
   * @param[out] Solout solution associated with the variable var
   */
  void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

protected:
  bool fIsAxisymmetric;
  bool fFourSpaces;
  TPZFNMatrix<3, REAL> fGravity;
};
