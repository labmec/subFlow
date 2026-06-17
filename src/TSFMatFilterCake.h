#pragma once

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include "TSFProblemData.h"

template <class TMEM>
class TSFMatFilterCake : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>, TPZMatErrorCombinedSpaces<STATE>> {

  using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TMEM>, TPZMatErrorCombinedSpaces<STATE>>;

public:
  /// @brief Default constructor
  TSFMatFilterCake();

  /// @brief Constructor with material id and dimension
  TSFMatFilterCake(int matid, int dimension);

  /// @brief Constructor with material id, dimension and filtercake properties
  TSFMatFilterCake(int matid, int dimension, TSFProblemData *simData);

  /// @brief Copy constructor
  TSFMatFilterCake(const TSFMatFilterCake &other);

  /// @brief Assignment operator
  /// @param other
  /// @return
  TSFMatFilterCake &operator=(const TSFMatFilterCake &other);

  /// @brief Destructor
  ~TSFMatFilterCake();

  /// @brief Returns an unique class identifier
  [[nodiscard]] int ClassId() const override;

  /// @brief Returns the name of the material
  std::string Name() const override { return "TSFMatFilterCake"; }

  /// @brief Print out the data associated with the material
  void Print(std::ostream &out = std::cout) const override;

  /// @brief Contribute Method
  void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

  /// @brief ContributeBC Method
  void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

  /// @brief Solution Method
  void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

  /// @brief Error Method
  void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

  /// @brief Set the required data at each integration point
  void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

  /// @brief Returns the dimension of the material
  int Dimension() const override { return fDimension; }

  /// @brief Returns the number of state variables associated with the material
  int NStateVariables() const override { return 1; }

  /// @brief Creates a new material of the same type
  /// @return Pointer to the new material
  virtual TPZMaterial *NewMaterial() const override { return new TSFMatFilterCake(*this); }

  /// @brief Returns the variable index associated with the name
  int VariableIndex(const std::string &name) const override;

  /// @brief Returns the number of variables associated with the variable indexed by var.
  int NSolutionVariables(int var) const override;

  /// @brief @ brief Set the alpha parameter (alpha = C / (1-phi)rho)
  void SetAlpha(REAL alpha) { fAlpha = alpha; }

  /// @brief @ brief Get the alpha parameter
  REAL GetAlpha() const { return fAlpha; }

  /// @brief @ brief Set the permeability according to the Carman-Kozeny equation
  void SetPerm(REAL perm) { fPerm = perm; }

  /// @brief @ brief Get the permeability according to the Carman-Kozeny equation
  REAL GetPerm() const { return fPerm; }

  /// @brief Set the axisymmetry flag
  void SetAxisymmetry(bool IsAxisymmetric) { fIsAxisymmetric = IsAxisymmetric; }

  /// @brief Simulation data
  void SetProblemData(TSFProblemData *simData) { fSimData = simData; }

protected:
  int fDimension;
  TSFProblemData *fSimData; // Simulation data
  REAL fAlpha;              // Coefficient related to the filter cake buildup, defined as alpha = C / (1-phi)rho
  REAL fPerm;               // Permeability according to the Carman-Kozeny equation
  bool fIsAxisymmetric;     // Whether the problem is axisymmetric
};