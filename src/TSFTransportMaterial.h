#pragma once

#include "TPZBndCondT.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatInterfaceCombinedSpaces.h"
#include "TPZMatInterfaceSingleSpace.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatWithMem.h"
#include "TPZMaterialDataT.h"
#include <stdio.h>

class TSFTransportMaterial : public TPZMatBase<STATE, TPZMatInterfaceSingleSpace<STATE>, TPZMatCombinedSpacesT<STATE>, TPZMatSingleSpaceT<STATE>, TPZMatInterfaceCombinedSpaces<STATE>> {

  using TBase = TPZMatBase<STATE, TPZMatInterfaceSingleSpace<STATE>, TPZMatCombinedSpacesT<STATE>, TPZMatSingleSpaceT<STATE>, TPZMatInterfaceCombinedSpaces<STATE>>;

private:
  /** @brief material dimension */ // TODO:: Candidate for deletion
  int m_dimension;

  /** @brief material dimension */
  int m_mat_id;

  /** @brief Directive that stands for Mass matrix assembly  */
  bool m_mass_matrix_Q;

  /** @brief Regular time step size  */
  REAL m_dt;

  /** @brief Porosity  */
  REAL m_phi;

  REAL m_fracture_epsilon;

public:
  /** @brief Default constructor */
  TSFTransportMaterial();

  /** @brief Constructor based on a material id */
  TSFTransportMaterial(int matid, int dimension);

  /** @brief Constructor based on a TSFTransportMaterial object */
  TSFTransportMaterial(const TSFTransportMaterial &other);

  /** @brief Assignment operator */
  TSFTransportMaterial &operator=(const TSFTransportMaterial &other);

  /** @brief Default destructor */
  ~TSFTransportMaterial();

  /** @brief Set the required data at each integration point */
  virtual void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

  /** @brief Set the required data at each integration point */
  virtual void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

  virtual void FillDataRequirementsInterface(TPZMaterialData &data) const override;

  virtual void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec_left, std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;

  /** @brief Returns the name of the material */
  virtual std::string Name() const override {
    return "TSFTransportMaterial";
  }

  /** @brief Returns the integrable dimension of the material */
  int Dimension() const override { return m_dimension; }

  /** @brief Sets material dimension */
  void SetDimension(int dim) { m_dimension = dim; }

  /** @brief Returns the number of state variables associated with the material */
  int NStateVariables() const override { return 1; } // Deprecated, must to be removed

  /** @brief Returns material copied form this object */
  virtual TPZMaterial *NewMaterial() const override {
    return new TSFTransportMaterial(*this);
  }
  virtual void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override{
  }
  virtual void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override {
  }
  /** @brief Print out the data associated with the material */
  void Print(std::ostream &out = std::cout) const override;

  /** @brief Returns the variable index associated with the name */
  int VariableIndex(const std::string &name) const override;

  /** @brief Returns the number of variables associated with varindex */
  int NSolutionVariables(int var) const override;

  // Contribute Methods being used

  /** @brief Returns the solution associated with the var index */
  void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) override;

  void Solution(const TPZMaterialDataT<REAL> &data, int var, TPZVec<REAL> &sol) override;

  void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

  void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;

  virtual void ContributeInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &datavecleft, const std::map<int, TPZMaterialDataT<STATE>> &datavecright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

  void ContributeBCInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

  virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override {

    DebugStop();
  }

  /**
   * Unique identifier for serialization purposes
   */
  int ClassId() const override;

  /**
   * Save the element data to a stream
   */
  void Write(TPZStream &buf, int withclassid);

  /**
   * Read the element data from a stream
   */
  void Read(TPZStream &buf, void *context) override;

  /** @brief Set directive that stands for Mass matrix assembly  */
  void SetMassMatrixAssembly(bool mass_matrix_Q) {
    m_mass_matrix_Q = mass_matrix_Q;
  }

  /** @brief Set directive that stands for Mass matrix assembly  */
  bool GetMassMatrixAssembly() {
    return m_mass_matrix_Q;
  }

  /** @brief Set regular time step size  */
  void SetTimeStep(REAL dt) {
    m_dt = dt;
  }

  /** @brief Get regular time step size  */
  REAL GetTimeStep() {
    return m_dt;
  }

  /** @brief Set porosity  */
  void SetPorosity(REAL phi) {
    m_phi = phi;
  }

  /** @brief Get porosity  */
  REAL GetPorosity() {
    return m_phi;
  }

  /** @brief Set fracture cross length  */
  void SetFractureCrossLength(REAL fracture_epsilon) {
    m_fracture_epsilon = fracture_epsilon;
  }

  /** @brief Get fracture cross length  */
  REAL GetFractureCrossLength() {
    return m_fracture_epsilon;
  }

  REAL FractureFactor(TPZMaterialData &data);

  virtual void SolutionInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec, const std::map<int, TPZMaterialDataT<STATE>> &datarightvec, int var, TPZVec<STATE> &Solout) override {
    DebugStop();
  }

  virtual void SolutionInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec, const std::map<int, TPZMaterialDataT<STATE>> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl *left, TPZCompEl *right) override {
    DebugStop();
  }

  virtual void
  SolutionInterface(const TPZMaterialDataT<STATE> &data, const TPZMaterialDataT<STATE> &dataleft, const TPZMaterialDataT<STATE> &dataright, const int var, TPZVec<STATE> &Solout) override {
    DebugStop();
  }

  virtual void
  ContributeInterface(const TPZMaterialDataT<STATE> &data, const TPZMaterialDataT<STATE> &dataleft, const TPZMaterialDataT<STATE> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
  }

  virtual void
  ContributeBCInterface(const TPZMaterialDataT<STATE> &data, const TPZMaterialDataT<STATE> &dataleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override {
  }
};