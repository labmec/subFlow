#include "TSFMatFilterCake.h"
#include "TSFFilterCakeMemory.h"

template <class TMEM>
TSFMatFilterCake<TMEM>::TSFMatFilterCake() : TPZRegisterClassId(&TSFMatFilterCake<TMEM>::ClassId), TBase(), fDimension(0), fSimData(nullptr), fAlpha(0.0), fCoef(0.0), fIsAxisymmetric(false) {}

template <class TMEM>
TSFMatFilterCake<TMEM>::TSFMatFilterCake(int matid, int dimension) : TPZRegisterClassId(&TSFMatFilterCake<TMEM>::ClassId), TBase(matid), fDimension(dimension), fSimData(nullptr), fAlpha(0.0), fCoef(0.0), fIsAxisymmetric(false) {}

template <class TMEM>
TSFMatFilterCake<TMEM>::TSFMatFilterCake(int matid, int dimension, TSFProblemData *simData) : TPZRegisterClassId(&TSFMatFilterCake<TMEM>::ClassId), TBase(matid), fDimension(dimension), fSimData(simData), fAlpha(0.0), fCoef(0.0), fIsAxisymmetric(false) {
  if (!fSimData) DebugStop();
  fIsAxisymmetric = fSimData->fTNumerics.fIsAxisymmetric;
  REAL C = fSimData->fTFilterCakeProperties.fC[matid];
  REAL Porosity = fSimData->fTFilterCakeProperties.fPorosity[matid];
  REAL Density = fSimData->fTFilterCakeProperties.fDensity[matid];
  REAL ParticleDiameter = fSimData->fTFilterCakeProperties.fParticleDiameter[matid];
  fAlpha = C / ((1.0 - Porosity) * Density);
  fCoef = (ParticleDiameter * ParticleDiameter) * Porosity * Porosity * Porosity / (180.0 * (1.0 - Porosity) * (1.0 - Porosity));
}

template <class TMEM>
TSFMatFilterCake<TMEM>::TSFMatFilterCake(const TSFMatFilterCake &other) : TPZRegisterClassId(&TSFMatFilterCake<TMEM>::ClassId), TBase(other) {
  fDimension = other.fDimension;
  fSimData = other.fSimData;
  fAlpha = other.fAlpha;
  fCoef = other.fCoef;
  fIsAxisymmetric = other.fIsAxisymmetric;
}

template <class TMEM>
TSFMatFilterCake<TMEM> &TSFMatFilterCake<TMEM>::operator=(const TSFMatFilterCake &other) {
  if (this == &other) return *this;
  fDimension = other.fDimension;
  fSimData = other.fSimData;
  fAlpha = other.fAlpha;
  fCoef = other.fCoef;
  fIsAxisymmetric = other.fIsAxisymmetric;
  return *this;
}

template <class TMEM>
TSFMatFilterCake<TMEM>::~TSFMatFilterCake() = default;

template <class TMEM>
[[nodiscard]] int TSFMatFilterCake<TMEM>::ClassId() const {
  return Hash("TSFMatFilterCake") ^ TBase::ClassId();
}

template <class TMEM>
void TSFMatFilterCake<TMEM>::Print(std::ostream &out) const {
  out << "Material Name: " << Name() << std::endl;
  out << "Material ID: " << TBase::Id() << std::endl;
  out << "Dimension: " << fDimension << std::endl;
  out << "Alpha (C / (1-phi)rho): " << fAlpha << std::endl;
  out << "Coefficient (Carman-Kozeny): " << fCoef << std::endl;
  out << "Axisymmetric: " << (fIsAxisymmetric ? "Yes" : "No") << std::endl;
}

template <class TMEM>
void TSFMatFilterCake<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
  TMEM &memory = TPZMatWithMem<TMEM>::GetMemory()->operator[](0); // the same memory item is used for all integration points
  REAL accumulatedVolume = memory.GetAccumulatedVolume();
  REAL thickness = accumulatedVolume / fAlpha;
  REAL absperm = fCoef / thickness;
  REAL invperm = fSimData->fTFluidProperties.fWaterViscosity * accumulatedVolume / (absperm * fAlpha);

  TPZFNMatrix<20, REAL> &phiU = datavec[0].phi;
  int nPhiU = phiU.Rows();

  TPZManVector<STATE, 3> usol = datavec[0].sol[0];

  REAL axiFactor = 1.0;
  if (fIsAxisymmetric) // Axisymmetric: assuming radius is aligned with the x axis
  {
    REAL r = datavec[0].x[0];
    axiFactor = 1.0 / (2.0 * M_PI * r);
    usol[0] *= axiFactor;
  }

  for (int i = 0; i < nPhiU; i++) {
    ef(i) += -1.0 * invperm * phiU(i, 0) * usol[0] * weight;
    for (int j = 0; j < nPhiU; j++) {
      ek(i, j) += invperm * phiU(i, 0) * phiU(j, 0) * weight * axiFactor;
    }
  }
}

template <class TMEM>
void TSFMatFilterCake<TMEM>::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
  int ndata = datavec.size();
  for (int i = 0; i < ndata; i++) {
    datavec[i].SetAllRequirements(false);
  }
  datavec[0].fNeedsSol = true;
}

template <class TMEM>
int TSFMatFilterCake<TMEM>::VariableIndex(const std::string &name) const {
  if (!strcmp("Flux", name.c_str())) return 1;
  return TPZMaterial::VariableIndex(name);
}

template <class TMEM>
int TSFMatFilterCake<TMEM>::NSolutionVariables(int var) const {
  if (var == 1) return 1;
  return TPZMaterial::NSolutionVariables(var);
}

template <class TMEM>
void TSFMatFilterCake<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
  // No contribution for the filter cake material since it is an internal interface
}

template <class TMEM>
void TSFMatFilterCake<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {}

template <class TMEM>
void TSFMatFilterCake<TMEM>::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {}

template class TSFMatFilterCake<TSFFilterCakeMemory>;