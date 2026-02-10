//
//  Created by Giovane Avancini on 05/09/25.
//

#pragma once

#include "Material/TPZMatTypes.h"
#include "pzreal.h"
#include "pzvec.h"
#include <functional>
#include <random>

class TSFFunctionsGenerator {
public:
  /// Enumerate defining the function type
  enum class EP0FunctionType { ENone = 0,
                               EConstantFunction = 1 };

  enum class ES0FunctionType { ENone = 0,
                               EConstantFunction = 1,
                               EPiecewiseFunction = 2,
                               ECircleLevelSetFunction = 4,
                               ERandomCirclesFunction = 5 };

  enum class EDarcyBCFunctionType { ENone = 0,
                                    EHydrostaticPressure = 1,
                                    EInferiorSaturatedPressure = 2,
                                    ESuperiorSaturatedPressure = 3 };

  enum class ETransportBCFunctionType { ENone = 0 };

  TSFFunctionsGenerator() {
    fS0FuncType = ES0FunctionType::ENone;
    fDarcyBCFuncType = EDarcyBCFunctionType::ENone;
    fTransportBCFuncType = ETransportBCFunctionType::ENone;
    fVal = 0.0;
  };

  TSFFunctionsGenerator &operator=(const TSFFunctionsGenerator &other) {
    fS0FuncType = other.fS0FuncType;
    fDarcyBCFuncType = other.fDarcyBCFuncType;
    fTransportBCFuncType = other.fTransportBCFuncType;
    fVal = other.fVal;
    return *this;
  };

  TSFFunctionsGenerator(const TSFFunctionsGenerator &other) {
    fS0FuncType = other.fS0FuncType;
    fDarcyBCFuncType = other.fDarcyBCFuncType;
    fTransportBCFuncType = other.fTransportBCFuncType;
    fVal = other.fVal;
  };

  ~TSFFunctionsGenerator() {}

  /// @brief
  /// @param funcType - function identifier
  /// @param value - optional to set a constant value
  void SetS0FuncType(ES0FunctionType funcType, REAL value = 0.0) {
    fS0FuncType = funcType;
    fVal = value;
  }

  /// @brief
  /// @param funcType - function identifier
  /// @param value - optional to set a constant value
  void SetP0FuncType(EP0FunctionType funcType, REAL value = 0.0) {
    fP0FuncType = funcType;
    fVal = value;
  }

  /// @brief
  /// @param funcType - function identifier
  /// @param value - optional to set a constant value
  void SetDarcyBCFuncType(EDarcyBCFunctionType funcType, REAL value = 0.0) {
    fDarcyBCFuncType = funcType;
    fVal = value;
  }

  /// @brief
  /// @param funcType - function identifier
  /// @param value - optional to set a constant value
  void SetTransportBCFuncType(ETransportBCFunctionType funcType, REAL value = 0.0) {
    fTransportBCFuncType = funcType;
    fVal = value;
  }

  std::function<REAL(const TPZVec<REAL> &)> CreateP0() {
    switch (fP0FuncType) {
    case EP0FunctionType::EConstantFunction: {
      REAL p0 = fVal;
      return [p0](const TPZVec<REAL> &pt) -> REAL {
        REAL p = p0;
        return p;
      };
    } break;
    case EP0FunctionType::ENone: {
      return nullptr;
    } break;
    default: {
      std::cout << " Function not implemented " << std::endl;
      DebugStop();
      return [](const TPZVec<REAL> &pt) -> REAL {
        return 0;
      };
    } break;
    }
  }

  std::function<REAL(const TPZVec<REAL> &)> CreateS0() {
    switch (fS0FuncType) {
    case ES0FunctionType::EConstantFunction: {
      REAL s0 = fVal;
      return [s0](const TPZVec<REAL> &pt) -> REAL {
        REAL s = s0;
        return s;
      };
    } break;
    case ES0FunctionType::ECircleLevelSetFunction: {
      return [](const TPZVec<REAL> &pt) -> REAL {
        REAL x, y, s, r, c, f;
        r = 0.15; // mm
        c = 0;
        x = pt[0];
        y = pt[1];
        f = -r * r + (-202.65 + x) * (-202.65 + x) + (-800 + y) * (-800 + y) - c;
        if (f < 0) {
          s = 0.0;
        } else {
          s = 1.0;
        }
        return s;
      };
    } break;
    case ES0FunctionType::EPiecewiseFunction: {
      REAL h = fVal;
      return [h](const TPZVec<REAL> &pt) -> REAL {
        REAL y, s;
        y = pt[1];
        if (y < h) {
          s = 0.0;
        } else {
          s = 1.0;
        }
        return s;
      };
    } break;
    case ES0FunctionType::ERandomCirclesFunction: {
      // Parameters for the domain and circles
      const REAL inner_radius = 30.15;                // mm
      const REAL outer_radius = inner_radius + 345.0; // mm
      const REAL domain_height = 1000.0;              // mm
      const int min_circles = 50;
      const int max_circles = 100;
      const REAL max_radius = 10.0;

      // Generate random circles only once
      static bool initialized = false;
      static std::vector<std::tuple<REAL, REAL, REAL>> circles; // (center_x, center_y, radius)
      if (!initialized) {
        initialized = true;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> n_dist(min_circles, max_circles);
        std::uniform_real_distribution<REAL> x_dist(inner_radius, outer_radius);
        std::uniform_real_distribution<REAL> y_dist(0.0, domain_height);
        std::uniform_real_distribution<REAL> r_dist(1.0, max_radius);

        int n_circles = n_dist(gen);
        circles.reserve(n_circles);
        for (int i = 0; i < n_circles; ++i) {
          REAL cx = x_dist(gen);
          REAL cy = y_dist(gen);
          REAL r = r_dist(gen);
          circles.emplace_back(cx, cy, r);
        }
      }

      return [](const TPZVec<REAL> &pt) -> REAL {
        REAL x = pt[0];
        REAL y = pt[1];
        for (const auto &circle : circles) {
          REAL cx = std::get<0>(circle);
          REAL cy = std::get<1>(circle);
          REAL r = std::get<2>(circle);
          REAL dx = x - cx;
          REAL dy = y - cy;
          if (dx * dx + dy * dy <= r * r) {
            return 0.0;
          }
        }
        return 1.0;
      };
    } break;
    case ES0FunctionType::ENone: {
      return nullptr;
    } break;
    default: {
      std::cout << " Function not implemented " << std::endl;
      DebugStop();
      return [](const TPZVec<REAL> &pt) -> REAL {
        return 0;
      };
    } break;
    }
  }

  ForcingFunctionBCType<REAL> CreateDarcyBC() {
    switch (fDarcyBCFuncType) {
    case EDarcyBCFunctionType::EHydrostaticPressure: {
      REAL p0 = fVal;
      const REAL rho = 1000.0;    // kg/m3
      const REAL g = 9.81;        // m/s2
      const REAL ref_level = 0.0; // m
      return [p0, rho, g, ref_level](const TPZVec<REAL> &loc, TPZVec<REAL> &rhsVal, TPZFMatrix<REAL> &matVal) {
        REAL z = loc[1]; // Assuming y-coordinate is vertical
        REAL p = p0 + rho * g * (z - ref_level);
        rhsVal[0] = p;
      };
    } break;
    case EDarcyBCFunctionType::EInferiorSaturatedPressure: {
      // This function is a piece-wise function where the pressure decreases linearly from p0 to pf from t0 to tf,
      // and then remains constant at pf.
      REAL p0 = 100.0e3;   // Initial pressure in Pa
      REAL pf = 0.01 * p0; // Final pressure in Pa
      REAL tf = 1.0;       // Time at which pressure reaches pf in seconds
      return [p0, pf, tf](const TPZVec<REAL> &loc, TPZVec<REAL> &rhsVal, TPZFMatrix<REAL> &matVal) {
        REAL t = loc[3]; // Assuming t = loc[3]
        REAL p = (t < tf) ? p0 * (1.0 - 0.99 * t) : pf;
        rhsVal[0] = p;
      };
    } break;
    case EDarcyBCFunctionType::ESuperiorSaturatedPressure: {
      // This function is a piece-wise function where the pressure increases linearly from p0 to pf from t0 to tf,
      // and then remains constant at pf.
      REAL p0 = 100.0;     // Initial pressure in Pa
      REAL pf = 10.0 * p0; // Final pressure in Pa
      REAL tf = 1.0;       // Time at which pressure reaches pf in seconds
      return [p0, pf, tf](const TPZVec<REAL> &loc, TPZVec<REAL> &rhsVal, TPZFMatrix<REAL> &matVal) {
        REAL t = loc[3]; // Assuming t = loc[3]
        REAL p = (t < tf) ? p0 * 10.0 * t : pf;
        rhsVal[0] = p;
      };
    } break;
    case EDarcyBCFunctionType::ENone: {
      return nullptr;
    } break;
    default: {
      std::cout << " Function not implemented " << std::endl;
      DebugStop();
      return [](const TPZVec<REAL> &loc, TPZVec<REAL> &rhsVal, TPZFMatrix<REAL> &matVal) {};
    } break;
    }
  }

  ForcingFunctionBCType<REAL> CreateTransportBC() {
    switch (fTransportBCFuncType) {
    case ETransportBCFunctionType::ENone: {
      return nullptr;
    } break;
    default: {
      std::cout << " Function not implemented " << std::endl;
      DebugStop();
      return [](const TPZVec<REAL> &loc, TPZVec<REAL> &rhsVal, TPZFMatrix<REAL> &matVal) {};
    } break;
    }
  }

private:
  /// @brief function type for initial pressure
  EP0FunctionType fP0FuncType;
  /// @brief function type for initial saturation
  ES0FunctionType fS0FuncType;
  /// @brief function type for Darcy boundary conditions
  EDarcyBCFunctionType fDarcyBCFuncType;
  /// @brief function type for Transport boundary conditions
  ETransportBCFunctionType fTransportBCFuncType;
  /// @brief Constant value used in some functions. Depending on the context, it may represent different physical quantities.
  REAL fVal;
};