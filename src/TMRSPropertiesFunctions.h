//
//  TMRSPropertiesFunctions.h
//  Buckley_Levertt
//
//  Created by Omar Durán on 6/16/20.
//

#ifndef TMRSPropertiesFunctions_h
#define TMRSPropertiesFunctions_h

#include "pzreal.h"
#include "TRMSpatialPropertiesMap.h"
#include <random>

class TMRSPropertiesFunctions
{
    public:
    
    /// Enumerate defining the function type
    enum EFunctionType { EConstantFunction = 0, EPiecewiseFunction = 1, ECircleLevelSetFunction = 2, EUNISIMFunction = 3, ESPECase10Function = 4, ERandomCirclesFunction = 5 };


    TMRSPropertiesFunctions(){
        m_function_type_kappa = EConstantFunction;
        m_function_type_phi = EConstantFunction;
        m_function_type_s0 = EConstantFunction;
        m_constant_val = 0.0;
    }
    TMRSPropertiesFunctions &operator=(TMRSPropertiesFunctions &other)
    {
        m_function_type_kappa = other.m_function_type_kappa;
        m_function_type_phi = other.m_function_type_phi;
        m_function_type_s0 = other.m_function_type_s0;
        return *this;
    };
    TMRSPropertiesFunctions(TMRSPropertiesFunctions &other)
    {
        m_function_type_kappa = other.m_function_type_kappa;
        m_function_type_phi = other.m_function_type_phi;
        m_function_type_s0 = other.m_function_type_s0;
    };
    ~TMRSPropertiesFunctions(){
        
    }
    
    void set_function_type_kappa(EFunctionType function_type_kappa){
        m_function_type_kappa   = function_type_kappa;
    }
    
    void set_function_type_phi(EFunctionType function_type_phi){
        m_function_type_phi     = function_type_phi;
    }
    
    /// @brief 
    /// @param function_type_s0 - function identifier
    /// @param value - optional to set a constant value if EConstantFunction is selected
    void set_function_type_s0(EFunctionType function_type_s0, REAL value = 0.0){
        m_function_type_s0      = function_type_s0;
        m_constant_val = value;
    }
    
    std::function<std::vector<REAL>(const TPZVec<REAL> & )> Create_Kappa_Phi(TRMSpatialPropertiesMap & map){
        
        return [& map] (const TPZVec<REAL> & pt) -> std::vector<REAL> {
            std::vector<REAL> kappa_and_phi;
            TPZManVector<REAL,3> x(pt);
            map.SampleKappaAndPhi(x,kappa_and_phi);
            return kappa_and_phi;
        };
        
    }
    
    std::function<std::vector<REAL>(const TPZVec<REAL> & )> Create_Kappa_Phi(){
        
        return [] (const TPZVec<REAL> & pt) -> std::vector<REAL> {
            std::vector<REAL> kappa_and_phi;
            REAL k_c,phi_c;
            k_c = 1.0e-7;
            phi_c = 0.1;
            REAL x = pt[0];
            REAL y = pt[1];
            REAL z = pt[2];
            REAL kx = k_c ;//* fabs(std::cos(0.2*x)*std::sin(0.1*y)*std::sin(0.1*z)) + k_c;
            REAL ky = 10.0 * k_c ;//* fabs(std::sin(0.1*x)*std::cos(0.2*y)*std::sin(0.1*z)) + k_c;
            REAL kz = 10.0 * k_c ;//* fabs(std::sin(0.1*x)*std::sin(0.1*y)*std::cos(0.2*z)) + k_c;
            REAL phi = phi_c ;//* fabs(std::cos(0.2*x)*std::cos(0.3*y)*std::cos(0.1*z)) + phi_c;
            kappa_and_phi.push_back(kx);
            kappa_and_phi.push_back(ky);
            kappa_and_phi.push_back(kz);
            kappa_and_phi.push_back(phi);
            return kappa_and_phi;
        };
        
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Kx(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kx;
                            x = pt[0];
                            y = pt[1];
                            kx = 1.0e-7;
                            return kx;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kx,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                kx = 1.0e-10;
                            }else{
                                kx = 1.0e-7;
                            }
                            return kx;
                        };
                }
                break;
            case EPiecewiseFunction:
            {
                return [] (const TPZVec<REAL> & pt) -> REAL {
                        REAL y,kx;
                        y = pt[1];
                        if(y>5){
                             kx = 1.0e-10;
                         }else{
                             kx = 1.0e-7;
                         }
                         return kx;
                    };
            }
            break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Ky(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,ky;
                            x = pt[0];
                            y = pt[0];
                            ky = 1.0e-7;
                            return ky ;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,ky,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                 ky = 1.0e-10;
                             }else{
                                 ky = 1.0e-7;
                             }
                             return ky;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,ky;
                            y = pt[1];
                            if(y>5){
                                 ky = 1.0e-10;
                             }else{
                                 ky = 1.0e-7;
                             }
                             return ky;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_Kz(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kz;
                            x = pt[0];
                            y = pt[1];
                            kz = 1.0e-6;
                            return kz ;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kz,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                 kz = 1.0e-10;
                             }else{
                                 kz = 1.0e-7;
                             }
                             return kz;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,kz;
                            y = pt[1];
                            if(y>5){
                                 kz = 1.0e-10;
                             }else{
                                 kz = 1.0e-7;
                             }
                             return kz;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_phi(){
        
        switch (m_function_type_phi) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,phi;
                            x = pt[0];
                            y = pt[1];
                            phi = 0.1;
                            return phi;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,phi,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                phi = 0.01;
                            }else{
                                phi = 0.1;
                            }
                            return phi;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,phi;
                            y = pt[1];
                            if(y<5){
                                 phi = 0.01;
                             }else{
                                 phi = 0.1;
                             }
                             return phi;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_s0(){
        
        switch (m_function_type_s0) {
            case EConstantFunction:
                {
                    REAL s0 = m_constant_val;
                    return [s0] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s;
                            x = pt[0];
                            y = pt[1];
                            s = s0;
                            return s;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s,r,c,f;
                            r = 0.15; // mm
                            c = 0;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-202.65 + x)*(-202.65 + x) + (-800 + y)*(-800 + y) - c;
                            if(f<0){
                                s = 0.0;
                            }else{
                                s = 1.0;
                            }
                            return s;
                        };
                }
                break;
            case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,s;
                            y = pt[1];
                            if(y<0.3){
                                s = 1.0;
                            }else if (y>=0.3 && y<0.7){
                                s = 0.0;
                            }else{
                                s = 1.0;
                             }
                             return s;
                        };
                }
                break;
            case ERandomCirclesFunction:
                {
                    // Parameters for the domain and circles
                    const REAL inner_radius = 30.15; // mm
                    const REAL outer_radius = inner_radius + 345.0; // mm
                    const REAL domain_height = 1000.0; // mm
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

                    return [] (const TPZVec<REAL> & pt) -> REAL {
                        REAL x = pt[0];
                        REAL y = pt[1];
                        for (const auto& circle : circles) {
                            REAL cx = std::get<0>(circle);
                            REAL cy = std::get<1>(circle);
                            REAL r = std::get<2>(circle);
                            REAL dx = x - cx;
                            REAL dy = y - cy;
                            if (dx*dx + dy*dy <= r*r) {
                                return 0.0;
                            }
                        }
                        return 1.0;
                    };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    private:
    
    EFunctionType m_function_type_kappa;
    EFunctionType m_function_type_phi;
    EFunctionType m_function_type_s0;
    REAL m_constant_val;
};

#endif /* TMRSPropertiesFunctions_h */
