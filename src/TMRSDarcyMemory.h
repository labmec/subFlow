//
//  TMRSDarcyMemory.hpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSDarcyMemory_h
#define TMRSDarcyMemory_h

#include <stdio.h>
#include "TMRSSavable.h"
#include "pzfmatrix.h"

class TMRSDarcyMemory : public TMRSSavable {
    
public:
    
    /// Absolute permeability
    TPZFNMatrix<9,REAL> m_kappa;
    
    /// Inver of absolute permeability
    TPZFNMatrix<9,REAL> m_kappa_inv;
    
    /// normal permeability
    REAL m_kappa_normal;
    
    /// Frac indexes
    std::pair<int, int> m_fracindexes;
    
    /// lagrangian porosity
    REAL m_phi;
    
    /// weight pressure at last state
    REAL m_p;
    
    /// weight pressure at current state
    REAL m_p_n;
    
    /// mass flux
    TPZManVector<REAL,3> m_flux;
    
    /// Fracture opening
    REAL m_ad;
    
    /// Default constructor
    TMRSDarcyMemory();
    
    /// Copy constructor
    TMRSDarcyMemory(const TMRSDarcyMemory & other);
    
    /// Assignement constructor
    const TMRSDarcyMemory & operator=(const TMRSDarcyMemory & other);
    
    /// Desconstructor
    virtual ~TMRSDarcyMemory();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TMRSDarcyMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;

    
    /// Set a permeability tensor
    void SetPermeability(TPZFNMatrix<9,REAL> kappa){
        m_kappa =  kappa;
        m_kappa_inv.Zero();
        DebugStop(); // Compute inverse please
    }
    
    /// Set an isotropic permeability tensor
    void SetPermeability(REAL kappa){
        m_kappa.Zero();
        m_kappa_inv.Zero();
        for (int i = 0; i < 3; i++) {
            m_kappa(i,i) = kappa;
            m_kappa_inv(i,i) = 1.0/kappa;
        }
    }
    
    /// Get an permeability tensor
    void SetPermeability(TPZFNMatrix<9,REAL> & kappa){
        
        if (kappa.Rows() != 3 && kappa.Cols() != 3) {
            DebugStop();
        }
        
        m_kappa = kappa;
        m_kappa.Inverse(m_kappa_inv, ELU);
    }
    
    /// Dimensional factor associated to multi-dimensional mixed operators
    void SetDimensionalFactor(REAL d){
        m_ad = d;
    }
    
    /// Set pore pressure at last state
    void Setp(REAL p)
    {
        m_p = p;
    }
    
    /// Get pore pressure at last state
    REAL p()
    {
        return m_p;
    }
    
    /// Set pore pressure at current state
    void Setp_n(REAL p_n)
    {
        m_p_n = p_n;
    }
    
    /// Get pore pressure at current state
    REAL p_n()
    {
        return m_p_n;
    }

    
    /// Set lagrangian porosity
    void Setphi(REAL phi)
    {
        m_phi = phi;
    }
    
    /// Get lagrangian porosity
    REAL phi()
    {
        return m_phi;
    }
    
    /// Set the flux
    void SetFlux(TPZVec<REAL> &flux)
    {
        m_flux = flux;
    }
    
    /// Set the flux
    void Flux(TPZVec<REAL> &flux)
    {
        flux = m_flux;
    }
    
    /// Set the normal permeability
    void SetNormalPermeability(REAL &normal_kappa)
    {
        m_kappa_normal = normal_kappa;
    }
    
    /// Get the normal permeability
    REAL GetNormalPermeability()
    {
       return m_kappa_normal;
    }
    
    /// Set the normal permeability
    void SetFracindexes(std::pair<int,int> &fracindexes)
    {
        m_fracindexes = fracindexes;
    }
    
    /// Get the normal permeability
    std::pair<int,int> GetFracindexes()
    {
       return m_fracindexes;
    }
    
};

#endif /* TMRSDarcyMemory_h */
