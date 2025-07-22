//
//  TMRSTransportMemory.h
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSTransportMemory_h
#define TMRSTransportMemory_h

#include <stdio.h>
#include "TMRSSavable.h"
#include "pzfmatrix.h"

class TMRSTransportMemory : public TMRSSavable {
    
public:
    
    /// Water saturation at last state
    REAL m_sw;
    
    /// Oil saturation at last state
    REAL m_so;
    
    /// Gas saturation at last state
    REAL m_sg;
    
    /// Water saturation at current state
    REAL m_sw_n;
    
    /// Oil saturation at current state
    REAL m_so_n;
    
    /// Gas saturation at current state
    REAL m_sg_n;
    
    /// Default constructor
    TMRSTransportMemory();
    
    /// Copy constructor
    TMRSTransportMemory(const TMRSTransportMemory & other);
    
    /// Assignement constructor
    const TMRSTransportMemory & operator=(const TMRSTransportMemory & other);
    
    /// Desconstructor
    virtual ~TMRSTransportMemory();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TMRSTransportMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    void Setsw(REAL sw)
    {
        m_sw = sw;
    }
    
    REAL sw()
    {
        return m_sw;
    }
    
    void Setso(REAL so)
    {
        m_so = so;
    }
    
    REAL so()
    {
        return m_so;
    }
    
    void Setsg(REAL sg)
    {
        m_sg = sg;
    }
    
    REAL sg()
    {
        return m_sg;
    }
    
    
    void Setsw_n(REAL sw_n)
    {
        m_sw_n = sw_n;
    }
    
    REAL sw_n()
    {
        return m_sw_n;
    }
    
    void Setso_n(REAL so_n)
    {
        m_so_n = so_n;
    }
    
    REAL so_n()
    {
        return m_so_n;
    }
    
    void Setsg_n(REAL sg_n)
    {
        m_sg = sg_n;
    }
    
    REAL sg_n()
    {
        return m_sg_n;
    }
    
};

#endif /* TMRSTransportMemory_h */
