//
//  TMRSMemory.h
//
//
//  Created by Omar Durán on 10/14/19.
//

#ifndef TMRSMemory_h
#define TMRSMemory_h

#include <stdio.h>
#include "TMRSDarcyMemory.h"
#include "TMRSTransportMemory.h"

class TMRSMemory : public TMRSDarcyMemory, public TMRSTransportMemory {
    
public:
    
    /// Default constructor
    TMRSMemory();
    
    /// Copy constructor
    TMRSMemory(const TMRSMemory & other);
    
    /// Assignement constructor
    const TMRSMemory & operator=(const TMRSMemory & other);
    
    /// Desconstructor
    virtual ~TMRSMemory();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TMRSMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
};

#endif /* TMRSMemory_h */
