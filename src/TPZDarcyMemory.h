//
//  TPZDarcyMemory.hpp
//
//  Created by José Villegas on 07/07/20.
//

#ifndef TPZDarcyMemory_h
#define TPZDarcyMemory_h

#include <stdio.h>
#include "TMRSSavable.h"
#include "pzfmatrix.h"

class TPZDarcyMemory : public TMRSSavable {
    
public:
    
    int fTransportCellIndex;
    
    /// Default constructor
    TPZDarcyMemory();
    
    /// Copy constructor
    TPZDarcyMemory(const TPZDarcyMemory & other);
    
    /// Assignement constructor
    const TPZDarcyMemory & operator=(const TPZDarcyMemory & other);
    
    /// Desconstructor
    virtual ~TPZDarcyMemory();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZDarcyMemory & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;


    
    
};

#endif /* TPZDarcyMemory_h */
