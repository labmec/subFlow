//
//  TMRSSavable.h
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#ifndef TMRSSavable_h
#define TMRSSavable_h

#include <stdio.h>
#include "TPZSavable.h"


class TMRSSavable : public virtual TPZSavable {
public:
    TMRSSavable();
    TMRSSavable(const TMRSSavable& orig);
    virtual int ClassId() const;
    virtual ~TMRSSavable();
private:
    
};


#endif /* TMRSSavable_h */
