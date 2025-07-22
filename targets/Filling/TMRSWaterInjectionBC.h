#ifndef TMRSWATERINJECTIONBC_H
#define TMRSWATERINJECTIONBC_H

#include "pzreal.h"

class TMRSWaterInjectionBC {
public:
    // Constructor
    TMRSWaterInjectionBC();

    // Destructor
    ~TMRSWaterInjectionBC();

    // Copy Constructor
    TMRSWaterInjectionBC(const TMRSWaterInjectionBC& other);

    // Assignment Operator
    TMRSWaterInjectionBC& operator=(const TMRSWaterInjectionBC& other);

private:
    REAL fMeanPressure;
    REAL fInflow;
    REAL fOutflow;

    
};

#endif // TMRSWATERINJECTIONBC_H