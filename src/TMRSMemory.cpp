//
//  TMRSMemory.cpp
//
//
//  Created by Omar Durán on 10/14/19.
//

#include "TMRSMemory.h"

TMRSMemory::TMRSMemory(){
    
}

TMRSMemory::TMRSMemory(const TMRSMemory & other){
    
}

const TMRSMemory & TMRSMemory::operator=(const TMRSMemory & other){
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }

    return *this;
    
}

TMRSMemory::~TMRSMemory(){
    
}

const std::string TMRSMemory::Name() const{
    return "TMRSMemory";
}

void TMRSMemory::Write(TPZStream &buf, int withclassid) const{
    
}

void TMRSMemory::Read(TPZStream &buf, void *context){
    
}

void TMRSMemory::Print(std::ostream &out) const{
    out << "TMRSMemory " << std::endl;
    TMRSDarcyMemory::Print(out);
    TMRSTransportMemory::Print(out);
}

int TMRSMemory::ClassId() const{
    return Hash("TMRSMemory") ^ TMRSDarcyMemory::ClassId() << 1 ^ TMRSTransportMemory::ClassId() << 2;
}
