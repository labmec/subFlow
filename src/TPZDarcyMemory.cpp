//
//  TPZDarcyMemory.cpp
//
//  Created by José Villegas on 07/07/20.
//

#include "TPZDarcyMemory.h"

TPZDarcyMemory::TPZDarcyMemory() {
    fTransportCellIndex = -1;
}

TPZDarcyMemory::TPZDarcyMemory(const TPZDarcyMemory & other){
    fTransportCellIndex     =   other.fTransportCellIndex;
}

const TPZDarcyMemory & TPZDarcyMemory::operator=(const TPZDarcyMemory & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    fTransportCellIndex     =   other.fTransportCellIndex;
    
    return *this;
}

TPZDarcyMemory::~TPZDarcyMemory(){

}

const std::string TPZDarcyMemory::Name() const{
    return "TPZDarcyMemory";
}

void TPZDarcyMemory::Write(TPZStream &buf, int withclassid) const{
    buf.Write(fTransportCellIndex);
   
}

void TPZDarcyMemory::Read(TPZStream &buf, void *context){

    buf.Read(&fTransportCellIndex);

}

void TPZDarcyMemory::Print(std::ostream &out) const{
    out << Name();
    out << "\fTransportCellIndex = " << fTransportCellIndex;

}

int TPZDarcyMemory::ClassId() const{
    return Hash("TPZDarcyMemory");
}
