//
//  TMRSTransportMemory.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#include "TMRSTransportMemory.h"


TMRSTransportMemory::TMRSTransportMemory(){
    m_sw    = 0.0;
    m_so    = 0.0;
    m_sg    = 0.0;
    m_sw_n  = 0.0;
    m_so_n  = 0.0;
    m_sg_n  = 0.0;
}

TMRSTransportMemory::TMRSTransportMemory(const TMRSTransportMemory & other){
    m_sw    = other.m_sw;
    m_so    = other.m_so;
    m_sg    = other.m_sg;
    m_sw_n  = other.m_sw_n;
    m_so_n  = other.m_so_n;
    m_sg_n  = other.m_sg_n;
}

const TMRSTransportMemory & TMRSTransportMemory::operator=(const TMRSTransportMemory & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    m_sw    = other.m_sw;
    m_so    = other.m_so;
    m_sg    = other.m_sg;
    m_sw_n  = other.m_sw_n;
    m_so_n  = other.m_so_n;
    m_sg_n  = other.m_sg_n;
    return *this;
}

TMRSTransportMemory::~TMRSTransportMemory(){
    
}

const std::string TMRSTransportMemory::Name() const{
    return "TMRSTransportMemory";
}

void TMRSTransportMemory::Write(TPZStream &buf, int withclassid) const{
    buf.Write(&m_sw);
    buf.Write(&m_so);
    buf.Write(&m_sg);
    buf.Write(&m_sw_n);
    buf.Write(&m_so_n);
    buf.Write(&m_sg_n);
}

void TMRSTransportMemory::Read(TPZStream &buf, void *context){
    buf.Read(&m_sw);
    buf.Read(&m_so);
    buf.Read(&m_sg);
    buf.Read(&m_sw_n);
    buf.Read(&m_so_n);
    buf.Read(&m_sg_n);
}

void TMRSTransportMemory::Print(std::ostream &out) const{
    out << "\n m_sw     = " << m_sw;
    out << "\n m_so     = " << m_so;
    out << "\n m_sg     = " << m_sg;
    out << "\n m_sw_n   = " << m_sw_n;
    out << "\n m_so_n   = " << m_so_n;
    out << "\n m_sg_n   = " << m_sg_n;
}

int TMRSTransportMemory::ClassId() const{
    return Hash("TMRSTransportMemory");
}
