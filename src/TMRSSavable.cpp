//
//  TMRSSavable.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#include "TMRSSavable.h"
#include "Hash/TPZHash.h"

TMRSSavable::TMRSSavable() {
}

TMRSSavable::TMRSSavable(const TMRSSavable& orig) {
}

TMRSSavable::~TMRSSavable() {
}

int TMRSSavable::ClassId() const {
    return Hash("TMRSSavable");
}
