#include "TSFFilterCakeMemory.h"

TSFFilterCakeMemory::TSFFilterCakeMemory() : fAccumulatedVolume(0.) {
  
}

TSFFilterCakeMemory::TSFFilterCakeMemory(const TSFFilterCakeMemory &other) : TSFSavable(other), fAccumulatedVolume(other.fAccumulatedVolume) {}

const TSFFilterCakeMemory &TSFFilterCakeMemory::operator=(const TSFFilterCakeMemory &other) {
  if (this != &other) {
    TSFSavable::operator=(other);
    fAccumulatedVolume = other.fAccumulatedVolume;
  }
  return *this;
}

int TSFFilterCakeMemory::ClassId() const { return Hash("TSFFilterCakeMemory"); }

TSFFilterCakeMemory::~TSFFilterCakeMemory() {}

const std::string TSFFilterCakeMemory::Name() const { return "TSFFilterCakeMemory"; }

void TSFFilterCakeMemory::Write(TPZStream &buf, int withclassid) const {
  buf.Write(&fAccumulatedVolume);
}

void TSFFilterCakeMemory::Read(TPZStream &buf, void *context) {
  buf.Read(&fAccumulatedVolume);
}

void TSFFilterCakeMemory::Print(std::ostream &out) const {
  out << "TSFFilterCakeMemory: " << std::endl;
  out << "Accumulated Volume: " << fAccumulatedVolume << std::endl;
}