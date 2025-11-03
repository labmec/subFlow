//
//  Created by Giovane Avancini on 02/09/25.
//

#include "TSFSavable.h"

TSFSavable::TSFSavable() {}

TSFSavable::TSFSavable(const TSFSavable &orig) {}

TSFSavable::~TSFSavable() {}

int TSFSavable::ClassId() const {
  return Hash("TSFSavable");
}
