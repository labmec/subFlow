//
//  Created by Giovane Avancini on 02/09/25.
//

#pragma once

#include "Hash/TPZHash.h"
#include "TPZSavable.h"
#include <stdio.h>

class TSFSavable : public virtual TPZSavable {
public:
  TSFSavable();
  TSFSavable(const TSFSavable &orig);
  virtual int ClassId() const;
  virtual ~TSFSavable();
};