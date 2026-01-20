//
//  Created by Giovane Avancini on 10/01/26.
//  Copied from iMRS and is evolving to fit subFlow needs
//  PLEASE ADAPT THIS FILE
//

#include "TSFAlgebraicTransport.h"

/// Default constructor
TSFAlgebraicTransport::TSFAlgebraicTransport() {}

/// Copy constructor
TSFAlgebraicTransport::TSFAlgebraicTransport(const TSFAlgebraicTransport &other) {
  fNFluxCoefficients = other.fNFluxCoefficients;
  fNVolumesTransport = other.fNVolumesTransport;
  fCellsData = other.fCellsData;
  fInterfaceData = other.fInterfaceData;
  fNFluxCoefficients = other.fNFluxCoefficients;
  fboundaryCMatVal = other.fboundaryCMatVal;
}

/// Assignment operator
const TSFAlgebraicTransport &TSFAlgebraicTransport::operator=(const TSFAlgebraicTransport &other) {
  fNFluxCoefficients = other.fNFluxCoefficients;
  fNVolumesTransport = other.fNVolumesTransport;
  fCellsData = other.fCellsData;
  fInterfaceData = other.fInterfaceData;
  fboundaryCMatVal = other.fboundaryCMatVal;
  return *this;
}

/// Default destructor
TSFAlgebraicTransport::~TSFAlgebraicTransport() {}