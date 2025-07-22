//
//  TPZMHMixedMeshControl.hpp
//  PZ
//


#ifndef TPZMHMixedMeshWithTransportControl_hpp
#define TPZMHMixedMeshWithTransportControl_hpp



#include <stdio.h>

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "TPZVTKGeoMesh.h"
//#include "TPZNullMaterial.h"
#include "TMRSApproxSpaceGenerator.h"

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedMeshWithTransportControl : public TPZMHMixedMeshControl
{
    
    TPZAutoPointer<TPZCompMesh> fcmeshTransport;
    TMRSApproxSpaceGenerator *faproxspace;
    
public:
    
    TPZMHMixedMeshWithTransportControl() : TPZMHMixedMeshControl()
    {
      
    }
    void SetApproxSpaceGenerator(TMRSApproxSpaceGenerator *aproxspace){
        faproxspace =aproxspace;
    }
    TMRSApproxSpaceGenerator *GetApproxSpaceGenerator(){
        return faproxspace;
    }
 
    TPZMHMixedMeshWithTransportControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices):TPZMHMixedMeshControl( gmesh, coarseindices){
        fcmeshTransport = new TPZCompMesh(gmesh);
    }
    TPZMHMixedMeshWithTransportControl(const TPZMHMixedMeshWithTransportControl &copy){
        TPZMHMixedMeshControl::operator=(copy);
        this->operator=(copy);
    }
    
    TPZMHMixedMeshWithTransportControl operator=(const TPZMHMixedMeshWithTransportControl &cp) {
        
        fcmeshTransport = cp.fcmeshTransport;
        TPZMHMixedMeshControl::operator=(cp);
        return *this;
    }
    
    TPZMHMixedMeshWithTransportControl(TPZAutoPointer<TPZGeoMesh> gmesh):TPZMHMixedMeshControl(gmesh)
    {
        
    }
   
    void BuildComputationalMesh(bool usersubstructure);
    void CreateHDivPressureMHMMesh();
    void CreateTransport();
    void BuildMultiPhysicsMesh();

};

#endif /* TPZMHMixedMeshChannelControl_hpp */

