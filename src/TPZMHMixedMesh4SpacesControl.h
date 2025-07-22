//
//  TPZMHMixedMeshControl.hpp
//  Numeric Methods
//


#ifndef TPZMHMixedMesh4SpacesControl_hpp
#define TPZMHMixedMesh4SpacesControl_hpp

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
#include "TPZNullMaterial.h"

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedMesh4SpacesControl : public TPZMHMixedMeshControl
{
    
    TPZAutoPointer<TPZCompMesh> fcmeshPressureAverg;
    TPZAutoPointer<TPZCompMesh> fcmeshFluxAverg ;
    
public:
    
    TPZMHMixedMesh4SpacesControl() : TPZMHMixedMeshControl()
    {
    }
    
//    TPZMHMixedMesh4SpacesControl(int dimension):TPZMHMixedMeshControl(dimension){
//        
//    }
//    
//
    TPZVec<TPZAutoPointer<TPZCompMesh> > GetMeshes()
    {
        TPZManVector<TPZAutoPointer<TPZCompMesh>,4> result(4);
        result[0] = fFluxMesh;
        result[1] = fPressureFineMesh;
        result[2] = fcmeshFluxAverg;
        result[3] = fcmeshPressureAverg;
        return result;
    }
    
    TPZMHMixedMesh4SpacesControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices):TPZMHMixedMeshControl( gmesh, coarseindices){
        fcmeshPressureAverg = new TPZCompMesh(gmesh);
        fcmeshFluxAverg = new TPZCompMesh(gmesh);
    }
    TPZMHMixedMesh4SpacesControl(const TPZMHMixedMesh4SpacesControl &copy){
        TPZMHMixedMeshControl::operator=(copy);
        this->operator=(copy);
    }
    
    TPZMHMixedMesh4SpacesControl operator=(const TPZMHMixedMesh4SpacesControl &cp) {
        
        fcmeshPressureAverg = cp.fcmeshPressureAverg;
        fcmeshFluxAverg =cp.fcmeshFluxAverg;
        TPZMHMixedMeshControl::operator=(cp);
        return *this;
    }
//    
    TPZMHMixedMesh4SpacesControl(TPZAutoPointer<TPZGeoMesh> gmesh):TPZMHMixedMeshControl(gmesh)
    {
        
    }
//
//    
//    TPZMHMixedMesh4SpacesControl(const TPZMHMixedMesh4SpacesControl &copy) : TPZMHMixedMeshControl(copy)
//    {
//        
//        fFluxMesh = copy.fFluxMesh;
//    }
//    
//    TPZMHMixedMesh4SpacesControl &operator=(const TPZMHMixedMesh4SpacesControl &cp)
//    {
//        fFluxMesh = cp.fFluxMesh;
//        TPZMHMixedMeshControl::operator=(cp);
//        return *this;
//    }
//    
    void BuildComputationalMesh(bool usersubstructure);
    void CreateHDivPressureMHMMesh();
    void CreateAverageFlux();
    void CreateAveragePressure();
    void BuildMultiPhysicsMesh();
    void HideTheElements();
    
    
    void PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, bool KeepOneLagrangian);
    void GroupandCondenseElements();
    void GroupandCondenseElementsEigen();
//    void HideTheElements();
//    
//    int64_t WhichSubdomain(TPZCompEl *cel);
    
};

#endif /* TPZMHMixedMeshChannelControl_hpp */

