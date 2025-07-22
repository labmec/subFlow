//  TPZMFSolutionTransfer.hpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#ifndef TPZMFSolutionTransfer_hpp
#define TPZMFSolutionTransfer_hpp

#include <stdio.h>
#include "TPZMultiphysicsCompMesh.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

/// auxiliary class that permits to transfer information from a substructured multi physics mesh to and from the atomic meshes
class TPZMFSolutionTransfer
{
   
    /**
     * @brief Structure that relates the origin blocks indexes in the substructured multiphysics mesh and target indexes in the atomic
     * meshes. Is used to transfer the solution.
     */
    struct Match{
        
        /**
         * @brief Block index in the substructured multiphysics mesh.
         */
        int64_t fblocknumber;
        /**
         * @brief Block pointer and block index in the target atomic mesh
         */
        std::pair<TPZBlock*, int64_t> fblockTarget;// block target
        
        /**
         * @brief Empty Constructor
         */
        Match(){
            
        }
        /**
         * @brief Copy Constructor
         */
        Match(const Match & copy){
            fblocknumber=copy.fblocknumber;
            fblockTarget = copy.fblockTarget;
        }
        
        /**
         * @brief operator =
         */
        
        Match &operator=(const Match &other){
            fblocknumber=other.fblocknumber;
            fblockTarget = other.fblockTarget;
            return *this;
        }
        
        /**
         * @brief Default Destructor
         */
        ~Match(){
        }
        
        bool operator<(const Match &other) const {
            return fblocknumber < other.fblocknumber;
        }
        
        /**
         * @brief Transfer the solution from the multiphysic mesh to the atomic meshes.
         * @param mfmesh target mesh.
         */
        
        void TransferFromMultiphysics(TPZCompMesh * mfmesh);
        /**
         * @brief Transfer the solution from the atomic meshes to the multiphysic mesh.
         * @param cmesh is the multiphysics mesh.
         * @note Internaly are taken the corresponding blocks.
         */
        void TransferToMultiphysics(TPZCompMesh * cmesh);
        
    };
    
    /**
     * @brief Structure that stores the match objects relating the meshes in the substructured multi physics mesh and the atomic meshes.
     */
    struct MeshTransferData{
        /**
         * @brief Original mesh (either multiphysics of TPZSubCompMesh).
         */
        TPZCompMesh * fcmesh_ori;
        /**
         * @brief Match objects relating the blocks of fcmesh_ori to the atomic meshes
         */
        std::set<Match> fconnecttransfer;
        
        /**
         * @brief Empty constructor
         */
        MeshTransferData(){
            
        }
        /**
         * @brief Copy constructor
         */
        MeshTransferData(const MeshTransferData & copy){
            fcmesh_ori=copy.fcmesh_ori;
            fconnecttransfer=copy.fconnecttransfer;
           
        }
        /**
         * @brief operator =
         */
        MeshTransferData &operator=(const MeshTransferData &other){
            fcmesh_ori=other.fcmesh_ori;
            fconnecttransfer=other.fconnecttransfer;
            return *this;
        }
        /**
         * @brief Default destructor
         */
        ~MeshTransferData(){
            
        }
        /**
         * @brief Build all teh matches between the multiphysic and the atomic meshes.
           @param cmesh Origin mesh
         */
        void BuildTransferData(TPZCompMesh* cmesh);
        
        /// Insert the connect matches of an arbirtary computational element
        // this method will cast the element to either a condensed element, element group or multiphysics element
        void InsertMatches(TPZCompEl *cel);
        
        /// Insert the connect matches from the multiphysics element
        void InsertMatches(TPZMultiphysicsElement *cel);
        
        /// Insert the connect matches for an element group
        void InsertMatches(TPZElementGroup *celgr);
        
        /// Insert the connect matches for a condensed computational element
        void InsertMatches(TPZCondensedCompEl *cond);
        
        /**
         * @brief Transfer the solution from the mesh (either multiphysics or subcmesh) to the atomic meshes for every match stored in fconnecttransfer
         */
        void TransferFromMultiphysics();
        
        /**
         * @brief Transfer the solution from the atomic meshes to the mesh (either multiphysics or subcmesh) for every match stored in fconnecttransfer
         */
        void TransferToMultiphysics();
       
        
    };
    public:
        /**
         * @brief Objects vector MeshTransferData, the transference should be done for the multiphysics mesh and their substructures.
         */
        TPZStack<MeshTransferData> fmeshTransfers;
        /**
         * @brief Build all the MeshTransferData for the multiphysic and their substructure.
         @param mfcmesh Origin mesh.
         */
        void BuildTransferData(TPZCompMesh* mfcmesh);
        /**
         * @brief Transfer the solution from the multiphysics mesh to the atomica meshes for every MeshTransferData stored in fmeshTransfers.
         */
        void TransferFromMultiphysics();
        /**
         * @brief Transfer teh solution from the atomic meshes to the multiphysic mesh for every MeshTransferData stored in fmeshTransfers.
         */
        void TransferToMultiphysics();
    
        
};

#endif /* TPZMFSolutionTransfer_hpp */
