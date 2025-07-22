//
//  TPZMFSolutionTransfer.cpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#include "TPZMFSolutionTransfer.h"


void TPZMFSolutionTransfer::Match::TransferFromMultiphysics(TPZCompMesh * cmesh){
    
    TPZBlock &blockMF = cmesh->Block();
    TPZBlock* blockAtomic = fblockTarget.first;
    TPZFMatrix<STATE> &matMF = *blockMF.Matrix<STATE>();
    TPZFMatrix<STATE> &matAT = *blockAtomic->Matrix<STATE>();
    int64_t seqMF = fblocknumber;
    int64_t seqAto = fblockTarget.second;
    int blocksizeAto = blockAtomic->Size(seqAto);
    int blocksizeMF = blockMF.Size(seqMF);
    if(blocksizeAto!=blocksizeMF){
        DebugStop();
    }
    for (int idf=0; idf<blocksizeAto; idf++) {
        int indAT = blockAtomic->Index(seqAto, idf);
        int indMF = blockMF.Index(seqMF, idf);
        matAT(indAT,0) = matMF(indMF,0);
//        blockAtomic->Put(seqAto, idf, 0, blockMF.Get(seqMF, idf, 0));
    }
}
void TPZMFSolutionTransfer::Match::TransferToMultiphysics(TPZCompMesh * cmesh){
//    cmesh->InitializeBlock();
    TPZBlock &blockMF = cmesh->Block();
    TPZBlock* blockAtomic = fblockTarget.first;
    TPZFMatrix<STATE> &matMF = *blockMF.Matrix<STATE>();
    TPZFMatrix<STATE> &matAT = *blockAtomic->Matrix<STATE>();
    int seqMF = fblocknumber;
    int seqAto = fblockTarget.second;
    
    int blocksizeAto = blockAtomic->Size(seqAto);
    int blocksizeMF = blockMF.Size(seqMF);
    if(blocksizeAto!=blocksizeMF){
        DebugStop();
    }
    for (int idf=0; idf<blocksizeAto; idf++) {
        int indAT = blockAtomic->Index(seqAto, idf);
        int indMF = blockMF.Index(seqMF, idf);
        matMF(indMF,0) = matAT(indAT,0);

//        blockMF.Put(seqMF, idf, 0, blockAtomic->Get(seqAto, idf, 0));
    }
}

void TPZMFSolutionTransfer::MeshTransferData::TransferToMultiphysics(){
   
    fcmesh_ori->InitializeBlock();
    for (auto match : fconnecttransfer) {
        match.TransferToMultiphysics(fcmesh_ori);
    }
}
void TPZMFSolutionTransfer::MeshTransferData::TransferFromMultiphysics(){
    for (auto match : fconnecttransfer) {
        match.TransferFromMultiphysics(fcmesh_ori);
    }
}

/// Insert the connect matches of an arbirtary computational element
// this method will cast the element to either a condensed element, element group or multiphysics element
void TPZMFSolutionTransfer::MeshTransferData::InsertMatches(TPZCompEl *cel)
{
    TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if(mfcel) return InsertMatches(mfcel);
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    if(elgr) return InsertMatches(elgr);
    TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
    if(cond) return InsertMatches(cond);
    TPZMultiphysicsInterfaceElement *intf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
    if(intf) return;
    return;
    DebugStop();
}

/// Insert the connect matches from the multiphysics element
void TPZMFSolutionTransfer::MeshTransferData::InsertMatches(TPZMultiphysicsElement *mul)
{
    TPZVec<int> act_spacesEl = mul->GetActiveApproxSpaces();
    int initialconnect = 0;
    for (int iespacetes=0; iespacetes<act_spacesEl.size(); iespacetes++) {
        if (act_spacesEl[iespacetes]==0) {
            continue;
        }
        TPZCompEl *Atomic_celfrom = mul->Element(iespacetes);
        if (!Atomic_celfrom) {
            continue;
        }
        int nconnects = Atomic_celfrom->NConnects();
        for (int icon = 0 ; icon < nconnects; icon++){
            TPZConnect &connectAtomic = Atomic_celfrom->Connect(icon);
            int64_t seqnumberAtomic = connectAtomic.SequenceNumber();
            TPZConnect &connectMF = mul->Connect(initialconnect + icon);
            int64_t seqnumberMF = connectMF.SequenceNumber();
            Match currentmatch;
            currentmatch.fblocknumber = seqnumberMF;
            TPZBlock *blockAto = &Atomic_celfrom->Mesh()->Block();
            std::pair<TPZBlock *, int64_t> target = std::make_pair(blockAto,seqnumberAtomic);
            currentmatch.fblockTarget = target;
            fconnecttransfer.insert(currentmatch);
        }
        initialconnect += nconnects;
    }
}


/// Insert the connect matches for an element group
void TPZMFSolutionTransfer::MeshTransferData::InsertMatches(TPZElementGroup *celgr)
{
    const TPZVec<TPZCompEl *> &elgr = celgr->GetElGroup();
    for(auto cel : elgr)
    {
        InsertMatches(cel);
    }
    
}

/// Insert the connect matches for a condensed computational element
void TPZMFSolutionTransfer::MeshTransferData::InsertMatches(TPZCondensedCompEl *cond)
{
    TPZCompEl *ref = cond->ReferenceCompEl();
    InsertMatches(ref);
}


void TPZMFSolutionTransfer::MeshTransferData::BuildTransferData(TPZCompMesh* cmesh){

    TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh);
    
    TPZCompMesh * targetmesh = 0;
    
    if(subcmesh){
        targetmesh = subcmesh;
    }
    
    if(multcmesh){
        targetmesh = multcmesh;
    }

    int nels = targetmesh->NElements();
    for (int iel =0; iel <nels; iel++)
    {
        TPZCompEl *celtarget = cmesh->Element(iel);
        if(!celtarget) continue;
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(celtarget);
        if(sub && subcmesh)
        {
            std::cout << "Nested submeshes are not supported at this point\n";
            DebugStop();
        }
        InsertMatches(celtarget);
    }
}


void TPZMFSolutionTransfer::TransferFromMultiphysics(){
    
    int nsolutionstransfers = fmeshTransfers.size();
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        fmeshTransfers[isoltrans].TransferFromMultiphysics();
    }
    
}
void TPZMFSolutionTransfer::TransferToMultiphysics(){

    int nsolutionstransfers = fmeshTransfers.size();
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        fmeshTransfers[isoltrans].TransferToMultiphysics();
    }
}

void TPZMFSolutionTransfer::BuildTransferData(TPZCompMesh* cmesh){
    
    MeshTransferData transdata;
    transdata.fcmesh_ori = cmesh;
    transdata.BuildTransferData(cmesh);
    fmeshTransfers.push_back(transdata);
    int nels = cmesh->NElements();
    for (int iel=0; iel<nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        TPZSubCompMesh * subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            continue;
        }
        MeshTransferData transdatasub;
        transdatasub.fcmesh_ori =subcmesh;
        transdatasub.BuildTransferData(subcmesh);
        fmeshTransfers.push_back(transdatasub);
    }
}
