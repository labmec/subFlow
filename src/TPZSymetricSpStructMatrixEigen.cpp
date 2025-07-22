//
//  TPZSymetricSpStructMatrixEigenEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#include "TPZSymetricSpStructMatrixEigen.h"
// #include "pz_config.h"
// #include "threads.h"


#include "pzgeoelbc.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

#include "TPZLinearAnalysis.h"
//#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzdxmesh.h"
#include <fstream>

#include "pzelmat.h"

#include "TPZSYSMPMatrix.h"
//#include "pzsysmp.h"
#include "TPZBndCondT.h"
#include "TPZTimer.h"
#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

#include <typeinfo>
using namespace std;

TPZStructMatrix * TPZSymetricSpStructMatrixEigen::Clone(){
    return new TPZSymetricSpStructMatrixEigen(*this);
}
TPZMatrix<STATE> * TPZSymetricSpStructMatrixEigen::Create(){
    
    /**
     *Longhin implementation
     */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    TPZMatrix<STATE> * mat = SetupMatrixData(elgraph, elgraphindex);
    return mat;
}

TPZMatrix<STATE> * TPZSymetricSpStructMatrixEigen::SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex){
    
    int64_t neq = fEquationFilter.NActiveEquations();
    TPZSYsmpMatrixEigen<STATE> * mat = new TPZSYsmpMatrixEigen<STATE>(neq,neq);
    
    /**Creates a element graph*/
    TPZRenumbering renum;
    renum.SetElementsNodes(elgraphindex.NElements() -1 ,fMesh->NIndependentConnects());
    renum.SetElementGraph(elgraph,elgraphindex);
    
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    renum.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    int64_t i;
    int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(i=0;i<nblock;i++){
        int64_t iblsize = fMesh->Block().Size(i);
        int64_t iblpos = fMesh->Block().Position(i);
        int64_t numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        totaleq += iblsize;
        int64_t icfirst = nodegraphindex[i];
        int64_t iclast = nodegraphindex[i+1];
        int64_t j;
        //longhin
        totalvar+=(iblsize*(iblsize+1))/2;
        for(j=icfirst;j<iclast;j++) {
            int64_t col = nodegraph[j];
            if (col < i) {
                continue;
            }
            
            if (col == i) {
                DebugStop();
            }
            
            int64_t colsize = fMesh->Block().Size(col);
            int64_t colpos = fMesh->Block().Position(col);
            int64_t numactive = fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
            totalvar += iblsize*colsize;
        }
    }
    
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
    
    nblock=fMesh->NIndependentConnects();
    
    // number of non-zeros
    int64_t nnzeros = 0;
    TPZVec<int64_t> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<STATE> EqValue(totalvar,0.);
    for(i=0;i<nblock;i++){
        int64_t iblsize = fMesh->Block().Size(i);
        int64_t iblpos = fMesh->Block().Position(i);
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t i=0; i<iblsize; i++) {
            rowdestindices[i] = iblpos+i;
        }
        fEquationFilter.Filter(rowdestindices);
        
        int64_t ibleq;
        // working equation by equation
        for(ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
            Eq[ieq] = pos;
            int64_t colsize,colpos,jbleq;
            int64_t diagonalinsert = 0;
            int64_t icfirst = nodegraphindex[i];
            int64_t iclast = nodegraphindex[i+1];
            int64_t j;
            for(j=icfirst;j<iclast;j++)
            {
                int64_t col = nodegraph[j];
                if (col < i) {
                    continue;
                }
                // force the diagonal block to be inserted
                // the nodegraph does not contain the pointer to itself
                if(!diagonalinsert && col > i)
                {
                    diagonalinsert = 1;
                    int64_t colsize = fMesh->Block().Size(i);
                    int64_t colpos = fMesh->Block().Position(i);
                    TPZManVector<int64_t> destindices(colsize);
                    for (int64_t i=0; i<colsize; i++) {
                        destindices[i] = colpos+i;
                    }
                    fEquationFilter.Filter(destindices);
                    int64_t jbleq;
                    for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                        //             if(colpos+jbleq == ieq) continue;
                        int64_t jeq = destindices[jbleq];
                        if (jeq < ieq) {
                            nnzeros++;
                            continue;
                        }
                        EqCol[pos] = destindices[jbleq];
                        EqValue[pos] = 0.;
                        //            colpos++;
                        pos++;
                        nnzeros++;
                    }
                }
                colsize = fMesh->Block().Size(col);
                colpos = fMesh->Block().Position(col);
                if (fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
                for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                    int64_t jeq = destindices[jbleq];
                    if (jeq < ieq) {
                        nnzeros++;
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    colpos++;
                    pos++;
                    nnzeros++;
                }
            }
            // all elements are below (last block certainly)
            if(!diagonalinsert)
            {
                diagonalinsert = 1;
                int64_t colsize = fMesh->Block().Size(i);
                int64_t colpos = fMesh->Block().Position(i);
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
                int64_t jbleq;
                for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                    //             if(colpos+jbleq == ieq) continue;
                    int64_t jeq = destindices[jbleq];
                    if (jeq < ieq) {
                        nnzeros++;
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    //            colpos++;
                    nnzeros++;
                    pos++;
                }
            }
            ieq++;
        }
    }
    
    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
//    m_triplets.resize(nnzeros);
    m_triplets.clear();
    
    return mat;
}


void TPZSymetricSpStructMatrixEigen::Serial_Assemble(TPZBaseMatrix & stiffnessB, TPZBaseMatrix & rhsB) {
    
    
    TPZMatrix<STATE> &stiffness = dynamic_cast<TPZMatrix<STATE>&>(stiffnessB);
    TPZFMatrix<STATE> &rhs = dynamic_cast<TPZFMatrix<STATE>&>(rhsB);
    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matRed = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (&stiffness);
    if (mat) {
        Serial_AssembleGlob(stiffness,rhs);
    }
    if (matRed) {
         Serial_AssembleSub(stiffness,rhs);
        
//        boost::posix_time::ptime timeactual=boost::posix_time::microsec_clock::local_time();
//        std::cout<<"TotalSolveSub: "<<timetotalSolve<<std::endl;
//        std::cout<<"TotalSubMeshes: "<<timetotal<<std::endl;
//        std::cout<<"TotalCalcStiff: "<<timeTotalCalcStiff<<std::endl;
//        std::cout<<"TotalAssemble: "<<timeTotalAssemble<<std::endl;
//        std::cout<<"TotalSetFromTriplets: "<<timeTotalSetFromTrriplets<<std::endl;
//        timetotal=timeactual-timeactual;
//        timetotalSolve=timeactual-timeactual;
//        timeTotalCalcStiff=timeactual-timeactual;
//        timeTotalAssemble=timeactual-timeactual;
//        timeTotalSetFromTrriplets=timeactual-timeactual;

        
    }
    
}
void TPZSymetricSpStructMatrixEigen::Serial_AssembleSub(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs) {
    
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matRed = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (&stiffness);
//    matRed->Zero();
    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrixT<STATE> ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
    
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
  //  TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh*>(fMesh);
  
    int64_t count = 0;
    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        TPZFastCondensedElement *fastcondEl = dynamic_cast<TPZFastCondensedElement *>(el);
        if (!el) continue;
        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = fMaterialIds.size();
        if(matidsize){
            if(!el->NeedsComputing(fMaterialIds)) continue;
        }
        
        count++;
        if (!(count % 1000)) {
            std::cout << '*';
            std::cout.flush();
        }
        if (!(count % 20000)) {
            std::cout << "\n";
        }
     
        calcstiff.start();
        ek.Reset();
        ef.Reset();
        el->CalcStiff(ek, ef);
        calcstiff.stop();
      
//        subcmesh->fTimeTotalCalcStiff += endtimeCalc-initimeCalc;
        
        TPZMatrix<STATE> * matpzmat=matRed->K00().operator->();
        TPZSYsmpMatrixEigen<STATE> *matk00eigen =dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmat);
        TPZFMatrix<REAL>* matk01Eigen = &matRed->K01();
        TPZFMatrix<REAL> *matk10Eigen = &matRed->K10();
        TPZFMatrix<STATE> *matpzmat11= &matRed->K11();
//        matk00eigen->Zero();
//        matk01Eigen->Zero();
//        matk10Eigen->Zero();
//        matk10Eigen->Zero();
//        matpzmat11->Zero();
        
        
        assemble.start();
     
        bool hasIndex=false;
        if (!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            int64_t nelems = ek.fSourceIndex.NElements();
            int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
            REAL prevval;
            
            for(icoef=0; icoef<nelems; icoef++) {
                ieq = ek.fDestinationIndex[icoef];
                ieqs = ek.fSourceIndex[icoef];
                for(jcoef=icoef; jcoef<nelems; jcoef++) {
                    jeq = ek.fDestinationIndex[jcoef];
                    jeqs = ek.fSourceIndex[jcoef];
                    prevval = matRed->GetVal(ieq,jeq);
                    prevval += ek.fMat(ieqs,jeqs);
                    int r = ieq;
                    int c = jeq;
                    REAL value =prevval;
                    int64_t row(r),col(c);
                    int fDim0 = matRed->Dim0();
                   
                    
                    if ( row > col ) Swap( &row, &col );
                    if (row<fDim0 &&  col<fDim0)  {
//                        if (fastcondEl) {
//                            if (!fastcondEl->hasIndexes) {
//                                int posfa =0;
//                                matk00eigen->PutVal(row,col,value, posfa);
//                                if (fastcondEl) {
//                                    fastcondEl->faValK00.push_back(posfa);
//                                    hasIndex=true;
//                                }
//                            }
//                            else{
//                                int posfa =fastcondEl->faValK00[count];
//                                matk00eigen->PutVal(posfa,value);
//                                countpos++;
//                            }
//                        }
//                        else{
//                            matk00eigen->PutVal(row,col,value);
//                        }
                        matk00eigen->PutVal(row,col,value);
                    }
                    if (row<fDim0 &&  col>=fDim0)  {
                        matk01Eigen->PutVal(row,col-fDim0,value);
                    }
                    if (row>=fDim0 &&  col<fDim0)  {
                        matk10Eigen->PutVal(row-fDim0,col,value);
                    }
                    if (row>=fDim0 &&  col>=fDim0) {
                        matpzmat11->PutVal(row-fDim0,col-fDim0,value);
                    }
                }
            }
            
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ek.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//            matRed->AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            int64_t nelem = ek.fSourceIndex.NElements();
            int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
            
            REAL prevval;
                for(icoef=0; icoef<nelem; icoef++) {
                    ieq = ek.fDestinationIndex[icoef];
                    ieqs = ek.fSourceIndex[icoef];
                    for(jcoef=icoef; jcoef<nelem; jcoef++) {
                        jeq = ek.fDestinationIndex[jcoef];
                        jeqs = ek.fSourceIndex[jcoef];
                        prevval = matRed->GetVal(ieq,jeq);
                        prevval += ek.fConstrMat(ieqs,jeqs);
//                        matRed->PutVal(ieq,jeq,prevval);
                        //
                        int r = ieq;
                        int c = jeq;
                        REAL value =prevval;
                        int64_t row(r),col(c);
                        int fDim0 = matRed->Dim0();
                        TPZMatrix<STATE> * matpzmat=matRed->K00().operator->();
                        TPZSYsmpMatrixEigen<STATE> *matk00eigen =dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmat);
                        TPZFMatrix<REAL>* matk01Eigen = &matRed->K01();
                        TPZFMatrix<REAL> *matk10Eigen = &matRed->K10();
                        TPZFMatrix<STATE> *matpzmat11= &matRed->K11();
                        
                        if ( row > col ) Swap( &row, &col );
                        if (row<fDim0 &&  col<fDim0)  {
//                            if (!fastcondEl->hasIndexes) {
//                                int posfa =0;
//                                matk00eigen->PutVal(row,col,value, posfa);
//                                if (fastcondEl) {
//                                    fastcondEl->faValK00.push_back(posfa);
//                                    hasIndex=true;
//                                }
//                            }
//                            else{
//                                int posfa =fastcondEl->faValK00[count];
//                                matk00eigen->PutVal(posfa,value);
//                                countpos++;
//                            }
                            matk00eigen->PutVal(row,col,value);
                        }
                        if (row<fDim0 &&  col>=fDim0)  {
                            matk01Eigen->PutVal(row,col-fDim0,value);
                        }
                        if (row>=fDim0 &&  col<fDim0)  {
                            matk10Eigen->PutVal(row-fDim0,col,value);
                        }
                        if (row>=fDim0 &&  col>=fDim0) {
                            matpzmat11->PutVal(row-fDim0,col-fDim0,value);
                        }
                    }
                }
            rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
        }
        if (hasIndex==true) {
            fastcondEl->hasIndexes=true;
        }
        
        assemble.stop();
        

    }//fim for iel
}


//
void  TPZSymetricSpStructMatrixEigen::Serial_AssembleGlob(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs) {
    

    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrixT<STATE> ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
    boost::posix_time::ptime actime = boost::posix_time::microsec_clock::local_time();
    fAsTotalAdkelsSub =actime-actime;
    fAsTotalCalcStifSub =actime-actime;
    fAsTotaCondensedSub =actime-actime;;
    fAsTotaAssembleSub =actime-actime;
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

    int64_t count = 0;
    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
    mat->Zero();
    std::ofstream filep("Mats.txt");
    std::ofstream filep2("DestIndexes.txt");
        for (iel = 0; iel < nelem; iel++) {
            TPZCompEl *el = elementvec[iel];
//            TPZFastCondensedElement *cond = dynamic_cast<TPZFastCondensedElement *>(el);
//            if(cond){
//                TPZCompEl *celm = cond->GetMultiphysics();
//                calcstiff.start();
//                ek.Reset();
//                ef.Reset();
//                celm->CalcStiff(ek, ef);
//            }
            TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(el);
            if (!el) continue;
            int matid = 0;
            TPZGeoEl *gel = el->Reference();
            
            
            if (gel) {
                matid = gel->MaterialId();
              
            }
            int matidsize = fMaterialIds.size();
            if(matidsize){
                if(!el->NeedsComputing(fMaterialIds)) continue;
            }
            
            count++;
            if (!(count % 1000)) {
                std::cout << '*';
                std::cout.flush();
            }
            if (!(count % 20000)) {
                std::cout << "\n";
            }
            calcstiff.start();
            ek.Reset();
            ef.Reset();
            
//            if(gel->Index()==61){
//                int ok=0;
//            }
//            if(iel == 70486)
//            {
//                std::ofstream sout("submesh.txt");
//                el->Print(sout);
//            }
//            std::cout << "Computing element " << iel <<  typeid(el).name() << std::endl;
            el->CalcStiff(ek, ef);
//            filep<<"iel: "<<gel->Index()<< " MatId: "<<gel->MaterialId()<<std::endl;
//                ek.fMat.Print("Ekmat= ", filep, EMathematicaInput);
//            
            
            
            calcstiff.stop();
            assemble.start();
            
            if (!ek.HasDependency()) {
                ek.ComputeDestinationIndices();
//                filep2<<"iel: "<< gel->Index()<<" "<<ek.fDestinationIndex<<std::endl;
                
                
                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
                        int64_t nelem = ek.fSourceIndex.NElements();
                        int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
                        REAL prevval;
                        for(icoef=0; icoef<nelem; icoef++) {
                            ieq = ek.fDestinationIndex[icoef];
                            ieqs = ek.fSourceIndex[icoef];
                            for(jcoef=icoef; jcoef<nelem; jcoef++) {
                                jeq = ek.fDestinationIndex[jcoef];
                                jeqs = ek.fSourceIndex[jcoef];
                                prevval = mat->GetVal(ieq,jeq);
                                prevval += ek.fMat(ieqs,jeqs);
                                mat->PutVal(ieq,jeq,prevval);
                                
                                
                            }
                        }
                 rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            assemble.stop();
        }//fim for iel
            else {
                // the element has dependent nodes
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
                        int64_t nelem = ek.fSourceIndex.NElements();
                        int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
                        REAL prevval;
                        for(icoef=0; icoef<nelem; icoef++) {
                            ieq = ek.fDestinationIndex[icoef];
                            ieqs = ek.fSourceIndex[icoef];
                            for(jcoef=icoef; jcoef<nelem; jcoef++) {
                                jeq = ek.fDestinationIndex[jcoef];
                                jeqs = ek.fSourceIndex[jcoef];
                                prevval = mat->GetVal(ieq,jeq);
                                prevval += ek.fConstrMat(ieqs,jeqs);
                                mat->PutVal(ieq,jeq,prevval);
                            }
                        }
                rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            }
//        if (count > 1000) std::cout << std::endl;
            if (subcmesh) {
//                fAsTotalAdkelsSub +=subcmesh->fTimeTotalAddKels;
//                fAsTotalCalcStifSub += subcmesh->fTimeTotalCalcStiff;
//                fAsTotaAssembleSub += subcmesh->fTimeAssemble;
//                fAsTotaCondensedSub += subcmesh->fTimeCondensed;
            }
    }
//    std::cout<<"SubAssembleTime: "<<fAsTotaAssembleSub<<std::endl;
//    std::cout<<"SubCondensedTime: "<<fAsTotaCondensedSub<<std::endl;
//    std::cout<<"SubCalcStiffTime: "<<fAsTotalCalcStifSub<<std::endl;
//    std::cout<<"SubAddKelTime: "<<fAsTotalAdkelsSub<<std::endl;
    
//    std::ofstream fileg("kg.txt");
//    stiffness.Print(fileg);
}


int TPZSymetricSpStructMatrixEigen::ClassId() const{
    return Hash("TPZSymetricSpStructMatrixEigen") ^
        TPZStructMatrixT<STATE>::ClassId() << 1 ^
        TPZStructMatrixOR<STATE>::ClassId() << 2;
}

void TPZSymetricSpStructMatrixEigen::Read(TPZStream& buf, void* context){
    TPZStructMatrixT<STATE>::Read(buf,context);
    TPZStructMatrixOR<STATE>::Read(buf,context);
}

void TPZSymetricSpStructMatrixEigen::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrixT<STATE>::Write(buf,withclassid);
    TPZStructMatrixOR<STATE>::Write(buf,withclassid);
}
//
//#ifndef STATE_COMPLEX
//#include "pzmat2dlin.h"
//
//int TPZSymetricSpStructMatrixEigen::main() {
//    
//    int refine=5;
//    int order=5;
//    
//    TPZGeoMesh gmesh;
//    TPZCompMesh cmesh(&gmesh);
//    double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},
//        {0.,1.,0.}};
//    
//    int i,j;
//    TPZVec<REAL> coord(3,0.);
//    for(i=0; i<4; i++) {
//        // initializar as coordenadas do no em um vetor
//        for (j=0; j<3; j++) coord[j] = coordstore[i][j];
//        
//        // identificar um espa� no vetor onde podemos armazenar
//        // este vetor
//        gmesh.NodeVec ().AllocateNewElement ();
//        
//        // initializar os dados do n�       gmesh.NodeVec ()[nodeindex].Initialize (i,coord,gmesh);
//    }
//    int el;
//    TPZGeoEl *gel;
//    for(el=0; el<1; el++) {
//        
//        // initializar os indices dos nos
//        TPZVec<int64_t> indices(4);
//        for(i=0; i<4; i++) indices[i] = i;
//        // O proprio construtor vai inserir o elemento na malha
//        //       gel = new TPZGeoElQ2d(el,indices,1,gmesh);
//        int64_t index;
//        gel = gmesh.CreateGeoElement(EQuadrilateral,indices,1,index);
//    }
//    gmesh.BuildConnectivity ();
//    
//    TPZVec<TPZGeoEl *> subel;
//    
//    cout << "Refinement ";
//    cin >> refine;
//    cout << endl;
//    DebugStop();
//    //    UniformRefine(refine,gmesh);
//    
//    TPZGeoElBC gelbc(gel,4,-4);
//    TPZMat2dLin *meumat = new TPZMat2dLin(1);
//    TPZFMatrix<STATE> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
//    meumat->SetMaterial (xk,xc,xf);
//    TPZMaterial * meumatptr(meumat);
//    cmesh.InsertMaterialObject(meumatptr);
//    
//    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
//    TPZMaterial * bnd = meumat->CreateBC (meumatptr,-4,0,val1,val2);
//    cmesh.InsertMaterialObject(bnd);
//    
//    cout << "Interpolation order ";
//    cin >> order;
//    cout << endl;
//    
//    //TPZCompEl::gOrder = order;
//    cmesh.SetDefaultOrder(order);
//    
//    cmesh.AutoBuild();
//    //    cmesh.AdjustBoundaryElements();
//    cmesh.InitializeBlock();
//    
//    ofstream output("outputPar.dat");
//    TPZAnalysis an(&cmesh,true,output);
//    
//    //TPZVec<int> numelconnected(cmesh.NEquations(),0);
//    TPZSymetricSpStructMatrixEigen mat(&cmesh);
//    
//    an.SetStructuralMatrix(mat);
//    
//    TPZStepSolver<STATE> sol;
//    sol.SetJacobi(100,1.e-5,0);
//    
//    
//    an.SetSolver(sol);
//    an.Run(output);
//    output.flush();
//    cout.flush();
//    return 0;
//    
//}
//#endif
