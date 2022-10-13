
#ifdef PETSC
#include "PetscUtil.h"
#include <petsc.h>
#include <petscviewerhdf5.h>

#include "Home.h"
#include "Comm.h"

//#define GHOST 

PetscErrorCode FormVKFunction(SNES snes,Vec velocities,Vec F,void *ptr);
void CaculateTensorLengthXBij(Home_t *home, Node_t *nodei, Node_t *nodej, real8 (*lXB)[3]);
PetscErrorCode FormVKJacobianNolinear(SNES snes,Vec X,Mat J,Mat jac,void *ptr);
PetscErrorCode FormVKJacobianLinear(Home_t *home, PetscVKCtx_t *vk);

/**
 * PetscVelocityKernelInitialize: initializes the velocity kernel context
 * @param vk a pointer to the PetscVKCtx_t structure that contains all the parameters for the velocity
 * kernel.
 */
PetscErrorCode PetscVelocityKernelInitialize(PetscVKCtx_t *vk)
{
    PetscErrorCode ierr;
    
    PetscFunctionBegin;

    PetscCall(PetscPrintf(MPI_COMM_WORLD,"Using Petsc velocity kernel to solve the dislocation velocities\n"));

    vk->nolinear      = PETSC_FALSE;
    vk->ltog          = PETSC_FALSE;
    vk->matFDColoring = PETSC_FALSE;

    /**
     * Defines the velocity kernel option -vk_nolinear with a default value of no, i. e., a linear solver is used.
     */
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "vk_" ,"options for velocity kernel",__FILE__); PetscCall(ierr);
    PetscCall(PetscOptionsBool("-nolinear", "Using nolinear solver (SNES)",__func__,vk->nolinear,&vk->nolinear,NULL));
    PetscCall(PetscOptionsBool("-ltog", "Assemble Jacobi matrix by using local index mapping global index",
                                __func__,vk->ltog,&vk->ltog,NULL));
    PetscCall(PetscOptionsBool("-jac_coloring", "Set Jacobi matrix coloring",
                                __func__,vk->matFDColoring,&vk->matFDColoring,NULL));
    ierr = PetscOptionsEnd(); PetscCall(ierr);

    PetscFunctionReturn(0);
}

/**
 * > This function assembles the global velocity kernel matrix and the global velocity vector
 * 
 * @param home the home structure
 */
PetscErrorCode AssembleGlobalVelocityKernel(Home_t *home)
{

    const int       bs = 3;
    int             i, j, nbrID, rows[3], cols[6], rowID, colID;
    int             localID, cellNum, cellID, maxNodes, totNodes;  
    real8           lXB[3][3], values[18];
    Node_t          *node, *nbr;
    int             *nNbrsNativeList, *nNbrsGhostList, *nnz;
    Cell_t          *cell;
#ifdef GHOST
    PetscInt     *ifrom;
#else
    IS           islocal, isglobal;
#endif
#ifdef  PETSC_DEBUG
    FILE  *fp;
#endif


    PetscVKCtx_t       *vk = home->vk;
    const int    &nNatives = vk->nNatives;
    const int    &nGlobals = vk->nGlobals;
    int          &nGhosts  = vk->nGhosts;
    Vec        &velocities = vk->velocities;
    Vec            &forces = vk->forces;
    int    homeNativeCells = home->nativeCellCount;
    int           myDomain = home->myDomain;
    ISLocalToGlobalMapping  &isl2g = vk->isl2g;


    PetscFunctionBegin;
    
#ifdef PETSC_DEBUG
    PetscCall(PetscViewerASCIIGetPointer(vk->viewer,&fp));
#endif
    
    /*
        1. Count the maximum number of nodes required for the current domain 
           (maxNodes), i.e. assuming that all ghost nodes are not the same
        2. Count the neighboring node components (native, ghost and total) of each 
           native node, which is used to assemble the adjacency matrix.
        3. Check the existence of adjacent nodes of all native nodes
    */
    PetscCall(PetscCalloc3(bs*nNatives,&nNbrsNativeList,bs*nNatives,&nNbrsGhostList,bs*nNatives,&nnz));
    for(i=0;i<nNatives*bs;i++){
        nnz[i] = bs; nNbrsNativeList[i++] = bs;
        nnz[i] = bs; nNbrsNativeList[i++] = bs;
        nnz[i] = bs; nNbrsNativeList[i]   = bs;
    }
    maxNodes = nNatives; localID = 0;
    for(cellNum=0;cellNum<homeNativeCells;cellNum++){
        cellID = home->cellList[cellNum];
        cell   = home->cellKeys[cellID];
        for(node=cell->nodeQ;node!=(Node_t *)NULL; node=node->nextInCell){
            if(node->myTag.domainID!=myDomain)continue;
            for(i=0;i<node->numNbrs;i++){
                nbr = GetNeighborNode(home, node, i);
                if(nbr == (Node_t *)NULL){
                    PrintNode(home,node);
                    Fatal("%s: can not find the nbr node %d\n",__func__,i);
                }
                if(nbr->myTag.domainID==myDomain){
                    nNbrsNativeList[3*localID] += bs;
                    nNbrsNativeList[3*localID+1] += bs;
                    nNbrsNativeList[3*localID+2] += bs;
                }else{
                    nNbrsGhostList[3*localID] += bs;
                    nNbrsGhostList[3*localID+1] += bs;
                    nNbrsGhostList[3*localID+2] += bs;
                    maxNodes += 1;
                }
                nnz[3*localID]   += bs;
                nnz[3*localID+1] += bs;
                nnz[3*localID+2] += bs;
            }
            localID++;
        }
    }

    /*
       1. Count the number of all nodes (natives and ghosts) in this domain
           and compress the node list of global index.
       2. Record the pointers to all nodes and assemble to vk->nodeList.
    */
    PetscCall(PetscMalloc4(maxNodes,&vk->gloIDList,maxNodes,&vk->nodeList,nNatives,
                &vk->nCstList,bs*nNatives,&vk->cstList));
    localID = 0;
    for(cellNum=0;cellNum<homeNativeCells;cellNum++){
        cellID = home->cellList[cellNum];
        cell   = home->cellKeys[cellID];
        for(node=cell->nodeQ;node!=(Node_t *)NULL; node=node->nextInCell){
            if(node->myTag.domainID!=myDomain){continue;}
            node->locID = localID;
            vk->nodeList[localID]    = node;
            vk->gloIDList[localID++] = node->gloID;
        }
    }
    for(cellNum=0;cellNum<homeNativeCells;cellNum++){
        cellID = home->cellList[cellNum];
        cell   = home->cellKeys[cellID];
        for(node=cell->nodeQ;node!=(Node_t *)NULL; node=node->nextInCell){
            if(node->myTag.domainID!=myDomain){continue;}
            for(i=0;i<node->numNbrs;i++){
                nbr = GetNeighborNode(home, node, i);
                if(nbr->myTag.domainID!=myDomain){
                    for(j=nNatives;j<localID;j++)if(nbr->gloID==vk->gloIDList[j])break;
                    if(j==localID){
                        nbr->locID = localID;
                        vk->nodeList[localID]    = nbr;
                        vk->gloIDList[localID++] = nbr->gloID;
                    }
                }
            }
        }
    }
    totNodes = localID;
    nGhosts  = totNodes-nNatives;
    
#ifdef GHOST
    printf("[%d] 3: nNatives=%d tot=%d  nGhost=%d\n",home->myDomain,nNatives,totNodes,nGhosts);
    PetscCall(PetscCalloc1(bs*nGhosts,&ifrom));
    for(i=nNatives;i<totNodes;i++){
        node    = vk->nodeList[i];
        localID = i-nNatives;
        ifrom[localID*bs]   = bs*node->gloID;
        ifrom[localID*bs+1] = bs*node->gloID+1;
        ifrom[localID*bs+2] = bs*node->gloID+2;
    }

    /*
      Create vector data structures
    */
    PetscCall(VecCreateGhost(PETSC_COMM_WORLD,bs*nNatives,bs*nGlobals,bs*nGhosts,ifrom,&velocities));
    PetscCall(VecSetOptionsPrefix(velocities,"vk_velocities_"));
    PetscCall(PetscObjectSetName((PetscObject) velocities,"velocities"));
    PetscCall(PetscFree(ifrom));    
#else
    PetscCall(VecCreate(MPI_COMM_WORLD,&velocities));
    PetscCall(VecSetOptionsPrefix(velocities,"vk_velocities_"));
    PetscCall(PetscObjectSetName((PetscObject) velocities,"velocities"));
    PetscCall(VecSetFromOptions(velocities));
    PetscCall(VecSetSizes(velocities,bs*nNatives,bs*nGlobals));
    PetscCall(VecSetBlockSize(velocities, bs));
#endif
    PetscCall(VecDuplicate(velocities,&forces));

#ifndef GHOST
    PetscCall(VecCreateSeq(MPI_COMM_SELF,bs*totNodes,&vk->locVel));
    
    /*
      Create the scatter between the global representation and the
      local representation
    */
    PetscCall(ISCreateStride(MPI_COMM_SELF,bs*totNodes,0,1,&islocal));
    PetscCall(ISCreateBlock(MPI_COMM_SELF,bs,totNodes,vk->gloIDList,PETSC_COPY_VALUES,&isglobal));
    PetscCall(VecScatterCreate(velocities,isglobal,vk->locVel,islocal,&vk->scatter));
    PetscCall(ISDestroy(&isglobal));
    PetscCall(ISDestroy(&islocal));
#endif

#ifdef PETSC_DEBUG
    /*
      Print all nodes of each domain in order, including ghost nodes.
    */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"cycle=%d, nGlobals=%d, global indexs:\n",home->cycle,nGlobals));
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] # %d|%d:",myDomain,nNatives,nGhosts));
    for(i=0;i<nNatives;i++){
        PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d",vk->gloIDList[i]));
    }
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD," |"));
    for(;i<totNodes;i++){
        PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d",vk->gloIDList[i]));
    }
    PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n"));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,NULL));
#endif

    if(vk->nolinear){
        /*
            Set up the velocity coeffient matrix (linear)
        */
        PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF,nNatives*bs,totNodes*bs,0,nnz, &vk->velCoeMat));
        PetscCall(MatSetOptionsPrefix(vk->velCoeMat,"vk_coef_"));
        PetscCall(PetscObjectSetName((PetscObject) vk->velCoeMat,"coefficient matrix"));
        PetscCall(MatSetFromOptions(vk->velCoeMat));
        PetscCall(MatSetBlockSize(vk->velCoeMat,bs));
        PetscCall(MatSetUp(vk->velCoeMat));

        for(i=0;i<nNatives;i++){
            rowID = bs*i;
            node = vk->nodeList[i];
            for(j=0;j<node->numNbrs;j++){
                nbr = GetNeighborNode(home, node, j);
                PetscCheck(nbr != (Node_t *)NULL,PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
                "nbr arm %d of node (%d,%d) is null pointer",
                j,node->myTag.domainID,node->myTag.index);
                nbrID = nbr->locID; colID=3*nbrID;
                if(nbrID<i)continue;
                CaculateTensorLengthXBij(home,node,nbr,lXB);
                values[0] = 2.0*lXB[0][0]; values[1] = 2.0*lXB[0][1]; values[2] = 2.0*lXB[0][2];
                values[3] = lXB[0][0];     values[4] = lXB[0][1];     values[5] = lXB[0][2];
       
                values[6] = 2.0*lXB[1][0]; values[7]  = 2.0*lXB[1][1]; values[8]  = 2.0*lXB[1][2];
                values[9] = lXB[1][0];     values[10] = lXB[1][1];     values[11] = lXB[1][2];
                
                values[12] = 2.0*lXB[2][0]; values[13] = 2.0*lXB[2][1]; values[14] = 2.0*lXB[2][2];
                values[15] = lXB[2][0];     values[16] = lXB[2][1];     values[17] = lXB[2][2];
       
                rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
                cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
                cols[3] = colID; cols[4] = colID+1; cols[5] = colID+2;
                PetscCall(MatSetValues(vk->velCoeMat,3,rows,6,cols,values,ADD_VALUES));

                rows[0] = colID; rows[1] = colID+1; rows[2] = colID+2;
                cols[0] = colID; cols[1] = colID+1; cols[2] = colID+2;
                cols[3] = rowID; cols[4] = rowID+1; cols[5] = rowID+2;
                if(nbrID < nNatives){
                    PetscCall(MatSetValues(vk->velCoeMat,3,rows,6,cols,values,ADD_VALUES));
                }
            }
            rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
            cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
            NodalConstraintTensor(home, vk, i, values);
            PetscCall(MatSetValues(vk->velCoeMat,3,rows,3,cols,values,ADD_VALUES));
        }
        PetscCall(MatAssemblyBegin(vk->velCoeMat,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(vk->velCoeMat,MAT_FINAL_ASSEMBLY));
        #ifdef PETSC_DEBUG
            PetscCall(MatView(vk->velCoeMat,vk->viewer));
        #endif
        /*
        Tell the matrix we will never add a new nonzero location to the
        matrix. If we do, it will generate an error.
        */
        PetscCall(MatSetOption(vk->velCoeMat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE));

    }

    PetscCall(MatCreate(MPI_COMM_WORLD,&vk->velJacobi));
    PetscCall(MatSetOptionsPrefix(vk->velJacobi,"vk_jacobi_"));
    PetscCall(PetscObjectSetName((PetscObject) vk->velJacobi,"Jacobi matrix"));
    PetscCall(MatSetFromOptions(vk->velJacobi));
    PetscCall(MatSetSizes(vk->velJacobi,bs*nNatives,bs*nNatives,bs*nGlobals,bs*nGlobals));
    PetscCall(MatSetBlockSize(vk->velJacobi,bs));
    PetscCall(MatMPIAIJSetPreallocation(vk->velJacobi,0,nNbrsNativeList,0,nNbrsGhostList));
    PetscCall(MatSeqAIJSetPreallocation(vk->velJacobi,0,nNbrsNativeList));
    PetscCall(PetscFree3(nNbrsNativeList,nNbrsGhostList,nnz));
    PetscCall(MatSetUp(vk->velJacobi));

    if(vk->ltog){
        /*
        The following routine allows us to set the matrix values in local ordering
        */
        PetscCall(ISLocalToGlobalMappingCreate(MPI_COMM_SELF,bs,totNodes,vk->gloIDList,PETSC_COPY_VALUES,&isl2g));
        PetscCall(MatSetLocalToGlobalMapping(vk->velJacobi,isl2g,isl2g));
        #ifdef PETSC_DEBUG
            PetscCall(ISLocalToGlobalMappingView(isl2g,PETSC_VIEWER_STDOUT_SELF));
        #endif
    }
    PetscFunctionReturn(0);
}

/**
 * This function solves the velocity kernel equations using PETSc
 */
PetscErrorCode SolveGlobalVelocityKernel(Home_t *home)
{
    // variables for restore nodal velocities
    int                 i, totNodes;
    Node_t              *node, **nodeList;
    PetscScalar         **vArray, vel;

    //  Variables for SNES 
    SNES                    snes;
    SNESLineSearch          ls;
    Vec                     r;
    ISColoring              iscoloring;
    MatColoring             mc;

#ifdef PETSC_DEBUG
    PetscInt                nfails, its;
    Param_t                 *param = home->param;
    FILE                    *file;
#endif

    //  Variables for KSP
    KSP                     ksp;
    PC                      pc;
    KSPConvergedReason      reason;

    PetscScalar   vshr           = home->param->vshr;
    PetscVKCtx_t  *vk            = home->vk;
    Mat           &velJacobi     = vk->velJacobi;
    Vec           &velocities    = vk->velocities;
    Vec           &forces        = vk->forces;
    MatFDColoring &matfdcoloring = vk->matJacColor;

    PetscFunctionBegin;

    /*
        Update the initial velocities and forces 
    */
    PetscCall(UpdateVelocityForce(home, vk));

    if(vk->nolinear){
        PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
        PetscCall(SNESSetType(snes,SNESNEWTONLS));
        PetscCall(SNESGetLineSearch(snes,&ls));
//      PetscCall(SNESLineSearchSetType(ls,SNESLINESEARCHCP));
        PetscCall(SNESLineSearchSetType(ls,SNESLINESEARCHBASIC));
        PetscCall(SNESSetFromOptions(snes));

        PetscCall(VecDuplicate(velocities,&r));
        PetscCall(SNESSetFunction(snes,r,FormVKFunction,(void *)home));

        if(vk->matFDColoring){
            matfdcoloring = 0;
            /* Get the data structure of Jac */
            PetscCall(FormVKJacobianNolinear(snes,velocities,velJacobi,velJacobi,home));
            /* Create coloring context */
            PetscCall(MatColoringCreate(velJacobi,&mc));
            PetscCall(MatColoringSetType(mc,MATCOLORINGSL));
            PetscCall(MatColoringSetFromOptions(mc));
            PetscCall(MatColoringApply(mc,&iscoloring));
            PetscCall(MatColoringDestroy(&mc));
            PetscCall(MatFDColoringCreate(velJacobi,iscoloring,&matfdcoloring));
            PetscCall(MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))FormVKFunction,home));
            PetscCall(MatFDColoringSetFromOptions(matfdcoloring));
            PetscCall(MatFDColoringSetUp(velJacobi,iscoloring,matfdcoloring));
            /* PetscCall(MatFDColoringView(matfdcoloring,PETSC_VIEWER_STDOUT_WORLD)); */
            PetscCall(SNESSetJacobian(snes,velJacobi,velJacobi,SNESComputeJacobianDefaultColor,matfdcoloring));
            PetscCall(ISColoringDestroy(&iscoloring));
        }else{
            PetscCall(SNESSetJacobian(snes,velJacobi,velJacobi,FormVKJacobianNolinear,home));
        }

        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Now solve the nonlinear system
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        PetscCall(SNESSolve(snes,forces,velocities));
        #ifdef PETSC_DEBUG
            PetscCall(SNESGetIterationNumber(snes,&its));
            PetscCall(SNESGetNonlinearStepFailures(snes,&nfails));
            PetscCall(PetscPrintf(MPI_COMM_WORLD,"number of SNES iterations = %D, ",its));
            PetscCall(PetscPrintf(MPI_COMM_WORLD,"number of unsuccessful steps = %D\n",nfails));
        #endif
    
        /*
            Destory
        */
        PetscCall(SNESDestroy(&snes));
        PetscCall(VecDestroy(&r));

    }else{

        PetscCall(FormVKJacobianLinear(home,vk));
        PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
        PetscCall(KSPSetType(ksp,KSPGMRES));
        PetscCall(KSPSetOperators(ksp,vk->velJacobi,vk->velJacobi));
        
        /*  
        Set linear solver defaults for this problem (optional).
        - By extracting the KSP and PC contexts from the KSP context,
            we can then directly call any KSP and PC routines to set
            various options.
        - The following four statements are optional; all of these
            parameters could alternatively be specified at runtime via
            KSPSetFromOptions();
        */
        PetscCall(KSPGetPC(ksp,&pc));
        PetscCall(PCSetType(pc,PCJACOBI));
        PetscCall(KSPSetTolerances(ksp,1.0E-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT));

        PetscCall(KSPSetFromOptions(ksp));
        PetscCall(KSPSolve(ksp,forces,velocities));    
        
        PetscCall(KSPGetConvergedReason(ksp,&reason));
#ifdef PETSC_DEBUG
        PetscCheck(reason > 0,PetscObjectComm((PetscObject)ksp),PETSC_ERR_PLIB,"KSPConvergedReason() = %s [%d]",KSPConvergedReasons[reason],reason);
#else
        if(reason<0){
            PetscPrintf(PETSC_COMM_SELF,"KSPConvergedReason() = %s [%d]\n",KSPConvergedReasons[reason],reason);
        }
#endif
        PetscCall(KSPDestroy(&ksp));
    }

#ifdef PETSC_DEBUG
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"cycle=%d, velocity solution (b/s):\n",home->cycle));
    PetscCall(VecView(velocities,PETSC_VIEWER_STDOUT_WORLD));
#endif

    PetscCall(VecScatterBegin(vk->scatter, velocities, vk->locVel, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(vk->scatter, velocities, vk->locVel, INSERT_VALUES, SCATTER_FORWARD));

    nodeList = vk->nodeList;
    totNodes = vk->nNatives +vk->nGhosts;
    PetscCall(VecGetArray2d(vk->locVel,totNodes,3,0,0,&vArray));

#ifdef PETSC_DEBUG
    real8 deltaVec[3], deltaMag, dot1, vec1[3], vec2[3], dot2, vAmag,ratio;
    bool stop = false;
    Node_t *stopNode=(Node_t *)NULL;
    PetscCall(PetscViewerASCIIGetPointer(vk->viewer,&file));
    fprintf(file, "cycle = %d:\n", home->cycle);    
#endif

    for(i=0; i<vk->nNatives;i++){
        node = nodeList[i];

#ifdef PETSC_DEBUG
        node->vX = vArray[i][0];
        node->vY = vArray[i][1];
        node->vZ = vArray[i][2];

        AdjustNodeDirection(home, node, vArray[i]);
        deltaVec[0] = (vArray[i][0] - node->vX)*param->burgMag;
        deltaVec[1] = (vArray[i][1] - node->vY)*param->burgMag;
        deltaVec[2] = (vArray[i][2] - node->vZ)*param->burgMag;

        VectorCopy(vec1,node->vX,node->vY,node->vZ);
        vAmag = Normal(vec1)*param->burgMag;
        VECTOR_COPY(vec2, vArray[i]);
        NormalizeVec(vec1);
        NormalizeVec(vec2);
        dot1 = DotProduct(vec1,vec2);
        dot2 = dot_product(vec2[0],vec2[1],vec2[2],node->nx[0],node->ny[0],node->nz[0]);

        deltaMag  = Normal(deltaVec);
        ratio = (vAmag>1E-1)?deltaMag/vAmag:0.0;
        fprintf(file, "|%-5.2e|{%-5.4f,%-5.4f,%-5.4f}(%.4e,%.4e,%d) ", ratio, vAmag*vec1[0], vAmag*vec1[1],vAmag*vec1[2],dot1,dot2,node->constraint);
        if(i%5==4)fprintf(file,"\n");
        if(ratio>0.1){
            fprintf(file,"\n");
            PrintNode(home, node, file);
            PrintNbr(home, node, file);
            fprintf(file,"\n");
        }

        if(dot_product(node->fX,node->fY,node->fZ,node->vX,node->vY,node->vZ)<-1){
            stop = true;
            stopNode = node; 
        }
#endif
        if(vk->nCstList[i]==1){
            NormalComponent(vArray[i],&vk->cstList[3*i],vArray[i]);
        }else if(vk->nCstList[i]==2){
            ParallelComponent(vArray[i],&vk->cstList[3*i],vArray[i]);
        }else if(vk->nCstList[i]==3){
            VECTOR_ZERO(vArray[i]);
        }
        vel = Normal(vArray[i]);
        if(vel > vshr){
            vel = vshr/vel; 
        }else{
            vel = 1.0;
        }
        node->vX = vArray[i][0]*vel;
        node->vY = vArray[i][1]*vel;
        node->vZ = vArray[i][2]*vel;
    }

    PetscCall(VecRestoreArray2d(vk->locVel,totNodes,3,0,0,&vArray));

#ifdef PETSC_DEBUG
    fprintf(file,"\n");
    if(stop){
#if 0
        PetscCall(MatView(vk->velJacobi,PETSC_VIEWER_STDOUT_WORLD));
        PetscCall(VecScale(velocities,param->burgMag));
        PetscCall(VecView(velocities, PETSC_VIEWER_STDOUT_WORLD));
        PetscCall(VecScale(forces,param->burgMag));
        PetscCall(VecView(forces, PETSC_VIEWER_STDOUT_WORLD));
        PetscCall(VecScale(vk->locVel,param->burgMag));
        PetscCall(VecView(vk->locVel, PETSC_VIEWER_STDOUT_SELF));
#endif
        PrintNode(home,stopNode);
        printf("Inverse!\n");
    }
#endif

    PetscFunctionReturn(0);
}

/**
 * > This function destroys the velocity kernel context
 * \param[out] vk the velocity kernel context
 */
PetscErrorCode DestoryGlobalVelocityKernel(PetscVKCtx_t *vk)
{
    PetscFunctionBegin;
    PetscCall(PetscFree4(vk->gloIDList,vk->nodeList,vk->nCstList,vk->cstList));

    if(vk->nolinear){
        PetscCall(MatDestroy(&vk->velCoeMat));
        if(vk->matFDColoring)PetscCall(MatFDColoringDestroy(&vk->matJacColor));
    }

#ifndef GHOST
        PetscCall(VecDestroy(&vk->locVel));
        PetscCall(VecScatterDestroy(&vk->scatter));
#endif
    PetscCall(MatDestroy(&vk->velJacobi));
    if(vk->ltog)PetscCall(ISLocalToGlobalMappingDestroy(&vk->isl2g));
    
    PetscCall(VecDestroy(&vk->velocities));
    PetscCall(VecDestroy(&vk->forces));
    PetscFunctionReturn(0);
}

/**
 * The function `UpdateVelocityForce` is used to update the initial velocities and forces of the nodes
 * 
 * \param[in] home the home structure
 * \param[in] vk the velocity kernel
 */
PetscErrorCode UpdateVelocityForce(Home_t *home, PetscVKCtx_t *vk){

    Node_t          *node, *nbr;
    PetscReal       **velMat, **forceMat;

    int             ncs;
    real8           gc[3][3], len, ratio;
    bool            zeroLenJunc;

    const   int      bs = 3;
    int     *nCstList   = vk->nCstList;
    real8   *cstList    = vk->cstList;
    Node_t  **nodeList  = vk->nodeList;
    real8   minSeg      = home->param->minSeg;
#ifndef GHOST
    Vec     &velocities = vk->velocities;
    Vec     &forces     = vk->forces;

    PetscFunctionBegin;

    PetscCall(VecGetArray2dWrite(velocities,vk->nNatives,bs,0,0,&velMat));
    PetscCall(VecGetArray2dWrite(forces,vk->nNatives,bs,0,0,&forceMat));
#else
    Vec     velocities, forces;
    int     vecSize, nLocalNodes;

    PetscFunctionBegin;

    PetscCall(VecGhostGetLocalForm(vk->velocities,&velocities));
    PetscCall(VecGhostGetLocalForm(vk->forces,&forces));

    PetscCall(VecGetSize(velocities, &vecSize));
    nLocalNodes = vecSize/bs;
#ifdef PETSC_DEBUG
    if(nLocalNodes == vk->nNatives)PetscCall(PetscPrintf(PETSC_COMM_WORLD,"nLocalNodes == nNatives"));
    else if(nLocalNodes == vk->nGlobals)PetscCall(PetscPrintf(PETSC_COMM_WORLD,"nLocalNodes == nGlobals"));
    else{
        Fatal("nLocalNodes wrong");
    }
#endif
    PetscCall(VecGetArray2dWrite(velocities,nLocalNodes,bs,0,0,&velMat));
    PetscCall(VecGetArray2dWrite(forces,nLocalNodes,bs,0,0,&forceMat));
#endif

    /**
     * Currently only the initial velocities of the nodes are extracted directly. 
     * When a nonlinear mobility law is used, the initial velocities of the iteration can be 
     * predicted based on time-approximation and space-approximation principles.
     */
    for(int i=0; i<vk->nNatives; i++){
        node = nodeList[i];
        if ((node->constraint == PINNED_NODE) ||
            (NumNonZeroLenSegs(home,node)==0)){
            VECTOR_ZERO(velMat[i]);
            VECTOR_ZERO(forceMat[i]);
            nCstList[i] = 3; 
            VECTOR_ZERO(cstList+bs*i);
            continue;
        }else{
            VectorCopy(velMat[i],node->vX,node->vY,node->vZ);
            VectorCopy(forceMat[i],node->fX,node->fY,node->fZ);

            GetGlideConstraintList(home, node, &ncs, gc, zeroLenJunc);
            nCstList[i] = ncs;
            if(ncs==1){
                VECTOR_COPY(cstList+bs*i, gc[0]);
                NormalComponent(forceMat[i],gc[0],forceMat[i]);
            }else if (ncs==2){
                VECTOR_COPY(cstList+bs*i, gc[2]);
                ParallelComponent(forceMat[i],gc[2],forceMat[i]);
            }else{
                VECTOR_ZERO(cstList+bs*i);
                VECTOR_ZERO(forceMat[i]);
            }

            ratio = 1.0;
            for(int j=0;j<node->numNbrs;j++){
                nbr = GetNeighborNode(home,node,j);
                if(nbr==(Node_t *)NULL)continue;
                len = NodeDistance(home, node, nbr);
                if(len<0.5*minSeg)ratio *= len/minSeg;    
            }
            VectorTimesSelf(forceMat[i],ratio);
        }

#ifndef IMAGE_FORCE
        real8 imageForce[3], imag, iDir[3], ndot;
        NodeImageForce(home, node, imageForce); // nodal image force parallel to sement glide dirs
        imag = Normal(imageForce);
        if(imag > 1.0){ // if has image force
            iDir[0] = imageForce[0]/imag;
            iDir[1] = imageForce[1]/imag;
            iDir[2] = imageForce[2]/imag;
   
        ndot = DotProduct(forceMat[i], iDir); // nodal force magitude along image force dir
   
        if(ndot < -imag){ // if the dirs are oppsite and abs(ndot) > imag, add the image force
            CoordAdd(forceMat[i], imageForce, forceMat[i]);
        }else if(ndot < 0){ // if the dirs are oppsite and image > abs(ndot), animate the force conponent along the image force dir
            forceMat[i][0] -= ndot*iDir[0];
            forceMat[i][1] -= ndot*iDir[1];
            forceMat[i][2] -= ndot*iDir[2];
        }
    }
#endif
    }

#ifndef GHOST
    PetscCall(VecRestoreArray2dWrite(velocities,vk->nNatives,bs,0,0,&velMat));
    PetscCall(VecRestoreArray2dWrite(forces,vk->nNatives,bs,0,0,&forceMat));
#else 
    PetscCall(VecRestoreArray2dWrite(velocities,nLocalNodes,bs,0,0,&velMat));
    PetscCall(VecRestoreArray2dWrite(forces,nLocalNodes,bs,0,0,&forceMat));
    PetscCall(VecGhostRestoreLocalForm(vk->velocities,&velocities));
    PetscCall(VecGhostRestoreLocalForm(vk->forces,&forces));

    PetscCall(VecGhostUpdateBegin(vk->velocities,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGhostUpdateEnd(vk->velocities,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGhostUpdateBegin(vk->forces,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGhostUpdateEnd(vk->forces,INSERT_VALUES,SCATTER_FORWARD));

#ifdef PETSC_DEBUG
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"initial velocities"));
    PetscCall(VecView(vk->velocities, PETSC_VIEWER_STDOUT_WORLD));
    Fatal("Stop");
#endif

#endif
    PetscFunctionReturn(0);
}

void CaculateTensorLengthXBij_FCC(Home_t *home, Node_t *nodei, Node_t *nodej, real8 (*lXB)[3])
{
    Param_t         *param = home->param;
    int             armID, j, k;
    real8           d[3], n[3], m[3], b[3], mag2, length, oneSixthLen;
    real8           costheta, costheta2, sintheta2;
    real8           Bline, Beclimb, BlmBecl, Bscrew, Bedge, Bglide;
    real8           lineTensor[3][3], glideTensor[3][3], climbTensor[3][3];
    bool            tdKink;

    Bscrew  = 1.0/param->MobScrew;
    Bedge   = 1.0/param->MobEdge;
/*
 *  Climb is set this way to help with convergence and the fact that this
 *  is a glide restricted mobility
 */
    Beclimb = 10000.0 * Bedge;
    Beclimb = 0.0;
    Bline   = 1.0 * MIN(Bscrew, Bedge);
    BlmBecl = Bline - Beclimb;

    d[0] = nodej->x - nodei->x;
    d[1] = nodej->y - nodei->y;
    d[2] = nodej->z - nodei->z;
    ZImage(param,d,d+1,d+2);
    mag2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    if(mag2<1E-12){ Matrix33Init(lXB); return;}
    length = sqrt(mag2);
    d[0] /= length; d[1] /= length; d[2] /= length;
    oneSixthLen = length/6.0e0;

    armID = GetArmID(nodei,nodej);

    n[0] = nodei->nx[armID];
    n[1] = nodei->ny[armID];
    n[2] = nodei->nz[armID];

    b[0] = nodei->burgX[armID];
    b[1] = nodei->burgY[armID];
    b[2] = nodei->burgZ[armID];

    /*          
     *  Calculate how close to screw the arm is
     */ 
    costheta = (d[0]*b[0] + d[1]*b[1] + d[2]*b[2]);
    costheta2 = (costheta*costheta);
    sintheta2 = 1.0 - costheta2;
    if(sintheta2<0)sintheta2 = 0.0;

    Bglide = Bscrew*costheta2 + Bedge*sintheta2;

    cross(n, d, m);
    NormalizeVec(m);
    tdKink = (nodei->onGB > -1 && nodej->onGB > -1 && fabs(DotProduct(b, param->TBNormal))>5E-2);

    if((DotProduct(m,m) < 1.0E-2) || ((GetPlaneType_FCC(home, nodei, nodej)& GLIDE_PLANE)==0) || tdKink==true){
        lXB[0][0] = oneSixthLen * (d[0]*d[0]*BlmBecl + Beclimb);
        lXB[0][1] = oneSixthLen * (d[0]*d[1]*BlmBecl);
        lXB[0][2] = oneSixthLen * (d[0]*d[2]*BlmBecl);
        lXB[1][1] = oneSixthLen * (d[1]*d[1]*BlmBecl + Beclimb);
        lXB[1][2] = oneSixthLen * (d[1]*d[2]*BlmBecl);
        lXB[2][2] = oneSixthLen * (d[2]*d[2]*BlmBecl + Beclimb);
        lXB[1][0] = lXB[0][1];
        lXB[2][0] = lXB[0][2];
        lXB[2][1] = lXB[1][2];
    }else{
        Vec3TransposeAndMult(d, lineTensor);
        Vec3TransposeAndMult(n, climbTensor);
        Vec3TransposeAndMult(m, glideTensor);
        for(j=0;j<3;j++){
            for(k=0;k<3;k++){
                lXB[j][k] = oneSixthLen*(Bglide*glideTensor[j][k] + Beclimb*climbTensor[j][k] + Bline*lineTensor[j][k]);
            }
        }
    }

    return;
}
void CaculateTensorLengthXBij_BCC(Home_t *home, Node_t *nodei, Node_t *nodej, real8 (*lXB)[3])
{
    Param_t         *param = home->param;
    int             armID;
    real8           d[3], n[3], m[3], b[3], mag2, length, oneSixthLen;
    real8           costheta, costheta2, sintheta2;
    real8           Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
    real8           Bscrew2, Beclimb2;
    real8           invBscrew2, invBedge2;
    real8           BlmBsc, BclmBsc, BglmBsc, BlmBecl;
    bool            tdKink;

    Bscrew = 1.0 / param->MobScrew;
    Bedge  = 1.0 / param->MobEdge;

/*
 *  Climb is set this way to help with convergence and the fact that this
 *  is a glide restricted mobility
 */
    Beclimb  = 10000.0 * Bedge;
    Beclimb  = 0.0;
    Bscrew2  = Bscrew * Bscrew;
    Beclimb2 = Beclimb * Beclimb;
    
    Bline      = 1.0 * MIN(Bscrew, Bedge);
    BlmBsc     = Bline - Bscrew;
    BlmBecl    = Bline - Beclimb; 
    
    invBscrew2 = 1.0 / (Bscrew*Bscrew);
    invBedge2  = 1.0 / (Bedge*Bedge);
    
    d[0] = nodej->x - nodei->x;
    d[1] = nodej->y - nodei->y;
    d[2] = nodej->z - nodei->z;
    ZImage(param,d,d+1,d+2);
    mag2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    if(mag2<1E-12){
        Matrix33Init(lXB);
        return;
    }
    length = sqrt(mag2);
    d[0] /= length; d[1] /= length; d[2] /= length;
    oneSixthLen = length/6.0e0;

    armID = GetArmID(nodei,nodej);

    n[0] = nodei->nx[armID];
    n[1] = nodei->ny[armID];
    n[2] = nodei->nz[armID];

    b[0] = nodei->burgX[armID];
    b[1] = nodei->burgY[armID];
    b[2] = nodei->burgZ[armID];

/*          
 *  Calculate how close to screw the arm is
 */ 
    costheta = (d[0]*b[0] + d[1]*b[1] + d[2]*b[2]);
    costheta2 = (costheta*costheta);
    sintheta2 = 1.0 - costheta2;
    if(sintheta2<0)sintheta2=0.0;

    cross(n, d, m);
    NormalizeVec(m);

    if(dot_product(m[0],m[1],m[2],nodei->fX+nodej->fX,nodei->fY+nodej->fY,nodei->fZ+nodej->fZ)<0){
        m[0]=-m[0]; m[1]=-m[1]; m[2]=-m[2];        
    }
    if(dot_product(d[0],d[1],d[2],nodei->fX+nodej->fX,nodei->fY+nodej->fY,nodei->fZ+nodej->fZ)<0){
        d[0]=-d[0]; d[1]=-d[1]; d[2]=-d[2];        
    }

    tdKink = (nodei->onGB > -1 && nodej->onGB > -1 && fabs(DotProduct(b, param->TBNormal))>5E-2);

    if((DotProduct(m,m) < 1.0E-2) || ((GetPlaneType_FCC(home, nodei, nodej)& GLIDE_PLANE)==0) || tdKink==true){
        lXB[0][0] = oneSixthLen * (d[0]*d[0]*BlmBecl + Beclimb);
        lXB[0][1] = oneSixthLen * (d[0]*d[1]*BlmBecl);
        lXB[0][2] = oneSixthLen * (d[0]*d[2]*BlmBecl);
        lXB[1][1] = oneSixthLen * (d[1]*d[1]*BlmBecl + Beclimb);
        lXB[1][2] = oneSixthLen * (d[1]*d[2]*BlmBecl);
        lXB[2][2] = oneSixthLen * (d[2]*d[2]*BlmBecl + Beclimb);
    }else{
/*
 *      Arm is a regular glide arm, so build the drag matrix
 *      assuming the dislocation is screw type
 */
        lXB[0][0] = oneSixthLen * (d[0]*d[0] * BlmBsc + Bscrew);
        lXB[0][1] = oneSixthLen * (d[0]*d[1] * BlmBsc);
        lXB[0][2] = oneSixthLen * (d[0]*d[2] * BlmBsc);
        lXB[1][1] = oneSixthLen * (d[1]*d[1] * BlmBsc + Bscrew);
        lXB[1][2] = oneSixthLen * (d[1]*d[2] * BlmBsc);
        lXB[2][2] = oneSixthLen * (d[2]*d[2] * BlmBsc + Bscrew);

/*
 *      Now correct the drag matrix for dislocations that are
 *      not screw
 */
        if ((1.0 - costheta2) > 1.01-12) {

            Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) * costheta2);
            Bglide = 1.0 / Bglide;
            Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) * costheta2);

            BclmBsc = Bclimb - Bscrew;
            BglmBsc = Bglide - Bscrew;
        
            lXB[0][0] += oneSixthLen * (n[0]*n[0] * BclmBsc + m[0]*m[0] * BglmBsc);
            lXB[0][1] += oneSixthLen * (n[0]*n[1] * BclmBsc + m[0]*m[1] * BglmBsc);
            lXB[0][2] += oneSixthLen * (n[0]*n[2] * BclmBsc + m[0]*m[2] * BglmBsc);
            lXB[1][1] += oneSixthLen * (n[1]*n[1] * BclmBsc + m[1]*m[1] * BglmBsc);
            lXB[1][2] += oneSixthLen * (n[1]*n[2] * BclmBsc + m[1]*m[2] * BglmBsc);
            lXB[2][2] += oneSixthLen * (n[2]*n[2] * BclmBsc + m[2]*m[2] * BglmBsc);
        }
    }

    lXB[1][0] = lXB[0][1];
    lXB[2][0] = lXB[0][2];
    lXB[2][1] = lXB[1][2];
    return;
}

void CaculateTensorLengthXBij(Home_t *home, Node_t *nodei, Node_t *nodej, real8 (*lXB)[3])
{
    Param_t *param =  home->param;
    if(param->materialType == MAT_TYPE_FCC){
        CaculateTensorLengthXBij_FCC(home, nodei, nodej, lXB);
    }else if(param->materialType == MAT_TYPE_BCC){
        CaculateTensorLengthXBij_BCC(home, nodei, nodej, lXB);
    }else{
        Fatal("%s: Not support material type", __func__);
    }
}

PetscErrorCode FormVKFunction(SNES snes,Vec velocities,Vec F,void *ptr)
{
#ifdef GHOST
    Vec                 lv;
    Vec                 locVel, locF;
    PetscInt            totNodes;
    const PetscScalar   *lvArray;
    const PetscScalar   *locFArray;
    PetscScalar         *locVelArray;
    PetscScalar         *lfArray;
#endif
    Vec                 lf;
    const int           bs = 3;
    Home_t              *home = (Home_t *)ptr;
    PetscVKCtx_t        *vk = home->vk; 

#ifndef GHOST 
    VecScatter          scatter = vk->scatter;
#endif

#ifdef PETSC_DEBUG
    static int  iter=0;
#endif
    PetscFunctionBegin;
#ifdef GHOST

    totNodes = vk->nNatives+vk->nGhosts;

    PetscCall(VecGhostGetLocalForm(velocities,&lv));
    PetscCall(VecGetArrayRead(lv,&lvArray));
    PetscCall(VecGetArray(locVel,&locVelArray));
    for(int i=0;i<totNodes;i++){
        locVelArray[bs*i]   = lvArray[bs*i];
        locVelArray[bs*i+1] = lvArray[bs*i+1];
        locVelArray[bs*i+2] = lvArray[bs*i+2];
    }
    PetscCall(VecRestoreArray(locVel,&locVelArray));
    PetscCall(VecRestoreArrayRead(lv,&lvArray));
    PetscCall(VecGhostRestoreLocalForm(velocities,&lv));

    PetscCall(VecCreateSeq(PETSC_COMM_SELF,bs*vk->nNatives,&locF));
    PetscCall(MatMult(vk->velCoeMat,locVel,locF));

    PetscCall(VecGhostGetLocalForm(F, &lf));
    PetscCall(VecGetArray(lf,&lfArray));
    PetscCall(VecGetArrayRead(locF,&locFArray));
    for(int i=0;i<totNodes;i++){
        lfArray[bs*i]   = locFArray[bs*i];
        lfArray[bs*i+1] = locFArray[bs*i+2];
        lfArray[bs*i+2] = locFArray[bs*i+2];
    }
    PetscCall(VecRestoreArrayRead(locF,&locFArray));
    PetscCall(VecRestoreArray(lf,&lfArray));
    PetscCall(VecGhostRestoreLocalForm(F, &lf));

    PetscCall(VecDestroy(&locVel));
    PetscCall(VecDestroy(&locF));
    
    PetscCall(VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD));
#else
    PetscCall(VecScatterBegin(scatter, velocities, vk->locVel, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, velocities, vk->locVel, INSERT_VALUES, SCATTER_FORWARD));

    PetscCall(VecCreateSeq(MPI_COMM_SELF,vk->nNatives*bs,&lf));
    PetscCall(VecSetBlockSize(lf,bs));

    PetscCall(VecGetLocalVector(F,lf));
    PetscCall(MatMult(vk->velCoeMat,vk->locVel,lf));
    PetscCall(VecRestoreLocalVector(F,lf));

    PetscCall(VecDestroy(&lf));
#endif
#ifdef PETSC_DEBUG
    iter++;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"VEC F cycle=%d iter=%d:\n",home->cycle, iter));
    PetscCall(VecView(F,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"forces:\n"));
    PetscCall(VecView(vk->forces,PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Mat coe:\n"));
    PetscCall(MatView(vk->velCoeMat,PETSC_VIEWER_STDOUT_SELF));
#ifdef GHOST
    PetscCall(VecView(vk->velocities,PETSC_VIEWER_STDOUT_WORLD));
#else
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"local velocity:\n"));
    PetscCall(VecView(vk->locVel,PETSC_VIEWER_STDOUT_SELF));
#endif
#endif
    PetscFunctionReturn(0);

}

PetscErrorCode FormVKJacobianNolinear(SNES snes,Vec X,Mat J,Mat jac,void *ptr)
{
    PetscInt        i,j,rows[3],cols[6],rowID,colID;
    real8           lXB[3][3], values[18];
    Node_t          *node, *nbr;
    Home_t          *home       = (Home_t *)ptr;
    PetscVKCtx_t    *vk         = home->vk;
    PetscInt        nNatives    = vk->nNatives;
    Node_t          **nodeList  = vk->nodeList;
    const PetscInt  bs          = 3;

    PetscFunctionBegin;
    
    PetscCall(MatZeroEntries(jac));

    for(i=0;i<nNatives;i++){
        node = nodeList[i];
        rowID = node->gloID*bs;
        for(j=0;j<node->numNbrs;j++){
            nbr = GetNeighborNode(home,node,j);
            if(nbr->gloID < node->gloID)continue;
            colID = nbr->gloID*bs;
            CaculateTensorLengthXBij(home,node,nbr,lXB);

            values[0] = 2.0*lXB[0][0]; values[1] = 2.0*lXB[0][1]; values[2] = 2.0*lXB[0][2];
            values[3] = lXB[0][0];     values[4] = lXB[0][1];     values[5] = lXB[0][2];

            values[6] = 2.0*lXB[1][0]; values[7]  = 2.0*lXB[1][1]; values[8]  = 2.0*lXB[1][2];
            values[9] = lXB[1][0];     values[10] = lXB[1][1];     values[11] = lXB[1][2];
            
            values[12] = 2.0*lXB[2][0]; values[13] = 2.0*lXB[2][1]; values[14] = 2.0*lXB[2][2];
            values[15] = lXB[2][0];     values[16] = lXB[2][1];     values[17] = lXB[2][2];

            rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
            cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
            cols[3] = colID; cols[4] = colID+1; cols[5] = colID+2;
            PetscCall(MatSetValues(jac,3,rows,6,cols,values,ADD_VALUES));

            rows[0] = colID; rows[1] = colID+1; rows[2] = colID+2;
            cols[0] = colID; cols[1] = colID+1; cols[2] = colID+2;
            cols[3] = rowID; cols[4] = rowID+1; cols[5] = rowID+2;

            PetscCall(MatSetValues(jac,3,rows,6,cols,values,ADD_VALUES));
        }
        rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
        cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
        NodalConstraintTensor(home, vk, i, values);
        PetscCall(MatSetValues(jac,3,rows,3,cols,values,ADD_VALUES));
    }
    PetscCall(MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY));

    /*
        Tell the matrix we will never add a new nonzero location to the
        matrix. If we do, it will generate an error.
     */
    PetscCall(MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE));

    /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
    PetscCall(MatSetOption(jac,MAT_SYMMETRIC,PETSC_TRUE));

#ifdef PETSC_DEBUG
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"View Jacobi\n"));
    PetscCall(MatView(jac,PETSC_VIEWER_STDOUT_WORLD));
#endif

    PetscFunctionReturn(0);
}

/**
 * > FormVKJacobianLinear() is a function that forms the Jacobian matrix 
 * 
 * @param home the home structure
 * @param vk the PetscVKCtx_t structure
 */
PetscErrorCode FormVKJacobianLinear(Home_t *home, PetscVKCtx_t *vk)
{
    int             i,j,rows[3],cols[6],rowID,colID,nbrID;
    real8           lXB[3][3], values[18];
    Node_t          *node, *nbr;
    PetscInt        nNatives    = vk->nNatives;
    Node_t          **nodeList  = vk->nodeList;
    const int       bs          = 3;

    PetscFunctionBegin;

    PetscCall(MatZeroEntries(vk->velJacobi));
    if(vk->ltog){
        for(i=0;i<nNatives;i++){
            node = nodeList[i];
            rowID = node->locID*bs;
            for(j=0;j<node->numNbrs;j++){
                nbr = GetNeighborNode(home,node,j);
                nbrID = nbr->locID; colID=bs*nbrID;
                if(nbrID<node->locID)continue;
                CaculateTensorLengthXBij(home,node,nbr,lXB);
                values[0] = 2.0*lXB[0][0]; values[1] = 2.0*lXB[0][1]; values[2] = 2.0*lXB[0][2];
                values[3] = lXB[0][0];     values[4] = lXB[0][1];     values[5] = lXB[0][2];
       
                values[6] = 2.0*lXB[1][0]; values[7]  = 2.0*lXB[1][1]; values[8]  = 2.0*lXB[1][2];
                values[9] = lXB[1][0];     values[10] = lXB[1][1];     values[11] = lXB[1][2];
                
                values[12] = 2.0*lXB[2][0]; values[13] = 2.0*lXB[2][1]; values[14] = 2.0*lXB[2][2];
                values[15] = lXB[2][0];     values[16] = lXB[2][1];     values[17] = lXB[2][2];
       
                rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
                cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
                cols[3] = colID; cols[4] = colID+1; cols[5] = colID+2;
                PetscCall(MatSetValuesLocal(vk->velJacobi,3,rows,6,cols,values,ADD_VALUES));

                if(nbr->locID<nNatives){
                    rows[0] = colID; rows[1] = colID+1; rows[2] = colID+2;
                    cols[0] = colID; cols[1] = colID+1; cols[2] = colID+2;
                    cols[3] = rowID; cols[4] = rowID+1; cols[5] = rowID+2;
                    PetscCall(MatSetValuesLocal(vk->velJacobi,3,rows,6,cols,values,ADD_VALUES));
                }
            }
            rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
            cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
            NodalConstraintTensor(home, vk, i, values);
            PetscCall(MatSetValuesLocal(vk->velJacobi,3,rows,3,cols,values,ADD_VALUES));
        }
    }else{
        for(i=0;i<nNatives;i++){
            node = nodeList[i];
            rowID = node->gloID*bs;
            for(j=0;j<node->numNbrs;j++){
                nbr = GetNeighborNode(home,node,j);
                nbrID = nbr->gloID; colID=bs*nbrID;
                if(nbrID<node->gloID)continue;
                CaculateTensorLengthXBij(home,node,nbr,lXB);
                values[0] = 2.0*lXB[0][0]; values[1] = 2.0*lXB[0][1]; values[2] = 2.0*lXB[0][2];
                values[3] = lXB[0][0];     values[4] = lXB[0][1];     values[5] = lXB[0][2];
       
                values[6] = 2.0*lXB[1][0]; values[7]  = 2.0*lXB[1][1]; values[8]  = 2.0*lXB[1][2];
                values[9] = lXB[1][0];     values[10] = lXB[1][1];     values[11] = lXB[1][2];
                
                values[12] = 2.0*lXB[2][0]; values[13] = 2.0*lXB[2][1]; values[14] = 2.0*lXB[2][2];
                values[15] = lXB[2][0];     values[16] = lXB[2][1];     values[17] = lXB[2][2];
       
                rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
                cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
                cols[3] = colID; cols[4] = colID+1; cols[5] = colID+2;
                PetscCall(MatSetValues(vk->velJacobi,3,rows,6,cols,values,ADD_VALUES));

                rows[0] = colID; rows[1] = colID+1; rows[2] = colID+2;
                cols[0] = colID; cols[1] = colID+1; cols[2] = colID+2;
                cols[3] = rowID; cols[4] = rowID+1; cols[5] = rowID+2;
                PetscCall(MatSetValues(vk->velJacobi,3,rows,6,cols,values,ADD_VALUES));
            }
            rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
            cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
            NodalConstraintTensor(home, vk, i, values);
            PetscCall(MatSetValues(vk->velJacobi,3,rows,3,cols,values,ADD_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(vk->velJacobi,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(vk->velJacobi,MAT_FINAL_ASSEMBLY));

     /*
       Tell the matrix we will never add a new nonzero location to the
       matrix. If we do, it will generate an error.
    */
    PetscCall(MatSetOption(vk->velJacobi,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE));

    /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
    PetscCall(MatSetOption(vk->velJacobi,MAT_SYMMETRIC,PETSC_TRUE));

#ifdef PETSC_DEBUG
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"View Jacobi\n"));
    PetscCall(MatView(vk->velJacobi,PETSC_VIEWER_STDOUT_WORLD));
#endif
    PetscFunctionReturn(0);
}

/**
 * > This function is called by the PETSc library to get the local kernel function
 * velocity calculation
 * 
 * @param home The home structure.
 * @param vk The context for the velocity kernel.
 */
PetscErrorCode InitilizeLocalVelocityKernel(Home_t *home, PetscVKCtx_t *vk)
{
    PetscFunctionBegin;
    vk->nNatives = 0;
    vk->nolinear = home->vk->nolinear;
    PetscFunctionReturn(0);
}

PetscErrorCode AddNodeToLocalVelocityKernel(Home_t *home, PetscVKCtx_t *vk, Node_t *node)
{
    int id,i;

    PetscFunctionBegin;
    for(i=0;i<vk->nNatives;i++){
        if(node == vk->nodeList[i]){
            PetscFunctionReturn(0);
        }
    }
    
    id = vk->nNatives;
    vk->nNatives++;
    if(id==0){
        PetscCall(PetscMalloc1(1, &vk->nodeList)); 
        vk->nodeList[id]=node; node->locID=id;
        PetscCall(PetscMalloc1(1, &vk->nCstList));
        PetscCall(PetscMalloc1(3, &vk->cstList));
    }else{
        PetscCall(PetscRealloc(vk->nNatives*sizeof(Node_t *),&vk->nodeList));
        vk->nodeList[id] = node; node->locID = id;
        PetscCall(PetscRealloc(vk->nNatives*sizeof(PetscInt), &vk->nCstList));
        PetscCall(PetscRealloc(3*vk->nNatives*sizeof(real8),&vk->cstList));

    }


    PetscFunctionReturn(0);
}

PetscErrorCode SolveLocalVelocityKernel(Home_t *home, PetscVKCtx_t *vk)
{
        Node_t      *node, *nbr;
        PetscInt    *nnz, i, j, k, nNatives;
        PetscInt    rows[3], cols[6], rowID, colID;
        real8       lXB[3][3], values[18], nbrVel[3], dragForce[3];
        PetscReal   **forceMat, **vArray, vel, vshr;
        KSP         ksp;
        PC          pc;
        KSPConvergedReason      reason;

        PetscFunctionBegin;
        
        vshr = home->param->vshr;
        nNatives = vk->nNatives;
        PetscCall(PetscMalloc1(3*nNatives, &nnz));

        PetscCall(VecCreateSeq(PETSC_COMM_SELF,3*nNatives,&vk->velocities));
        PetscCall(VecDuplicate(vk->velocities,&vk->forces));
        PetscCall(UpdateVelocityForce(home,vk));

        for(i=0;i<nNatives;i++){
            node = vk->nodeList[i];
            nnz[3*i] = 3; nnz[3*i+1] = 3; nnz[3*i+2] = 3;
            for(j=0;j<node->numNbrs;j++){
                nbr = GetNeighborNode(home, node, j);
                if(nbr == (Node_t *)NULL)continue;
                for(k=0; k<nNatives;k++){
                    if(nbr == vk->nodeList[k])break;
                }
                if(k!=nNatives){
                    nnz[3*i] +=3; nnz[3*i+1] += 3; nnz[3*i+2] += 3;
                }else{
                    nbr->locID = -1;
                }
            }
        }
        
        PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF,3*nNatives, 3*nNatives, 0, nnz, &vk->velJacobi));
        PetscCall(MatZeroEntries(vk->velJacobi));

        PetscCall(VecGetArray2dWrite(vk->forces,nNatives,3,0,0,&forceMat));
        
        for(i=0;i<nNatives;i++){
            node = vk->nodeList[i];
            rowID = i*3;
            for(j=0;j<node->numNbrs;j++){
                nbr = GetNeighborNode(home,node,j); 
                if(nbr->locID >= 0  &&  nbr->locID < node->locID)continue;
                colID = 3*nbr->locID;

                CaculateTensorLengthXBij(home, node, nbr, lXB);

                if(colID >=0 ){
                    values[0] = 2.0*lXB[0][0]; values[1] = 2.0*lXB[0][1]; values[2] = 2.0*lXB[0][2];
                    values[3] = lXB[0][0];     values[4] = lXB[0][1];     values[5] = lXB[0][2];
                   
                    values[6] = 2.0*lXB[1][0]; values[7]  = 2.0*lXB[1][1]; values[8]  = 2.0*lXB[1][2];
                    values[9] = lXB[1][0];     values[10] = lXB[1][1];     values[11] = lXB[1][2];
                    
                    values[12] = 2.0*lXB[2][0]; values[13] = 2.0*lXB[2][1]; values[14] = 2.0*lXB[2][2];
                    values[15] = lXB[2][0];     values[16] = lXB[2][1];     values[17] = lXB[2][2];
                   
                    rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
                    cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
                    cols[3] = colID; cols[4] = colID+1; cols[5] = colID+2;
                    PetscCall(MatSetValues(vk->velJacobi,3,rows,6,cols,values,ADD_VALUES));

                    rows[0] = colID; rows[1] = colID+1; rows[2] = colID+2;
                    cols[0] = colID; cols[1] = colID+1; cols[2] = colID+2;
                    cols[3] = rowID; cols[4] = rowID+1; cols[5] = rowID+2;
                    PetscCall(MatSetValues(vk->velJacobi,3,rows,6,cols,values,ADD_VALUES));
                }else{
                    values[0] = 2.0*lXB[0][0]; values[1] = 2.0*lXB[0][1]; values[2] = 2.0*lXB[0][2];
                    values[3] = 2.0*lXB[1][0]; values[4] = 2.0*lXB[1][1]; values[5] = 2.0*lXB[1][2];
                    values[6] = 2.0*lXB[2][0]; values[7] = 2.0*lXB[2][1]; values[8] = 2.0*lXB[2][2];
                    rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
                    cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
                    PetscCall(MatSetValues(vk->velJacobi,3,rows,3,cols,values,ADD_VALUES));

                    VectorCopy(nbrVel, nbr->vX, nbr->vY, nbr->vZ);
                    Matrix33Vector3Multiply(lXB,nbrVel,dragForce);
                    forceMat[i][0] -= dragForce[0];
                    forceMat[i][1] -= dragForce[1];
                    forceMat[i][2] -= dragForce[2];
                }
            }
            rows[0] = rowID; rows[1] = rowID+1; rows[2] = rowID+2;
            cols[0] = rowID; cols[1] = rowID+1; cols[2] = rowID+2;
            NodalConstraintTensor(home, vk, i, values);
            PetscCall(MatSetValues(vk->velJacobi,3,rows,3,cols,values,ADD_VALUES));
        }

        PetscCall(VecRestoreArray2dWrite(vk->forces,nNatives,3,0,0,&forceMat));

        PetscCall(MatAssemblyBegin(vk->velJacobi,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(vk->velJacobi,MAT_FINAL_ASSEMBLY));
       
         /*
           Tell the matrix we will never add a new nonzero location to the
           matrix. If we do, it will generate an error.
        */
        PetscCall(MatSetOption(vk->velJacobi,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE));
       
        /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
        PetscCall(MatSetOption(vk->velJacobi,MAT_SYMMETRIC,PETSC_TRUE));

        PetscCall(KSPCreate(PETSC_COMM_SELF,&ksp));
        PetscCall(KSPSetType(ksp,KSPGMRES));
        PetscCall(KSPSetOperators(ksp,vk->velJacobi,vk->velJacobi));
        
        /*  
        Set linear solver defaults for this problem (optional).
        - By extracting the KSP and PC contexts from the KSP context,
            we can then directly call any KSP and PC routines to set
            various options.
        - The following four statements are optional; all of these
            parameters could alternatively be specified at runtime via
            KSPSetFromOptions();
        */
        PetscCall(KSPGetPC(ksp,&pc));
        PetscCall(PCSetType(pc,PCJACOBI));
        PetscCall(KSPSetTolerances(ksp,1.0E-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT));

        PetscCall(KSPSetFromOptions(ksp));
        PetscCall(KSPSolve(ksp,vk->forces,vk->velocities));    
        
        PetscCall(KSPGetConvergedReason(ksp,&reason));
#ifdef PETSC_DEBUG
        PetscCheck(reason > 0,PetscObjectComm((PetscObject)ksp),PETSC_ERR_PLIB,"KSPConvergedReason() = %s [%d]",KSPConvergedReasons[reason],reason);
#else
        if(reason<0){
            PetscPrintf(PETSC_COMM_SELF,"KSPConvergedReason() = %s [%d]\n",KSPConvergedReasons[reason],reason);
        }
#endif

        PetscCall(VecGetArray2d(vk->velocities,nNatives,3,0,0,&vArray));
        for(i=0;i<nNatives;i++){
            node = vk->nodeList[i];

            if(vk->nCstList[i]==1){
                NormalComponent(vArray[i],&vk->cstList[3*i],vArray[i]);
            }else if(vk->nCstList[i]==2){
                ParallelComponent(vArray[i],&vk->cstList[3*i],vArray[i]);
            }else if(vk->nCstList[i]==3){
                VECTOR_ZERO(vArray[i]);
            }

            vel = Normal(vArray[i]);
            if(vel > vshr){
                vel = vshr/vel; 
            }else{
                vel = 1.0;
            }
            node->vX = vArray[i][0]*vel;
            node->vY = vArray[i][1]*vel;
            node->vZ = vArray[i][2]*vel;
        }
        PetscCall(VecRestoreArray2d(vk->velocities,nNatives,3,0,0,&vArray));

        PetscCall(PetscFree(nnz));
        PetscCall(KSPDestroy(&ksp));
        PetscCall(PetscFree(vk->nodeList));
        PetscCall(PetscFree(vk->nCstList));
        PetscCall(PetscFree(vk->cstList));
        PetscCall(VecDestroy(&vk->velocities));
        PetscCall(VecDestroy(&vk->forces));
        PetscCall(MatDestroy(&vk->velJacobi));

        PetscFunctionReturn(0);
}

#endif