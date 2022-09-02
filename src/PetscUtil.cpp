
#ifdef PETSC
#include <petsc.h>
#include "PetscUtil.h"
#include "Home.h"

#define JACOBI

PetscErrorCode FormFunction(SNES, Vec, Vec, void*);
PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
/**
 * `PetscInitialize` initializes the PETSc library and reads the options file
 * 
 * @param home This is the home structure that contains all the information about the simulation.
 * 
 * @return The petsc error code.
 */
PetscErrorCode PetscInitialize(Home_t *home){
    Param_t *param = home->param;

    PetscFunctionBegin;

    PetscCall(PetscInitializeNoArguments()); 

    if(strcmp(param->petscOptionsFile, "0") == 0)PetscFunctionReturn(0);

/* Reading the options file and printing the options to the screen. */
    PetscCall(PetscOptionsInsertFileYAML(PETSC_COMM_WORLD,NULL,
                    param->petscOptionsFile,PETSC_TRUE));
    PetscCall(PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(PetscTest(home));

    PetscFunctionReturn(0);
}


typedef struct{
    PetscInt    nlocal, nghost, myDomain, numDomains, nvertices;
    PetscInt    *vertices;
    Vec         localx, localF;
    VecScatter  scatter;
    FILE        *fp;
}User_t;

/**
 * A function that is used to test the PETSc library.
 * This is a small demonstration for building a system of partial differential equations (PDEs) 
 * of dislocation velocity solutions.
 * 
 * @param home the home object
 */
PetscErrorCode PetscTest(Home_t *home) {

    PetscMPIInt         rank, size;
    PetscInt            n = 4, nvertices, bs=2, initID[2];
    PetscInt            nlocl, nghost = 0, i;
    SNES                snes;          // nonlinear solver
    Vec                 x, r, b;          // solution, residual and right-hand vectors
    PetscScalar const   *array;
    PetscScalar         initValue[2];
    IS                  is, isglobal;
    User_t              user;
    FILE                *fp;
    PetscInt            bid[4]={0,1,2,3};
    PetscScalar         bList[4]={21.0515,3.07944,-3.6212,8.06307};

#ifdef JACOBI
    Mat                     Jac;
    ISLocalToGlobalMapping  isl2g;
    MatFDColoring           matfdcoloring = 0;
    PetscBool               fd_jacobian_coloring = PETSC_FALSE;
#endif
    char            partName[256];

    PetscFunctionBegin;

    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Petsc Test\n"));
    sprintf(partName,"output.%d",rank);
    user.fp = fopen(partName,"w");
    fp      = user.fp;

    if(size>1){
        if(rank<2){
            nlocl = 1; nghost=1;
            nvertices = nlocl + nghost;
            PetscCall(PetscMalloc1(nvertices, &user.vertices));
            user.vertices[0] = rank;
            if(rank==0)user.vertices[1]=1;else user.vertices[1]=0;
        }else{
            nlocl = 0; nghost=(rank%3);
            nvertices = nlocl + nghost;
            if(nghost){
                PetscCall(PetscMalloc1(nghost, &user.vertices));
                if(nghost==1)user.vertices[0] = rank%2; else {user.vertices[0]=0; user.vertices[1]=1;}
            }
        }
    }else{
        nlocl = 2; nghost = 0; nvertices = 2;
        PetscCall(PetscMalloc1(nvertices,&user.vertices));
        user.vertices[0] = 0; user.vertices[1]=1;
    }
    user.nlocal = nlocl; user.nghost = nghost; user.nvertices = nvertices;
    user.myDomain = rank; user.numDomains = size;

    PetscCall(PetscSynchronizedPrintf(MPI_COMM_WORLD,"nlocal=%d, nghost=%d: ", nlocl, nghost));
    for(i=0;i<nlocl+nghost;i++) {PetscCall(PetscSynchronizedPrintf(MPI_COMM_WORLD,"%d ",user.vertices[i]));}
    PetscCall(PetscSynchronizedPrintf(MPI_COMM_WORLD,"\n"));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));

    PetscCall(VecCreateMPI(MPI_COMM_WORLD, nlocl*bs, n, &x));
    PetscCall(VecSetBlockSize(x, bs));
    PetscCall(VecCreateSeq(MPI_COMM_SELF, bs*nvertices, &user.localx));
    PetscCall(VecSetBlockSize(user.localx, bs));
    PetscCall(VecDuplicate(user.localx,&user.localF));

    PetscCall(ISCreateStride(MPI_COMM_SELF,bs*nvertices,0,1,&is));
    PetscCall(ISCreateBlock(MPI_COMM_SELF,bs,nvertices,user.vertices,PETSC_COPY_VALUES,&isglobal));
    PetscCall(PetscFPrintf(MPI_COMM_SELF,fp,"is and global is:\n"));
    PetscCall(ISView(isglobal,PETSC_VIEWER_STDOUT_SELF));

    PetscCall(VecScatterCreate(x,isglobal,user.localx,is,&user.scatter));
    PetscCall(VecScatterView(user.scatter,PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(ISDestroy(&is));
    PetscCall(ISDestroy(&isglobal));

#ifdef JACOBI
  /*
     Create matrix data structure; Just to keep the example simple, we have not done any
     preallocation of memory for the matrix. In real application code with big matrices,
     preallocation should always be done to expedite the matrix creation.
  */
    PetscCall(MatCreate(MPI_COMM_WORLD,&Jac));
    PetscCall(MatSetSizes(Jac,PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(Jac));
    PetscCall(MatSetBlockSize(Jac,bs));
    PetscCall(MatSetUp(Jac));

  /*
    The following routine allows us to set the matrix values in local ordering
  */
    PetscCall(ISLocalToGlobalMappingCreate(MPI_COMM_SELF,bs,nvertices,user.vertices,PETSC_COPY_VALUES,&isl2g));
    PetscCall(MatSetLocalToGlobalMapping(Jac,isl2g,isl2g));

#endif
    
    PetscCall(PetscFPrintf(MPI_COMM_SELF, fp, "-1 Vec Set\n"));
    if(size == 1){
        initID[0]=0; initID[1]=1; initValue[0]=1.0; initValue[1]=2.0;
        PetscCall(VecSetValues(x,2,initID,initValue,INSERT_VALUES));
        initID[1]=2; initID[1]=3; initValue[0]=3.0; initValue[1]=4.0;
        PetscCall(VecSetValues(x,2,initID,initValue,INSERT_VALUES));
    }else{
        if(rank<2){
            if(rank==0){
                initID[0]=0; initID[1]=1; initValue[0]=1.0; initValue[1]=2.0;
            }else if(rank==1){
                initID[0]=2; initID[1]=3; initValue[0]=3.0; initValue[1]=4.0;
            }
            PetscCall(VecSetValues(x,2,initID,initValue,INSERT_VALUES));
        }
    }
    PetscCall(VecAssemblyBegin(x));
    PetscCall(VecAssemblyEnd(x));
    PetscCall(VecSet(x,1.0));         // initial iterate
    PetscCall(PetscPrintf(MPI_COMM_WORLD, "Initial x\n"));
    PetscCall(VecView(x, PETSC_VIEWER_STDOUT_WORLD));

    /*
       Print out each vector, including the ghost padding region.
    */
    PetscCall(VecGetArrayRead(user.localx,&array));
    for (i=0; i<bs*(nlocl+nghost); i++) {
      PetscCall(PetscFPrintf(PETSC_COMM_SELF,fp,"[%d](%s) %D %g\n",rank, i<bs*nlocl?"native":"ghost",
                i,(double)PetscRealPart(array[i])));
    }
    PetscCall(VecRestoreArrayRead(user.localx,&array));

    PetscCall(PetscFPrintf(MPI_COMM_SELF, fp, "VecView\n"));
    PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(VecDuplicate(x,&r));
    PetscCall(VecDuplicate(x,&b));

    PetscCall(VecSetValues(b,4,bid,bList,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(b));
    PetscCall(VecAssemblyEnd(b));
    PetscCall(PetscPrintf(MPI_COMM_WORLD,"VecView b\n"));
    PetscCall(VecView(b,PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));
    PetscCall(SNESSetType(snes,SNESNEWTONLS));
    PetscCall(SNESSetFunction(snes,r,FormFunction,(void*)&user));

#ifdef JACOBI
    PetscCall(PetscOptionsGetBool(NULL,NULL,"-fd_jacobian_coloring",&fd_jacobian_coloring,0));
    if (!fd_jacobian_coloring) {
      PetscCall(SNESSetJacobian(snes,Jac,Jac,FormJacobian,(void*)&user));
    } else {  /* Use matfdcoloring */
      ISColoring   iscoloring;
      MatColoring  mc;
  
      /* Get the data structure of Jac */
      PetscCall(FormJacobian(snes,x,Jac,Jac,&user));
      /* Create coloring context */
      PetscCall(MatColoringCreate(Jac,&mc));
      PetscCall(MatColoringSetType(mc,MATCOLORINGSL));
      PetscCall(MatColoringSetFromOptions(mc));
      PetscCall(MatColoringApply(mc,&iscoloring));
      PetscCall(MatColoringDestroy(&mc));
      PetscCall(MatFDColoringCreate(Jac,iscoloring,&matfdcoloring));
      PetscCall(MatFDColoringSetFunction(matfdcoloring,(PetscErrorCode (*)(void))FormFunction,&user));
      PetscCall(MatFDColoringSetFromOptions(matfdcoloring));
      PetscCall(MatFDColoringSetUp(Jac,iscoloring,matfdcoloring));
      /* PetscCall(MatFDColoringView(matfdcoloring,PETSC_VIEWER_STDOUT_WORLD)); */
      PetscCall(SNESSetJacobian(snes,Jac,Jac,SNESComputeJacobianDefaultColor,matfdcoloring));
      PetscCall(ISColoringDestroy(&iscoloring));
    }
#endif
    PetscCall(SNESSetFromOptions(snes));
    PetscCall(SNESSolve(snes,b,x));
    PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));

    if(nvertices) {PetscCall(PetscFree(user.vertices));}
    PetscCall(VecDestroy(&user.localx));
    PetscCall(VecDestroy(&user.localF));
    PetscCall(VecScatterDestroy(&user.scatter));

    PetscCall(SNESDestroy(&snes));  
    PetscCall(VecDestroy(&x));  
    PetscCall(VecDestroy(&r));
    PetscCall(VecDestroy(&b));
    fclose(fp);
#ifdef JACOBI
    PetscCall(MatDestroy(&Jac));
#endif
    PetscFunctionReturn(0);
}

PetscErrorCode FormFunction(SNES snes, Vec x, Vec F, void *ctx) {

    const PetscReal     *ax;
    PetscReal           *aF;
    PetscInt            i;


/* A pointer to the user structure. */
    User_t     *user     = (User_t *)ctx;
    Vec        lx        = user->localx;
    VecScatter scatter   = user->scatter;
    PetscInt   nlocal    = user->nlocal;
    FILE       *fp       = user->fp;
    PetscInt   myDomain  = user->myDomain;
    PetscInt   numDomains  = user->numDomains;
    
    PetscFunctionBegin;
    PetscCall(PetscFPrintf(MPI_COMM_SELF,fp,"%s: 1, nlocal=%d\n",__func__,nlocal));

    PetscCall(VecScatterBegin(scatter,x,lx,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter,x,lx,INSERT_VALUES,SCATTER_FORWARD));

    /*
     * Get pointer to vector data
     */
    PetscCall(VecGetArrayRead(lx,&ax));
    PetscCall(VecGetArray(F,&aF));

    for(i=0;i<nlocal;i++){
        if(numDomains>1){
            if(!myDomain){
                aF[0] = ax[0]*PetscCosReal(M_PI/ax[1]/6.0) + ax[0]*PetscExpReal(ax[2]);
                aF[1] = ax[1]*PetscSinReal(M_PI/ax[0]/6.0) + ax[2]*PetscLogReal(ax[1]) ;
            }else{
                aF[0] = ax[0]*PetscLogReal(M_PI/ax[3]/6.0) + ax[1]*PetscSinReal(ax[2]/10);
                aF[1] = ax[1]*PetscExpReal(M_PI/ax[3]/6.0) + ax[0]*PetscCosReal(ax[0]/10); 
            }
            PetscPrintf(MPI_COMM_SELF,"[%d]: aF %g %g\n",myDomain, aF[0],aF[1]);
            PetscFPrintf(MPI_COMM_SELF,fp,"[%d]: aF %g %g\n",myDomain, aF[0],aF[1]);
        }else{
            aF[0] = ax[0]*PetscCosReal(M_PI/ax[1]/6.0) + ax[0]*PetscExpReal(ax[2]);
            aF[1] = ax[1]*PetscSinReal(M_PI/ax[0]/6.0) + ax[2]*PetscLogReal(ax[1]);

            aF[2] = ax[2]*PetscLogReal(M_PI/ax[1]/6.0) + ax[3]*PetscSinReal(ax[0]/10);
            aF[3] = ax[3]*PetscExpReal(M_PI/ax[1]/6.0) + ax[2]*PetscCosReal(ax[2]/10);
            PetscPrintf(MPI_COMM_SELF,"[%d]: aF %g %g %g %g\n",myDomain, aF[0],aF[1], aF[2],aF[3]);
            PetscFPrintf(MPI_COMM_SELF,fp,"[%d]: aF %g %g %g %g\n",myDomain, aF[0],aF[1],aF[2],aF[3]);
        }
    }

    PetscCall(VecRestoreArrayRead(lx,&ax));
    PetscCall(VecRestoreArray(F,&aF));
    PetscCall(PetscPrintf(MPI_COMM_WORLD,"Output F\n"));
    PetscCall(VecView(F,PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(PetscFPrintf(MPI_COMM_SELF,fp,"%s 1\n",__func__));
    
    PetscFunctionReturn(0);
}
/* --------------------  Evaluate Jacobian F'(x) -------------------- */
/*
   FormJacobian - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  X - input vector
.   - optional user-defined context, as set by SNESSetJacobian()

   Output Parameters:
.  A - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure

*/
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat J,Mat jac,void *ctx)
{
    const PetscReal         *ax;
    PetscScalar             af[4][4], v[4]={0.0,0.0,0.0,0.0}, value[4];
    PetscScalar             t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    PetscInt                i, j, id[4],col[4];

/* A pointer to the user structure. */
    User_t     *user      = (User_t *)ctx;
    Vec        lx         = user->localx;
    VecScatter scatter    = user->scatter;
    PetscInt   nlocal     = user->nlocal;
    FILE       *fp        = user->fp;
    PetscInt   myDomain   = user->myDomain;
    PetscInt   numDomains = user->numDomains;
    PetscInt   nvertices  = user->nvertices;

    PetscFunctionBegin;

    PetscCall(VecScatterBegin(scatter,x,lx,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter,x,lx,INSERT_VALUES,SCATTER_FORWARD));
   
    /*
       Get pointer to vector data
    */
    PetscCall(VecGetArrayRead(lx,&ax));

    if(!numDomains){
        v[0] = ax[0]; v[1]=ax[1]; v[2]=ax[2]; v[3]=ax[3];
        id[0]= 0; id[1]=1; id[2]=2; id[3]=3;
    }else{
        if(!myDomain){
            v[0] = ax[0]; v[1]=ax[1]; v[2]=ax[2]; v[3]=ax[3];
            id[0]= 0; id[1]=1; id[2]=2; id[3]=3;
        }else if(myDomain==1){
            v[0] = ax[2]; v[1]=ax[3]; v[2]=ax[0]; v[3]=ax[1];
            id[0]= 2; id[1]=3; id[2]=0; id[3]=1;
        }
    }
   
    if(nlocal){
        t1 = 0.1e1 / v[1];
        t2 = 0.3141592654e1 / 0.6e1;
        t3 = t2 * t1;
        t4 = PetscExpReal(v[2]);
        t5 = PetscPowReal(t1, 0.2e1);
        t6 = 0.1e1 / v[0];
        t7 = t2 * t6;
        t1 = v[2] * t1;
        t8 = v[0] / 0.10e2;
        t9 = PetscExpReal(t3);
        t10 = v[2] / 0.10e2;

        af[0][0] = PetscCosReal(t3) + t4;
        af[0][1] = t2 * v[0] * t5 * PetscSinReal(t3);
        af[0][2] = v[0] * t4;
        af[0][3] = 0;
        af[1][0] = -t2 * v[1] * PetscPowReal(t6, 0.2e1) * PetscCosReal(t7);
        af[1][1] = PetscSinReal(t7) + t1;
        af[1][2] = PetscLogReal(v[1]);
        af[1][3] = 0;
        af[2][0] = v[3] * PetscCosReal(t8) / 0.10e2;
        af[2][1] = -t1;
        af[2][2] = PetscLogReal(t3);
        af[2][3] = PetscSinReal(t8);
        af[3][0] = 0;
        af[3][1] = -t2 * v[3] * t5 * t9;
        af[3][2] = -t10 * PetscSinReal(t10) + PetscCosReal(t10);
        af[3][3] = t9;

        for (i=0; i < 2*nlocal; i++) {
            for (j = 0; j < 2*nvertices; j++) {
                col[j] = j;
                value[j] = af[id[i]][id[j]];
            }
            PetscFPrintf(MPI_COMM_SELF,fp,"jab[%d][0:4]=(%d)%g (%d)%g (%d)%g (%d)%g\n",
                            i,id[0],value[id[0]],id[1],value[id[1]],id[2],value[id[2]],id[3],value[id[3]]); 
            /*
              Set the matrix values in the local ordering. Note that in order to use this
              feature we must call the routine MatSetLocalToGlobalMapping() after the
              matrix has been created.
            */
            PetscCall(MatSetValuesLocal(jac,1,&i,2*nvertices,col,value,INSERT_VALUES));
       
        }
    }
   
    /*
       Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().
       Between these two calls, the pointer to vector data has been restored to
       demonstrate the use of overlapping communicationn with computation.
    */
    PetscCall(MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY));
    PetscCall(VecRestoreArrayRead(lx,&ax));
    PetscCall(MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY));
   
    /*
       Tell the matrix we will never add a new nonzero location to the
       matrix. If we do, it will generate an error.
    */
    PetscCall(MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE));
    /* MatView(jac,PETSC_VIEWER_STDOUT_SELF); */
    PetscFunctionReturn(0);
}
#endif