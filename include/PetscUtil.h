#ifndef PetscUtil_h
#define PetscUtil_h
#include "Typedefs.h"
#include "Home.h"
#include <petsc.h>

struct _petscVelocityKernelCtx{
    PetscBool   nolinear;       /* Whether to use a nonlinear solver */
    PetscBool   ltog;           /* Assemble Jacobi matrix by using local index mapping global index */
    PetscBool   matFDColoring;  /* Set Jacobi matrix coloring, default true*/
    int         nGlobals;       /* Number of gloal nodes to calculate velocity*/
    int         nNatives;       /* Number of native nodes to calculate velocity*/
    int         *gloIDList;     /* Global index list of all nodes (natives and ghosts) in this domain*/
    int         nGhosts;        /* Number of earth ghost nodes to calculate velocity*/
    Node_t      **nodeList;     /* Node ptr list of all nodes (natives and ghosts) in this domain */
    Vec         velocities;     /* Velocities vector to be solved*/
    Vec         forces;         /* Nodal forces */
    Vec         locVel;         /* Velocities of all nodes in this domain */
    Mat         velCoeMat;      /* Coefficient matrix of nodal velocities, local*/
    Mat         velJacobi;      /* Jocobi matrix of nadal velocities, MPI*/
    int         *nCstList;      /* Num constraits for native node list*/
    real8       *cstList;       /* constraint list for native nodes*/
    VecScatter  scatter;        /* Scatter for velocities*/
    ISLocalToGlobalMapping  isl2g; /* Mapping of local to global indexes*/
    MatFDColoring           matJacColor; /* Matrix Jacobi coloring */
#ifdef PETSC_DEBUG
    FILE                    *fp;
    PetscViewer             viewer;             
#endif
};

/**
 * `PetscInitialize` initializes the PETSc library and reads the options file
 */
PetscErrorCode PetscInitialize(Home_t *home);

/**
 * > GenerateGlobalIndex() assigns a global index to each native node in the domain
 */
PetscErrorCode GenerateGlobalIndex(Home_t *home);

/**
 * PetscVelocityKernelInitialize: initializes the velocity kernel context
 */
PetscErrorCode PetscVelocityKernelInitialize(PetscVKCtx_t *vk);

/**
 * The function `UpdateVelocityForce` is used to update the initial velocities and forces of the nodes
 * 
 * \param[in] home the home structure
 * \param[in,out] vk the velocity kernel
 */
PetscErrorCode UpdateVelocityForce(Home_t *home, PetscVKCtx_t *vk);

/**
 * > This function assembles the global velocity kernel matrix and the global velocity vector
 */
PetscErrorCode AssembleGlobalVelocityKernel(Home_t *home);

/**
 * This function solves the velocity kernel equations using PETSc
 */
PetscErrorCode SolveGlobalVelocityKernel(Home_t *home);

/**
 * > This function destroys the velocity kernel context
 * \param[out] vk the velocity kernel context
 */
PetscErrorCode DestoryGlobalVelocityKernel(PetscVKCtx_t *vk);

/**
 * The function takes in a node ID and returns a 3x3 matrix that is used to penalize the node's motion.
 * 
 * @param home the home structure
 * @param vk the PetscVKCtx_t structure
 * @param nativeID the native ID of the node
 * @param values a 9-element array that will be filled with the constraint tensor.
 */
void NodalConstraintTensor(Home_t *home, PetscVKCtx_t *vk, int nativeID, real8 *values);
#endif