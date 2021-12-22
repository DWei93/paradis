/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisStep
 *      Description:  This function controls everything needed for a
 *                    single step of a ParaDiS simulation including
 *                    force calculations, ghost cell communications,
 *                    node migration, dynamic load balance, output
 *                    generation, etc.
 *
 *-----------------------------------------------------------------------*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "Home.h"
#include "Util.h"
//#include "DisplayC.h"
#include "Comm.h"
#include "Mobility.h"
#include "Decomp.h"
#include "QueueOps.h"
#include <mcheck.h>
#define MALLOC_TRACE mem.log 

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#ifdef PARALLEL
#include "mpi.h"
#endif

/*
 *      By default, there are no runtime checks to see if all of
 *      the dislocations have annihilated themselves.  To enable
 *      a check with an abort if it happens, simply define the
 *      DEBUG_CHECK_FOR_ZERO_SEG value below to 1 rather than zero.
 */
#define DEBUG_CHECK_FOR_ZERO_SEG 0


/*
 *      For debugging only.  If DEBUG_STEP is not defined, all
 *      calls to Synchronize() will be replaced with an empty
 *      block of code, but if it is defined, the calls will
 *      be replaced with a call to syncronize the code and log
 *      a message.
 */
#ifdef DEBUG_STEP
#define Synchronize(a,b) _Synchronize((a),(b))
#else
#define Synchronize(a,b) {}
#endif

/*
 *      Explicitly synchronize parallel execution via an MPI
 *      barrier and print a message when all tasks have reached
 *      the barrier.  For debug only.
 */
void _Synchronize(Home_t *home, char *msg)
{

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if (home->myDomain == 0) {
            printf(" *** %s: All tasks synchronized\n", msg);
            fflush(NULL);
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ApplyDeltaStress
 *      Description: Increment the force/vel for all native nodes based
 *                   on the change in applied stress (deltaStress)
 *                   during a timestep.  Since the applied stress plays
 *                   no impact on the segment/segment interactions this
 *                   is much less expensive than recomputing all the n^2
 *                   segment interactions.
 *
 *------------------------------------------------------------------------*/
static void ApplyDeltaStress(Home_t *home, real8 deltaStress[3][3])
{
        int     i, j, nbrArm, nbrIsLocal;
        real8   x1, y1, z1, x2, y2, z2;
        real8   bx1, by1, bz1, dx, dy, dz;
        real8   f1[3], f2[3];
        Node_t  *node, *nbr;
        Param_t *param;
        
        param = home->param;
		

/*
 *      Loop over all native nodes
 */        
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
        
            x1=node->x;
            y1=node->y;
            z1=node->z;
        
/*
 *          For each node, recalculate forces for all the node's arms
 *          that are either owned by this node, or terminate non-locally.
 */
            for (j = 0; j < node->numNbrs; j++) {

                nbr = GetNeighborNode(home, node, j);

                if (nbr == (Node_t *)NULL) {
					PrintNode(node);
                    printf("WARNING: Neighbor %d (%d,%d) not found at %s line %d\n",
                           j,node->nbrTag[j].domainID, node->nbrTag[j].index,__FILE__, __LINE__);
                    continue;
                }

                nbrIsLocal = (nbr->myTag.domainID == home->myDomain);

                if (nbrIsLocal) {
                    if (OrderNodes(node, nbr) >= 0) continue;
                    nbrArm = GetArmID(home, nbr, node);
                }
        
                bx1 = node->burgX[j];
                by1 = node->burgY[j];
                bz1 = node->burgZ[j];
        
                dx=nbr->x-x1;
                dy=nbr->y-y1;
                dz=nbr->z-z1;
        
                ZImage(param, &dx, &dy, &dz) ;
        
                x2=x1+dx;
                y2=y1+dy;
                z2=z1+dz;
				
				
                ExtPKForce(deltaStress, bx1, by1, bz1, x1, y1, z1,
                        x2, y2, z2, f1, f2);
				
				
        
                AddtoNodeForce(node,f1);
                AddtoArmForce(node, j, f1);
				
                if (nbrIsLocal) {
                    AddtoNodeForce(nbr, f2);
                    AddtoArmForce(nbr, nbrArm, f2);
                }
            }

            (void)EvaluateMobility(home, node);
        }
		
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReevaluateForces
 *      Description: Look for any local nodes whose force/velocity data
 *                   have been flagged as obsolete.  For all such nodes
 *                   recompute the forces and velocities.
 *
 *------------------------------------------------------------------------*/
void ReevaluateForces(Home_t *home)
{
        int     i;
        Node_t  *node;
        Param_t *param;

        param = home->param;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
            if (node->flags & NODE_RESET_FORCES) {

                SetOneNodeForce(home, node);
                EvaluateMobility(home, node);

                node->flags &= (~NODE_RESET_FORCES);
            }
        }

        return;
}


void ParadisStep(Home_t *home)
{
        int        i;
        int        doAll = 1;
        real8      deltaStress[3][3];
        Param_t    *param;
        Node_t     *node;
        char       name[50];
        static int firstTime = 1;
   

        param = home->param;

/*
 *      If this step is an initial load-balance-only step just
 *      perform the minimal work needed to estimate per-process
 *      load, shift boundaries, and migrate nodes among processors.
 */

        if (param->numDLBCycles > 0) {
/*
 *          Note:  When param->numDLBCycles > 0, NodeForce() assumes
 *          the cycle is a DLB-only cycle and only counts the number
 *          of force calcs that would be done without actually calculating
 *          any forces.
 */
            NodeForce(home, FULL);
            Rebalance(home, DLB_USE_FORCECALC_COUNT);

/*
 *          Any time the boundaries are changed, we need to migrate
 *          nodes to their new owning domains and go through all the
 *          ghost node communication stuff.
 */
            Migrate(home);
            RecycleGhostNodes(home);
            SortNativeNodes(home);
            CommSendGhosts(home);

            home->cycleForceCalcCount = 0;

            return;
        }
		

#ifndef NO_XWINDOW
        while (WinIsPaused()) {
            sleep(1);
        }
#endif

        DebugGB(home,15);
/*
 *      Calculate the net charge tensor for each cell (includes global comm)
 */
        CellCharge(home);
		
        DebugGB(home,16);
/*
 *      Calculate new force and velocity data for all nodes or a selected
 *      subset and distribute the new data out to neighboring domains.
 *      The first cycle we'll recalculate all forces and velocities
 *      We do this to get an initial estimate of forces on the first cycle,
 *      After that, we only need to recompute values for nodes that were
 *      involved in topological changes the previous step.
 */
        if (firstTime) {
            NodeForce(home, FULL);
            CalcNodeVelocities(home, 1, doAll);
            CommSendVelocity(home);
            firstTime = 0;

        } else {
            NodeForce(home, PARTIAL);
            CalcNodeVelocities(home, 0, doAll);
            CommSendVelocity(home);

        }


        DebugGB(home,17);
/*
 *      Invoke the selected time step integration method.  The
 *      selected method will calculate the time step as well as
 *      move nodes to their correct locations and communicate 
 *      the new nodal force/velocity data to neighboring domains.
 */

		
        if (strcmp(param->timestepIntegrator, "forward-euler") == 0) {
            ForwardEulerIntegrator(home);
        } else if (strcmp(param->timestepIntegrator, "trapezoid") == 0) {
            TrapezoidIntegrator(home);
        } else {
/*
 *          Used to be specified as 'backard-euler', so if integration
 *          method is unrecognized, use trapezoid as default
 */
            TrapezoidIntegrator(home);
        }
        DebugGB(home,18);
//      DumpDeltaPlasticStrain(home, 0);


/*
 *      Before doing topological changes, set flags indicating any
 *      nodes exempt from topological changes.  These flags are used
 *      in both splitting multi-arm nodes and collisions, so this
 *      function should be invoked before either of those items are done.
 */
        InitTopologyExemptions(home);
        DebugGB(home,19);

#if CHECK_GB_TOPO
        GenerateOutput(home, STAGE_DEBUG, "-2Trap");
#endif

        DebugGB(home,20);

        if(param->xBoundType == Free ||
           param->yBoundType == Free ||
           param->zBoundType == Free ||
           param->gbStrength > 0){
            ClearOpList(home);

            AdjustSurfaceNode(home);

            TimerStart(home, SEND_FS);
            CommSendRemesh(home);
            TimerStop(home, SEND_FS);
	
            TimerStart(home, FIX_FS);
            FixRemesh(home, "Surface");
            TimerStop(home, FIX_FS);
#if CHECK_GB_TOPO
            GenerateOutput(home, STAGE_DEBUG, "-3FS");
#endif
        }

//      DumpDeltaPlasticStrain(home, 1);
        DebugGB(home,21);
        if(param->gbStrength>0 || param->polycrystal && param->nGrain>1){
            ClearOpList(home);
            HandleGrainBoundary(home);
            DebugGB(home,22);

            TimerStart(home, SEND_GB);
            CommSendRemesh(home);
            TimerStop(home, SEND_GB);
	
            TimerStart(home, FIX_GB);
            FixRemesh(home, "GB");
            TimerStop(home, FIX_GB);

        }

        DebugGB(home,23);

/*
 *      In some simulations, it is necessary to recalculate and distribute
 *      the glide plane infromation for segments after the nodes have been
 *      repositioned.  Do so now if needed.
 */
        ResetGlidePlanes(home);
        DebugGB(home,24);

/*
 *      Increment the per-burgers vector density gain/loss with
 *      changes for this cycle.  This must be done immediately
 *      after timestep integration!
 *
 *      Note: This is currently only applicable to BCC simulations.
 */
        GetDensityDelta(home);
        DebugGB(home,25);

        DebugGB(home,26);

/*
 *      The call to GenerateOutput will update the time and cycle counters,
 *      determine if any output needs to be generated at this stage, and
 *      call the appropriate I/O functions if necessary.
 */
//      DumpDeltaPlasticStrain(home, 2);
        DebugGB(home,1);

/*
 *      Now do all the topological changes from segment interactions
 *      (collisions, multinode splitting)...  Clear the list of local
 *      operations that will be sent to the remote domains for processsing,
 *      then split any multi-arm nodes that need splitting, cross slip
 *      nodes (as needed/allowed), handle all local collisions, then
 *      send remote nodes the list of ops needed to keep their data in sync.
 */
        ClearOpList(home);
		
        SortNodesForCollision(home);

        SplitMultiNodes(home);

        DebugGB(home,2);
#ifdef _FEM
        SplitSurfaceNodes(home);
#endif

/*
 *      Call a generic cross-slip dispatch function that will call
 *      (if necessary) the cross-slip function appropriate to the
 *      type of material in use.
 */
//      DumpDeltaPlasticStrain(home, 3);
        CrossSlip(home);
        DebugGB(home,3);
#if CHECK_GB_TOPO
        GenerateOutput(home, STAGE_DEBUG, "-4CrossSlip");
#endif
/*
 *      Search for dislocation segments in close proximity to each other
 *      and if necessary handle any collision between them.
 */
//      DumpDeltaPlasticStrain(home, 4);
        HandleCollisions(home);
		
//      DumpDeltaPlasticStrain(home, 5);
        DebugGB(home,4);

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, POST_COLLISION_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, POST_COLLISION_BARRIER);
#endif
#endif
		
        TimerStart(home, COL_SEND_REMESH);
        CommSendRemesh(home);
        TimerStop(home, COL_SEND_REMESH);
		
        DebugGB(home,5);
        TimerStart(home, COL_FIX_REMESH);
        FixRemesh(home, "Collision");
        TimerStop(home, COL_FIX_REMESH);
		
        DebugGB(home,6);

#ifdef _FEM
	AdjustNodePosition(home, 1); 
#endif 

/*
 *      Under certain circumstances, parallel topological changes can
 *      create double links between nodes; links which can not be detected
 *      until after FixRemesh() is called... so, a quick check has to be
 *      done to clean up these potential double-links here, or they will
 *      cause problems later on.  Should only have to check nodes local
 *      to this domain.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if(RemoveDoubleLinks(home, node, 0) == 2){
                snprintf(name,50,"RDL%d.1",home->cycle);
                Tecplot(home, name, 0, 0, 0 ,0, 0, 1);
                printf("Unusual RemoveDoubleLinks , 1\n");
            }
            node->flags &= ~NODE_CHK_DBL_LINK;
        }

#if CHECK_GB_TOPO
        GenerateOutput(home, STAGE_DEBUG, "-5Collision");
#endif
/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif
       
/*
 *      Invoke mesh coarsen/refine
 */
        DebugGB(home,7);
        Remesh(home);
        DebugGB(home,8);


#if CHECK_GB_TOPO
        GenerateOutput(home, STAGE_DEBUG, "-6Remesh");
#endif
/*
 *      Calculate the new plastic strain.
 */
        DeltaPlasticStrain(home);
 //     DumpDeltaPlasticStrain(home, 6);
        DebugGB(home,9);
        GenerateOutput(home, STAGE_CYCLE,"-1cycle");
        DebugGB(home,10);
/*

 *      Define load curve and calculate change in applied stress this cycle
 */
        LoadCurve(home, deltaStress);
		

/*
 *      This is only needed when we do force calcs for only
 *      a subset of the nodes at the beginning of the timestep.  It will
 *      adjust the nodal forces based on the current delta stress and
 *      recalculate the nodal velocities so we have more accurate values
 *      when we enter the timestep integrator at the beginning of the next
 *      cycle.
		 
 */
        ApplyDeltaStress(home, deltaStress);

/*
 *      If necessary, use the current load data to generate a new
 *      domain decomposition to rebalance the workload among the
 *      processors.
 */
        Rebalance(home, DLB_USE_WALLCLK_TIME);
//        DebugGB(home,10);

/*
 *      Send any nodes that have moved beyond the domain's
 *      boundaries to the domain the node now belongs to.
 */

        Migrate(home);
//        DebugGB(home,11);
/*
 *      Recycle all the ghost nodes: move them back to the free Queue
 */
        RecycleGhostNodes(home);
//        DebugGB(home,12);

/*
 *      Sort the native nodes into their proper subcells.  
 */
		
        SortNativeNodes(home);
//       DebugGB(home,13);

/*
 *      Communicate ghost cells to/from neighbor domains
 */
		
        CommSendGhosts(home);
//     DebugGB(home,14);

#ifdef NAN_CHECK
/*
 *      For debug only:  Abort if any of the nodes have position or
 *      velocity values that are NaNs or infinites.  Be sure to do this
 *      before we write the restart files so we don't get bad data
 *      in the restart.
 */
        CheckForNANS(home, 0);
#endif

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
 #ifdef DEBUG_MEM
         ParadisMemCheck();
#endif

        CheckMemUsage(home, "ParadisStep-complete");

/*
 *      Zero out the count of force calculations done this cycle
 *      so the load-balancing in the next step is based on accurate
 *      values.
 */
        home->cycleForceCalcCount = 0;

/*
 *      For Debugging and testing only...
 */
#if DEBUG_CHECK_FOR_ZERO_SEG
        CheckForEmptySimulation(home);
#endif

	
        return;
}
