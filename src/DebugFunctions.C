/*-------------------------------------------------------------------------
 *
 *      Module:       DebugFunctions.c
 *
 *      Description:  This module contains a variety of functions
 *                    that are intended for debug purposes only and
 *                    are not normally invoked by the code.  During
 *                    the debugging process, calls to these functions
 *                    can be inserted in the code as necessary.
 *
 *
 *      Included functions:
 *
 *          CheckSegLengths()
 *          CheckForces()
 *          CheckForNANS()
 *          CheckForEmptySimulation()
 *          CheckForUndefinedPlanes()
 *
 *-----------------------------------------------------------------------*/
#include "Home.h"

#include "sys/stat.h"
#include "sys/types.h"
#ifdef PARALLEL
#include "mpi.h"
#endif


/**************************************************************************
 *
 *      Function:     CheckForces
 *      Description:  Check the force and velocity components of each
 *                    node looking for any with values exceeding some
 *                    hard-code thresholds.  If such a node is found,
 *                    the code aborts with an error identifying the offending
 *                    node.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               a long segment is found.
 *
 *************************************************************************/
void CheckForces(Home_t *home, const char *msg)
{
        int    i, j;
        int    numNodes, numNbrs;
        real8  forceThreshold, velThreshold;
        Node_t *node;

        forceThreshold = 1.0e+14;
        velThreshold   = 1.0e+14;

        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node->numNbrs;

            if ((node->fX > forceThreshold) ||
                (node->fY > forceThreshold) ||
                (node->fZ > forceThreshold)) {
#if 1
                    Fatal("%s:  Node (%d,%d) force = %e %e %e\n", 
                           msg, node->myTag.domainID, node->myTag.index,
                           node->fX, node->fY, node->fZ);
#else
                    printf("%s:  Node (%d,%d) force = %e %e %e\n", 
                           msg, node->myTag.domainID, node->myTag.index,
                           node->fX, node->fY, node->fZ);
#endif
            }

            if ((node->vX > velThreshold) ||
                (node->vY > velThreshold) ||
                (node->vZ > velThreshold)) {
#if 1
                    Fatal("%s:  Node (%d,%d) velocity = %e %e %e\n", 
                           msg, node->myTag.domainID, node->myTag.index,
                           node->vX, node->vY, node->vZ);
#else
                    printf("%s:  Node (%d,%d) velocity = %e %e %e\n", 
                           msg, node->myTag.domainID, node->myTag.index,
                           node->vX, node->vY, node->vZ);
#endif
            }

        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckSegLengths
 *      Description:  Check the length of every segment owned by the 
 *                    current domain and print out an error message
 *                    any time a segment exceeding the maximum segment
 *                    lenth <maxSeg> by 10% or more is found.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               a long segment is found.
 *
 *************************************************************************/
void CheckSegLengths(Home_t *home, const char *msg)
{
        int    i, j;
        int    numNodes, numNbrs;
        real8  xLen, yLen, zLen, totLen;
        Node_t *node1, *node2;


        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node1->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                node2 = GetNeighborNode(home, node1, j);

                if (node2 == (Node_t *)NULL) {
                    printf("%s:  Segment (%d,%d)--(%d,%d) neighbor not found\n",
                           msg, node1->myTag.domainID, node1->myTag.index,
                           node1->nbrTag[j].domainID, node1->nbrTag[j].index);
                    fflush(NULL);
                    continue;
                }

                xLen = node2->x - node1->x;
                yLen = node2->y - node1->y;
                zLen = node2->z - node1->z;

                ZImage(home->param, &xLen, &yLen, &zLen);

                totLen = sqrt(xLen*xLen + yLen*yLen + zLen*zLen);

                if (totLen > (home->param->maxSeg * 1.1)) {
                    printf("%s:  Segment (%d,%d)--(%d,%d) length = %lf\n",
                           msg, node1->myTag.domainID, node1->myTag.index,
                           node2->myTag.domainID, node2->myTag.index, totLen);
                    fflush(NULL);
                }
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckForNaNs
 *      Description:  Check for NaNs (or Infs) in the position, force and 
 *                    velocity components for each local node.  If a NaN 
 *                    or Inf is found, the code will abort.
 *
 *************************************************************************/
void CheckForNANS(Home_t *home, real8 flag)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if (isnan(node->x) ||
                isnan(node->y) ||
                isnan(node->z)) {
                PrintNode(node);
                Fatal("Node position contains NaNs.  Aborting! %3.1f",flag);
            }
             
            if (isinf(node->x) ||
                isinf(node->y) ||
                isinf(node->z)) {
                PrintNode(node);
                Fatal("Node position contains Infs.  Aborting! %3.1f",flag);
            }
             
            if (isnan(node->vX) ||
                isnan(node->vY) ||
                isnan(node->vZ)) {
                PrintNode(node);
                Fatal("Node velocity contains NaNs.  Aborting! %3.1f",flag);
            }
             
            if (isinf(node->vX) ||
                isinf(node->vY) ||
                isinf(node->vZ)) {
                PrintNode(node);
                Fatal("Node velocity contains Infs.  Aborting! %3.1f",flag);
            }
             
            if (isnan(node->fX) ||
                isnan(node->fY) ||
                isnan(node->fZ)) {
                PrintNode(node);
                Fatal("Node force contains NaNs.  Aborting! %3.1f",flag);
            }
             
            if (isinf(node->fX) ||
                isinf(node->fY) ||
                isinf(node->fZ)) {
                PrintNode(node);
                Fatal("Node force contains Infs.  Aborting! %3.1f",flag);
            }
             
        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckForEmptySimulation
 *      Description:  Just a sanity check to kill the simulation if we
 *                    end up in a situation with no more dislocations in
 *                    the problem.
 *
 *************************************************************************/
void CheckForEmptySimulation(Home_t *home)
{
        int    i, localNodeCnt, globalNodeCnt;
        Node_t *node;

        localNodeCnt = 0;
        globalNodeCnt = 0;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            } else {
                localNodeCnt++;
            }
        }

#ifdef PARALLEL
        MPI_Reduce(&localNodeCnt, &globalNodeCnt, 1, MPI_INT, MPI_SUM,
                    0, MPI_COMM_WORLD);
#else
        globalNodeCnt = localNodeCnt;
#endif
        if ((home->myDomain == 0) && (globalNodeCnt == 0)) {
            Fatal("All dislocations in simulation have been "
                  "annihilated.\nSimulation terminating NOW!");
        }
}


/**************************************************************************
 *
 *      Function:     CheckForUndefinedPlanes
 *      Description:  Check for any local segment for which the glide
 *                    plane is undefined.  If found, print the nodal
 *                    information and an error message.
 *                    
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               a segment with an undefined glide plane is found.
 *
 *************************************************************************/
void CheckForUndefinedPlanes(Home_t *home, const char *msg)
{
        int    i, j;
        int    numNodes, numNbrs;
        Node_t *node;


        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                if ((fabs(node->nx[j]) < 1.0e-4) &&
                    (fabs(node->ny[j]) < 1.0e-4) &&
                    (fabs(node->nz[j]) < 1.0e-4)) {
                    printf("Error: Undefined glide plane\n");
                    PrintNode(node);
                }
            }
        }

        return;
}

void CheckAll(Home_t *home, const char *format, ...)
{
	return;
	char    msg[512];
	va_list args;

	va_start(args, format);
	vsnprintf(msg, sizeof(msg)-1, format, args);
	msg[sizeof(msg)-1] = 0;
	va_end(args);


	
	printf("Debug (%d): %s\n",home->myDomain, msg);
	
	Node_t 		*node;
	int			i;
	
	for(i=0; i<home->newNodeKeyPtr; i++){
		if((node = home->nodeKeys[i]) == (Node_t *)NULL)continue;
		if(node->numNbrs > 16 || node->numNbrs < 0){
			PrintTag(node);
			Fatal("Debug (%d): %s, native node->numNbrs > 16 || node->numNbrs < 0, 1", home->myDomain,msg);
		}
	}
	
	node = home->ghostNodeQ;
	while (node) {
		
		if(node->numNbrs > 16 || node->numNbrs < 0){
			PrintTag(node);
			Fatal("Debug (%d): %s, node->numNbrs > 16 || node->numNbrs < 0, 2", home->myDomain,msg);
		}
		
		node = node->next;
	}

}

typedef struct {
	Segment_t *seg1;
	Segment_t *seg2;
	int       setSeg1Forces;
	int       setSeg2Forces;
} SegmentPair_t;


typedef struct {
	Segment_t *seg;
	Cell_t    *cell;
	int       cellID;
} NativeSeg_t;

void CheckLocal(Home_t *home, int flag){
	return;
	int        i, cellNum, cellID, cellSegCnt;
	int        arm, armCount;
	int        homeDomain, homeCells, homeNativeCells;
	int        totalNativeSegs;
	int        sendDomCnt, allocSize, numDomains;
	int        setSeg1Forces, setSeg2Forces;
	int        *nativeSegCounts, *totalSegCounts;
	int        *sendDomList, *globalMsgCnts, *localMsgCnts;
	int        segPairListCnt = 0, segPairListSize = 0;
	int        nativeSegListCnt = 0;
	real8      MU, NU, a, Ecore, extstress[3][3];
	Node_t     *node1, *node2, *node3, *node4;
	Node_t     *node, *nbr;
	Cell_t     *cell;
	Param_t    *param;
	
	
	homeCells       = home->cellCount;
	homeDomain      = home->myDomain;
	param           = home->param;

	printf("CheckLocal: Task %d, flag %d \n",home->myDomain,flag);
	
	for (cellNum = 0 ; cellNum < homeCells; cellNum++) {
	
		cell = home->cellKeys[home->cellList[cellNum]];
		
	
		for (node = cell->nodeQ; node != (Node_t *)NULL; node = node->nextInCell) {
	
			if (node->myTag.domainID == homeDomain) continue;
			
			if(node->numNbrs < 0 || node->numNbrs > 16){
				PrintTag(node);
				Fatal("CheckLocal: Task %d,node->numNbrs < 0 || node->numNbrs > 16 at %s line %d, flag is %d\n",
					  home->myDomain, __FILE__, __LINE__, flag);
			}
	
			for (arm = 0; arm < node->numNbrs; arm++) {
	
				if(node->nbrTag[arm].domainID < 0 || 
					node->nbrTag[arm].domainID >= home->numDomains||
				    node->nbrTag[arm].index < 0){
					Fatal("CheckLocal: Task %d, Neighbor not found at %s line %d, flag is %d\n",
						  home->myDomain, __FILE__, __LINE__, flag);
					
				}
#if 0	
				if (nbr == (Node_t *)NULL) {
					PrintNode(node);
					Fatal("CheckLocal: Task %d, Neighbor not found at %s line %d, flag is %d\n",
						  home->myDomain, __FILE__, __LINE__, flag);
				
				}
#endif	
			}
		}
	
	}
	return;
}

void CheckNodeNANS(Node_t *node){
	int flag=0;

	if(node == (Node_t *)NULL)return;

            if (isnan(node->x) ||
                isnan(node->z)) {
                PrintNode(node);
                Fatal("Node position contains NaNs.  Aborting! %3.1f",flag);
            }
             
            if (isinf(node->x) ||
                isinf(node->y) ||
                isinf(node->z)) {
                PrintNode(node);
                Fatal("Node position contains Infs.  Aborting! %3.1f",flag);
            }
             
            if (isnan(node->vX) ||
                isnan(node->vY) ||
                isnan(node->vZ)) {
                PrintNode(node);
                Fatal("Node velocity contains NaNs.  Aborting! %3.1f",flag);
            }
             
            if (isinf(node->vX) ||
                isinf(node->vY) ||
                isinf(node->vZ)) {
                PrintNode(node);
                Fatal("Node velocity contains Infs.  Aborting! %3.1f",flag);
            }
             
            if (isnan(node->fX) ||
                isnan(node->fY) ||
                isnan(node->fZ)) {
                PrintNode(node);
                Fatal("Node force contains NaNs.  Aborting! %3.1f",flag);
            }
             
            if (isinf(node->fX) ||
                isinf(node->fY) ||
                isinf(node->fZ)) {
                PrintNode(node);
                Fatal("Node force contains Infs.  Aborting! %3.1f",flag);
            }
}



void DebugGB(Home_t *home, int flag){
	int 	i,j,k,l, crossSlipNode, sameNode;
	int	    normalID, inGrain, onGB;
	real8   dot,resiBurg[3], resiVec[3];
	Param_t *param;
	Grain_t	*grain;
	Node_t	*node,*nbr;
	int 	inGrainPos, corInGrain, corOnGB; 
	real8   vec1[3],vec2[3], vec[3], vecMag;
	real8	normalCryst[3];
	char    name[256];
        FILE    *fp;
	char    subdir[256], stepdir[256], filename[256];

	param = home->param;

//	printf("%d %d\n",home->myDomain, flag);

	for(i=0; i<home->newNodeKeyPtr; i++){
		if((node = home->nodeKeys[i]) == (Node_t *)NULL)continue;

        if(flag>22 || flag <= 9){
            if(!IsOnGB(node)){
                NodeLocation(home, node->x, node->y, node->z, inGrain, onGB);
                if(node->inGrain != inGrain){
                    PrintTag(node);
                    printf(" DebugGB - 1: %d %d\n", home->myDomain, flag);
                }
            }
        }

        if(IsOnGB(node))continue;
        for(j=0; j<node->numNbrs; j++){
            nbr = GetNeighborNode(home, node, j);
            if(nbr == (Node_t *)NULL)continue;
            if(!IsOnGB(nbr) && nbr->inGrain != node->inGrain && flag < 11){
                PrintTag(node);
                PrintTag(nbr);
                Fatal("DebugGB - 2: %d %d", home->myDomain, flag);
            }
        }
    }

    return;
	node = home->ghostNodeQ;

    while (node) {
        if(node->inGrain == -1 && node->onGB == -1){
            PrintNode(node);
            Fatal("DebugGB %d %d", home->myDomain, flag);
        }
    
        node = node->next;
    }
}

