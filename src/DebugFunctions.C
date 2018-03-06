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


#define DIR_DEBUG_GB "debugGB"

void TecplotCurrDomain(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
             int writePrologue, int writeEpilogue, int numArms)
{
        int      i, j, thisDomain;
        int      newNodeKeyPtr, narm;
        int      burgID, normalID; 
	int 	 arm;
        real8    x,y,z;
        real8    x2,y2,z2;
        real8    x2o,y2o,z2o;
        real8    prx,pry,prz;
        real8    vx,vy,vz;
        real8    bX, bY, bZ;
        real8    Lx, Ly, Lz;
        real8     normal[3], normalCryst[3];
        real8	 armInGrain;
        real8	 minDistance, pos[3];
        real8	 velocity, displacement, vecOldNew[3];
        char     fileName[256];
        Node_t   *node, *nbrNode;
        Param_t  *param;
        FILE     *fp;
        Plane_t	 plane;
        struct stat statbuf;


        param      = home->param;
        thisDomain = home->myDomain;

        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;
		
/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s.plt",
                     baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_DEBUG_GB, baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(fileName, &statbuf);
#endif
#endif
	numArms = 0;
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
	    for(j=0; j<node->numNbrs; j++){
                nbrNode = GetNeighborNode(home, node, j);
                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }
		numArms++;
	    }
	}

	while (node) {
		    	
	    if(node->numNbrs > 16 || node->numNbrs < 0){
		node = node->next;
		continue;
	    }

    	    for (arm = 0; arm < node->numNbrs; arm++) {
	        	
	        nbrNode = GetNeighborNode(home, node, arm);
	        if(nbrNode == (Node_t *)NULL){
	            continue;
	        }
		numArms++;
	    }
	    node = node->next;
	}

	firstInGroup = 1;
	writePrologue = 1;

/*
 *      If this process is the first member of its I/O group
 *      it needs to create the output file and do any intializations
 *      needed.
 */
        if (firstInGroup) {
/*
 *          First task in the I/O group must open the data file for writing
 *          to overwrite any existing file of the same name.
 */
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
            }

            if (writePrologue) {
		fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,diplacement,numNbr,ArmInGrain,burgID,normalID,constraint,domain,ID,native\n");	
		VECTOR_ZERO(plane.point);
				
		fprintf(fp, "zone i = %d  F=POINT\n", numArms); 
            }
        } else {
/*
 *          Any task NOT first in its I/O group must open the data file
 *          in an append mode so everything gets added to the end of
 *          the file.
 */
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
            }
        }
 

/*
 *      Generate the plot data for the segments associated with this
 *      domain's data
 */
        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {;
                continue;
            }
        
            x = node->x;
            y = node->y;
            z = node->z;
			
            vecOldNew[XC] = x - node->oldx;
            vecOldNew[YC] = y - node->oldy;
            vecOldNew[ZC] = z - node->oldz;
            	
            ZImage(param, vecOldNew, vecOldNew+1, vecOldNew+2);
            	
            displacement = Normal(vecOldNew);
          
/*
 *          Don't assume PBC is enabled...
 */
            if (param->xBoundType == Periodic) {
                x = x - Lx*rint(x/Lx);
            }

            if (param->yBoundType == Periodic) {
                y = y - Ly*rint(y/Ly);
            }

            if (param->zBoundType == Periodic) {
                z = z - Lz*rint(z/Lz);
            }
        
            for (j = 0; j < node->numNbrs; j++) {
#if 0
/*
 *              This check should be obsolete; should never
 *              have index < zero at this point in execution
 */
                if ((index = node->nbrTag[j].index) < 0) {
                    continue; 
                }
#endif
        
/*
 *  FIX ME? This will result in segments crossing domain boundaries
 *  to be added to the file by both domains!  Is this what is wanted?
 */
				
/*
                if ((node->nbrTag[j].domainID == thisDomain) && 
                    (node->nbrTag[j].index < i)) {
                    continue;
                }
 */
                nbrNode = GetNeighborNode(home, node, j);
				
                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if(home->param->materialType == MAT_TYPE_FCC){
                    burgID = GetBurgID_FCC(home, node, nbrNode);
                }else{
                    Fatal("The burgers ID of non-FCC material has not defined\n");
                }
				
                normal[XC] = node->nx[j];
                normal[YC] = node->ny[j];
                normal[ZC] = node->nz[j];

                if(home->param->useLabFrame){
	            ChangeOrientation(home, node, nbrNode,3);
                    Matrix33Vector3Multiply(home->rotMatrixInverse, normal, normalCryst);
                }else{
                    VECTOR_COPY(normalCryst, normal);
                }

                normalID = GlidePlaneID_FCC(normalCryst);

                x2o = nbrNode->x;
                y2o = nbrNode->y;
                z2o = nbrNode->z;
        
                prx = x2o-x;
                pry = y2o-y;
                prz = z2o-z;
        
                if (param->xBoundType == Periodic) {
                    prx = prx - Lx*rint(prx/Lx);
                }

                if (param->yBoundType == Periodic) {
                    pry = pry - Ly*rint(pry/Ly);
                }

                if (param->zBoundType == Periodic) {
                    prz = prz - Lz*rint(prz/Lz);
                }
        
                x2 = x + prx;
                y2 = y + pry;
                z2 = z + prz;

                vx =  node->vX * param->burgMag;
                vy =  node->vY * param->burgMag;
                vz =  node->vZ * param->burgMag;
	        velocity = sqrt(vx*vx + vy*vy + vz*vz);
                  
                narm = node->numNbrs;
				
                plane.normal[XC] = node->nx[j];
                plane.normal[YC] = node->ny[j];
                plane.normal[ZC] = node->nz[j];
                pos[XC] = x;
                pos[YC] = y;
                pos[ZC] = z;
                	
                minDistance = PointPlaneDistance(pos, &plane);
                	
                armInGrain = ChangeOrientation(home, node, nbrNode, 0);
        
/*
 *              This part is for output of junction investigation
 */
				
	        fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %7.4f %d %3.1f %d %d %d %d %d %d\n", 
	                x,y,z,prx/2,pry/2,prz/2,velocity,displacement,narm,armInGrain,burgID,normalID
	                ,node->constraint,node->myTag.domainID,node->myTag.index,1);
                
            }
        }
		
		
		
	node = home->ghostNodeQ;
		
	while (node) {
			
	    if(node->numNbrs > 16 || node->numNbrs < 0){
	        printf("node(%d,%d,%3.1f)-->position(%e,%e,%e)--->old(%e,%e,%e)",
		        node->myTag.domainID, node->myTag.index, node->inGrain,
			node->x,node->y,node->z,node->oldx,node->oldy,node->oldz);
	        node = node->next;
	        printf("Tecplot (%d), node->numNbrs > 16 || node->numNbrs < 0, 2", home->myDomain);
		continue;
	    }
			
	    x = node->x;
            y = node->y;
            z = node->z;
			
	    vecOldNew[XC] = x - node->oldx;
	    vecOldNew[YC] = y - node->oldy;
	    vecOldNew[ZC] = z - node->oldz;
	    	
	    ZImage(param, vecOldNew, vecOldNew+1, vecOldNew+2);
	    	
	    displacement = Normal(vecOldNew);
	    	
	    if (param->xBoundType == Periodic) {
	        x = x - Lx*rint(x/Lx);
	    }
	    	
	    if (param->yBoundType == Periodic) {
	        y = y - Ly*rint(y/Ly);
	    }
	    	
	    if (param->zBoundType == Periodic) {
	        z = z - Lz*rint(z/Lz);
	    }
	    
	    for (arm = 0; arm < node->numNbrs; arm++) {
	    	
	        nbrNode = GetNeighborNode(home, node, arm);
	        if(nbrNode == (Node_t *)NULL){
	            continue;
	        }
	    	
	        if(home->param->materialType == MAT_TYPE_FCC){
	            burgID = GetBurgID_FCC(home, node, nbrNode);
	        }else{
	            Fatal("The burgers ID of non-FCC material has not defined\n");
	        }
	    	
	        x2o = nbrNode->x;
	        y2o = nbrNode->y;
	        z2o = nbrNode->z;
	        	
	        prx = x2o-x;
	        pry = y2o-y;
	        prz = z2o-z;
	    	
	        if (param->xBoundType == Periodic) {
	            prx = prx - Lx*rint(prx/Lx);
	        }
	        	
	        if (param->yBoundType == Periodic) {
	            pry = pry - Ly*rint(pry/Ly);
	        }
	        	
	        if (param->zBoundType == Periodic) {
	            prz = prz - Lz*rint(prz/Lz);
	        }
	        	
	        x2 = x + prx;
	        y2 = y + pry;
	        z2 = z + prz;
	        	
	        vx =  node->vX * param->burgMag;
	        vy =  node->vY * param->burgMag;
	        vz =  node->vZ * param->burgMag;
	        velocity = sqrt(vx*vx + vy*vy + vz*vz);
	        	
	        narm = node->numNbrs;
	        	
	        plane.normal[XC] = node->nx[arm];
	        plane.normal[YC] = node->ny[arm];
	        plane.normal[ZC] = node->nz[arm];
	        pos[XC] = x;
	        pos[YC] = y;
	        pos[ZC] = z;
	        
                if(home->param->materialType == MAT_TYPE_FCC){
                    burgID = GetBurgID_FCC(home, node, nbrNode);
                }else{
                    Fatal("The burgers ID of non-FCC material has not defined\n");
                }

                normal[XC] = node->nx[j];
                normal[YC] = node->ny[j];
                normal[ZC] = node->nz[j];

                if(home->param->useLabFrame){
	            ChangeOrientation(home, node, nbrNode,3);
                    Matrix33Vector3Multiply(home->rotMatrixInverse, normal, normalCryst);
                }else{
                    VECTOR_COPY(normalCryst, normal);
                }

                normalID = GlidePlaneID_FCC(normalCryst);
	        
	        armInGrain = ChangeOrientation(home, node, nbrNode, 0);
	        
	        fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %7.4f %d %3.1f %d %d %d %d %d %d\n", 
	                x,y,z,prx/2,pry/2,prz/2,velocity,displacement,narm,armInGrain,burgID,normalID
	                ,node->constraint,node->myTag.domainID,node->myTag.index,0);
	        	
	        }
	    	
	        node = node->next;
	    }

 /*
  *     Handle anything that needs to be done after all data has been
  *     written to the file
  */
	    if (writeEpilogue) {
	        fprintf(fp, "\n");
	    }
	     
	    fclose(fp);
	     
	return;
}



void DebugGB(Home_t *home, int flag){
	return;
	int 	i,j,k,l, crossSlipNode, sameNode;
	int	normalID;
	real8   dot,resiBurg[3], resiVec[3];
	Param_t *param;
	Grain_t	*grain;
	Node_t	*node,*nbr;
	real8	inGrainPos, corInGrain; 
	real8   vec1[3],vec2[3], vec[3], vecMag;
	real8	normalCryst[3];
	char    name[256];
        FILE    *fp;
	char    subdir[256], stepdir[256], filename[256];

	param = home->param;
	
        snprintf(subdir, sizeof(subdir), "./%s", DIR_DEBUG_GB);

//	printf("%d %d\n",home->myDomain, flag);

#if SEG_FAULT_CHECK
        snprintf(stepdir, sizeof(stepdir), "%s/%06d", subdir, home->cycle);
        snprintf(filename, sizeof(filename), "%s/%06d", stepdir,
                 home->myDomain);
        snprintf(name, sizeof(name), "%s.%d.plt", filename,
                 flag);

        if (mkdir(subdir, S_IRWXU) != 0) {
            if (errno != EEXIST) {
                Fatal("Open error %d on directory %s",
                      errno, subdir);
            }
        }
        if (mkdir(stepdir, S_IRWXU) != 0) {
            if (errno != EEXIST) {
                Fatal("Open error %d on directory %s",
                      errno, stepdir);
            }
        }

        if ((fp = fopen(filename, "a")) == (FILE *)NULL) {
            Fatal("Open error on debugGB error data file %s", filename);
        }

	fprintf(fp, "%d\n", flag); 
	fclose(fp);

	if(flag == 6){
        //    TecplotCurrDomain(home, filename, 0, 0, 0, 0, 0);
	}
#endif
	for(i=0; i<home->newNodeKeyPtr; i++){
	    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
		continue;
	    }
            if ((node->numNbrs == 1) && (node->constraint == UNCONSTRAINED)) {
                PrintNode(node);
                snprintf(name,50,"DebuGBSingly%06d.%d.%d",home->cycle,home->myDomain,flag);
                Tecplot(home, name, 0, 0, 0, 0, 0, 1);
                Fatal("DebugGB: Node (%d,%d) singly linked!, flag %d",
                      node->myTag.domainID, node->myTag.index, flag);
            }
//	    if(node->inGrain -floor(node->inGrain) == 0)continue;
//
	    VECTOR_ZERO(resiBurg);
	    for(j=0; j<node->numNbrs; j++){
	        resiBurg[XC] += node->burgX[j];
	        resiBurg[YC] += node->burgY[j];
	        resiBurg[ZC] += node->burgZ[j];

		nbr = GetNeighborNode(home,node,j);
		if(nbr == (Node_t *)NULL)continue;
		
		sameNode = 0;
		VECTOR_ZERO(resiVec);
		for(k=0; k<nbr->numNbrs; k++){
		    resiVec[XC] += nbr->burgX[k];
		    resiVec[YC] += nbr->burgY[k];
		    resiVec[ZC] += nbr->burgZ[k];

		    if(nbr->nbrTag[k].domainID == node->myTag.domainID &&
		       nbr->nbrTag[k].index == node->myTag.index){
			sameNode = 1;
		    }
		}
#if 1
		if(Normal(resiVec) > 0.01 && nbr->constraint == 0){
                    snprintf(name,50,"NW%06d.%d.%d",home->cycle,home->myDomain,flag);
                    Tecplot(home, name, 0, 0, 0, 0, 0, 1);
		    PrintNode(node);
		    PrintNode(nbr);
		    printf("Neighbor Wrong\n");
		}	
#endif
		if(sameNode == 0){
                    snprintf(name,50,"WC%06d.%d.%d",home->cycle,home->myDomain,flag);
                    Tecplot(home, name, 0, 0, 0, 0, 0, 1);
		    PrintNode(node);
		    PrintNode(nbr);
		    printf("Wrong connection\n");
		}

#if 0
                if(nbr->myTag.domainID < node->myTag.domainID || 
                   (nbr->myTag.domainID == node->myTag.domainID &&
                    nbr->myTag.index < node->myTag.index))continue;
#endif
		    if(node->inGrain != nbr->inGrain &&
		       node->inGrain - floor(node->inGrain) == 0 &&
		       nbr->inGrain - floor(nbr->inGrain) == 0 &&
	               !SurfaceNode(node) && !SurfaceNode(nbr)){
			PrintNode(node);
			PrintNode(nbr);
		        Fatal("\nTask %d, flag %d, Wrong inGrain\n",home->myDomain,flag);
		    }

		vec2[XC] = nbr->x - node->x;
		vec2[YC] = nbr->y - node->y;
		vec2[ZC] = nbr->z - node->z;
		ZImage(param,vec2,vec2+1,vec2+2);
		if(Normal(vec2) < 1e-12){
		    PrintTag(node);
		    PrintTag(nbr);
		    printf("\nTask %d, flag %d, zero length segement\n",home->myDomain,flag);
		}	
		
		vec[XC] = node->nx[j];
		vec[YC] = node->ny[j];
		vec[ZC] = node->nz[j];
                
                if(home->param->useLabFrame){
	            ChangeOrientation(home, node, nbr,3);
                    Matrix33Vector3Multiply(home->rotMatrixInverse, vec, normalCryst);
                }else{
                    VECTOR_COPY(normalCryst, vec);
                }

                normalID = GlidePlaneID_FCC(normalCryst);

		if(normalID < 4 && node->numNbrs < 3){
		    NormalizeVec(vec2);
		    dot = DotProduct(vec2,vec);
		    if(fabs(dot) > 2.0e-1){
//			printf("DebugNormal, (%d,%d) dot is %f\n",node->myTag.domainID, node->myTag.index,dot);
		    }
		}
		
	        corInGrain = InWhichGrain(home, node->x, node->y, node->z);
		if(node->inGrain - floor(node->inGrain) ==0 && normalID == 3){
		    if(corInGrain != -1 && corInGrain != node->inGrain && flag > 8 && flag < 11 || 
	               nbr->inGrain - floor(nbr->inGrain) == 0.5){
                        snprintf(name,50,"DebuGBNormal%06d.%d.%d.2",home->cycle,home->myDomain,flag);
                        Tecplot(home, name, 0, 0, 0, 0, 0, 1);
			PrintNode(node);
			PrintNode(nbr);
			Fatal("segement which normal ID is 3 penetrated though TB\n");
			break;
		    }
		}

	    }

	    if(node->constraint == 0 && Normal(resiBurg) > 0.01){
                    snprintf(name,50,"RB%06d.%d.%d",home->cycle,home->myDomain,flag);
                    Tecplot(home, name, 0, 0, 0, 0, 0, 1);
		    PrintNode(node);
		    Fatal("resiBurg > 0\n");
	    }
	}
}
