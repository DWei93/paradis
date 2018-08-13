/*****************************************************************************
 *
 *      Module:         RemeshRule_2B.c
 *      Description:    This module contains functions to coarsen
 *                      or refine the mesh topology.
 *                      
 *      Included functions:
 *              MeshCoarsen()
 *              MeshRefine()
 *              RemeshRule_2()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "Home.h"
#include "Util.h"
#include "QueueOps.h"
#include "Mobility.h"

static int dbgDom;

/*-------------------------------------------------------------------------
 *
 *      Function:    MeshCoarsen
 *      Description: 
 *
 *------------------------------------------------------------------------*/
static void MeshCoarsen(Home_t *home)
{
        int     i,j, k, q, thisDomain, mergeDone, mergeStatus, globalOp;
        int     localCoarsenCnt, globalCoarsenCnt, hasRemoteNbr;
	    int	    breakLoop = 0, type;
        real8   tmpMag2, vec1[3], vec2[3];
        real8   cellLength, cutoffLength1, cutoffLength2;
        real8   vec1x, vec1y, vec1z;
        real8   vec2x, vec2y, vec2z;
        real8   vec3x, vec3y, vec3z;
        real8   r1, r2, r3;
        real8   s, area2, areaMin, areaMin2, delta;
        real8   dvec1xdt, dvec1ydt, dvec1zdt;
        real8   dvec2xdt, dvec2ydt, dvec2zdt;
        real8   dvec3xdt, dvec3ydt, dvec3zdt;
        real8   dr1dt, dr2dt, dr3dt, dsdt, darea2dt;
        real8   newPos[3], f0seg1[3], f1seg1[3];
        real8   gp0[3], gp1[3], tmp3[3];
        real8   burg[3],vec[3],vecStored;
        real8   dis;
	    real8	mergedGrain;
	    int	    riskySeg1, riskySeg2, armID;    
        int     armInGrain1, armInGrain2, armOnGB1, armOnGB2;   
        Tag_t   nbr1Tag, nbr2Tag, oldTag1, oldTag2, oldTag3;
        Node_t  *node, *nbr, *nbr1, *nbr2, *mergedNode;
        Param_t *param;
        Node_t  bkupNodeList[2];
	    char 	name[50];
        Plane_t GBPlane;
        bool    inhibitMerge,violateGlidePlanesOK = 0; 
        int     mergedBasic;
        int     planeType;

        thisDomain = home->myDomain;
        param      = home->param;

        cellLength = home->param->Lx / home->param->nXcells;
        cellLength = MIN(cellLength, home->param->Ly / home->param->nYcells);
        cellLength = MIN(cellLength, home->param->Lz / home->param->nZcells);

        cutoffLength1 = param->maxSeg;
        cutoffLength2 = MIN(cutoffLength1, 0.45 * cellLength);

        areaMin    = param->remeshAreaMin;
        areaMin2   = areaMin * areaMin;
        delta      = 1.0e-16;

        localCoarsenCnt = 0;
        globalCoarsenCnt = 0;
/*
 *      Loop through all the nodes native to this domain looking for
 *      nodes that should be coarsened out.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
        
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
        
/*
 *          Check for various conditions that will exempt a node from removal:
 *
 *          1) does not have  2 arms or the segement is not a GB segement
 *          2) node is a 'fixed' node or surface node
 *          3) node is flagged as exempt from coarsen operations
 *          4) current domain does not 'own' at least one of the  segments
 *          5) If the node's arms are on different glide planes we might
 *             not allow the node to be removed.
 */
#if 1

#else
            if (node->numNbrs < 2 || node->numNbrs > 3) continue;

            breakLoop = 0;
		    if(param->polycrystal && IsOnGB(node)){
                
                for(j=0; breakLoop == 0 && j<node->numNbrs; j++){
                    nbr1 = GetNeighborNode(home, node, j);
                    if(nbr1 == (Node_t *)NULL){
                        printf("Warnig: neighbor not found. (at %s, L%d)\n",
                                __FILE__, __LINE__);
                    }
                    if(!IsOnGB(nbr1))continue;
                    if(!DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr1->myTag))continue;
                    dis = NodeDistance(home, node, nbr1);
                    if(fabs(dis) > 5)continue;

                    burg[0] = node->burgX[j];
                    burg[1] = node->burgY[j];
                    burg[2] = node->burgZ[j];

                    if(BurgType_FCC(burg) != BASIC_BURG){
                        continue;
                    }
                    if(SurfaceNode(node) && SurfaceNode(nbr1))continue;

                    planeType = GetPlaneType_FCC(home, node, nbr1);

                    if(!((planeType & GLIDE_PLANE) > 0 && (planeType & GB_PLANE) == 0))continue;

                    globalOp = (nbr1->myTag.domainID != thisDomain);
                    for (q = 0; q < nbr1->numNbrs; q++) {
                        globalOp |= (nbr1->nbrTag[q].domainID != thisDomain);
                    }

                    for (q = 0; q < node->numNbrs; q++) {
                        globalOp |= (node->nbrTag[q].domainID != thisDomain);
                    }

                    vec[0] = nbr1->x - node->x;
                    vec[1] = nbr1->y - node->y;
                    vec[2] = nbr1->z - node->z;
                    ZImage(home->param, vec, vec+1, vec+2);

                    newPos[0] = node->x + 0.5*vec[0];
                    newPos[1] = node->y + 0.5*vec[1];
                    newPos[2] = node->z + 0.5*vec[2];

                    if(BeforeMergeNodes(home, node, nbr1, newPos, &dis)){
                        if(fabs(dis) > 5.0){
                            PrintNode(node);
                            PrintNode(nbr1);
                            Fatal("Warning in %s at %d, dis %f\n", __FILE__, __LINE__, dis);
                        }
                    }else{
                        continue;
                    }

                    FoldBox(home->param, newPos, newPos+1, newPos+2);

                    MergeNode(home, OPCLASS_REMESH, node, nbr1, newPos,
                              &mergedNode, &mergeStatus, globalOp);
            
                    mergeDone = mergeStatus & MERGE_SUCCESS;

                    if(mergedNode != (Node_t *)NULL){
                        SplitIntersectNode(home, mergedNode, OPCLASS_REMESH, globalOp);
                        mergedNode->flags |= NO_MESH_COARSEN;
                    }

                    breakLoop = 1;
		            break;
                }
		    }
            
            if(breakLoop)continue;
	        if(node ==(Node_t *)NULL)continue;
#endif
	        if(node->numNbrs != 2 )continue;

            if (node->constraint == PINNED_NODE) continue;

            nbr1 = GetNeighborNode(home, node, 0);
            nbr2 = GetNeighborNode(home, node, 1);
            if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

#if CROSS_SLIP_GB_MODEL
#else
            if(nbr1 != (Node_t *)NULL && nbr2 != (Node_t *)NULL){
                if(NodeOnGB(home, nbr1) & INTERSECT_NODE &&
                   NodeOnGB(home, nbr2) & INTERSECT_NODE){
                    continue;
                }
            }
#endif
            if (node->flags & NO_MESH_COARSEN) {
                continue;
            }
        
            if (!DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr1->myTag) &&
                !DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr2->myTag)) {
                continue;
            }

            hasRemoteNbr  = (node->myTag.domainID != nbr1->myTag.domainID);
            hasRemoteNbr |= (node->myTag.domainID != nbr2->myTag.domainID);

            if(home->param->polycrystal){
                GetArmGrain(home, node, nbr1, armInGrain1, armOnGB1);
                GetArmGrain(home, node, nbr2, armInGrain2, armOnGB2);
                              
                if(armOnGB1 != armOnGB2)continue;
                if(armOnGB1 == -1 && armInGrain1 != armInGrain2)continue;
            }
/*
 *          Calculate the lengths of the node's 2 arms plus
 *          the distance between the two neighbor nodes.
 *
 *          If periodic boundaries are enabled, the nodes may
 *          be on opposite side of the problem space, so adjust
 *          the lengths/distances accordingly.
 */
            vec1x = nbr1->x - node->x;
            vec1y = nbr1->y - node->y;
            vec1z = nbr1->z - node->z;
              
            vec2x = nbr2->x - node->x;
            vec2y = nbr2->y - node->y;
            vec2z = nbr2->z - node->z;
        
            vec3x = vec2x - vec1x;
            vec3y = vec2y - vec1y;
            vec3z = vec2z - vec1z;
              
            ZImage(param, &vec1x, &vec1y, &vec1z);
            ZImage(param, &vec2x, &vec2y, &vec2z);
            ZImage(param, &vec3x, &vec3y, &vec3z);
        
            r1 = sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);
            r2 = sqrt(vec2x*vec2x + vec2y*vec2y + vec2z*vec2z);
            r3 = sqrt(vec3x*vec3x + vec3y*vec3y + vec3z*vec3z);
        
/*
 *          If we are enforcing the use of segment glide planes, we do not
 *          want to coarsen out a node whose arms are on different glide planes.
 *
 *          However, when glide planes constraints are being enforced, we
 *          tend to get a lot of debris (i.e. tiny triangles or quadrangles,
 *          really short segments, etc) that oscillate quickly and adversely
 *          affect the timestep.  So, under the following circumstances,
 *          we'll allow the code to violate glide plane constraints.
 *
 *          1) If the node has an attached segment less than 1b in length
 *          2) If the node has two segments less than 20% of the minimum
 *             segment length AND both of the node's neighbors are connected
 *             to each other (i.e. 3 nodes form a triangle)
 *          3) If the glide planes are allowed to be 'fuzzy' add
 *             a few extra exceptions (see below)
 */
            inhibitMerge = 0;
            violateGlidePlanesOK = 0;
            if (param->enforceGlidePlanes) {
                int   connectionID; 
                if ((r1 < 1.0) ||
                    (r2 < 1.0) ||
                    ((MAX(r1, r2) < MAX(1.0, 0.20 * param->minSeg) &&
                     (Connected(nbr1, nbr2, &connectionID))))) {
                    violateGlidePlanesOK = 1;
                }

                gp0[0] = node->nx[0];
                gp0[1] = node->ny[0];
                gp0[2] = node->nz[0];

                gp1[0] = node->nx[1];
                gp1[1] = node->ny[1];
                gp1[2] = node->nz[1];

/*
 *              Glide planes are being used, but if we are allowing
 *              some 'fuzziness' in the planes there are a couple
 *              situations in which we allow glide plane constraints
 *              to be violated.
 *
 *                1)  the two segment glide planes are within a small
 *                    number of degrees
 *                2)  If either segment is shorter than the annihilation
 *                    distance
 *                3)  when each segment is mapped to the closest
 *                    precise glide plane, if the two precise planes are
 *                    the same
 */
                if (param->allowFuzzyGlidePlanes) {
                    if (fabs(DotProduct(gp0, gp1)) > 0.9555) {
                        violateGlidePlanesOK = 1;
                    } else if ((r1 < param->rann) || (r2 <  param->rann)) {
                        violateGlidePlanesOK = 1;
                    } else {
                        real8 burg1[3], burg2[3];
                        real8 lineDir1[3], lineDir2[3];
                        real8 testPlane1[3], testPlane2[3];

                        burg1[0] = node->burgX[0];
                        burg1[1] = node->burgY[0];
                        burg1[2] = node->burgZ[0];

                        burg2[0] = node->burgX[1];
                        burg2[1] = node->burgY[1];
                        burg2[2] = node->burgZ[1];
                        lineDir1[0] = vec1x;
                        lineDir1[1] = vec1y;
                        lineDir1[2] = vec1z;

                        lineDir2[0] = vec2x;
                        lineDir2[1] = vec2y;
                        lineDir2[2] = vec2z;

/*
 *                      FindGlidePlanes() just uses l cross b if fuzzy
 *                      planes are allowed, so we temporarily reset the value
 *                      as a quick kludge to find the closest precise plane.
 */
                        param->allowFuzzyGlidePlanes = 0;
                        FindGlidePlane(home, burg1, lineDir1, testPlane1, node, nbr1);
                        FindGlidePlane(home, burg2, lineDir2, testPlane2, node, nbr2);
                        param->allowFuzzyGlidePlanes = 1;

                        if (fabs(DotProduct(testPlane1, testPlane2)) > 0.99) {
                            violateGlidePlanesOK = 1;
                        }
                    }
                }

                if (!violateGlidePlanesOK) {
                    cross(gp0, gp1, tmp3);
                    if (fabs(DotProduct(tmp3, tmp3)) > 1.0e-3) {
                        continue;
                    }
                }
            }

            if(SurfaceNode(nbr1) && SurfaceNode(node) ||
               SurfaceNode(nbr2) && SurfaceNode(node))continue;

            if(IsOnGB(node)){
                if(nbr1->onGB != nbr2->onGB)continue;
                if(nbr1->inGrain != nbr2->inGrain)continue;
                if((NodeOnGB(home, node) & INTERSECT_NODE) == 0){
                    if((NodeOnGB(home, nbr1) & INTERSECT_NODE) == 1 &&
                       (NodeOnGB(home, nbr2) & INTERSECT_NODE) == 1)continue;
                }
            }else{
                if(IsOnGB(nbr1) && IsOnGB(nbr2))continue;
            }
/*
 *          If coarsening out a node would leave a segment longer
 *          than a defined length, the node should not be removed.
 *          This 'cutoff length' is the maximum segment length
 *          if all involved nodes are within the same domain; but
 *          if a remote node is involved, we need to set the
 *          cutoff length to (at most) 1/2 the cell length.  This
 *          is needed because the remote node may potentially be involved
 *          in a simultaneous mesh coarsening in the remote domain,
 *          and although the node would not be removed, it could
 *          be repositioned resulting in a segment that spanned
 *          more than 2 cells... this is a bad thing.
 */
            if (((hasRemoteNbr == 0) && (r3 > cutoffLength1)) ||
                ((hasRemoteNbr == 1) && (r3 > cutoffLength2))) {
                continue;
            }

/*
 *          Check if the area of the triangle defined by node
 *          and its two neighbors, plus determine if that area
 *          is increasing or decreasing.
 */

            s = 0.5 * (r1 + r2 + r3);
            area2 = (s * (s-r1) * (s-r2) * (s-r3));
        
            dvec1xdt = nbr1->vX - node->vX;
            dvec1ydt = nbr1->vY - node->vY;
            dvec1zdt = nbr1->vZ - node->vZ;
              
            dvec2xdt = nbr2->vX - node->vX;
            dvec2ydt = nbr2->vY - node->vY;
            dvec2zdt = nbr2->vZ - node->vZ;
              
            dvec3xdt = dvec2xdt - dvec1xdt;
            dvec3ydt = dvec2ydt - dvec1ydt;
            dvec3zdt = dvec2zdt - dvec1zdt;
              
            dr1dt = ((vec1x * dvec1xdt) + (vec1y * dvec1ydt) +
                     (vec1z * dvec1zdt)) / (r1 + delta);
        
            dr2dt = ((vec2x * dvec2xdt) + (vec2y * dvec2ydt) +
                     (vec2z * dvec2zdt)) / (r2 + delta);
        
            dr3dt = ((vec3x * dvec3xdt) + (vec3y * dvec3ydt) +
                     (vec3z * dvec3zdt)) / (r3 + delta);
        
        
            dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);
        
            darea2dt = (dsdt * (s-r1) * (s-r2) * (s-r3));
            darea2dt += s * (dsdt-dr1dt) * (s-r2) * (s-r3);
            darea2dt += s * (s-r1) * (dsdt-dr2dt) * (s-r3);
            darea2dt += s * (s-r1) * (s-r2) * (dsdt-dr3dt);
        
/*
 *          If the area is less than the specified minimum and shrinking,
 *          or one of the arms is less than the minimum segment length, the
 *          node should be removed.
 */

            if (((area2 < areaMin2) && (darea2dt < 0.0)) ||
                ((r1 < param->minSeg) || (r2 < param->minSeg))) {
        
                EstCoarsenForces(home, nbr1, node, nbr2, f0seg1, f1seg1);
        
                mergeDone = 0;

                nbr1Tag = nbr1->myTag;
                nbr2Tag = nbr2->myTag;
        
/*
 *              If either of the neighbor nodes (or any of their neighbors)
 *              is in a remote domain, the operation must be treated as global.
 *              This is necessary to prevent an inconsistent linkage problem.
 *              For example, given the segments A--B--C--D where nodes A, B
 *              and C are in domain 1, D is in domain 2, and segment C--D is
 *              owned by domain 1:  domain 1 could coarsen node B into A, then
 *              node C into D.  If the first operation was not communicated
 *              to the remote domain, an inconsitency would arise.
 *
 *              NOTE:  It is safe to not distribute the purely local coarsen
 *                     operations so long as no other topological operations
 *                     are done after Remesh() but before the ghost node
 *                     are redistributed.
 */
                globalOp = ((nbr1->myTag.domainID != thisDomain) ||
                            (nbr2->myTag.domainID != thisDomain));

                for (q = 0; q < nbr1->numNbrs; q++) {
                    globalOp |= (nbr1->nbrTag[q].domainID != thisDomain);
                }
        
                for (q = 0; q < nbr2->numNbrs; q++) {
                    globalOp |= (nbr2->nbrTag[q].domainID != thisDomain);
                }
        
                oldTag1 = nbr1->myTag;
                oldTag2 = node->myTag;
                oldTag3 = nbr2->myTag;
/*
 *              If the first neighbor is not exempt from a coarsen
 *              operation, attempt to merge the nodes.
 */
                if ((nbr1->flags & NO_MESH_COARSEN) == 0) {
                    if(violateGlidePlanesOK == 1 && r1 < 1){
                        newPos[0] = nbr1->x;
                        newPos[1] = nbr1->y;
                        newPos[2] = nbr1->z;
                    }else{
                        vec[0] = (node->x - nbr1->x);
                        vec[1] = (node->y - nbr1->y);
                        vec[2] = (node->z - nbr1->z);
                        ZImage(param, vec, vec+1, vec+2);
                        VECTOR_COPY(vec1, vec);
                        
                        if(AdjustNodeDirection(home, nbr1, vec) == 1){
                            VECTOR_COPY(vec, vec1);
                            newPos[0] = nbr1->x + 0.5*vec[0];
                            newPos[1] = nbr1->y + 0.5*vec[1];
                            newPos[2] = nbr1->z + 0.5*vec[2];
                        }else{
                            newPos[0] = nbr1->x;
                            newPos[1] = nbr1->y;
                            newPos[2] = nbr1->z;
                        }
                    }

                    FoldBox(home->param,newPos, newPos+1, newPos+2);
                    MergeNode(home, OPCLASS_REMESH, node, nbr1, newPos,
                              &mergedNode, &mergeStatus, globalOp);
                    mergeDone = mergeStatus & MERGE_SUCCESS;
                }
/*
 *              If the merge could not be done, try using
 *              the other neighbor.
 */
                if (mergeDone == 0) {
                    if ((nbr2->flags & NO_MESH_COARSEN) == 0) {
                        if(violateGlidePlanesOK == 1 && r2 < 1){
                            newPos[0] = nbr2->x;
                            newPos[1] = nbr2->y;
                            newPos[2] = nbr2->z;
                        }else{
                            vec[0] = (node->x - nbr2->x);
                            vec[1] = (node->y - nbr2->y);
                            vec[2] = (node->z - nbr2->z);
                            ZImage(param, vec, vec+1, vec+2);
                            VECTOR_COPY(vec1, vec);
                         
                            if(AdjustNodeDirection(home, nbr2, vec) == 1){
                                VECTOR_COPY(vec, vec1);
                                newPos[0] = nbr2->x + 0.5*vec[0];
                                newPos[1] = nbr2->y + 0.5*vec[1];
                                newPos[2] = nbr2->z + 0.5*vec[2];
                            }else{
                                newPos[0] = nbr2->x;
                                newPos[1] = nbr2->y;
                                newPos[2] = nbr2->z;
                            }
                        }
    
                        FoldBox(home->param,newPos, newPos+1, newPos+2);
                        MergeNode(home, OPCLASS_REMESH, node, nbr2, newPos,
                                  &mergedNode, &mergeStatus, globalOp);
                        mergeDone = mergeStatus & MERGE_SUCCESS;
                    }
                }

/*
 *              If the merge was successful, update the forces
 *              on the remaining nodes.   Otherwise go back and
 *              continue looking for more nodes to coarsen out.
 */
				
                if (mergeDone == 0) continue;
        
                localCoarsenCnt++;

#ifdef DEBUG_TOPOLOGY_CHANGES
                if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                    printf("Coarsen: (%d,%d)--(%d,%d)--(%d,%d)\n",
                           oldTag1.domainID, oldTag1.index,
                           oldTag2.domainID, oldTag2.index,
                           oldTag3.domainID, oldTag3.index);
                }
#endif
                nbr1 = GetNodeFromTag(home, nbr1Tag);
                nbr2 = GetNodeFromTag(home, nbr2Tag);
        
                if ((nbr1 == (Node_t *)NULL) && (nbr2 == (Node_t *)NULL)) {
                    continue;
                }

/*
 *              The merge will have placed the resultant node at the
 *              location of either nbr1 or nbr2, but the merge function
 *              determines which of the two specified nodes is deleted
 *              and which survives, which means that nbr1 or nbr2 may
 *              have been the deleted and the <node> repositioned to
 *              the correct location... so if one of the nbr nodes does
 *              not exist anymore, <mergedNode> (if it exists) should
 *              be the node that replaced the nbr.
 */
                if (nbr1 == (Node_t *)NULL) {
                    nbr1 = mergedNode;
                } else if (nbr2 == (Node_t *)NULL) {
                    nbr2 = mergedNode;
                }

/*
 *              At this point, if we don't have a node at the location
 *              of at least one of the original nbr nodes, looks like
 *              some nodes were orphaned and deleted, so the force
 *              estimates we made are not applicable.
 */
                if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                    continue;
                }
        
                mergedNode->flags |= NO_MESH_COARSEN;
/*
 *              Reset force/velocity for the two remaining nodes
 *              and mark the forces for those nodes and all their
 *              neighbors as obsolete.  This is done because the
 *              force estimates above are good enough for
 *              the remainder of this timestep, but we need to 
 *              recalculate more exact forces for these nodes
 *              before (or at the beginning of) the next timestep.
 */
                ResetSegForces(home, nbr1, &nbr2->myTag, f0seg1[0],
                               f0seg1[1], f0seg1[2], 1);

                ResetSegForces(home, nbr2, &nbr1->myTag, f1seg1[0],
                               f1seg1[1], f1seg1[2], 1);

/*
 *              Originally, we were resetting the velocity of the
 *              neighboring nodes here, but it appeared we were
 *              better off NOT doing that, so we've dropped that.
 *              for now...
 */
                MarkNodeForceObsolete(home, nbr1);
                MarkNodeForceObsolete(home, nbr2);
        
                for (q = 0; q < nbr1->numNbrs; q++) {
                    nbr = GetNodeFromTag(home, nbr1->nbrTag[q]);
                    if (nbr == (Node_t *)NULL) continue;
                    MarkNodeForceObsolete(home, nbr);
                }
        
                for (q = 0; q < nbr2->numNbrs; q++) {
                    nbr = GetNodeFromTag(home, nbr2->nbrTag[q]);
                    if (nbr == (Node_t *)NULL) continue;
                    MarkNodeForceObsolete(home, nbr);
                }
		
/*
 *              If we are enforcing glide planes but allowing them to be
 *              slightly fuzzy, we need to recalculate the glide plane for
 *              the new segment.
 */
                if (param->enforceGlidePlanes &&
                    param->allowFuzzyGlidePlanes) {

                    RecalcSegGlidePlane(home, nbr1, nbr2, 1);
                }
            }
        }
        
#ifdef DEBUG_LOG_MESH_COARSEN
#ifdef PARALLEL
        MPI_Reduce(&localCoarsenCnt, &globalCoarsenCnt, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
        globalCoarsenCnt = localCoarsenCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Remesh: coarsen count = %d\n", globalCoarsenCnt);
        }
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    TrySegBisect
 *      Description: Attempt to bisect a discretization segment during
 *                   mesh refinement by splitting the specified node.
 *                   Some restrictions apply and if any of the criteria
 *                   are not met, the segment will not be bisected.
 *
 *      Arguments:
 *          origNode  Pointer to location containing the pointer to
 *                    the node that we are considering splitting.  
 *                    On return to the caller this will contain a 
 *                    pointer to the node left at the position
 *                    <origNode> was in on entry to this function.
 *          nbr       Pointer to node at the end of the segment to
 *                    be bisected.
 *          armID     Index of the segment in <origNode>'s arm list
 *          segLen    Length (in units of b) of the segment to be
 *                    bisected.
 *          vec       Vector from <origNode> to <nbr>
 *          area2     area (squared) enclosed by <origNode>, <nbr> and
 *                    the other neighbor of <origNode> which is not 
 *                    provided to this function.
 *          areaMax2  Maximum area (squared) needed before the segment
 *                    may be bisected (unless segLen has exceeded the
 *                    maximum permitted segment length).
 *          darea2dt  Indicates if the area indicated by <area2> is
 *                    increasing or decreasing in size.
 *          splitIsOK Additional flag indicating if the segment is
 *                    permitted to be bisected.  (Caller may have 
 *                    determined the bisection should not occur).
 *                    0 == not permitted, 1 == permitted.
 *          didBisect Flag indicating if the other segment attached
 *                    to <origNode> was bisected this cycle.  Mostly for
 *                    possible future use...  0 == not bisected,
 *                    1 == bisected.
 *
 *      Returns:  1 if the segment was bisected, 0 if not
 *          
 *------------------------------------------------------------------------*/
static int TrySegBisect(Home_t *home, Node_t **origNode, Node_t *nbr,
                        int armID, real8 segLen, real8 vec[3], real8 area2,
                        real8 areaMax2, real8 darea2dt, int splitIsOK,
                        int didBisect)
{
        int     splitStatus, thisDomain, globalOp;
        int     *armList, armCount;
        real8   newVel[3], newPos[3], nodeVel[3], nodePos[3],posInGrain;
        real8   f0seg1[3], f0seg2[3], f1seg1[3], f1seg2[3], newInGrain;
        Tag_t   oldTag1, oldTag2;
        Node_t  *node, *splitNode1, *splitNode2;
        Param_t *param;

        thisDomain = home->myDomain;
        param = home->param;
        armList = &armID;
        armCount = 1;
        didBisect = 0;
        node = *origNode;

#if 0
/*
 *      The idea here was that if we already bisected 1 segment of
 *      the 2-node, we'd only bisect the second segment this cycle
 *      if the second segment exceeded the max seg length. I'm not
 *      sure yet if that is the best thing to do, so for now we
 *      won't prevent the second segment from being bisected because
 *      the first was.
 */
        if (!splitIsOK || (didBisect && (segLen < param->maxSeg))) {
            return(0);
        }
#else
        if (!splitIsOK) {
            return(0);
        }
#endif

/*
 *      If the segment is over the max segment length, or the area
 *      of the triangle defined by the node and its two neighbors
 *      is above the limit AND the area is increasing in size AND
 *      the segment is not too small to bisect, we cut it.
 */

		
        if ((segLen > param->maxSeg) ||
            ((area2 > areaMax2) &&
             (segLen >= param->minSeg * 2.0) &&
             (darea2dt >= 0.0))) {
/*
 *          When bisecting a segment, we always move exactly one arm
 *          from the original node to the node being created... in this
 *          case, the arm with index <armID>
 */
            newVel[0] = (node->vX + nbr->vX) * 0.5;
            newVel[1] = (node->vY + nbr->vY) * 0.5;
            newVel[2] = (node->vZ + nbr->vZ) * 0.5;

            newPos[0] = node->x + (vec[0] * 0.5);
            newPos[1] = node->y + (vec[1] * 0.5);
            newPos[2] = node->z + (vec[2] * 0.5);

            EstRefinementForces(home, node, nbr, newPos, vec,
                                f0seg1, f1seg1, f0seg2, f1seg2);

            FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);

/*
 *          This should be a global operation distributed
 *          out to remote domains only if the neighbor node
 *          is in another domain.
 */
            globalOp = (nbr->myTag.domainID != node->myTag.domainID);

            nodePos[0] = node->x;
            nodePos[1] = node->y;
            nodePos[2] = node->z;

            nodeVel[0] = node->vX;
            nodeVel[1] = node->vY;
            nodeVel[2] = node->vZ;

            oldTag1 = node->myTag;
            oldTag2 = nbr->myTag;

			
            splitStatus = SplitNode(home, OPCLASS_REMESH,
                                    node, nodePos, newPos,
                                    nodeVel, newVel,
                                    armCount, armList,
                                    globalOp, &splitNode1,
                                    &splitNode2, 0);

            if (splitStatus == SPLIT_SUCCESS) {

                didBisect = 1;
                *origNode = splitNode1;
/*
 *              The force estimates above are good enough for the
 *              remainder of this timestep, but mark the force and
 *              velocity data for some nodes as obsolete so that
 *              more accurate forces will be recalculated either at
 *              the end of this timestep, or the beginning of the next.
 */
                MarkNodeForceObsolete(home, splitNode1);
                MarkNodeForceObsolete(home, splitNode2);
                MarkNodeForceObsolete(home, nbr);

/*
 *              Reset nodal forces on all nodes involved
 *              in the split.
 */
                ResetSegForces(home, splitNode1, &splitNode2->myTag,
                               f0seg1[0], f0seg1[1], f0seg1[2], 1);

                ResetSegForces(home, splitNode2, &splitNode1->myTag,
                               f1seg1[0], f1seg1[1], f1seg1[2], 1);

                ResetSegForces(home, splitNode2, &nbr->myTag,
                               f0seg2[0], f0seg2[1], f0seg2[2], 1);

                ResetSegForces(home, nbr, &splitNode2->myTag,
                               f1seg2[0], f1seg2[1], f1seg2[2], 1);

                (void)EvaluateMobility(home, splitNode1);
                (void)EvaluateMobility(home, splitNode2);
                (void)EvaluateMobility(home, nbr);

/*
 *              When debugging, dump some info on
 *              topological changes taking place and
 *              the nodes involved
 */
#if 0
				if(globalOp == 1){
					if ((dbgDom < 0) || (dbgDom == thisDomain)) {
						printf("Remesh/refine1:  (%d,%d)--(%d,%d) ==> "
							   "(%d,%d)--(%d,%d)--(%d,%d)\n",
							   oldTag1.domainID, oldTag1.index,
							   oldTag2.domainID, oldTag2.index,
							   splitNode1->myTag.domainID,
							   splitNode1->myTag.index,
							   splitNode2->myTag.domainID,
							   splitNode2->myTag.index,
							   nbr->myTag.domainID,
							   nbr->myTag.index);
//						PrintNode(splitNode1);
//						PrintNode(splitNode2);
//						PrintNode(nbr);
					}
				}
#endif
            }  /* if (splitStatus == SPLIT_SUCCESS) */
        }

        return(didBisect);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    MeshRefine
 *      Description: 
 *
 *------------------------------------------------------------------------*/
static void MeshRefine(Home_t *home)
{
        int     i, thisDomain, didBisect, seg, splitIsOK;
        int     splitStatus, armIndex, armCount = 1, globalOp;
        int     armID, *armList;
        int     localRefineCnt, globalRefineCnt;
        int     splitOK[2], splitSegList[2];
        real8   areaMax, areaMax2, maxSeg2, newInGrain;
        real8   delta, r1, r2, r3, s, area2, segLen;
        real8   dvec1xdt, dvec1ydt, dvec1zdt,posInGrain;
        real8   dvec2xdt, dvec2ydt, dvec2zdt;
        real8   dvec3xdt, dvec3ydt, dvec3zdt;
        real8   dr1dt, dr2dt, dr3dt, dsdt, darea2dt;
        real8   *vec, vec1[3], vec2[3], vec3[3];
        real8   nodePos[3], nodeVel[3], newVel[3], newPos[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        real8   dis, burg[3];
        int     planeType;
        Tag_t   oldTag1, oldTag2;
        Node_t  *node, *nbr, *nbr1, *nbr2, *splitNode1, *splitNode2;
        Param_t *param;

/*
 *      Initialize some of the constants we'll need during this call.
 */
        thisDomain = home->myDomain;
        param      = home->param;

        areaMax  = param->remeshAreaMax;
        areaMax2 = areaMax * areaMax;
        maxSeg2  = param->maxSeg * param->maxSeg;
        delta    = 1.0e-16;

        localRefineCnt = 0;
        globalRefineCnt = 0;
        
/*
 *      Loop through all the  native nodes looking for segments
 *      that need to be refined.
 */
        armList = &armID;
        
        for (i = 0; i < home->newNodeKeyPtr; i++) {
        
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
	        if (SurfaceNode(node))continue;
        
/*
 *          If the node has only two arms, use both arm lengths and
 *          the area in the triangle defined by the node and its 
 *          neighbors to decide whether the node should be split or not.
 */
            if (node->numNbrs == 2) {
        
/*
 *              Calculate the lengths of the node's 2 arms plus
 *              the distance between the two neighbor nodes.
 */
                nbr1 = GetNeighborNode(home, node, 0);
                nbr2 = GetNeighborNode(home, node, 1);
        
                if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                vec1[0] = nbr1->x - node->x;
                vec1[1] = nbr1->y - node->y;
                vec1[2] = nbr1->z - node->z;
              
                vec2[0] = nbr2->x - node->x;
                vec2[1] = nbr2->y - node->y;
                vec2[2] = nbr2->z - node->z;
        
                vec3[0] = vec2[0] - vec1[0];
                vec3[1] = vec2[1] - vec1[1];
                vec3[2] = vec2[2] - vec1[2];
        
                ZImage(param, &vec1[0], &vec1[1], &vec1[2]);
                ZImage(param, &vec2[0], &vec2[1], &vec2[2]);
                ZImage(param, &vec3[0], &vec3[1], &vec3[2]);
        
                r1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
                r2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
                r3 = sqrt(vec3[0]*vec3[0] + vec3[1]*vec3[1] + vec3[2]*vec3[2]);
        
                s = 0.5 * (r1 + r2 + r3);
                area2 = (s * (s - r1) * (s - r2) * (s - r3));
        
/*
 *              Determine if the area of the triangle defined by the node
 *              and its two neighbors is increasing or decreasing.
 */
                dvec1xdt = nbr1->vX - node->vX;
                dvec1ydt = nbr1->vY - node->vY;
                dvec1zdt = nbr1->vZ - node->vZ;

                dvec2xdt = nbr2->vX - node->vX;
                dvec2ydt = nbr2->vY - node->vY;
                dvec2zdt = nbr2->vZ - node->vZ;

                dvec3xdt = dvec2xdt - dvec1xdt;
                dvec3ydt = dvec2ydt - dvec1ydt;
                dvec3zdt = dvec2zdt - dvec1zdt;

                dr1dt = ((vec1[0] * dvec1xdt) + (vec1[1] * dvec1ydt) +
                         (vec1[2] * dvec1zdt)) / (r1 + delta);

                dr2dt = ((vec2[0] * dvec2xdt) + (vec2[1] * dvec2ydt) +
                         (vec2[2] * dvec2zdt)) / (r2 + delta);

                dr3dt = ((vec3[0] * dvec3xdt) + (vec3[1] * dvec3ydt) +
                         (vec3[2] * dvec3zdt)) / (r3 + delta);

                dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);

                darea2dt = (dsdt * (s-r1) * (s-r2) * (s-r3));
                darea2dt += s * (dsdt-dr1dt) * (s-r2) * (s-r3);
                darea2dt += s * (s-r1) * (dsdt-dr2dt) * (s-r3);
                darea2dt += s * (s-r1) * (s-r2) * (dsdt-dr3dt);

/*
 *              Check if the current domain owns the segments.
 *              It may only split a segment it owns...
 */
                splitOK[0] = DomainOwnsSeg(home, OPCLASS_REMESH,
                                           thisDomain, &nbr1->myTag);

                splitOK[1] = DomainOwnsSeg(home, OPCLASS_REMESH,
                                           thisDomain, &nbr2->myTag);
/*
 *              If both nodes of the segment are flagged to have forces
 *              updated, forces and velocities on the node may not be
 *              good enough to accurately position the new node, so
 *              don't split this segment this cycle.
 */
                if (((node->flags & NODE_RESET_FORCES) != 0) &&
                    ((nbr1->flags & NODE_RESET_FORCES) != 0)) {
                    splitOK[0] = 0;
                }

                if (((node->flags & NODE_RESET_FORCES) != 0) &&
                    ((nbr2->flags & NODE_RESET_FORCES) != 0)) {
                    splitOK[1] = 0;
                }

/*
 *              It would be preferable to bisect the longest segment
 *              first, so check the segment lengths and set the
 *              preferred order.
 */
                if (r1 > r2) {
                    splitSegList[0] = 1;
                    splitSegList[1] = 2;
                } else {
                    splitSegList[0] = 2;
                    splitSegList[1] = 1;
                }
/*
 *              Loop over both segments in the preferred order and try to 
 *              bisect them.
 *
 *              Note: If we bisect a segment on the first loop iteration,
 *              the following segment will only be bisected if it exceeds
 *              the maximum seg length.
 */
                didBisect = 0;

                for (seg = 0; seg < 2; seg++) {
                    if (splitSegList[seg] == 1) {
                        nbr = nbr1;
                        vec = vec1;
                        segLen = r1;
                        splitIsOK = splitOK[0];
                    } else {
                        nbr = nbr2;
                        vec = vec2;
                        segLen = r2;
                        splitIsOK = splitOK[1];
                    }

                    armID = GetArmID(home, node, nbr);

			        if(home->param->polycrystal && IsOnGB(node) && IsOnGB(nbr)){

                        dis = NodeDistance(home, node, nbr);

                        burg[0] = node->burgX[armID];
                        burg[1] = node->burgY[armID];
                        burg[2] = node->burgZ[armID];

                        planeType = GetPlaneType_FCC(home, node, nbr);

                        if((planeType & GLIDE_PLANE) > 0 && (planeType & GB_PLANE) == 0 &&
                           fabs(dis) < 5 && BurgType_FCC(burg) == BASIC_BURG)continue;
				
			        }

                    didBisect = TrySegBisect(home, &node, nbr, armID, segLen,
                                             vec, area2, areaMax2, darea2dt,
                                             splitIsOK, didBisect);

                    localRefineCnt += didBisect;

                }  /* for (seg = 0; seg < 2; ...) */

            } else {
/*
 *              For nodes with other than exactly two arms, we
 *              just bisect any arm exceeding the max allowable
 *              segment length, but only do the bisection from the
 *              domain "owning" the segment.
 *
 *              A "split" is only considered a global operation
 *              to be sent to remote domains if the segment being
 *              split spans domains.
 */
                for (armIndex = 0; armIndex < node->numNbrs; ) {
        
                    nbr1 = GetNeighborNode(home, node, armIndex);
        
                    if (nbr1 == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        armIndex++;
                        continue;
                    }

                    if (!DomainOwnsSeg(home, OPCLASS_REMESH,
                                       thisDomain, &nbr1->myTag)) {
                        armIndex++;
                        continue;
                    }
        
                    vec1[0] = nbr1->x - node->x;
                    vec1[1] = nbr1->y - node->y;
                    vec1[2] = nbr1->z - node->z;
        
                    ZImage(param, &vec1[0], &vec1[1], &vec1[2]);
        
                    r1 = vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2];
        
                    if (r1 > maxSeg2) {
        
                        newPos[0] = node->x + (vec1[0] * 0.5);
                        newPos[1] = node->y + (vec1[1] * 0.5);
                        newPos[2] = node->z + (vec1[2] * 0.5);
        
                        newVel[0] = (node->vX + nbr1->vX) * 0.5;
                        newVel[1] = (node->vY + nbr1->vY) * 0.5;
                        newVel[2] = (node->vZ + nbr1->vZ) * 0.5;
        
                        EstRefinementForces(home, node, nbr1, newPos, vec1,
                                            f0seg1, f1seg1, f0seg2, f1seg2);

                        FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);

/*
 *                      When bisecting a segment, we always move
 *                      exactly one arm from the original node
 *                      to the node being created...
 */
                        *armList = armIndex;
        
/*
 *                      This should be a global operation
 *                      distributed out to remote domains only
 *                      if the segment spans domains.
 */
                        globalOp = (nbr1->myTag.domainID !=
                                    node->myTag.domainID);
        
                        nodePos[0] = node->x;
                        nodePos[1] = node->y;
                        nodePos[2] = node->z;
        
                        nodeVel[0] = node->vX;
                        nodeVel[1] = node->vY;
                        nodeVel[2] = node->vZ;
        
                        oldTag1 = node->myTag;
                        oldTag2 = nbr1->myTag;
						
			            if(home->param->polycrystal && IsOnGB(node) && IsOnGB(nbr1)){

                            dis = NodeDistance(home, node, nbr1);

                            burg[0] = node->burgX[armIndex];
                            burg[1] = node->burgY[armIndex];
                            burg[2] = node->burgZ[armIndex];

                            planeType = GetPlaneType_FCC(home, node, nbr1);

                            if((planeType & GLIDE_PLANE) > 0 && (planeType & GB_PLANE) == 0 &&
                               fabs(dis) < 5 && BurgType_FCC(burg) == BASIC_BURG){
                                armIndex++;
                                continue;
                            }
				
			            }

                        splitStatus = SplitNode(home, OPCLASS_REMESH,
                                                node, nodePos, newPos,
                                                nodeVel, newVel,
                                                armCount, armList,
                                                globalOp, &splitNode1,
                                                &splitNode2, 0);
        
                        if (splitStatus == SPLIT_SUCCESS) {

                            localRefineCnt++;

                            node = splitNode1;
/*
 *                          The force estimates above are good enough for
 *                          the remainder of this timestep, but mark the
 *                          force and velocity data for some nodes as
 *                          obsolete so that more accurate forces will be
 *                          recalculated either at the end of this timestep,
 *                          or the beginning of the next.
 */
                            MarkNodeForceObsolete(home, splitNode1);
                            MarkNodeForceObsolete(home, splitNode1);
                            MarkNodeForceObsolete(home, nbr1);
        
/*
 *                          Reset nodal forces on all nodes involved
 *                          in the split.
 */
                            ResetSegForces(home, splitNode1, &splitNode2->myTag,
                                           f0seg1[0], f0seg1[1], f0seg1[2], 1);
        
                            ResetSegForces(home, splitNode2, &splitNode1->myTag,
                                           f1seg1[0], f1seg1[1], f1seg1[2], 1);
        
                            ResetSegForces(home, splitNode2, &nbr1->myTag,
                                           f0seg2[0], f0seg2[1], f0seg2[2], 1);
        
                            ResetSegForces(home, nbr1, &splitNode2->myTag,
                                           f1seg2[0], f1seg2[1], f1seg2[2], 1);
        
                            (void)EvaluateMobility(home, splitNode1);
                            (void)EvaluateMobility(home, splitNode2);
                            (void)EvaluateMobility(home, nbr1);
        
/*
 *                          When debugging, dump some info on
 *                          topological changes taking place and
 *                          the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                                printf("Remesh/refine3:  (%d,%d)--(%d,%d) ==> "
                                       "(%d,%d)--(%d,%d)--(%d,%d)\n",
                                       oldTag1.domainID, oldTag1.index,
                                       oldTag2.domainID, oldTag2.index,
                                       splitNode1->myTag.domainID,
                                       splitNode1->myTag.index,
                                       splitNode2->myTag.domainID,
                                       splitNode2->myTag.index,
                                       nbr1->myTag.domainID,
                                       nbr1->myTag.index);
                                PrintNode(splitNode1);
                                PrintNode(splitNode2);
                                PrintNode(nbr1);
                            }
#endif
                        }
                    } else armIndex++;
                }
            }
        }
        
#ifdef DEBUG_LOG_MESH_REFINE
#ifdef PARALLEL
        MPI_Reduce(&localRefineCnt, &globalRefineCnt, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
        globalRefineCnt = localRefineCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Remesh: refine count = %d\n", globalRefineCnt);
        }
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RemeshRule_2
 *      Description: Base function which invokes all subroutines
 *                   needed for handling mesh operations specific
 *                   to the second remesh rule
 *
 *------------------------------------------------------------------------*/
void RemeshRule_2(Home_t *home)
{

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        MeshCoarsen(home);
        MeshRefine(home);
		

        return;
}
