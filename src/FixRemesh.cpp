/****************************************************************************
 *
 *	Function:	FixRemesh
 *	Description:	This function is called by each domain after
 *			a set of topology changes has been done (i.e.
 *			remesh, collision handling or SplitMultiNodes()).
 *			The function uses the list of operations provided
 *			by each domain to update its local topology with
 *			the necessary changes.
 *
 *			The opList received from each neighbor has been
 *			stored in remDom->inBuf, with length in
 *			remDom->inBufLen, by CommSendRemsh.
 *
 *****************************************************************************/

#include "Home.h"
#include "Util.h"


void FixRemesh(Home_t *home, const char flag[256])
{
	int			    idom, domIdx, rcvOpCount, opidx, i;
	real8			newPlane[3];
    real8   		cellSize[3], cellIndex[3];
	RemoteDomain_t	*remDom;
	Operate_t		*rcvOpList, *op;
	Node_t			*node1, *node2;
	Tag_t			tag1, tag2, tag3;
    Param_t     	*param;

#if 0
	Node_t 			*node3, *nbr;
#endif


    param = home->param;

#if CHECK_FIXREMESH
	printf("FR %s: Task %d, cycle=%d\n", flag, home->myDomain, home->cycle);

#endif
/*
 *	Loop through the neighbor domains
 */
	for (idom = 0; idom < home->remoteDomainCount; idom++) {

		domIdx = home->remoteDomains[idom];
		remDom = home->remoteDomainKeys[domIdx];

		rcvOpList = (Operate_t *) remDom->inBuf;
		rcvOpCount = remDom->inBufLen / sizeof(Operate_t);

/*
 *		Loop through the operations from this neighbor
 */
		for (opidx = 0; opidx < rcvOpCount; opidx++) {

			op = rcvOpList + opidx;
			switch(op->type) {

/*
 *			Change the node1->node2 link so it becomes
 *			a node1->node3 link instead.
 */
					
			
			case CHANGE_CONNECTION:
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
				
#if CHECK_FIXREMESH
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("CC (%d,%d), (%d,%d), (%d,%d). ",op->dom1,op->idx1,
					  op->dom2,op->idx2,op->dom3,op->idx3);

				node2 = GetNodeFromIndex(home, op->dom2,op->idx2);
				node3 = GetNodeFromIndex(home, op->dom3,op->idx3);
				PrintTag(node1);
				PrintTag(node2);
				PrintTag(node3);
				printf("\n");
#endif				
				if (node1 == NULL) break;
				
				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;
				tag3.domainID = op->dom3;
				tag3.index    = op->idx3;

				ChangeConnection (home, node1, &tag2, &tag3, 0);
				break;

/*
 *			Add an arm from node1 to node2
 */
			case INSERT_ARM:
				
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
				
#if CHECK_FIXREMESH		
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("IA (%d,%d), (%d,%d). ",op->dom1,op->idx1,
					  op->dom2,op->idx2);

				node2 = GetNodeFromIndex(home, op->dom2,
					                 op->idx2);
				PrintTag(node1);
				PrintTag(node2);
				printf("\n");
#endif
				if (node1 == NULL) break;
				
				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				InsertArm (home, node1, &tag2, op->bx, op->by,
					   op->bz, op->nx, op->ny, op->nz, 0);
				break;

/*
 *			Create a ghost node with the necessary domain
 *			and index, and set the coordinates and velocity
 *			appropriately.  Note: the node initially has
 *			no arms, they will be added later.  
 */
			case SPLIT_NODE:

#if 0
/*
 *                              Temporarily remove restriction that
 *                              domain must have knowledge of node
 *                              being split... if this does not cause
 *                              problems, remve it completely...
 */
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
				if (node1 == (Node_t *)NULL) break;
#endif
				
				node2 = GetNewGhostNode(home, op->dom2,
							op->idx2);
				FreeNodeArms(node2);
				
				node2->constraint = op->dom3;

				node2->x = op->x;
				node2->y = op->y;
				node2->z = op->z;

				node2->oldx = op->bx;
				node2->oldy = op->by;
				node2->oldz = op->bz;
				
#if CHECK_FIXREMESH
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("SN (%d,%d), (%d,%d) ",op->dom1,op->idx1,
					  op->dom2,op->idx2);

				node1 = GetNodeFromIndex(home, op->dom1,
	                                op->idx1);
				PrintTag(node1);
				PrintTag(node2);
				printf("\n");
#endif
/*
 *				For the node split operation, the glide plane
 *				values weren't required, but the node velocity
 *				was, so the glide plane values for the
 *				operation were comandeered to hold the
 *				velocity data.
 */
				node2->vX = op->nx;
				node2->vY = op->ny;
				node2->vZ = op->nz;

/*
 *              Insure the node is assigned to the proper cell
 */
                cellSize[0] = param->Lx / param->nXcells;
                cellSize[1] = param->Ly / param->nYcells;
                cellSize[2] = param->Lz / param->nZcells;

                cellIndex[0] = (int)(node2->x-param->minSideX) /
                                (int)cellSize[0];
                cellIndex[1] = (int)(node2->y-param->minSideY) /
                                (int)cellSize[1];
                cellIndex[2] = (int)(node2->z-param->minSideZ) /
                                (int)cellSize[2];

                cellIndex[0] = MAX(0, cellIndex[0]);
                cellIndex[1] = MAX(0, cellIndex[1]);
                cellIndex[2] = MAX(0, cellIndex[2]);

                cellIndex[0] = MIN(param->nXcells,cellIndex[0]);
                cellIndex[1] = MIN(param->nYcells,cellIndex[1]);
                cellIndex[2] = MIN(param->nZcells,cellIndex[2]);

                cellIndex[0]++;
                cellIndex[1]++;
                cellIndex[2]++;

                node2->cellIdx = EncodeCellIdx(home,
                                               cellIndex[0],
                                               cellIndex[1],
                                               cellIndex[2]);
	        break;

/*
 *			Mark a node so that new force and velocity
 *			values will be computed for the node.  Probably
 *			because topology changes in another domain have
 *			affected this node.
 */
			case(MARK_FORCES_OBSOLETE):
				

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				node1->flags |= NODE_RESET_FORCES;

				break;

			case(MARK_NO_MESH_COARSEN):
				

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				node1->flags |= NO_MESH_COARSEN;

				break;

			case(MARK_NO_COLLISIONS):
				

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				node1->flags |= NO_COLLISIONS;

				break;
/*
 *			Reset the forces on the specified segment of node1.
 *                      The op structure x,y,z values are loaded here with
 *			force values.
 */
			case(RESET_SEG_FORCES):
			

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				ResetSegForces(home, node1, &tag2, op->x,
					       op->y, op->z, 0);

				break;

			case(RESET_SEG_FORCES2):
			

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);

				if (node1 == (Node_t *)NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				ResetSegForces2(home, node1, &tag2,
                                               op->x, op->y, op->z,
                                               op->nx, op->ny, op->nz, 0);

				break;

/*
 *			Reset the coordinates of the specified node.
 */
			case(RESET_COORD):

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
#if CHECK_FIXREMESH				
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("RC (%d,%d), (%f,%f,%f) ",op->dom1,op->idx1,
					  op->x,op->y,op->z);

				PrintTag(node1);
				printf("\n");
#endif				
				if (node1 == (Node_t *)NULL) break;

				node1->x = op->x;
				node1->y = op->y;
				node1->z = op->z;

                if(fabs(op->nx)>1E-10 || fabs(op->ny)>1E-10 || fabs(op->nz)>1E-10){
                    node1->oldx = op->nx;
                    node1->oldy = op->ny;
                    node1->oldz = op->nz;
                }
				break;

/*
 * 			Update the grain of node. 
 */				
				case(UPDATE_GRAIN):

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
#if CHECK_FIXREMESH				
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("UG (%d,%d), %d %d ",op->dom1,op->idx1,
					  op->inGrain, op->onGB);

				PrintTag(node1);
				printf("\n");
#endif				
				if (node1 == (Node_t *)NULL) {
					break;
				}

				if (node1 == (Node_t *)NULL) break;

                if(node1->inGrain == op->inGrain &&
                   node1->onGB == op->onGB)break;
				

				if(op->dom1 == home->myDomain){
					PrintTag(node1);
					printf("Task %d, Warning , the native node change inGrain %d %d\n",
                           home->myDomain, op->inGrain, op->onGB);
				}
#if 0
                if(op->inGrain > -1){
                    for(i=0; i<node1->numNbrs; i++){
                        if((nbr = GetNeighborNode(home,node1,i)) == (Node_t *)NULL)continue;

                        if(!IsOnGB(nbr) && nbr->inGrain != op->inGrain){
                            printf("Warning: task %d, flag %s, the node (%d,%d) change grain from ", 
                                home->myDomain, flag, node1->myTag.domainID, node1->myTag.index);
                            printf("(%d,%d) to (%d,%d).\n", node1->inGrain, node1->onGB,
                                op->inGrain, op->onGB);
                            PrintTag(node1);
                            PrintTag(nbr);
                            printf("\n");
                        }
                    }

                }

#endif
				node1->inGrain = op->inGrain;
				node1->onGB = op->onGB;


				break;
				
			case(UPDATE_CONSTRAINT):
				node1 = GetNodeFromIndex(home, op->dom1,op->idx1);
#if CHECK_FIXREMESH				
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("UC (%d,%d), %d ",op->dom1,op->idx1,
					  op->inGrain);
				PrintTag(node1);
				printf("\n");
#endif				

				if (node1 == (Node_t *)NULL) {
					break;
				}
				if (node1 == (Node_t *)NULL) break;

				node1->constraint = op->inGrain;
				
				break;

/*
 *			Remove the specified node
 */
			case(REMOVE_NODE) :

				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
#if CHECK_FIXREMESH				
				printf("FR %s: Task %d, ", flag, home->myDomain);
				printf("RN (%d,%d) ",op->dom1,op->idx1);

				PrintTag(node1);
				printf("\n");
#endif
				if (node1 == NULL) break;

				RemoveNode(home, node1, 0);
				break;

				
/*
 *			Change the burgers vector of a  node's arm
 *			to the given burgers vector.  If the new
 *			burgers vector is zero, the arm will be
 *			removed.
 */
			case(CHANGE_ARM_BURG) :
				
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
				
#if CHECK_FIXREMESH	
				node2 = GetNodeFromIndex(home, op->dom2,
										 op->idx2);
				
				printf("FR %s: Task %d, ", flag, home->myDomain);
				if(op->bx == 0.0 && op->by == 0.0 && op->by == 0.0){
					printf("DC (%d,%d), (%d,%d).",op->dom1,op->idx1,
					       op->dom2,op->idx2);				
				}else{
					printf("CAB (%d,%d), (%d,%d)--> ",op->dom1,op->idx1,
					       op->dom2,op->idx2);				
					printf("[%f,%f,%f](%f,%f,%f). ",op->bx,op->by,op->bz,
	                                       op->nx, op->ny, op->nz);
				}

				PrintTag(node1);
				PrintTag(node2);
				printf("\n");
#endif
				if (node1 == (Node_t *)NULL) break;

				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				ChangeArmBurg(home, node1, &tag2,
					      op->bx, op->by, op->bz,
					      op->nx, op->ny, op->nz,
					      0, DEL_SEG_NONE);
				break;

/*
 *			Reset the glide plane for a specific segment.
 */
			case RESET_GLIDE_PLANE:

				tag1.domainID = op->dom1;
				tag1.index    = op->idx1;
				tag2.domainID = op->dom2;
				tag2.index    = op->idx2;

				newPlane[0] = op->nx;
				newPlane[1] = op->ny;
				newPlane[2] = op->nz;

				ResetGlidePlane(home, newPlane,
						&tag1, &tag2, 0);
				break;
            
            case CHANGE_ARM_SIGNSF:
				node1 = GetNodeFromIndex(home, op->dom1,
							 op->idx1);
                
                if(node1 == (Node_t *)NULL)break;

                for(i=0; i<node1->numNbrs; i++){
                    if(op->dom2 == node1->nbrTag[i].domainID &&
                       op->idx2 == node1->nbrTag[i].index){
                        node1->signSF[i] = op->bx;
                    }
                }                

                break;
			default:
				Fatal("FixRemesh : invalid Remesh operation");
			}
		}

/*
 *		We're done with this domain's opList; release it
 */
		free(remDom->inBuf);
	}
	return;
}
