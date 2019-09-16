/**************************************************************************
 *
 *  Function    : Mobility_FCC_0b
 *  Author      : Tom Arsenlis
 *  Description : This is alternate version of a generic mobility function
 *                for FCC materials.  It is very similar in function and
 *                structure to the BCC_glide mobility although glide planes
 *                rather than burgers vectors are used to identify junctions.
 *                Additionally, this function attempts to identify nodes
 *                that both have short segments and have reversed the 
 *                their direction of velocity.  When found, these nodes
 *                are artifically slowed for a single timestep in an effort
 *                to dampen motion of flickering/oscillating nodes.
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
 *
 ***************************************************************************/

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif


int Mobility_FCC_0b(Home_t *home, Node_t *node)
{
	int     i, numNbrs, numNorms,grainNode; 
	int     numNonZeroLenSegs = 0, hasShortSeg = 0;
	int	    numGlideConstraints;
	real8   dx, dy, dz;
	real8   nx, ny, nz, nx2, ny2, nz2; 
	real8   temp;
	real8   mag, mag2, halfMag, invMag;
	real8   invbMag2, bMag2, costheta, costheta2;
	real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
	real8   Bscrew2, Beclimb2;
	real8   invBscrew2, invBedge2;
	real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
	real8   shortSegCutoff;
	real8   eps = 1.0e-12;
	real8   tor=1.e-5;
	real8   shortSegRatio = 1.0;
	real8   nForce[3], nVel[3];
	real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
	                        {0.0, 0.0, 0.0},
	                        {0.0, 0.0, 0.0}};
	real8   (*invDragMatrix)[3], glideConstraints[3][3];
	real8   invDragpoArray[3][3];
	real8   invBtotal[3][3];
	real8   normCrystal[3][3];
	
	real8   velocity,directionV[3];
	real8	forceGB[3];
	Node_t  *nbrNode;
	int	nodeType;
	Param_t *param;
    bool    zeroLenJunc;

	param  = home->param;

	Bscrew = 1.0 / param->MobScrew;
	Bedge  = 1.0 / param->MobEdge;
	
	invDragMatrix       = invDragpoArray;

	Matrix33Init(invDragMatrix);
	
/*
 *  Climb is set this way to help with convergence and the fact that this
 *  is a glide restricted mobility
 */
	Beclimb    = 100.0 * Bedge;
	 
	Bscrew2    = Bscrew * Bscrew;
	Beclimb2   = Beclimb * Beclimb;
	
	Bline      = 1.0 * MIN(Bscrew, Bedge);
	BlmBsc     = Bline - Bscrew;
	BlmBecl    = Bline - Beclimb; 
	
	invBscrew2 = 1.0 / (Bscrew*Bscrew);
	invBedge2  = 1.0 / (Bedge*Bedge);
	
	numNbrs = node->numNbrs;
	shortSegCutoff = 0.5 * param->minSeg;

/*
 *  If node is 'pinned' in place by constraints, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
	if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
	    node->vX = 0.0;
	    node->vY = 0.0;
	    node->vZ = 0.0;
	    return(0);
	}

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
	for (i = 0; i < numNbrs; i++) {
	    if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
	        continue;
	    }

	    dx = node->x - nbrNode->x;
	    dy = node->y - nbrNode->y;
	    dz = node->z - nbrNode->z;
	    
	    ZImage(param, &dx, &dy, &dz);

/*
 *          This assumes that lengths are given in terms of burgers vector units
 */
	    if ((dx*dx + dy*dy + dz*dz) > eps) {
	        numNonZeroLenSegs++;
	    }
	}

	if (numNonZeroLenSegs == 0) {
	    node->vX = 0.0;
	    node->vY = 0.0;
	    node->vZ = 0.0;
	    return(0);
	}

/*
 *  Get the normal for the first segment in both the lab and crystal frames
 */
    GetGlideConstraintList(home, node, 
                           &numGlideConstraints, glideConstraints, zeroLenJunc);

/*
 *  Glide constraints have been set, now initialize the velocity
 */
	node->vX = 0.0e0;
	node->vY = 0.0e0;
	node->vZ = 0.0e0;

/*
 *  Begin construction of the node drag matrix
 *  Loop over all arms of the node, adding each arm's contribution
 *  to the drag matrix.
 */
		if (numNorms < 2) {
		    for (i = 0; i < numNbrs; i++) {
			real8 b[3], d[3], m[3], n[3], nCryst[3];

			b[0] = node->burgX[i];
			b[1] = node->burgY[i];
			b[2] = node->burgZ[i];

			bMag2 = (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
			invbMag2 = 1.0 / bMag2;

	/*
	 *              Calculate the length of the arm and its tangent line direction
	 */
			nbrNode = GetNeighborNode(home, node, i);

			if (nbrNode == (Node_t *)NULL) {
			    printf("WARNING: Neighbor not found at %s line %d\n",
				   __FILE__, __LINE__);
			    continue;
			}

			d[0] = nbrNode->x - node->x;
			d[1] = nbrNode->y - node->y;
			d[2] = nbrNode->z - node->z;

			ZImage(param, &d[0], &d[1], &d[2]);

			mag2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

	/*
	 *              If the segment is zero length (which can happen when
	 *              the mobility function is being called from SplitMultiNodes())
	 *              just skip the segment.
	 */
			if (mag2 < eps) {
			    continue;
			}

			mag     = sqrt(mag2);
			halfMag = mag/2.0;
			invMag  = 1.0 / mag;

			d[0] *= invMag;
			d[1] *= invMag;
			d[2] *= invMag;

	/*
	 *              If the node has a very short segment, we need to flag it
	 *              for later
	 */
	            if (mag < shortSegCutoff) {
	                hasShortSeg = 1;
                    shortSegRatio = MIN(shortSegRatio, mag/param->minSeg);
                }

/*
 *              Calculate how close to screw the arm is
 */
                costheta = (d[0]*b[0] + d[1]*b[1] + d[2]*b[2]);
                costheta2 = (costheta*costheta) * invbMag2;
    
/*
 *              arms not on [1 1 1] planes don't move as readily as
 *              other arms, so must be handled specially.
 *
 *              If needed, rotate a copy of the glide plane vector from the
 *              laboratory frame to the crystal frame.
 */
                n[0] = node->nx[i];
                n[1] = node->ny[i];
                n[2] = node->nz[i];

                if (param->useLabFrame) {
                   VectorLabToCrystal(home, node, nbrNode, n, nCryst);
                } else {
                    VECTOR_COPY(nCryst, n);
                }

                if (fabs(nCryst[0] * nCryst[1] * nCryst[2]) < eps) {
                    if (numNbrs == 2) {
                        Btotal[0][0] += halfMag * Beclimb;
                        Btotal[1][1] += halfMag * Beclimb;
                        Btotal[2][2] += halfMag * Beclimb;
                    } else {
                        Btotal[0][0] += halfMag * (d[0]*d[0]*BlmBecl + Beclimb);
                        Btotal[0][1] += halfMag * (d[0]*d[1]*BlmBecl);
                        Btotal[0][2] += halfMag * (d[0]*d[2]*BlmBecl);
                        Btotal[1][1] += halfMag * (d[1]*d[1]*BlmBecl + Beclimb);
                        Btotal[1][2] += halfMag * (d[1]*d[2]*BlmBecl);
                        Btotal[2][2] += halfMag * (d[2]*d[2]*BlmBecl + Beclimb);
                    }
                } else  {
/*
 *                  Arm is a regular glide arm, so build the drag matrix
 *                  assuming the dislocation is screw type
 */
                    Btotal[0][0] += halfMag * (d[0]*d[0] * BlmBsc + Bscrew);
                    Btotal[0][1] += halfMag * (d[0]*d[1] * BlmBsc);
                    Btotal[0][2] += halfMag * (d[0]*d[2] * BlmBsc);
                    Btotal[1][1] += halfMag * (d[1]*d[1] * BlmBsc + Bscrew);
                    Btotal[1][2] += halfMag * (d[1]*d[2] * BlmBsc);
                    Btotal[2][2] += halfMag * (d[2]*d[2] * BlmBsc + Bscrew);

/*
 *                  Now correct the drag matrix for dislocations that are
 *                  not screw
 */
                    if ((1.0 - costheta2) > eps) {

                        cross(n, d, m);

                        Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) *
                                      costheta2);
                        Bglide = 1.0 / Bglide;
                        Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) *
                                      costheta2);

#ifdef NAN_CHECK
                        if (isnan(Bglide) != 0) {
                            Fatal("Mobility_FCC_0b: Bglide = NaN\n"
                                "  Bglide = sqrt(invBedge2 + "
                                "(invBscrew2-invBedge2)*costheta2)\n"
                                "  where invBedge2 = %lf, invBscrew2 = %lf, "
                                "costheta2 = %lf", invBedge2, invBscrew2,
                                costheta2);
                        }

                        if (isnan(Bclimb) != 0) {
                            Fatal("Mobility_FCC_0b: Bclimb = NaN\n"
                                "  Bclimb = sqrt(Beclimb2 + "
                                "(Bscrew2-Beclimb2)*costheta2)\n"
                                "  where Beclimb2 = %lf, Bscrew2 = %lf, "
                                "costheta2 = %lf", Beclimb2, Bscrew2,
                                costheta2);
                        }
#endif
                        BclmBsc = Bclimb - Bscrew;
                        BglmBsc = Bglide - Bscrew;


                        Btotal[0][0] += halfMag * (n[0]*n[0] * BclmBsc +
                                        m[0]*m[0] * BglmBsc);
                        Btotal[0][1] += halfMag * (n[0]*n[1] * BclmBsc +
                                        m[0]*m[1] * BglmBsc);
                        Btotal[0][2] += halfMag * (n[0]*n[2] * BclmBsc +
                                        m[0]*m[2] * BglmBsc);
                        Btotal[1][1] += halfMag * (n[1]*n[1] * BclmBsc +
                                        m[1]*m[1] * BglmBsc);
                        Btotal[1][2] += halfMag * (n[1]*n[2] * BclmBsc +
                                        m[1]*m[2] * BglmBsc);
                        Btotal[2][2] += halfMag * (n[2]*n[2] * BclmBsc +
                                        m[2]*m[2] * BglmBsc);
                    }
                }  /* End non-[0 0 1] arm */
            }  /* End loop over arms */

            Btotal[1][0] = Btotal[0][1];
            Btotal[2][0] = Btotal[0][2];
            Btotal[2][1] = Btotal[1][2];

/*
 *          At this point we should check if the matrix is invertable and
 *          if it isn't, find the eigen values and eigen vectors of the drag
 *          matrix, then invert the drag matrix keeping zero eigen values
 *          as zero.
 *
 *          FIX ME!  For now, we're assuming the matrix is invertable.
 */

#if 0
	DoResistance(home, node, nForce);
#else
    nForce[0] = node->fX;
    nForce[1] = node->fY;
    nForce[2] = node->fZ;
#endif

            if (Matrix33Invert(Btotal, invBtotal) == 0) {
                Fatal("Mobility_FCC_0b: Cannot invert 3X3 matrix!");
            }

            Matrix33Vector3Multiply(invBtotal, nForce, nVel);

/*
 *          Now that the nVel has been set, apply glide corrections
 */
            if (numNorms == 0){

                temp = nVel[0] * glideConstraints[0][0] +
                       nVel[1] * glideConstraints[0][1] +
                       nVel[2] * glideConstraints[0][2];

                nVel[0] -= temp * glideConstraints[0][0];
                nVel[1] -= temp * glideConstraints[0][1];
                nVel[2] -= temp * glideConstraints[0][2];

            } else { /*numNorms==1*/

                temp = nVel[0] * glideConstraints[2][0] +
                       nVel[1] * glideConstraints[2][1] +
                       nVel[2] * glideConstraints[2][2];

                nVel[0] = temp * glideConstraints[2][0];
                nVel[1] = temp * glideConstraints[2][1];
                nVel[2] = temp * glideConstraints[2][2];
            }

            node->vX = nVel[0];
            node->vY = nVel[1];
            node->vZ = nVel[2];
        }

/*
 *      Sometimes nodes oscillate at high frequency causing timestep to drop.
 *      Usually, these are 2-nodes that have short segments attached.  So, if
 *      we find a 2-node that has a least 1 short segment, and the velocity
 *      of the node has changed direction since the previous cycle, slow
 *      the node down for this cycle.  Hopefully this will help dampen the
 *      oscillation and improve the timestep.
 */
        if ((numNbrs == 2) && (hasShortSeg)) {
            if (((node->vX * node->oldvX) < 0.0) ||
                ((node->vY * node->oldvY) < 0.0) ||
                ((node->vZ * node->oldvZ) < 0.0)) {
                node->vX *= shortSegRatio;
                node->vY *= shortSegRatio;
                node->vZ *= shortSegRatio;
            }
        }
		

#ifdef _FEM
/*
 *      The velocity of any surface node along the negative surface
 *      normal direction should be zero to prevent the node moving into
 *      the box.  Make a call to adjsut the velocity accordingly.
 *
 *      Note: If the node moves outside the box, it's allowed here, then
 *      position is adjusted in AdjustNodePosition().
 */
//        AdjustSurfaceNodeVel(home, node);
#endif

	velocity = sqrt(node->vX * node->vX+
			        node->vY * node->vY+
			        node->vZ * node->vZ);
	if(velocity > param->vshr && param->vshr > 0){
	    node->vX = node->vX * param->vshr /velocity ;
	    node->vY = node->vY * param->vshr /velocity ;
	    node->vZ = node->vZ * param->vshr /velocity ;

	}	

        return(0);
}
