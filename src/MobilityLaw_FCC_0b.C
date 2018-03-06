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
	int	    *numGlideConstraints, numPointer;
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
	real8   (*invDragMatrix)[3], (*glideConstraints)[3];
	real8   invDragpoArray[3][3], glideConsArray[3][3];
	real8   invBtotal[3][3];
	real8   normCrystal[3][3];
	
	real8   velocity,directionV[3];

	real8	forceGB[3];
	Node_t  *nbrNode;
	int	nodeType;
	Param_t *param;

	param  = home->param;

	nodeType = NodeOnGB(home, node);
	if((nodeType & TB_NODE) && !(nodeType & INTERSECT_NODE)){
	    Bscrew = 1.0 / param->MobScrew / 10.0;
	    Bedge  = 1.0 / param->MobEdge / 10.0;
	}else{
	    Bscrew = 1.0 / param->MobScrew;
	    Bedge  = 1.0 / param->MobEdge;
	}
	
	invDragMatrix       = invDragpoArray;
	glideConstraints    = glideConsArray;
	numGlideConstraints = &numPointer;

	Matrix33Init(invDragMatrix);
	*numGlideConstraints= 0;
	
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
	grainNode = (node->inGrain -floor(node->inGrain) == 0 ? 1 : 0);

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
	    numNorms = 0; 
		
	    if ((nbrNode = GetNeighborNode(home, node, 0)) == (Node_t *)NULL) {
		Fatal("WARNING: Neighbor not found at %s line %d\n",
		      __FILE__, __LINE__);
	    }
	    if(grainNode){	
		glideConstraints[0][XC] = node->nx[0];
		glideConstraints[0][YC] = node->ny[0];
		glideConstraints[0][ZC] = node->nz[0]; 
	    }else{
		glideConstraints[0][XC] = param->TBNormal[XC];
		glideConstraints[0][YC] = param->TBNormal[YC];
		glideConstraints[0][ZC] = param->TBNormal[ZC];
	    }

	    NormalizeVec(glideConstraints[0]);

	    if (param->useLabFrame) {
		if(grainNode){
		    ChangeOrientation(home, node, nbrNode, 3);
		}else{
		    ChangeOrientation(home, node, node, 3);
		}
		Matrix33Vector3Multiply(home->rotMatrixInverse, glideConstraints[0],
					normCrystal[0]);
	    } else {
		VECTOR_COPY(normCrystal[0], glideConstraints[0]);
	    }

	/*
	 *  if this test passes then the dislocaiton segment is a junction
	 *  segment which will be constrained to grow and shrink but not glide.
	 *  Here glide planes are assumed to be of 111 type and any zero in the
	 *  plane data is assumed to be associated with a non-glide plane
	 */
	    temp = normCrystal[0][XC] *
		   normCrystal[0][YC] *
		   normCrystal[0][ZC];

	    if(fabs(temp) < eps && !grainNode){
		FormatVector(normCrystal[0],"TB direction in Crystal");
		Fatal("the TB plane is a junction plane");
	    }

	    if (fabs(temp) < eps) {			// Junction plane

		numNorms = 1;

		glideConstraints[2][XC] = node->x - nbrNode->x;
		glideConstraints[2][YC] = node->y - nbrNode->y;
		glideConstraints[2][ZC] = node->z - nbrNode->z;

		ZImage(param, &glideConstraints[2][XC], &glideConstraints[2][YC], 
				   &glideConstraints[2][ZC]);

		NormalizeVec(glideConstraints[2]);
		cross(glideConstraints[2], glideConstraints[0], glideConstraints[1]);
			NormalizeVec(glideConstraints[1]);
	    }
	    

	    if (numNbrs > 1) { 

		i = grainNode;

		while ((numNorms < 2) && (i < numNbrs)){

		    nx = node->nx[i];
		    ny = node->ny[i];
		    nz = node->nz[i];

		    Normalize(&nx, &ny, &nz);
	    
		    if (numNorms == 0) {

			temp = fabs(glideConstraints[0][XC] * nx +
				    glideConstraints[0][YC] * ny +
				    glideConstraints[0][ZC] * nz);

			if (fabs(1.0e0-temp) > tor) {		// second effective glide plane !

			    numNorms = 1;

			    nx -= temp * glideConstraints[0][XC];
			    ny -= temp * glideConstraints[0][YC];
			    nz -= temp * glideConstraints[0][ZC];

			    Normalize(&nx, &ny, &nz);

			    glideConstraints[1][XC] = nx;
			    glideConstraints[1][YC] = ny; 
			    glideConstraints[1][ZC] = nz;

			    cross(glideConstraints[0], 
							  glideConstraints[1], 
							  glideConstraints[2]);
			}

		    } else {/* numNorms==1*/

			temp = glideConstraints[2][XC] * nx +
			       glideConstraints[2][YC] * ny +
			       glideConstraints[2][ZC] * nz;

			if (fabs(temp) > tor ) {		// third effective glide plane
			    numNorms = 2;
			}
		    }

	/*
	 *          check to see if the normal is a non-glide plane and then add to
	 *          the checks so that it constrains the junction dislocation to
	 *          only move along its line
	 */
		    temp = normCrystal[0][XC] *
			   normCrystal[0][YC] *
			   normCrystal[0][ZC];

		    if (fabs(temp) < eps) {

			nbrNode = GetNeighborNode(home, node, i);

			dx = node->x - nbrNode->x;
			dy = node->y - nbrNode->y;
			dx = node->z - nbrNode->z;

			ZImage(param, &dx, &dy, &dz);
			Normalize(&dx, &dy, &dz);
			xvector(nx, ny, nz, dx, dy, dz, &nx2, &ny2, &nz2);

			if (numNorms == 0) {

			    temp = fabs(glideConstraints[0][XC] * nx2 +
					glideConstraints[0][YC] * ny2 +
					glideConstraints[0][ZC] * nz2);

			    if (fabs(1.0e0 - temp) > tor ) {

				numNorms = 1;

				nx2 -= temp * glideConstraints[0][XC];
				ny2 -= temp * glideConstraints[0][YC];
				nz2 -= temp * glideConstraints[0][ZC];

				Normalize(&nx2, &ny2, &nz2);

				glideConstraints[1][XC] = nx2;
				glideConstraints[1][YC] = ny2; 
				glideConstraints[1][ZC] = nz2;

				cross(glideConstraints[0], glideConstraints[1], glideConstraints[2]);
			    }

			} else {/* numNorms==1*/

			    temp = glideConstraints[2][XC] * nx2 +
				   glideConstraints[2][YC] * ny2 +
				   glideConstraints[2][ZC] * nz2;

			    if (fabs(temp) > tor) {
				numNorms = 2;
			    }
			}
		    }

		    i++;

		} /* end while */

	    } /* end if (numNbrs > 1) */   

	/*
	 *  Glide constraints have been set, now initialize the velocity
	 */
		node->vX = 0.0e0;
		node->vY = 0.0e0;
		node->vZ = 0.0e0;

	/*
	 *      Begin construction of the node drag matrix
	 *
	 *      Loop over all arms of the node, adding each arm's contribution
	 *      to the drag matrix.
	 */
		if (numNorms < 2) {
		    for (i = 0; i < numNbrs; i++) {
			real8 b[3], d[3], m[3], n[3], nCryst[3];

			b[XC] = node->burgX[i];
			b[YC] = node->burgY[i];
			b[ZC] = node->burgZ[i];

			bMag2 = (b[XC]*b[XC] + b[YC]*b[YC] + b[ZC]*b[ZC]);
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

			d[XC] = nbrNode->x - node->x;
			d[YC] = nbrNode->y - node->y;
			d[ZC] = nbrNode->z - node->z;

			ZImage(param, &d[XC], &d[YC], &d[ZC]);

			mag2 = d[XC]*d[XC] + d[YC]*d[YC] + d[ZC]*d[ZC];

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

			d[XC] *= invMag;
			d[YC] *= invMag;
			d[ZC] *= invMag;

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
                costheta = (d[XC]*b[XC] + d[YC]*b[YC] + d[ZC]*b[ZC]);
                costheta2 = (costheta*costheta) * invbMag2;
    
/*
 *              arms not on [1 1 1] planes don't move as readily as
 *              other arms, so must be handled specially.
 *
 *              If needed, rotate a copy of the glide plane vector from the
 *              laboratory frame to the crystal frame.
 */
                n[XC] = node->nx[i];
                n[YC] = node->ny[i];
                n[ZC] = node->nz[i];

                if (param->useLabFrame) {
	            ChangeOrientation(home, node, nbrNode, 3);
                    Matrix33Vector3Multiply(home->rotMatrixInverse, n, nCryst);
                } else {
                    VECTOR_COPY(nCryst, n);
                }

                if (fabs(nCryst[XC] * nCryst[YC] * nCryst[ZC]) < eps) {
                    if (numNbrs == 2) {
                        Btotal[0][0] += halfMag * Beclimb;
                        Btotal[1][1] += halfMag * Beclimb;
                        Btotal[2][2] += halfMag * Beclimb;
                    } else {
                        Btotal[0][0] += halfMag * (d[XC]*d[XC]*BlmBecl + Beclimb);
                        Btotal[0][1] += halfMag * (d[XC]*d[YC]*BlmBecl);
                        Btotal[0][2] += halfMag * (d[XC]*d[ZC]*BlmBecl);
                        Btotal[1][1] += halfMag * (d[YC]*d[YC]*BlmBecl + Beclimb);
                        Btotal[1][2] += halfMag * (d[YC]*d[ZC]*BlmBecl);
                        Btotal[2][2] += halfMag * (d[ZC]*d[ZC]*BlmBecl + Beclimb);
                    }
                } else  {
/*
 *                  Arm is a regular glide arm, so build the drag matrix
 *                  assuming the dislocation is screw type
 */
                    Btotal[0][0] += halfMag * (d[XC]*d[XC] * BlmBsc + Bscrew);
                    Btotal[0][1] += halfMag * (d[XC]*d[YC] * BlmBsc);
                    Btotal[0][2] += halfMag * (d[XC]*d[ZC] * BlmBsc);
                    Btotal[1][1] += halfMag * (d[YC]*d[YC] * BlmBsc + Bscrew);
                    Btotal[1][2] += halfMag * (d[YC]*d[ZC] * BlmBsc);
                    Btotal[2][2] += halfMag * (d[ZC]*d[ZC] * BlmBsc + Bscrew);

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


                        Btotal[0][0] += halfMag * (n[XC]*n[XC] * BclmBsc +
                                        m[XC]*m[XC] * BglmBsc);
                        Btotal[0][1] += halfMag * (n[XC]*n[YC] * BclmBsc +
                                        m[XC]*m[YC] * BglmBsc);
                        Btotal[0][2] += halfMag * (n[XC]*n[ZC] * BclmBsc +
                                        m[XC]*m[ZC] * BglmBsc);
                        Btotal[1][1] += halfMag * (n[YC]*n[YC] * BclmBsc +
                                        m[YC]*m[YC] * BglmBsc);
                        Btotal[1][2] += halfMag * (n[YC]*n[ZC] * BclmBsc +
                                        m[YC]*m[ZC] * BglmBsc);
                        Btotal[2][2] += halfMag * (n[ZC]*n[ZC] * BclmBsc +
                                        m[ZC]*m[ZC] * BglmBsc);
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

	    DoResistance(home, node, nForce);

            if (Matrix33Invert(Btotal, invBtotal) == 0) {
                Fatal("Mobility_FCC_0b: Cannot invert 3X3 matrix!");
            }

            Matrix33Vector3Multiply(invBtotal, nForce, nVel);

/*
 *          Now that the nVel has been set, apply glide corrections
 */
            if (numNorms == 0){

                temp = nVel[0] * glideConstraints[0][XC] +
                       nVel[1] * glideConstraints[0][YC] +
                       nVel[2] * glideConstraints[0][ZC];

                nVel[0] -= temp * glideConstraints[0][XC];
                nVel[1] -= temp * glideConstraints[0][YC];
                nVel[2] -= temp * glideConstraints[0][ZC];

            } else { /*numNorms==1*/

                temp = nVel[0] * glideConstraints[2][XC] +
                       nVel[1] * glideConstraints[2][YC] +
                       nVel[2] * glideConstraints[2][ZC];

                nVel[0] = temp * glideConstraints[2][XC];
                nVel[1] = temp * glideConstraints[2][YC];
                nVel[2] = temp * glideConstraints[2][ZC];
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
		

        AdjustSurfaceNodeVel(home->param, node);

#ifdef _FEM
/*
 *      The velocity of any surface node along the negative surface
 *      normal direction should be zero to prevent the node moving into
 *      the box.  Make a call to adjsut the velocity accordingly.
 *
 *      Note: If the node moves outside the box, it's allowed here, then
 *      position is adjusted in AdjustNodePosition().
 */
        AdjustSurfaceNodeVel(home, node);
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
