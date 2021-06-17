/****************************************************************************
 *
 *      Module:         CrossSlipFCC.c
 *
 *      Author:         Converted/adapted from CrossSlipBCC.c
 *
 *      Description:    This module contains functions for allowing
 *                      dislocations in FCC materials to cross slip to
 *                      a glide plane other than its original plane
 *
 *                      NOTE:  See comments at the top of the CrossSlip.c
 *                             module for a general description of the
 *                             cross-slip mechanism.
 *
 *      Includes public functions:
 *          CrossSlipFCC()
 *
 ***************************************************************************/
#include <math.h>
#include "Home.h"
#include "Util.h"

static int dbgDom;


/*---------------------------------------------------------------------------
 *
 *      Function:       CrossSlipFCC
 *
 *      Description:    Examines all nodes local to the domain, determines
 *                      if the node should cross slip, if the node is
 *                      permitted to cross slip, and if so, adjusts node
 *                      positions and segment glide planes as necessary.
 *
 *-------------------------------------------------------------------------*/
void CrossSlipFCC(Home_t *home)
{
        int    i, j, n;
        int    plane1, plane2, fplane;
        int    pinned1, pinned2;
        int    opClass, thisDom;
        int    resetNodePos, resetNbr1Pos, resetNbr2Pos;
        int    nbr1ArmID, nbr2ArmID;
        int    seg1_is_screw, seg2_is_screw, bothseg_are_screw;
        int     burgID;
        real8  tmp, eps, thetacrit, sthetacrit, s2thetacrit, areamin;
        real8  test1, test2, test3, testmax1, testmax2, testmax3;
        real8  vec1dotb, vec2dotb, vec3dotb, fdotglide;
        real8  f1dotplane1, f1dotplane2, f1dotplanef, newfdotglide;
        real8  nodep[3], nbr1p[3], nbr2p[3];
        real8  vec1[3], vec2[3], vec3[3];
        real8  segplane1[3], segplane2[3], newplane[3];
        real8  tmp3[3], tmp3B[3], tmp3C[3];
        real8  glideDirLab[3][3], glideDirCrystal[3][3];
        real8  fLab[3], fCrystal[3];
        real8  burgLab[3], burgCrystal[3];
        real8  nodePosOrig[3], nbr1PosOrig[3], nbr2PosOrig[3];
        real8  newForce[3], newSegForce[3];
        real8  segForceOrig[4][3];
        real8  noiseFactor, weightFactor, fnodeThreshold;
        real8  L1, L2, shearModulus, burgSize;
        real8  zipperThreshold;
        Node_t *node, *nbr1, *nbr2;
        Param_t *param;
        real8   rotMatrix[3][3], rotMatrixInverse[3][3];
        BurgInfo_t  *burgData;
        GB_t *gb;
        real8  lamp = 1.0, temVec[3], cof;
        int armInGrain1, armInGrain2, armOnGB1, armOnGB2, curInGrain, curOnGB;
        int iter, iterMax=8;

#ifdef DEBUG_CROSSSLIP_DOMAIN
        dbgDom = DEBUG_CROSSSLIP_DOMAIN;
#else
        dbgDom = -1;
#endif

        param = home->param;
        thisDom = home->myDomain;
        shearModulus = param->shearModulus;

        eps = 1.0e-06;
        thetacrit = 2.0 / 180.0 * M_PI;
        sthetacrit = sin(thetacrit);
        s2thetacrit = sthetacrit * sthetacrit;
        areamin = param->remeshAreaMin;

/*
 *      In order to cross-slip a node, the force on the cross-slip
 *      plane must exceed a value based on the force on the current
 *      plane plus an amount that should guarantee the force on
 *      the cross-slip plane is not just in the noise.  The next
 *      two values are used in calculating that needed force.
 */
        noiseFactor=1e-05;
        weightFactor=1.0;

/*
 *      For purposes of determining segment ownership, we need to treat
 *      cross slip operations the same as we do during collision handling.
 */
        opClass = OPCLASS_COLLISION;

/*
 *      Loop through all 2-nodes native to this domain
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if (node->numNbrs != 2) {
                continue;
            }

            if (node->constraint != UNCONSTRAINED) {
                continue;
            }

            nbr1 = GetNeighborNode(home, node, 0);
            nbr2 = GetNeighborNode(home, node, 1);

            if ((nbr1 == (Node_t *)NULL) ||
                (nbr2 == (Node_t *)NULL)) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }


            burgLab[0] = node->burgX[0];
            burgLab[1] = node->burgY[0];
            burgLab[2] = node->burgZ[0];
            VECTOR_COPY(burgCrystal, burgLab);

/*
 *          If needed, rotate a copy of the burgers vector from the lab
 *          frame to the crystal frame.  Otherwise, lab and crystal frames
 *          are identical.
 */
            if (param->useLabFrame) {
                VectorLabToCrystal(home, node, nbr1, burgLab, burgCrystal);
            }

            burgSize = sqrt(DotProduct(burgLab, burgLab));

            NormalizeVec(burgLab);
            NormalizeVec(burgCrystal);

/*
 *          We must not change the glide plane for any segment
 *          that is not 'owned' by the node in the current domain.
 *          Otherwise, neighboring domains *could* change the
 *          glide plane for the same segment... this would be bad.
 */
            if ((!DomainOwnsSeg(home, opClass, thisDom, &node->nbrTag[0])) ||
                (!DomainOwnsSeg(home, opClass, thisDom, &node->nbrTag[1]))) {
                continue;
            }

/*
 *          Only consider glide dislocations.  If burgers vector is not
 *          a [1 1 0] type, ignore it.
 */
            burgData = GetBurgData(home, node, nbr1);
            burgID = BurgID_FCC(burgData, burgLab);
            if(burgID >= burgData->numBurgVectors)continue;
            if((burgData->burgType[burgID] & EDGE_BURG) == 0)continue;
/*
 *          Only consider non-GB dislocations. If any node of the two
 *          segements is on GB, ignore it.
 *          FIXME: perhaps the neighors can be on GB, we need consider
 *          it with dislocation emission on GB.
 *                                                              D. Wei
 */
            if( param->polycrystal >0 ){
                if(IsOnGB(nbr1) || IsOnGB(nbr2) || IsOnGB(node))continue;
#if 1
                if(fabs(ClosedGB(home, node->x, node->y, node->z, &gb))<param->rc)continue;
                if(fabs(ClosedGB(home, nbr1->x, nbr1->y, nbr1->z, &gb))<param->rc)continue;
                if(fabs(ClosedGB(home, nbr2->x, nbr2->y, nbr2->z, &gb))<param->rc)continue;
#endif
            }
/*
 *          Before we commit to doing a cross-slip event, we may shift
 *          nodes and recalculate some forces to determine if it is
 *          better to cross-slip the node or leave it where it is.
 *          in order to do that, we need to save some of the nodal info
 *          so we can restore the original state if it is determined
 *          that the cross-slip is not wise.
 */
            nbr1ArmID = GetArmID(home, nbr1, node);
            nbr2ArmID = GetArmID(home, nbr2, node);

            SaveCrossSlipInfo(node, nbr1, nbr2, nbr1ArmID, nbr2ArmID,
                              segForceOrig, nodePosOrig, nbr1PosOrig,
                              nbr2PosOrig);

            GetArmGrain(home, nbr1, node, armInGrain1, armOnGB1);
            GetArmGrain(home, nbr2, node, armInGrain2, armOnGB2);

            lamp = MIN(LengthAmplification(home, armInGrain1, armOnGB1), LengthAmplification(home, armInGrain2, armOnGB2));
/*
 *          Get the positions of all three nodes and convert the neighbor
 *          node positions to the PBC coordinates relative to the main
 *          node.  These adjusted positions will be used during the
 *          cross slip algorithm, and if updated during the process, will
 *          be copied back into the respective nodal data structures and
 *          shifted back into the primary image space.
 */
            nodep[0] = node->x; nodep[1] = node->y; nodep[2] = node->z;
            nbr1p[0] = nbr1->x; nbr1p[1] = nbr1->y; nbr1p[2] = nbr1->z;
            nbr2p[0] = nbr2->x; nbr2p[1] = nbr2->y; nbr2p[2] = nbr2->z;

            PBCPOSITION(param, nodep[0], nodep[1], nodep[2],
                        &nbr1p[0], &nbr1p[1], &nbr1p[2]);
            PBCPOSITION(param, nodep[0], nodep[1], nodep[2],
                        &nbr2p[0], &nbr2p[1], &nbr2p[2]);

            resetNodePos = 0;
            resetNbr1Pos = 0;
            resetNbr2Pos = 0;

/*
 *          If the node is a point on a long screw then we can consider
 *          it for possible cross slip.
 */
            vec1[0] = nbr1p[0] - nbr2p[0];
            vec1[1] = nbr1p[1] - nbr2p[1];
            vec1[2] = nbr1p[2] - nbr2p[2];

            vec2[0] = nodep[0] - nbr1p[0];
            vec2[1] = nodep[1] - nbr1p[1];
            vec2[2] = nodep[2] - nbr1p[2];

            vec3[0] = nodep[0] - nbr2p[0];
            vec3[1] = nodep[1] - nbr2p[1];
            vec3[2] = nodep[2] - nbr2p[2];

/*
 *          Force vector (initially in laboratory frame)
 */
            fLab[0] = node->fX;
            fLab[1] = node->fY;
            fLab[2] = node->fZ;

/*
 *          Calculate some test conditions
 */
            test1 = DotProduct(vec1, burgLab);
            test2 = DotProduct(vec2, burgLab);
            test3 = DotProduct(vec3, burgLab);

            test1 = test1 * test1;
            test2 = test2 * test2;
            test3 = test3 * test3;

            testmax1 = DotProduct(vec1, vec1);
            testmax2 = DotProduct(vec2, vec2);
            testmax3 = DotProduct(vec3, vec3);

/*
 *          Set up the tests to see if this dislocation is close enough to
 *          screw to be considered for cross slip.  For a segment to be close
 *          to screw it must be within 2*thetacrit defined above
 */
            seg1_is_screw = ((testmax2 - test2) < (testmax2 * s2thetacrit));
            seg2_is_screw = ((testmax3 - test3) < (testmax3 * s2thetacrit));
            bothseg_are_screw =
                    (((testmax2 - test2) < (4.0 * testmax2 * s2thetacrit)) &&
                     ((testmax3 - test3) < (4.0 * testmax3 * s2thetacrit)) &&
                     ((testmax1 - test1) < (testmax1 * s2thetacrit)));

            if (seg1_is_screw || seg2_is_screw || bothseg_are_screw) {
/*
 *              Set the force threshold for noise level within the code
 */
                L1 = sqrt(testmax2);
                L2 = sqrt(testmax3);
                fnodeThreshold = noiseFactor * shearModulus * burgSize *
                                 0.5 * (L1 + L2);

/*
 *              Find which glide planes the segments are on
 *              e.g. for burg = [ 1  1  0 ], the two glide directions are
 *                              [ 1 -1  2 ] and
 *                              [ 1 -1 -2 ] 
 *
 *              Use burgers vectors in crystal frame to generate initial glide
 *              planes in crystal frame.
 */
                tmp = 1.0; 

                for (n = 0; n < 3; n++) {
                   if ( fabs(burgCrystal[n]) > eps ) {
                      glideDirCrystal[0][n] =
                              (burgCrystal[n]*tmp > 0) ? 1.0 : -1.0;
                      glideDirCrystal[1][n] =
                              (burgCrystal[n]*tmp > 0) ? 1.0 : -1.0;
                      tmp = -1.0;
                   } 
                   else {
                      glideDirCrystal[0][n] =  2.0;
                      glideDirCrystal[1][n] = -2.0;
                   }
                }

/*
 *              For FCC there are only two slip planes for screw dislocation
 */
                VECTOR_ZERO(glideDirCrystal[2]);

/*
 *              Normalization
 */
                for (n = 0; n < 3; n++) {
                    glideDirCrystal[0][n] *= sqrt(1.0/6.0);
                    glideDirCrystal[1][n] *= sqrt(1.0/6.0);
                }

                segplane1[0] = node->nx[0];
                segplane1[1] = node->ny[0];
                segplane1[2] = node->nz[0];

                segplane2[0] = node->nx[1];
                segplane2[1] = node->ny[1];
                segplane2[2] = node->nz[1];

/*
 *              Need copies of the glideDir matrix and the force vector
 *              in both the lab frame and the crystal frame for later use.
 *              Also need to rotate the segment planes into the crystal frame.
 */
                if (param->useLabFrame) {
                    int     row;
                    real8   tmpPlane[3];

                    if(param->polycrystal){
                        Grain_t *grain;
                        grain = home->grainKeys[node->inGrain];
                        memcpy(rotMatrix, grain->orientation, sizeof(real8 [3][3]));
                        memcpy(rotMatrixInverse, grain->orientationInv, sizeof(real8 [3][3]));
                    }else{
                        memcpy(rotMatrix, home->rotMatrix, sizeof(real8 [3][3]));
                        memcpy(rotMatrixInverse, home->rotMatrixInverse, sizeof(real8 [3][3]));
                    }   
/*
 *                  Rotate initial glide dir from crystal to lab frame
 */
                    for (row = 0; row < 3; row++) {
                        Matrix33Vector3Multiply(rotMatrix,
                                                glideDirCrystal[row],
                                                glideDirLab[row]);
                    }
/*
 *                  Rotate segment planes to crystal frame
 */
                    Matrix33Vector3Multiply(rotMatrixInverse, segplane1,
                                            tmpPlane);
                    VECTOR_COPY(segplane1, tmpPlane);

                    Matrix33Vector3Multiply(rotMatrixInverse, segplane2,
                                            tmpPlane);
                    VECTOR_COPY(segplane2, tmpPlane);

/*
 *                  Rotate force vector to crystal frame
 */
                    Matrix33Vector3Multiply(rotMatrixInverse, fLab,
                                            fCrystal);
                } else {
/*
 *                  Lab and crystal frames are identical, so just copy
 *                  the vectors
 */
                    VECTOR_COPY(glideDirLab[0], glideDirCrystal[0]);
                    VECTOR_COPY(glideDirLab[1], glideDirCrystal[1]);
                    VECTOR_COPY(glideDirLab[2], glideDirCrystal[2]);
                    VECTOR_COPY(fCrystal, fLab);
                }

                Matrix33Vector3Multiply(glideDirCrystal, segplane1, tmp3);
                Matrix33Vector3Multiply(glideDirCrystal, segplane2, tmp3B);
                Matrix33Vector3Multiply(glideDirCrystal, fCrystal, tmp3C);

                plane1 = 0;
                plane2 = 0;
                fplane = 0;

/*
 *              For FCC there are only two slip planes for screw dislocation
 */
                for (j = 1; j < 2; j++) {
                    plane1 = (fabs(tmp3[j])  < fabs(tmp3[plane1]) ) ? j:plane1;
                    plane2 = (fabs(tmp3B[j]) < fabs(tmp3B[plane2])) ? j:plane2;
                    fplane = (fabs(tmp3C[j]) > fabs(tmp3C[fplane])) ? j:fplane;
                }

/*
 *              Calculate the new plane in the lab frame since it will only be
 *              used from here on out to update the existing nodal data which
 *              is in the lab frame.
 */
                NormalizedCrossVector(burgLab, glideDirLab[fplane], newplane);

                if ((bothseg_are_screw) && (plane1 == plane2) &&
                    (plane1 != fplane) &&
                    (fabs(tmp3C[fplane]) >
                     (weightFactor*fabs(tmp3C[plane1])+fnodeThreshold))) {
/*
 *                  Both segments are close to screw and the average direction
 *                  is close to screw.
 *
 *                  Determine if the neighbor nodes should be considered
 *                  immobile.
 */
                    pinned1 = NodePinned(home, nbr1, plane1, glideDirLab);
                    pinned2 = NodePinned(home, nbr2, plane2, glideDirLab);

                if (pinned1) {
                    if ((!pinned2) ||
                        ((testmax1-test1) < (eps*eps*burgSize*burgSize))) {

                        vec1dotb = DotProduct(vec1, burgLab);
                        vec2dotb = DotProduct(vec2, burgLab);

                        if (!pinned2) {
/*
 *                          Neighbor 2 can be moved, so proceed with the
 *                          cross-slip operation.
 */
                            nbr2p[0] = nbr1p[0] - (vec1dotb * burgLab[0]);
                            nbr2p[1] = nbr1p[1] - (vec1dotb * burgLab[1]);
                            nbr2p[2] = nbr1p[2] - (vec1dotb * burgLab[2]);

                            NodeLocation(home, nbr2p[0], nbr2p[1], nbr2p[2], curInGrain, curOnGB);
                            if(nbr2->inGrain != curInGrain || nbr2->onGB != curOnGB){
                                continue;
                            }
                        }

/*
 *                      If neighbor2 is pinned, it is already perfectly
 *                      aligned with neighbor1 in the screw direction
 *                      so there is no need to move it.
 */
                        nodep[0] = nbr1p[0] + (vec2dotb * burgLab[0]);
                        nodep[1] = nbr1p[1] + (vec2dotb * burgLab[1]);
                        nodep[2] = nbr1p[2] + (vec2dotb * burgLab[2]);

                        NodeLocation(home, nodep[0], nodep[1], nodep[2], curInGrain, curOnGB);
                        if(node->inGrain != curInGrain || node->onGB != curOnGB){
                            continue;
                        }

                        fdotglide = DotProduct(fLab, glideDirLab[fplane]);
                        tmp = areamin * lamp / fabs(vec1dotb) * 2.0 *
                               (1.0 + eps) * Sign(fdotglide);

                        VECTOR_COPY(temVec, nodep);
                        for(iter=iterMax;iter>=0; iter--){
                            cof = real8(iter)/real8(iterMax);
                            CoordAdd(1.0, temVec, tmp*cof, glideDirLab[fplane], nodep);
                            NodeLocation(home, nodep[0], nodep[1], nodep[2], curInGrain, curOnGB);
                            if(node->inGrain == curInGrain && node->onGB == curOnGB){
                                break;
                            }
                        }

/*
 *                      It looks like we should do the cross-slip, but to
 *                      be sure, we need to move the nodes and evaluate
 *                      the force on node in the new configuration.  If
 *                      it appears the node will not continue to move out
 *                      on the new plane, skip the cross-slip event and
 *                      restore the old configuration.
 */
                        ResetPosition(param, node, nodep);
                        ResetPosition(param, nbr2, nbr2p);

                        SetOneNodeForce(home, node);

                        newForce[0] = node->fX;
                        newForce[1] = node->fY;
                        newForce[2] = node->fZ;

                        newfdotglide = DotProduct(newForce,glideDirLab[fplane]);

                        if ((Sign(newfdotglide) * Sign(fdotglide)) < 0.0) {
                            ResetPosition(param, node, nodePosOrig);
                            ResetPosition(param, nbr2, nbr2PosOrig);
                            RestoreCrossSlipForce(node, nbr1, nbr2, nbr1ArmID,
                                                  nbr2ArmID, segForceOrig);
                            continue;
                        }

                        resetNbr2Pos = 1;
                        resetNodePos = 1;

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "2nd neighbor is movable");
                        }
#endif
/*
 *                      Now reset the glide plane for both segments
 */
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr1->myTag, 1);
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr2->myTag, 1);
                    }
                } else {
/*
 *                  Neighbor 1 can be moved, so proceed with the cross-slip
 *                  operation.
 */
                        vec1dotb = DotProduct(vec1, burgLab);

                        nbr1p[0] = nbr2p[0] + (vec1dotb * burgLab[0]);
                        nbr1p[1] = nbr2p[1] + (vec1dotb * burgLab[1]);
                        nbr1p[2] = nbr2p[2] + (vec1dotb * burgLab[2]);

                        NodeLocation(home, nbr1p[0], nbr1p[1], nbr1p[2], curInGrain, curOnGB);
                        if(curInGrain != nbr1->inGrain || curOnGB != nbr1->onGB){
                            continue;
                        }

                        vec3dotb = DotProduct(vec3, burgLab);

                        nodep[0] = nbr2p[0] + (vec3dotb * burgLab[0]);
                        nodep[1] = nbr2p[1] + (vec3dotb * burgLab[1]);
                        nodep[2] = nbr2p[2] + (vec3dotb * burgLab[2]);

                        NodeLocation(home, nodep[0], nodep[1], nodep[2], curInGrain, curOnGB);
                        if(node->inGrain != curInGrain || node->onGB != curOnGB){
                            continue;
                        }

                        fdotglide = DotProduct(fLab, glideDirLab[fplane]);
                        tmp = areamin * lamp / fabs(vec1dotb) * 2.0 *
                               (1.0 + eps) * Sign(fdotglide);

                        VECTOR_COPY(temVec, nodep);
                        for(iter=iterMax;iter>=0; iter--){
                            cof = real8(iter)/real8(iterMax);
                            CoordAdd(1.0, temVec, tmp*cof, glideDirLab[fplane], nodep);
                            NodeLocation(home, nodep[0], nodep[1], nodep[2], curInGrain, curOnGB);
                            if(node->inGrain == curInGrain && node->onGB == curOnGB){
                                break;
                            }
                        }
/*
 *                      It looks like we should do the cross-slip, but to
 *                      be sure, we need to move the nodes and evaluate
 *                      the force on node in the new configuration.  If
 *                      it appears the node will not continue to move out
 *                      on the new plane, skip the cross-slip event and
 *                      restore the old configuration.
 */
                        ResetPosition(param, node, nodep);
                        ResetPosition(param, nbr1, nbr1p);

                        SetOneNodeForce(home, node);

                        newForce[0] = node->fX;
                        newForce[1] = node->fY;
                        newForce[2] = node->fZ;

                        newfdotglide = DotProduct(newForce,glideDirLab[fplane]);

                        if ((Sign(newfdotglide) * Sign(fdotglide)) < 0.0) {
                            ResetPosition(param, node, nodePosOrig);
                            ResetPosition(param, nbr1, nbr1PosOrig);
                            RestoreCrossSlipForce(node, nbr1, nbr2, nbr1ArmID,
                                                  nbr2ArmID, segForceOrig);
                            continue;
                        }

                        resetNbr1Pos = 1;
                        resetNodePos = 1;

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "1st neighbor is movable");
                        }
#endif
/*
 *                      Now reset the glide plane for both segments
 */
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr1->myTag, 1);
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr2->myTag, 1);
                }

                } else if ((seg1_is_screw) &&
                           (plane1 != plane2) &&
                           (plane2 == fplane) &&
                           (fabs(tmp3C[fplane]) >
                            (weightFactor*fabs(tmp3C[plane1])+fnodeThreshold))){

/*
 *                  Zipper condition met for first segment.  If the first
 *                  neighbor is either not pinned or pinned but already
 *                  sufficiently aligned, proceed with the cross-slip event
 */
                    pinned1 = NodePinned(home, nbr1, plane1, glideDirLab);

                    if ((!pinned1) ||
                        ((testmax2-test2) < (eps*eps*burgSize*burgSize))) {
/*
 *                      Before 'zippering' a segment, try a quick sanity check
 *                      to see if it makes sense.  If the force on the segment
 *                      to be 'zippered' is less than 5% larger on the new
 *                      plane than the old plane, leave the segment alone.
 */
                        newSegForce[0] = node->armfx[0]+nbr1->armfx[nbr1ArmID];
                        newSegForce[1] = node->armfy[0]+nbr1->armfy[nbr1ArmID];
                        newSegForce[2] = node->armfz[0]+nbr1->armfz[nbr1ArmID];

                        zipperThreshold = noiseFactor * shearModulus *
                                          burgSize *  L1;

                        f1dotplane1 = fabs(DotProduct(newSegForce,
                                                      glideDirLab[plane1]));
                        f1dotplanef = fabs(DotProduct(newSegForce, newplane));

                        if (f1dotplanef < zipperThreshold + f1dotplane1) {
                            continue;
                        }

                        if (!pinned1) {
                            vec2dotb = DotProduct(vec2, burgLab);

                            nbr1p[0] = nodep[0] - (vec2dotb * burgLab[0]);
                            nbr1p[1] = nodep[1] - (vec2dotb * burgLab[1]);
                            nbr1p[2] = nodep[2] - (vec2dotb * burgLab[2]);

                            NodeLocation(home, nbr1p[0], nbr1p[1], nbr1p[2], curInGrain, curOnGB);
                            if(curInGrain != nbr1->inGrain || curOnGB != nbr1->onGB){
                                continue;
                            }

                            resetNbr1Pos = 1;
                        }

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "zipper conditions met for 1st neighbor");
                        }
#endif
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr1->myTag, 1);
                    }

                } else if ((seg2_is_screw) &&
                           (plane1 != plane2) &&
                           (plane1 == fplane) &&
                           (fabs(tmp3C[fplane]) >
                            (weightFactor*fabs(tmp3C[plane2])+fnodeThreshold))){
/*
 *                  Zipper condition met for second segment
 */
                    pinned2 = NodePinned(home, nbr2, plane2, glideDirLab);

                    if ((!pinned2) ||
                        ((testmax2-test2) < (eps*eps*burgSize*burgSize))) {
/*
 *                      Before 'zippering' a segment, try a quick sanity check
 *                      to see if it makes sense.  If the force on the segment
 *                      to be 'zippered' is less than 5% larger on the new
 *                      plane than the old plane, leave the segment alone.
 */
                        newSegForce[0] = node->armfx[1]+nbr2->armfx[nbr2ArmID];
                        newSegForce[1] = node->armfy[1]+nbr2->armfy[nbr2ArmID];
                        newSegForce[2] = node->armfz[1]+nbr2->armfz[nbr2ArmID];

                        zipperThreshold = noiseFactor * shearModulus *
                                          burgSize *  L2;

                        f1dotplane2 = fabs(DotProduct(newSegForce,
                                                      glideDirLab[plane2]));
                        f1dotplanef = fabs(DotProduct(newSegForce, newplane));

                        if (f1dotplanef < zipperThreshold + f1dotplane2) {
                            continue;
                        }

                        if (!pinned2) {
                            vec3dotb = DotProduct(vec3, burgLab);

                            nbr2p[0] = nodep[0] - (vec3dotb * burgLab[0]);
                            nbr2p[1] = nodep[1] - (vec3dotb * burgLab[1]);
                            nbr2p[2] = nodep[2] - (vec3dotb * burgLab[2]);

                            NodeLocation(home, nbr2p[0], nbr2p[1], nbr2p[2], curInGrain, curOnGB);
                            if(nbr2->inGrain != curInGrain || nbr2->onGB != curOnGB){
                                continue;
                            }
                            resetNbr2Pos = 1;
                        }

#ifdef DEBUG_CROSSSLIP_EVENTS
                        if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                            DumpCrossSlipEvent(node, newplane,
                                    "zipper conditions met for 2nd neighbor");
                        }
#endif
                        ResetGlidePlane(home, newplane, &node->myTag,
                                        &nbr2->myTag, 1);
                    }
                }
/*
 *              If any of the nodes were repositioned, shift the new coordinates
 *              back into the primary image space and update the corresponding
 *              node structures.  The operation will also be sent to the
 *              remote domains for processing.
 */
                if (resetNodePos) {
                    FoldBox(param, &nodep[0], &nodep[1], &nodep[2]);
                    RepositionNode(home, nodep, &node->myTag, 1);
                    NodeLocation(home, nodep[0], nodep[1], nodep[2], curInGrain, curOnGB);
                    if(curInGrain != node->inGrain || curOnGB != node->onGB){
                        PrintNode(home, node);
                        Fatal("%s: inGrain=%d onGB=%d",__func__, curInGrain, curOnGB);
                    }
                }

                if (resetNbr1Pos) {
                    FoldBox(param, &nbr1p[0], &nbr1p[1], &nbr1p[2]);
                    RepositionNode(home, nbr1p, &nbr1->myTag, 1);
                    NodeLocation(home, nbr1p[0], nbr1p[1], nbr1p[2], curInGrain, curOnGB);
                    if(curInGrain != nbr1->inGrain || curOnGB != nbr1->onGB){
                        PrintNode(home, nbr1);
                        Fatal("%s: inGrain=%d onGB=%d",__func__, curInGrain, curOnGB);
                    }
                }

                if (resetNbr2Pos) {
                    FoldBox(param, &nbr2p[0], &nbr2p[1], &nbr2p[2]);
                    RepositionNode(home, nbr2p, &nbr2->myTag, 1);
                    NodeLocation(home, nbr2p[0], nbr2p[1], nbr2p[2], curInGrain, curOnGB);
                    if(curInGrain != nbr2->inGrain || curOnGB != nbr2->onGB){
                        PrintNode(home, nbr2);
                        Fatal("%s: inGrain=%d onGB=%d",__func__, curInGrain, curOnGB);
                    }
                }

            }  /* end of screw check */

        }  /* end loop over all nodes */

        return;
}
