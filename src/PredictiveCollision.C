/*****************************************************************************
 *
 *      Module:         PredictiveCollision.c
 *      Description:    This module contains various functions used
 *                      for detecting various types of collisions
 *                      (segment/segment, node/node, zipping) and
 *                      dealing with those collisions.  These functions
 *                      are specific to the type 2 collision handling
 *                      which attempts to predict if there is a future
 *                      time at which segments will physicaly intersect,
 *                      and uses this data as collision criteria rather
 *                      than the simple proximity criteria of the old
 *                      mechanism
 *
 *                      Each domain handles local collisions and
 *                      distributes the necessary topological changes
 *                      to remote domains via the same mechanism
 *                      employed by remesh.
 *
 *                      NOTE: Certain strict rules govern what topological
 *                      changes are permitted based on noda and segment
 *                      ownership.  See comments at the beginning of
 *                      the module Topology.c for details of the
 *                      rule.  Additional restrictions may be implemented
 *                      for the collision handling; see code below
 *                      for details of any such restrictions.
 *
 *      Included functions:
 *
 *          FindCollisionPoint()
 *          FindCollisionPointAndTime()
 *          PredictiveCollisions()
 *
 *****************************************************************************/
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"
#include "Math.h"

static int dbgDom;


/*---------------------------------------------------------------------------
 *
 *      Function:       FindCollisionPoint
 *      Description:    This function attempts to select a collision
 *                      point on a plane common to the two nodes.
 *
 *      Arguments:
 *          node1, node2   pointers to the two node structures
 *          x, y, z        pointers to locations in which to return
 *                         the coordinates of the point at which the
 *                         two nodes should be collided.
 *
 *-------------------------------------------------------------------------*/
static void FindCollisionPoint(Home_t *home, Node_t *node1, Node_t *node2,
                               real8 *x, real8 *y, real8 *z)
{
        int     i, j, m, n;
        int     conditionsmet, Nsize;
        int     planeDefined;
        real8   L, invL, tmp;
        real8   norm, invnorm;
        real8   n1mag2, n2mag2, eps;
        real8   dx, dy, dz;
        real8   n1x, n1y, n1z;
        real8   n2x, n2y, n2z;
        real8   dirx, diry, dirz;
        real8   p1[3], p2[3];
        real8   plane[3], vector[3];
        real8   newplanecond, npc2, onemnpc4, detN;
        real8   Nmat[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
        real8   Matrix[6][6], invMatrix[6][6];
        real8   V[6], result[6];
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

        eps = 1.0e-12;

        p1[0] = node1->x;
        p1[1] = node1->y;
        p1[2] = node1->z;

        p2[0] = node2->x;
        p2[1] = node2->y;
        p2[2] = node2->z;

        PBCPOSITION(param, p1[0], p1[1], p1[2], &p2[0], &p2[1], &p2[2]);

/*
 *      If a node is a 'fixed' node it can't be relocated, so use
 *      that node's coordinates as the collision point
 */
        if (node1->constraint == PINNED_NODE) {
            *x = p1[0];
            *y = p1[1];
            *z = p1[2];
            return;
        } else if (node2->constraint == PINNED_NODE) {
            *x = p2[0];
            *y = p2[1];
            *z = p2[2];
            return;
        }


       newplanecond = 0.875;
       npc2         = newplanecond * newplanecond;

       tmp          = 1.0 - newplanecond;
       onemnpc4     = tmp * tmp * tmp * tmp;

       vector[0] = 0.0;
       vector[1] = 0.0;
       vector[2] = 0.0;

       Nsize = 0;

       for (i = 0; i < node1->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node1, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p1[0] - nbrNode->x;
               dy = p1[1] - nbrNode->y;
               dz = p1[2] - nbrNode->z;

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node1->burgX[i], node1->burgY[i],
                       node1->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node1->vX, node1->vY, node1->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } else if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node1->nx[i];
                   plane[1] = node1->ny[i];
                   plane[2] = node1->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }


               if (planeDefined) {

                   switch (Nsize) {
                   case 0:
                       conditionsmet = 1;
                       break;

                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;

                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];

                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }

                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p1);
                       Nsize++;
                   }
               }
           }
       }

       for (i = 0; i < node2->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node2, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p2[0] - nbrNode->x;
               dy = p2[1] - nbrNode->y;
               dz = p2[2] - nbrNode->z;

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node2->burgX[i], node2->burgY[i],
                       node2->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node2->vX, node2->vY, node2->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } else if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node2->nx[i];
                   plane[1] = node2->ny[i];
                   plane[2] = node2->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }

               if ((n1mag2 > eps) || (n2mag2 > eps)) {
                   switch (Nsize) {
                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;
                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];
                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }
 
                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p2);
                       Nsize++;
                   }
               }
           }
       }

/*
 *         Upper left 3X3 of Matrix is identity matrix.
 *         Matrix rows 3 thru 3+(Nsize-1) colums 0 thru 2 are Nmat.
 *         Matrix columns 3 thru 3+(Nsize-1) rows 0 thru 2 are transpose of Nmat.
 *         All remaining elements are zeroed.
 */
       for (i = 0; i < 6; i++) {
           for (j = 0; j < 6; j++) {
               Matrix[i][j] = 0.0;
           }
       }

       Matrix[0][0] = 1.0;
       Matrix[1][1] = 1.0;
       Matrix[2][2] = 1.0;

       for (i = 0; i < Nsize; i++) {
           for (j = 0; j < 3; j++) {
               Matrix[3+i][j] = Nmat[i][j];
               Matrix[j][3+i] = Nmat[i][j];
           }
       }

       V[0] = 0.5 * (p1[0] + p2[0]);
       V[1] = 0.5 * (p1[1] + p2[1]);
       V[2] = 0.5 * (p1[2] + p2[2]);
       V[3] = vector[0];
       V[4] = vector[1];
       V[5] = vector[2];

       Nsize += 3;

       MatrixInvert((real8 *)Matrix, (real8 *)invMatrix, Nsize, 6);
       MatrixMult((real8 *)invMatrix, Nsize, Nsize, 6,
                  (real8 *)V, 1, 1,
                  (real8 *)result, 1);

       *x = result[0];
       *y = result[1];
       *z = result[2];

       return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FindCollisionPointAndTime
 *      Description:    Given two segments and their velocities, determine
 *                      if the segments will physically intersect at any
 *                      given time in the future.  If so, return to the
 *                      caller the position at which the collision will 
 *                      occur, the time (in units of current deltaTT)
 *                      at which the collision would occur, etc.
 *                      
 *
 *      Arguments:
 *          p1,p2       Endpoint nodes of segment 1
 *          p3,p4       Endpoint nodes of segment 2
 *          v1,v2 v3,v4 Velocity of nodes p1 thru p4 respectively
 *          cPoint      Will be set to the coordinates of the collision
 *                      point if the two segments will collide.
 *          cTime       Set to the time (in units of current timestep
 *                      length (deltaTT)) at which the two segments will
 *                      collide, or negative if not colliding.
 *          minDist2    Min distance between segments, squared.
 *          L1, L2
 *
 *-------------------------------------------------------------------------*/
static void FindCollisionPointAndTime(Home_t *home, real8 p1[3], real8 p2[3],
                                  real8 p3[3], real8 p4[3], real8 v1[3],
                                  real8 v2[3], real8 v3[3], real8 v4[3],
                                  real8 cPoint[3], real8 *cTime,
                                  real8 *minDist2, real8 *L1, real8 *L2)
{
        int   i;
        real8 A, B, C, D, E, F, G, H, J;
        real8 dt, matDet, eps;
        real8 L1a, L2a, L3, L3a;
        real8 dist2, ddist2dt, minDist2Tol;
        real8 tmpMat22[2][2], tmpMat22Inv[2][2];
        real8 tmpVec2a[2], tmpVec2b[2];
        real8 tmpVec3[3];
        real8 R0[3], seg1[3], seg2[3], vDiff[3];
        real8 uBar10[3], vBar10[3];
        real8 rhs[3];
        real8 mat[3][3], matInv[3][3];


        eps = 1.0e-12;
        dt = home->param->deltaTT;

        for (i = 0; i < 3; i++) {
            R0[i] = p1[i] - p3[i];
            seg1[i] = p2[i] - p1[i];
            seg2[i] = p4[i] - p3[i];
            vBar10[i] = (v1[i] + v2[i]) * 5.0 * dt;
            uBar10[i] = (v3[i] + v4[i]) * 5.0 * dt;
            vDiff[i] = vBar10[i] - uBar10[i];
        }

        A =   DotProduct(seg1, seg1);
        B = -(DotProduct(seg1, seg2));
        C =   DotProduct(seg1, vDiff);
        D =   DotProduct(seg2, seg2);
        E = -(DotProduct(seg2, vDiff));
        F =   DotProduct(vDiff, vDiff);
        G =   DotProduct(R0, seg1);
        H = -(DotProduct(R0, seg2));
        J =   DotProduct(R0, vDiff);

        mat[0][0] = A;  mat[0][1] = C;  mat[0][2] = B;
        mat[1][0] = C;  mat[1][1] = F;  mat[1][2] = E;
        mat[2][0] = B;  mat[2][1] = E;  mat[2][2] = D;

        rhs[0] = -G;
        rhs[1] = -J;
        rhs[2] = -H;

/*
 *      Check the solution of all the faces, edges and corners.
 */
        L3 = 0.0;  /* Start with zero time */

        GetMinDist(p1[0], p1[1], p1[2], v1[0], v1[1], v1[2],
                   p2[0], p2[1], p2[2], v2[0], v2[1], v2[2],
                   p3[0], p3[1], p3[2], v3[0], v3[1], v3[2],
                   p4[0], p4[1], p4[2], v4[0], v4[1], v4[2],
                   minDist2, &ddist2dt, L1, L2);

/*
 *      Check the end node pairs of the segments
 */
        if (F > eps) {

            L1a = 0.0;  /* Check for this node 1 on seg1 */
            L2a = 0.0;  /* Check for this node 1 on seg2 */
            L3a = J / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);
              
                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 1.0;  /* Check for this node 2 on seg1 */
            L2a = 0.0;  /* Check for this node 1 on seg2 */
            L3a = (J + C) / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + seg1[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 0.0;  /* Check for this node 1 on seg1 */
            L2a = 1.0;  /* Check for this node 2 on seg2 */
            L3a = (J + E) / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] - seg2[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 1.0;  /* Check for this node 2 on seg1 */
            L2a = 1.0;  /* Check for this node 2 on seg2 */
            L3a = (J + C + E) / -F;

            if (L3a > 0.0) {

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + seg1[i] - seg2[i] + vDiff[i] * L3a;
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

        } /* end if (F > eps) */

/*
 *      Check single end nodes on seg1
 */
        if (fabs(D*F - E*E) > eps) {
            L1a = 0.0;
            tmpMat22[0][0] = mat[1][1];
            tmpMat22[0][1] = mat[1][2];
            tmpMat22[1][0] = mat[2][1];
            tmpMat22[1][1] = mat[2][2];

            Matrix22Invert(tmpMat22, tmpMat22Inv);
            Matrix22Vector2Mult(tmpMat22Inv, &rhs[1], tmpVec2a);

            if ((tmpVec2a[1] >= 0.0) &&
                (tmpVec2a[1] <= 1.0) &&
                (tmpVec2a[0] >= 0.0)) {

                L2a = tmpVec2a[1];
                L3a = tmpVec2a[0];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] - (seg2[i] * L2a) + (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

            L1a = 1.0;

            tmpVec2a[0] = rhs[1] - mat[1][0];
            tmpVec2a[1] = rhs[2] - mat[2][0];
 
            /* note: tmpMat22Inv calculated above */
            Matrix22Vector2Mult(tmpMat22Inv, tmpVec2a, tmpVec2b);

            if ((tmpVec2b[1] >= 0.0) &&
                (tmpVec2b[1] <= 1.0) &&
                (tmpVec2b[0] >= 0.0)) {

                L2a = tmpVec2b[1];
                L3a = tmpVec2b[0];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] - (seg2[i] * L2a) + (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

        } /* end if (fabs(D*F - E*E) > eps) */

/*
 *      Check single end nodes on seg2
 */
        if (fabs(A*F - C*C) > eps) {

            L2a = 0.0;

            tmpMat22[0][0] = mat[0][0];
            tmpMat22[0][1] = mat[0][1];
            tmpMat22[1][0] = mat[1][0];
            tmpMat22[1][1] = mat[1][1];

            Matrix22Invert(tmpMat22, tmpMat22Inv);
            Matrix22Vector2Mult(tmpMat22Inv, rhs, tmpVec2a);

            if ((tmpVec2a[0] >= 0.0) &&
                (tmpVec2a[0] <= 1.0) &&
                (tmpVec2a[1] >= 0.0)) {

                L1a = tmpVec2a[0];
                L3a = tmpVec2a[1];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + (seg1[i] * L1a) + (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }
            
            L2a = 1.0;

            tmpVec2a[0] = rhs[0] - mat[0][2];
            tmpVec2a[1] = rhs[1] - mat[1][2];
 
            /* note: tmpMat22Inv calculated above */
            Matrix22Vector2Mult(tmpMat22Inv, tmpVec2a, tmpVec2b);

            if ((tmpVec2b[0] >= 0.0) &&
                (tmpVec2b[0] <= 1.0) &&
                (tmpVec2b[1] >= 0.0)) {

                L1a = tmpVec2b[0];
                L3a = tmpVec2b[1];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + (seg1[i] * L1a) - seg2[i] +
                                 (vDiff[i] * L3a);
                }

                dist2 = DotProduct(tmpVec3, tmpVec3);

                if (dist2 < *minDist2) {
                    *L1 = L1a;
                    *L2 = L2a;
                    L3 = L3a;
                    *minDist2 = dist2;
                }
            }

        }  /* end if (fabs(A*F - C*C) > eps) */

        matDet = Matrix33Det(mat);

        if (fabs(matDet) > eps) {

            if (!Matrix33Invert(mat, matInv)) {
#if 0
                Fatal("Unable to invert matrix in FindCollisionPointAndTime");
#endif
                *L1 = 0.0;
                *L2 = 0.0;
                *cTime = -1;

                return;
            }

            Matrix33Vector3Multiply(matInv, rhs, tmpVec3);

            if ((tmpVec3[0] >= 0.0) && (tmpVec3[0] <= 1.0) &&
                (tmpVec3[2] >= 0.0) && (tmpVec3[2] <= 1.0) &&
                (tmpVec3[1] >= 0.0)) {
/*
 *              An interior point has been found
 */
                *L1 = tmpVec3[0];
                *L2 = tmpVec3[2];
                L3 = tmpVec3[1];

                for (i = 0; i < 3; i++) {
                    tmpVec3[i] = R0[i] + (seg1[i] * *L1) - (seg2[i] * *L2) +
                                 (vDiff[i] * L3);
                }

                *minDist2 = DotProduct(tmpVec3, tmpVec3);
            }

        } /* end if (fabs(matDet) > eps) */


/*
 *      In general, if the dislocation segments on their current
 *      trajectories will not eventually intersect we don't want
 *      to do a collision and set the collision time negative to
 *      indicate no collision.  First, determine the tolerance for
 *      deciding if the segments will be within collision distance in
 *      the future.  Use a tight tolerance if strict glide planes
 *      are enforced, loosen up the tolerance a little in all
 *      other cases.
 */
        if (home->param->enforceGlidePlanes &&
            (home->param->allowFuzzyGlidePlanes == 0)) {
            minDist2Tol = eps;
        } else {
            minDist2Tol = 1.0;
        }


/*
 *      If collision criteria is not met, set collision time to a
 *      negative value.  Otherwise calculate the collision point and
 *      time.
 */
        if (*minDist2 > minDist2Tol) {

            *cTime = -1.0;

        } else {

            *cTime = 10.0 * L3;

            for (i = 0; i < 3; i++) {
                cPoint[i] = 0.5 * ((p1[i] + p3[i]) +
                                   (*L1 * seg1[i]) +
                                   (*L2 * seg2[i]) +
                                   (vBar10[i] + uBar10[i]) * L3);
            }
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       PredictiveCollisions
 *      Description:    Loop though all native nodes, identify segments
 *                      or nodes that should be collided (or zipped) by
 *                      the current domain and handle all such collisions.
 *
 *                      Note the following restrictions on collisions:
 *
 *                      - A 'pinned' node be deleted during a collision
 *                      - A node which has been the target of a merge
 *                        (this timestep) may not be deleted during any 
 *                        other type of collision in this timestep.
 *
 *      All nodes with arms owned by another domain have previously
 *      been marked as non-deletable nodes for this cycle.
 *
 *-------------------------------------------------------------------------*/

inline bool ParallelGB(Node_t *node, real8 n[3]){
    for(int i=0; i<node->numNbrs; i++){
        if(fabs(dot_product(n[0],n[1],n[2],node->nx[i],node->ny[i],node->nz[i])) > 0.95)return true;
    }    
    return false;
}

void PredictiveCollisions(Home_t *home)
{
        int     i, j, k, q, arm12, arm21, arm34, arm43;
        int     thisDomain, splitStatus, mergeStatus;
        int     armAB, armBA;
        int     globalOp = 1, didCollision;
        int     close2node1, close2node2, close2node3, close2node4;
        int     splitSeg1, splitSeg2;
        int     cell2Index, nbrCell2Index, nextIndex;
        int     cell2X, cell2Y, cell2Z, cx, cy, cz;
        int     localCollisionCnt, globalCollisionCnt;
        int     collisionConditionIsMet, adjustCollisionPoint;
        int     inGrain1=-1, inGrain2=-1, onGB1=-1, onGB2=-1;
        real8   mindist2, dist2, ddist2dt, L1, L2, eps, half;
        real8   cTime, cPoint[3];
        real8   x1, y1, z1, vx1, vy1, vz1;
        real8   x2, y2, z2, vx2, vy2, vz2;
        real8   x3, y3, z3, vx3, vy3, vz3;
        real8   x4, y4, z4, vx4, vy4, vz4;
        real8   seg1Lx, seg1Ly, seg1Lz;
        real8   seg2Lx, seg2Ly, seg2Lz;
        real8   newx, newy, newz;
        real8   newvx, newvy, newvz;
        real8   pnew[3], burg1[3], burg2[3];
        real8   oldfp0[3], oldfp1[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        real8   nodeVel[3], newNodeVel[3];
        real8   newPos[3], newVel[3];
        real8   vec1[3], vec2[3], vec3[3];
        real8   p1[3], p2[3], p3[3], p4[3];
        real8   v1[3], v2[3], v3[3], v4[3];
        real8   pos[3];
        Tag_t   oldTag1, oldTag2;
        Node_t  *node1, *node2, *node3, *node4, *tmpNbr;
        Node_t  *mergenode1, *mergenode2, *targetNode;
        Node_t  *splitNode1, *splitNode2;
        GB_t    *gb;
        Grain_t *grain;
        int     gbIndex;
        Plane_t plane;
        Param_t *param;
        int newInGrain, newOnGB, type; real8 dis; bool canMerge;

        int     gbCollision;
        int     splitInGrain, splitOnGB;

#ifdef _FEM
        int     resetSurfaceProperties, femSurface[2];
        real8   femSurfaceNorm[3];
#endif

        thisDomain = home->myDomain;
        param      = home->param;

/* QUESTION!  What is reasonable for mindist2? */
        mindist2 = param->rann * param->rann;
        eps      = 1.0e-12;
        half     = 0.5;

        localCollisionCnt = 0;
        globalCollisionCnt = 0;

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        TimerStart(home, COLLISION_HANDLING);

        DebugGB(home, 401);
/*
 *      Start looping through native nodes looking for segments to collide...
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            didCollision = 0;

/*
 *          Loop through all cell2s neighboring the node.  Only
 *          nodes in these neighboring cell2s are candidates for
 *          collisions.
 *
 *          NOTE: All nodes are assigned cell2 membership prior to
 *          entering collision handling.  However, nodes added during
 *          node separation and collision handling are not given
 *          cell2 membership for this cycle.  In most cases this
 *          is not a problem since those nodes have usually been
 *          exempted from subsequent collisions this cycle.  There
 *          are certain circumstances, though, in which the exemptions
 *          have not been set.  If we encounter one such node that
 *          has no cell2 membership, skip it for this cycle, or bad
 *          things happen.
 */
            cell2Index = node1->cell2Idx;
            if (cell2Index < 0) {
                continue;
            }

            DecodeCell2Idx(home, cell2Index, &cell2X, &cell2Y, &cell2Z);

            for (cx = cell2X - 1; cx <= cell2X + 1; cx++) {
             for (cy = cell2Y - 1; cy <= cell2Y + 1; cy++) {
              for (cz = cell2Z - 1; cz <= cell2Z + 1; cz++) {

                nbrCell2Index = EncodeCell2Idx(home, cx, cy, cz);

/*
 *              Loop though all nodes in the neighbor cell2
 */
                nextIndex = home->cell2[nbrCell2Index];

                while (nextIndex >= 0) {

                    if (didCollision) break;

                    node3 = home->cell2QentArray[nextIndex].node;
                    nextIndex = home->cell2QentArray[nextIndex].next;

                    if (node3 == (Node_t *)NULL) continue;
                    if (node3->flags & NO_COLLISIONS) continue;

                    if (CollisionNodeOrder(home, &node1->myTag,
                                           &node3->myTag) >= 0) {
                        continue;
                    }

/*
 *                  Loop over all arms of node1.  Skip any arms that
 *                  terminate at node3 (those hinge arms will be dealt
 *                  with later)
 */
                    for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

                        if (didCollision) break;

                        node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

                        if (node2 == (Node_t *)NULL) continue;
                        if (node2->flags & NO_COLLISIONS) continue;

                        if (CollisionNodeOrder(home, &node1->myTag,
                                               &node2->myTag) > 0) {
                            continue;
                        }

                        if ((node2->myTag.domainID == node3->myTag.domainID) &&
                            (node2->myTag.index    == node3->myTag.index   )) {
                            continue;
                        }

/*
 *                      Segment node1/node2 may only be used in a collision
 *                      if the segment is owned by the current domain.
 */
                        if (!DomainOwnsSeg(home, OPCLASS_COLLISION,
                                           thisDomain, &node2->myTag)) {
                            continue;
                        }

                        if(SurfaceNode(node1) && SurfaceNode(node2))continue;

/*
 *                      Loop over all arms of node3.  Skip any segments that
 *                      terminate at node2 (those hinge arms will be dealt
 *                      with later), and skip any segments not owned by
 *                      node3 -- we'll deal with those segments when we
 *                      hit the node owning the segment
 *                      
 */
                        for (arm34 = 0; arm34 < node3->numNbrs; arm34++) {
                            real8 dx, dy, dz;

                            if (didCollision) break;

                            node4 = GetNodeFromTag(home, node3->nbrTag[arm34]);

                            if (node4 == (Node_t *)NULL) continue;
                            if (node4->flags & NO_COLLISIONS) continue;

                            if ((node4->myTag.domainID==node2->myTag.domainID)&&
                                (node4->myTag.index   ==node2->myTag.index)) {
                                continue;
                            }

                            if(SurfaceNode(node3) && SurfaceNode(node3))continue;
            
                            if (CollisionNodeOrder(home, &node3->myTag,
                                                   &node4->myTag) > 0) {
                                continue;
                            }

                            gbCollision = 0;
                            if(home->param->polycrystal){
                                GetArmGrain(home, node1, node2, inGrain1, onGB1);
                                GetArmGrain(home, node3, node4, inGrain2, onGB2);

                                splitInGrain = -1;
                                splitOnGB = -1;
                                if(onGB1 > -1 || onGB2 > -1){
                                    if(onGB1 == onGB2){
                                        splitOnGB = onGB1;
                                    }else if(onGB1 > -1 && onGB2 > -1){
                                        continue;
                                    }else if(onGB1 > -1){
                                        splitOnGB = onGB1;
                                        gbCollision = 2;
                                    }else if(onGB2 > -1){
                                        splitOnGB = onGB2;
                                        gbCollision = 2;
                                    }else{
                                        continue;
                                    }
                                }else{
                                    if(inGrain1 == inGrain2){
                                        splitInGrain = inGrain1;
                                    }else{
                                        if(inGrain1 == -1 || inGrain2 == -1){
                                            Fatal("in %s at %d", __FILE__, __LINE__);
                                        }
                                        real8 refPos[3]={node1->x, node1->y, node1->z};
                                        if(GetPublicGB(home, home->grainKeys[inGrain1],
                                           home->grainKeys[inGrain2], &gb, refPos)){
                                           splitOnGB = gb->GBID;
                                            gbCollision = 1;
                                        }else{
                                            continue;
                                        }
                                    }
                                }
                            }else{
                                splitInGrain = -1;
                                splitOnGB = -1;
                            }
/*
 *                          At this point, segment node3/node4 is owned by
 *                          node3.  If node3 is not native to this domain,
 *                          the segment may not be used in a collision since
 *                          the domain doing to collision must own both
 *                          segments.
 */
                            if (node3->myTag.domainID != thisDomain) {
                                continue;
                            }

                            x1 = node1->x; y1 = node1->y; z1 = node1->z;
                            x2 = node2->x; y2 = node2->y; z2 = node2->z;
                            x3 = node3->x; y3 = node3->y; z3 = node3->z;
                            x4 = node4->x; y4 = node4->y; z4 = node4->z;
    
                            vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
                            vx2 = node2->vX; vy2 = node2->vY; vz2 = node2->vZ;
                            vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
                            vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;
    
                            PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
                            PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                            PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);
    
                            p1[0] = x1;  p1[1] = y1;  p1[2] = z1;
                            p2[0] = x2;  p2[1] = y2;  p2[2] = z2;
                            p3[0] = x3;  p3[1] = y3;  p3[2] = z3;
                            p4[0] = x4;  p4[1] = y4;  p4[2] = z4;

                            v1[0] = vx1;  v1[1] = vy1;  v1[2] = vz1;
                            v2[0] = vx2;  v2[1] = vy2;  v2[2] = vz2;
                            v3[0] = vx3;  v3[1] = vy3;  v3[2] = vz3;
                            v4[0] = vx4;  v4[1] = vy4;  v4[2] = vz4;

/*
 *                          It is possible to have a zero-length segment
 *                          (created by a previous collision).  If we find
 *                          such a segment, do not try to use it in any
 *                          subsequent collisions.
 */
                            vec1[0] = x2 - x1;
                            vec1[1] = y2 - y1;
                            vec1[2] = z2 - z1;

                            if (DotProduct(vec1, vec1) < 1.0e-20) {
                                continue;
                            }

                            vec2[0] = x4 - x3;
                            vec2[1] = y4 - y3;
                            vec2[2] = z4 - z3;

                            if (DotProduct(vec2, vec2) < 1.0e-20) {
                                continue;
                            }

/*
 *                          Find the minimum distance between the two segments
 *                          and determine if they should be collided.
 */
                            GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                                       x2, y2, z2, vx2, vy2, vz2,
                                       x3, y3, z3, vx3, vy3, vz3,
                                       x4, y4, z4, vx4, vy4, vz4,
                                       &dist2, &ddist2dt, &L1, &L2);
    
/*
 *                          First check of segments already intersect.  If
 *                          not, find out if they will collide in the future.
 *
 *                          Note: If the separation between the segments is
 *                          less than 1b, treat them as if they are already
 *                          intersecting
 */
                            if (dist2 < 1.0) {
                                cPoint[0] = x1 + vec1[0] * L1;
                                cPoint[1] = y1 + vec1[1] * L1;
                                cPoint[2] = z1 + vec1[2] * L1;
                                collisionConditionIsMet = 1;
                            } else if (dist2 < mindist2) {
/*
 *                              FIX ME!  If segments are to far away but
 *                              moving fast, they can pass right through
 *                              each other with no collision
 *
 *                              Only do a rigorous treatment of points 
 *                              within the distance filter.  Find the 
 *                              collision point, collision time and the
 *                              points on the segments where the two segments
 *                              will be colliding.
 */
                                FindCollisionPointAndTime(home, p1, p2, p3, p4,
                                                      v1, v2, v3, v4, cPoint,
                                                      &cTime, &dist2, &L1, &L2);

                                collisionConditionIsMet = ((cTime > 0.0) &&
                                                           (cTime < 10.0));
                            } else {
                                collisionConditionIsMet = 0;
                            }

                            if (collisionConditionIsMet) {
/*
 *                              Segments are unconnected and colliding.
 *                              Identify the first node to be merged.  If the
 *                              collision point is close to one of the nodal
 *                              endpoints, use that node, otherwise insert a
 *                              new node in the segment.
 *
 *                              NOTE: The current domain owns node1 but may
 *                              not own node2.  If it does not own node2, we 
 *                              cannot allow the collision to use node2
 *                              even if the collision point is close to 
 *                              that node.
 */
                                adjustCollisionPoint = 0;
                                vec1[0] = cPoint[0] - p1[0];
                                vec1[1] = cPoint[1] - p1[1];
                                vec1[2] = cPoint[2] - p1[2];
    
                                vec2[0] = cPoint[0] - p2[0];
                                vec2[1] = cPoint[1] - p2[1];
                                vec2[2] = cPoint[2] - p2[2];

                                close2node1 = (DotProduct(vec1,vec1)<mindist2);
                                close2node2 = (DotProduct(vec2,vec2)<mindist2);
    
                                if ((node2->myTag.domainID != thisDomain) &&
                                    close2node2) {
                                    continue;
                                }

                                if (close2node1) {
                                    mergenode1 = node1;
                                    splitSeg1 = 0;
                                    adjustCollisionPoint = 1;
                                } else if (close2node2) {
                                    mergenode1 = node2;
                                    splitSeg1 = 0;
                                    adjustCollisionPoint = 1;
                                } else {
                                    splitSeg1 = 1;
                                }
#if 0
                                if(gbCollision){
                                    printf("Task %d, gb collision %d on (%d,%d): (%d,%d,%d)--(%d,%d,%d) * (%d,%d,%d)--(%d,%d,%d)\n",
                                            home->myDomain, gbCollision, splitInGrain, splitOnGB,
                                            node1->myTag.domainID,node1->myTag.index,node1->inGrain,
                                            node2->myTag.domainID,node2->myTag.index,node2->inGrain,
                                            node3->myTag.domainID,node3->myTag.index,node3->inGrain,
                                            node4->myTag.domainID,node4->myTag.index,node4->inGrain);
                                }
#endif                       
/*
 *                              If we need to add a new node to the first
 *                              segment, do it now.
 */
                                if (splitSeg1) {
                                     real8 pos0[3], pos1[3];

                                     newx = x1 * (1.0-L1) + x2*L1;
                                     newy = y1 * (1.0-L1) + y2*L1;
                                     newz = z1 * (1.0-L1) + z2*L1;
    
                                     newvx = vx1 * (1.0-L1) + vx2*L1;
                                     newvy = vy1 * (1.0-L1) + vy2*L1;
                                     newvz = vz1 * (1.0-L1) + vz2*L1;
    
/*
 *                                   Estimate resulting forces on all segments
 *                                   involved in the split.
 */
                                     arm21 = GetArmID(home, node2, node1);
    
                                     oldfp0[0] = node1->armfx[arm12];
                                     oldfp0[1] = node1->armfy[arm12];
                                     oldfp0[2] = node1->armfz[arm12];
    
                                     oldfp1[0] = node2->armfx[arm21];
                                     oldfp1[1] = node2->armfy[arm21];
                                     oldfp1[2] = node2->armfz[arm21];
    
                                     pos0[0] = x1;   pos1[0] = x2;
                                     pos0[1] = y1;   pos1[1] = y2;
                                     pos0[2] = z1;   pos1[2] = z2;
    
                                     pnew[0] = newx;
                                     pnew[1] = newy;
                                     pnew[2] = newz;
    
                                     newNodeVel[0] = newvx;
                                     newNodeVel[1] = newvy;
                                     newNodeVel[2] = newvz;

                                     nodeVel[0] = vx1;
                                     nodeVel[1] = vy1;
                                     nodeVel[2] = vz1;

                                     burg1[0] = node1->burgX[arm12];
                                     burg1[1] = node1->burgY[arm12];
                                     burg1[2] = node1->burgZ[arm12];

                                     FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
                                                 oldfp1, pnew, f0seg1, f1seg1,
                                                 f0seg2, f1seg2);
    
                                     oldTag1 = node1->myTag;
                                     oldTag2 = node2->myTag;

                                     FoldBox(param, &pnew[0], &pnew[1], &pnew[2]);

                                     splitStatus = SplitNode(home,
                                                             OPCLASS_COLLISION,
                                                             node1, pos0, pnew,
                                                             nodeVel,
                                                             newNodeVel, 1,
                                                             &arm12, globalOp,
                                                             &splitNode1,
                                                             &splitNode2, 0
                                                             ,node1->inGrain, node1->onGB
                                                             ,inGrain1, onGB1
                                                                            );
/*
 *                                   If we were unable to split the node
 *                                   go back to looking for more collision
 *                                   candidates.
 */
                                     if (splitStatus == SPLIT_FAILED) {
                                         continue;
                                     }

/*
 *                                   The force estimates above are good enough
 *                                   for the remainder of this timestep, but
 *                                   mark the force and velocity data for
 *                                   some nodes as obsolete so that more
 *                                   accurate forces will be recalculated
 *                                   either at the end of this timestep, or
 *                                   the beginning of the next.
 */
                                     mergenode1 = splitNode2;

                                     MarkNodeForceObsolete(home, splitNode2);
    
                                     for (q = 0; q < splitNode2->numNbrs; q++) {
                                         tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
                                         if (tmpNbr != (Node_t *)NULL) {
                                             tmpNbr->flags |= NODE_RESET_FORCES;
                                         }
                                     }
    
/*
 *                                   Reset nodal forces on nodes involved in the
 *                                   split.
 */
                                     ResetSegForces(home, splitNode1,
                                                    &splitNode2->myTag,
                                                    f0seg1[0], f0seg1[1],
                                                    f0seg1[2], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &splitNode1->myTag,
                                                    f1seg1[0], f1seg1[1],
                                                    f1seg1[2], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &node2->myTag,
                                                    f0seg2[0], f0seg2[1],
                                                    f0seg2[2], 1);
    
                                     ResetSegForces(home, node2,
                                                    &splitNode2->myTag,
                                                    f1seg2[0], f1seg2[1],
                                                    f1seg2[2], 1);
    
                                     (void)EvaluateMobility(home, splitNode1);
                                     (void)EvaluateMobility(home, splitNode2);
                                     (void)EvaluateMobility(home, node2);

/*
 *                                  When debugging, dump some info on
 *                                  topological changes taking place and
 *                                  the nodes involved
 */
#if 0
                                    if(gbCollision){
                                        printf("%d  Split-1(SegCollision): "
                                               "(%d,%d)--(%d,%d) ==> "
                                               "(%d,%d,%d)--(%d,%d,%d)--(%d,%d,%d)\n",home->myDomain,
                                               oldTag1.domainID, oldTag1.index,
                                               oldTag2.domainID, oldTag2.index,
                                               splitNode1->myTag.domainID,
                                               splitNode1->myTag.index,
                                               splitNode1->inGrain,
                                               splitNode2->myTag.domainID,
                                               splitNode2->myTag.index,
                                               splitNode2->inGrain,
                                               node2->myTag.domainID,
                                               node2->myTag.index, node2->inGrain);
                                     }
#endif
/*
 *                                   When we actually do a collision, we flag
 *                                   it so we don't attempt additional changes
 *                                   on either segment.  It's possible, however
 *                                   that we will split the segment but later
 *                                   determine a collision will not be done.
 *                                   Since the original node1/node2 segment
 *                                   has been bisected, we must treat it as
 *                                   if a collision had occurred so we do not
 *                                   attempt to collide the now non-existent
 *                                   node1/node2 segment.
 */
                                     didCollision = 1;

                                }          /* if (splitSeg1) */
/*
 *                              Identify the second node to be merged
 *
 *                              Note: The current domain owns node3 but may not
 *                              own node4.  If it does not own node4, we 
 *                              cannot allow the collision to use node4
 *                              even if the collision point is close to 
 *                              that node.
 */
                                vec1[0] = cPoint[0] - p3[0];
                                vec1[1] = cPoint[1] - p3[1];
                                vec1[2] = cPoint[2] - p3[2];
    
                                vec2[0] = cPoint[0] - p4[0];
                                vec2[1] = cPoint[1] - p4[1];
                                vec2[2] = cPoint[2] - p4[2];

                                close2node3 = (DotProduct(vec1,vec1)<mindist2);
                                close2node4 = (DotProduct(vec2,vec2)<mindist2);
    
                                if ((node4->myTag.domainID != thisDomain) &&
                                    close2node4) {
                                    continue;
                                }

                                if (close2node3) {
                                    mergenode2 = node3;
                                    splitSeg2 = 0;
                                    adjustCollisionPoint = 1;
                                } else if (close2node4) {
                                    mergenode2 = node4;
                                    splitSeg2 = 0;
                                    adjustCollisionPoint = 1;
                                } else {
                                    splitSeg2 = 1;
                                }


/*
 *                              If we need to add a new node to the second
 *                              segment, do it now.
 */
                                if (splitSeg2) {
                                     real8 pos0[3], pos1[3];

                                     newx = x3 * (1.0-L2) + x4*L2;
                                     newy = y3 * (1.0-L2) + y4*L2;
                                     newz = z3 * (1.0-L2) + z4*L2;
    
                                     newvx = vx3 * (1.0-L2) + vx4*L2;
                                     newvy = vy3 * (1.0-L2) + vy4*L2;
                                     newvz = vz3 * (1.0-L2) + vz4*L2;
    
/*
 *                                   Estimate resulting forces on all segments
 *                                   involved in the split.
 */
                                     arm43 = GetArmID(home, node4, node3);
    
                                     burg1[0] = node3->burgX[arm34];
                                     burg1[1] = node3->burgY[arm34];
                                     burg1[2] = node3->burgZ[arm34];
    
                                     oldfp0[0] = node3->armfx[arm34];
                                     oldfp0[1] = node3->armfy[arm34];
                                     oldfp0[2] = node3->armfz[arm34];
    
                                     oldfp1[0] = node4->armfx[arm43];
                                     oldfp1[1] = node4->armfy[arm43];
                                     oldfp1[2] = node4->armfz[arm43];
    
                                     pos0[0] = x3;   pos1[0] = x4;
                                     pos0[1] = y3;   pos1[1] = y4;
                                     pos0[2] = z3;   pos1[2] = z4;
    
                                     pnew[0] = newx;
                                     pnew[1] = newy;
                                     pnew[2] = newz;
    
                                     newNodeVel[0] = newvx;
                                     newNodeVel[1] = newvy;
                                     newNodeVel[2] = newvz;

                                     nodeVel[0] = vx3;
                                     nodeVel[1] = vy3;
                                     nodeVel[2] = vz3;

                                     FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
                                                 oldfp1, pnew, f0seg1, f1seg1,
                                                 f0seg2, f1seg2);
    
                                     oldTag1 = node3->myTag;
                                     oldTag2 = node4->myTag;

                                     FoldBox(param, &pnew[0], &pnew[1], &pnew[2]);
                                     splitStatus = SplitNode(home,
                                                            OPCLASS_COLLISION,
                                                            node3, pos0,
                                                            pnew, nodeVel,
                                                            newNodeVel, 1,
                                                            &arm34, globalOp,
                                                            &splitNode1,
                                                            &splitNode2, 0
                                                            ,node3->inGrain, node3->onGB
                                                            ,inGrain2, onGB2);
/*
 *                                   If we were unable to split the node
 *                                   go back to looking for more collision
 *                                   candidates.
 */
                                     if (splitStatus == SPLIT_FAILED) {
                                         continue;
                                     }
/*
 *                                   The force estimates above are good enough
 *                                   for the remainder of this timestep, but
 *                                   mark the force and velocity data for some
 *                                   nodes as obsolete so that more accurate
 *                                   forces will be recalculated either at the
 *                                   end of this timestep, or the beginning of
 *                                   the next.
 */
                                     mergenode2 = splitNode2;

                                     MarkNodeForceObsolete(home, splitNode2);
    
                                     for (q = 0; q < splitNode2->numNbrs; q++) {
                                         tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
                                         if (tmpNbr != (Node_t *)NULL) {
                                             tmpNbr->flags |= NODE_RESET_FORCES;
                                         }
                                     }
    
/*
 *                                   Reset nodal forces on nodes involved in the
 *                                   split.
 */
                                     ResetSegForces(home, splitNode1,
                                                    &splitNode2->myTag,
                                                    f0seg1[0], f0seg1[1],
                                                    f0seg1[2], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &splitNode1->myTag,
                                                    f1seg1[0], f1seg1[1],
                                                    f1seg1[2], 1);
    
                                     ResetSegForces(home, splitNode2,
                                                    &node4->myTag,
                                                    f0seg2[0], f0seg2[1],
                                                    f0seg2[2], 1);
    
                                     ResetSegForces(home, node4,
                                                    &splitNode2->myTag,
                                                    f1seg2[0], f1seg2[1],
                                                    f1seg2[2], 1);
 
                                     (void)EvaluateMobility(home, splitNode1);
                                     (void)EvaluateMobility(home, splitNode2);
                                     (void)EvaluateMobility(home, node4);
#if 0
/*
 *                                   When debugging, dump some info on
 *                                   topological changes taking place and
 *                                   the nodes involved
 */
                                     if (gbCollision) {
                                         printf("%d  Split-2(SegCollision): "
                                               "(%d,%d)--(%d,%d) ==> "
                                               "(%d,%d,%d)--(%d,%d,%d)--(%d,%d,%d)\n", home->myDomain,
                                               oldTag1.domainID,
                                               oldTag1.index,
                                               oldTag2.domainID,
                                               oldTag2.index,
                                               splitNode1->myTag.domainID,
                                               splitNode1->myTag.index, splitNode1->inGrain,
                                               splitNode2->myTag.domainID,
                                               splitNode2->myTag.index, splitNode2->inGrain,
                                               node4->myTag.domainID,
                                               node4->myTag.index,node3->inGrain);
                                    }
#endif
                                }  /* if (splitSeg2) */
    
                                if(SurfaceNode(mergenode1) && SurfaceNode(mergenode2))continue;
                                if(splitInGrain >-1 && (mergenode1->onGB>-1 || mergenode2->onGB>-1)){
                                    splitInGrain = -1; splitOnGB = (mergenode1->onGB>-1?mergenode1->onGB:mergenode2->onGB);
                                }
//                              if(splitOnGB > -1 && (mergenode1->numNbrs+mergenode2->numNbrs)>6)continue;
//                              if(mergenode1->onGB != onGB1 || mergenode2->onGB != onGB2)continue;
//                              if(mergenode1->inGrain != inGrain1 || mergenode2->inGrain != inGrain2)continue;
    
/*
 *                             If the initially selected collision point
 *                             was close enough to one of the segment
 *                             endpoints that we decided to use the
 *                             endpoint node for the merge, we need to
 *                             recalculate the final collision point.
 */
                               if (adjustCollisionPoint) {
                                   FindCollisionPoint(home, mergenode1,
                                                      mergenode2, &newPos[0],
                                                      &newPos[1], &newPos[2]);
                               } else {
                                   newPos[0] = cPoint[0];
                                   newPos[1] = cPoint[1];
                                   newPos[2] = cPoint[2];
                               }
                                
                               FoldBox(param, newPos, newPos+1, newPos+2);
#if 0
                               if(gbCollision){
                                    printf("%d ",home->myDomain); PrintTag(mergenode1); PrintTag(mergenode2);
                                    printf("plan to merge at {%g,%g,%g}\n",newPos[0],newPos[1],newPos[2]);
                                    
                               }
#endif
#ifdef _FEM
/*
 *                             If colliding 2 surface nodes, we may have to
 *                             adjust the collision point so it too is on the
 *                             surface.
 */
                               resetSurfaceProperties = 0;

                               if ((mergenode1->constraint == SURFACE_NODE) &&
                                   (mergenode2->constraint == SURFACE_NODE)) {
                                   Node_t *seg1Node2, *seg2Node2;

                                   seg1Node2 = (mergenode1 == node1) ?
                                               node2 : node1;
                                   seg2Node2 = (mergenode2 == node3) ?
                                               node4 : node3;

                                   FEM_AdjustCollisionPoint(mergenode1,
                                                            seg1Node2,
                                                            mergenode2,
                                                            seg2Node2,
                                                            newPos, femSurface,
                                                            femSurfaceNorm);
                                   resetSurfaceProperties = 1;
                               }
#endif

                               newVel[0] = half * (mergenode1->vX +
                                                   mergenode2->vX);
                               newVel[1] = half * (mergenode1->vY +
                                                   mergenode2->vY);
                               newVel[2] = half * (mergenode1->vZ +
                                                   mergenode2->vZ);
    
                               oldTag1 = mergenode1->myTag;
                               oldTag2 = mergenode2->myTag;
                    
                               if(param->polycrystal){
                                   if(splitOnGB>-1){
                                       gb = home->GBKeys[splitOnGB];
                                       PointPlaneDistance(newPos, gb->n, gb->p, dis, newPos);
                                       if(fabs(dis)>param->rTol){
                                           continue;
                                       }
                                   }
                                   
                                   canMerge = BeforeMergeNodes(home, mergenode1,mergenode2,newPos,&dis,type);
                                   if(!canMerge){
                                         mergenode1->flags |= NO_COLLISIONS;
                                         mergenode2->flags |= NO_COLLISIONS;
                                         continue;
                                   }
                                   NodeLocation(home, newPos[0],newPos[1],newPos[2],newInGrain,newOnGB);
                                   if(splitOnGB>-1){if(splitOnGB != newOnGB)canMerge=false;}
                                   else{if(splitInGrain!=newInGrain)canMerge=false;}                            
                                   
                                   if(!canMerge){
                                        if(splitOnGB>-1){
                                            gb = home->GBKeys[splitOnGB];
                                            PointPlaneDistance(newPos, gb->n, gb->p, dis, newPos);
                                            if(fabs(dis)>param->rTol || 
                                               (mergenode1->inGrain>-1 && ParallelGB(mergenode1,gb->n)) ||
                                               (mergenode2->inGrain>-1 && ParallelGB(mergenode2,gb->n)) ){
                                                mergenode1->flags |= NO_COLLISIONS;
                                                mergenode2->flags |= NO_COLLISIONS;
                                                continue;
                                            }
                                            
                                            NewPosToGB(home, mergenode1, -1, splitOnGB, newPos, vec1);
                                            NewPosToGB(home, mergenode2, -1, splitOnGB, newPos, vec2);
                                            vec1[0] -= vec2[0];
                                            vec1[1] -= vec2[1];
                                            vec1[2] -= vec2[2];
                                            ZImage(param, vec1, vec1+1, vec1+2);
                                        
                                            if(Normal(vec1) > 2.0*param->rTol){
                                                PrintNode(mergenode1);
                                                PrintNode(mergenode2);
                                                printf("Warning: can not merge the two nodes in PredictiveCollision. %f\n",
                                                       Normal(vec1));
                                                continue;
                                            }
                                        
                                            newPos[0] = vec2[0] + 0.5*vec1[0];
                                            newPos[1] = vec2[1] + 0.5*vec1[1];
                                            newPos[2] = vec2[2] + 0.5*vec1[2];
                                        }else{
                                            if(newInGrain != splitInGrain){
                                        
                                                AdjustSplitedNodePosition(home, mergenode1, newPos, vec1);
                                                AdjustSplitedNodePosition(home, mergenode2, newPos, vec2);
                                        
                                        if(mergenode1->onGB>-1 && mergenode1->constraint==7){
                                            printf("Task %d %s-2: attempt move GB node ",home->myDomain,__func__);
                                            PrintTag(home, mergenode1); 
                                            printf(" form {%g,%g,%g} to {%g,%g,%g}\n", newPos[0],newPos[1],newPos[2],vec1[0],vec1[1],vec1[2]);
                                        }
                                        if(mergenode2->onGB>-2 && mergenode2->constraint==7){
                                            printf("Task %d %s-1: attempt move GB node ",home->myDomain,__func__);
                                            PrintTag(home, mergenode2); 
                                            printf(" form {%g,%g,%g} to {%g,%g,%g}\n", newPos[0],newPos[1],newPos[2],vec2[0],vec2[1],vec2[2]);
                                        }
                                                vec1[0] -= vec2[0];
                                                vec1[1] -= vec2[1];
                                                vec1[2] -= vec2[2];
                                                ZImage(param, vec1, vec1+1, vec1+2);
                                        
                                                newPos[0] = vec2[0] + 0.5*vec1[0];
                                                newPos[1] = vec2[1] + 0.5*vec1[1];
                                                newPos[2] = vec2[2] + 0.5*vec1[2];
                                            }
                                        }
                                    }
                                }
                                FoldBox(param, newPos, newPos+1, newPos+2);
//                              printf("Task %d %d %s-1: (%d,%d,%d)-(%d,%d,%d) plan to merge at (%g,%g,%g)\n", home->myDomain,home->cycle, __func__,
//                                     mergenode1->myTag.domainID, mergenode1->myTag.index, mergenode1->inGrain,
//                                     mergenode2->myTag.domainID, mergenode2->myTag.index, mergenode2->inGrain,
//                                     newPos[0], newPos[1], newPos[2]);
                               MergeNode(home, OPCLASS_COLLISION, mergenode1,
                                         mergenode2, newPos, &targetNode,
                                         &mergeStatus, globalOp,
                                         splitInGrain, splitOnGB);
//                              if((mergeStatus & MERGE_SUCCESS) == 0){
//                                  printf("Task %d %s-1: merge failed, status %d\n", home->myDomain, __func__,mergeStatus);
//                              }else{
//                                  printf("Task %d %s-1: merge succeed, ", home->myDomain, __func__);
//                                  PrintTag(targetNode);
//                                  printf("\n");
//                              }
    
/*
 *                             If the merge did not succeed, go back and
 *                             continue looking for collision candidates.
 */
                               if ((mergeStatus & MERGE_SUCCESS) == 0) {
                                   continue;
                               }
#ifdef _FEM
/*
 *                             Need to explicitly reset surface properties
 *                             after colliding 2 surface nodes.
 */
                               if (resetSurfaceProperties) {
                                   targetNode->fem_Surface[0] = femSurface[0];
                                   targetNode->fem_Surface[1] = femSurface[1];
                                   targetNode->fem_Surface_Norm[0] =
                                           femSurfaceNorm[0];
                                   targetNode->fem_Surface_Norm[1] =
                                           femSurfaceNorm[1];
                                   targetNode->fem_Surface_Norm[2] =
                                           femSurfaceNorm[2];
                               }
#endif

/*
 *                             When debugging, dump some info on topological
 *                             changes taking place and the nodes involved
 */
#if 0
                               printf("  Merge(SegCollision): "
                                      "(%d,%d) and (%d,%d) at (%d,%d,%d,%d)\n",
                                      oldTag1.domainID, oldTag1.index,
                                      oldTag2.domainID, oldTag2.index,
                                      targetNode->myTag.domainID,
                                      targetNode->myTag.index, targetNode->inGrain, targetNode->onGB);
#endif

/*
 *                             If the target node exists after the merge,
 *                             reset its velocity, and reset the topological
 *                             change exemptions for the target node; use
 *                             NodeTopologyExemptions() to get the basic
 *                             exemptions, then exempt all the node's arms
 *                             from additional segment collisions this cycle.
 */
                               if (targetNode != (Node_t *)NULL) {
/*
 *                                 If we are enforcing glide planes but
 *                                 allowing some fuzziness in the planes, we
 *                                 also need to recalculate the glide 
 *                                 planes for the segments attched to the
 *                                 collision node.
 */
//                                  if(IsOnGB(targetNode))printf("Task %d: PredictiveCollsion merged node(%d,%d) arms %d\n", home->myDomain
//                                      ,targetNode->myTag.domainID,targetNode->myTag.index, targetNode->numNbrs);
                                   if(!IsOnGB(targetNode) && home->param->polycrystal){
                                        NodeLocation(home, targetNode->x, targetNode->y,
                                                     targetNode->z, inGrain1, onGB1);
                                        if(inGrain1 != targetNode->inGrain){
                                            PrintNode(targetNode);
                                            Fatal("Task %d There is something wrong in %s at %d %d - 1",home->myDomain,
                                                   __FILE__, __LINE__, mergeStatus);
                                        } 
                                   }

                                   if (param->enforceGlidePlanes &&
                                       param->allowFuzzyGlidePlanes) {
                                       int n;
                                       for (n=0; n<targetNode->numNbrs; n++) {
                                           tmpNbr = GetNodeFromTag(home,
                                                   targetNode->nbrTag[n]);
                                           RecalcSegGlidePlane(home,
                                                               targetNode,
                                                               tmpNbr, 1);
                                       }
                                   }

/*
 *                                 Estimate velocity so mobility function
 *                                 has a reasonable starting point
 */
                                   targetNode->vX = newVel[0];
                                   targetNode->vY = newVel[1];
                                   targetNode->vZ = newVel[2];

                                   (void)EvaluateMobility(home, targetNode);

                                   targetNode->flags |= NO_COLLISIONS;
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                                   targetNode->multiNodeLife = 0;
#endif
                               }
    
                               didCollision = 1;
                               localCollisionCnt++;

                            }  /* If conditions for collision are met */
                        }  /* Loop over node3 arms */
                    }  /* Loop over node1 arms */
                }  /* while(nextIndex >= 0) */
            }  /* Loop over neighboring cell2s */
           }
          }
        }  /* for (i = 0 ...) */
        DebugGB(home, 402);

/*
 *      Now we have to loop for collisions on hinge joints (i.e zipping)
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            for (j = 0; j < node1->numNbrs; j++) {

                if (node1->myTag.domainID != node1->nbrTag[j].domainID)
                    continue;

                for (k = j + 1; k < node1->numNbrs; k++) {
                    real8 dx, dy, dz;

                    if (node1->myTag.domainID != node1->nbrTag[k].domainID)
                        continue;

                    node3 = GetNodeFromTag(home, node1->nbrTag[j]);
                    node4 = GetNodeFromTag(home, node1->nbrTag[k]);

                    splitInGrain = -1;
                    splitOnGB = -1;
                    gbCollision = 0;
                    if(home->param->polycrystal){
                        GetArmGrain(home, node1, node3, inGrain1, onGB1);
                        GetArmGrain(home, node1, node4, inGrain2, onGB2);
                                
                        if(onGB1 > -1 || onGB2 > -1){
                            if(onGB1 == onGB2){
                                splitOnGB = onGB1;
                            }else if(onGB1 > -1 && onGB2 > -1){
                                continue;
                            }else if(onGB1 > -1){
                                splitOnGB = onGB1;
                                gbCollision = 2;
                            }else if(onGB2 > -1){
                                splitOnGB = onGB2;
                                gbCollision = 2;
                            }else{
                                continue;
                            }
                        }else{
                            if(inGrain1 == inGrain2){
                                splitInGrain = inGrain1;
                            }else{
                                if(inGrain1 == -1 || inGrain2 == -1){
                                    Fatal("in %s at %d", __FILE__, __LINE__);
                                }
                                real8 refPos[3]={node1->x, node1->y, node1->z};
                                if(GetPublicGB(home,home->grainKeys[inGrain1],
                                   home->grainKeys[inGrain2], &gb,refPos)){
                                   splitOnGB = gb->GBID;
                                    gbCollision =1;
                                }else{
                                    continue;
                                }
                            }
                        }
                    }

                    x1 = node1->x; y1 = node1->y; z1 = node1->z;
                    x3 = node3->x; y3 = node3->y; z3 = node3->z;
                    x4 = node4->x; y4 = node4->y; z4 = node4->z;

                    vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
                    vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
                    vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;

                    PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                    PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

                    x2  = x1;  y2  = y1;  z2  = z1;
                    vx2 = vx1; vy2 = vy1; vz2 = vz1;

/*
 *                  It is possible to have a zero-length segment
 *                  (created by a previous collision).  If we find
 *                  such a segment, do not try to use it in any
 *                  subsequent collisions.
 */
                    dx = x1 - x3;
                    dy = y1 - y3;
                    dz = z1 - z3;

                    if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                        continue;
                    }

                    dx = x1 - x4;
                    dy = y1 - y4;
                    dz = z1 - z4;

                    if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                        continue;
                    }

/*
 *                  Find the minimum distance between the the node1/node4
 *                  segment and the point at node3 to determine if they
 *                  should be collided.
 */
                    GetMinDist(x1, y1, z1, vx1, vy1, vz1,
                               x4, y4, z4, vx4, vy4, vz4,
                               x3, y3, z3, vx3, vy3, vz3,
                               x3, y3, z3, vx3, vy3, vz3,
                               &dist2, &ddist2dt, &L1, &L2);

                    if (dist2 < mindist2) {

                        FindCollisionPointAndTime(home, p1, p4, p3, p3,
                                              v1, v4, v3, v3, cPoint, &cTime,
                                              &dist2, &L1, &L2);

                        collisionConditionIsMet = ((cTime > 0.0) &&
                                                   (cTime < 10.0));
                    } else {
                        collisionConditionIsMet = 0;
                    }

                    if (collisionConditionIsMet) {

/*
 *                      A collision should occur, so calculate the ideal
 *                      collision point between the node1/node4 segment and
 *                      point node3, using the 'closest' point identified above
 */

/*
 *                      Node3 is used as one of the collision points, but see
 *                      if the collision point is close to one of the existing
 *                      nodes on the node1/node4 segment
 */
                        mergenode1 = node3;

                        vec1[0] = cPoint[0] - x1;
                        vec1[1] = cPoint[1] - y1;
                        vec1[2] = cPoint[2] - z1;

                        vec2[0] = cPoint[0] - x4;
                        vec2[1] = cPoint[1] - y4;
                        vec2[2] = cPoint[2] - z4;

                        close2node1 = (DotProduct(vec1,vec1)<mindist2);
                        close2node4 = (DotProduct(vec2,vec2)<mindist2);

/*
 *                      If the collision point is close to one of the
 *                      endpoints of the node1/node4 segment, but the
 *                      node in question is not owned by the current
 *                      domain, skip the collision this time around.
 */
                        if ((close2node1&&(node1->myTag.domainID!=thisDomain))||
                            (close2node4&&(node4->myTag.domainID!=thisDomain))){
                            continue;
                        }

                        if (close2node1) {
                             mergenode2 = node1;
                             adjustCollisionPoint = 1;
                        } else if (close2node4) {
                             mergenode2 = node4;
                             adjustCollisionPoint = 1;
                        } else { 
                             real8 pos0[3], pos1[3];

                             adjustCollisionPoint = 0;
/*
 *                           Collision point is not close enough to one of
 *                           the segment endpoints, so bisect the segment
 *                           and add a new node.
 */
                             newPos[0] = x1*(1.0-L1) + x4*L1;
                             newPos[1] = y1*(1.0-L1) + y4*L1;
                             newPos[2] = z1*(1.0-L1) + z4*L1;

                             newVel[0] = vx1*(1.0-L1) + vx4*L1;
                             newVel[1] = vy1*(1.0-L1) + vy4*L1;
                             newVel[2] = vz1*(1.0-L1) + vz4*L1;

                             nodeVel[0] = node3->vX;
                             nodeVel[1] = node3->vY;
                             nodeVel[2] = node3->vZ;

/*
 *                           Estimate resulting forces on all segments
 *                           involved in the split.
 */
                             armAB = GetArmID(home, node1, node4);
                             armBA = GetArmID(home, node4, node1);

                             burg1[0] = node1->burgX[armAB];
                             burg1[1] = node1->burgY[armAB];
                             burg1[2] = node1->burgZ[armAB];

                             oldfp0[0] = node1->armfx[armAB];
                             oldfp0[1] = node1->armfy[armAB];
                             oldfp0[2] = node1->armfz[armAB];

                             oldfp1[0] = node4->armfx[armBA];
                             oldfp1[1] = node4->armfy[armBA];
                             oldfp1[2] = node4->armfz[armBA];

                             pos0[0] = x1;   pos1[0] = x4;
                             pos0[1] = y1;   pos1[1] = y4;
                             pos0[2] = z1;   pos1[2] = z4;

                             FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
                                         oldfp1, newPos, f0seg1, f1seg1,
                                         f0seg2, f1seg2);

                             FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);

/*
 *                           Attempt to split the segment.  If the split
 *                           fails, we can't continue with this collision.
 */
                             oldTag1 = node1->myTag;
                             oldTag2 = node4->myTag;

                             splitStatus = SplitNode(home, OPCLASS_COLLISION,
                                                     node1, pos0, newPos,
                                                     nodeVel, newVel, 1,
                                                     &k, globalOp,
                                                     &splitNode1,
                                                     &splitNode2, 0
                                                     ,node1->inGrain, node1->onGB
                                                     ,inGrain2, onGB2
                                                                     );

                             if (splitStatus != SPLIT_SUCCESS) {
                                 continue;
                             }

/*
 *                           The force estimates above are good enough for
 *                           the remainder of this timestep, but mark the
 *                           force and velocity data for some nodes as obsolete
 *                           so that more accurate forces will be recalculated
 *                           either at the end of this timestep, or the
 *                           beginning of the next.
 */
                             mergenode2 = splitNode2;

                             MarkNodeForceObsolete(home, splitNode1);
                             MarkNodeForceObsolete(home, splitNode2);
                             MarkNodeForceObsolete(home, node4);

/*
 *                           Reset nodal forces for nodes involved in the
 *                           split.
 */
                             ResetSegForces(home,splitNode1,&splitNode2->myTag,
                                            f0seg1[0], f0seg1[1], f0seg1[2], 1);

                             ResetSegForces(home,splitNode2,&splitNode1->myTag,
                                            f1seg1[0], f1seg1[1], f1seg1[2], 1);

                             ResetSegForces(home, splitNode2, &node4->myTag,
                                            f0seg2[0], f0seg2[1], f0seg2[2], 1);

                             ResetSegForces(home, node3, &splitNode2->myTag,
                                            f1seg2[0], f1seg2[1], f1seg2[2], 1);

                             (void)EvaluateMobility(home, splitNode1);
                             (void)EvaluateMobility(home, splitNode2);
                             (void)EvaluateMobility(home, node4);
#if 0
/*
 *                           When debugging, dump some info on topological
 *                           changes taking place and the nodes involved
 */
                             if (gbCollision) {
                                 printf("%d,%d Split-1(Hinge): "
                                        "(%d,%d)--(%d,%d) ==> "
                                        "(%d,%d,%d)--(%d,%d,%d)--(%d,%d,%d)\n", home->myDomain,gbCollision,
                                        oldTag1.domainID, oldTag1.index,
                                        oldTag2.domainID, oldTag2.index,
                                        splitNode1->myTag.domainID,
                                        splitNode1->myTag.index, splitNode1->inGrain,
                                        splitNode2->myTag.domainID, 
                                        splitNode2->myTag.index,splitNode2->inGrain,
                                        node4->myTag.domainID,
                                        node4->myTag.index,node4->inGrain);
                             }
#endif
                        }

                        if(splitInGrain >-1 && (mergenode1->onGB>-1 || mergenode2->onGB>-1)){
                            splitInGrain = -1; splitOnGB = (mergenode1->onGB>-1?mergenode1->onGB:mergenode2->onGB);
                        }
//                      if(splitOnGB > -1 && (mergenode1->numNbrs+mergenode2->numNbrs)>6)continue;
/*
 *                      If the initially selected collision point
 *                      was close enough to one of the segment
 *                      endpoints that we decided to use the
 *                      endpoint node for the merge, we need to
 *                      recalculate the final collision point.
 */
                        if (adjustCollisionPoint) {
                            FindCollisionPoint(home, mergenode1,
                                               mergenode2, &newPos[0],
                                               &newPos[1], &newPos[2]);
                        }

                        FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);
#if 0
                        if(gbCollision){
                             printf("%d ",home->myDomain); PrintTag(mergenode1); PrintTag(mergenode2);
                             printf("plan to merge at {%g,%g,%15.7e}\n",newPos[0],newPos[1],newPos[2]);
                        }
#endif
/*
 *                      The merge is going to happen, so there's a few
 *                      more nodes we'll need to mark for force/velocity
 *                      re-evaluation.
 */
                        MarkNodeForceObsolete(home, node1);
                        MarkNodeForceObsolete(home, node3);
                        MarkNodeForceObsolete(home, node4);

                        for (q = 0; q < node1->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node1->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

                        for (q = 0; q < node3->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node3->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

                        for (q = 0; q < node4->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node4->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

/*
 *                      Try the collision
 */
                        oldTag1 = mergenode1->myTag;
                        oldTag2 = mergenode2->myTag;

                        if(param->polycrystal){
                            if(splitOnGB>-1){
                                gb = home->GBKeys[splitOnGB];
                                PointPlaneDistance(newPos, gb->n, gb->p, dis, newPos);
                                if(fabs(dis)>2.0*param->rTol){
                                    continue;
                                }
                            }

                            canMerge = BeforeMergeNodes(home, mergenode1,mergenode2,newPos,&dis,type);
                            if(!canMerge){
                                mergenode1->flags |= NO_COLLISIONS;
                                mergenode2->flags |= NO_COLLISIONS;
                                continue;
                            }
                    
                            NodeLocation(home, newPos[0],newPos[1],newPos[2],newInGrain,newOnGB);
                            if(splitOnGB>-1){if(splitOnGB != newOnGB)canMerge=false;}
                            else{if(splitInGrain!=newInGrain)canMerge=false;}                            

                            if(!canMerge){
                                if(splitOnGB>-1){
                                    gb = home->GBKeys[splitOnGB];
                                    PointPlaneDistance(newPos, gb->n, gb->p, dis, newPos);
                                    if(fabs(dis)>param->rTol || 
                                       (mergenode1->inGrain>-1 && ParallelGB(mergenode1,gb->n)) ||
                                       (mergenode2->inGrain>-1 && ParallelGB(mergenode2,gb->n)) ){
                                        mergenode1->flags |= NO_COLLISIONS;
                                        mergenode2->flags |= NO_COLLISIONS;
                                        continue;
                                    }
                                    
                                    NewPosToGB(home, mergenode1, -1, splitOnGB, newPos, vec1);
                                    NewPosToGB(home, mergenode2, -1, splitOnGB, newPos, vec2);
                                    vec1[0] -= vec2[0];
                                    vec1[1] -= vec2[1];
                                    vec1[2] -= vec2[2];
                                    ZImage(param, vec1, vec1+1, vec1+2);

                                    if(Normal(vec1) > param->rTol){
                                        PrintNode(mergenode1);
                                        PrintNode(mergenode2);
                                        printf("Warning: can not merge the two nodes in PredictiveCollision. %f\n",
                                               Normal(vec1));
                                        continue;
                                    }

                                    newPos[0] = vec2[0] + 0.5*vec1[0];
                                    newPos[1] = vec2[1] + 0.5*vec1[1];
                                    newPos[2] = vec2[2] + 0.5*vec1[2];
                                }else{
                                    if(newInGrain != splitInGrain){

                                        AdjustSplitedNodePosition(home, mergenode1, newPos, vec1);
                                        AdjustSplitedNodePosition(home, mergenode2, newPos, vec2);
                                        if(mergenode1->onGB>-1 && mergenode1->constraint==7){
                                            printf("Task %d %s-1: attempt move GB node ",home->myDomain,__func__);
                                            PrintTag(home, mergenode1); 
                                            printf(" form {%g,%g,%g} to {%g,%g,%g}\n", newPos[0],newPos[1],newPos[2],vec1[0],vec1[1],vec1[2]);
                                        }
                                        if(mergenode2->onGB>-1 && mergenode2->constraint==7){
                                            printf("Task %d %s-1: attempt move GB node ",home->myDomain,__func__);
                                            PrintTag(home, mergenode2); 
                                            printf(" form {%g,%g,%g} to {%g,%g,%g}\n", newPos[0],newPos[1],newPos[2],vec2[0],vec2[1],vec2[2]);
                                        }

                                        vec1[0] -= vec2[0];
                                        vec1[1] -= vec2[1];
                                        vec1[2] -= vec2[2];
                                        ZImage(param, vec1, vec1+1, vec1+2);

                                        newPos[0] = vec2[0] + 0.5*vec1[0];
                                        newPos[1] = vec2[1] + 0.5*vec1[1];
                                        newPos[2] = vec2[2] + 0.5*vec1[2];
                                    }
                                }
                            }
                        }

                        FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);
//                              printf("Task %d %d %s-2: (%d,%d,%d)-(%d,%d,%d) plan to merge at (%g,%g,%g)\n", home->myDomain,home->cycle, __func__,
//                                     mergenode1->myTag.domainID, mergenode1->myTag.index, mergenode1->inGrain,
//                                     mergenode2->myTag.domainID, mergenode2->myTag.index, mergenode2->inGrain,
//                                     newPos[0], newPos[1], newPos[2]);
                        MergeNode(home, OPCLASS_COLLISION, mergenode1,
                                  mergenode2, newPos, &targetNode,
                                  &mergeStatus, globalOp,
                                  splitInGrain, splitOnGB);

//                              if((mergeStatus & MERGE_SUCCESS) == 0){
//                                  printf("Task %d %s-2: merge failed, status %d\n", home->myDomain, __func__,mergeStatus);
//                              }else{
//                                  printf("Task %d %s-2: merge succeed, ", home->myDomain, __func__);
//                                  PrintTag(targetNode);
//                                  printf("\n");
//                              }
                        localCollisionCnt++;

/*
 *                      If the target node exists after the merge, reevaulate
 *                      the node velocity, mark the forces as obsolete,
 *                      and reset the topological change exemptions for
 *                      the target node and exempt all the node's arms
 *                      from subsequent segment collisions this cycle.
 */
                        if (targetNode != (Node_t *)NULL) {
//                          if(IsOnGB(targetNode))printf("Task %d: PredictiveCollsion hinge merged node(%d,%d) arms %d\n", home->myDomain
//                              ,targetNode->myTag.domainID,targetNode->myTag.index, targetNode->numNbrs);
                            if(!IsOnGB(targetNode) && home->param->polycrystal){
                                 NodeLocation(home, targetNode->x, targetNode->y,
                                              targetNode->z, inGrain1, onGB1);
                                 if(inGrain1 != targetNode->inGrain){
                                     PrintNode(targetNode);
                                     Fatal("There is something wrong in %s at %d - 2",
                                            __FILE__, __LINE__);
                                 } 
                            }
/*
 *                          If we are enforcing glide planes but
 *                          allowing some fuzziness in the planes, we
 *                          also need to recalculate the glide 
 *                          planes for the segments attched to the
 *                          collision node.
 */
                            if (param->enforceGlidePlanes &&
                                param->allowFuzzyGlidePlanes) {
                                int n;
                                for (n = 0; n < targetNode->numNbrs; n++) {
                                    tmpNbr = GetNodeFromTag(home,
                                            targetNode->nbrTag[n]);
                                    RecalcSegGlidePlane(home, targetNode,
                                                        tmpNbr, 1);
                                }
                            }

                            (void)EvaluateMobility(home, targetNode);

                            targetNode->flags |= NODE_RESET_FORCES;
                            targetNode->flags |= NO_COLLISIONS;
#if 0
                            if (gbCollision) {
                                printf("Hinge(%d,%d)  Merge(HingeCollision): "
                                       "(%d,%d) and (%d,%d) at (%d,%d)\n",home->myDomain,gbCollision,
                                       oldTag1.domainID, oldTag1.index,
                                       oldTag2.domainID, oldTag2.index,
                                       targetNode->myTag.domainID,
                                       targetNode->myTag.index);
                            }
#endif
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                            targetNode->multiNodeLife = 0;
#endif
                        }

                    }  /* conditions for zip were met */
                }  /* for (k = 0...) */
            }  /* for (j = 0...) */
        }  /* for (i = 0...) */
        DebugGB(home, 403);


#ifdef DEBUG_LOG_COLLISIONS
#ifdef PARALLEL
        MPI_Reduce(&localCollisionCnt, &globalCollisionCnt, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
        globalCollisionCnt = localCollisionCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Collision count = %d\n", globalCollisionCnt);
        }
#endif

        TimerStop(home, COLLISION_HANDLING);

        return;
}
