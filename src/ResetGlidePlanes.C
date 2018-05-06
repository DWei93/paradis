#include "Home.h"
#include "Comm.h"
/*---------------------------------------------------------------------------
 *
 *      Function:     ResetGlidePlanes
 *      Description:  If needed, this function will recalculate the glide
 *                    plane info for all locally owned segments, then 
 *                    distribute the plane information to the remote
 *                    domains which have these segments as ghosts.
 *
 *-------------------------------------------------------------------------*/
void ResetGlidePlanes(Home_t *home)
{
        int     i, j;
        int     nbrSegID, numNodes, numNbrs;
        real8   lDir[3], burg[3], newPlane[3];
        Node_t  *node, *nbrNode;
        Param_t *param;

        param = home->param;
        numNodes = home->newNodeKeyPtr;

/*
 *      Currently, the only situation in which we need to go through this
 *      process is when glide planes are being enforced but allowed some
 *      'fuzziness'.
 */
        if (!(param->enforceGlidePlanes && param->allowFuzzyGlidePlanes)) {
            return;
        }
 
/*
 *      Loop through all the nodes looking for segments owned
 *      by the local domain.
 */
        for (i = 0; i < numNodes; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    continue;
                }

/*
 *              Only check the segment from one end so there's no
 *              duplication of work done.
 */
                if (OrderNodes(node, nbrNode) < 0) {
                    continue;
                }

                burg[0] = node->burgX[j];
                burg[1] = node->burgY[j];
                burg[2] = node->burgZ[j];

                lDir[0] = nbrNode->x - node->x;
                lDir[1] = nbrNode->y - node->y;
                lDir[2] = nbrNode->z - node->z;

                ZImage(param, &lDir[0], &lDir[1], &lDir[2]);
                NormalizeVec(lDir);

/*
 *              Calculate the new glide plane.  If it is undefined,
 *              the segment is screw and should simply maintain the
 *              current glide plane.
 */
                FindPreciseGlidePlane(home, burg, lDir, newPlane, node, nbrNode);

                if (fabs(DotProduct(newPlane, newPlane)) < 1.0e-03) {
                    continue;
                }

/*
 *              Segment is not screw, so reset the segment' plane
 *              to the newly calculated one.
 */
                node->nx[j] = newPlane[0];
                node->ny[j] = newPlane[1];
                node->nz[j] = newPlane[2];

                nbrSegID = GetArmID(home, nbrNode, node);

                nbrNode->nx[nbrSegID] = newPlane[0];
                nbrNode->ny[nbrSegID] = newPlane[1];
                nbrNode->nz[nbrSegID] = newPlane[2];
            }
        }

/*
 *      Now send the glide plane info to all the domains that use it in the
 *      ghost data
 */
        CommSendGhostPlanes(home);
                     
        return;
}
