/*****************************************************************************
 *
 *      Module:       CreateConfig.c
 *      Description:  Contains the majority of the functions to generate
 *                    nodal data for various types of initial dislocation
 *                    structures.
 *
 *      Includes functions:
 *
 *              CreateEdges()
 *              CreatePrismaticLoop()
 *              CreateScrewConfig()
 *              CreateFiniteMixedConfig()
 *              CreateFCCConfig()
 *              CreateFCCIrradConfig()
 *              CreateFCCPerfectLoop()  >>> Not yet fully implemented <<<
 *
 *****************************************************************************/
#include "Home.h"
#include "InData.h"
#include "Tag.h"
#include "Util.h"
#include "ParadisGen.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

#include <vector>

void ExtPKForce(real8 str[3][3],
                real8 bx, real8 by, real8 bz,
                real8 x1, real8 y1, real8 z1,
                real8 x2, real8 y2, real8 z2,
                real8 f1[3], real8 f2[3])
{
    real8 strb[3], xi[3], ft[3];
    int j;
	
	
    xi[0] = x2-x1;
    xi[1] = y2-y1;
    xi[2] = z2-z1;

    for (j = 0; j < 3; j++) {
        strb[j] = str[j][0]*bx + str[j][1]*by + str[j][2]*bz;
    }

    cross(strb, xi, ft);
	
	
    for (j = 0; j < 3; j++) {
        f1[j] = ft[j]*0.5;
        f2[j] = ft[j]*0.5;
    }
}


static void IncDislocationDensity(InData_t *inData, real8 *totDislocLen)
{
        int      i, armID, nbrIndex;
        real8    xLen, yLen, zLen;
        Param_t  *param;
        Node_t   *node, *nbrNode;

        param = inData->param;

        for (i = 0; i < inData->nodeCount; i++) {

            node = &inData->node[i];

            for (armID = 0; armID < node->numNbrs; armID++) {

                nbrIndex = node->nbrTag[armID].index;
                if (nbrIndex < i) continue;

                nbrNode = &inData->node[nbrIndex];

                xLen = nbrNode->x - node->x;   
                yLen = nbrNode->y - node->y;   
                zLen = nbrNode->z - node->z;   

                ZImage(param, &xLen, &yLen, &zLen);

                *totDislocLen += sqrt(xLen*xLen + yLen*yLen + zLen*zLen);
            }
        }

        return;
}

static void IncDislocationDensity(InData_t *inData, int startID, real8 *totDislocLen)
{
        int      i, armID, nbrIndex;
        real8    xLen, yLen, zLen;
        Param_t  *param;
        Node_t   *node, *nbrNode;

        param = inData->param;

        for (i = startID; i < inData->nodeCount; i++) {

            node = &inData->node[i];

            for (armID = 0; armID < node->numNbrs; armID++) {

                nbrIndex = node->nbrTag[armID].index;
                if (nbrIndex < i) continue;

                nbrNode = &inData->node[nbrIndex];

                xLen = nbrNode->x - node->x;   
                yLen = nbrNode->y - node->y;   
                zLen = nbrNode->z - node->z;   

                ZImage(param, &xLen, &yLen, &zLen);

                *totDislocLen += sqrt(xLen*xLen + yLen*yLen + zLen*zLen);
            }
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     CreatePrismaticLoop
 *      Description:  Generate nodal data for one or more prismatic loops.
 *                    Loops will be either interstitial or vacancy loops;
 *                    default behavior is to generate interstitials.
 *
 *                    Note: [1 1 1] type loops will be hexagonal loops
 *                    
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          loopType    Type of loops to create (i.e. 111, 100)
 *                        0 == mixture (default)
 *                        1 == [1 1 1] all types
 *                        2 == [1 0 0] [0 1 0] [0 0 1] types
 *          useInterstitial If set 1, all loops will be interstitial loops.
 *                      Otherwise loops will be vacancy loops.
 *          numLoops    Number of prismatic loops to create
 *          radius      radius (in units of b) of the loop?
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreatePrismaticLoop(Home_t *home, InData_t *inData, int cubeLength,
                         int loopType, int useInterstitial, int numLoops,
                         real8 radius, int seed, real8 *totDislocLen,
                         int dislocType)
{
        int     id, loopIndex, burgIndex, dIndex, nextNode, newNodeIndex;
        int     minBurgIndex, maxBurgIndex, lastBlock, numSegs;
        int     startRemeshIndex = 0;
        int     nbr1Index;
        real8   cubeSize;
        real8   x, y, z, ux, uy, uz;
        real8   rx, ry, rz;
        real8   loopCtrX,loopCtrY,loopCtrZ;
        real8   invSqrt3, invSqrt6, twoInvSqrt6;
        real8   burg[7][3];
        real8   glidePlane[7][3], xp[7][3], zp[7][3];
        real8   xpp[3], zpp[3], ypp[3], rotMatrix[3][3]; 
        real8   vec1[3], vec2[3], tmpVec[3];
        real8   tr[4][6][3];
        Param_t *param;
        Node_t  *node, *nbr1Node;

        param = inData->param;
        cubeSize = (real8)cubeLength;

/*
 *      All loops are now hexagonal and of [1 1 1] type burgers vectors for now
 */
        numSegs = 6;
        minBurgIndex = 0;
        maxBurgIndex = 3;
      
/*
 *      Define the sets of burger's vectors that may be used.  For each
 *      burgers vector we'll choose a glide plane normal vector which is
 *      perpendicular to the loop and along the burgers vector direction.
 */
        invSqrt3 = 1.0 / sqrt(3.0);
        invSqrt6 = 1.0 / sqrt(6.0);
        twoInvSqrt6 = 2.0 * invSqrt6;

        /*  [1 1 1]  */
        burg[0][0] =  invSqrt3;
        burg[0][1] =  invSqrt3;
        burg[0][2] =  invSqrt3;

        /*  [-1 1 1]  */
        burg[1][0] = -invSqrt3;
        burg[1][1] =  invSqrt3;
        burg[1][2] =  invSqrt3;

        /*  [1 -1 1]  */
        burg[2][0] =  invSqrt3;
        burg[2][1] = -invSqrt3;
        burg[2][2] =  invSqrt3;

        /*  [1 1 -1]  */
        burg[3][0] =  invSqrt3;
        burg[3][1] =  invSqrt3;
        burg[3][2] = -invSqrt3;

/*
 *      Define the translation vectors as 4 X 6 X 3 array 
 *
 *          tr[ ][ ][ ] - translation vector
 *             |  |  |
 *             |  |  x,y,z
 *             |  |
 *             |  6 directions
 *             |
 *             4 [111] type burgers vectors
 */
        if (useInterstitial == 1) {
/*
 *          Interstitial loops for [1 1 1] burgers vector
 */
            tr[0][0][0]=  invSqrt6;
            tr[0][0][1]=  invSqrt6;
            tr[0][0][2]= -twoInvSqrt6;

            /* [-1  2 -1] */
            tr[0][1][0]= -invSqrt6;
            tr[0][1][1]=  twoInvSqrt6;
            tr[0][1][2]= -invSqrt6;

            /* [-2  1  1] */
            tr[0][2][0]= -twoInvSqrt6;
            tr[0][2][1]=  invSqrt6;
            tr[0][2][2]=  invSqrt6;

            /* [-1 -1  2] */
            tr[0][3][0]= -invSqrt6;
            tr[0][3][1]= -invSqrt6;
            tr[0][3][2]=  twoInvSqrt6;

            /* [ 1 -2  1] */
            tr[0][4][0]=  invSqrt6;
            tr[0][4][1]= -twoInvSqrt6;
            tr[0][4][2]=  invSqrt6;

            /* [ 2 -1 -1] */
            tr[0][5][0]=  twoInvSqrt6;
            tr[0][5][1]= -invSqrt6;
            tr[0][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [-1 1 1] burgers vector
 */
            /* [-1  1 -2] */
            tr[1][0][0]= -invSqrt6;
            tr[1][0][1]=  invSqrt6;
            tr[1][0][2]= -twoInvSqrt6;

            /* [-2 -1 -1] */
            tr[1][1][0]= -twoInvSqrt6;
            tr[1][1][1]= -invSqrt6;
            tr[1][1][2]= -invSqrt6;

            /* [-1 -2  1] */
            tr[1][2][0]= -invSqrt6;
            tr[1][2][1]= -twoInvSqrt6;
            tr[1][2][2]=  invSqrt6;

            /* [ 1 -1  2] */
            tr[1][3][0]=  invSqrt6;
            tr[1][3][1]= -invSqrt6;
            tr[1][3][2]=  twoInvSqrt6;

            /* [ 2  1  1] */
            tr[1][4][0]=  twoInvSqrt6;
            tr[1][4][1]=  invSqrt6;
            tr[1][4][2]=  invSqrt6;

            /* [ 1  2 -1] */
            tr[1][5][0]=  invSqrt6;
            tr[1][5][1]=  twoInvSqrt6;
            tr[1][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [1 -1 1] burgers vector
 */
            /* [ 1 -1 -2] */
            tr[2][0][0]=  invSqrt6;
            tr[2][0][1]= -invSqrt6;
            tr[2][0][2]= -twoInvSqrt6;

            /* [ 2  1 -1] */
            tr[2][1][0]=  twoInvSqrt6;
            tr[2][1][1]=  invSqrt6;
            tr[2][1][2]= -invSqrt6;

            /* [ 1  2  1] */
            tr[2][2][0]=  invSqrt6;
            tr[2][2][1]=  twoInvSqrt6;
            tr[2][2][2]=  invSqrt6;

            /* [-1  1  2] */
            tr[2][3][0]= -invSqrt6;
            tr[2][3][1]=  invSqrt6;
            tr[2][3][2]=  twoInvSqrt6;

            /* [-2 -1  1] */
            tr[2][4][0]= -twoInvSqrt6;
            tr[2][4][1]= -invSqrt6;
            tr[2][4][2]=  invSqrt6;

            /* [-1 -2 -1] */
            tr[2][5][0]= -invSqrt6;
            tr[2][5][1]= -twoInvSqrt6;
            tr[2][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [1 1 -1] burgers vector
 */
            /* [ 1  1  2] */
            tr[3][0][0]=  invSqrt6;
            tr[3][0][1]=  invSqrt6;
            tr[3][0][2]=  twoInvSqrt6;

            /* [ 2 -1  1] */
            tr[3][1][0]=  twoInvSqrt6;
            tr[3][1][1]= -invSqrt6;
            tr[3][1][2]=  invSqrt6;

            /* [ 1 -2 -1] */
            tr[3][2][0]=  invSqrt6;
            tr[3][2][1]= -twoInvSqrt6;
            tr[3][2][2]= -invSqrt6;

            /* [-1 -1 -2] */
            tr[3][3][0]= -invSqrt6;
            tr[3][3][1]= -invSqrt6;
            tr[3][3][2]= -twoInvSqrt6;

            /* [-2  1 -1] */
            tr[3][4][0]= -twoInvSqrt6;
            tr[3][4][1]=  invSqrt6;
            tr[3][4][2]= -invSqrt6;

            /* [-1  2  1] */
            tr[3][5][0]= -invSqrt6;
            tr[3][5][1]=  twoInvSqrt6;
            tr[3][5][2]=  invSqrt6;
        } else{
/*
 *          Vacancy [111] loops
 */
            /* [ 1  1 -2] */ tr[0][0][0]=  invSqrt6;    tr[0][0][1]=  invSqrt6;    tr[0][0][2]= -twoInvSqrt6;
            /* [-1  2 -1] */ tr[0][5][0]= -invSqrt6;    tr[0][5][1]=  twoInvSqrt6; tr[0][5][2]= -invSqrt6;
            /* [-2  1  1] */ tr[0][4][0]= -twoInvSqrt6; tr[0][4][1]=  invSqrt6;    tr[0][4][2]=  invSqrt6; 
            /* [-1 -1  2] */ tr[0][3][0]= -invSqrt6;    tr[0][3][1]= -invSqrt6;    tr[0][3][2]=  twoInvSqrt6;
            /* [ 1 -2  1] */ tr[0][2][0]=  invSqrt6;    tr[0][2][1]= -twoInvSqrt6; tr[0][2][2]=  invSqrt6;
            /* [ 2 -1 -1] */ tr[0][1][0]=  twoInvSqrt6; tr[0][1][1]= -invSqrt6;    tr[0][1][2]= -invSqrt6;

/*
 *          Vacancy [-1 1 1] loops
 */
            /* [-1  1 -2] */
            tr[1][0][0]= -invSqrt6;
            tr[1][0][1]=  invSqrt6;
            tr[1][0][2]= -twoInvSqrt6;

            /* [ 1  2 -1] */
            tr[1][1][0]=  invSqrt6;
            tr[1][1][1]=  twoInvSqrt6;
            tr[1][1][2]= -invSqrt6;

            /* [ 2  1  1] */
            tr[1][2][0]=  twoInvSqrt6;
            tr[1][2][1]=  invSqrt6;
            tr[1][2][2]=  invSqrt6;

            /* [ 1 -1  2] */
            tr[1][3][0]=  invSqrt6;
            tr[1][3][1]= -invSqrt6;
            tr[1][3][2]=  twoInvSqrt6;

            /* [-1 -2  1] */
            tr[1][4][0]= -invSqrt6;
            tr[1][4][1]= -twoInvSqrt6;
            tr[1][4][2]=  invSqrt6;

            /* [-2 -1 -1] */
            tr[1][5][0]= -twoInvSqrt6;
            tr[1][5][1]= -invSqrt6;
            tr[1][5][2]= -invSqrt6;

/*
 *          Vacancy [1 -1 1] loops
 */
            /* [ 1 -1 -2] */
            tr[2][0][0]=  invSqrt6;
            tr[2][0][1]= -invSqrt6;
            tr[2][0][2]= -twoInvSqrt6;

            /* [-1 -2 -1] */
            tr[2][1][0]= -invSqrt6;
            tr[2][1][1]= -twoInvSqrt6;
            tr[2][1][2]= -invSqrt6;

            /* [-2 -1 1] */
            tr[2][2][0]= -twoInvSqrt6;
            tr[2][2][1]= -invSqrt6;
            tr[2][2][2]=  invSqrt6;

            /* [-1  1  2] */
            tr[2][3][0]= -invSqrt6;
            tr[2][3][1]=  invSqrt6;
            tr[2][3][2]=  twoInvSqrt6;

            /* [ 1  2  1] */
            tr[2][4][0]=  invSqrt6;
            tr[2][4][1]=  twoInvSqrt6;
            tr[2][4][2]=  invSqrt6;

            /* [ 2  1 -1] */
            tr[2][5][0]=  twoInvSqrt6;
            tr[2][5][1]=  invSqrt6;
            tr[2][5][2]= -invSqrt6;

/*
 *          Vacancy [1 1 -1] loops
 */
            /* [ 1  1  2] */
            tr[3][0][0]=  invSqrt6;
            tr[3][0][1]=  invSqrt6;
            tr[3][0][2]=  twoInvSqrt6;

            /* [-1  2  1] */
            tr[3][1][0]= -invSqrt6;
            tr[3][1][1]=  twoInvSqrt6;
            tr[3][1][2]=  invSqrt6;

            /* [-2  1 -1] */
            tr[3][2][0]= -twoInvSqrt6;
            tr[3][2][1]=  invSqrt6;
            tr[3][2][2]= -invSqrt6;

            /* [-1 -1 -2] */
            tr[3][3][0]= -invSqrt6;
            tr[3][3][1]= -invSqrt6;
            tr[3][3][2]= -twoInvSqrt6;

            /* [ 1 -2 -1] */
            tr[3][4][0]=  invSqrt6;
            tr[3][4][1]= -twoInvSqrt6;
            tr[3][4][2]= -invSqrt6;

            /* [ 2 -1  1] */
            tr[3][5][0]=  twoInvSqrt6;
            tr[3][5][1]= -invSqrt6;
            tr[3][5][2]=  invSqrt6;
        }

/*
 *      FIX ME!  Need to modify the code to deal with <100> type burgers
 *               vector later
 */

        inData->nodeCount = 0;
        nextNode = 0;
        burgIndex = maxBurgIndex;

/*
 *      Create one loop at a time, cycling through burgers vectors as we go.
 */
        for (loopIndex = 0; loopIndex < numLoops; loopIndex++) {

            if (++burgIndex > maxBurgIndex) {
                burgIndex = minBurgIndex;
            }

/*
 *          Increase the size of the node array enough to hold
 *          all the new nodes created for this loop.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numSegs;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numSegs);

            loopCtrX = (randm(&seed)-0.5) * cubeSize;
            loopCtrY = (randm(&seed)-0.5) * cubeSize;
            loopCtrZ = (randm(&seed)-0.5) * cubeSize;

            for (id = nextNode; id < (nextNode+numSegs); id++) {

                node =  &inData->node[id];

/*
 *              Pick the index for the direction component of the tr array
 */
                dIndex = (id - nextNode) % numSegs;

                node->x = loopCtrX + radius * tr[burgIndex][dIndex][0];
                node->y = loopCtrY + radius * tr[burgIndex][dIndex][1];
                node->z = loopCtrZ + radius * tr[burgIndex][dIndex][2];

/*
 *              Set up the node
 */
                node->constraint = 0;
                node->myTag.domainID = dislocType;
                node->myTag.index = id;

                AllocNodeArms(node, 2);

                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];

                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = nextNode + ((id+1)%numSegs);

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = nextNode + ((id+numSegs-1)%numSegs);
            }

/*
 *          Need to set the glide plane normal for each segment.
 *          Couldn't do it above because we didn't have all the
 *          neighbor node positions until all the nodes in the
 *          loop were created.
 */
            for (id = nextNode; id < (nextNode+numSegs); id++) {

                node = &inData->node[id];

                nbr1Index = node->nbrTag[1].index;
                nbr1Node = &inData->node[nbr1Index];

                ux = nbr1Node->x - node->x;
                uy = nbr1Node->y - node->y;
                uz = nbr1Node->z - node->z;

                Normalize(&ux,&uy,&uz);

/* 
 *              l cross b gives normal vector for each segment, they may
 *              not be (110)
 */
                vec1[0] = ux;
                vec1[1] = uy;
                vec1[2] = uz;

                vec2[0] =  node->burgX[1];
                vec2[1] =  node->burgY[1];
                vec2[2] =  node->burgZ[1];

                NormalizedCrossVector(vec1, vec2, tmpVec);

                node->nx[1] = tmpVec[0];
                node->ny[1] = tmpVec[1];
                node->nz[1] = tmpVec[2];

                nbr1Node->nx[0] = tmpVec[0];
                nbr1Node->ny[0] = tmpVec[1];
                nbr1Node->nz[0] = tmpVec[2];

            }  /*  for (id = nextNode; ...) */

            nextNode += numSegs; 
            lastBlock = (loopIndex == (numLoops-1));
            
/*
 *          When we've generated the nodal data for the final loop,
 *          write all the data out to disk.
 */
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                InitRemesh(inData, dislocType, startRemeshIndex);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
                nextNode = 0;
            }

        }  /* for (loopIndex = 0; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateScrewConfig
 *      Description:  The function name says it all
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateScrewConfig(Home_t *home, InData_t *inData, int cubeLength,
                       int numChains, int seed, real8 *totDislocLen,
                       int dislocType)
{
        int      ic, np, ip, id0;
        int      burgIndex, gpIndex, newNodeIndex, lastBlock;
        int      nbr1Index, nbr2Index;
        int      startRemeshIndex = 0;
        real8    xp[3], yp[3], zp[3], cubeSize;
        real8    burg[4][3], glidePlane[4][6][3];
        Param_t  *param;
        Node_t   *node;

        param = inData->param;
        cubeSize = (real8)cubeLength;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateScrewConfig", numChains);
        }

/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */

/*
 *      Set up an array of valid burgers vectors
 */
        inData->nburg = 4;

        burg[0][0] =  0.5773503;
        burg[0][1] =  0.5773503;
        burg[0][2] =  0.5773503; 

        burg[1][0] = -0.5773503;
        burg[1][1] =  0.5773503;
        burg[1][2] =  0.5773503; 

        burg[2][0] =  0.5773503;
        burg[2][1] = -0.5773503;
        burg[2][2] =  0.5773503; 

        burg[3][0] =  0.5773503;
        burg[3][1] =  0.5773503;
        burg[3][2] = -0.5773503; 

/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

        id0 = 0;
        inData->nodeCount = 0;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

            burgIndex = ic%4;
            gpIndex = (ic / 4) % 6;

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] - (cubeSize * burg[burgIndex][0]);
            yp[0] = yp[1] - (cubeSize * burg[burgIndex][1]);
            zp[0] = zp[1] - (cubeSize * burg[burgIndex][2]);

            xp[2] = xp[1] + (cubeSize * burg[burgIndex][0]);
            yp[2] = yp[1] + (cubeSize * burg[burgIndex][1]);
            zp[2] = zp[1] + (cubeSize * burg[burgIndex][2]);

/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = id0 + nbr1Index;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[burgIndex][gpIndex][0];
                node->ny[0] = glidePlane[burgIndex][gpIndex][1];
                node->nz[0] = glidePlane[burgIndex][gpIndex][2];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = id0 + nbr2Index;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[burgIndex][gpIndex][0];
                node->ny[1] = glidePlane[burgIndex][gpIndex][1];
                node->nz[1] = glidePlane[burgIndex][gpIndex][2];

            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the current block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == numChains - 1);
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCIrradConfig
 *      Description:  Generate an FCC problem space fill with a combination
 *                    of line dislocations and Frank Sessile Hexagonal
 *                    Interstitial Loops.
 *                    (Masato Hiratani)
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *          numLoops    Number of hexagonal loops to create
 *          hexl
 *
 *-------------------------------------------------------------------------*/
void CreateFCCIrradConfig(Home_t *home, InData_t *inData, int cubeLength,
                          int numChains, int seed, int numLoops, real8 hexl,
                          real8 *totDislocLen, int dislocType)
{
        int      ic, np, ip, id0, i;
        int      nplane, indp, indb, indf, inds, indr;
        int      newNodeIndex, lastBlock;
        int      startRemeshIndex = 0;
        real8    inv3, inv6, invsq2, sq2over3;
        real8    xp[6], yp[6], zp[6], cubeSize;
        real8    tnx[4], tny[4], tnz[4], burg[12][3];
        Node_t   *node;
        Param_t  *param;

        if (numLoops <= 0) {
            Fatal("%s: numLoops is %d, but must be > 0.\n",
                  "CreateFCCIrradConfig", numLoops);
        }

        param = inData->param;
        cubeSize = (real8)cubeLength;

/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */

        inv3     = 0.3333333;
        inv6     = 0.1666666;
        invsq2   = 0.70710678;
        sq2over3 = 0.81649658;

        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta plane normals
 */
        nplane = 4;
 
        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV is parallel to the plane normal i.e. FS 1/3<111>
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0;
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i];

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i];

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0;

            Normalize(&tnx[i], &tny[i], &tnz[i]);
        }

        id0 = 0;
        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {

            np = 3;    /* number of points along 1 chain */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;                /* index of plane */
            indb = indp*3+ic%3;         /* index of BV */
            indf = ((ic-ic%12)/12)%2;   /* index of alternative line sense */
            inds = indp*3+(ic+indf+1)%3;        /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2);     /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] + indr * cubeSize * burg[inds][0];
            yp[1] = yp[0] + indr * cubeSize * burg[inds][1];
            zp[1] = zp[0] + indr * cubeSize * burg[inds][2];

            xp[2] = xp[1] + indr * cubeSize * burg[inds][0];
            yp[2] = yp[1] + indr * cubeSize * burg[inds][1];
            zp[2] = zp[1] + indr * cubeSize * burg[inds][2];

/*
 *          PBC
 */        
            for (i = 0; i < 3; i++) {
                if (xp[i] < -cubeSize/2) xp[i] += cubeSize;
                if (yp[i] < -cubeSize/2) yp[i] += cubeSize;
                if (zp[i] < -cubeSize/2) zp[i] += cubeSize;
            }
/*
            printf("ic indp indb indf inds indr %d %d %d %d %d %d\n",
                   ic,indp,indb,indf,inds,indr);
*/
            for (ip = 0; ip < np; ip++) {
                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = PINNED_NODE;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;
                node->burgX[0] = burg[indb][0];
                node->burgY[0] = burg[indb][1];
                node->burgZ[0] = burg[indb][2];
                node->nx[0] = tnx[indp];
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;
                node->burgX[1] = -burg[indb][0];
                node->burgY[1] = -burg[indb][1];
                node->burgZ[1] = -burg[indb][2];
                node->nx[1] = tnx[indp];
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Then, if the count of nodes currently contained in memoru
 *          exceeds the threshhold write the current block of nodal data
 *          to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

/*
 *          We won't write nodal data out here, but rather after the
 *          hexagonal loops are added below.
 */
            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* end of chains */

/*
 *      Place hexagonal loops
 */
        for (ic = 0; ic < numLoops; ic++) {
            np = 6;     /* creation of hexagonal loop */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;              /* index of plane */
            indb = indp*3;            /* index of BV */
            indf = ((ic-ic%12)/12)%2; /* index of alternativeline sense */
            inds = indp*3+(ic+indf+1)%3;    /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2); /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] - inv6*hexl*burg[indb][0];
            yp[1] = yp[0] - inv6*hexl*burg[indb][1];
            zp[1] = zp[0] - inv6*hexl*burg[indb][2];
            xp[2] = xp[1] + inv6*hexl*burg[indb+2][0];
            yp[2] = yp[1] + inv6*hexl*burg[indb+2][1];
            zp[2] = zp[1] + inv6*hexl*burg[indb+2][2];
            xp[3] = xp[2] - inv6*hexl*burg[indb+1][0];
            yp[3] = yp[2] - inv6*hexl*burg[indb+1][1];
            zp[3] = zp[2] - inv6*hexl*burg[indb+1][2];
            xp[4] = xp[3] + inv6*hexl*burg[indb][0];
            yp[4] = yp[3] + inv6*hexl*burg[indb][1];
            zp[4] = zp[3] + inv6*hexl*burg[indb][2];
            xp[5] = xp[4] - inv6*hexl*burg[indb+2][0];
            yp[5] = yp[4] - inv6*hexl*burg[indb+2][1];
            zp[5] = zp[4] - inv6*hexl*burg[indb+2][2];
      
/*
 *          PBC
 */
            for (i = 0; i < np; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }

            for (ip = 0; ip < np; ip++) {
                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = PINNED_NODE;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;
                node->burgX[0] = sq2over3*tnx[indp];
                node->burgY[0] = sq2over3*tny[indp];
                node->burgZ[0] = sq2over3*tnz[indp];
                node->nx[0] = tnx[indp];
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;
                node->burgX[1] = -sq2over3*tnx[indp];
                node->burgY[1] = -sq2over3*tny[indp];
                node->burgZ[1] = -sq2over3*tnz[indp];
                node->nx[1] = tnx[indp];
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
/*
                printf("lnode(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index,
                       node->burgX[0],node->burgY[0],node->burgZ[0],
                       node->burgX[1],node->burgY[1],node->burgZ[1]);
                printf("lnode(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index,
                       node->nx[0],node->ny[0],node->nz[0],
                       node->nx[1],node->ny[1],node->nz[1]);
 */     
            }   
/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == (numLoops - 1));
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* end of loops */

        return; 
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCConfig
 *      Description:  Generates random line configurations compatible
 *                    with MobilityRule_FCC1.  (Masato Hiratani)
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateFCCConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
						 real8 *totDislocLen)
{
        int      i, lineIndex, chain, baseNodeID, numPoints,ranNum;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      x0, d, nLayers, xn, layer;
        int      startRemeshIndex = 0;
        real8    srcLen, thickness;
        real8    p0[3], p1[3], p2[3], p3[3], randomLine[3],k0[3];
        real8    p4[3], p5[3];
        real8	 randomAngle,srcLenMin, srcLenMax,lenRange;
        int	     cubeLength, numSources,  seed, dislocType;
        int	     polycrystal, nGrain, finish;
        int      firstBurgIndex, srcLenType;
        int      inGrain1, onGB1, inGrain2, onGB2, inGrain0, onGB0;
        int      inGrain3, onGB3;
        real8	 twinRatio, maxSeg, maxDensity, burgMag = 2.556E-10;
        real8	 volume, glideDir[3];
        real8	 tnx[4], tny[4], tnz[4], burg[12][3], cubeSize,heightSize;
        real8	 n0[3], b[3], n1[3];	/* plane normal vector, burgers, normal vector*/
        Node_t   *node;
        Param_t  *param;
        real8	 invsq2, basePos[3], b1[3], b2[3], f1[3], f2[3], f3[3], f4[3];
        real8	 rotMatrix[3][3],rotLine[3],rotBurg[3],rotPlane[3], btemp[3];
        int      intID, inGrain4, onGB4, inGrain5, onGB5;
        real8    pi = 3.1415926535897932385, totSrcLen;
        real8    separateLength = 40.0;

        bool        finished;
        BurgInfo_t  *bd;
        	
        real8                                       str[3][3];
        int                                         appliedDir;
        real8                                       sign;
        std::vector<std::vector<Dislocation_t> >    disArray;
        std::vector<std::vector<real8> >            burgList;

        real8    normalABC[3];
        normalABC[0] = 1;
        normalABC[1] = 1;
        normalABC[2] = 1;
        VECTOR_ZERO(str[0]);
        VECTOR_ZERO(str[1]);
        VECTOR_ZERO(str[2]);

        for(i=0;i<3;i++){
            str[i][i] = inArgs->loadDir[i];
        }

/*
 *      Grap the command-line argument values needed for this function
 */
        home->param->materialType = MAT_TYPE_FCC;
        	
        cubeLength  = inArgs->cubeLength;
	    srcLenType  = inArgs->srcLenType;
        nLayers = inArgs->nLayers;
       	maxDensity  = inArgs->maxDensity;
        numSources  = inArgs->numFRSrcs;
        srcLenMin   = inArgs->frLenMin;
        srcLenMax   = inArgs->frLenMax;
        seed        = inArgs->seed;
        dislocType  = inArgs->type;
        polycrystal = inArgs->polycrystal;
        nGrain	    = inArgs->nGrain;
        twinRatio   = inArgs->twinRatio;
        maxSeg      = (real8)inArgs->maxSegLen;
        	
        d = inArgs->deltaNums;
		
        lenRange = srcLenMax - srcLenMin;
		
        param = inData->param;
        cubeSize	=	(real8)cubeLength;
        	
        n1[0] = 0.0;	n1[1] = 0.0;	n1[2] = 0.0;

        srcLen = cubeSize;		
	    printf("the max density is %e \n", maxDensity);
/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        baseNodeID = 0;
        inData->nodeCount = 0;
		
        heightSize = 2 * param->zBoundMax;
	    volume = heightSize * cubeSize * cubeSize;
        thickness = heightSize / ((real8)nLayers);
        	
        x0 = (d*nLayers*(1 - nLayers) + 2*numSources)/(2 * nLayers);

	    totSrcLen = 0.0;
	    finish = 0;
        for (layer = 0; layer < nLayers; layer++) {
	    if(finish)break;
            xn = layer * d + x0;
            if(xn <= 0){
                xn = 0;
                continue;
            }
            	
            printf("the num of  sources in layers %d is %d\n",layer, xn);
            	
            for(chain = 0; chain < xn; ){
             
                numPoints = (inArgs->extended == 0) ? 3 : 6;
   
                finished = 0;
                while(!finished){
                    p0[0] = (randm(&seed)-0.5) * cubeSize;
                    p0[1] = (randm(&seed)-0.5) * cubeSize;			
                    p0[2] = randm(&seed)* thickness + ((real8)layer) * thickness - 0.5*heightSize;
                    NodeLocation(home, p0[0],  p0[1], p0[2], inGrain0, onGB0,15);
                    bd = GetBurgData(home, inGrain0, onGB0);
                    if(onGB0 > -1)continue;

                    ranNum = chain;
                    gpBurgIndex = (ranNum % 12) / 2;
                    if((bd->burgType[gpBurgIndex] & EDGE_BURG) == 0)continue;

                    gpIndex = bd->burgOnPlaneID[gpBurgIndex][ranNum % 2];
                    if((bd->planeType[gpIndex] & GLIDE_PLANE) == 0 )continue;
   
                    VECTOR_COPY(n0, bd->planeList[gpIndex]);
                    VECTOR_COPY(b , bd->burgList[gpBurgIndex]);
                    if(!PerpendicularDirection(b,param->TBNormal)){break;}
                    VECTOR_COPY(randomLine, b);

                    for (i = 0; i < 3; i++) {
                        p1[i] = p0[i] + (srcLen * randomLine[i]*0.5);
                        p2[i] = p0[i] - (srcLen * randomLine[i]*0.5);
                    } 

                    NodeLocation(home, p1[0],  p1[1], p1[2], inGrain1, onGB1,15);
                    NodeLocation(home, p2[0],  p2[1], p2[2], inGrain2, onGB2,15);
             
                    if(onGB1 > -1 || inGrain0 != inGrain1 || onGB2 > -1 || inGrain2 != inGrain1){
                        continue;
                    }	                

                    if(inArgs->extended != 0){
                        VECTOR_COPY(basePos, p0);
                        cross(n0, randomLine, glideDir);
                        NormalizeVec(glideDir);
                        p0[0] -= (separateLength*glideDir[0]*0.5);
                        p0[1] -= (separateLength*glideDir[1]*0.5);
                        p0[2] -= (separateLength*glideDir[2]*0.5);

                        p3[0] = p0[0] + separateLength*glideDir[0];
                        p3[1] = p0[1] + separateLength*glideDir[1];
                        p3[2] = p0[2] + separateLength*glideDir[2];

                        p1[0] -= (separateLength*glideDir[0]*0.5);
                        p1[1] -= (separateLength*glideDir[1]*0.5);
                        p1[2] -= (separateLength*glideDir[2]*0.5);

                        p4[0] = p1[0] + separateLength*glideDir[0];
                        p4[1] = p1[1] + separateLength*glideDir[1];
                        p4[2] = p1[2] + separateLength*glideDir[2];

                        p2[0] -= (separateLength*glideDir[0]*0.5);
                        p2[1] -= (separateLength*glideDir[1]*0.5);
                        p2[2] -= (separateLength*glideDir[2]*0.5);

                        p5[0] = p2[0] + separateLength*glideDir[0];
                        p5[1] = p2[1] + separateLength*glideDir[1];
                        p5[2] = p2[2] + separateLength*glideDir[2];
#if 1
                        NodeLocation(home, p0[0],  p0[1], p0[2], inGrain0, onGB0,15);
                        NodeLocation(home, p1[0],  p1[1], p1[2], inGrain1, onGB1,15);
                        NodeLocation(home, p2[0],  p2[1], p2[2], inGrain2, onGB2,15);
                        NodeLocation(home, p3[0],  p3[1], p3[2], inGrain3, onGB3,15);
                        NodeLocation(home, p4[0],  p4[1], p4[2], inGrain4, onGB4,15);
                        NodeLocation(home, p5[0],  p3[1], p5[2], inGrain5, onGB5,15);
                    
                        if(onGB0 > -1 || onGB3 > -1 || onGB1 > -1 || onGB2 > -1 ||
                           onGB4 > -1 || onGB5 > -1 || inGrain4 != inGrain0 || inGrain5 != inGrain0 ||
                          inGrain0 != inGrain2 || inGrain3 != inGrain2 || inGrain0 != inGrain1)continue;
#endif
                    }

                    finished = 1;
                }

                if(inArgs->extended != 0){
                    if(!GeometricCriterion(bd, b, randomLine, n0, 2,
                       burgList, disArray))continue;

                    if(disArray.size() != 1)printf("Warning: the extened dislcoation has some wrong\n");
                    if(disArray[0].size() != 2)printf("Warning: the extened dislcoation has some wrong\n");
                    
                    VECTOR_COPY(b1, burgList[disArray[0][0].burgID]);
                    VECTOR_COPY(b2, burgList[disArray[0][1].burgID]);

                    appliedDir = abs(inArgs->extended)-1; 
                    sign = (inArgs->extended > 0) ? 1.0 : -1.0;
                    VECTOR_ZERO(str[0]);
                    VECTOR_ZERO(str[1]);
                    VECTOR_ZERO(str[2]);
                    if(appliedDir < 3){
                        str[appliedDir][appliedDir] = sign * 1.0E8;
                    }else{
                        switch(appliedDir){
                            case 3:
                                str[1][2] = sign * 1.0E8;
                                str[2][1] = str[1][2];
                                break;
                            case 4:
                                str[0][2] = sign * 1.0E8;
                                str[2][0] = str[0][2];
                                break;
                            case 5:
                                str[0][1] = sign * 1.0E8;
                                str[1][0] = str[0][1];
                                break;
                            default:
                                Fatal("Wrong extend type");
                                break;
                        }
                    }

                    ExtPKForce(str, b[0], b[1], b[2], 
                               basePos[0], basePos[1], basePos[2],
                               p1[0], p1[1], p1[2], f1, f2);
                    ExtPKForce(str, b1[0], b1[1], b1[2], 
                               basePos[0], basePos[1], basePos[2],
                               p1[0], p1[1], p1[2], f3, f4);
                    if(DotProduct(f1, f3) < 0){
                        VECTOR_COPY(btemp, b1);
                        VECTOR_COPY(b1, b2);
                        VECTOR_COPY(b2, btemp);
                        ExtPKForce(str, b1[0], b1[1], b1[2], 
                                   basePos[0], basePos[1], basePos[2],
                                   p1[0], p1[1], p1[2], f3, f4);
                        if(DotProduct(f1, f3) < 0)printf("Warning: in %s at %d", __FILE__, __LINE__);
                    }
                }
                if(finished == 0){chain++;continue;}

	            totSrcLen += srcLen;				
/*           
 *              Reallocate the node array with sufficient size to add
 *              all the new nodes defining this chain.
 */          
                newNodeIndex = inData->nodeCount;
                inData->nodeCount += numPoints;
                inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                                   * sizeof(Node_t));
                memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
		 	   
                if(inArgs->extended == 0){
/*           
 *                   Set up the 3 nodes we're using to define the dislocation line.
 *                   First do point p0.
 */               
                     node = &inData->node[baseNodeID];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID;
                     node->x = p0[0];
                     node->y = p0[1];
                     node->z = p0[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain0;
                     node->onGB = onGB0;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = b[0];
                     node->burgY[0] = b[1];
                     node->burgZ[0] = b[2];
                  
                     node->burgX[1] = -b[0];
                     node->burgY[1] = -b[1];
                     node->burgZ[1] = -b[2];
                  
/*                
 *                   For screw dislocations, use the glide plane from the table.
 *                   For edge dislocations, glide plane is the cross product
 *                   of the burgers vector and line direction.
 */               
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1;
                  
                     node->x = p1[0];
                     node->y = p1[1];
                     node->z = p1[2];
                  
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain1;
                     node->onGB = onGB1;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = -b[0];
                     node->burgY[0] = -b[1];
                     node->burgZ[0] = -b[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                  
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+2;
                  
                     node->burgX[1] =  b[0];
                     node->burgY[1] =  b[1];
                     node->burgZ[1] =  b[2];
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2;
                  
                     node->x = p2[0];
                     node->y = p2[1];
                     node->z = p2[2];
                  
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain2;
                     node->onGB = onGB2;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = b[0];
                     node->burgY[0] = b[1];
                     node->burgZ[0] = b[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];

                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+1;
                  
                     node->burgX[1] =-b[0];
                     node->burgY[1] =-b[1];
                     node->burgZ[1] =-b[2];
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                }else{
                     node = &inData->node[baseNodeID];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID;
                     node->x = p0[0];
                     node->y = p0[1];
                     node->z = p0[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain0;
                     node->onGB = onGB0;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = b1[0];
                     node->burgY[0] = b1[1];
                     node->burgZ[0] = b1[2];
                     node->signSF[0] = 1.0;
                  
                     node->burgX[1] = -b1[0];
                     node->burgY[1] = -b1[1];
                     node->burgZ[1] = -b1[2];
                     node->signSF[1] = -1.0;
                  
/*                
 *                   For screw dislocations, use the glide plane from the table.
 *                   For edge dislocations, glide plane is the cross product
 *                   of the burgers vector and line direction.
 */               
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1;
                  
                     node->x = p1[0];
                     node->y = p1[1];
                     node->z = p1[2];
                  
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain1;
                     node->onGB = onGB1;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = -b1[0];
                     node->burgY[0] = -b1[1];
                     node->burgZ[0] = -b1[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     node->signSF[0] = -1.0;
                  
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+2;
                  
                     node->burgX[1] =  b1[0];
                     node->burgY[1] =  b1[1];
                     node->burgZ[1] =  b1[2];
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                     node->signSF[1] = 1.0;
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2;
                  
                     node->x = p2[0];
                     node->y = p2[1];
                     node->z = p2[2];
                  
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain2;
                     node->onGB = onGB2;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = b1[0];
                     node->burgY[0] = b1[1];
                     node->burgZ[0] = b1[2];
                     node->signSF[0] = 1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];

                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+1;
                  
                     node->burgX[1] =-b1[0];
                     node->burgY[1] =-b1[1];
                     node->burgZ[1] =-b1[2];
                     node->signSF[1] = -1.0;
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
//
                     node = &inData->node[baseNodeID+3];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+3;
                     node->x = p3[0];
                     node->y = p3[1];
                     node->z = p3[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain3;
                     node->onGB = onGB3;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1 + 3;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2 + 3;
                  
                     node->burgX[0] = b2[0];
                     node->burgY[0] = b2[1];
                     node->burgZ[0] = b2[2];
                     node->signSF[0] = -1.0;
                  
                     node->burgX[1] = -b2[0];
                     node->burgY[1] = -b2[1];
                     node->burgZ[1] = -b2[2];
                     node->signSF[1] = 1.0;
                  
/*                
 *                   For screw dislocations, use the glide plane from the table.
 *                   For edge dislocations, glide plane is the cross product
 *                   of the burgers vector and line direction.
 */               
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1+3];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1+3;
                  
                     node->x = p4[0];
                     node->y = p4[1];
                     node->z = p4[2];
                  
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain4;
                     node->onGB = onGB4;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID+3;
                  
                     node->burgX[0] = -b2[0];
                     node->burgY[0] = -b2[1];
                     node->burgZ[0] = -b2[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     node->signSF[0] = 1.0;
                  
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+2+3;
                  
                     node->burgX[1] =  b2[0];
                     node->burgY[1] =  b2[1];
                     node->burgZ[1] =  b2[2];
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                     node->signSF[1] = -1.0;
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2+3];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2+3;
                  
                     node->x = p5[0];
                     node->y = p5[1];
                     node->z = p5[2];
                  
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain5;
                     node->onGB = onGB5;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID+3;
                  
                     node->burgX[0] = b2[0];
                     node->burgY[0] = b2[1];
                     node->burgZ[0] = b2[2];
                     node->signSF[0] = -1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];

                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+1+3;
                  
                     node->burgX[1] =-b2[0];
                     node->burgY[1] =-b2[1];
                     node->burgZ[1] =-b2[2];
                     node->signSF[1] = 1.0;
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                }
/*           
 *              The initial segments created are not necessarily limited to
 *              param->maxSegLen, so a call to InitRemesh() is needed to
 *              chop any excessively long segments into proper lengths.
 *              When we've generated the nodal data for the final chain,
 *              write the block of nodal data to the file.
 */	         
                InitRemesh(inData, dislocType, startRemeshIndex);
//                lastBlock = (chain == (numSources - 1));
             
               
	            if(totSrcLen > maxDensity*volume*burgMag*burgMag){
		            finish = 1;
		            break;
	            }else{
		            finish = 0;
	            }
                chain++;
                baseNodeID = inData->nodeCount;
                startRemeshIndex = baseNodeID;

            }  /* for (chain = 0; chain < numChains; ) */
		
	    }/* for(layer=0; layer<nLayers;...)*/
	
	    IncDislocationDensity(inData, totDislocLen);
	    param->nodeCount = inData->nodeCount;
	    WriteInitialNodeData(home, inData, 1);
	    FreeInNodeArray(inData, inData->nodeCount);
	    inData->nodeCount = 0;
	
        return;
}


/*---------------------------------------------------------------------------
 *
 *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *            @                                    @
 *            @ THIS FUNCTION IS NOT YET COMPLETE! @
 *            @                                    @
 *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *
 *      Function:     CreateFCCPerfectLoop
 *      Description:  Generate an initial set of nodal data consisting
 *                    of Glissile perfect loops for FCC.  (Masato Hiratani)
 *
 *      Arguments:
 *          cubeLength Length of cubic problem space in a single
 *                     dimension (units of b)
 *          numChains  Number of chains to create
 *          seed       Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateFCCPerfectLoop(Home_t *home, InData_t *inData, int cubeLength,
                          int numChains, int seed, real8 *totDislocLen,
                          int dislocType)
{
        int      nplane, indexplane, ic, np, ip, id0, i;
        int      inds11, inds12, inds21, inds22;
        int      indp, indp1, indp2, indb1, indb2, indr;
        int      newNodeIndex, lastBlock;
        int      startRemeshIndex = 0;
        real8    prisml, invsq2, cubeSize;
        real8    xp[4], yp[4], zp[4];
        real8    tempbx,tempby,tempbz,tnx[4], tny[4], tnz[4];
        real8    sxp,syp,szp, bvx, bvy, bvz, siz;
        real8    burg[12][3];
        real8    vec1[3], vec2[3], tmpVec[3];
        Param_t  *param;
        Node_t   *node;

        param = inData->param;
        invsq2 = 0.70710678;
        prisml = 0.3535533*cubeLength;        /* size of loop */
        cubeSize = (real8)cubeLength;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be a multiple of 12.\n",
                  "CreateFCCPerfectLoop", numChains);
        }
    
/*
 *      The nodal data is generated in a single block.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes. The InitRemesh() calls also accept the dislocation
 *      type to set the tag.domainID value for any new nodes it adds to the
 *      node array.
 */
        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta planes
 */
        nplane =  4;

        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0; 
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i]; 

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i]; 

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0; 

            Normalize(&tnx[i],&tny[i],&tnz[i]);   
        }
    
        id0 = 0;
        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {

            np = 4;    /* number of points along 1 chain */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *) realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          BV cycle 6, determination of zone axis
 */
            indp = ic%6;  /* zone axis */
            if (indp == 0) { indp1 = 0; indp2 = 1;}
            if (indp == 1) { indp1 = 0; indp2 = 2;}
            if (indp == 2) { indp1 = 0; indp2 = 3;}
            if (indp == 3) { indp1 = 1; indp2 = 2;}
            if (indp == 4) { indp1 = 1; indp2 = 3;}
            if (indp == 5) { indp1 = 2; indp2 = 3;}
/*
            printf("n1= %e %e %e\n", tnx[indp1], tny[indp1], tnz[indp1]);
            printf("n2= %e %e %e\n", tnx[indp2], tny[indp2], tnz[indp2]);
*/
            vec1[0] = tnx[indp1];
            vec1[1] = tny[indp1];
            vec1[2] = tnz[indp1];

            vec2[0] = tnx[indp2];
            vec2[1] = tny[indp2];
            vec2[2] = tnz[indp2];

            NormalizedCrossVector(vec1, vec2, tmpVec);

            tempbx = tmpVec[0];
            tempby = tmpVec[1];
            tempbz = tmpVec[2];

            if (tempbx == 0) {
                indb1  = indp1*3;
                inds11 = indb1+1;
                inds12 = indb1+2;
                indb2  = indp2*3;
                inds21 = indb2+1;
                inds22 = indb2+2;
            } else if (tempby == 0) {
                indb1  = indp1*3+1;
                inds11 = indb1+1;
                inds12 = indb1-1;
                indb2  = indp2*3+1;
                inds21 = indb2+1;
                inds22 = indb2-1;
            } else if (tempbz == 0) {
                indb1  = indp1*3+2;
                inds11 = indb1-2;
                inds12 = indb1-1;
                indb2  = indp2*3+2;
                inds21 = indb2-2;
                inds22 = indb2-1;
            } else {    
                printf("n1xn2= %e %e %e\n", tempbx, tempby, tempbz);
                Fatal("CreateFCCPerfectLoop: wrong zone BV");
            }

            indr = 1-2*(((ic-ic%6)/6)%2); /* sign change: I-loop & V-loop */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            if (indr > 0) {
                xp[1] = xp[0] + indr*prisml*burg[inds11][0];
                yp[1] = yp[0] + indr*prisml*burg[inds11][1];
                zp[1] = zp[0] + indr*prisml*burg[inds11][2];
                xp[2] = xp[1] + indr*prisml*burg[inds12][0];
                yp[2] = yp[1] + indr*prisml*burg[inds12][1];
                zp[2] = zp[1] + indr*prisml*burg[inds12][2];
                sxp   = xp[2] + indr*prisml*burg[inds21][0];
                syp   = yp[2] + indr*prisml*burg[inds21][1];
                szp   = zp[2] + indr*prisml*burg[inds21][2];
                xp[3] = xp[2] - indr*prisml*burg[inds11][0];
                yp[3] = yp[2] - indr*prisml*burg[inds11][1];
                zp[3] = zp[2] - indr*prisml*burg[inds11][2];
                indexplane = indp1;
             } else {
                xp[1] = xp[0] + indr*prisml*burg[inds21][0];
                yp[1] = yp[0] + indr*prisml*burg[inds21][1];
                zp[1] = zp[0] + indr*prisml*burg[inds21][2];
                xp[2] = xp[1] + indr*prisml*burg[inds22][0];
                yp[2] = yp[1] + indr*prisml*burg[inds22][1];
                zp[2] = zp[1] + indr*prisml*burg[inds22][2];
                sxp   = xp[2] + indr*prisml*burg[inds11][0];
                syp   = yp[2] + indr*prisml*burg[inds11][1];
                szp   = zp[2] + indr*prisml*burg[inds11][2];
                xp[3] = xp[2] - indr*prisml*burg[inds21][0];
                yp[3] = yp[2] - indr*prisml*burg[inds21][1];
                zp[3] = zp[2] - indr*prisml*burg[inds21][2];
                indexplane = indp2;
            }

            GetUnitVector(1, sxp, syp, szp, xp[1], yp[1], zp[1],
                          &bvx, &bvy, &bvz, &siz); 

            printf("ic,indb1 indr indexplane=\n");
            printf("%d %d %d %d \n", ic, indb1, indr, indexplane);    
            printf("n= %e %e %e\n", tnx[indexplane],
                   tny[indexplane], tnz[indexplane]);    
            printf("b= %e %e %e\n", bvx,bvy,bvz);

/*
 *          PBC
 */        
            for (i = 0; i < np; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }
/*
            printf("ic indp indr %d %d %d\n",ic,indp, indr);
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+id0;

                node->burgX[0] = bvx;
                node->burgY[0] = bvy;
                node->burgZ[0] = bvz;

                node->nx[0] = 0.0;
                node->ny[0] = 0.0;
                node->nz[0] = 0.0;


                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+id0;

                node->burgX[1] = -bvx;
                node->burgY[1] = -bvy;
                node->burgZ[1] = -bvz;

                node->nx[1] = 0.0;
                node->ny[1] = 0.0;
                node->nz[1] = 0.0;

/*
                printf("node(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index, 
                       node->burgX[0],node->burgY[0],node->burgZ[0],
                       node->burgX[1],node->burgY[1],node->burgZ[1]);
                printf("node(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
                       node->myTag.domainID, node->myTag.index, 
                       node->nx[0],node->ny[0],node->nz[0],
                       node->nx[1],node->ny[1],node->nz[1]);
 */   
            }
/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (ic == (numChains - 1));
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;

        } /* loop over chains */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFiniteMixedConfig
 *      Description:  Creates random screw and edge dislocations,
 *                    but assumes a finite problem space in at
 *                    least 1 dimension and terminates the dislocations
 *                    at the free surface(s) 
 *
 *      Arguments:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numChains     Number of chains to create
 *          seed          Seed value for random number generator
 *          totDislocLen  Pointer to location at which to return
 *                        to caller the total combined length of
 *                        all created dislocations.
 *
 *-------------------------------------------------------------------------*/
void CreateFiniteMixedConfig(Home_t *home, InData_t *inData, int cubeLength,
                             int numChains, int seed, real8 *totDislocLen,
                             int dislocType)
{
        int      i, lineIndex, chain, baseNodeID, numPoints;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      startRemeshIndex = 0;
        int      signFact; 
        int      nbr1Index, nbr2Index, numConnections, numSegs;
        int      intersectIndex1, intersectIndex2;
        int      pbc[3];
        real8    burgSign, vecLen;
        real8    surfCoord, len, minLen1, minLen2;
        real8    intersectSurfCoord1, intersectSurfCoord2;
        real8    p0[3], newPos[3];
        real8    range[3], rangeCntr[3];
        real8    minCoord[3], maxCoord[3];
        real8    intersectPos1[3], intersectPos2[3], vector[3];
        real8    normalizedBurg[3];
        real8    linedir[16][3], glidePlane[4][6][3], plane[3];
        Node_t   *node;
        Param_t  *param;
        
        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateFiniteMixedConfig", numChains);
        }
        
        param = inData->param;

        pbc[0] = (param->xBoundType == Periodic);
        pbc[1] = (param->yBoundType == Periodic);
        pbc[2] = (param->zBoundType == Periodic);

/*
 *      Just a sanity check...
 */
        if (pbc[0]*pbc[1]*pbc[2] != 0) {
            Fatal("CreateFiniteMixedConfig: Requires free surfaces in at"
                  "least one dimension\n    but all boundaries are periodic!");
        }

/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Type  [1 1 1] burgers vector
 */
        linedir[0][0] =  0.5773503;  
        linedir[0][1] =  0.5773503;
        linedir[0][2] =  0.5773503; 
        
        linedir[1][0] = -0.8164966;
        linedir[1][1] =  0.4082483;
        linedir[1][2] =  0.4082483; 
        
        linedir[2][0] =  0.4082483;
        linedir[2][1] = -0.8164966;
        linedir[2][2] =  0.4082483; 
        
        linedir[3][0] =  0.4082483;
        linedir[3][1] =  0.4082483; 
        linedir[3][2] = -0.8164966;
        
/*
 *      Type [-1 1 1] burgers vector
 */
        linedir[4][0] = -0.5773503;
        linedir[4][1] =  0.5773503;
        linedir[4][2] =  0.5773503; 
        
        linedir[5][0] =  0.8164966;
        linedir[5][1] =  0.4082483;
        linedir[5][2] =  0.4082483; 
        
        linedir[6][0] =  0.4082483;
        linedir[6][1] =  0.8164966;
        linedir[6][2] = -0.4082483; 
        
        linedir[7][0] =  0.4082483;
        linedir[7][1] = -0.4082483; 
        linedir[7][2] =  0.8164966;
        
/*
 *      Type [1 -1 1] burgers vector
 */
        linedir[8][0] =  0.5773503;
        linedir[8][1] = -0.5773503;
        linedir[8][2] =  0.5773503; 
        
        linedir[9][0] =  0.4082483;
        linedir[9][1] =  0.8164966;
        linedir[9][2] =  0.4082483; 
        
        linedir[10][0] =  0.8164966;
        linedir[10][1] =  0.4082483;
        linedir[10][2] = -0.4082483; 
        
        linedir[11][0] = -0.4082483;
        linedir[11][1] =  0.4082483; 
        linedir[11][2] =  0.8164966;
        
/*
 *      Type [1 1 -1] burgers vector
 */
        linedir[12][0] =  0.5773503;
        linedir[12][1] =  0.5773503;
        linedir[12][2] = -0.5773503; 
        
        linedir[13][0] = -0.4082483;
        linedir[13][1] =  0.8164966;
        linedir[13][2] =  0.4082483; 
        
        linedir[14][0] =  0.8164966;
        linedir[14][1] = -0.4082483;
        linedir[14][2] =  0.4082483; 
        
        linedir[15][0] =  0.4082483;
        linedir[15][1] =  0.4082483; 
        linedir[15][2] =  0.8164966;
        
/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.  For edges, glide plane
 *      will simply be cross product between burgers vector and 
 *      linedir.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

        minCoord[0] = param->xBoundMin;
        maxCoord[0] = param->xBoundMax;
        minCoord[1] = param->yBoundMin;
        maxCoord[1] = param->yBoundMax;
        minCoord[2] = param->zBoundMin;
        maxCoord[2] = param->zBoundMax;

        for (i = 0; i < 3; i++) {
            range[i] = maxCoord[i] - minCoord[i];
            rangeCntr[i] = 0.5 * (minCoord[i] + maxCoord[i]);
        }

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        baseNodeID = 0;
        inData->nodeCount = 0;

        for (chain = 0; chain < numChains; chain++) {
        
            numPoints = 3;
            lineIndex = chain % 16; 
            burgIndex = 4 * (lineIndex / 4);
            gpBurgIndex = lineIndex / 4;
            gpIndex = (chain / 16) % 6;
            isScrew = (chain % 4) == 0;

            normalizedBurg[0] = linedir[burgIndex][0];
            normalizedBurg[1] = linedir[burgIndex][1];
            normalizedBurg[2] = linedir[burgIndex][2];

            Normalize(&normalizedBurg[0], &normalizedBurg[1],
                      &normalizedBurg[2]);

/*
 *          Select an initial point (p0) that is within the boundaries
 *          of the simulation and then calculate the positions at which
 *          a dislocation line with the given line direction would
 *          interesect the nearest free surface in each direction.
 */
            for (i = 0; i < 3; i++) {
                p0[i] = (randm(&seed)-0.5) * (0.5 * range[i]) + rangeCntr[i];
            }

            minLen1 = 1.0e+20;
            minLen2 = 1.0e+20;

            for (i = 0; i < 3; i++) {
                if (pbc[i] == 0) {
                    signFact = (linedir[burgIndex][i] < 0.0 ? -1 : 1);
                    surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
                    len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
                    if (len < minLen1) {
                        minLen1 = len;
                        intersectIndex1 = i;
                        intersectSurfCoord1 = surfCoord;
                    }

                    signFact = -signFact;
                    surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
                    len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
                    if (len < minLen2) {
                        minLen2 = len;
                        intersectIndex2 = i;
                        intersectSurfCoord2 = surfCoord;
                    }
                }
            }

/*
 *          We know how far the dislocation can extend in each direction, now
 *          calculate the exact intersect point in both directions
 */
            for (i = 0; i < 3; i++) {

                if (i == intersectIndex1) {
                    intersectPos1[i] = intersectSurfCoord1;
                } else {
                    intersectPos1[i] = p0[i] + (minLen1 * normalizedBurg[i]);
                }

                if (i == intersectIndex2) {
                    intersectPos2[i] = intersectSurfCoord2;
                } else {
                    intersectPos2[i] = p0[i] - (minLen2 * normalizedBurg[i]);
                }
            }

/*
 *          Find a vector from the first intersection point to the second,
 *          calculate how many segments the line should be broken into based
 *          on the <maxSeg> value.
 */
            for (i = 0; i < 3; i++) {
                vector[i] = intersectPos2[i] - intersectPos1[i];
            }

            vecLen = sqrt(vector[0]*vector[0] +
                          vector[1]*vector[1] +
                          vector[2]*vector[2]);

            numSegs = (int)(vecLen / (.95 * param->maxSeg)) + 1;
            numPoints = numSegs + 1;

            for (i = 0; i < 3; i++) {
                vector[i] /= (real8)numSegs;
            }

/*
 *          Reallocate the node array with sufficient size to add
 *          all the new nodes defining this chain.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numPoints;
            inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                             * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
        

/*
 *          Starting with the first intersection point, create a
 *          series of dislocation segments ending at the second point.
 */
            newPos[0] = intersectPos1[0];
            newPos[1] = intersectPos1[1];
            newPos[2] = intersectPos1[2];

            FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);

            for (i = 0; i < numPoints; i++) {
                if (i == 0) {
                    numConnections = 1;
                    burgSign = -1.0;
                    nbr1Index = 1;
                } else if (i == (numPoints - 1)) {
                    numConnections = 1;
                    burgSign = 1.0;
                    nbr1Index = i - 1;
/*
 *                  Make sure the final node is at the surface
 *                  intersection point
 */
                    newPos[0] = intersectPos2[0];
                    newPos[1] = intersectPos2[1];
                    newPos[2] = intersectPos2[2];
                    FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);
                } else {
                    numConnections = 2;
                    burgSign = 1.0;
                    nbr1Index = i - 1;
                    nbr2Index = i + 1;
                }

                node = &inData->node[baseNodeID+i];
                node->myTag.domainID = dislocType;
                node->myTag.index    = baseNodeID+i;

                node->x = newPos[0];
                node->y = newPos[1];
                node->z = newPos[2];

                node->constraint = UNCONSTRAINED;

                AllocNodeArms(node, numConnections);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index    = baseNodeID + nbr1Index;

                node->burgX[0] = burgSign * linedir[burgIndex][0];
                node->burgY[0] = burgSign * linedir[burgIndex][1];
                node->burgZ[0] = burgSign * linedir[burgIndex][2];

/*
 *              For screw dislocations, use the glide plane from the table.
 *              For edge dislocations, glide plane is the cross product
 *              of the burgers vector and line direction.
 */
                if (isScrew) {
                    node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
                    node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
                    node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
                    if (numConnections == 2) {
                        node->nx[1] = glidePlane[gpBurgIndex][gpIndex][0];
                        node->ny[1] = glidePlane[gpBurgIndex][gpIndex][1];
                        node->nz[1] = glidePlane[gpBurgIndex][gpIndex][2];
                    }
                } else {
                    cross(linedir[burgIndex], linedir[lineIndex], plane);
                    plane[0] = (floor(plane[0] * 1.0e+07)) * 1.0e-07;
                    plane[1] = (floor(plane[1] * 1.0e+07)) * 1.0e-07;
                    plane[2] = (floor(plane[2] * 1.0e+07)) * 1.0e-07;
                    node->nx[0] = plane[0];
                    node->ny[0] = plane[1];
                    node->nz[0] = plane[2];
                    if (numConnections == 2) {
                        node->nx[1] = plane[0];
                        node->ny[1] = plane[1];
                        node->nz[1] = plane[2];
                    }
                }

/*
 *              Calculate the next node's position relative to this one.
 */
                newPos[0] += vector[0];
                newPos[1] += vector[1];
                newPos[2] += vector[2];

                FoldBox(param, &newPos[0], &newPos[1], &newPos[2]);

                if (numConnections == 1) {
                    node->constraint = SURFACE_NODE;
                    continue;
                }

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index    = baseNodeID + nbr2Index;

                node->burgX[1] = -burgSign * linedir[burgIndex][0];
                node->burgY[1] = -burgSign * linedir[burgIndex][1];
                node->burgZ[1] = -burgSign * linedir[burgIndex][2];
            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);

            lastBlock = (chain == (numChains - 1));

            if (lastBlock) {
                    IncDislocationDensity(inData, totDislocLen);
                    param->nodeCount = inData->nodeCount;
                    WriteInitialNodeData(home, inData, lastBlock);
                    FreeInNodeArray(inData, inData->nodeCount);
                    inData->nodeCount = 0;
            }

            baseNodeID = inData->nodeCount;
            startRemeshIndex = baseNodeID;

        }  /* for (chain = 0; chain < numChains; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFRSource
 *      Description:  Generate a configuration consisting of one or more
 *                    frank-read sources of random lengths between the
 *                    specified minimum and maximum lengths.  
 *
 *                    NOTE:  This function currently assumes periodic
 *                    boundary conditions!
 *
 *      Arguments:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numSources    Number of frank read sources to create
 *          srcLenMin     Minimum length of any created frank-read source
 *          srcLenMax     Maximum length of any created frank-read source
 *          seed          Seed value for random number generator
 *          totDislocLen  Pointer to location at which to return
 *                        to caller the total combined length of
 *                        all created dislocations.
 *
 *-------------------------------------------------------------------------*/

#define INTERACTION_TEST 0

void GetFreePath(Param_t *param, real8 pos1[3], real8 pos2[3], real8 n[0], Grain_t *grain, 
                 real8 *interDir, real8 &lambda_SAS1, real8 &lambda_GB1, real8 &lambda_SAS2, real8 &lambda_GB2)
{
    int         i, j, k, m, th, type; 
    real8       interPos[3], poss[2][3], vec1[3], vec2[3], vec3[3], dis;
    Plane_t     plane;

    lambda_GB1 = 0.0;
    lambda_GB2 = 0.0;
    lambda_SAS1 = 1E10;
    lambda_SAS2 = 1E10;

    VECTOR_COPY(plane.point, pos1);
    VECTOR_COPY(plane.normal, n);
    VECTOR_ZERO(interDir);

    for(k=0; k<grain->surface.size(); k++){
        type = PlanePlaneIntersection(
                  plane.normal, plane.point,
                  grain->surface[k].normal, grain->surface[k].point, vec1, vec2); 
        
        if(type != LINE_INTERSECTION)continue;

        th = 0;
        for(m=0; m<grain->surface.size(); m++){
            if(m==k)continue;
            type = LinePlaneIntersection(vec1, vec2, grain->surface[m].normal,
                         grain->surface[m].point, vec3, interPos);

            if(type != POINT_INTERSECTION)continue;

            if(fabs(interPos[0]) > param->xBoundMax+1.0){
                continue;
            }
            if(fabs(interPos[1]) > param->yBoundMax+1.0){
                continue;
            }
            if(fabs(interPos[2]) > param->zBoundMax+1.0){
                continue;
            }

            VECTOR_COPY(poss[th], interPos);

            if(th == 1){
                th = 0;
                break;
            }else{
                th++;
            }
        }

        if(k == 0){
            VECTOR_COPY(interDir, vec1);
            PointLineIntersection(pos1, vec1, vec2, vec3, &dis);
            lambda_GB1 = fabs(dis);
            PointLineIntersection(pos2, vec1, vec2, vec3, &dis);
            lambda_GB2 = fabs(dis);
        }else{
            PointLineIntersection(pos1, vec1, vec2, vec3, &dis);
            lambda_SAS1 = MIN(fabs(dis), lambda_SAS1);
            PointLineIntersection(pos2, vec1, vec2, vec3, &dis);
            lambda_SAS2 = MIN(fabs(dis), lambda_SAS2);
        }
    
    }

    return;
}

void CreateFRSource(Home_t *home, InData_t *inData, InArgs_t *inArgs,
						 real8 *totDislocLen)
{
        int      i, lineIndex, chain, baseNodeID, numPoints=3,ranNum;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      x0, d, nLayers, xn, layer, time;
        int      startRemeshIndex = 0;
        real8    srcLen, thickness,lenRange;
        real8    p0[3], p1[3], p2[3], p3[3], randomLine[3];
        real8	 randomAngle, edge, srcLenMin, srcLenMax;
        int	     cubeLength, numSources,  seed, dislocType;
        int	     polycrystal, nGrain, finish;
        int      firstBurgIndex, srcLenType;
        int      inGrain1, onGB1, inGrain2, onGB2, inGrain0, onGB0;
        int      inGrain3, onGB3;
        real8	 twinRatio, maxSeg, maxDensity, burgMag = 2.556E-10;
        real8	 volume, glideDir[3];
        real8	 tnx[4], tny[4], tnz[4], burg[12][3], cubeSize,heightSize;
        real8	 n0[3], b[3], n1[3];	/* plane normal vector, burgers, normal vector*/
        Node_t   *node;
        Param_t  *param;
        real8	 invsq2, basePos[3], b1[3], b2[3], f1[3], f2[3], f3[3], f4[3];
        real8	 rotMatrix[3][3],rotLine[3],rotBurg[3],rotPlane[3], btemp[3];
        int      intID, burgID, planeID;
        real8    pi = 3.1415926535897932385, totSrcLen;
        real8    separateLength = 40.0;
        bool        finished, singleSrcTest = 0;
        BurgInfo_t  *bd, *gbd;
        GB_t        *gb;
        Grain_t     *grain;
        	
        real8                                       str[3][3], gbNormal[3];
        int                                         appliedDir, j, k;
        real8                                       sign;
        std::vector<std::vector<Dislocation_t> >    disArray;
        std::vector<std::vector<real8> >            burgList;
        std::vector<std::vector<Dislocation_t> >::iterator      itDisList;

        real8   normalABC[3], lambda_FR, lambda_SAS1, lambda_GB1, lambda_SAS2, lambda_GB2, interDir[3];
        real8   barrierFactor, crossProd[3], vec[3];
        real8   vec1[3], vec2[3], s0,  paras[4][4], zaxis[3] = {0,0,1}, schmid;
        int l; real8 dis;

        normalABC[0] = 1;
        normalABC[1] = 1;
        normalABC[2] = 1;

#ifdef SAS_ANALYZE
        FILE    *fp;
        
        if((fp = fopen(inArgs->outFile, "w")) == (FILE *)NULL){
            Fatal("can not open file %s", inArgs->outFile);
        }
#endif

/*
 *      Grap the command-line argument values needed for this function
 */
        home->param->materialType = MAT_TYPE_FCC;
        	
        cubeLength  = inArgs->cubeLength;
	    srcLenType  = inArgs->srcLenType;
        nLayers     = inArgs->nLayers;
       	maxDensity  = inArgs->maxDensity;
        numSources  = inArgs->numFRSrcs;
        srcLenMin   = inArgs->frLenMin;
        srcLenMax   = inArgs->frLenMax;
        seed        = inArgs->seed;
        dislocType  = inArgs->type;
        polycrystal = inArgs->polycrystal;
        nGrain	    = inArgs->nGrain;
        twinRatio   = inArgs->twinRatio;
        maxSeg      = (real8)inArgs->maxSegLen;
        burgID      = inArgs->burgID;
        planeID     = inArgs->planeID;
        	
        param = inData->param;
        d = inArgs->deltaNums;
        cubeSize	=	(real8)cubeLength;
        if(srcLenMax<0) {srcLenMax = -srcLenMax; srcLenType = 2;} 
        if(srcLenMin<0) {srcLenMin = -srcLenMin; srcLenType = 2;}
        lenRange = srcLenMax - srcLenMin;

        bool freeLine=(param->xBoundType != Periodic && param->yBoundType != Periodic && param->zBoundType != Periodic);

        real8 n3[3];

        if(numSources <= 0){
            singleSrcTest = 1;
        }
        	
        if(polycrystal){
            gb = home->GBKeys[0];
            gbd = &(gb->burgData);
        }
        n1[0] = 0.0;	n1[1] = 0.0;	n1[2] = 0.0;
        real8 maximalLen, a1, a2, randmVal;
        bool partialDislocation;

/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        baseNodeID = 0;
        inData->nodeCount = 0;
		
        heightSize = 2 * param->zBoundMax;
	    volume = heightSize * cubeSize * cubeSize;
        thickness = heightSize / ((real8)nLayers);
        	
        x0 = (d*nLayers*(1 - nLayers) + 2*numSources)/(2 * nLayers);
	    totSrcLen = 0.0;
        real8 primarySrcLen = 0.0;
	    finish = 0;

        int nSlips;
        if(param->polycrystal>0){ if(inArgs->nTDs>0)nSlips=N_SS; else nSlips=N_LD; }else nSlips=N_SG;
        real8 *proba, *primary, *Proba, *percen, maxSchmid=0;
        proba = (real8 *)malloc(nSlips*sizeof(real8));
        Proba = (real8 *)malloc(nSlips*sizeof(real8));
        percen = (real8 *)malloc(nSlips*sizeof(real8));
        primary = (real8 *)malloc(nSlips*sizeof(real8));

        for(i=0; i<nSlips; i++){
            j = slip_sys[i][2]>-1?slip_sys[i][2]:0;
            bd=GetBurgData(home, j, -1);
            proba[i] = fabs(DotProduct(inArgs->loadDir, bd->burgList[slip_sys[i][0]])*
                            DotProduct(inArgs->loadDir,bd->planeList[slip_sys[i][1]])/Normal(bd->burgList[slip_sys[i][0]]));
            if(proba[i]>maxSchmid){
                maxSchmid = proba[i];
            }
            if(home->nGBs>1)Proba[i]=proba[i]; else Proba[i]=1.0;
            percen[i]=0.0;
        }
        for(i=0;i<nSlips;i++)if(fabs(proba[i]-maxSchmid)<1E-4)primary[i]=1.0; else primary[i]=MAX(proba[i]/maxSchmid, 0.33333);
        for(i=1;i<nSlips;i++)Proba[i]+=Proba[i-1];
        for(i=0;i<nSlips;i++)Proba[i]/=Proba[nSlips-1];
        
        for (layer = 0; layer < nLayers; layer++) {

	        if(finish)break;

            xn = layer * d + x0;
            if(xn <= 0){
                xn = 0;
                continue;
            }
            	
            for(chain = 0; chain < xn; ){
             
                if(xn > nSlips && chain < nSlips && (2*maxDensity*volume*burgMag*burgMag/srcLenMax > nSlips)){
                    ranNum = chain%nSlips;
                }else{
                    randmVal= randm(&seed);
                    for(i=0;i<nSlips;i++)if(randmVal<Proba[i]&& (randmVal>((i==0)?0:Proba[i-1])) )ranNum=i;
                }
                gpBurgIndex = slip_sys[ranNum][0]; gpIndex=slip_sys[ranNum][1];
    
                finished = 0;
                while(!finished){
#if 1
                    p1[0] = (randm(&seed)-0.5) * cubeSize;
                    p1[1] = (randm(&seed)-0.5) * cubeSize;			
                    p1[2] = randm(&seed)* thickness + ((real8)layer) * thickness - 0.5*heightSize;

#else
                    VECTOR_ZERO(p0);
                    p0[1] = -1000;
                    p0[2] = -200;
#endif

                    NodeLocation(home, p1[0],  p1[1], p1[2], inGrain1, onGB1,15);
                    bd = GetBurgData(home, inGrain1, onGB1);
                    if(onGB1 > -1)continue;

#if 1
                    if(gpBurgIndex>=6){
                        ClosedGB(home, p1, &gb);
                        bd = &(gb->burgData);
                        PointPlaneDistance(p1, gb->n, gb->p, dis, p1);
                        inGrain1 = -1; onGB1 = gb->GBID;
                    }else if(gpIndex!=3){
                        if(inArgs->polycrystal ==1 && (inGrain1%2)!=slip_sys[ranNum][2]) continue;
                    }
#else
                    gpIndex = 1;
                    gpBurgIndex = 5;
//                  p1[2] =  home->grainKeys[inGrain1]->cz;
                    VECTOR_ZERO(p0);
#endif   
                    VECTOR_COPY(n0, bd->planeList[gpIndex]);
                    VECTOR_COPY(b , bd->burgList[gpBurgIndex]);
                    VECTOR_COPY(n3, bd->planeList[3]);
                    if(!PerpendicularDirection(n0,b))Fatal("burg is not perpendicular to glide plane");

                    time = 0;
                    while(true){
                        if(time>100)break;
                        time++;
                        randomAngle	= acos((randm(&seed) - 0.5) * 2.0);
                        randomAngle = ((randm(&seed) > 0.5) ? randomAngle+pi : randomAngle);
                        
                        cross(n0, b, n1);
                        NormalizeVec(n1);
                	    
                        for(i = 0; i < 3; i++){
                            randomLine[i] = n1[i] * sin(randomAngle) + b[i] * cos(randomAngle);
                        }
                        if(polycrystal>0 && gpIndex != 3){
                            edge = home->grainKeys[inGrain1]->thick/0.942809;
                        }else edge = cubeLength;
                        NormalizeVec(randomLine);
#if 0
                        VECTOR_ZERO(p1); VECTOR_ZERO(randomLine); randomLine[0]=1;
#endif
/*           
 *                     Length of the frank read source should be a random length
 *                     between srcLenMin and srcLenMax.
 */          
	                    switch(srcLenType){
	                        case 0:
                                if (lenRange > 0) {
                                    srcLen = (srcLenMin + (randm(&seed) * lenRange));
                                    if(srcLen > srcLenMax){
                                       srcLen = srcLenMax;
                                    }
                                } else {
                                    srcLen = srcLenMin;
                                }
	 	                        break;
                        
	 	                    case 1:
	 	                        srcLen = -log(randm(&seed))/sqrt(maxDensity)/burgMag;
                                printf("chain=%d length=%f\n",chain,srcLen);
	 	                        break;
                            case 2:
                                if (lenRange > 0) {
                                    srcLen = (srcLenMin + (randm(&seed) * lenRange))*edge;
                                    if(srcLen > srcLenMax*edge){
                                        srcLen = srcLenMax*edge;
                                    }
                                } else {
                                    srcLen = srcLenMin*edge;
                                }
	 	                        break;
                            default:
                                break;
	                    }

                        if(home->nGBs>1)srcLen *= primary[ranNum];

                        if(home->nGBs<=1 || gpIndex==3){
                            for (i = 0; i < 3; i++) {
                                p0[i] = p1[i] + (srcLen * randomLine[i])/2.0;
                                p2[i] = p1[i] + (srcLen * randomLine[i]);
                            } 
                            NodeLocation(home, p0[0],  p0[1], p0[2], inGrain0, onGB0,15);
                            NodeLocation(home, p2[0],  p2[1], p2[2], inGrain2, onGB2,15);
                            if(inGrain1 != inGrain0 || inGrain1 != inGrain2 || onGB1 != onGB0 || onGB1 != onGB2)continue;
                        }else{
                            
                            VECTOR_COPY(p0, p1); inGrain0=inGrain1; onGB0 = onGB1;
                            for (i = 0; i < 3; i++) {
                                p1[i] = p0[i] - (srcLen * randomLine[i])/2.0;
                                p2[i] = p0[i] + (srcLen * randomLine[i])/2.0;
                            } 
                            NodeLocation(home, p1[0],  p1[1], p1[2], inGrain1, onGB1, 15);
                            NodeLocation(home, p2[0],  p2[1], p2[2], inGrain2, onGB2, 15);
                            grain = home->grainKeys[inGrain0];

                            if(!(onGB1<0 && inGrain1==inGrain0) && !(onGB2<0 && inGrain2==inGrain0)){
                                a1 = atan(grain->thick/0.942809/(srcLen));
                                cross(param->TBNormal, n0, vec1); NormalizeVec(vec1);
                                cross(vec1, n0, vec2); NormalizeVec(vec2);
                                CoordAdd(((randm(&seed)>0.5)?-1.0:1.0)*cos(a1),vec1,((randm(&seed)>0.5)?-1.0:1.0)*sin(a1),vec2, randomLine);
                                NormalizeVec(randomLine);
                
                                for(i=0;i<grain->nBounds;i++){
                                    gb = grain->GB[i];
                                    gbNormal[0] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[0];
                                    gbNormal[1] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[1];
                                    gbNormal[2] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[2];
                                
                                    if(DotProduct(gbNormal,randomLine)<0){
                                        VECTOR_COPY(vec2, gb->p);
                                        PBCPOSITION(param, p0[0], p0[1], p0[2], vec2, vec2+1, vec2+2);
                                        LinePlaneIntersection(randomLine,p0,gbNormal,vec2,vec1,p1);
                                        inGrain1=-1; onGB1 = gb->GBID;
                                        break;
                                    }
                                }
                                if(i==grain->nBounds)Fatal("can not find gb in %s, %d",__FILE__,__LINE__);

                                for(i=0;i<grain->nBounds;i++){
                                    gb = grain->GB[i];
                                    gbNormal[0] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[0];
                                    gbNormal[1] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[1];
                                    gbNormal[2] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[2];
                                
                                    if(DotProduct(gbNormal,randomLine)>0){
                                        VECTOR_COPY(vec2, gb->p);
                                        PBCPOSITION(param, p0[0], p0[1], p0[2], vec2, vec2+1, vec2+2);
                                        LinePlaneIntersection(randomLine,p0,gbNormal,vec2,vec1,p2);
                                        inGrain2=-1; onGB2 = gb->GBID;
                                        break;
                                    }
                                }
                                if(i==grain->nBounds)Fatal("can not find gb in %s, %d",__FILE__,__LINE__);
                            }else if(!(onGB1<0 && inGrain1==inGrain0)){
                                for(i=0;i<grain->nBounds;i++){
                                    gb = grain->GB[i];
                                    gbNormal[0] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[0];
                                    gbNormal[1] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[1];
                                    gbNormal[2] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[2];
                                
                                    if(DotProduct(gbNormal,randomLine)<0){
                                        VECTOR_COPY(vec2, gb->p);
                                        PBCPOSITION(param, p0[0], p0[1], p0[2], vec2, vec2+1, vec2+2);
                                        LinePlaneIntersection(randomLine,p0,gbNormal,vec2,vec1,p1);
                                        inGrain1=-1; onGB1 = gb->GBID;
                                        break;
                                    }
                                }
                                if(i==grain->nBounds)Fatal("can not find gb in %s, %d",__FILE__,__LINE__);
                            }else if(!(onGB2<0 && inGrain2==inGrain0)){
                                for(i=0;i<grain->nBounds;i++){
                                    gb = grain->GB[i];
                                    gbNormal[0] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[0];
                                    gbNormal[1] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[1];
                                    gbNormal[2] = (grain->type==FOLDED_GRAIN?-1:1)*grain->sign[i]*gb->n[2];
                                
                                    if(DotProduct(gbNormal,randomLine)>0){
                                        VECTOR_COPY(vec2, gb->p);
                                        PBCPOSITION(param, p0[0], p0[1], p0[2], vec2, vec2+1, vec2+2);
                                        LinePlaneIntersection(randomLine,p0,gbNormal,vec2,vec1,p2);
                                        inGrain2=-1; onGB2 = gb->GBID;
                                        break;
                                    }
                                }
                                if(i==grain->nBounds)Fatal("can not find gb in %s, %d",__FILE__,__LINE__);
                            }

                            srcLen = PointPointDistance(p1,p2);
                        }

                        if(inArgs->extended != 0){
                            VECTOR_COPY(basePos, p0);
                            cross(n0, randomLine, glideDir);
                            NormalizeVec(glideDir);
                            p0[0] -= (separateLength*glideDir[0]*0.5); // drag partial dislocation
                            p0[1] -= (separateLength*glideDir[1]*0.5);
                            p0[2] -= (separateLength*glideDir[2]*0.5);

                            p3[0] = p0[0] + separateLength*glideDir[0]; // leading partial dislocation
                            p3[1] = p0[1] + separateLength*glideDir[1];
                            p3[2] = p0[2] + separateLength*glideDir[2];
                            NodeLocation(home, p0[0],  p0[1], p0[2], inGrain0, onGB0,15);
                            NodeLocation(home, p3[0],  p3[1], p3[2], inGrain3, onGB3,15);
                    
                            if(onGB0 > -1 || onGB3 > -1 || inGrain0 != inGrain2 || inGrain3 != inGrain2)continue;
                        }   

                        if(srcLenType == 0 || inArgs->gb==true){
                            if(!InCuboid(param, p0, 0) || !InCuboid(param, p1, 0) || !InCuboid(param, p2, 0)){
                                continue;
                            }
	                    }

                        break;
                    }
                    if(time>100)continue;
                    break;
                }

                if(inArgs->extended != 0){
                    if(!PartialDislocations(bd, b, n0, b1, b2)){
                        partialDislocation = true;
                    }else partialDislocation = false;
                }else partialDislocation= false;

                if(primary[ranNum]>0.99) primarySrcLen += srcLen;
                percen[ranNum] += srcLen;
	            totSrcLen += srcLen;				
#ifdef SAS_ANALYZE
                if(polycrystal){
                    lambda_FR = srcLen/2;
                    GetFreePath(param, p1, p2, n0, home->grainKeys[inGrain0], 
                                interDir, lambda_SAS1, lambda_GB1, lambda_SAS2, lambda_GB2);

                    ExtPKForce(str, b[0], b[1], b[2], 
                               p1[0], p1[1], p1[2],
                               p2[0], p2[1], p2[2], f1, f2);
                   
                    real8 dis; 
                    PointPlaneIntersection(p0, gb->n, gb->p, vec, &dis);
                    vec[0] -= p0[0];
                    vec[1] -= p0[1];
                    vec[2] -= p0[2];
                    if(DotProduct(vec, f1) < 0){
                        lambda_GB1 = -lambda_GB1;
                        lambda_GB2 = -lambda_GB2;
                    }

                    s0 = DotProduct(zaxis, n0)*DotProduct(zaxis, b);
                    VECTOR_ZERO(paras[0]);
                    VECTOR_ZERO(paras[1]);
                    VECTOR_ZERO(paras[2]);
                    VECTOR_ZERO(paras[3]);

                    real8   schmidList[2] = {0,0}, maxSchmid = -1;
                    int     maxSet = 0;
                
                    if(lambda_GB1 != 0.0 && fabs(s0) > 1.0E-2){
                        std::vector<std::vector< Dislocation_t> >().swap(disArray);
                        std::vector<std::vector<double> >().swap(burgList);

                        if(!GeometricCriterion(gbd, b, interDir, n0, 1,
                            burgList, disArray)){
                            Fatal("something wrong with geometric criterion");
                        }

/*
 *                      Now, we don't need the resultant configure including a lattice burg on GB,
 *                      or a partial burg in Grain.
 */
                        for(j=0; j<disArray.size(); ){
                            barrierFactor = -1.0;
                            for(k=0; k<disArray[j].size(); k++){
                                barrierFactor += disArray[j][k].energyFactor;
                                cross(gb->n, gbd->planeList[disArray[j][k].planeID], crossProd);
                                VECTOR_COPY(vec, burgList[disArray[j][k].burgID]);
                                if(Normal(crossProd) < 1.0E-2){
                                    if(BurgType_FCC(vec) == BASIC_BURG)break;
                                }else{
                                    if(BurgType_FCC(vec) != BASIC_BURG)break;
                                }
                            }
  
                            if(k != disArray[j].size() || j > 0 && 
                               barrierFactor > 0){
                                itDisList = disArray.begin() + j;
                                disArray.erase(itDisList);
                            }else{
                                j++;
                            }
                        }
                        if(disArray.size()>2)Fatal("someting wrong with disArray");
                    
                        for(j=0; j<disArray.size(); j++){
                            if(disArray[j].size()>2)Fatal("someting wrong with disArray[%d]", j);     

                            for(k=0; k<disArray[j].size(); k++){
                                paras[2*j+k][0] = disArray[j][k].energyFactor;
                                VECTOR_COPY(vec, burgList[disArray[j][k].burgID]);
                                paras[2*j+k][1] = Normal(vec);
                                NormalizeVec(vec);
                                schmid = DotProduct(zaxis, vec)*DotProduct(zaxis, gbd->planeList[disArray[j][k].planeID]);

	                            if(fabs(paras[2*j+k][1] - 1.0E0) < 1.0E-2){
                                    paras[2*j+k][2] = (schmid*s0 > 0) ? schmid : 0;
	                            }else{
                                    paras[2*j+k][2] = schmid;
                                }
                                paras[2*j+k][2] = fabs(paras[2*j+k][2]);
                                
                            }
                            schmidList[j] = DotProduct(zaxis,b)*DotProduct(zaxis, gbd->planeList[disArray[j][0].planeID]);
                        }
                    }
                    
                    for(j=0; j<2; j++){
                        if(fabs(schmidList[j]) > maxSchmid){
                            maxSet = j;
                            maxSchmid = fabs(schmidList[j]);
                        }
                    }

                    fprintf(fp,"%e %e %e %e   ", lambda_FR, lambda_SAS1, lambda_GB1, fabs(s0));
                    for(j=2*maxSet; j<2*maxSet+2; j++){
                        fprintf(fp, "%e %e %e ", paras[j][0], paras[j][1], paras[j][2]);
                    }
                    fprintf(fp, "\n");

                    fprintf(fp,"%e %e %e %e ", lambda_FR, lambda_SAS2, lambda_GB2, fabs(s0));
                    for(j=2*maxSet; j<2*maxSet+2; j++){
                        fprintf(fp, "%e %e %e ", paras[j][0], paras[j][1], paras[j][2]);
                    }
                    fprintf(fp, "\n");

                }
#endif 
/*           
 *              Reallocate the node array with sufficient size to add
 *              all the new nodes defining this chain.
 */          
                numPoints = (inArgs->extended == 0 || partialDislocation) ? 3 : 4;
                newNodeIndex = inData->nodeCount;
                inData->nodeCount += numPoints;
                inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                                   * sizeof(Node_t));
                memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
		 	   
                if(inArgs->extended == 0 || partialDislocation){
/*           
 *                   Set up the 3 nodes we're using to define the dislocation line.
 *                   First do point p0.
 */               
                     node = &inData->node[baseNodeID];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID;
                     node->x = p0[0];
                     node->y = p0[1];
                     node->z = p0[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain0;
                     node->onGB = onGB0;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = -b[0];
                     node->burgY[0] = -b[1];
                     node->burgZ[0] = -b[2];
                  
                     node->burgX[1] = b[0];
                     node->burgY[1] = b[1];
                     node->burgZ[1] = b[2];
                  
/*                
 *                   For screw dislocations, use the glide plane from the table.
 *                   For edge dislocations, glide plane is the cross product
 *                   of the burgers vector and line direction.
 */               
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1;
                  
                     node->x = p1[0];
                     node->y = p1[1];
                     node->z = p1[2];
                  
                     node->constraint = PINNED_NODE;
                     node->inGrain = inGrain1;
                     node->onGB = onGB1;
                     AllocNodeArms(node, 1);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = b[0];
                     node->burgY[0] = b[1];
                     node->burgZ[0] = b[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                  
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2;
                  
                     node->x = p2[0];
                     node->y = p2[1];
                     node->z = p2[2];
                  
                     node->constraint = PINNED_NODE;
                     node->inGrain = inGrain2;
                     node->onGB = onGB2;
                     AllocNodeArms(node, 1);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = -b[0];
                     node->burgY[0] = -b[1];
                     node->burgZ[0] = -b[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                }else{
                     node = &inData->node[baseNodeID];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID;
                     node->x = p0[0];
                     node->y = p0[1];
                     node->z = p0[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain0;
                     node->onGB = onGB0;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = -b1[0];
                     node->burgY[0] = -b1[1];
                     node->burgZ[0] = -b1[2];
                     node->signSF[0] = -1.0;
                  
                     node->burgX[1] = b1[0];
                     node->burgY[1] = b1[1];
                     node->burgZ[1] = b1[2];
                     node->signSF[1] = 1.0;
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1;
                  
                     node->x = p1[0];
                     node->y = p1[1];
                     node->z = p1[2];
                  
                     node->constraint = PINNED_NODE;
                     node->inGrain = inGrain1;
                     node->onGB = onGB1;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+3;
                  
                     node->burgX[0] = b1[0];
                     node->burgY[0] = b1[1];
                     node->burgZ[0] = b1[2];
                     node->signSF[0] = 1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                  
                     node->burgX[1] = b2[0];
                     node->burgY[1] = b2[1];
                     node->burgZ[1] = b2[2];
                     node->signSF[1] = -1.0;
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2;
                  
                     node->x = p2[0];
                     node->y = p2[1];
                     node->z = p2[2];
                  
                     node->constraint = PINNED_NODE;
                     node->inGrain = inGrain2;
                     node->onGB = onGB2;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+3;
                  
                     node->burgX[0] = -b1[0];
                     node->burgY[0] = -b1[1];
                     node->burgZ[0] = -b1[2];
                     node->signSF[0] = -1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];

                     node->burgX[1] = -b2[0];
                     node->burgY[1] = -b2[1];
                     node->burgZ[1] = -b2[2];
                     node->signSF[1] = 1.0;
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
/*                
 *                   Now point p3...
 */               
                     node = &inData->node[baseNodeID+3];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+3;
                     node->x = p3[0];
                     node->y = p3[1];
                     node->z = p3[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain3;
                     node->onGB = onGB3;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = -b2[0];
                     node->burgY[0] = -b2[1];
                     node->burgZ[0] = -b2[2];
                     node->signSF[0] = 1.0;
                  
                     node->burgX[1] = b2[0];
                     node->burgY[1] = b2[1];
                     node->burgZ[1] = b2[2];
                     node->signSF[1] = -1.0;
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];

                }
/*           
 *              The initial segments created are not necessarily limited to
 *              param->maxSegLen, so a call to InitRemesh() is needed to
 *              chop any excessively long segments into proper lengths.
 *              When we've generated the nodal data for the final chain,
 *              write the block of nodal data to the file.
 */	         
                InitRemesh(inData, dislocType, startRemeshIndex);
//                lastBlock = (chain == (numSources - 1));
             
               
                chain++;
                baseNodeID = inData->nodeCount;
                startRemeshIndex = baseNodeID;

	            if(totSrcLen> maxDensity*volume*burgMag*burgMag){
		            finish = 1;
		            break;
	            }else{
		            finish = 0;
	            }
            }  /* for (chain = 0; chain < numChains; ) */
            printf("layer %d: xn=%d, chain=%d\n",layer,xn,chain+1);
	    }/* for(layer=0; layer<nLayers;...)*/
        printf("Slip systems | location  | burgID | planeID  | Schmid Factor | percentage\n"); 
        for(i=0; i<nSlips; i++){
            printf("    %3d      |   %3s     |   %3d  |   %3d    |    %5.3f      |  %5.2f  \n",
                i+1, (slip_sys[i][0]>=6)?"TB":((slip_sys[i][1]==3)?"M/T":(slip_sys[i][2]==0?"M":"T")), slip_sys[i][0], slip_sys[i][1], proba[i], percen[i]*100/totSrcLen);
        }

        free(proba);
        free(Proba);
        free(percen);
        free(primary);
        printf("%d chains has been distributed in matrix/twin\n"
               "Averaged number of chains over grain is %d\n" 
               "%4.2f%% of srcs are of primary slip systems\n",
                chain, chain/home->nGrains, primarySrcLen*100/totSrcLen);

#ifdef SAS_ANALYZE
        fclose(fp);
#endif
        
	    IncDislocationDensity(inData, totDislocLen);
	    param->nodeCount = inData->nodeCount;
	    WriteInitialNodeData(home, inData, 1);
	    FreeInNodeArray(inData, inData->nodeCount);
	    inData->nodeCount = 0;
	
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateEdges
 *      Description:  Creates "edge" dislocations in the [100] [010] and [001]
 *                    directions.  Each of the three line senses has 4
 *                    combinations of burgers vector and normals.  With the
 *                    opposite sign of the line sense vectors considered,
 *                    the total number of types of dislocations is 24,
 *                    so number of chains specified should be a multiple
 *                    of 24 to include all types of.  (Okay, they're
 *                    not pure edge, but they are not screw either)
 *
 *      Arguments:
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *
 *-------------------------------------------------------------------------*/
void CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
                 int numChains, int seed, real8 *totDislocLen,
                 int dislocType)
{
        int           ic, ip, np, id0, lastBlock, signFactor;
        int           newNodeIndex, startRemeshIndex;
        int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
        real8         posFactor, cubeSize;
        real8         xp[3], yp[3], zp[3];
        Param_t       *param;
        Node_t        *node;
        static real8  lineDir[3][3] = {
                      {0.0, 0.0, 1.0},
                      {0.0, 1.0, 0.0},
                      {1.0, 0.0, 0.0}};
        static real8  burg[4][3] = {
                      { 0.5773503,  0.5773503,  0.5773503},
                      { 0.5773503,  0.5773503, -0.5773503},
                      { 0.5773503, -0.5773503,  0.5773503},
                      {-0.5773503,  0.5773503,  0.5773503}};
        static real8  glidePlane[12][3] = {
                      {  0.7071068, -0.7071068,  0},  /* ldir [001] b [111]  */
                      {  0.7071068, -0.7071068,  0},  /* ldir [001] b [11-1] */
                      {  0.7071068,  0.7071068,  0},  /* ldir [001] b [1-11] */
                      {  0.7071068,  0.7071068,  0},  /* ldir [001] b [-111] */
                      {  0.7071068,  0, -0.7071068},  /* ldir [010] b [111]  */
                      {  0.7071068,  0,  0.7071068},  /* ldir [010] b [11-1] */
                      {  0.7071068,  0, -0.7071068},  /* ldir [010] b [1-11] */
                      {  0.7071068,  0,  0.7071068},  /* ldir [010] b [-111] */
                      {  0, -0.7071068,  0.7071068},  /* ldir [100] b [111]  */
                      {  0,  0.7071068,  0.7071068},  /* ldir [100] b [11-1] */
                      {  0,  0.7071068,  0.7071068},  /* ldir [100] b [1-11] */
                      {  0, -0.7071068,  0.7071068}}; /* ldir [100] b [-111] */


        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateEdges", numChains);
        }

        param = inData->param;
        cubeSize = (real8)cubeLength;
        id0 = 0;
        inData->nodeCount = 0;
        startRemeshIndex = 0;
        posFactor = 0.333 * cubeSize;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            gpIndex = ic % 12; 
            burgIndex = gpIndex % 4;
            ldIndex = gpIndex / 4;
            
/*
 *          First 12 burgers vector/normal sets use positive line
 *          sense, next set of 12 uses opposite line sense, and
 *          so on.
 */
            signFactor = ((ic / 12) & 0x01) ? -1 : 1;

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;

            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][0]);
            yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][1]);
            zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][2]);

            xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][0]);
            yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][1]);
            zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][2]);

/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                node->constraint = UNCONSTRAINED;
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+id0;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = id0 + nbr1Index;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[gpIndex][0];
                node->ny[0] = glidePlane[gpIndex][1];
                node->nz[0] = glidePlane[gpIndex][2];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = id0 + nbr2Index;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[gpIndex][0];
                node->ny[1] = glidePlane[gpIndex][1];
                node->nz[1] = glidePlane[gpIndex][2];

            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the current block of nodal data to the file.
 */
            InitRemesh(inData, dislocType, startRemeshIndex);
			
            lastBlock = (ic == numChains - 1);
            if (lastBlock) {
                IncDislocationDensity(inData, totDislocLen);
                param->nodeCount = inData->nodeCount;
                WriteInitialNodeData(home, inData, lastBlock);
                FreeInNodeArray(inData, inData->nodeCount);
                inData->nodeCount = 0;
            }

            id0 = inData->nodeCount;
            startRemeshIndex = id0;
        }
        return;
}

void CreateLoopSource(Home_t *home, InData_t *inData, InArgs_t *inArgs, real8 *totDislocLen)
{
    int     nLoops = inArgs->numLoops;
    int     seed = inArgs->seed;
    real8   cubeLength = real8(inArgs->cubeLength);
    const int     radius = inArgs->radius<=sqrt(3)/2?inArgs->radius*cubeLength:inArgs->radius;

    real8   cenPos[3], n[3], b[3], r, rRes, g[3], gRot[3], rot[3][3], segDir[3];
    int     inGrain0=-1, onGB0=-1, burgID, planeID;
    BurgInfo_t *bd;
    int     ranNum;

    GB_t    *gb;
    int baseNodeID=0, nSegs,nodeID, nbrID1, nbrID2;
    inData->nodeCount = 0;
    Param_t *param = inData->param;
    real8   maxSeg = param->maxSeg;
    real8   theta;
    Node_t  *node;
    
    home->param->materialType = MAT_TYPE_FCC;

    for(int i=0; i<nLoops;i++){
        bool finished=false; 
        rRes = radius;
        int j=-1;
        while(finished==false ){
            
            cenPos[0] = (randm(&seed)-0.5)*cubeLength;
            cenPos[1] = (randm(&seed)-0.5)*cubeLength;
            cenPos[2] = (randm(&seed)-0.5)*cubeLength;
            
            if(home->nGBs>0){
                NodeLocation(home, cenPos[0], cenPos[1], cenPos[2], inGrain0, onGB0,15);
                if(onGB0 > -1) continue;
                ranNum = (i<N_SS && nLoops>N_SS)?i:int(1.0E5*randm(&seed))%N_SS;
                if(slip_sys[ranNum][2]==0 || slip_sys[ranNum][2]==1){
                    if((inGrain0%2)!=slip_sys[ranNum][2])continue;
                    bd = GetBurgData(home, inGrain0, onGB0);
                    planeID = slip_sys[ranNum][1];
                    if(planeID != 3){
                        r = (ClosedGB(home, cenPos, &gb)-1.0)*3.0/sqrt(8.0);
                        if(r*sqrt(2)<maxSeg/2.0)continue;
                        rRes -=r; if(rRes<0)finished=true;
                    }else{
                        r = radius; if(r<maxSeg/sqrt(2.0))Fatal("radius %g < %g", r, maxSeg/sqrt(2));
                        finished = true;
                    }
                }else if (slip_sys[ranNum][2] == 3){
                    bd = GetBurgData(home, 0, -1);
                    r = radius; if(r<maxSeg/sqrt(2.0))Fatal("radius %g < %g", r, maxSeg/sqrt(2));
                    finished = true;
                }else{
                    real8 dis; j++;
                    bd = GetBurgData(home, -1, j);
                    gb = home->GBKeys[j];
                    PointPlaneDistance(cenPos, gb->n, gb->p, dis, cenPos);
                    inGrain0=-1; onGB0 = j;
                    r = radius; if(r<maxSeg/sqrt(2.0))Fatal("radius %g < %g", r, maxSeg/sqrt(2));
                    if(j>=home->nGBs-1)finished=true;
                }
            }else{
                ranNum = (i<12&&nLoops>12)?i:int(1.0E5*randm(&seed))%12; 
                bd = &(home->burgData);
                r = radius; if(r<maxSeg/sqrt(2.0))Fatal("radius %g < %g", r, maxSeg/sqrt(2));
                finished = true;
            }
            burgID = slip_sys[ranNum][0]; planeID = slip_sys[ranNum][1];
            VECTOR_COPY(b, bd->burgList[burgID]);
            VECTOR_COPY(n, bd->planeList[planeID]);
            if(!PerpendicularDirection(n,b))Fatal("burg is not perpendicular to glide plane");
            if(randm(&seed)<0.5){b[0]=-b[0]; b[1]=-b[1]; b[2]=-b[2];}

            theta = 2*asin(maxSeg/2.0/r); nSegs = MAX(4,int(2*M_PI/theta)); theta = 2*M_PI/real8(nSegs);

            inData->nodeCount += nSegs;
            inData->node = (Node_t *)realloc(inData->node, inData->nodeCount*sizeof(Node_t));
            memset(&inData->node[baseNodeID], 0, sizeof(Node_t) * nSegs);

            cross(n,b,g); NormalizeVec(g);
            g[0] *= r; g[1]*=r; g[2]*=r;

            for(int k=0; k<nSegs; k++){
                nodeID = baseNodeID+k;
                nbrID1 = k>0?nodeID-1:baseNodeID+nSegs-1;
                nbrID2 = k<nSegs-1?nodeID+1:baseNodeID;

                RotateWithAxis(n, theta*k, rot);
                Matrix33Vector3Multiply(rot, g, gRot);
                CoordMinus(g,gRot,segDir);
                *totDislocLen += Normal(segDir);

                node = &inData->node[nodeID];
                
                node->myTag.domainID = 0;
                node->myTag.index    = nodeID;
                node->x = cenPos[0]+gRot[0];
                node->y = cenPos[1]+gRot[1];
                node->z = cenPos[2]+gRot[2];
		 	    
                node->constraint = UNCONSTRAINED;
                node->inGrain = inGrain0;
                node->onGB = onGB0;
                
                AllocNodeArms(node, 2);
		 	    
                node->nbrTag[0].domainID = 0;
                node->nbrTag[0].index    = nbrID1;
                node->nbrTag[1].domainID = 0;
                node->nbrTag[1].index    = nbrID2;
                
                node->burgX[0] = b[0];
                node->burgY[0] = b[1];
                node->burgZ[0] = b[2];
                
                node->burgX[1] = -b[0];
                node->burgY[1] = -b[1];
                node->burgZ[1] = -b[2];
                
                node->nx[0] = n[0];
                node->ny[0] = n[1];
                node->nz[0] = n[2];
                	
                node->nx[1] = n[0];
                node->ny[1] = n[1];
                node->nz[1] = n[2];
            }
            baseNodeID += nSegs; 
        }
    }

	IncDislocationDensity(inData, totDislocLen);
	param->nodeCount = inData->nodeCount;
	WriteInitialNodeData(home, inData, 1);
	FreeInNodeArray(inData, inData->nodeCount);
	inData->nodeCount = 0;
	
    return;
}

void SliptArm(InData_t *inData, int id1, int id2, real8 burg[3], real8 n[3], int &newNodeID, int inGrain, int onGB){
   
    inData->nodeCount +=2;
    inData->node = (Node_t *)realloc(inData->node, inData->nodeCount*sizeof(Node_t));
    Node_t *node1 = &inData->node[id1];
    Node_t *node2 = &inData->node[id2];
    memset(&inData->node[inData->nodeCount-1], 0, sizeof(Node_t) * 1);
    Node_t *newNode1 = &inData->node[inData->nodeCount-2];
    Node_t *newNode2 = &inData->node[inData->nodeCount-1];

    newNode1->myTag.domainID = 0;
    newNode1->myTag.index = newNodeID;
    newNode2->myTag.domainID = 0;
    newNode2->myTag.index = newNodeID+1;

    real8 pos1[3]={node1->x,node1->y,node1->z};
    real8 pos2[3]={node2->x,node2->y,node2->z};
    real8 newPos[3], t1=1.0/3.0, t2=2.0/3.0;

    CoordAdd(1.0-t1,pos1,t1,pos2,newPos);
    newNode1->x=newPos[0]; newNode1->y=newPos[1]; newNode1->z=newPos[2];
    CoordAdd(1.0-t2,pos1,t2,pos2,newPos);
    newNode2->x=newPos[0]; newNode2->y=newPos[1]; newNode2->z=newPos[2];

    newNode1->constraint = UNCONSTRAINED;
    newNode2->constraint = UNCONSTRAINED;
    newNode1->inGrain = inGrain; newNode1->onGB = onGB;
    newNode2->inGrain = inGrain; newNode2->onGB = onGB;

    AllocNodeArms(newNode1, 2);
    AllocNodeArms(newNode2, 2);
	
    newNode1->nbrTag[0].domainID = 0;
    newNode1->nbrTag[0].index    = node1->myTag.index;
    newNode1->nbrTag[1].domainID = 0;
    newNode1->nbrTag[1].index    = newNode2->myTag.index;

    newNode2->nbrTag[0].domainID = 0;
    newNode2->nbrTag[0].index    = newNode1->myTag.index;
    newNode2->nbrTag[1].domainID = 0;
    newNode2->nbrTag[1].index    = node2->myTag.index;
    
    newNode1->burgX[0] = -burg[0];
    newNode1->burgY[0] = -burg[1];
    newNode1->burgZ[0] = -burg[2];
    newNode1->burgX[1] =  burg[0];
    newNode1->burgY[1] =  burg[1];
    newNode1->burgZ[1] =  burg[2];
    newNode1->nx[0] = n[0];
    newNode1->ny[0] = n[1];
    newNode1->nz[0] = n[2];
    newNode1->nx[1] = n[0];
    newNode1->ny[1] = n[1];
    newNode1->nz[1] = n[2];

    newNode2->burgX[0] = -burg[0];
    newNode2->burgY[0] = -burg[1];
    newNode2->burgZ[0] = -burg[2];
    newNode2->burgX[1] =  burg[0];
    newNode2->burgY[1] =  burg[1];
    newNode2->burgZ[1] =  burg[2];
    newNode2->nx[0] = n[0];
    newNode2->ny[0] = n[1];
    newNode2->nz[0] = n[2];
    newNode2->nx[1] = n[0];
    newNode2->ny[1] = n[1];
    newNode2->nz[1] = n[2];

    int id12 = GetArmID(node1, node2);
    int id21 = GetArmID(node2, node1); 
    if(id12 == id21 && id12 == -1){
        id12 = node1->numNbrs;
        ReallocNodeArms(node1, id12+1);
        id21 = node2->numNbrs;
        ReallocNodeArms(node2, id21+1);
    }else{
        if(id12 == -1 || id21 == -1)Fatal("can not find the relationship");
    }
    
    node1->nbrTag[id12].domainID = 0;
    node1->nbrTag[id12].index = newNode1->myTag.index;
    node1->burgX[id12] = burg[0];
    node1->burgY[id12] = burg[1];
    node1->burgZ[id12] = burg[2];
    node1->nx[id12] = n[0];
    node1->ny[id12] = n[1];
    node1->nz[id12] = n[2];

    node2->nbrTag[id21].domainID = 0;
    node2->nbrTag[id21].index = newNode2->myTag.index;
    node2->burgX[id21] = -burg[0];
    node2->burgY[id21] = -burg[1];
    node2->burgZ[id21] = -burg[2];
    node2->nx[id21] = n[0];
    node2->ny[id21] = n[1];
    node2->nz[id21] = n[2];

    newNodeID +=2;
    return;
} 
inline bool Perpendicular(Node_t *node1, Node_t *node2){
    int id = GetArmID(node1, node2);
    real8 lDir[3] = {node1->x-node2->x,node1->y-node2->y,node1->z-node2->z}; normalize_vec(lDir);
    real8 n[3] = {node1->nx[id],node1->ny[id],node1->nz[id]};
    return PerpendicularDirection(lDir,n);
}

void CreateFCCPerfectLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs, real8 *totDislocLen)
{
    int         nLoops = inArgs->numLoops;
    int         seed = inArgs->seed;
    const real8 cubeLength = real8(inArgs->cubeLength);
    const real8  radius = inArgs->radius<=sqrt(3)/2?inArgs->radius*cubeLength:inArgs->radius;
    const real8  radius2 = inArgs->radius2;
    const real8 radius3 = inArgs->radius3;
    const real8 invVol=1.0E0/(pow(cubeLength,3.0)*pow(2.556e-10,2));
    const real8 maxDensity = inArgs->maxDensity;
    BurgInfo_t  *bd; Grain_t *grain; GB_t *gb;
    Node_t *node,  *nextNode;
    real8   rList[4];

    int baseNodeID=0, randNum, grainID, burgID, newNodeIndex=0, startRemeshIndex=0;
    inData->nodeCount = 0;
    Param_t *param = inData->param;
    real8   maxSeg = param->maxSeg;
    real8   n0[3], n[4][3], p[4][3], cenPos[3], r, dis, lDirList[4][3],lPoiList[4][3];
    real8   burg[3], gbPos[3], l1[3], l2[3], xi[3], t;
    int     inGrain=-1, onGB=-1, ig, og, onGBList[4];
    const int loopIDs[3][2]={{0,1},{0,1},{1,0}};
    int      lID, pid1, pid2, nodeGBList[4];
    
    home->param->materialType = MAT_TYPE_FCC;

    int i;for(i=0; (i<nLoops) && (maxDensity>invVol*(*totDislocLen)); i++){
        randNum = int(1.0E5*randm(&seed));
        grainID=(nLoops>6*home->nGrains)&&(i<6*home->nGrains)?i/6:(randNum%(home->nGrains*6))/6;
        burgID=(nLoops>6*home->nGrains)&&(i<6*home->nGrains)?i%6:randNum%6;
        grain = home->grainKeys[grainID];
        bd = &grain->burgData;
        VECTOR_COPY(burg,bd->burgList[burgID]);
        if(randm(&seed)>0.5){burg[0]=-burg[0];burg[1]=-burg[1];burg[2]=-burg[2];}
        
        while(true){
            cenPos[0] = (randm(&seed)-0.5)*cubeLength;
            cenPos[1] = (randm(&seed)-0.5)*cubeLength;
            cenPos[2] = (randm(&seed)-0.5)*cubeLength;
            NodeLocation(home, cenPos[0], cenPos[1], cenPos[2], inGrain, onGB,15);
            if((inArgs->polycrystal==1 && inGrain == grainID && onGB <0) || (inArgs->polycrystal==0))break;
        }

        lID = int(1.0E5*randm(&seed))%3; 
        pid1=bd->burgOnPlaneID[burgID][loopIDs[lID][0]];
        pid2=bd->burgOnPlaneID[burgID][loopIDs[lID][1]];
        VECTOR_COPY(n[0],bd->planeList[pid1]);
        VECTOR_COPY(n[1],bd->planeList[pid2]);
        VECTOR_COPY(n[2],bd->planeList[pid1]);
        VECTOR_COPY(n[3],bd->planeList[pid2]);

        rList[0] = (param->polycrystal>0 && pid1<3)?radius3:radius+randm(&seed)*(radius2-radius);
        rList[1] = (param->polycrystal>0 && pid2<3)?radius3:radius+randm(&seed)*(radius2-radius);
        rList[2] = (param->polycrystal>0 && pid1<3)?radius3:radius+randm(&seed)*(radius2-radius);
        rList[3] = (param->polycrystal>0 && pid2<3)?radius3:radius+randm(&seed)*(radius2-radius);
        
#if 0
        if(inArgs->polycrystal==1 && (burgID%2==1)){
            while(true){
                int loopPlane = int(randm(&seed)*1.0E5)%7;
                if(loopPlane == bd->burgOnPlaneID[burgID][0] || loopPlane == bd->burgOnPlaneID[burgID][1] ||
                   loopPlane == bd->burgOnPlaneID[burgID][2])continue; 
                VECTOR_COPY(n0, bd->planeList[loopPlane]);break;
            }
        }else{
            if(inArgs->polycrystal==1){
                VECTOR_COPY(n0, bd->planeList[3]);        
            }else{
                while(true){
                    int loopPlane = int(randm(&seed)*1.0E5)%7;
                    if(loopPlane == bd->burgOnPlaneID[burgID][0] || loopPlane == bd->burgOnPlaneID[burgID][1] ||
                       loopPlane == bd->burgOnPlaneID[burgID][2])continue; 
                    VECTOR_COPY(n0, bd->planeList[loopPlane]);break;
                }
            }
        }
#else
        cross(n[0],burg,l1);
        cross(n[1],burg,l2);

        NormalizeVec(l1); NormalizeVec(l2);
        DefineN0:
        t=randm(&seed);
        CoordAdd(t, l1, 1.0-t, l2, xi);
        NormalizeVec(xi);
        t=randm(&seed)/2;

        CoordAdd(t, xi, 1.0-t, burg, n0);
        NormalizeVec(n0);


        if(param->polycrystal>0 && (burgID%2)==0 && fabs(Dot(param->TBNormal,n0))<0.8) goto DefineN0;
#endif
        CoordAdd(1.0, cenPos, -rList[0], &n[0][0], &p[0][0]);
        CoordAdd(1.0, cenPos, -rList[1], &n[1][0], &p[1][0]);
        CoordAdd(1.0, cenPos,  rList[2], &n[2][0], &p[2][0]);
        CoordAdd(1.0, cenPos,  rList[3], &n[3][0], &p[3][0]);

        for(int j=0; j<4; j++){
            onGBList[j] = -1;
            if(LINE_INTERSECTION!=PlanePlaneIntersection(n0,cenPos,n[j],p[j],lDirList[j],lPoiList[j])){
                Fatal("can not find line intersection");
            }
            CoordMinus(lPoiList[j],cenPos,lDirList[j]); normalize_vec(lDirList[j]);
            CoordAdd(1.0,cenPos,rList[j],lDirList[j],lPoiList[j]);
            NodeLocation(home, lPoiList[j][0],lPoiList[j][1],lPoiList[j][2], ig, og,15);
            if((og<0 && inGrain==ig) || (!AlignedDirection(bd->planeList[3],n[j]))){
                cross(n0,n[j],lDirList[j]); normalize_vec(lDirList[j]);
                continue;
            }
        
            r=10.0;
            while(true){
                CoordAdd(1.0,cenPos,r,lDirList[j],lPoiList[j]);
                NodeLocation(home, lPoiList[j][0],lPoiList[j][1],lPoiList[j][2], ig, og,15);
                if(ig!=inGrain)break;
                r+=10.0;
            }
            ClosedGB(home, lPoiList[j], &gb); VECTOR_COPY(gbPos, gb->p);
            if(lPoiList[j][2]<param->zBoundMin)gbPos[2]+=(2.0*param->zBoundMin);
            else if(lPoiList[j][2]>param->zBoundMax)gbPos[2]+=(2.0*param->zBoundMax);

            real8 vec[3]; if(POINT_INTERSECTION!=LinePlaneIntersection(lDirList[j],cenPos,gb->n,gbPos,vec,lPoiList[j])){
                Fatal("can not find point intersection");
            }
            onGBList[j] = gb->GBID;
            cross(n0,n[j],lDirList[j]); normalize_vec(lDirList[j]);
        }
        
        for(int j=0; j<4; j++){
            real8 vec[3], dis; int type;
            type =LineLineIntersection(lDirList[j],lPoiList[j],lDirList[(j==3?0:j+1)],lPoiList[(j==3?0:j+1)],vec,p[j],&dis,1.0E-3,1.0E-4);
            if(type != POINT_INTERSECTION){
                Fatal("can not find point intersection!, dis=%g, type =%d, j=%d",dis,type,j);
            }
            int nextIndex = j==3?0:j+1;
            if(onGBList[j]==-1 && onGBList[nextIndex]==-1){
                nodeGBList[j]=-1;
            }else{
                nodeGBList[j] = (onGBList[j]>-1?onGBList[j]:onGBList[nextIndex]);
            }
        }        
        bool dangerLoop = (pid1!=3 && pid2!=3 && !AlignedDirection(n0, bd->planeList[3]))?true:false;
        if(dangerLoop){
            int j; for(j=0; j<4; j++){
                NodeLocation(home, p[j][0],p[j][1],p[j][2], ig, og,15); if(ig!=inGrain)break;
            }
            if(j<4){
                real8 ldir1[3]; CoordMinus(p[2],p[0],ldir1); normalize_vec(ldir1);
                real8 ldir2[3]; CoordMinus(p[3],p[1],ldir2); normalize_vec(ldir2);
                real8 vec1[3], vec2[3], lList[4][3];
                int k,l; real8 factor=1.0;
                if(LineLineIntersection(ldir1,p[0],ldir2,p[1],vec1,vec2,&dis)!=POINT_INTERSECTION){
                    Fatal("%d: can not find intersection",__LINE__);
                }
                CoordMinus(cenPos,vec2,vec1); 
                for(k=0; k<4; k++){CoordAdd(p[k],vec1,p[k]); CoordMinus(p[k],cenPos,lList[k]);}
                for(k=0;k<4;k++){
                    NodeLocation(home, p[k][0],p[k][1],p[k][2], ig, og,15); if(ig!=inGrain)break;
                }
                if(k<4){
                    while(ig!=inGrain||og>-1){
                        factor *= 0.98;
                        for(l=0; l<4; l++){
                            CoordAdd(1.0,cenPos,factor,lList[l],p[l]);
                            NodeLocation(home, p[l][0],p[l][1],p[l][2], ig, og,15);
                            if(ig!=inGrain||og>-1)break;
                        }
                    }
//                  printf("interation succeed. factor=%g\n",factor);
                }
            }
            if(inArgs->polycrystal !=0){
                for(int gbID=0; gbID<grain->nBounds; gbID++){
                    gb = grain->GB[gbID];
                    real8 dis1=PointPlaneDistance(cenPos, gb->n, gb->p);
                    for(j=0; j<4; j++){if(dis1*PointPlaneDistance(p[j],gb->n,gb->p) < 0)break;}
                    if(j<4){
                        real8 ldir1[3]; CoordMinus(p[2],p[0],ldir1); normalize_vec(ldir1);
                        real8 ldir2[3]; CoordMinus(p[3],p[1],ldir2); normalize_vec(ldir2);
                        real8 vec1[3], vec2[3], lList[4][3];
                        int k,l; real8 factor=1.0E-10;
                        if(LineLineIntersection(ldir1,p[0],ldir2,p[1],vec1,vec2,&dis)!=POINT_INTERSECTION){
                            Fatal("%d: can not find intersection",__LINE__);
                        }
                        CoordMinus(cenPos,vec2,vec1); 
                        for(k=0; k<4; k++){
                            CoordAdd(p[k],vec1,p[k]); 
                            CoordMinus(p[k],cenPos,lList[k]); 
                        }
                        l=4; ig=inGrain; og=-1;
                        while(ig==inGrain && og==-1 && l==4){
                            factor *= 1.1;
                            for(l=0; l<4; l++){
                                CoordAdd(1.0,cenPos,factor,lList[l],p[l]);
                                NodeLocation(home, p[l][0],p[l][1],p[l][2], ig, og,15);
                                if(ig!=inGrain||og>-1)break;
                            }
                        }
                        factor /= 1.1;
                        for(l=0; l<4; l++){
                            CoordAdd(1.0,cenPos,factor,lList[l],p[l]);
                            NodeLocation(home, p[l][0],p[l][1],p[l][2], ig, og,15);
                            if(ig!=inGrain || onGB>-1)Fatal("shit {%g,%g,%g}",p[l][0],p[l][1],p[l][2]);
                        }
                    
                    }
                }
            }
        }
            if(inArgs->polycrystal !=0){
                int j;
                for(int gbID=0; gbID<grain->nBounds; gbID++){
                    gb = grain->GB[gbID];
                    real8 dis1=PointPlaneDistance(cenPos, gb->n, gb->p);
                    for(j=0; j<4; j++){
                        if(nodeGBList[j]>-1)continue;
                        if(dis1*PointPlaneDistance(p[j],gb->n,gb->p) < 0.0)
                        break;
                    }

                    if(j<4){
                        real8 ldir1[3]; CoordMinus(p[2],p[0],ldir1); normalize_vec(ldir1);
                        real8 ldir2[3]; CoordMinus(p[3],p[1],ldir2); normalize_vec(ldir2);
                        real8 vec1[3], vec2[3], lList[4][3];
                        int k,l; real8 factor=1.0E-10;
                        if(LineLineIntersection(ldir1,p[0],ldir2,p[1],vec1,vec2,&dis)!=POINT_INTERSECTION){
                            Fatal("%d: can not find intersection",__LINE__);
                        }
                        CoordMinus(cenPos,vec2,vec1); 
                        for(k=0; k<4; k++){
                            CoordAdd(p[k],vec1,p[k]); 
                            CoordMinus(p[k],cenPos,lList[k]); 
                        }
                        l=4; ig=inGrain; og=-1;
                        while(ig==inGrain && og==-1 && l==4){
                            factor *= 1.1;
                            for(l=0; l<4; l++){
                                CoordAdd(1.0,cenPos,factor,lList[l],p[l]);
                                NodeLocation(home, p[l][0],p[l][1],p[l][2], ig, og,15);
                                if(ig!=inGrain||og>-1)break;
                            }
                        }
                        factor /= 1.1;
                        for(l=0; l<4; l++){
                            CoordAdd(1.0,cenPos,factor,lList[l],p[l]);
                            NodeLocation(home, p[l][0],p[l][1],p[l][2], ig, og,15);
                            nodeGBList[l] = -1;
                            if(ig!=inGrain || onGB>-1)Fatal("shit {%g,%g,%g}",p[l][0],p[l][1],p[l][2]);
                        }
                    }
                }
            }

        newNodeIndex = inData->nodeCount;
        inData->nodeCount += 4;
        inData->node = (Node_t *)realloc(inData->node, inData->nodeCount*sizeof(Node_t));
        memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * 4);

        for(int j=0; j<4; j++){
            if(!PerpendicularDirection(n[j],lDirList[j]))Fatal("n+_l");
            if(!PerpendicularDirection(n[j],burg))Fatal("n+_burg");
            node = &inData->node[baseNodeID+j];
            int nextIndex = j==3?0:j+1;
            int lastIndex = j==0?3:j-1;

            node->myTag.domainID = 0;
            node->myTag.index = baseNodeID+j;
            node->x = p[j][0];
            node->y = p[j][1];
            node->z = p[j][2];

            real8 ls[3]; CoordMinus(p[j],p[nextIndex],ls); normalize_vec(ls);
            if(!PerpendicularDirection(ls, n[nextIndex]))Fatal("kkk");

            node->constraint = UNCONSTRAINED;
            if(nodeGBList[j] == -1){
                node->inGrain = inGrain; node->onGB = -1;
            }else{
                node->inGrain =-1; node->onGB = nodeGBList[j];
            }
            NodeLocation(home, node->x, node->y, node->z, ig,og);
            if(node->onGB != og){
                FormatVector(p[j],"p[j]");
                FormatVector(n[j],"n[j]");
                FormatVector(n[nextIndex],"n[nt]");
                Fatal("j=%d: [%d,%d]!=[%d,%d] onGBList[%d]=%d, onGBList[%d]=%d",j,node->inGrain,node->onGB,ig,og,j,onGBList[j],nextIndex,onGBList[nextIndex]);
            }
            
            AllocNodeArms(node, 2);
		 	
            node->nbrTag[0].domainID = 0;
            node->nbrTag[0].index    = lastIndex+baseNodeID;
            node->nbrTag[1].domainID = 0;
            node->nbrTag[1].index    = nextIndex+baseNodeID;
            
            node->burgX[0] = -burg[0];
            node->burgY[0] = -burg[1];
            node->burgZ[0] = -burg[2];
            
            node->burgX[1] =  burg[0];
            node->burgY[1] =  burg[1];
            node->burgZ[1] =  burg[2];
            
            node->nx[0] = n[j][0];
            node->ny[0] = n[j][1];
            node->nz[0] = n[j][2];
            	
            node->nx[1] = n[nextIndex][0];
            node->ny[1] = n[nextIndex][1];
            node->nz[1] = n[nextIndex][2];
        }

        int newNodeID = inData->nodeCount;
        for(int j=0; j<4; j++){
            int nextIndex = j==3?0:j+1;
            node = &inData->node[baseNodeID+j];
            nextNode = &inData->node[baseNodeID+nextIndex];
            if(!Perpendicular(nextNode, node)) Fatal("%s: haha",__func__);
   
            if(node->onGB>-1 && nextNode->onGB>-1){
                if(node->onGB!=nextNode->onGB){
                    SliptArm(inData, baseNodeID+j, baseNodeID+nextIndex, burg, n[nextIndex], newNodeID, inGrain, onGB);
                }else{
                    real8 b1[3], b2[3];
                    std::vector<std::vector<Dislocation_t> >    disArray;
                    std::vector<std::vector<real8> >            burgList;

                    if(!GeometricCriterion(bd, burg, lDirList[nextIndex], n[nextIndex], 2,
                       burgList, disArray)){
                        Fatal("something wrong with geometric criterion");
                    }
                    if(disArray.size() != 1)printf("Warning: the extened dislcoation has some wrong\n");
                    if(disArray[0].size() != 2)printf("Warning: the extened dislcoation has some wrong\n");
                    
                    VECTOR_COPY(b1, burgList[disArray[0][0].burgID]);
                    VECTOR_COPY(b2, burgList[disArray[0][1].burgID]);
                    SliptArm(inData, baseNodeID+j, baseNodeID+nextIndex, b1, n[nextIndex], newNodeID, node->inGrain, node->onGB);
                    node = &inData->node[baseNodeID+j];
                    nextNode = &inData->node[baseNodeID+nextIndex];
                    SliptArm(inData, baseNodeID+j, baseNodeID+nextIndex, b2, n[nextIndex], newNodeID, node->inGrain, node->onGB);
                }
            }else{
                SliptArm(inData, baseNodeID+j, baseNodeID+nextIndex, burg, n[nextIndex], newNodeID, inGrain, onGB);
            }
        }

	    IncDislocationDensity(inData, startRemeshIndex, totDislocLen);
        InitRemesh(inData, 0, startRemeshIndex);
        baseNodeID = inData->nodeCount;
        startRemeshIndex = baseNodeID;
    } 
    printf("%d have beend distributed\n",i);

	param->nodeCount = inData->nodeCount;
	WriteInitialNodeData(home, inData, 1);
	FreeInNodeArray(inData, inData->nodeCount);
	inData->nodeCount = 0;
    return;   
}


void CreateSingleArmSource(Home_t *home, InData_t *inData, InArgs_t *inArgs, real8 *totDislocLen)
{
        int      i, lineIndex, chain, baseNodeID, numPoints=3,ranNum;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      x0, d, nLayers, xn, layer, time;
        int      startRemeshIndex = 0;
        real8    srcLen, thickness,lenRange;
        real8    p0[3], p1[3], p2[3], p3[3], p4[3], randomLine[3];
        real8	 randomAngle, edge, srcLenMin, srcLenMax;
        int	     cubeLength, numSources,  seed, dislocType;
        int	     polycrystal, nGrain, finish;
        int      firstBurgIndex, srcLenType;
        int      inGrain1, onGB1, inGrain2, onGB2, inGrain0, onGB0;
        int      inGrain3, onGB3, inGrain4, onGB4;
        real8	 twinRatio, maxSeg, maxDensity, burgMag = 2.556E-10;
        real8	 volume, glideDir[3];
        real8	 tnx[4], tny[4], tnz[4], burg[12][3], cubeSize,heightSize;
        real8	 n0[3], b[3], n1[3];	/* plane normal vector, burgers, normal vector*/
        Node_t   *node;
        Param_t  *param;
        real8	 invsq2, basePos[3], b1[3], b2[3], f1[3], f2[3], f3[3], f4[3];
        real8	 rotMatrix[3][3],rotLine[3],rotBurg[3],rotPlane[3], btemp[3];
        int      intID, burgID, planeID;
        real8    pi = 3.1415926535897932385, totSrcLen;
        real8    separateLength = 40.0;
        bool        finished, singleSrcTest = 0;
        BurgInfo_t  *bd, *gbd;
        GB_t        *gb;
        Grain_t     *grain;
        int constraint2, constraint4;
        	
        real8                                       str[3][3], gbNormal[3];
        int                                         appliedDir, j, k;
        real8                                       sign;
        std::vector<std::vector<Dislocation_t> >    disArray;
        std::vector<std::vector<real8> >            burgList;
        std::vector<std::vector<Dislocation_t> >::iterator      itDisList;

        real8   normalABC[3], lambda_FR, lambda_SAS1, lambda_GB1, lambda_SAS2, lambda_GB2, interDir[3];
        real8   barrierFactor, crossProd[3], vec[3];
        real8   vec1[3], vec2[3], s0,  paras[4][4], zaxis[3] = {0,0,1}, schmid;
        int l; real8 dis;

        normalABC[0] = 1;
        normalABC[1] = 1;
        normalABC[2] = 1;

#ifdef SAS_ANALYZE
        FILE    *fp;
        
        if((fp = fopen(inArgs->outFile, "w")) == (FILE *)NULL){
            Fatal("can not open file %s", inArgs->outFile);
        }
#endif

/*
 *      Grap the command-line argument values needed for this function
 */
        home->param->materialType = MAT_TYPE_FCC;
        	
        cubeLength  = inArgs->cubeLength;
	    srcLenType  = inArgs->srcLenType;
        nLayers     = inArgs->nLayers;
       	maxDensity  = inArgs->maxDensity;
        numSources  = inArgs->numFRSrcs;
        srcLenMin   = inArgs->frLenMin;
        srcLenMax   = inArgs->frLenMax;
        seed        = inArgs->seed;
        dislocType  = inArgs->type;
        polycrystal = inArgs->polycrystal;
        nGrain	    = inArgs->nGrain;
        twinRatio   = inArgs->twinRatio;
        maxSeg      = (real8)inArgs->maxSegLen;
        burgID      = inArgs->burgID;
        planeID     = inArgs->planeID;
        	
        param = inData->param;
        d = inArgs->deltaNums;
        cubeSize	=	(real8)cubeLength;
        if(srcLenMax<0) {srcLenMax = -srcLenMax; srcLenType = 2;} 
        if(srcLenMin<0) {srcLenMin = -srcLenMin; srcLenType = 2;}
        lenRange = srcLenMax - srcLenMin;

        bool freeLine=(param->xBoundType != Periodic && param->yBoundType != Periodic && param->zBoundType != Periodic);
        if(freeLine==false)Fatal("Singe arm source can only be distributed in finite volume");

        real8 n3[3];

        if(numSources <= 0){
            singleSrcTest = 1;
        }
        	
        if(polycrystal){
            gb = home->GBKeys[0];
            gbd = &(gb->burgData);
        }
        n1[0] = 0.0;	n1[1] = 0.0;	n1[2] = 0.0;
        real8 maximalLen, a1, a2, randmVal;
        bool partialDislocation;

/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        baseNodeID = 0;
        inData->nodeCount = 0;
		
        heightSize = 2 * param->zBoundMax;
	    volume = heightSize * cubeSize * cubeSize;
        thickness = heightSize / ((real8)nLayers);
        	
        x0 = (d*nLayers*(1 - nLayers) + 2*numSources)/(2 * nLayers);
	    totSrcLen = 0.0;
        real8 primarySrcLen = 0.0;
	    finish = 0;

        int nSlips;
        if(param->polycrystal>0){ if(inArgs->nTDs>0)nSlips=N_SS; else nSlips=N_LD; }else nSlips=N_SG;
        real8 *proba, *primary, *Proba, *percen, maxSchmid=0;
        proba = (real8 *)malloc(nSlips*sizeof(real8));
        Proba = (real8 *)malloc(nSlips*sizeof(real8));
        percen = (real8 *)malloc(nSlips*sizeof(real8));
        primary = (real8 *)malloc(nSlips*sizeof(real8));

        for(i=0; i<nSlips; i++){
            j = slip_sys[i][2]>-1?slip_sys[i][2]:0;
            bd=GetBurgData(home, j, -1);
            proba[i] = fabs(DotProduct(inArgs->loadDir, bd->burgList[slip_sys[i][0]])*
                            DotProduct(inArgs->loadDir,bd->planeList[slip_sys[i][1]])/Normal(bd->burgList[slip_sys[i][0]]));
            if(proba[i]>maxSchmid){
                maxSchmid = proba[i];
            }
            if(home->nGBs>1)Proba[i]=proba[i]; else Proba[i]=1.0;
            percen[i]=0.0;
        }
        for(i=0;i<nSlips;i++)if(fabs(proba[i]-maxSchmid)<1E-4)primary[i]=1.0; else primary[i]=MAX(proba[i]/maxSchmid, 0.33333);
        for(i=1;i<nSlips;i++)Proba[i]+=Proba[i-1];
        for(i=0;i<nSlips;i++)Proba[i]/=Proba[nSlips-1];
        
        for (layer = 0; layer < nLayers; layer++) {

	        if(finish)break;

            xn = layer * d + x0;
            if(xn <= 0){
                xn = 0;
                continue;
            }
            	
            for(chain = 0; chain < xn; ){
 

                if(xn > nSlips && chain < nSlips && (2*maxDensity*volume*burgMag*burgMag/srcLenMax > nSlips)){
                    ranNum = chain%nSlips;
                }else{
                    randmVal= randm(&seed);
                    for(i=0;i<nSlips;i++)if(randmVal<Proba[i]&& (randmVal>((i==0)?0:Proba[i-1])) )ranNum=i;
                }
                gpBurgIndex = slip_sys[ranNum][0]; gpIndex=slip_sys[ranNum][1];
    
                finished = 0;
                while(!finished){
                    p1[0] = (randm(&seed)-0.5) * cubeSize*0.99;
                    p1[1] = (randm(&seed)-0.5) * cubeSize*0.99;			
                    p1[2] = 0.99*randm(&seed)* thickness + ((real8)layer) * thickness - 0.5*heightSize;

                    NodeLocation(home, p1[0],  p1[1], p1[2], inGrain1, onGB1,15);
                    bd = GetBurgData(home, inGrain1, onGB1);
                    if(onGB1 > -1)continue;

#if 1
                    if(gpBurgIndex>=6){
                        ClosedGB(home, p1, &gb);
                        bd = &(gb->burgData);
                        PointPlaneDistance(p1, gb->n, gb->p, dis, p1);
                        inGrain1 = -1; onGB1 = gb->GBID;
                    }else if(gpIndex!=3){
                        if(inArgs->polycrystal ==1 && (inGrain1%2)!=slip_sys[ranNum][2]) continue;
                    }
#else
                    gpIndex = 1;
                    gpBurgIndex = 5;
//                  p1[2] =  home->grainKeys[inGrain1]->cz;
                    VECTOR_ZERO(p0);
#endif   
                    VECTOR_COPY(n0, bd->planeList[gpIndex]);
                    VECTOR_COPY(b , bd->burgList[gpBurgIndex]);
                    VECTOR_COPY(n3, bd->planeList[3]);
                    if(!PerpendicularDirection(n0,b))Fatal("burg is not perpendicular to glide plane");

                    time = 0;
                    while(true){
                        if(time>100)break;
                        time++;
                        randomAngle	= acos((randm(&seed) - 0.5) * 2.0);
                        randomAngle = ((randm(&seed) > 0.5) ? randomAngle+pi : randomAngle);
                        
                        cross(n0, b, n1);
                        NormalizeVec(n1);
                	    
                        for(i = 0; i < 3; i++){
                            randomLine[i] = n1[i] * sin(randomAngle) + b[i] * cos(randomAngle);
                        }
                        NormalizeVec(randomLine);

                        if(!IntersectionOfRayBox(param, randomLine, p1, p2, constraint2))continue;
                        p0[0]=0.5*(p1[0]+p2[0]);
                        p0[1]=0.5*(p1[1]+p2[1]);
                        p0[2]=0.5*(p1[2]+p2[2]);
                        NodeLocation(home, p0[0],  p0[1], p0[2], inGrain0, onGB0,15);
                        NodeLocation(home, p2[0],  p2[1], p2[2], inGrain2, onGB2,15);
                        if(inGrain1 != inGrain0 || inGrain1 != inGrain2 || onGB1 != onGB0 || onGB1 != onGB2)continue;

                        if(inArgs->extended != 0){
                            VECTOR_COPY(basePos, p0);
                            cross(n0, randomLine, glideDir);
                            NormalizeVec(glideDir);
                            p0[0] -= (separateLength*glideDir[0]*0.5); // drag partial dislocation
                            p0[1] -= (separateLength*glideDir[1]*0.5);
                            p0[2] -= (separateLength*glideDir[2]*0.5);

                            p3[0] = p0[0] + separateLength*glideDir[0]; // leading partial dislocation
                            p3[1] = p0[1] + separateLength*glideDir[1];
                            p3[2] = p0[2] + separateLength*glideDir[2];


                            if(!IntersectionOfRayBox(param, randomLine, p3, p4, constraint4))continue;

                            NodeLocation(home, p0[0],  p0[1], p0[2], inGrain0, onGB0,15);
                            NodeLocation(home, p3[0],  p3[1], p3[2], inGrain3, onGB3,15);
                            NodeLocation(home, p4[0],  p4[1], p4[2], inGrain4, onGB4,15);
                    
                            if(onGB0 > -1 || onGB3 > -1 || inGrain0 != inGrain2 || inGrain3 != inGrain2 || inGrain3!=inGrain4 || onGB4!=onGB3)continue;
                        }   

                        break;
                    }
                    if(time>100)continue;
                    break;
                }
                srcLen = PointPointDistance(p1,p2);
                if(inArgs->extended != 0){
                    if(!PartialDislocations(bd, b, n0, b1, b2)){
                        partialDislocation = true;
                    }else {
                        partialDislocation = false;
	                    totSrcLen += srcLen;
                    }
                }else partialDislocation= false;

                if(primary[ranNum]>0.99) primarySrcLen += srcLen;
                percen[ranNum] += srcLen;
	            totSrcLen += srcLen;

#ifdef SAS_ANALYZE
                if(polycrystal){
                    lambda_FR = srcLen/2;
                    GetFreePath(param, p1, p2, n0, home->grainKeys[inGrain0], 
                                interDir, lambda_SAS1, lambda_GB1, lambda_SAS2, lambda_GB2);

                    ExtPKForce(str, b[0], b[1], b[2], 
                               p1[0], p1[1], p1[2],
                               p2[0], p2[1], p2[2], f1, f2);
                   
                    real8 dis; 
                    PointPlaneIntersection(p0, gb->n, gb->p, vec, &dis);
                    vec[0] -= p0[0];
                    vec[1] -= p0[1];
                    vec[2] -= p0[2];
                    if(DotProduct(vec, f1) < 0){
                        lambda_GB1 = -lambda_GB1;
                        lambda_GB2 = -lambda_GB2;
                    }

                    s0 = DotProduct(zaxis, n0)*DotProduct(zaxis, b);
                    VECTOR_ZERO(paras[0]);
                    VECTOR_ZERO(paras[1]);
                    VECTOR_ZERO(paras[2]);
                    VECTOR_ZERO(paras[3]);

                    real8   schmidList[2] = {0,0}, maxSchmid = -1;
                    int     maxSet = 0;
                
                    if(lambda_GB1 != 0.0 && fabs(s0) > 1.0E-2){
                        std::vector<std::vector< Dislocation_t> >().swap(disArray);
                        std::vector<std::vector<double> >().swap(burgList);

                        if(!GeometricCriterion(gbd, b, interDir, n0, 1,
                            burgList, disArray)){
                            Fatal("something wrong with geometric criterion");
                        }


/*
 *                      Now, we don't need the resultant configure including a lattice burg on GB,
 *                      or a partial burg in Grain.
 */
                        for(j=0; j<disArray.size(); ){
                            barrierFactor = -1.0;
                            for(k=0; k<disArray[j].size(); k++){
                                barrierFactor += disArray[j][k].energyFactor;
                                cross(gb->n, gbd->planeList[disArray[j][k].planeID], crossProd);
                                VECTOR_COPY(vec, burgList[disArray[j][k].burgID]);
                                if(Normal(crossProd) < 1.0E-2){
                                    if(BurgType_FCC(vec) == BASIC_BURG)break;
                                }else{
                                    if(BurgType_FCC(vec) != BASIC_BURG)break;
                                }
                            }
  
                            if(k != disArray[j].size() || j > 0 && 
                               barrierFactor > 0){
                                itDisList = disArray.begin() + j;
                                disArray.erase(itDisList);
                            }else{
                                j++;
                            }
                        }
                        if(disArray.size()>2)Fatal("someting wrong with disArray");
                    
                        for(j=0; j<disArray.size(); j++){
                            if(disArray[j].size()>2)Fatal("someting wrong with disArray[%d]", j);     

                            for(k=0; k<disArray[j].size(); k++){
                                paras[2*j+k][0] = disArray[j][k].energyFactor;
                                VECTOR_COPY(vec, burgList[disArray[j][k].burgID]);
                                paras[2*j+k][1] = Normal(vec);
                                NormalizeVec(vec);
                                schmid = DotProduct(zaxis, vec)*DotProduct(zaxis, gbd->planeList[disArray[j][k].planeID]);

	                            if(fabs(paras[2*j+k][1] - 1.0E0) < 1.0E-2){
                                    paras[2*j+k][2] = (schmid*s0 > 0) ? schmid : 0;
	                            }else{
                                    paras[2*j+k][2] = schmid;
                                }
                                paras[2*j+k][2] = fabs(paras[2*j+k][2]);
                                
                            }
                            schmidList[j] = DotProduct(zaxis,b)*DotProduct(zaxis, gbd->planeList[disArray[j][0].planeID]);
                        }
                    }
                    
                    for(j=0; j<2; j++){
                        if(fabs(schmidList[j]) > maxSchmid){
                            maxSet = j;
                            maxSchmid = fabs(schmidList[j]);
                        }
                    }

                    fprintf(fp,"%e %e %e %e   ", lambda_FR, lambda_SAS1, lambda_GB1, fabs(s0));
                    for(j=2*maxSet; j<2*maxSet+2; j++){
                        fprintf(fp, "%e %e %e ", paras[j][0], paras[j][1], paras[j][2]);
                    }
                    fprintf(fp, "\n");

                    fprintf(fp,"%e %e %e %e ", lambda_FR, lambda_SAS2, lambda_GB2, fabs(s0));
                    for(j=2*maxSet; j<2*maxSet+2; j++){
                        fprintf(fp, "%e %e %e ", paras[j][0], paras[j][1], paras[j][2]);
                    }
                    fprintf(fp, "\n");

                }
#endif 
/*           
 *              Reallocate the node array with sufficient size to add
 *              all the new nodes defining this chain.
 */          
                numPoints = (inArgs->extended == 0 || partialDislocation) ? 3 : 5;
                newNodeIndex = inData->nodeCount;
                inData->nodeCount += numPoints;
                inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                                   * sizeof(Node_t));
                memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
		 	   
                if(inArgs->extended == 0 || partialDislocation){
/*           
 *                   Set up the 3 nodes we're using to define the dislocation line.
 *                   First do point p0.
 */               
                     node = &inData->node[baseNodeID];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID;
                     node->x = p0[0];
                     node->y = p0[1];
                     node->z = p0[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain0;
                     node->onGB = onGB0;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = -b[0];
                     node->burgY[0] = -b[1];
                     node->burgZ[0] = -b[2];
                  
                     node->burgX[1] = b[0];
                     node->burgY[1] = b[1];
                     node->burgZ[1] = b[2];
                  
/*                
 *                   For screw dislocations, use the glide plane from the table.
 *                   For edge dislocations, glide plane is the cross product
 *                   of the burgers vector and line direction.
 */               
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1;
                  
                     node->x = p1[0];
                     node->y = p1[1];
                     node->z = p1[2];
                  
                     node->constraint = PINNED_NODE;
                     node->inGrain = inGrain1;
                     node->onGB = onGB1;
                     AllocNodeArms(node, 1);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = b[0];
                     node->burgY[0] = b[1];
                     node->burgZ[0] = b[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                  
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2;
                  
                     node->x = p2[0];
                     node->y = p2[1];
                     node->z = p2[2];
                  
                     node->constraint = constraint2;
                     node->inGrain = inGrain2;
                     node->onGB = onGB2;
                     AllocNodeArms(node, 1);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = -b[0];
                     node->burgY[0] = -b[1];
                     node->burgZ[0] = -b[2];
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                }else{
                     node = &inData->node[baseNodeID];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID;
                     node->x = p0[0];
                     node->y = p0[1];
                     node->z = p0[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain0;
                     node->onGB = onGB0;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 2;
                  
                     node->burgX[0] = -b1[0];
                     node->burgY[0] = -b1[1];
                     node->burgZ[0] = -b1[2];
                     node->signSF[0] = -1.0;
                  
                     node->burgX[1] = b1[0];
                     node->burgY[1] = b1[1];
                     node->burgZ[1] = b1[2];
                     node->signSF[1] = 1.0;
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
                  
/*                
 *                   Now point p1...
 */               
                     node = &inData->node[baseNodeID+1];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+1;
                  
                     node->x = p1[0];
                     node->y = p1[1];
                     node->z = p1[2];
                  
                     node->constraint = PINNED_NODE;
                     node->inGrain = inGrain1;
                     node->onGB = onGB1;
                     AllocNodeArms(node, 2);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID+3;
                  
                     node->burgX[0] = b1[0];
                     node->burgY[0] = b1[1];
                     node->burgZ[0] = b1[2];
                     node->signSF[0] = 1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                  
                     node->burgX[1] = b2[0];
                     node->burgY[1] = b2[1];
                     node->burgZ[1] = b2[2];
                     node->signSF[1] = -1.0;
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];
/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+2];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+2;
                  
                     node->x = p2[0];
                     node->y = p2[1];
                     node->z = p2[2];
                  
                     node->constraint = constraint2;
                     node->inGrain = inGrain2;
                     node->onGB = onGB2;
                     AllocNodeArms(node, 1);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID;
                  
                     node->burgX[0] = -b1[0];
                     node->burgY[0] = -b1[1];
                     node->burgZ[0] = -b1[2];
                     node->signSF[0] = -1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];

/*                
 *                   Now point p3...
 */               
                     node = &inData->node[baseNodeID+3];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+3;
                     node->x = p3[0];
                     node->y = p3[1];
                     node->z = p3[2];
		 	        
                     node->constraint = UNCONSTRAINED;
                     node->inGrain = inGrain3;
                     node->onGB = onGB3;
                  
                     AllocNodeArms(node, 2);
		 	        
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID + 1;
                     node->nbrTag[1].domainID = dislocType;
                     node->nbrTag[1].index    = baseNodeID + 4;
                  
                     node->burgX[0] = -b2[0];
                     node->burgY[0] = -b2[1];
                     node->burgZ[0] = -b2[2];
                     node->signSF[0] = 1.0;
                  
                     node->burgX[1] = b2[0];
                     node->burgY[1] = b2[1];
                     node->burgZ[1] = b2[2];
                     node->signSF[1] = -1.0;
                  
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];
                     	
                     node->nx[1] = n0[0];
                     node->ny[1] = n0[1];
                     node->nz[1] = n0[2];

/*                
 *                   Now point p2...
 */               
                     node = &inData->node[baseNodeID+4];
                  
                     node->myTag.domainID = dislocType;
                     node->myTag.index    = baseNodeID+4;
                  
                     node->x = p4[0];
                     node->y = p4[1];
                     node->z = p4[2];
                  
                     node->constraint = constraint4;
                     node->inGrain = inGrain4;
                     node->onGB = onGB4;
                     AllocNodeArms(node, 1);
                  
                     node->nbrTag[0].domainID = dislocType;
                     node->nbrTag[0].index    = baseNodeID+3;
                  
                     node->burgX[0] = -b2[0];
                     node->burgY[0] = -b2[1];
                     node->burgZ[0] = -b2[2];
                     node->signSF[0] = -1.0;
                     node->nx[0] = n0[0];
                     node->ny[0] = n0[1];
                     node->nz[0] = n0[2];

                }
/*           
 *              The initial segments created are not necessarily limited to
 *              param->maxSegLen, so a call to InitRemesh() is needed to
 *              chop any excessively long segments into proper lengths.
 *              When we've generated the nodal data for the final chain,
 *              write the block of nodal data to the file.
 */	         
                InitRemesh(inData, dislocType, startRemeshIndex);
//                lastBlock = (chain == (numSources - 1));
             
               
                chain++;
                baseNodeID = inData->nodeCount;
                startRemeshIndex = baseNodeID;

	            if(totSrcLen> maxDensity*volume*burgMag*burgMag){
		            finish = 1;
		            break;
	            }else{
		            finish = 0;
	            }
            }  /* for (chain = 0; chain < numChains; ) */
		
	    }/* for(layer=0; layer<nLayers;...)*/
        printf("Slip systems | location  | burgID | planeID  | Schmid Factor | percentage\n"); 
        for(i=0; i<nSlips; i++){
            printf("    %3d      |   %3s     |   %3d  |   %3d    |    %5.3f      |  %5.2f  \n",
                i+1, (slip_sys[i][0]>=6)?"TB":((slip_sys[i][1]==3)?"M/T":(slip_sys[i][2]==0?"M":"T")), slip_sys[i][0], slip_sys[i][1], proba[i], percen[i]*100/totSrcLen);
        }

        free(proba);
        free(Proba);
        free(percen);
        free(primary);
        printf("%d chains has been distributed in matrix/twin\n"
               "Averaged number of chains over grain is %d\n" 
               "%4.2f%% of srcs are of primary slip systems\n",
                chain, chain/home->nGrains, primarySrcLen*100/totSrcLen);

#ifdef SAS_ANALYZE
        fclose(fp);
#endif
        
	    IncDislocationDensity(inData, totDislocLen);
	    param->nodeCount = inData->nodeCount;
	    WriteInitialNodeData(home, inData, 1);
	    FreeInNodeArray(inData, inData->nodeCount);
	    inData->nodeCount = 0;
	
        return;
}
