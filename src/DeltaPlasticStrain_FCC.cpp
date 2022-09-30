/***************************************************************************
 *   
 *      Author:       Sylvie Aubry (based on DeltaPlasticStrain_BCC.c)
 *
 *      Function:     DeltaPlasticStrain_FCC
 *
 *      Description:  Calculate the plastic strain increment
 *                    by dislocation motion in FCC materials.
 *
 *      Public functions:
 *          DeltaPlasticStrain_FCC()
 *
 *      Private functions:
 *          InitDeltaPlasticStrain_FCC()
 *
 ***************************************************************************/
#include "Home.h"
#include "Node.h"
#include "Util.h"
#include "Mobility.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*
 *      We only deal with specific burgers vectors in this function.
 *      Define the number of main burgers vectors we handle and the
 *      number of planes for each burgers vector.
 */
#define NUM_BURG    6
#define NUM_PLANES  3

#ifdef VELOCITY_STATISTICS
void GetBinID(real8 a0, real8 q, int n, real8 x, real8 w, int binID[2], real8 coeff[2])
{
    binID[0] = 0;
    binID[1] = 0;
    coeff[0] = 0.0;
    coeff[1] = 0.0;
#if 1
    real8 upper, invLower;
    if(x < a0*pow(q,0.5)){
        binID[0] = 0;
//        coeff[0] = x*w/(a0*pow(q,0.5));
        coeff[0] = w;
    }else if(x >= a0*pow(q,((real8)n))-0.5){
        binID[0] = n-1;
        coeff[0] = w;
//        coeff[0] = x*w/(a0*pow(q,((real8)n)-0.5));
    }else{
        for(int i=0; i<n-1; i++){
            if(x >= a0*pow(q,((real8)i)+0.5)){
                binID[0] = i;
                binID[1] = i+1;
                upper = -a0*pow(q,1.5+((real8)i))+x;
                invLower = 1.0/(a0*(pow(q,((real8)i)+0.5)-pow(q,((real8)i)+1.5))); 
                coeff[0] = w*upper*invLower;
                if(isnan(coeff[0]))coeff[0] = 0.0;
                coeff[1] = w*(1.0 - coeff[0]);
            }
        }
    }

    return;
#else
    if(x <= a0){
        coeff[0] = w;
        return;
    }

    for(int i=0; i<n; i++){
        if(x >= a0*pow(q, (real8)i) && x < a0*pow(q, ((real8)i) + 1.0)){
            binID[0] = i;
            coeff[0] = w;
        }
    }
#endif

    return;
}
#endif
/***************************************************************************
 *
 *      Function:    InitDeltaPlasticStrain_FCC
 *
 *      Description: Several of the arrays used during the delta plastic
 *                   strain calculations are dynamically calculated, but
 *                   once set, never need to be changed again for the
 *                   duration of the simulation.  This function should be
 *                   invoked during the first call to DeltaPlasticStrain_FCC()
 *                   to initialize those arrays, and never used again
 *                   after that.
 *
 ***************************************************************************/
static void InitDeltaPlasticStrain_FCC(int useLabFrame, real8 rotMatrix[3][3],
                                       real8 nanv[NUM_PLANES][3][NUM_BURG],
                                       real8 tanv[NUM_PLANES][3][NUM_BURG],
                                       real8 burgv[NUM_BURG][3],
                                       real8 bnorm[NUM_BURG][3])
{
	
        int bIndex;


        nanv[0][0][0]=0.5773502692; nanv[0][1][0]=-0.5773502692; nanv[0][2][0]=0.5773502692;
        nanv[1][0][0]=-0.5773502692; nanv[1][1][0]=0.5773502692; nanv[1][2][0]=0.5773502692;
        nanv[2][0][0]=0.0000000000; nanv[2][1][0]=0.0000000000; nanv[2][2][0]=1.0000000000;

        nanv[0][0][1]=0.5773502692; nanv[0][1][1]=0.5773502692; nanv[0][2][1]=0.5773502692;
        nanv[1][0][1]=0.5773502692; nanv[1][1][1]=0.5773502692; nanv[1][2][1]=-0.5773502692;
        nanv[2][0][1]=0.0000000000; nanv[2][1][1]=0.0000000000; nanv[2][2][1]=1.0000000000;

        nanv[0][0][2]=0.5773502692; nanv[0][1][2]=0.5773502692; nanv[0][2][2]=-0.5773502692;
        nanv[1][0][2]=-0.5773502692; nanv[1][1][2]=0.5773502692; nanv[1][2][2]=0.5773502692;
        nanv[2][0][2]=0.0000000000; nanv[2][1][2]=1.0000000000; nanv[2][2][2]=0.0000000000;

        nanv[0][0][3]=0.5773502692; nanv[0][1][3]=0.5773502692; nanv[0][2][3]=0.5773502692;
        nanv[1][0][3]=0.5773502692; nanv[1][1][3]=-0.5773502692; nanv[1][2][3]=0.5773502692;
        nanv[2][0][3]=0.0000000000; nanv[2][1][3]=1.0000000000; nanv[2][2][3]=0.0000000000;

        nanv[0][0][4]=0.5773502692; nanv[0][1][4]=0.5773502692; nanv[0][2][4]=-0.5773502692;
        nanv[1][0][4]=0.5773502692; nanv[1][1][4]=-0.5773502692; nanv[1][2][4]=0.5773502692;
        nanv[2][0][4]=1.0000000000; nanv[2][1][4]=0.0000000000; nanv[2][2][4]=0.0000000000;

        nanv[0][0][5]=0.5773502692; nanv[0][1][5]=0.5773502692; nanv[0][2][5]=0.5773502692;
        nanv[1][0][5]=-0.5773502692; nanv[1][1][5]=0.5773502692; nanv[1][2][5]=0.5773502692;
        nanv[2][0][5]=1.0000000000; nanv[2][1][5]=0.0000000000; nanv[2][2][5]=0.0000000000;

        tanv[0][0][0]=-0.4082482905; tanv[0][1][0]=0.4082482905; tanv[0][2][0]=0.8164965809;
        tanv[1][0][0]=-0.4082482905; tanv[1][1][0]=0.4082482905; tanv[1][2][0]=-0.8164965809;
        tanv[2][0][0]=-0.7071067812; tanv[2][1][0]=0.7071067812; tanv[2][2][0]=0.0000000000;

        tanv[0][0][1]=0.4082482905; tanv[0][1][1]=0.4082482905; tanv[0][2][1]=-0.8164965809;
        tanv[1][0][1]=-0.4082482905; tanv[1][1][1]=-0.4082482905; tanv[1][2][1]=-0.8164965809;
        tanv[2][0][1]=0.7071067812; tanv[2][1][1]=0.7071067812; tanv[2][2][1]=-0.0000000000;

        tanv[0][0][2]=0.4082482905; tanv[0][1][2]=-0.8164965809; tanv[0][2][2]=-0.4082482905;
        tanv[1][0][2]=0.4082482905; tanv[1][1][2]=0.8164965809; tanv[1][2][2]=-0.4082482905;
        tanv[2][0][2]=0.7071067812; tanv[2][1][2]=0.0000000000; tanv[2][2][2]=-0.7071067812;

        tanv[0][0][3]=-0.4082482905; tanv[0][1][3]=0.8164965809; tanv[0][2][3]=-0.4082482905;
        tanv[1][0][3]=0.4082482905; tanv[1][1][3]=0.8164965809; tanv[1][2][3]=0.4082482905;
        tanv[2][0][3]=-0.7071067812; tanv[2][1][3]=0.0000000000; tanv[2][2][3]=-0.7071067812;

        tanv[0][0][4]=0.8164965809; tanv[0][1][4]=-0.4082482905; tanv[0][2][4]=0.4082482905;
        tanv[1][0][4]=-0.8164965809; tanv[1][1][4]=-0.4082482905; tanv[1][2][4]=0.4082482905;
        tanv[2][0][4]=0.0000000000; tanv[2][1][4]=-0.7071067812; tanv[2][2][4]=0.7071067812;

        tanv[0][0][5]=-0.8164965809; tanv[0][1][5]=0.4082482905; tanv[0][2][5]=0.4082482905;
        tanv[1][0][5]=-0.8164965809; tanv[1][1][5]=-0.4082482905; tanv[1][2][5]=-0.4082482905;
        tanv[2][0][5]=-0.0000000000; tanv[2][1][5]=0.7071067812; tanv[2][2][5]=0.7071067812;

        bnorm[0][0]=0.7071067812; bnorm[0][1]=0.7071067812; bnorm[0][2]=0.0000000000;
        bnorm[1][0]=0.7071067812; bnorm[1][1]=-0.7071067812; bnorm[1][2]=0.0000000000;
        bnorm[2][0]=0.7071067812; bnorm[2][1]=0.0000000000; bnorm[2][2]=0.7071067812;
        bnorm[3][0]=0.7071067812; bnorm[3][1]=0.0000000000; bnorm[3][2]=-0.7071067812;
        bnorm[4][0]=0.0000000000; bnorm[4][1]=0.7071067812; bnorm[4][2]=0.7071067812;
        bnorm[5][0]=0.0000000000; bnorm[5][1]=0.7071067812; bnorm[5][2]=-0.7071067812;

        burgv[0][0]=0.7071067812; burgv[0][1]=0.7071067812; burgv[0][2]=0.0000000000;
        burgv[1][0]=0.7071067812; burgv[1][1]=-0.7071067812; burgv[1][2]=0.0000000000;
        burgv[2][0]=0.7071067812; burgv[2][1]=0.0000000000; burgv[2][2]=0.7071067812;
        burgv[3][0]=0.7071067812; burgv[3][1]=0.0000000000; burgv[3][2]=-0.7071067812;
        burgv[4][0]=0.0000000000; burgv[4][1]=0.7071067812; burgv[4][2]=0.7071067812;
        burgv[5][0]=0.0000000000; burgv[5][1]=0.7071067812; burgv[5][2]=-0.7071067812;

/*
 *      If necessary, rotate the nanv and tanv arrays from the crystal frame
 *      to the lab frame
 */
        if (useLabFrame) {
            int 		plane;

            for (plane = 0; plane < NUM_PLANES; plane++) {

                for (bIndex = 0; bIndex < NUM_BURG; bIndex++) {
                    real8 nanvCryst[3], nanvLab[3], tanvCryst[3], tanvLab[3];

                    nanvCryst[0] = nanv[plane][0][bIndex];
                    nanvCryst[1] = nanv[plane][1][bIndex];
                    nanvCryst[2] = nanv[plane][2][bIndex];

                    tanvCryst[0] = tanv[plane][0][bIndex];
                    tanvCryst[1] = tanv[plane][1][bIndex];
                    tanvCryst[2] = tanv[plane][2][bIndex];

                    Matrix33Vector3Multiply(rotMatrix,nanvCryst,nanvLab);
                    Matrix33Vector3Multiply(rotMatrix,tanvCryst,tanvLab);

                    nanv[plane][0][bIndex] = nanvLab[0];
                    nanv[plane][1][bIndex] = nanvLab[1];
                    nanv[plane][2][bIndex] = nanvLab[2];

                    tanv[plane][0][bIndex] = tanvLab[0];
                    tanv[plane][1][bIndex] = tanvLab[1];
                    tanv[plane][2][bIndex] = tanvLab[2];
                }
            }
        }

/*
 *      If necessary, rotate the bnorm array from the crystal frame to
 *      the lab frame
 */
        if (useLabFrame) {

            for (bIndex = 0; bIndex < NUM_BURG; bIndex++) {
                real8 bnormLab[3];

                Matrix33Vector3Multiply(rotMatrix,bnorm[bIndex],bnormLab);
                VECTOR_COPY(bnorm[bIndex], bnormLab);
            }
        }

#if 0
/*
 *      Some stuff for debugging
 */
        for (i = 0; i < NUM_BURG; i++) {
            printf("bnorm[%d]=%f %f %f\n", i,
                   bnorm[i][0], bnorm[i][1], bnorm[i][2]);
        }
        printf("\n");

        for (i = 0; i < NUM_BURG; i++) {
            for (j = 0; j < NUM_PLANES; j++) {
                printf("nanv[%d][%d]=%f %f %f\n", i, j,
                       nanv[j][0][i], nanv[j][1][i], nanv[j][2][i]);
            }
            printf("\n");
        }

        for (i = 0; i < NUM_BURG; i++) {
            for (j = 0; j < NUM_PLANES; j++) {
                printf("tanv[%d][%d]=%f %f %f\n", i, j,
                       tanv[j][0][i], tanv[j][1][i], tanv[j][2][i]);
            }
            printf("\n");
        }

#endif
        return;
}


void DeltaPlasticStrain_FCC(Home_t *home)
{

        int          i, j;
        real8        factor0, bmagsq;
        real8        gdspn[6], gdstn[6], gdtstn[6];
        real8        dstn[3][3], dspn[3][3], dtstn[3][3];
        real8        Ltot[NUM_BURG][4], areaSwept[NUM_BURG][7];
        real8        rotMatrix[3][3];
        Param_t      *param;
        static int   initDone = 0;
        static real8 burgv[NUM_BURG][3], bnorm[NUM_BURG][3];
        static real8 nanv[NUM_PLANES][3][NUM_BURG];
        static real8 tanv[NUM_PLANES][3][NUM_BURG];
        
#ifdef VELOCITY_STATISTICS
        real8   velocity, coeff[2];
        real8   delAccumulation[8], disLength[3];
        delAccumulation[0] = 0.0;
        delAccumulation[1] = 0.0;
        delAccumulation[2] = 0.0;
        delAccumulation[3] = 0.0;
        delAccumulation[4] = 0.0;
        delAccumulation[5] = 0.0;
        delAccumulation[6] = 0.0;
        delAccumulation[7] = 0.0;
        disLength[0] = 0.0;
        disLength[1] = 0.0;
        disLength[2] = 0.0;
#endif
/*
 *      Ltot : density of 
 *         [burg][0] : screw
 *         [burg][1] : edge 1
 *         [burg][2] : edge_2
 *         [burg][3] : edge_3
 *      for each burgers vector.
 */
        VECTOR_ZERO(rotMatrix[0]);
        VECTOR_ZERO(rotMatrix[1]);
        VECTOR_ZERO(rotMatrix[2]);

        rotMatrix[0][0] = 1.0;
        rotMatrix[1][1] = 1.0;
        rotMatrix[2][2] = 1.0;

        param = home->param;

        real8   extstress[3][3];

        extstress[0][0] = param->appliedStress[0];
        extstress[1][1] = param->appliedStress[1];
        extstress[2][2] = param->appliedStress[2];
        extstress[1][2] = param->appliedStress[3];
        extstress[2][0] = param->appliedStress[4];
        extstress[0][1] = param->appliedStress[5];
        extstress[2][1] = extstress[1][2];
        extstress[0][2] = extstress[2][0];
        extstress[1][0] = extstress[0][1];

        bmagsq = param->burgMag * param->burgMag;
        param->disloDensity = 0.0;
        param->disloDensity_TB = 0.0;

        memset(Ltot, 0, sizeof(Ltot));
        memset(areaSwept, 0, sizeof(areaSwept));

        Matrix33Init(dstn);
        Matrix33Init(dspn);
        Matrix33Init(dtstn);


        for (i = 0; i < param->numBurgGroups; i++) {
            param->partialDisloDensity[i] = 0.0;
        }

        for (i=0; i < 6; i++) {
            param->delpStrain[i] = 0.0;                                         
            param->delpTwinStrain[i] = 0.0;
        }                                                                


/*
 *      The first time we enter this function, initialize a few
 *      arrays that will remain static for the duration of the
 *      simulation once they are set.
 */
		
        if (!param->useLabFrame && !initDone++) {
            InitDeltaPlasticStrain_FCC(0, rotMatrix, nanv, tanv, burgv, bnorm);
        }
	
/*
 *      Loop through all the nodes calculating the flux contribution
 *      from each of the node's segments.  We want the loop threaded,
 *      but we also need to sum thread-specific totals when we're
 *      done, so we can't just automatically distribute the loop
 *      iterations.  Instead, we'll set up a threaded execution block
 *      where each thread calculates it's own loop indices and then
 *      safely adds its contributions to the summed totals in a 'critical'
 *      section after its loop iterations have been completed.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
   {
        int    k, l;
        int    nodeIndex;
        int    threadID, threadIterStart, threadIterEnd;
        real8  thread_disloDensity=0.0, thread_disloDensity_TB=0.0;
        real8  *thread_partialDensity;
        real8  thread_dstn[3][3], thread_dspn[3][3], thread_dyad[3][3];
        real8  thread_dtstn[3][3];
        real8  thread_Ltot[NUM_BURG][4], thread_areaSwept[NUM_BURG][7];

#ifdef VELOCITY_STATISTICS
        real8  thread_delAccumulation[8], thread_disLength[3], f1[3], f2[3];

        thread_disloDensity = 0.0;
        thread_disloDensity_TB = 0.0;
        thread_delAccumulation[0] = 0.0;
        thread_delAccumulation[1] = 0.0;
        thread_delAccumulation[2] = 0.0;
        thread_delAccumulation[3] = 0.0;
        thread_delAccumulation[4] = 0.0;
        thread_delAccumulation[5] = 0.0;
        thread_delAccumulation[6] = 0.0;
        thread_delAccumulation[7] = 0.0;
        thread_disLength[0] = 0.0;
        thread_disLength[1] = 0.0;
        thread_disLength[2] = 0.0;
#endif
        thread_partialDensity = (real8 *)calloc(1, param->numBurgGroups *
                                               sizeof(real8));

        Matrix33Init(thread_dstn);
        Matrix33Init(thread_dspn);
        Matrix33Init(thread_dtstn);

        memset(thread_Ltot, 0, sizeof(thread_Ltot));
        memset(thread_areaSwept, 0, sizeof(thread_areaSwept));

        GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                  &threadIterStart, &threadIterEnd);

/*
 *      Each thread now loops over its portion of the node list
 */
        for (nodeIndex=threadIterStart; nodeIndex<threadIterEnd; nodeIndex++) {
            int    m;
            int    numSeg, segIndex;
            real8  localstrain[3][3];
            Node_t *node;

            if ((node = home->nodeKeys[nodeIndex]) == (Node_t *)NULL) {
                continue;
            }

            numSeg = node->numNbrs;

            Matrix33Init(localstrain);
        
/*
 *          Loop over every segment attached to the node
 */        
            for (segIndex = 0; segIndex < numSeg; segIndex++) {
               real8  bx, by, bz;
               real8  ex, ey, ez;
               real8  nx, ny, nz;
               real8  hhx, hhy, hhz;
               real8  tmpx, tmpy, tmpz;
               real8  deltaDensity;
               real8  size;
               real8  bCryst[3];
               real8  zeta[3];
			   real8  pp[3],p1[3],p2[3];
               real8  tmpDstn;
               Node_t *nbr;

               int  burgType;
               int  meshNx, meshNy, meshNz, eleID[27];
               real8    lenPer[27], cen1[3], cen2[3];

#ifdef VELOCITY_STATISTICS
               real8  disLen, burgMag, vtau;
               int  binID[2];
               real8  vec1[3], vec2[3], lDir[3], nDir[3];
               real8 dot, gDir[3], burg[3], normal[3]; 
               int  nbrArmID;
               bool active, positive;
#endif
               nbr = GetNeighborNode(home, node, segIndex);

               if (nbr == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               /* Avoid double counting */

               if (OrderNodes(node, nbr) != -1) {
                   continue;
               }
			   
			   if(param->useLabFrame){
                    if(param->polycrystal){

                        int     inGrain, onGB;
                        Grain_t *grain;
                        GB_t    *GB;

                        GetArmGrain(home, node, nbr, inGrain, onGB);
                        if(onGB > -1){
                            GB = home->GBKeys[onGB];
                            grain = home->grainKeys[GB->nbrGrain[0]];
                        }else{
                            grain = home->grainKeys[inGrain];    
                        }
                        
                        memcpy(rotMatrix, grain->orientation, sizeof(int [3][3]));
                    }else{
                        memcpy(rotMatrix, home->rotMatrix, sizeof(int [3][3]));
                    }

					InitDeltaPlasticStrain_FCC(1, rotMatrix, nanv, tanv, burgv, bnorm);
			   }
            
#ifdef VELOCITY_STATISTICS
                lDir[0] = nbr->x - node->x;
                lDir[1] = nbr->y - node->y;
                lDir[2] = nbr->z - node->z;
                ZImage(param, lDir, lDir+1, lDir+2);

                ExtPKForce(extstress, node->burgX[segIndex], node->burgY[segIndex], node->burgZ[segIndex],
                            node->x, node->y, node->z, node->x+lDir[0], node->y+lDir[1], node->z+lDir[2], f1, f2);

                disLen = Normal(lDir)/param->maxSeg; // units maxSegLen;
                NormalizeVec(lDir);

                vec1[0] = (node->vX + nbr->vX)*0.5*param->burgMag; 
                vec1[1] = (node->vY + nbr->vY)*0.5*param->burgMag; 
                vec1[2] = (node->vZ + nbr->vZ)*0.5*param->burgMag; 
                velocity = Normal(vec1);
                if(isnan(velocity)) {
                    velocity = 0.0;
                    vec1[0] = 0.0;
                    vec1[1] = 0.0;
                    vec1[2] = 0.0;
                }

                nbrArmID = GetArmID(home, nbr, node);
                vec2[0] = node->armfx[segIndex] + nbr->armfx[nbrArmID];
                vec2[1] = node->armfy[segIndex] + nbr->armfy[nbrArmID];
                vec2[2] = node->armfz[segIndex] + nbr->armfz[nbrArmID];
                burgMag = sqrt(node->burgX[segIndex]*node->burgX[segIndex]+
                               node->burgY[segIndex]*node->burgY[segIndex]+
                               node->burgZ[segIndex]*node->burgZ[segIndex]);

                nDir[0] = node->nx[segIndex];
                nDir[1] = node->ny[segIndex];
                nDir[2] = node->nz[segIndex];
                dot = DotProduct(nDir, vec2);
                
                vec2[0] -= (dot*nDir[0]); 
                vec2[1] -= (dot*nDir[1]); 
                vec2[2] -= (dot*nDir[2]); 

                cross(lDir, nDir, gDir);
                
                NormalizeVec(gDir);
                dot = DotProduct(vec2, gDir);
                vec2[0] = dot*gDir[0];
                vec2[1] = dot*gDir[1];
                vec2[2] = dot*gDir[2];
                
                if(burgMag <1.0E-10){
                    Fatal("somthing wrong in file %s at %d\n", __FILE__, __LINE__);
                }

                vtau = ((disLen < 1.0E-5) ? 0.0 : fabs(DotProduct(vec1, vec2)/1E6/param->maxSeg)); // units MPa.m/s

                if((GetPlaneType_FCC(home, node, nbr) & GLIDE_PLANE) > 0 && fabs(burgMag - 1.0) < 1.0E-2){
                    real8 zaxis[3] = {0.0, 0.0, 1.0};
                    
                    burg[0] = node->burgX[segIndex];
                    burg[1] = node->burgY[segIndex];
                    burg[2] = node->burgZ[segIndex];
                    normal[0] = node->nx[segIndex];
                    normal[1] = node->ny[segIndex];
                    normal[2] = node->nz[segIndex];
                    NormalizeVec(burg);
                    NormalizeVec(normal);
                    dot = fabs(DotProduct(burg, zaxis)*DotProduct(normal, zaxis));
                    
                    if(fabs(dot-param->maxSchmid) < 1.0E-2){
                        active = 1;
                    }else active = 0;
                }else{
                    active = 0;
                }

                dot = 2*DotProduct(f1, nDir);
                f1[0] -= (dot*nDir[0]); 
                f1[1] -= (dot*nDir[1]); 
                f1[2] -= (dot*nDir[2]); 
                dot = DotProduct(f1, vec1);
                positive = ((dot >= 0) ? 1 : 0);           

//                    printf("%d %d %d: %f %f %f\n", (int)active, BurgID(home, node, nbr),
//                            PlaneID(home, node, nbr), Normal(vec1),Normal(vec2), vtau);
                if(active){
                    thread_delAccumulation[2] += (DotProduct(vec1, vec1)*disLen);
                    thread_delAccumulation[6] += vtau;
                    if(!positive){
                        thread_disLength[0] += disLen;
                        thread_delAccumulation[0] += (DotProduct(vec1, vec1)*disLen);
                        thread_delAccumulation[4] += vtau;
                        GetBinID(param->velMin, param->baseVel, param->nVelBins, velocity, disLen, binID, coeff);
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                        {
                            param->velTable[binID[0]][1] += coeff[0];
                            param->velTable[binID[1]][1] += coeff[1];
                        }
                        
                        GetBinID(param->dissMin, param->baseDiss, param->nDissBins, vtau, disLen, binID, coeff);
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                        {
                            param->dissTable[binID[0]][1] += coeff[0];
                            param->dissTable[binID[1]][1] += coeff[1];
                        }
                    }else{
                        thread_disLength[1] += disLen;
                        thread_delAccumulation[1] += (DotProduct(vec1, vec1)*disLen);
                        thread_delAccumulation[5] += vtau;
                        GetBinID(param->velMin, param->baseVel, param->nVelBins, velocity, disLen, binID, coeff);
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                        {
                            param->velTable[binID[0]][2] += coeff[0];
                            param->velTable[binID[1]][2] += coeff[1];
                        }

                        GetBinID(param->dissMin, param->baseDiss, param->nDissBins, vtau, disLen, binID, coeff);
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                        {
                            param->dissTable[binID[0]][2] += coeff[0];
                            param->dissTable[binID[1]][2] += coeff[1];
                        }
                    }


                }else{
                    thread_disLength[2] += disLen;
                    thread_delAccumulation[3] += (DotProduct(vec1, vec1)*disLen);
                    thread_delAccumulation[7] += vtau;

                    GetBinID(param->velMin, param->baseVel, param->nVelBins, velocity, disLen, binID, coeff);
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                    {
                        param->velTable[binID[0]][3] += coeff[0];
                        param->velTable[binID[1]][3] += coeff[1];
                    }

                    GetBinID(param->dissMin, param->baseDiss, param->nDissBins, vtau, disLen, binID, coeff);
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                    {
                        param->dissTable[binID[0]][3] += coeff[0];
                        param->dissTable[binID[1]][3] += coeff[1];
                    }
                }

#endif
/*
 *             This check is in support of some local site-specific
 *             post-processing only and all *normal* ParaDiS simulations
 *             will simply continue on through this check.
 */
  /*             if (param->ignoreNonLoopFlux) {
                   if (!HAS_ANY_OF_CONSTRAINTS(node->constraint, LOOP_NODE) ||
                       !HAS_ANY_OF_CONSTRAINTS(nbr->constraint,  LOOP_NODE)) {
                       continue;
                   }
               }
 */
/*
 *             For later, we need a copy of the burgers vector guaranteed
 *             to be in the crystal frame
 */
               bx = node->burgX[segIndex];
               by = node->burgY[segIndex];
               bz = node->burgZ[segIndex];

               if (param->useLabFrame) {
                   real8 bLab[3] = {bx, by, bz};
                   VectorLabToCrystal(home, node, nbr, bLab, bCryst);
               } else {
                   bCryst[0] = bx;
                   bCryst[1] = by;
                   bCryst[2] = bz;
               }


                burgType = BurgType_FCC(bCryst);

/*
 * 		In the output, we would like to generate the plastic strain
 * 		field, so we should record all the plastic strain in each 
 * 		element. Here, we get the eleNum of the current area position. 
 */
	        p1[0] = node->x;
	        p1[1] = node->y;
	        p1[2] = node->z;
	         
	        p2[0] = node->oldx;
	        p2[1] = node->oldy;
	        p2[2] = node->oldz;
	         
	        PBCPOSITION(param,p1[0],p1[1],p1[2],&p2[0],&p2[1],&p2[2]);
	         
	        pp[0] = p1[0] + p2[0];
	        pp[1] = p1[1] + p2[1];
	        pp[2] = p1[2] + p2[2];

            cen1[0] = 0.5*pp[0];
            cen1[1] = 0.5*pp[1];
            cen1[2] = 0.5*pp[2];
	         
	        p2[0] = nbr->oldx;
	        p2[1] = nbr->oldy;
	        p2[2] = nbr->oldz;
	         
	        PBCPOSITION(param,p1[0],p1[1],p1[2],&p2[0],&p2[1],&p2[2]);

            VECTOR_COPY(cen2, p2);
	         
	        pp[0] +=  p2[0];
	        pp[1] +=  p2[1];
	        pp[2] +=  p2[2];
	         
	        p2[0] = nbr->x;
	        p2[1] = nbr->y;
	        p2[2] = nbr->z;
	         
	        PBCPOSITION(param,p1[0],p1[1],p1[2],&p2[0],&p2[1],&p2[2]);

            cen2[0] = 0.5*(cen2[0]+p2[0]);
            cen2[1] = 0.5*(cen2[1]+p2[1]);
            cen2[2] = 0.5*(cen2[2]+p2[2]);

	        pp[0] +=  p2[0];
	        pp[1] +=  p2[1];
	        pp[2] +=  p2[2];
			   
	        VECTOR_MULTIPLY(pp,0.25);

/*
 *          Some bugs can cause trouble here, such as the nodal position
 *          is nan.                                             D. Wei
 */
            if(isnan(pp[0]) || isnan(pp[1]) || isnan(pp[2])){
                PrintNode(node);
                PrintNode(nbr);
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
	         
	        Coor2EleNum(home,pp,meshNx,meshNy,meshNz);
	        FoldBox(param,&pp[0],&pp[1],&pp[2]);
            DecomposeSegmentInRubik(home, meshNx,meshNy,meshNz, cen1, cen2, lenPer, eleID);
			   
	        ex = nbr->x - node->oldx; 
	        ey = nbr->y - node->oldy; 
	        ez = nbr->z - node->oldz; 

	        hhx = node->x - nbr->oldx;
	        hhy = node->y - nbr->oldy;
	        hhz = node->z - nbr->oldz;

/*
 *          Some bugs can cause trouble here, such as the nodal old 
 *          psoition is nan.                                  D. Wei
 */
            if(isnan(ex) || isnan(ey) || isnan(ez) ||
               isnan(hhx) || isnan(hhy) || isnan(hhz)){
                PrintNode(node);
                PrintNode(nbr);
                Fatal("in %s at %d", __FILE__, __LINE__);
            }

	        ZImage(param, &ex, &ey, &ez);
	        ZImage(param, &hhx, &hhy, &hhz);
	        
	        nx = (ey*hhz - ez*hhy) * 0.5;
	        ny = (ez*hhx - ex*hhz) * 0.5;
	        nz = (ex*hhy - ey*hhx) * 0.5;
	         
	        thread_dyad[0][0] = nx*bx;
	        thread_dyad[0][1] = nx*by;
	        thread_dyad[0][2] = nx*bz;
	        thread_dyad[1][0] = ny*bx;
	        thread_dyad[1][1] = ny*by;
	        thread_dyad[1][2] = ny*bz;
	        thread_dyad[2][0] = nz*bx;
	        thread_dyad[2][1] = nz*by;
	        thread_dyad[2][2] = nz*bz;

/*
 *          Caculate total plastci strain, twinning deformation and elemental 
 *          plastic strain.                                         D. Wei
 */
            for (k = 0; k < 3; k++){
                for (l = 0; l < 3; l++){
                    tmpDstn = (thread_dyad[l][k]+thread_dyad[k][l])*0.5;
                    thread_dstn[l][k] += tmpDstn;

                    if(burgType == SHOCKLEY_BURG && (node->onGB>-1 && nbr->onGB>-1)){
                        thread_dtstn[l][k] += tmpDstn;
                    }
			 	  
	                if(l==k){
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                        {
                            for(m=0;m<27;m++)home->pstnEle[eleID[m]][l] += tmpDstn*lenPer[m]; // element
                        }
                    }else if(l<k){
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
                        {
                            for(m=0;m<27;m++)home->pstnEle[eleID[m]][6-l-k] += tmpDstn*lenPer[m]; // element
                        }
                    }

                    thread_dspn[l][k] += (thread_dyad[l][k] - thread_dyad[k][l])*0.5;
                    localstrain[l][k] += tmpDstn;
                }
            }

            tmpx = nbr->x - node->x;
            tmpy = nbr->y - node->y;
            tmpz = nbr->z - node->z;

            ZImage(param, &tmpx, &tmpy, &tmpz);

            zeta[0] = tmpx;
            zeta[1] = tmpy;
            zeta[2] = tmpz;

            size=sqrt(zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2]);

/*
 *          deltaDensity = size/param->simVol/bmagsq;
 */
            deltaDensity = size * param->burgVolFactor;
            thread_disloDensity += deltaDensity;
            if(burgType == SHOCKLEY_BURG)thread_disloDensity_TB+=deltaDensity;

	    }  /* for (segIndex = 0; segIndex < nc; ...) */
	   

/**********************************************************************/

/*
 *          computed localstrain[3][3] for the node
 */
            if ((param->loadType == 1)) {
                real8  localeRate;
                localeRate =
                   (param->edotdir[0]*localstrain[0][0]*param->edotdir[0] +
                    param->edotdir[0]*localstrain[0][1]*param->edotdir[1] +
                    param->edotdir[0]*localstrain[0][2]*param->edotdir[2] +
                    param->edotdir[1]*localstrain[1][0]*param->edotdir[0] +
                    param->edotdir[1]*localstrain[1][1]*param->edotdir[1] +
                    param->edotdir[1]*localstrain[1][2]*param->edotdir[2] +
                    param->edotdir[2]*localstrain[2][0]*param->edotdir[0] +
                    param->edotdir[2]*localstrain[2][1]*param->edotdir[1] +
                    param->edotdir[2]*localstrain[2][2]*param->edotdir[2])/param->simVol;

                if (localeRate > 0.0) {
                    node->sgnv = 1;
                } else if (localeRate<0.0) {
                    node->sgnv = -1;
                } else {
                    node->sgnv = 0;
                }
            }
        }  /* for (nodeIndex = 0; nodeIndex < home->newNodeKeyPtr; ...) */

/*
 *      At this point, each thread has a portion of several arrays.
 *      The threads now update the global arrays with their contributions
 *      in a 'critical' block so there's no conflicts among the threads.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_DPS_FCC_A)
#endif
       {
            for (k = 0; k < 3; k++) {
                for (l = 0; l < 3; l++) {
                    dstn[k][l] += thread_dstn[k][l];
                    dspn[k][l] += thread_dspn[k][l];
                    dtstn[k][l] += thread_dtstn[k][l];
                }
            }
#ifdef VELOCITY_STATISTICS
            disLength[0] += thread_disLength[0];
            disLength[1] += thread_disLength[1];
            disLength[2] += thread_disLength[2];
            delAccumulation[0] += thread_delAccumulation[0];
            delAccumulation[1] += thread_delAccumulation[1];
            delAccumulation[2] += thread_delAccumulation[2];
            delAccumulation[3] += thread_delAccumulation[3];
            delAccumulation[4] += thread_delAccumulation[4];
            delAccumulation[5] += thread_delAccumulation[5];
            delAccumulation[6] += thread_delAccumulation[6];
            delAccumulation[7] += thread_delAccumulation[7];
#endif
#if 0
            for (k = 0; k < NUM_BURG; k++) {
                for (l = 0; l < 4; l++) {
                    Ltot[k][l] += thread_Ltot[k][l];
                    areaSwept[k][l] += thread_areaSwept[k][l];
                }
                for (l = 4; l < 7; l++) {
                    areaSwept[k][l] += thread_areaSwept[k][l];
                }
            }
#endif
            for (k = 0; k < param->numBurgGroups; k++) {
                param->partialDisloDensity[k] += thread_partialDensity[k];
            }

            param->disloDensity += thread_disloDensity;
            param->disloDensity_TB += thread_disloDensity_TB;

        }  /* end "omp critical" section */

        if (thread_partialDensity != NULL) {
            free(thread_partialDensity);
        }

    }  /* end "omp parallel" section */

#if 0
/*
 *      Some stuff for debugging
 */
        printf("Ltot ---------------- \n");
        for (i = 0; i < NUM_BURG; i++) {
            printf("%e %e %e %e \n", Ltot[i][0], Ltot[i][1],
                   Ltot[i][2], Ltot[i][3]);
        }
        printf("areaSwept printing ---- \n");
        for (i = 0; i < NUM_BURG; i++) {
            printf(" %e %e %e %e %e %e %e\n",
                   areaSwept[i][0], areaSwept[i][1],
                   areaSwept[i][2], areaSwept[i][3],
                   areaSwept[i][4], areaSwept[i][5],
                   areaSwept[i][6]);
        }
#endif

/*
 *      Density decomposition
 */
        factor0 = param->simVol*bmagsq;

        for (i = 0; i < NUM_BURG; i++) {
            for (j = 0; j < 4; j++) {
                param->FCC_dLtot[i][j] = Ltot[i][j] / factor0; 
            }
            for (j = 0; j < 7; j++) {
                param->FCC_dfluxtot[i][j] = areaSwept[i][j] / param->simVol / 
                                            param->realdt;
            }
        }
 
#if 0
        printf("Ltot ---------------- \n");
        for (i = 0; i < NUM_BURG; i++) {
            printf("%e %e %e %e \n", Ltot[i][0], Ltot[i][1],
                   Ltot[i][2], Ltot[i][3]);
        }
#endif


/*
 *      accumulate delta strain this subcycle into delta strain this cycle
 */
        param->delpStrain[0] = dstn[0][0]/param->simVol;
        param->delpStrain[1] = dstn[1][1]/param->simVol;
        param->delpStrain[2] = dstn[2][2]/param->simVol;
        param->delpStrain[3] = dstn[1][2]/param->simVol;
        param->delpStrain[4] = dstn[0][2]/param->simVol;
        param->delpStrain[5] = dstn[0][1]/param->simVol;

        param->delpTwinStrain[0] = dtstn[0][0]/param->simVol;
        param->delpTwinStrain[1] = dtstn[1][1]/param->simVol;
        param->delpTwinStrain[2] = dtstn[2][2]/param->simVol;
        param->delpTwinStrain[3] = dtstn[1][2]/param->simVol;
        param->delpTwinStrain[4] = dtstn[0][2]/param->simVol;
        param->delpTwinStrain[5] = dtstn[0][1]/param->simVol;
        
        param->delpSpin[0] = dspn[0][0]/param->simVol;
        param->delpSpin[1] = dspn[1][1]/param->simVol;
        param->delpSpin[2] = dspn[2][2]/param->simVol;
        param->delpSpin[3] = dspn[1][2]/param->simVol;
        param->delpSpin[4] = dspn[0][2]/param->simVol;
        param->delpSpin[5] = dspn[0][1]/param->simVol;

#ifdef VELOCITY_STATISTICS                                                                                
        param->delAccumulation[0] = delAccumulation[0]; 
        param->delAccumulation[1] = delAccumulation[1]; 
        param->delAccumulation[2] = delAccumulation[2]; 
        param->delAccumulation[3] = delAccumulation[3]; 
        param->delAccumulation[4] = delAccumulation[4]; 
        param->delAccumulation[5] = delAccumulation[5]; 
        param->delAccumulation[6] = delAccumulation[6]; 
        param->delAccumulation[7] = delAccumulation[7]; 
        param->disLength[0] = disLength[0];
        param->disLength[1] = disLength[1];
        param->disLength[2] = disLength[2];
#endif
#ifdef PARALLEL
/*
 *      We've calculated processor specific values, now sum the
 *      delta strain from all processors, and accumulate into net strain
 */
        MPI_Allreduce(param->delpStrain, gdstn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(param->delpSpin, gdspn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(param->delpTwinStrain, gdtstn, 6, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

#ifdef VELOCITY_STATISTICS                                                                                
        MPI_Allreduce(param->delAccumulation, delAccumulation, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(param->disLength, disLength, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        for (i = 0; i < 6; i++) {
            param->delpStrain[i] = gdstn[i];
            param->delpSpin[i] = gdspn[i];
            param->delpTwinStrain[i] = gdtstn[i];
        }

#ifdef VELOCITY_STATISTICS                                                                                
        param->disLength[0] = disLength[0];
        param->disLength[1] = disLength[1];
        param->disLength[2] = disLength[2];
        param->delAccumulation[0] = delAccumulation[0];
        param->delAccumulation[1] = delAccumulation[1];
        param->delAccumulation[2] = delAccumulation[2];
        param->delAccumulation[3] = delAccumulation[3];
        param->delAccumulation[4] = delAccumulation[4]; 
        param->delAccumulation[5] = delAccumulation[5]; 
        param->delAccumulation[6] = delAccumulation[6]; 
        param->delAccumulation[7] = delAccumulation[7]; 
#endif
/*
 *      Flux decomposition
 */
        MPI_Allreduce(param->FCC_dLtot,    param->FCC_Ltot,    NUM_BURG*4,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(param->FCC_dfluxtot, param->FCC_fluxtot, NUM_BURG*7,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
/*
 *      For serial compilation, no need to accumulate values from remote
 *      processors, just copy the data into the param struct.
 */
        for (i = 0; i < NUM_BURG; i++) {
            for (j = 0; j < 4; j++) {
                param->FCC_Ltot[i][j] = param->FCC_dLtot[i][j];
            }
            for (j = 0; j < 7; j++) {
                param->FCC_fluxtot[i][j] = param->FCC_dfluxtot[i][j];
            }
        } 

        for (i = 0; i < 6; i++) {
            gdstn[i] = param->delpStrain[i];
            gdspn[i] = param->delpSpin[i];
            gdtstn[i] = param->delpTwinStrain[i];
        }
#ifdef VELOCITY_STATISTICS                                                                                
        param->disLength[0] = disLength[0];
        param->disLength[1] = disLength[1];
        param->disLength[2] = disLength[2];
        param->delAccumulation[0] = delAccumulation[0];
        param->delAccumulation[1] = delAccumulation[1];
        param->delAccumulation[2] = delAccumulation[2];
        param->delAccumulation[3] = delAccumulation[3];
        param->delAccumulation[4] = delAccumulation[4]; 
        param->delAccumulation[5] = delAccumulation[5]; 
        param->delAccumulation[6] = delAccumulation[6]; 
        param->delAccumulation[7] = delAccumulation[7]; 
#endif
#endif /* if PARALLEL */


#if 0
/*
 *      For debugging only.  See comments further on
 */
        if (home->myDomain == 0) {
            int   k, l;
            real8 pstn[3][3];
            real8 val[3][3];
            real8 dyad[3][3], dyad0[3][3], dyad1[3][3], dyad2[3][3];

/*
 *          Strain increment using decomposition here
 */
            for (i = 0; i < NUM_BURG; i++) {

                /* climb */
                dyad[0][0] = bnorm[i][0]*bnorm[i][0];
                dyad[0][1] = bnorm[i][0]*bnorm[i][1];
                dyad[0][2] = bnorm[i][0]*bnorm[i][2];

                dyad[1][0] = bnorm[i][1]*bnorm[i][0];
                dyad[1][1] = bnorm[i][1]*bnorm[i][1];
                dyad[1][2] = bnorm[i][1]*bnorm[i][2];

                dyad[2][0] = bnorm[i][2]*bnorm[i][0];
                dyad[2][1] = bnorm[i][2]*bnorm[i][1];
                dyad[2][2] = bnorm[i][2]*bnorm[i][2];

                for (k = 0; k < 3; k++){
                    for (l = 0; l < 3; l++){
                        pstn[k][l] += param->FCC_fluxtot[i][0] *
                                      (dyad[l][k]+dyad[k][l]) * 0.5 /
                                      param->simVol;
                    }
                }
              
                /* Three planes */
                dyad0[0][0] = bnorm[i][0]*nanv[0][0][i];
                dyad0[0][1] = bnorm[i][0]*nanv[0][1][i];
                dyad0[0][2] = bnorm[i][0]*nanv[0][2][i]; 

                dyad0[1][0] = bnorm[i][1]*nanv[0][0][i];
                dyad0[1][1] = bnorm[i][1]*nanv[0][1][i];
                dyad0[1][2] = bnorm[i][1]*nanv[0][2][i]; 

                dyad0[2][0] = bnorm[i][2]*nanv[0][0][i];
                dyad0[2][1] = bnorm[i][2]*nanv[0][1][i];
                dyad0[2][2] = bnorm[i][2]*nanv[0][2][i]; 

                dyad1[0][0] = bnorm[i][0]*nanv[1][0][i];
                dyad1[0][1] = bnorm[i][0]*nanv[1][1][i];
                dyad1[0][2] = bnorm[i][0]*nanv[1][2][i]; 

                dyad1[1][0] = bnorm[i][1]*nanv[1][0][i];
                dyad1[1][1] = bnorm[i][1]*nanv[1][1][i];
                dyad1[1][2] = bnorm[i][1]*nanv[1][2][i]; 

                dyad1[2][0] = bnorm[i][2]*nanv[1][0][i];
                dyad1[2][1] = bnorm[i][2]*nanv[1][1][i];
                dyad1[2][2] = bnorm[i][2]*nanv[1][2][i]; 

                dyad2[0][0] = bnorm[i][0]*nanv[2][0][i];
                dyad2[0][1] = bnorm[i][0]*nanv[2][1][i];
                dyad2[0][2] = bnorm[i][0]*nanv[2][2][i]; 

                dyad2[1][0] = bnorm[i][1]*nanv[2][0][i];
                dyad2[1][1] = bnorm[i][1]*nanv[2][1][i];
                dyad2[1][2] = bnorm[i][1]*nanv[2][2][i]; 

                dyad2[2][0] = bnorm[i][2]*nanv[2][0][i];
                dyad2[2][1] = bnorm[i][2]*nanv[2][1][i];
                dyad2[2][2] = bnorm[i][2]*nanv[2][2][i]; 

                for (k = 0; k < 3; k++) {
                    for(l = 0; l < 3; l++) {
                  
                        pstn[k][l] += ((param->FCC_fluxtot[i][1]+param->FCC_fluxtot[i][4]) * (dyad0[l][k]+dyad0[k][l]) +
                                       (param->FCC_fluxtot[i][2]+param->FCC_fluxtot[i][5]) * (dyad1[l][k]+dyad1[k][l]) +
                                       (param->FCC_fluxtot[i][3]+param->FCC_fluxtot[i][6]) * (dyad2[l][k]+dyad2[k][l])) * 0.5 / param->simVol;
                    }
                }
            }  /* for (i=0;i<NUM_BURG; etc. */

/* 
 *          dstn2 tensor the total sum of all decomposed strain flux components
 *          The sum should be identical to delpStrain if the decomposition is
 *          done correctly .. To check if both methods produce the same
 *          result ..  Do not delete the print statements.
 */
#if 0
            printf("--simVoll = %e -----\n", param->simVol);
            printf("--realdt = %e -----\n", param->realdt);

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    val[i][j] = pstn[i][j]*(param->simVol * param->realdt);
                }
            }

            printf("gdstn= %.8e %.8e %.8e %.8e %.8e %.8e \n",
                   gdstn[0], gdstn[1], gdstn[2],
                   gdstn[3], gdstn[4], gdstn[5]);
            printf("pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",
                   val[0][0], val[1][1], val[2][2],
                   val[1][2], val[0][2], val[0][1]);

            printf("gdstn/pstn = %.8e %.8e %.8e %.8e %.8e %.8e \n",
                   gdstn[0]/val[0][0], gdstn[1]/val[1][1], gdstn[2]/val[2][2],
                   gdstn[3]/val[1][2], gdstn[4]/val[0][2], gdstn[5]/val[0][1]);

#endif
        } /*if masternode */
#endif


#if 0
        if (home->myDomain == 0) {
            printf("param->delpStrain[0]= %e\n",param->delpStrain[0]);
            printf("param->delpStrain[1]= %e\n",param->delpStrain[1]);
            printf("param->delpStrain[2]= %e\n",param->delpStrain[2]);
            printf("param->delpStrain[3]= %e\n",param->delpStrain[3]);
            printf("param->delpStrain[4]= %e\n",param->delpStrain[4]);
            printf("param->delpStrain[5]= %e\n",param->delpStrain[5]);
        }
#endif                                                                  
#if 0
        if (home->myDomain == 0) {
            printf("param->delpSpin[0]= %e\n",param->delpSpin[0]);
            printf("param->delpSpin[1]= %e\n",param->delpSpin[1]);
            printf("param->delpSpin[2]= %e\n",param->delpSpin[2]);
            printf("param->delpSpin[3]= %e\n",param->delpSpin[3]);
            printf("param->delpSpin[4]= %e\n",param->delpSpin[4]);
            printf("param->delpSpin[5]= %e\n",param->delpSpin[5]);
        }

#endif                                                                  

    return;
}
