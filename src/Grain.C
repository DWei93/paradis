/****************************************************************************
* 
*	Module:		Grain.c
*
*	Description:	This module contains all the subroutines needed in 
*                  	grain/twin boundary model.
*
*	07/06/2016:	created by H. Fan, D. Wei 
*
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "Home.h" 
#include "Grain.h"
#include "Util.h" 

bool IsOnGB(Node_t *node)
{
    return(node->onGB > -1);
}

void PrintNeighbor(Home_t *home, Node_t *node){
	int 	i;
	Node_t	*nbr;
	if(node == (Node_t *)NULL){
	    printf("node is NULL\n");
	    return;
	}
	for(i=0; i<node->numNbrs; i++){
	    printf("node (%d,%d) arm %d is (%d,%d):\n",
	           node->myTag.domainID, node->myTag.index,
	           i, node->nbrTag[i].domainID, node->nbrTag[i].index);

	    nbr = GetNeighborNode(home, node, i);
	    PrintNode(nbr);
	}
	return;
}

void FormatMatrix(real8 M[3][3], const char *msg){
	int i;
	printf("%s matrix:\n",msg);
	for (i = 0; i < 3; i++){
	    printf("%20.15e %20.15e %20.15e\n", M[i][0], M[i][1], M[i][2]);
	}

	return;
}

void FormatMatrix66(real8 M[6][6], const char *msg){
	int i;
	printf("%s matrix:\n",msg);
	for (i = 0; i < 6; i++){
		printf("%18.7e %18.7e %18.7e %18.7e %18.7e %18.7e\t\n", M[i][0], M[i][1], M[i][2],
			  M[i][3], M[i][4], M[i][5]);
	}
	return;
}

void FormatVector(real8 vec[3], const char *msg){
	printf("%s vector:\n",msg);
	printf("{%.15f,%.15f,%.15f}\n", vec[0], vec[1], vec[2]);
}


void InitGrainStructure(Grain_t *grain, int numBounds){
	grain->nBounds = numBounds;
    grain->GB = (GB_t **)malloc(numBounds*sizeof(GB_t *));	
    grain->sign = (real8 *)malloc(numBounds*sizeof(real8));
    return;
}

void InitGrain(Home_t *home, Param_t *param){
	int		    i,j;
	Grain_t		*grain;
    GB_t        *GB;
	real8		HalfMatrixThick, maxSeg;
	real8		TwinOrientation[3][3], temR[3][3];
	int		    nGrain, type, GBID;
	real8		twinRatio;
    real8       I[3][3] = {
                            {1.0e+00, 0.0e+00, 0.0e+00},
                            {0.0e+00, 1.0e+00, 0.0e+00},
                            {0.0e+00, 0.0e+00, 1.0e+00}
                          };

	/*the following variables are placed in param.h*/
	int		    numTGrains;
	real8		halfTwinThick;
	int		    invFlag;
	real8		lx=0, ly=0, lz=0;	/*half cube size in param.h*/

	lx = param->xBoundMax;
	ly = param->yBoundMax;
	lz = param->zBoundMax;
	maxSeg = param->maxSeg;

	if(param->polycrystal == 0)param->nGrain = 0;
	nGrain = param->nGrain;

	for(i=0; i<3; i++){
	    VECTOR_COPY(home->rotMatrix[i], param->matrixOrientation[i]);
	}
    GenerateFCCBurgInfo(&(home->burgData), param->matrixOrientation);
    GenerateFCCBurgInfo(&(home->crystBurgData), I);
	Matrix33Invert(home->rotMatrix, home->rotMatrixInverse);

	if(fabs(home->rotMatrix[0][0]) > 0.9999999999 &&  
	   fabs(home->rotMatrix[1][1]) > 0.9999999999 &&
	   fabs(home->rotMatrix[2][2]) > 0.9999999999){
	    param->useLabFrame = 0;
	}else{
	    param->useLabFrame = 1;
	}

    if(param->polycrystal == 0){
        if(param->xBoundType == Periodic && 
           param->yBoundType == Periodic && 
           param->zBoundType == Periodic ){
            type = 1;
        }else type = 2;
    }else{
        if(param->nGrain == 1){
            if(param->xBoundType == Periodic && 
               param->yBoundType == Periodic && 
               param->zBoundType == Periodic ){
                type = 3;
            }else{
                Fatal("at file %s line %d", __FILE__, __LINE__);
            }
        }else{
            if(param->xBoundType == Periodic && 
               param->yBoundType == Periodic && 
               param->zBoundType == Periodic ){
                type = 4;
            }else{
                type = 5;
            }

        }
    }

    if(type == 3)Fatal("Not supposed this type now");
/*
 *	Type I: polycrystal = 0, nGrain = 0, PBC
 * 	origin PBC condition.
 */	

/*
 *	Type II: polycrystal = 0, nGrain = 0, Free surface.
 *	single crystal, free surface.
 */
	if(param->polycrystal ==0 ){
        home->nGrains = 0;
        home->nGBs = 0;
	    if(home->myDomain == 0){
	        FormatMatrix(param->matrixOrientation,"Orientation");
	    }
	    return;
	}

/*
 *	Type III: polycrystal = 1, nGrain =1, PBC 
 *	only have grain boundary
 */
	if(nGrain == 0)Fatal("polycrstal type, but nGrain = 0");

	home->nGrains = param->nGrain;				/*num of grain and twins */

	home->grainKeys	= (Grain_t **)malloc(home->nGrains * sizeof(Grain_t *));
	param->useLabFrame = 1;
	

/*
 *	Initialize twin grain.
 */
	twinRatio = param->twinRatio;
	numTGrains = nGrain - 1;

	real8	CDeltaLab[3] = {-1,-1,2}, CDelta[3];
	real8	planeABC[3] = {1,1,1};
	NormalizeVec(CDeltaLab);
	NormalizeVec(planeABC);

	Matrix33Vector3Multiply(param->matrixOrientation, CDeltaLab, CDelta);
    RotateWithAxis(CDelta, PI, temR);
	Matrix33Mult33(temR, param->matrixOrientation, TwinOrientation);

	Matrix33Vector3Multiply(param->matrixOrientation, planeABC, param->TBNormal);
    RotateWithAxis(param->TBNormal, PI, param->rotateNormal);


/*
 *	    Type V: polycrystal = 1, nGrain >= 2, PBC.
 *	    polycrystal, twin, RVE, XOZ plane is twin plane
 */
    if(type == 4){
	    halfTwinThick = twinRatio * ly / (real8)numTGrains;
	    HalfMatrixThick = ly / (real8)numTGrains - halfTwinThick;
	    	
	    if (HalfMatrixThick <= 0){
	        printf("the thickness of matrix is %f\n ", HalfMatrixThick);
	        Fatal("the thickness of matrix is negative value at %s line %d\n", __FILE__, __LINE__);
	    }

        home->nGBs = home->nGrains;
        home->GBKeys = (GB_t **)malloc(home->nGrains*sizeof(GB_t *));   
        for(i=0; i<home->nGrains; i++){
            home->GBKeys[i] = (GB_t *)calloc(1, sizeof(GB_t));
            home->GBKeys[i]->GBID = i;
            VECTOR_ZERO(home->GBKeys[i]->n);
            VECTOR_ZERO(home->GBKeys[i]->p);
            home->GBKeys[i]->n[1] = 1.0;
            home->GBKeys[i]->p[1] = -ly + HalfMatrixThick + 2*((real8)i)*halfTwinThick;
            home->GBKeys[i]->nbrGrain[0] = i;
            home->GBKeys[i]->nbrGrain[1] = ((i == home->nGrains -1) ? 0 : i+1);
        }

	    for(i=0; i<nGrain; i++){
	        home->grainKeys[i] = (Grain_t *)calloc(1, sizeof(Grain_t));
	        grain = home->grainKeys[i];
	        InitGrainStructure(grain, 2);
    
	        grain->type = (i==0) ? FOLDED_GRAIN : ODINARY_GRAIN;
	        grain->cy = -ly + (2 * ((real8)i) - 1) * (HalfMatrixThick + halfTwinThick);
            grain->cx = 0.0;
            grain->cz = 0.0;
            grain->grainID = i;

            GBID = i -1;
            if(GBID < 0) GBID = nGrain - 1;
            grain->GB[0] = home->GBKeys[GBID];
            grain->sign[0] = -1.0;
            grain->GB[1] = home->GBKeys[i];
            grain->sign[1] =  1.0;

            if(grain->type == ODINARY_GRAIN || i%2 == 0){
                memcpy(grain->orientation, TwinOrientation, sizeof(real8 [3][3]));
            }else{
                memcpy(grain->orientation, param->matrixOrientation, sizeof(real8 [3][3]));
            }

            GenerateFCCBurgInfo(&(grain->burgData), grain->orientation);
	        invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	        if(!invFlag){
	            Fatal("Invert twin orientation faild!\n");
	        }
        }

        for(i=0; i<home->nGBs; i++){
            GenerateGBBurgInfo(home, home->GBKeys[i]);    
        }
        return;
    }

	if(nGrain == 2){
/*
 *	    Type VI: polycrystal = 1, nGrain =2, Free surface.
 *	    simulat twin sample.
 */
	    real8   hMatrix, p[4][3];
	    real8   interPos[3], bound, rho;
	    int	    inner = 0;
	    Plane_t plane;
	    Line_t  edge;

	    VECTOR_COPY(plane.normal, param->TBNormal);

        
	    for(i=0; i<3; i++){
	        VECTOR_ZERO(edge.direction);
	        VECTOR_ZERO(plane.point);

	        edge.direction[i] = 1;
	        if(DotProduct(param->TBNormal, edge.direction) > 0){
	            rho = 1.0 - twinRatio; }
	        else{
	            rho = twinRatio;  }

	        if(i == 0){
	            bound = lx;
	            hMatrix = 2 * lx * rho - lx; 
	            p[0][0] = 0; p[0][1] = -ly; p[0][2] = -lz;
	            p[1][0] = 0; p[1][1] = -ly; p[1][2] =  lz;
	            p[2][0] = 0; p[2][1] =  ly; p[2][2] = -lz;
	            p[3][0] = 0; p[3][1] =  ly; p[3][2] =  lz;
	        }else if(i == 1){
	            bound = ly;
	            hMatrix = 2 * ly * rho - ly; 
	            p[0][0] = -lx; p[0][1] = 0; p[0][2] = -lz;
	            p[1][0] =  lx; p[1][1] = 0; p[1][2] = -lz;
	            p[2][0] = -lx; p[2][1] = 0; p[2][2] =  lz;
	            p[3][0] =  lx; p[3][1] = 0; p[3][2] =  lz;
	        }else if(i == 2){
	            bound = lz;
	            hMatrix = 2 * lz * rho - lz; 
	            p[0][0] = -lx; p[0][1] = -ly; p[0][2] = 0;
	            p[1][0] =  lx; p[1][1] = -ly; p[1][2] = 0;
	            p[2][0] = -lx; p[2][1] =  ly; p[2][2] = 0;
	            p[3][0] =  lx; p[3][1] =  ly; p[3][2] = 0;
	        }
	        plane.point[i] = hMatrix;	        

	        inner = 1;
	        for(j=0; j<4; j++){
	            VECTOR_COPY(edge.point, p[j]);

	            if(!IntersectionOfLinePlane(&edge, &plane, interPos)){
	                inner = 0; break; }
	            
	            if(fabs(interPos[i]) > bound){
	                inner = 0; break; }
	        }/* for (j=0; j ...)*/
	        if(inner == 1)break;

	    }/* for(i=0; i ...)*/
	    
	    if(inner == 0)Fatal("InitGrain: can not find TB plane");
	        
        home->nGBs = 1;
        home->GBKeys = (GB_t **)malloc(sizeof(GB_t *));
        i = 0;
        home->GBKeys[i] = (GB_t *)calloc(1, sizeof(GB_t));
        home->GBKeys[i]->GBID = i;
        VECTOR_COPY(home->GBKeys[i]->n, param->TBNormal);
        VECTOR_COPY(home->GBKeys[i]->p, plane.point);
        home->GBKeys[i]->nbrGrain[0] = i;
        home->GBKeys[i]->nbrGrain[1] = i+1;

	    home->grainKeys[0] = (Grain_t *)calloc(1, sizeof(Grain_t));
        grain = home->grainKeys[0];
	    InitGrainStructure(grain, 1);
	    grain->type = ODINARY_GRAIN;
        grain->grainID = 0;
        grain->sign[0] = 1.0;
        grain->GB[0] = home->GBKeys[0];
        memcpy(grain->orientation, param->matrixOrientation, sizeof(real8 [3][3]));
        GenerateFCCBurgInfo(&(grain->burgData), grain->orientation);
	    invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	    if(!invFlag){
	        Fatal("Invert twin orientation faild!\n");
	    }

	    home->grainKeys[1] = (Grain_t *)calloc(1, sizeof(Grain_t));
	    grain = home->grainKeys[1];
	    InitGrainStructure(grain, 1);
	    grain->type = ODINARY_GRAIN;
        grain->grainID = 1;
        grain->sign[0] = -1.0;
        grain->GB[0] = home->GBKeys[0];
    
        memcpy(grain->orientation, TwinOrientation, sizeof(real8 [3][3]));
	    invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	    if(!invFlag){
	        Fatal("Invert twin orientation faild!\n");
	    }

	    if(home->myDomain == 0){
	        FormatMatrix(param->matrixOrientation,"basis orientation");
	        printf("\n");
	        FormatMatrix(TwinOrientation,"twin orientation");
	        printf("\n");
	        FormatVector(param->TBNormal, "twin boundary normal");
	        printf("\n");
	        FormatVector(CDelta, "CDelta (use the axis to rotate basis to twin)");
	    }
	}else{
	    Fatal("Not suppose free surface polycrystal problem");
	}

	return;
}

void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int& inGrain, int& onGB){
    int         i, j, nGrains, type, breakLoop = 0;
    real8       vec[3], dis, pos[3], sign;
	Grain_t		*grain;
	Param_t		*param;

    inGrain = -1;
    onGB = -1;
	param = home->param;

	if(param->polycrystal == 0)return;

	nGrains = home->nGrains;
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    FoldBox(param, pos, pos+1, pos+2);

    for(i=0; i<nGrains; i++){
        if(breakLoop)break;
        grain = home->grainKeys[i];
        for(j=0; j<grain->nBounds; j++){
            type = PointPlaneIntersection(pos, grain->GB[j]->n, 
                   grain->GB[j]->p, vec, &dis);
            dis *= grain->sign[j];

            if(type == POINT_INTERSECTION){
                if(dis < 0.0){
                    inGrain = grain->grainID;
                }else if (dis == 0){
                    inGrain = -1;
                }else{
                    inGrain = ((grain->GB[j]->nbrGrain[0] == grain->grainID) ?
                              grain->GB[j]->nbrGrain[1] : grain->GB[j]->nbrGrain[0]);
                }
                onGB = grain->GB[j]->GBID;
                breakLoop = 1; 
                break;
            }

            if(grain->type == FOLDED_GRAIN){
                if(dis < 0.0){
                    inGrain = grain->grainID;
                    onGB = -1;
                    breakLoop = 1; 
                    break;
                }else{
                    continue;
                }
            }else{
                if(dis > 0.0){
                    break;
                }else{
                    if(j == grain->nBounds - 1){
                        inGrain = grain->grainID;
                        onGB = -1;
                        breakLoop = 1;
                        break;
                    }
                }
            }
        }   /*for(j=0; j<grain->nBounds; j++) */
    }   /*for(i=0; i<home->nGrains; i++) */

    if(breakLoop == 0)Fatal("at file %s at line %d", __FILE__, __LINE__);

	return;
}

Node_t *GetinDataNodeFromTag(InData_t *inData, Tag_t tag){
	Node_t 	*node;
	int  	i;
	
	for (i = 0; i < inData->nodeCount; i++){
		node = &inData->node[i];
		if(node == (Node_t *)NULL)continue;
		
		if(node->myTag.domainID == tag.domainID &&
		   node->myTag.index  == tag.index)return(node);
	}
	return((Node_t *)NULL);
}

int JudgeNodeInCuboid(real8 *x1, real8 *y1, real8 *z1, real8 x2, real8 y2, real8 z2){
	return(fabs(*x1) <= x2 && fabs(*y1) <= y2 && fabs(*z1) <= z2);
}

void GetTBBurg_FCC(int nID, real8 burgTB[3]){
    if(nID == 0){
        /* Adelta*/
        burgTB[0] = -1;
        burgTB[1] =  2;
        burgTB[2] = -1;
    }else if(nID == 1){
        /* Bdelta*/
        burgTB[0] =  2;
        burgTB[1] = -1;
        burgTB[2] = -1;
    }else if (nID == 2){
        /* Cdelta*/
        burgTB[0] = -1;
        burgTB[1] = -1;
        burgTB[2] =  2;
    }else Fatal("GetTBBurg_FCC: normal ID is %d\n", nID);
    NormalizeVec(burgTB);
    burgTB[0] /= sqrt(3.0);
    burgTB[1] /= sqrt(3.0);
    burgTB[2] /= sqrt(3.0);
}

void VectorLabToCrystal(Home_t *home, Node_t *node1, Node_t *node2, 
                        real8 lab[3], real8 *cryst)
{
        if(home->param->polycrystal){
            int     inGrain, onGB;
            Grain_t *grain;
            GB_t    *GB;

            GetArmGrain(home, node1, node2, inGrain, onGB);
            if(onGB > -1){
                GB = home->GBKeys[onGB];
                grain = home->grainKeys[GB->nbrGrain[0]];
            }else{
                grain = home->grainKeys[inGrain];    
            }
            Matrix33Vector3Multiply(grain->orientationInv, lab, cryst);
        }else{
            if(!home->param->useLabFrame){
                VECTOR_COPY(cryst, lab);
            }else{
                Matrix33Vector3Multiply(home->rotMatrixInverse, 
                                        lab, cryst);
            }
        }
        return;
}

void VectorCrystalToLab(Home_t *home, Node_t *node1, Node_t *node2, 
                        real8 cryst[3], real8 *lab)
{
        if(home->param->polycrystal){
            int     inGrain, onGB;
            Grain_t *grain;
            GB_t    *GB;

            GetArmGrain(home, node1, node2, inGrain, onGB);
            if(onGB > -1){
                GB = home->GBKeys[onGB];
                grain = home->grainKeys[GB->nbrGrain[0]];
            }else{
                grain = home->grainKeys[inGrain];    
            }
            Matrix33Vector3Multiply(grain->orientation, cryst, lab);
        }else{
            if(!home->param->useLabFrame){
                VECTOR_COPY(lab, cryst);
            }else{
                Matrix33Vector3Multiply(home->rotMatrix, 
                                        cryst, lab);
            }
        }
        return;
}

int BurgID_FCC(BurgInfo_t *burgData, real8 burg[3])
{
        int     i;
        real8   crossProd[3];
        
        for(i=0; i<burgData->numBurgVectors; i++){
            cross(burgData->burgList[i], burg, crossProd);
            if(Normal(crossProd) < 1.0E-3)break;
        }
        
        return(i);
}

int GlidePlaneID_FCC(BurgInfo_t *burgData, real8 plane[3])
{
        int     i;
        real8   crossProd[3];
        
        for(i=0; i<burgData->numPlanes; i++){
            cross(burgData->planeList[i], plane, crossProd);
            if(Normal(crossProd) < 1.0E-3)break;
        }
        
        return(i);
}

int	BurgType_FCC(real8 burg[3]){

	real8 eps = 1e-2;
	real8 burgMag;

	burgMag = Normal(burg);

/*	
 *	0: zero burgers
 */
    if(burgMag < eps)return(ZERO_BURG);

/*
 * 	1: basic Burgers vector
 */		
	if(fabs(burgMag - 1) < eps){
		return(BASIC_BURG);
	}
/*	
 *	2: TB burgers
 */
	
	if(fabs(burgMag - 0.5773502691896258) < eps) {
		return(SHOCKLEY_BURG);
	}

/*	
 *	3: junction 1 burgers, like [Sqrt(2),0,0]
 */
	if(fabs(burgMag - 1.414213562373095) < eps) {
		return(JUNCTION1_BURG);
	}

/*	
 *	3: junction 2 burgers, like DA+DB
 */
	if(fabs(burgMag - 1.732050807568877) < eps) {
		return(JUNCTION2_BURG);
	}
	return(UNKNOWN_BURG);

}


int GlidePlaneType_FCC(BurgInfo_t *burgData, real8 plane[3])
{
    int planeID;

    planeID = GlidePlaneID_FCC(burgData, plane);

    if(planeID >= NUM_PLANES_IN_FCC) return(UNKOWN_PLANE);
 
    return(burgData->planeType[planeID]);
}

int GetBurgID_FCC(Home_t *home, Node_t *node, Node_t *nbr){

        int         inGrain, onGB, segID;
        real8       burg[3];
        BurgInfo_t  *burgData;

        if(home->param->polycrystal){
            GetArmGrain(home, node, node, inGrain, onGB);
            if(onGB > -1){
                burgData = &(home->GBKeys[onGB]->burgData);
            }else{
                burgData = &(home->grainKeys[inGrain]->burgData);
            }
        }else{
            burgData = &(home->burgData);
        }
        
	    segID  = GetArmID(home, node, nbr);
        burg[0] = node->burgX[segID];
        burg[1] = node->burgY[segID];
        burg[2] = node->burgZ[segID];        
    
        return(BurgID_FCC(burgData, burg)); 
}


/*-------------------------------------------------------------------------
 *
 *      Function:    HasShortSeg
 *
 *      Description: Check all segments attached to the specified node
 *                   to see if any of the segments is considered 'short'
 *                   according to the provided criteria.
 *
 *      Parameters:
 *          in:  node          Pointer to node structure
 *          in:  shortSegLenSq Square of the length under which any segment
 *                             should be considered a short segment.
 *
 *      Returns: 1 if the node has any short segments
 *               0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int HasShortSeg(Home_t *home, Node_t *node, real8 shortSegLenSq)
{
        int segID;
		
		if(node == (Node_t *)NULL){
			printf("Task %d: %s, line %d -- Node = NULL",
				   home->myDomain, __FILE__, __LINE__);
			return(0);
		}

        for (segID = 0; segID < node->numNbrs; segID++) {
            Node_t *nbrNode;

            nbrNode = GetNeighborNode(home, node, segID);

            if (nbrNode != (Node_t *)NULL) {
                real8  segLenSq;
                real8  segVector[3];

                segVector[0] = node->x - nbrNode->x;
                segVector[1] = node->y - nbrNode->y;
                segVector[2] = node->z - nbrNode->z;

                ZImage(home->param, &segVector[0], &segVector[1],
                       &segVector[2]);

                segLenSq = DotProduct(segVector, segVector);

                if (segLenSq < shortSegLenSq) {
                    return(1);
                }
            }
        }

        return(0);
}

/****************************************************************************  
*
*	Assign grain number to node then add this operation to 
*	Oplist to send to other domains. 
*
****************************************************************************/  
void AssignGrain(Home_t *home, Node_t *node, int inGrain, int onGB) 
{
	if(node == (Node_t *)NULL){
		printf("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		return;
	}

	node->inGrain = inGrain;
	node->onGB = onGB;

	Operate_t *op;

/*
 *   Make sure the buffer allocated for the operation list is
 *   large enough to contain another operation.  If not, increase
 *   the size of the buffer.
 */
    if (home->OpCount >= home->OpListLen) {
		ExtendOpList (home);
	}
        
	op = &home->opList[home->OpCount];

	op->type = UPDATE_GRAIN;
	op->dom1 = node->myTag.domainID;
	op->idx1 = node->myTag.index;
	op->inGrain = inGrain;
	op->onGB = onGB;

	home->OpCount++;

	return;
}


void AssignConstraint(Home_t *home, Node_t *node, int constraint) 
{
	if(node == (Node_t *)NULL){
		printf("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		return;
	}

	node->constraint = constraint;

	Operate_t *op;

/*
 *   Make sure the buffer allocated for the operation list is
 *   large enough to contain another operation.  If not, increase
 *   the size of the buffer.
 */
    if (home->OpCount >= home->OpListLen) {
		ExtendOpList (home);
	}
        
	op = &home->opList[home->OpCount];

	op->type = UPDATE_CONSTRAINT;
	op->dom1 = node->myTag.domainID;
	op->idx1 = node->myTag.index;
	op->inGrain = constraint;

	home->OpCount++;
	return;
}


int NodeOnGB(Home_t *home, Node_t *node){
	int 	i, numNbrs, burgType, type = 0;
	real8	burg[3];
	Node_t 	*nbr;
	
	if(node == (Node_t *)NULL){
		Fatal("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);

	}

	if(!IsOnGB(node))return(INGRAIN_NODE);

	numNbrs = node->numNbrs;
	
	for(i=0; i<numNbrs; i++){
		burg[0] = node->burgX[i];
	    burg[1] = node->burgY[i];
	    burg[2] = node->burgZ[i];
		
	    burgType = BurgType_FCC(burg);
		
        switch(burgType){
            case SHOCKLEY_BURG:
                type |= TB_NODE;
                break;
            case BASIC_BURG:
                type |= INTERSECT_NODE;
                break;
            default:
                type |= JUNCTION_NODE;
                break;
        }
	}
	return(type);
}

int AdjustNodeDirection(Home_t *home, Node_t *node, real8 *vec){

    int 	NM;
    real8 	GC[3][3];

	NM = 0;

	GetGlideConstraintList(home, node, &NM, GC);
	ApplyConstraintsToVelocity(NM, GC, vec);

	return(NM);
}


real8 ClosedGB(Home_t *home, real8 *pos, GB_t **GB){
	real8	    tempDistance, vec[3];
	real8 	    minDistance = 1e14;
	GB_t 	    *temGB;
	    	
    if(!home->param->polycrystal)return(0.0);
    int     i;
    for(i=0; i<home->nGBs; i++){
        temGB = home->GBKeys[i];
        PointPlaneIntersection(pos, temGB->n, temGB->p, vec, &tempDistance);
        if(fabs(tempDistance) < fabs(minDistance)){
            *GB = temGB;
            minDistance = tempDistance;
        }
    }
	return(minDistance);
}

int GetGBIndexOfGrain(Grain_t *grain, int inGB){
    int i;
    for(i=0; i<grain->nBounds; i++){
         if(grain->GB[i]->GBID == inGB)return(i);
    }
    return(-1);
}

bool GetPublicGB(Grain_t *grain1, Grain_t *grain2, GB_t **GB)
{
    int id, i;
    for(i=0; i<grain1->nBounds; i++){
        id = GetGBIndexOfGrain(grain1, grain2->GB[i]->GBID);
        if(id > -1){
            *GB = grain1->GB[id];
            return(1);
        } 
    }
    return(0);
}

bool NewPosToGB(Home_t *home, Node_t *node, int inGrain, int onGB, real8 *newPos){
    
    real8       dis, vec[3], vec1[3], vec2[3];
    real8       pos[3], interPos[3], vec3[3]; 
    int         currInGrain, currOnGB, type, index;
    Param_t     *param;
    GB_t        *GB;
    Grain_t     *grain;
    Plane_t     *plane;

    param = home->param;

	if(node == (Node_t *)NULL){
        Fatal("Task %d: %s, line %d -- Node = NULL",
               home->myDomain, __FILE__, __LINE__);
	}

	pos[0] = node->x; pos[1] = node->y; pos[2] = node->z;
	FoldBox(param, pos, pos+1, pos+2);

    if(onGB < 0)Fatal("in %s at %d", __FILE__, __LINE__);
    GB = home->GBKeys[onGB];

    VECTOR_COPY(vec1, GB->n);
    VECTOR_COPY(vec2, GB->p);
    if(inGrain > -1){
        grain = home->grainKeys[inGrain];
        index = GetGBIndexOfGrain(grain, onGB);
        if(index < 0)Fatal("in %s at %d", __FILE__, __LINE__);
        vec2[0] -= (grain->sign[index] * param->rann * GB->n[0]);
        vec2[1] -= (grain->sign[index] * param->rann * GB->n[1]);
        vec2[2] -= (grain->sign[index] * param->rann * GB->n[2]);
    }
    
    VECTOR_COPY(vec, vec1);
    AdjustNodeDirection(home, node, vec);
    if(Normal(vec) < 1E-12){
        printf("Warning: in %s at %d", __FILE__, __LINE__);
        VECTOR_ZERO(newPos);
        return(0);
    }
    NormalizeVec(vec);
    
    type = LinePlaneIntersection(vec, pos, vec1, vec2, vec3, newPos); 
    if(type != POINT_INTERSECTION){
        printf("Warning: in %s at %d", __FILE__, __LINE__);
        VECTOR_ZERO(newPos);
        return(0);
    }
    
    FoldBox(param, newPos, newPos+1, newPos+2);
    return(1);
}


/****************************************************************************  
*
*	In order to avoid numerical error, we must ensure the node is pulled 
*	back into its grain.
*
*****************************************************************************/
void PullBack(Home_t *home, Node_t *node){

    int     currInGrain, currOnGB;
	real8	newPos[3];
	
	if(node == (Node_t *)NULL){
        printf("Task %d: %s, line %d -- Node = NULL",
		       home->myDomain, __FILE__, __LINE__);
        return;
	}

	newPos[0] = node->oldx;
	newPos[1] = node->oldy;
	newPos[2] = node->oldz;
	
	NodeLocation(home,newPos[0],newPos[1],newPos[2], currInGrain, currOnGB);

    if(IsOnGB(node)){
        if(node->onGB != currOnGB){
            if(!NewPosToGB(home, node, -1, node->onGB, newPos)){
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
        }
    }else{
        if(currOnGB == -1){
            if(node->inGrain != currInGrain){
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
        }else{
            if(!NewPosToGB(home, node, node->inGrain, currOnGB, newPos)){
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
        }
    }

	FoldBox(home->param, newPos, newPos+1, newPos+2);
    RepositionNode(home,newPos,&node->myTag,1);

	return;
}

int SameTag(Tag_t *tag1, Tag_t *tag2){

	return(tag1->domainID == tag2->domainID &&
	       tag1->index == tag2->index);	
}

real8 Distance(Home_t *home,real8 *p1,real8 *p2){
	
	real8 vec[3];
	
	FoldBox(home->param, p1, p1+1, p1+2);
	FoldBox(home->param, p2, p2+1, p2+2);
	
	vec[0] = p1[0]-p2[0];
	vec[1] = p1[1]-p2[1];
	vec[2] = p1[2]-p2[2];

	ZImage(home->param,&vec[0],&vec[1],&vec[2]);

	return(Normal(vec));
}

real8 PointPlaneDistance(real8 *point, Plane_t *plane){
	real8 vec[3];
	
	vec[0] = point[0] - plane->point[0];
	vec[1] = point[1] - plane->point[1];
	vec[2] = point[2] - plane->point[2];
	
	return(DotProduct(vec, plane->normal));
}

int IntersectionOfLinePlane(Line_t *line, Plane_t *plane, real8 *p){
	real8	t,vec1[3],vec2[3],dot1,dot2;

	vec1[0] = plane->point[0] - line->point[0];
	vec1[1] = plane->point[1] - line->point[1];
	vec1[2] = plane->point[2] - line->point[2];
	
	
	dot1 = DotProduct(vec1, plane->normal);
	dot2 = DotProduct(plane->normal, line->direction);
	
	if(fabs(dot2) > 1.0E-6){
	    t = dot1/dot2;
	    p[0] = line->point[0] + t * line->direction[0];
	    p[1] = line->point[1] + t * line->direction[1];
	    p[2] = line->point[2] + t * line->direction[2];
	    return(1);
	}else{
	    return(0);
	}
}

int NodeOnLine(real8 *p, Line_t *line){
	real8	eps = 1e-6,vec[3],dot;
	int		onLine = 1;
	
	onLine &= (fabs(p[0] - line->point[0] ) < eps);
	onLine &= (fabs(p[1] - line->point[1] ) < eps);
	onLine &= (fabs(p[2] - line->point[2] ) < eps);
	
	if(onLine)return(onLine);
	
	vec[0] = p[0] - line->point[0];
	vec[1] = p[1] - line->point[1];
	vec[2] = p[2] - line->point[2];
	
	NormalizeVec(vec);
	dot = DotProduct(vec,line->direction);
	
	if(1-fabs(dot) < eps)return(1);
	
	
	return(0);
}

int LineOnPlane(Line_t *line, Plane_t *plane){
	real8	dot1, vec1[3], dot2;
	real8 	eps = 1e-3;
	
	dot1 = DotProduct(line->direction, plane->normal);
	vec1[0] = line->point[0] - plane->point[0];
	vec1[1] = line->point[1] - plane->point[1];
	vec1[2] = line->point[2] - plane->point[2];
	NormalizeVec(vec1);
	dot2 = DotProduct(vec1, plane->normal);
	
	if(fabs(dot1) < eps && fabs(dot2) <eps){
		return(1);
	}
	return(0);
}


void GetArmGrain(Home_t *home, Node_t *node1, Node_t *node2, int &inGrain, int &onGB){

    inGrain = -1;
    onGB = -1;
    if(home->param->polycrystal == 0)return;

    if(IsOnGB(node1)){
        if(IsOnGB(node2)){
            if(node1->onGB == node2->onGB){
                onGB = node1->onGB;
            }else{
                int     i, j;
                bool    found = 0;
                GB_t    *GB1, *GB2;
                GB1 = home->GBKeys[node1->onGB];
                GB2 = home->GBKeys[node2->onGB];

                for(i=0; i<2; i++){
                    if(found)break;
                    for(j=0; j<2; j++){
                        if(GB1->nbrGrain[i] == GB2->nbrGrain[j]){
                            found = 1;
                            break;
                        }
                    }
                }
                inGrain = GB1->nbrGrain[i];
                if(!found)Fatal("file %s at line %d", __FILE__, __LINE__);
            }
        }else{
            inGrain = node2->inGrain;
        }
    }else {
        inGrain = node1->inGrain;
        if(!IsOnGB(node2) && node1->inGrain != node2->inGrain){
            printf("Warning in GetArmGrain: node1 in %d grain, but node2 in %d grain",
                   node1->inGrain, node2->inGrain);
        }
    }
    return;
}


int	IsZeroVector(real8 *vec){
	return(fabs(vec[0]) < 1.0e-10 && 
		   fabs(vec[1]) < 1.0e-10 && 
		   fabs(vec[2]) < 1.0e-10);
}

/****************************************************************************  
 *
 *	Get the element number the corrent position belongs to.
 *
 *****************************************************************************/
int Coor2EleNum(Home_t *home,real8 position[3]){

	int	nx,ny,nz;
	real8	eleLx,eleLy,eleLz;
    Param_t *param;

    param = home->param;
	eleLx = param->meshSizeX;
	eleLy = param->meshSizeY;
	eleLz = param->meshSizeZ;

	nx = (position[0] + param->Lx/2)/eleLx;
	ny = (position[1] + param->Ly/2)/eleLy;
	nz = (position[2] + param->Lz/2)/eleLz;

	if(nx == param->meshNumX) nx--;
	if(ny == param->meshNumY) ny--;
	if(nz == param->meshNumZ) nz--;

	return(nz*param->meshNumX*param->meshNumY + ny*param->meshNumX +nx);

}

void PrintTag(Node_t *node){

    if(node == (Node_t *)NULL){
        printf("(NULL)");
    }else{
        printf("[%d %d %d %d %d %d]-{%e,%e,%e} ",node->myTag.domainID,
               node->myTag.index,node->numNbrs,node->inGrain,
               node->onGB, node->constraint, node->x,node->y,node->z);
	}
    return;
}

void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]){
    
    NormalizeVec(a);
    R[0][0] = (1.0 - cos(theta))*a[0]*a[0] + cos(theta);
    R[0][1] = a[0]*a[1]*(1.0 - cos(theta)) - a[2]*sin(theta);
    R[0][2] = a[0]*a[2]*(1.0 - cos(theta)) + a[1]*sin(theta);
    R[1][0] = a[0]*a[1]*(1.0 - cos(theta)) + a[2]*sin(theta);
    R[1][1] = (1.0 - cos(theta))*a[1]*a[1] + cos(theta);
    R[1][2] = a[1]*a[2]*(1.0 - cos(theta)) - a[0]*sin(theta);
    R[2][0] = a[0]*a[2]*(1.0 - cos(theta)) - a[1]*sin(theta);
    R[2][1] = a[1]*a[2]*(1.0 - cos(theta)) + a[0]*sin(theta);
    R[2][2] = (1.0 - cos(theta))*a[2]*a[2] + cos(theta);
    return;
}

int SurfaceNode(Node_t *node){
	if(node->constraint >= 1 && node->constraint <= 3)return(1);
	return(0);
}

int ExceedBox(Param_t *param, real8 *pos){
	int flags = 0x0000;
        if (param->xBoundMax + 1.0E-6 < fabs(pos[0]))flags |= EXCEEDX;
        if (param->yBoundMax + 1.0E-6 < fabs(pos[1]))flags |= EXCEEDY;
        if (param->zBoundMax + 1.0E-6 < fabs(pos[2]))flags |= EXCEEDZ;
	return(flags);
}

int IntersectionOfLineBox(Param_t *param, Line_t *line, real8 interPos[2][3]){
	int	i,j,k;
	real8	cof;
	Plane_t plane;
	real8	box[3];

	box[0] = param->xBoundMax;
	box[1] = param->yBoundMax;
	box[2] = param->zBoundMax;

	if(ExceedBox(param, line->point)){
	    Fatal("IntersectionOfLineBox: the line exceed the box\n");
	}

	j = 0;
	k = 0;
	VECTOR_ZERO(interPos[0]);
	VECTOR_ZERO(interPos[1]);
	for(i=0; i<3; i++){
	    VECTOR_ZERO(plane.normal);
	    plane.normal[i] = 1.0;
	    for(j=0; j<2; j++){
	        VECTOR_ZERO(plane.point);
	        cof = (j == 0) ? -1.0 : 1.0;
	        plane.point[i] = cof * box[i]; 
	        if(!IntersectionOfLinePlane(line, &plane, interPos[k]))continue;
	        if(ExceedBox(param, interPos[k]))continue;
	        k++;
	        if(k == 2)return(1);
	    }
	}

	return(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     ImageForce
 *      Description:  
 *
 *      Arguments:
 *          x1,y1,z1  Coordinates of segment's first endpoint
 *          x2,y2,z2  Coordinates of segment's second endpoint
 *          bx,by,bz  Components of burgers vector from first endpoint
 *                    to second endpoint
 *          f1        Array in which to return image force at pos1
 *          f2        Array in which to return image force at pos2
 *
 *-------------------------------------------------------------------------*/
void ImageForce(Home_t *home, real8 armGrain, real8 x1, real8 y1, real8 z1,
                  real8 x2, real8 y2, real8 z2, real8 bx, real8 by, real8 bz,
                  real8 f1[3], real8 f2[3])
{

    int         grainID;
    Grain_t    *grain;
    Param_t     *param;

    param = home->param;

    VECTOR_ZERO(f1);
    VECTOR_ZERO(f2);


    return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       IntersectionOfNodes
 *      Description:  
 *
 *      Arguments:
 *          return  0:  No any intersection.
 *                  1:  Point   vec1 is used to be position
 *                  2:  Line    vec1 and vec2 are used to be direction and
 *                              a point on the line, respectively.
 *                  3:  Plane   vec1 and vec2 are used to be normal and
 *                              a point on the plane, respectively.
 *
 *-------------------------------------------------------------------------*/
int NodesIntersection(Home_t *home, Node_t *node1, Node_t *node2,
                             real8 *vec1, real8 *vec2)
{
    int     type;
    int     numGlideConstraintsNode1, numGlideConstraintsNode2;
    real8   glideConstraintsNode1[3][3], glideConstraintsNode2[3][3]; 
    real8   p1[3], p2[3], dis;
    real8   direction1[3], direction2[3];
    
    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    p1[0] = node1->x;
    p1[1] = node1->y;
    p1[2] = node1->z;
    p2[0] = node2->x;
    p2[1] = node2->y;
    p2[2] = node2->z;

    PBCPOSITION(home->param, p1[0], p1[1], p1[2], p2, p2+1, p2+2);

    GetGlideConstraintList(home, node1, 
                           &numGlideConstraintsNode1, glideConstraintsNode1);
    GetGlideConstraintList(home, node2, 
                           &numGlideConstraintsNode2, glideConstraintsNode2);
    
    switch (numGlideConstraintsNode1){
        case 1:
            switch (numGlideConstraintsNode2){
                case 1:
                    type = PlanePlaneIntersection(glideConstraintsNode1[0], p1,
                           glideConstraintsNode2[0], p2, vec1, vec2);
                    break;
                case 2:
                    type = LinePlaneIntersection(glideConstraintsNode2[2], p2, glideConstraintsNode1[0], p1,
                           vec1, vec2);
                    break;
                case 3:
                    type = PointPlaneIntersection(p2, glideConstraintsNode1[0], p1, vec2, &dis); 
                    break;
                default:
                    Fatal("NodesIntersectio: node2 %d",numGlideConstraintsNode2);
                    break;
            }
            break;
        case 2:
            switch (numGlideConstraintsNode2){
                case 1:
                    type = LinePlaneIntersection(glideConstraintsNode1[2], p1, glideConstraintsNode2[0], p2,
                           vec1, vec2); 
                    break;
                case 2:
                    cross(glideConstraintsNode2[0],glideConstraintsNode2[1],direction2); 
                    NormalizeVec(direction2);
                    type = LineLineIntersection(glideConstraintsNode1[2], p1, direction2, p2,
                           vec1, vec2, &dis); 
                    break;
                case 3:
                    type = PointLineIntersection(p2, glideConstraintsNode1[2], p1, vec2, &dis); 
                    break;
                default:
                    Fatal("NodesIntersectio: node2 %d",numGlideConstraintsNode2);
                    break;
            }

            break;
        case 3:
            switch (numGlideConstraintsNode2){
                case 1:
                    type = PointPlaneIntersection(p1, glideConstraintsNode2[0], p2,
                           vec2, &dis); 
                    break;
                case 2:
                    type = PointLineIntersection(p1, glideConstraintsNode2[2], p2, vec2, &dis); 
                    break;
                case 3:
                    type = PointPointIntersection(p1, p2, vec2, &dis); 
                    break;
                default:
                    Fatal("NodesIntersectio: node2 %d",numGlideConstraintsNode2);
                    break;
            }
            break;
        default:
            Fatal("NodesIntersection: node1 %d",numGlideConstraintsNode1);
            break;
    }
    return(type);
}

bool BeforeMergeArm(Home_t *home, Node_t *node1, Node_t *node2, real8 *vec1, real8 *vec2, int *type)
{
    bool    inhibitMerge = 0;
    real8   vec[3];

    *type = NodesIntersection(home, node1, node2, vec1, vec2); 
    switch (*type) {
        case NO_INTERSECTION:
            inhibitMerge = 1;
            break;
        case POINT_INTERSECTION:
            break;
        case LINE_INTERSECTION:
            break;
        case PLANE_INTERSECTION:
            vec[0] = 0.5*(node2->x - node1->x);
            vec[1] = 0.5*(node2->y - node1->y);
            vec[2] = 0.5*(node2->z - node1->z);
            ZImage(home->param, vec, vec+1, vec+2);
            vec2[0] = node1->x + vec[0];
            vec2[1] = node1->y + vec[1];
            vec2[2] = node1->z + vec[2];
            break;
        default:
            Fatal("BeforeMergeArm: wrong type %d", __LINE__);
            break;
    }
    return(inhibitMerge);
}

void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3])
{
        int         i;

        
        burgData->numBurgVectors = NUM_BURGS_IN_FCC;
        burgData->numPlanes   = NUM_PLANES_IN_FCC;
        burgData->numKeyLines = NUM_KEY_LINES_IN_FCC;

        burgData->burgList = (real8 (*)[3])calloc(3*burgData->numBurgVectors, sizeof(real8));
        burgData->planeList = (real8 (*)[3])calloc(3*burgData->numPlanes, sizeof(real8));
        burgData->planeType = (int *)calloc(burgData->numPlanes, sizeof(int));
        burgData->keyLineList = (real8 (*)[3])calloc(3*burgData->numKeyLines, sizeof(real8));

        burgData->numPlanesPerBurg = (int *)calloc(burgData->numBurgVectors, sizeof(int));        
        burgData->burgFirstPlaneIndex = (int *)calloc(burgData->numBurgVectors, sizeof(int));

        burgData->numBurgsPerPlane = (int *)calloc(burgData->numPlanes, sizeof(int));        
        burgData->planeFirstBurgIndex = (int *)calloc(burgData->numPlanes, sizeof(int));

        burgData->numPlanesPerKeyLine = (int *)calloc(burgData->numKeyLines, sizeof(int));        
        burgData->keyLineFirstPlaneIndex = (int *)calloc(burgData->numKeyLines, sizeof(int));
        

        static real8 burgList[NUM_BURGS_IN_FCC][3] = {
/* 
 *          basic Roman-Roman burgers vectors       (6)
 *          magnitude: 1.0
 */
           { 0.7071067811865475,    0.7071067811865475,     0.0e+00   },   // DC -- 0
           { 0.7071067811865475,    -0.7071067811865475,    0.0e+00   },   // BA -- 1
           { 0.7071067811865475,    0.0e+00   ,    0.7071067811865475 },   // DA -- 2
           { 0.7071067811865475,    0.0e+00   ,    -0.7071067811865475},   // BC -- 3
           { 0.0e+00   ,   0.7071067811865475,     0.7071067811865475 },   // DB -- 4
           { 0.0e+00   ,   0.7071067811865475,     -0.7071067811865475},   // AC -- 5
/* 
 *          Shocklay Roman-Greek burgers vectors    (12)
 *          magnitude: 1/Sqrt(3)
 */
           { 0.4714045207910317,       0.2357022603955158,        0.2357022603955158 },      // Db -- 6
           {-0.4714045207910317,       0.2357022603955158,        0.2357022603955158 },      // dB -- 7
           { 0.4714045207910317,       -0.2357022603955158,       0.2357022603955158 },      // gA -- 8
           { 0.4714045207910317,       0.2357022603955158,        -0.2357022603955158},      // aC -- 9
           { 0.2357022603955158,       0.4714045207910317,        0.2357022603955158 },      // Da -- 10
           {-0.2357022603955158,       0.4714045207910317,        0.2357022603955158 },      // gB -- 11
           { 0.2357022603955158,       -0.4714045207910317,       0.2357022603955158 },      // dA -- 12
           { 0.2357022603955158,       0.4714045207910317,        -0.2357022603955158},      // bC -- 13
           { 0.2357022603955158,       0.2357022603955158,        0.4714045207910317 },      // Dg -- 14
           {-0.2357022603955158,       0.2357022603955158,        0.4714045207910317 },      // aB -- 15
           { 0.2357022603955158,       -0.2357022603955158,       0.4714045207910317 },      // bA -- 16
           { 0.2357022603955158,       0.2357022603955158,        -0.4714045207910317},      // dC -- 17
/* 
 *          compound Roman-Roman burgers vectors    (3) 
 *          magnitude: Sqrt(2)
 */
           { 1.414213562373095,         0.0e+00,                    0.0e+00           },   // DC+BA -- 18
           { 0.0e+00,                   1.414213562373095,          0.0e+00           },   // DC+AB -- 19
           { 0.0e+00,                   0.0e+00,                    1.414213562373095 }    // CA+DB -- 20
        };

        for(i=0; i<NUM_BURGS_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, burgList[i], burgData->burgList[i]);
        }

        static int numPlanesPerBurg[NUM_BURGS_IN_FCC] = { 3, 3, 3, 3, 3, 3, 
                                                          2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                          4, 4, 4 }; 
        memcpy(burgData->numPlanesPerBurg,numPlanesPerBurg,sizeof(int [NUM_BURGS_IN_FCC]));

        static int burgFirstPlaneIndex[NUM_BURGS_IN_FCC] = { 0, 3, 6, 9, 12, 15, 
                                                            17,19,21,23,25,27,29,31,33,35,37,39,
                                                            43,47,51};
        memcpy(burgData->burgFirstPlaneIndex,burgFirstPlaneIndex,sizeof(int [NUM_BURGS_IN_FCC]));

        real8 planeList[NUM_PLANES_IN_FCC][3] = {
            {-0.5773502691896258,      0.5773502691896258,       -0.5773502691896258 },          // a plane  ID: 0
            {0.5773502691896258,       -0.5773502691896258,       -0.5773502691896258},          // b plane  ID: 1
            {-0.5773502691896258,      -0.5773502691896258,       0.5773502691896258 },          // g plane  ID: 2
            {0.5773502691896258,       0.5773502691896258,        0.5773502691896258 },          // d plane  ID: 3

            {1.0e+00,                   0.0e+00,                           0.0e+00   },          // (100) junction plane     ID: 4
            {0.0e+00,                   1.0e+00,                           0.0e+00   },          // (010) junction plane     ID: 5
            {0.0e+00,                   0.0e+00,                           1.0e+00   },          // (001) junction plane     ID: 6

            {0.0e+00,                   0.7071067811865475,        0.7071067811865475 },          // (011) junction plane     ID: 7
            {0.0e+00,                   -0.7071067811865475,       0.7071067811865475 },          // (0-11) junction plane    ID: 8

            {0.7071067811865475,        0.0e+00,                    0.7071067811865475 },          // (101) junction plane     ID: 9
            {-0.7071067811865475,       0.0e+00,                    0.7071067811865475 },          // (-101) junction plane    ID: 10

            {0.7071067811865475,        0.7071067811865475,          0.0e+00           },          // (110) junction plane     ID: 11
            {-0.7071067811865475,       0.7071067811865475,          0.0e+00           }          // (-110) junction plane    ID: 12
        };

        for(i=0; i<NUM_PLANES_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, planeList[i], burgData->planeList[i]);
        }
        static int planeType[NUM_PLANES_IN_FCC] = { 0, 0, 0, 0,                           // basic glide planes
                                                    1, 1, 1, 1, 1, 1, 1, 1, 1 };          // junction planes
        memcpy(burgData->planeType,planeType,NUM_PLANES_IN_FCC*sizeof(int));

        static int numBurgsPerPlane[NUM_PLANES_IN_FCC] = { 6, 6, 6, 6,                   // basic glide planes
                                                          4, 4, 4,                      // Junction plane 1 
                                                          3, 3, 3, 3, 3, 3              // Junction plane 2
                                                        }; 

        memcpy(burgData->numBurgsPerPlane, numBurgsPerPlane,sizeof(int [NUM_PLANES_IN_FCC]));

        static int planeFirstBurgIndex[NUM_PLANES_IN_FCC] = { 0, 6, 12, 18,
                                                              24,28,32,
                                                              36,39,42,45,48,52 };  

        memcpy(burgData->planeFirstBurgIndex,planeFirstBurgIndex,sizeof(int [NUM_PLANES_IN_FCC]));

        static real8 keyLineList[NUM_KEY_LINES_IN_FCC][3] = {
           { 0.7071067811865475,    0.7071067811865475,     0.0e+00   },   // DC -- 0
           { 0.7071067811865475,    -0.7071067811865475,    0.0e+00   },   // BA -- 1
           { 0.7071067811865475,    0.0e+00   ,    0.7071067811865475 },   // DA -- 2
           { 0.7071067811865475,    0.0e+00   ,    -0.7071067811865475},   // BC -- 3
           { 0.0e+00   ,   0.7071067811865475,     0.7071067811865475 },   // DB -- 4
           { 0.0e+00   ,   0.7071067811865475,     -0.7071067811865475}    // AC -- 5
        };

        for(i=0; i<NUM_KEY_LINES_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, keyLineList[i], burgData->keyLineList[i]);
        }

        static int numPlanesPerKeyLine[NUM_KEY_LINES_IN_FCC] = { 3, 3, 3, 3, 3, 3};
 
        memcpy(burgData->numPlanesPerKeyLine,numPlanesPerKeyLine,sizeof(int [NUM_KEY_LINES_IN_FCC]));

        static int keyLineFirstPlaneIndex[NUM_BURGS_IN_FCC] = { 0, 3, 6, 9, 12, 15 }; 
        memcpy(burgData->keyLineFirstPlaneIndex,keyLineFirstPlaneIndex,sizeof(int [NUM_KEY_LINES_IN_FCC]));

        static int burgOnPlaneID[54] = {
                                          /* DC on plane (0) */
                                         0,  // a plane 
                                         1,  // b plane
                                         6, 
                                         
                                          /* BA on plane (1) */
                                         2,  // g plane
                                         3,  // d plane
                                         6, 
                                         
                                          /* DA on plane (2) */
                                         1,  // b plane
                                         2,  // g plane
                                         5, 
                                         
                                          /* BC on plane (3) */
                                         0,  // a plane
                                         3,  // d plane
                                         5, 
                                         
                                          /* DB on plane (4) */
                                         0,  // a plane
                                         2,  // g plane
                                         4, 
                                         
                                          /* AC on plane (5) */
                                         1,  // b plane
                                         3,  // d plane
                                         4, 
                                         
                                          /* Db (6) */
                                         1,
                                         8,
                                         
                                          /* dB (7) */
                                         3,
                                         8,
                                         
                                          /* gA (8) */
                                         2,
                                         7,
                                         
                                          /* aC (9) */
                                         0,
                                         7,
                                         
                                          /* Da (10) */
                                         0,
                                         10,
                                         
                                          /* gB (11) */
                                         2,
                                         9,
                                         
                                          /* dA (12) */
                                         3,
                                         10,
                                         
                                          /* bC (13) */
                                         1,
                                         9,
                                         
                                          /* Dg (14) */
                                         2,
                                         12,
                                         
                                          /* aB (15) */
                                         0,
                                         11,
                                         
                                          /* bA (16) */
                                         1,
                                         11,
                                         
                                          /* dC (17) */
                                         3,
                                         12,

                                          /* [100] (18) */
                                         5, /* junction plane */
                                         6, /* junction plane */
                                         7, /* junction plane */
                                         8, /* junction plane */
                                         
                                          /* [010] (19) */
                                         4, /* junction plane */
                                         6, /* junction plane */
                                         9, /* junction plane */
                                         10, /* junction plane */
                                         
                                           /* [001] (20) */
                                         4, /* junction plane */
                                         5, /* junction plane */
                                         11, /* junction plane */
                                         12 /* junction plane */
                                       };
        burgData->burgOnPlaneID = (int *)calloc(54, sizeof(int));        
        memcpy(burgData->burgOnPlaneID,burgOnPlaneID,sizeof(int [54]));
        
        static int planeOnBurgID[54] = { 0,  3,  4,  9, 10, 15,     // plane (0)
                                         0,  2,  5,  6, 13, 16,     // plane (1)
                                         1,  2,  4,  8, 11, 14,     //    .
                                         1,  3,  5,  7, 12, 17,     //    .
                                         4,  5, 19, 20,             //  
                                         2,  3, 18, 20,             //
                                         0,  1, 18, 19,             //
                                         8,  9, 18,                 //
                                         6,  7, 18,                 //
                                        11, 13, 19,                 //
                                        10, 12, 19,                 //
                                        15, 16, 20,                 //
                                        14, 17, 20                  //
                                       };
        burgData->planeOnBurgID = (int *)calloc(54, sizeof(int));        
        memcpy(burgData->planeOnBurgID,planeOnBurgID,sizeof(int [54]));

        static int keyLineOnPlaneID[18] = {
                                            0, 1, 6,
                                            2, 3, 6,
                                            1, 2, 5,
                                            0, 3, 5,
                                            0, 2, 4,
                                            1, 3, 4
                                          };
        burgData->keyLineOnPlaneID = (int *)calloc(18, sizeof(int));        
        memcpy(burgData->keyLineOnPlaneID, keyLineOnPlaneID,sizeof(int [54]));

        return;
}

void FreeBurgInfoStr(BurgInfo_t *burgData)
{
    free(burgData->numPlanesPerBurg); 
    burgData->numPlanesPerBurg = NULL;

    free(burgData->burgFirstPlaneIndex); 
    burgData->burgFirstPlaneIndex = NULL;

    free(burgData->burgList); 
    burgData->burgList = NULL;

    free(burgData->planeList); 
    burgData->planeList = NULL;

    free(burgData->planeType); 
    burgData->planeType = NULL;

    free(burgData->burgOnPlaneID); 
    burgData->burgOnPlaneID = NULL;

    return;
}

int CompressDirectionList(real8 **dirList, int nRows)
{
    int     i, j;
    real8   crossProdt[3];
    int     newNRows = 0;

    for(i=0; i<nRows; i++){
        for(j=0; j<newNRows; j++){
            cross(dirList[i], dirList[j], crossProdt);
            if(Normal(crossProdt) < 1.0E-2)break;

            if(j == newNRows){
                VECTOR_COPY(dirList[newNRows], dirList[i])
                newNRows++;
            }
        }
    }
    return(newNRows);
}

void GenerateGBBurgInfo(Home_t *home, GB_t *GB)
{

    Grain_t     *grain1, *grain2;
    BurgInfo_t  *bd1, *bd2;
    
    grain1 = home->grainKeys[GB->nbrGrain[0]];
    grain2 = home->grainKeys[GB->nbrGrain[1]];

    bd1 = &(grain1->burgData);
    bd2 = &(grain2->burgData);


#if 1
    
#else
    GenerateFCCBurgInfo(&(GB->burgData), grain1->orientation);
#endif

    return;
}

