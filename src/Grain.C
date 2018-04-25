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
	printf("%.15f %.15f %.15f\n", vec[0], vec[1], vec[2]);
}

void CopyMatrix33(real8 M1[3][3], real8 M2[3][3]){
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++){
			M2[i][j] = M1[i][j];
		}
	}
	return;
}

void InitOrientation(real8 MatrixOrientation[3][3], real8 TwinOrientation[3][3]){
	real8 	newz[3];
	int	i;
	
	real8 	newx[3] = {1, -1, 0};
	real8 	newy[3] = {1,  1, 1};
	
	cross(newx, newy, newz);
	NormalizeVec(newx);
	NormalizeVec(newy);
	NormalizeVec(newz);
	
	for(i = 0; i < 3; i++){
		MatrixOrientation[0][i] = newx[i];
		MatrixOrientation[1][i] = newy[i];
		MatrixOrientation[2][i] = newz[i];
		
		TwinOrientation[0][i] = - newx[i];
		TwinOrientation[1][i] = - newy[i];
		TwinOrientation[2][i] =   newz[i];
	}
//	FormatMatrix(MatrixOrientation,"MatrixOrientation");
	return;
}

void InitGrainStructure(Grain_t *grain, int numBound){
	grain->numBoundaries = numBound;
	
	grain->BoundaryIndex = (int *)calloc(grain->numBoundaries, sizeof(int));
	grain->BoundaryType = (int *)calloc(grain->numBoundaries, sizeof(int));
	
	grain->nx = (real8 *)calloc(grain->numBoundaries, sizeof(real8));
	grain->ny = (real8 *)calloc(grain->numBoundaries, sizeof(real8));
	grain->nz = (real8 *)calloc(grain->numBoundaries, sizeof(real8));
	
	grain->x0 = (real8 *)calloc(grain->numBoundaries, sizeof(real8));
	grain->y0 = (real8 *)calloc(grain->numBoundaries, sizeof(real8));
	grain->z0 = (real8 *)calloc(grain->numBoundaries, sizeof(real8));

}

void InitGrain(Home_t *home, Param_t *param){
	int		i,j;
	Grain_t		*grain;
	real8		HalfMatrixThick, maxSeg;
	real8		TwinOrientation[3][3], temR[3][3];
	int		nGrain;
	real8		twinRatio;

	/*the following variables are placed in param.h*/
	int		numTGrains;
	real8		halfTwinThick;
	int		invFlag;
	real8		lx=0, ly=0, lz=0;	/*half cube size in param.h*/

	lx = param->xBoundMax;
	ly = param->yBoundMax;
	lz = param->zBoundMax;
	maxSeg = param->maxSeg;

	grain = (Grain_t *)calloc(1, sizeof(Grain_t));

	if(param->polycrystal == 0)param->nGrain = 0;
	nGrain = param->nGrain;

	for(i=0; i<3; i++){
	    VECTOR_COPY(home->rotMatrix[i], param->matrixOrientation[i]);
	}
	Matrix33Invert(home->rotMatrix, home->rotMatrixInverse);

	if(fabs(home->rotMatrix[0][0]) > 0.9999999999 &&  
	   fabs(home->rotMatrix[1][1]) > 0.9999999999 &&
	   fabs(home->rotMatrix[2][2]) > 0.9999999999){
	    param->useLabFrame = 0;
	}else{
	    param->useLabFrame = 1;
	}


/*
 *	Type I: polycrystal = 0, nGrain = 0, PBC
 * 	origin PBC condition.
 */	

/*
 *	Type II: polycrystal = 0, nGrain = 0, Free surface.
 *	single crystal, free surface.
 */
	if(param->polycrystal ==0 ){
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

	home->newGrainKeyPtr = nGrain;				/*num of grain and twins */
	home->grainKeys	= (Grain_t **)malloc(home->newGrainKeyPtr * sizeof(Grain_t *));
	param->useLabFrame = 1;
	
/*
 *	Initilize first grain.
 */
	home->grainKeys[0] = (Grain_t *)calloc(1, sizeof(Grain_t));
	grain = home->grainKeys[0];
	InitGrainStructure(grain, 6);

	grain->type = MATRIX_GRAIN;
	grain->cy = 0;

	grain->nx[0] = -1;
	grain->x0[0] = -lx;
	grain->nx[1] =  1;
	grain->x0[1] =  lx;

	grain->ny[2] = -1;
	grain->y0[2] = -ly;
	grain->ny[3] =  1;
	grain->y0[3] =  ly;

	grain->nz[4] = -1;
	grain->z0[4] = -lz;
	grain->nz[5] =  1;
	grain->z0[5] =  lz;

	CopyMatrix33(param->matrixOrientation, grain->orientation);
	invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	if(!invFlag){
	    Fatal("Invert matrix orientation faild!\n");
	}

	if(nGrain == 1)return;

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
//	FormatVector(CDelta, "CDelta");

	Matrix33Vector3Multiply(param->matrixOrientation, planeABC, param->TBNormal);
        RotateWithAxis(param->TBNormal, PI, param->rotateNormal);
//	FormatVector(param->TBNormal, "TBNormal");

	if(param->xBoundType == 0 && param->yBoundType == 0 && 
	   param->zBoundType == 0){	
/*
 *	    Type V: polycrystal = 1, nGrain >= 2, PBC.
 *	    polycrystal, twin, RVE, XOZ plane is twin plane
 */
	    halfTwinThick = twinRatio * ly / (real8)numTGrains;
	    HalfMatrixThick = ly / (real8)numTGrains - halfTwinThick;
	    	
	    if (2.0 * halfTwinThick < maxSeg || HalfMatrixThick < maxSeg){
#if 0
	        printf("a segment may penetrate through a grain at %s line %d, max seg is %f (%f,%f)\n", 
	              __FILE__, __LINE__,maxSeg,halfTwinThick,HalfMatrixThick);
#endif
	    }
	    if (HalfMatrixThick <= 0){
	        printf("the thickness of matrix is %f\n ", HalfMatrixThick);
	        Fatal("the thickness of matrix is negative value at %s line %d\n", __FILE__, __LINE__);
	    }
    
	    for(i=1; i<nGrain; i++){
	        home->grainKeys[i] = (Grain_t *)calloc(1, sizeof(Grain_t));
	        grain = home->grainKeys[i];
	        InitGrainStructure(grain, 2);
    
	        grain->type = TWIN_GRAIN;
	        grain->cy = -ly + (2 * ((real8)i) - 1) * (HalfMatrixThick + halfTwinThick);
	        grain->ny[0] = -1;
	        grain->ny[1] =  1;
	        grain->y0[0] = grain->cy - halfTwinThick;
	        grain->y0[1] = grain->cy + halfTwinThick;
    
	        CopyMatrix33(TwinOrientation, grain->orientation);
	        invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	        if(!invFlag){
	            Fatal("Invert twin orientation faild!\n");
	        }
	        
	    }
	    return;

	}else if(nGrain == 2){
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

	        if(i == XC){
	            bound = lx;
	            hMatrix = 2 * lx * rho - lx; 
	            p[0][XC] = 0; p[0][YC] = -ly; p[0][ZC] = -lz;
	            p[1][XC] = 0; p[1][YC] = -ly; p[1][ZC] =  lz;
	            p[2][XC] = 0; p[2][YC] =  ly; p[2][ZC] = -lz;
	            p[3][XC] = 0; p[3][YC] =  ly; p[3][ZC] =  lz;
	        }else if(i == YC){
	            bound = ly;
	            hMatrix = 2 * ly * rho - ly; 
	            p[0][XC] = -lx; p[0][YC] = 0; p[0][ZC] = -lz;
	            p[1][XC] =  lx; p[1][YC] = 0; p[1][ZC] = -lz;
	            p[2][XC] = -lx; p[2][YC] = 0; p[2][ZC] =  lz;
	            p[3][XC] =  lx; p[3][YC] = 0; p[3][ZC] =  lz;
	        }else if(i == ZC){
	            bound = lz;
	            hMatrix = 2 * lz * rho - lz; 
	            p[0][XC] = -lx; p[0][YC] = -ly; p[0][ZC] = 0;
	            p[1][XC] =  lx; p[1][YC] = -ly; p[1][ZC] = 0;
	            p[2][XC] = -lx; p[2][YC] =  ly; p[2][ZC] = 0;
	            p[3][XC] =  lx; p[3][YC] =  ly; p[3][ZC] = 0;
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
	        
	    home->grainKeys[1] = (Grain_t *)calloc(1, sizeof(Grain_t));
	    grain = home->grainKeys[1];
	    InitGrainStructure(grain, 1);
	     
	    grain->type = TWIN_GRAIN;
    
	    CopyMatrix33(TwinOrientation, grain->orientation);
	    invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	    if(!invFlag){
	        Fatal("Invert twin orientation faild!\n");
	    }

	    grain->nx[0] = -param->TBNormal[XC];
	    grain->ny[0] = -param->TBNormal[YC];
	    grain->nz[0] = -param->TBNormal[ZC];

	    grain->x0[0] = plane.point[XC];
	    grain->y0[0] = plane.point[YC];
	    grain->z0[0] = plane.point[ZC];
	    
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

real8 GetInGrainFromPos(Home_t *home, real8 x, real8 y, real8 z){
	Grain_t		*grain;
	real8		ID;
	real8		nx, ny, nz, x0, y0, z0;
	real8		eps, yy, distance;
	int		i, j, numGrains;
	int		grainFlag, boundaryFlag;
	int		existTwin = 0;
	Param_t		*param;

	if(home->param->polycrystal == 0)return(0.0);
	param = home->param;
	numGrains = home->newGrainKeyPtr;

	ID = 0;
	eps = home->param->rann;
	
	if(numGrains > 1){
	    existTwin = 1;
	}

	for (i = existTwin; i < numGrains; i++){
	    if ((grain = home->grainKeys[i])==(Grain_t *)NULL){
	        continue;
	    }

	    grainFlag = 0;
	    boundaryFlag = 0;

	    for (j = 0; j < grain->numBoundaries; j++){
	        nx = grain->nx[j]; ny = grain->ny[j]; nz = grain->nz[j];
	        x0 = grain->x0[j]; y0 = grain->y0[j]; z0 = grain->z0[j];

	        yy = nx*(x - x0) + ny*(y - y0) + nz*(z - z0);
	        distance = fabs(yy);

	        if (yy <= eps){
	            grainFlag = 1;	  /* inner or on the boundary*/
	            if (distance <= eps){
	                boundaryFlag = 1; /* On the boundary : default is twin boundary */
	                                  /* some matix nodes or twin nodes are part of the category*/
	                break;
	            }
	        }else{
	            grainFlag = 0;	 /* the node is not here */
	            break;
	        }
	    }

	    if (grainFlag == 1 && boundaryFlag == 0){
	        ID = (real8)i;
	    }
	    else if (grainFlag == 1 && boundaryFlag == 1){
	        ID = ((real8)i) + BOUNDARY_NODE;
	        break;
	    }
	}

	if(ID == 0.5 && param->xBoundType == 1 &&
           param->yBoundType == 1 && param->zBoundType == 1){
	    ID = 0;
	} 
	return(ID);
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

int	BurgType_FCC(Home_t *home, real8 *burg){

/*	
 *	0: zero burgers
 */
	real8 eps = 1e-2;

	if(fabs(burg[0]) < eps &&
	   fabs(burg[1]) < eps &&
	   fabs(burg[2]) < eps){
		return(ZERO_BURG);
	}


/*
 * 	1: basic Burgers vector
 */		
	real8 burgMag;

	burgMag = Normal(burg);

	if(fabs(burgMag - 1) < eps){
		return(BASIC_BURG);
	}
/*	
 *	2: TB burgers
 */
	
	if(fabs(burgMag - 1.0/sqrt(3.0)) < eps) {
		return(TB_BURG);
	}

	return(UNKNOWN_BURG);


}

int GlidePlaneID_FCC(real8 nCryst[3]){
    int    i, nID;
    real8  crossProd[3];
    static int normalIDMap[5] = {
                                 0,           // alpha plane (BCD): [-1  1 -1]
                                 1,           // beta  plane (ADC): [ 1 -1 -1]
                                 2,           // gama  plane (BDA): [-1 -1  1]
                                 3,           // delta plane (ABC): [ 1  1  1]
                                 4
                                };

    static real8 normalList[4][3] = {
               {-1,  1, -1},
               { 1, -1, -1},
               {-1, -1,  1},
               { 1,  1,  1}

                                  };

    for(i=0; i<4; ){
        cross(normalList[i], nCryst,crossProd);
        if (fabs(DotProduct(crossProd, crossProd)) < 1.0e-3) {
            break;
        }
        i++;
    }
    nID = normalIDMap[i];
    return(nID);
}

void GetTBBurg_FCC(int nID, real8 burgTB[3]){
    if(nID == 0){
        /* Adelta*/
        burgTB[XC] = -1;
        burgTB[YC] =  2;
        burgTB[ZC] = -1;
    }else if(nID == 1){
        /* Bdelta*/
        burgTB[XC] =  2;
        burgTB[YC] = -1;
        burgTB[ZC] = -1;
    }else if (nID == 2){
        /* Cdelta*/
        burgTB[XC] = -1;
        burgTB[YC] = -1;
        burgTB[ZC] =  2;
    }else Fatal("GetTBBurg_FCC: normal ID is %d\n", nID);
    NormalizeVec(burgTB);
    burgTB[XC] /= sqrt(3.0);
    burgTB[YC] /= sqrt(3.0);
    burgTB[ZC] /= sqrt(3.0);
}

int BurgID_FCC(real8 bCryst[3]){

        int          burgIndex, burgID;
        real8        bCrystMag, crossProd[3], bEr;
        static int   numBurgs = 21;
/*
 *      Define the IDs for each type (or grouping) of burgers vector
 *      we're interested in (plus a value with which to associate any
 *      burgers vector for which we don't explicitly look.
 */
        static int   burg2IDMap[21] =
                { 0,                		/* [ 1  1  0] */
                  1,                		/* [ 1 -1  0] */
                  2,                		/* [ 1  0  1] */
                  3,                		/* [ 1  0 -1] */
                  4,                		/* [ 0  1  1] */
                  5,               		/* [ 0  1 -1] */
                  6, 7,  8,  9, 10, 11, 	/* all 12 [2 1 1] types */
                  12, 13, 14, 15, 16, 17, 
		  18, 19, 20,			/* [100] types*/
                  };             		
/*
 *      Define the list of burgers vectors for which we are explicitly
 *      looking.
 */
        static real8 burgList[21][3] = {
										/* glide burgers vectors */
                { 1.0,  1.0,  0.0},
                { 1.0, -1.0,  0.0},
                { 1.0,  0.0,  1.0},
                { 1.0,  0.0, -1.0},
                { 0.0,  1.0,  1.0},
                { 0.0,  1.0, -1.0},
										/* all 12[211] types*/
                { 2.0,  1.0,  1.0},		// Dbeta	6
                {-2.0,  1.0,  1.0},		// deltaB	7
                { 2.0, -1.0,  1.0},		// gamaA	8
                { 2.0,  1.0, -1.0},		// alphaC	9
                { 1.0,  2.0,  1.0},		// Dalpha	10	
                {-1.0,  2.0,  1.0},		// gamaB	11
                { 1.0, -2.0,  1.0},		// deltaA	12
                { 1.0,  2.0, -1.0},		// betaC	13
                { 1.0,  1.0,  2.0},		// Dgama	14
                {-1.0,  1.0,  2.0},		// alphaB	15
                { 1.0, -1.0,  2.0},		// betaA	16
		{ 1.0,  1.0, -2.0},		// deltaC	17

						/* other Junctiontypes*/
		{ 1.0,  0.0,  0.0},		//			18
                { 0.0,  1.0,  0.0},		//			19
                { 0.0,  0.0,  1.0}		//			20
	                            };	
		
/*
 *      Look through the list of reference burgers vectors for one
 *      that matches the segment burgers vector.  If we find a match
 *      return the ID associated with that burgers vector.  If no
 *      match is found, the ID defaults to the final one in the list
 *      to indicate it was some *other* burgers vector type.
 */

        for (burgIndex = 0; burgIndex < numBurgs; ) {
            cross(bCryst, burgList[burgIndex], crossProd);
            if (fabs(DotProduct(crossProd, crossProd)) < 1.0e-3) {
                burgID = burg2IDMap[burgIndex];
                break;
            }
            burgIndex++;
        }
		
        if(burgIndex == numBurgs){					/* unknow types */
            burgID = 33;				
        }
		
        if(burgID >= 6 && burgID <= 17){			/* [112] Junction types */
            bCrystMag = Normal(bCryst);
            bEr = fabs(bCrystMag - 1.0 / sqrt(3.0));
            if(bEr > 1.0e-3) burgID += 15;
        }
        
        if(burgID < 6){
            bCrystMag = Normal(bCryst);
            bEr = fabs(bCrystMag - 1.0);
            if(bEr > 1.0e-2) burgID = -burgID;
        }
 
        return(burgID);
}

int GetBurgID_FCC(Home_t *home, Node_t *node, Node_t *nbr){

        int          burgID, segID;
        int	     intGID;
        real8        bLab[3], bCryst[3];
/*
 *      If needed, rotate the provided burgers vector from the lab frame
 *      into the crystal frame.
 */
				
	segID  = GetArmID(home, node, nbr);
	ChangeOrientation(home, node, nbr,3);
        if (home->param->useLabFrame) {
            bLab[XC] = node->burgX[segID];
            bLab[YC] = node->burgY[segID];
            bLab[ZC] = node->burgZ[segID];
            Matrix33Vector3Multiply(home->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[XC] = node->burgX[segID];
            bCryst[YC] = node->burgY[segID];
            bCryst[ZC] = node->burgZ[segID];
        }
		
        burgID = BurgID_FCC(bCryst);
               
        return(burgID);
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

                segVector[XC] = node->x - nbrNode->x;
                segVector[YC] = node->y - nbrNode->y;
                segVector[ZC] = node->z - nbrNode->z;

                ZImage(home->param, &segVector[XC], &segVector[YC],
                       &segVector[ZC]);

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
void AssignGrain(Home_t *home, Node_t *node, real8 inGrain) 
{
	if(node == (Node_t *)NULL){
		printf("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		return;
	}

	node->inGrain = inGrain;

	AddOpGB(home,UPDATE_GRAIN,
		node->myTag.domainID,node->myTag.index,
		inGrain);

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

	AddOpGB(home,UPDATE_CONSTRAINT,
		node->myTag.domainID,node->myTag.index,
		(real8)constraint);

	return;
}

/****************************************************************************  
*
*	Add a topological operation to the list to be sent to 
*	remote domains for processing. Note this is only used 
*	for sending grain number of node.
*
****************************************************************************/  
void AddOpGB (Home_t *home, OpType_t type, int dom1, int idx1, 
			  real8 inGrain) 
{
	Operate_t *op;

		/*
		 *          Make sure the buffer allocated for the operation list is
		 *          large enough to contain another operation.  If not, increase
		 *          the size of the buffer.
		 */
		if (home->OpCount >= home->OpListLen) {
			ExtendOpList (home);
		}
        
		op = &home->opList[home->OpCount];

		op->type = type;
		op->dom1 = dom1;
		op->idx1 = idx1;
		op->inGrain = inGrain;

		home->OpCount++;

	return;
}

real8 InWhichGrain(Home_t *home, real8 x, real8 y, real8 z){
	Grain_t		*grain;
	real8		ID, eps = 0;
	real8		nx, ny, nz, x0, y0, z0;
	real8		yy, distance;
	int		i, j, numGrains, numBoundaries;
	int		grainFlag;

	numGrains = home->newGrainKeyPtr;
	
	if(home->param->polycrystal == 0)return(0);

	ID = 0;
	
//	FoldBox(home->param, &x, &y, &z);
	for (i = 1; i < numGrains; i++){
	    if ((grain = home->grainKeys[i])==(Grain_t *)NULL){
	        continue;
	    }

	    grainFlag = 0;
	    numBoundaries = grain->numBoundaries;

	    for (j = 0; j < numBoundaries; j++){
	        nx = grain->nx[j]; ny = grain->ny[j]; nz = grain->nz[j];
	        x0 = grain->x0[j]; y0 = grain->y0[j]; z0 = grain->z0[j];
	        	
	        yy = nx*(x - x0) + ny*(y - y0) + nz*(z - z0);
	        distance = fabs(yy);

	        if (yy < -eps){
	            grainFlag = 1;				
	        }
	        else if(yy > eps){
	            grainFlag = 0;						/* the node is not here */
	            break;
	        }else{
	            return(-1);
	        }
	    }

	    if (grainFlag == 1){
	        ID = (real8)i;
	        return(ID);
	    }
	}
	return(0);
}


int NodeOnGB(Home_t *home, Node_t *node){
	int 	i, numNbrs, burgType, type = 0;
	real8	burg[3];
	Node_t 	*nbr;
	
	if(node == (Node_t *)NULL){
		Fatal("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);

	}

	if(node->inGrain - floor(node->inGrain) == 0)return(INGRAIN_NODE);

	numNbrs = node->numNbrs;
	
	for(i=0; i<numNbrs; i++){
		burg[0] = node->burgX[i];
	    burg[1] = node->burgY[i];
	    burg[2] = node->burgZ[i];
		
	    burgType = BurgType_FCC(home,burg);
		
	    if(burgType == TB_BURG) type |= TB_NODE;
	    if(burgType == BASIC_BURG) type |= INTERSECT_NODE;
	    if(burgType == UNKNOWN_BURG) type |= RESIDUAL_NODE;
	}
	return(type);
}

real8 CorrectInGrain(Home_t *home, Node_t *node, real8 inGrainPos){
	real8 inGrain;
	
	if(node == (Node_t *)NULL){
		Fatal("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		
	}
	
	if(inGrainPos == 0.5) return(MATRIX_NODE);
/*
 *	Before CorrectInGrain, we must ensure that the node which inGrain equal	to
 *	n is a matrix or twin node. 
 */		
	if(inGrainPos - floor(inGrainPos) == 0){
		return(inGrainPos);
	}
	
	if(NodeOnGB(home, node) & TB_NODE) {
		return(inGrainPos);
	}else if(NodeOnGB(home, node) & INTERSECT_NODE
			 && node->numNbrs == 2){
		Node_t *nbr1,*nbr2;
		nbr1 = GetNeighborNode(home, node, 0);
		nbr2 = GetNeighborNode(home, node, 1);
		if(nbr1->inGrain != nbr2->inGrain &&
		   nbr1->inGrain - floor(nbr1->inGrain) ==0 &&
		   nbr2->inGrain - floor(nbr2->inGrain) ==0){
			return(inGrainPos);
		}else{
			Fatal("the node is not a INTERSECT_NODE");
		}
	}else{
		real8 		vec[3], plane[3], vecDotPlane;
		Grain_t		*grain;
		int			i, numBoundaries;
		
		grain = home->grainKeys[(int)inGrainPos];
		
		numBoundaries = grain->numBoundaries;
		for(i=0; i<numBoundaries; i++){
			vec[XC] = node->x - grain->x0[i];
			vec[YC] = node->y - grain->y0[i];
			vec[ZC] = node->z - grain->z0[i];
			plane[XC] = grain->nx[i];
			plane[YC] = grain->ny[i];
			plane[ZC] = grain->nz[i];
			vecDotPlane = DotProduct(vec, plane);
			if(fabs(vecDotPlane) > TWIN_LAYER_RADIUS) continue;
			if(vecDotPlane < 0){
				return(floor(inGrainPos));
			}else if(vecDotPlane > 0){
				return(MATRIX_NODE);
			}else{
				return(node->inGrain);
			}
		}
	}
	
	Fatal(" Boundary node (%d,%d) inGrain (%f) is wrong",
		  node->myTag.domainID,node->myTag.index,node->inGrain);
	
}

real8 ChangeOrientation(Home_t *home, Node_t *node, Node_t *nbr, int flag){
	
	static int currInGrain = -999;
	real8	inGrain1, inGrain2, armInGrain;
	real8	inGrainPos;
	real8	cuurInGrain;
	
	if(!home->param->polycrystal) {
            return(0.0);
	}
	
	
	if(node == (Node_t *)NULL){
            Fatal("Task %d: %s, line %d -- Node = NULL",
                  home->myDomain, __FILE__, __LINE__);
	}
	
	
	inGrain1 = node->inGrain;
	inGrain2 = nbr->inGrain;	
	
	if(inGrain1 == inGrain2){
            armInGrain = inGrain1;
	}else{
            if(inGrain1 - floor(inGrain1) == 0 && inGrain2 - floor(inGrain2) == 0.5){
                armInGrain = inGrain1;
            }else if(inGrain2 - floor(inGrain2) == 0 && inGrain1 - floor(inGrain1) == 0.5){
                armInGrain = inGrain2;
            }else if(inGrain1 - floor(inGrain1) == 0 && inGrain2 - floor(inGrain2) == 0){	
	        if(SurfaceNode(node) && SurfaceNode(nbr)){
	            armInGrain = node->inGrain;
//		    printf("Surface segment: (%d,%d)--(%d,%d)\n",node->myTag.domainID,node->myTag.index,
//	                   nbr->myTag.domainID, nbr->myTag.index);
	        }else{
		    PrintNode(node);
		    PrintNode(nbr);
                    Fatal("ChangeOrienatation: Wrong inGrain\n");
	        }
	    }else{
                Fatal(" ChangeOrientation: wrong inGrain, Node %d's inGrain is %f, neighbor %d's inGrain is %f.",
                      node->myTag.index,node->inGrain,nbr->myTag.index,nbr->inGrain);
	    }
	}
	if(flag == 0)return(armInGrain);
	
	
	
	Grain_t	*grain;
	grain = home->grainKeys[(int)armInGrain];
	
	if(flag == 1){
	    CopyMatrix33(grain->orientationInv,home->rotMatrixInverse);
	}
	else if(flag == 2){
	    CopyMatrix33(grain->orientation, home->rotMatrix);
	}
	else if(flag == 3){
	    if(((int)armInGrain) == currInGrain)return(armInGrain);
	    currInGrain = (int)armInGrain;
	    	
	    CopyMatrix33(grain->orientationInv,home->rotMatrixInverse);
	    CopyMatrix33(grain->orientation, home->rotMatrix);
	}
	
	return(armInGrain);
}


real8 ClosedGB(Home_t *home, real8 *pos, Plane_t *GBPlane){
	real8	tempDistance;
	int	existTB = 0;
	real8 	minDistance = 1e14;
	
	FoldBox(home->param, pos, pos+1, pos+2);
	
	if(home->param->nGrain > 1)existTB = 1;
	if(Normal(GBPlane->normal) < 1.0e-3){	// the GB was not defined, so should define the nearest GB. 
	    int			i,j;
	    Grain_t 	*grain;
	    Plane_t 	plane;
	    	
	    for(i=existTB; i<home->newGrainKeyPtr; i++){
	        grain = home->grainKeys[i];
	        if(grain->type == MATRIX_GRAIN)continue;	// skip matrix grain.
	        for(j=0; j<grain->numBoundaries; j++){
	            plane.normal[XC] = grain->nx[j];
	            plane.normal[YC] = grain->ny[j];
	            plane.normal[ZC] = grain->nz[j];
	            plane.point[XC] = grain->x0[j];
	            plane.point[YC] = grain->y0[j];
	            plane.point[ZC] = grain->z0[j];
	            	
	            tempDistance = PointPlaneDistance(pos, &plane);
	            if(fabs(tempDistance) < fabs(minDistance)){
	                minDistance = tempDistance;
	                VECTOR_COPY(GBPlane->normal,plane.normal);
	                VECTOR_COPY(GBPlane->point,plane.point);
	            }
	        }
	    }
		
	}else{
	    minDistance = PointPlaneDistance(pos, GBPlane);
	}
	return(minDistance);
}

/****************************************************************************  
*
*	The GB nodes may glide out of GB due to dislocation climb and 
*	computational error, so here we modify the coordinate to ensure node
*	on GB.
*
*****************************************************************************/
void PullToGB(Home_t *home, Node_t *node, real8 oldInGrain){

	
	real8		eps, pos[3], minDistance, newPos[3],vec[3],oldPos[3];
	real8 		sb, currInGrain, changeGrain,n1[3],n2[3];
	int		crossSlipNode = 0;
	Plane_t	 	GBPlane,oldPlane;
	Line_t		path;
	Grain_t 	*oldGrain;
	
	if(node == (Node_t *)NULL){
            printf("Task %d: %s, line %d -- Node = NULL",
                   home->myDomain, __FILE__, __LINE__);
            return;
	}
	
	eps = home->param->rann;
	VECTOR_ZERO(GBPlane.normal);
	VECTOR_ZERO(GBPlane.point);
	VECTOR_ZERO(path.point);
	VECTOR_ZERO(path.direction);
	
	pos[XC] = node->x; pos[YC] = node->y; pos[ZC] = node->z;
	FoldBox(home->param, pos, pos+1, pos+2);
	currInGrain = InWhichGrain(home, pos[XC], pos[YC], pos[ZC]);

	minDistance = ClosedGB(home, pos, &GBPlane);
	
	if(oldInGrain - floor(oldInGrain) == 0.5){
            if(fabs(minDistance) > 2.0*eps){
                printf("warning: the GB node (%d,%d) diviate far from GB, minDistance is %f\n"
                       ,node->myTag.domainID,node->myTag.index,
                       minDistance);
            }
            newPos[XC] = pos[XC] - GBPlane.normal[XC] * minDistance;
            newPos[YC] = pos[YC] - GBPlane.normal[YC] * minDistance;
            newPos[ZC] = pos[ZC] - GBPlane.normal[ZC] * minDistance;
	
	    if(!(NodeOnGB(home, node) & TB_NODE)){
	        newPos[XC] =  node->oldx;
	        newPos[YC] =  node->oldy;
	        newPos[ZC] =  node->oldz;
	    }

	}else{
            if(fabs(minDistance) > 1.0e1*eps){
                printf("warning: the grain node (%d,%d,%3.1f) diviate far from GB, minDistance is %f\n"
                       ,node->myTag.domainID,node->myTag.index,
                       node->inGrain,minDistance);
            }
		
            if(currInGrain != -1.0){
                if(currInGrain != oldInGrain){
                    changeGrain = 1.0;
                }else{
                    changeGrain = -1.0;
                }
			
                sb = ((minDistance < 0) ? -1.0 : 1.0);
                sb *= changeGrain;
            }else{
                oldGrain = home->grainKeys[(int)oldInGrain];
                sb = ((oldGrain->type == TWIN_GRAIN) ? 1.0 : -1.0);
            }
		
            GBPlane.point[XC] -= GBPlane.normal[XC] * eps * sb;
            GBPlane.point[YC] -= GBPlane.normal[YC] * eps * sb;
            GBPlane.point[ZC] -= GBPlane.normal[ZC] * eps * sb;
		
	    if(node->numNbrs == 2){
	        if(fabs(node->nx[0] * node->nx[1] +
	                node->ny[0] * node->ny[1] +
	                node->nz[0] * node->nz[1]) < 0.995)crossSlipNode = 1; 
	    } 

            if(crossSlipNode){
	        n1[XC] = node->nx[0];
	        n1[YC] = node->ny[0];
	        n1[ZC] = node->nz[0];

	        n2[XC] = node->nx[1];
	        n2[YC] = node->ny[1];
	        n2[ZC] = node->nz[1];

	        cross(n1,n2,path.direction);
	        NormalizeVec(path.direction);
	        VECTOR_COPY(path.point, pos);

	        if(!IntersectionOfLinePlane(&path, &GBPlane, newPos)){
                    Fatal("PullToGB: IntersectionOfLinPlane faild");
	        }
            }else{
                VECTOR_COPY(oldPlane.point, pos);
                oldPlane.normal[XC] = node->nx[0];
                oldPlane.normal[YC] = node->ny[0];
                oldPlane.normal[ZC] = node->nz[0];
                	
                if(!IntersectionOfPlanePlane(&oldPlane,&GBPlane,&path)){
		    PrintNode(node);
               //     Tecplot(home, "PullToGB", 0, 0, 0, 0, 0);
                    Fatal("PullToGB: IntersectionOfPlanePlane faild");
                }
			
                VECTOR_COPY(newPos, path.point);
            }

            if(oldInGrain != InWhichGrain(home, newPos[XC],newPos[YC],newPos[ZC]) ){
                PrintNode(node);
                Fatal("PullToGB: the new position is not in grain %3.1f",oldInGrain);
            }
	}

	if(Normal(newPos) > 1E5){
	    PrintNode(node);
	    PrintNeighbor(home,node);
            Fatal("PullToGB: the new posion is (%f,%f,%f)",newPos[XC], newPos[YC], newPos[ZC]);
	}
	
	FoldBox(home->param,&newPos[0],&newPos[1],&newPos[2]);
	RepositionNode(home,newPos,&node->myTag,1);

}

/****************************************************************************  
*
*	In order to avoid numerical error, we must ensure the node is pulled 
*	back into its grain.
*
*****************************************************************************/
void PullBack(Home_t *home, Node_t *node, real8 oldInGrain){

	real8	inGrainPos;
	real8	newPos[3];
	
	if(node == (Node_t *)NULL){
            printf("Task %d: %s, line %d -- Node = NULL",
		   home->myDomain, __FILE__, __LINE__);
            return;
	}
	
	newPos[0] = node->oldx;
	newPos[1] = node->oldy;
	newPos[2] = node->oldz;
	
	FoldBox(home->param, newPos, newPos+1, newPos+2);
	
#if 0
/*
 *	ensure that the old position is in old grain. 
 */
	real8 	checkInGrain;
	
	if(oldInGrain - floor(oldInGrain) == 0.5){
		checkInGrain = GetInGrainFromPos(home, newPos[XC], newPos[YC], newPos[ZC]);
		
	}else{
		checkInGrain = InWhichGrain(home, newPos[XC], newPos[YC], newPos[ZC]);
	}
	
	if(checkInGrain != node->inGrain){
		PrintNode(node);
		Fatal("CHECK_PULLBACK: checkInGrain (%3.1f) != oldInGrain.",checkInGrain);
	}
#endif
	
	
	inGrainPos = GetInGrainFromPos(home,newPos[0],newPos[1],newPos[2]);

	if(inGrainPos != oldInGrain && node->numNbrs == 2){
            PullToGB(home, node, oldInGrain);
	}else{
            RepositionNode(home,newPos,&node->myTag,1);
	}

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
	
	vec[XC] = point[XC] - plane->point[XC];
	vec[YC] = point[YC] - plane->point[YC];
	vec[ZC] = point[ZC] - plane->point[ZC];
	
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

int ClosestPointOfPointLine(real8 *p, Line_t *line, real8 *outP){
	if(NodeOnLine(p, line)){
		VECTOR_COPY(outP, p);
		return(0);
	}
	
	real8 t,dotPLD, dotLPLD;
	
	dotPLD = DotProduct(p, line->direction);
	dotLPLD = DotProduct(line->point, line->direction);
	
	t = dotPLD - dotLPLD;
	
	outP[XC] = line->point[XC] + line->direction[XC] * t;
	outP[YC] = line->point[YC] + line->direction[YC] * t;
	outP[ZC] = line->point[ZC] + line->direction[ZC] * t;
	
	return(1);

}


int IntersectionOfPlanePlane(Plane_t *plane1, Plane_t *plane2, Line_t *line){
	real8 dotNormal;
	
	dotNormal = DotProduct(plane1->normal, plane2->normal);
	if(fabs(1.0-fabs(dotNormal)) < 1.0e-2){
	    printf("warning: the two planes has a same normal vector. at %s, line %d \n",__FILE__,__LINE__);
	    return(0);
	}
	
	NormalizedCrossVector(plane1->normal, plane2->normal, line->direction);
	
	real8	p1x,p1y,p1z,n1x,n1y,n1z;
	real8	p2x,p2y,p2z,n2x,n2y,n2z;
	real8	dx,dy,dz, denominator;
	
	p1x = plane1->point[XC]; p1y = plane1->point[YC]; p1z = plane1->point[ZC]; 
	p2x = plane2->point[XC]; p2y = plane2->point[YC]; p2z = plane2->point[ZC]; 
	n1x = plane1->normal[XC]; n1y = plane1->normal[YC]; n1z = plane1->normal[ZC]; 
	n2x = plane2->normal[XC]; n2y = plane2->normal[YC]; n2z = plane2->normal[ZC]; 
	dx = line->direction[XC]; dy = line->direction[YC]; dz = line->direction[ZC];
	
	denominator = (dz*n1y*n2x - dy*n1z*n2x - dz*n1x*n2y + dx*n1z*n2y + dy*n1x*n2z - 
				   dx*n1y*n2z);
	if(denominator == 0){
		FormatVector(plane1->normal,"plane1->normal");
		FormatVector(plane1->point,"plane1->point");
		FormatVector(plane2->normal,"plane2->normal");
		FormatVector(plane2->point,"plane2->point");
		FormatVector(line->direction,"line->direction");
		
		printf("warning: denominator is zero. at %s, line %d \n",__FILE__,__LINE__);
		VECTOR_COPY(line->direction,plane1->normal);
		IntersectionOfLinePlane(line, plane2, line->point);
		return(0);
	}
	
	line->point[XC] = -(dz*n1x*n2y*p1x - dx*n1z*n2y*p1x - dy*n1x*n2z*p1x + dx*n1y*n2z*p1x + 
					   dz*n1y*n2y*p1y - dy*n1z*n2y*p1y + dz*n1y*n2z*p1z - 
					   dy*n1z*n2z*p1z - dz*n1y*n2x*p2x + dy*n1z*n2x*p2x - 
					   dz*n1y*n2y*p2y + dy*n1z*n2y*p2y - dz*n1y*n2z*p2z + dy*n1z*n2z*p2z)/denominator;
	
	line->point[YC] = -((-dz)*n1x*n2x*p1x + dx*n1z*n2x*p1x - dz*n1y*n2x*p1y + 
					   dy*n1z*n2x*p1y - dy*n1x*n2z*p1y + dx*n1y*n2z*p1y - dz*n1x*n2z*p1z + 
					   dx*n1z*n2z*p1z + dz*n1x*n2x*p2x - dx*n1z*n2x*p2x + 
					   dz*n1x*n2y*p2y - dx*n1z*n2y*p2y + dz*n1x*n2z*p2z - dx*n1z*n2z*p2z)/denominator;
	
	line->point[ZC] = -(dy*n1x*n2x*p1x - dx*n1y*n2x*p1x + dy*n1x*n2y*p1y - dx*n1y*n2y*p1y - 
					   dz*n1y*n2x*p1z + dy*n1z*n2x*p1z + dz*n1x*n2y*p1z - 
					   dx*n1z*n2y*p1z - dy*n1x*n2x*p2x + dx*n1y*n2x*p2x - 
					   dy*n1x*n2y*p2y + dx*n1y*n2y*p2y - dy*n1x*n2z*p2z + dx*n1y*n2z*p2z)/denominator;
	
	return(1);
}

real8 GetArmInGrain(real8 inGrain1, real8 inGrain2){
	if(inGrain1 == inGrain2)return(inGrain1);
	if(inGrain1 - floor(inGrain1) == 0 && 
			inGrain2 - floor(inGrain2) == 0.5){
		return(inGrain1);
	}else if(inGrain2 - floor(inGrain2) == 0 && 
			  inGrain1 - floor(inGrain1) == 0.5){
		return(inGrain2);
	}else{
		Fatal("GetArmInGrain: the arm penetrate through GB.");
	}
}

void ClosetPointToGB(Param_t *param, Plane_t *plane, Plane_t *GBPlane, real8 *p){
	 Line_t  line;

	VECTOR_ZERO(line.direction);
	VECTOR_ZERO(line.point);
	
	IntersectionOfPlanePlane(plane, GBPlane, &line);

	VECTOR_COPY(p,line.point);
	FoldBox(param,p,p+1,p+2);
	
}

int	IsZeroVector(real8 *vec){
	return(fabs(vec[XC]) < 1.0e-10 && 
		   fabs(vec[YC]) < 1.0e-10 && 
		   fabs(vec[ZC]) < 1.0e-10);
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

	nx = (position[0] + home->param->Lx/2)/eleLx;
	ny = (position[1] + home->param->Ly/2)/eleLy;
	nz = (position[2] + home->param->Lz/2)/eleLz;

	if(nx == param->meshNumX) nx--;
	if(ny == param->meshNumY) ny--;
	if(nz == param->meshNumX) nz--;

	return(nz*param->meshNumX*param->meshNumY + ny*param->meshNumX +nx);

}

void AdjustSurfaceNodeVel(Param_t *param, Node_t *node)       //adjust surface node volecity
{
	if(node->constraint > 3 || node->constraint < 1) return;

	real8 nSfc[3],nPln[3],dir[3],resLen;

	nSfc[0] = node->constraint==1;
	nSfc[1] = node->constraint==2;
	nSfc[2] = node->constraint==3;

	nPln[0] = node->nx[0];
	nPln[1] = node->ny[0];
	nPln[2] = node->nz[0];


	cross(nSfc,nPln,dir);

	if(Normal(dir) < 1e-5)return;

	NormalizeVec(dir);

	resLen = node->vX * dir[0] +
		node->vY * dir[1] +
		node->vZ * dir[2] ;

	node->vX = resLen * dir[0];
	node->vY = resLen * dir[1];
	node->vZ = resLen * dir[2];
	
	if(node->inGrain - floor(node->inGrain) == 0.5){
	    if(fabs(DotProduct(dir, param->TBNormal)) > 0.001){
	        node->vX = 0.0;
	        node->vY = 0.0;
	        node->vZ = 0.0;
	        return;
	    }

	}

	if(node->numNbrs > 1){
	    int 	i;
	    real8	nDir[3];

	    for(i=1; i<node->numNbrs; i++){
	        nDir[XC] = node->nx[i];
	        nDir[YC] = node->ny[i];
	        nDir[ZC] = node->nz[i];
	        if(fabs(DotProduct(dir, nDir) > 0.001)){
	            node->vX = 0.0;
	            node->vY = 0.0;
	            node->vZ = 0.0;
	            return;
	        }
	    }
	}

//	printf("surface node velocitiy is %f, origin is %f\n",resLen * param->burgMag);

//if(fabs(fabs(node->x) - param->Lx/2.0) < 20) node->vX = 0.0;  //gb surface nodes
}

#if 0
void   AdjustSurfaceNode(Home_t *home)       //adjust surface node volecity
	{


	int 	iii,i,j;
	real8	currx,curry,currz, nodePos[3], nbrPos[3];
	real8	interPos[3], oldPos[3], distance[3];
	real8	maxDis;
	int	constraintSurface;
 	Param_t *param;
	Node_t  *node, *nbr, *newNode;
	Line_t	line;
	Plane_t	plane;
	param = home->param;
  
	for (iii = 0; iii < home->newNodeKeyPtr; iii++){

        if ((node = home->nodeKeys[iii]) == (Node_t *)NULL) continue;

	    currx=node->x; curry=node->y; currz=node->z;

	    VECTOR_ZERO(distance);
	    if(param->xBoundType == Free &&  fabs(currx)>param->xBoundMax) { 
				 
	        distance[XC] = fabs(currx) - param->xBoundMax;
//		node->x=node->oldx; node->y=node->oldy; node->z=node->oldz;

	    }

	    if(param->yBoundType == Free &&  fabs(curry)>param->yBoundMax) { 
			 
	        distance[YC] = fabs(curry) - param->yBoundMax;
//		node->x=node->oldx; node->y=node->oldy; node->z=node->oldz;

	    }

	    if(param->zBoundType == Free &&  fabs(currz)>param->zBoundMax) { 
			
	        distance[ZC] = fabs(currz) - param->zBoundMax;
//		node->x=node->oldx; node->y=node->oldy; node->z=node->oldz;

	    }
	    
	    maxDis = 0;
	    for(i=0; i<3; i++){
	        if(distance[i] > maxDis){
	            maxDis = distance[i];
	            constraintSurface = i;
	        }
	    }

	    if(maxDis > 0)node->constraint = constraintSurface;
/*
	    if(node->constraint >= 1 && node->constraint <= 3) continue;

	    if(node->numNbrs != 2) continue;

	    int jjj,surfaceNbr;
	    Node_t *nbr;

	    surfaceNbr = 0;

	    for(jjj=0;jjj<node->numNbrs;jjj++){

		if((nbr=GetNeighborNode(home,node,jjj))==(Node_t *)NULL) continue;

		if(nbr->constraint >= 1 && nbr->constraint <= 3) surfaceNbr++;

		if(surfaceNbr == 2) {
		    node->constraint = nbr->constraint;
		    break;
		}

	    }
*/

	}

	node = home->ghostNodeQ;

        while (node) {

	    currx=node->x; curry=node->y; currz=node->z;

	    VECTOR_ZERO(distance);
	    if(param->xBoundType == Free &&  fabs(currx)>param->xBoundMax) { 
				 
	        distance[XC] = fabs(currx) - param->xBoundMax;
//		node->x=node->oldx; node->y=node->oldy; node->z=node->oldz;

	    }

	    if(param->yBoundType == Free &&  fabs(curry)>param->yBoundMax) { 
			 
	        distance[YC] = fabs(curry) - param->yBoundMax;
//		node->x=node->oldx; node->y=node->oldy; node->z=node->oldz;

	    }

	    if(param->zBoundType == Free &&  fabs(currz)>param->zBoundMax) { 
			
	        distance[ZC] = fabs(currz) - param->zBoundMax;
//		node->x=node->oldx; node->y=node->oldy; node->z=node->oldz;

	    }
	    
	    maxDis = 0;
	    for(i=0; i<3; i++){
	        if(distance[i] > maxDis){
	            maxDis = distance[i];
	            constraintSurface = i;
	        }
	    }

	    if(maxDis > 0)node->constraint = constraintSurface;

            node = node->next;
	}

	for (i = 0;i < home->newNodeKeyPtr; i++){

        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
	    if(node->constraint < 1 || node->consraint > 3)continue;
	    nodePos[XC] = node->x;
	    nodePos[YC] = node->y;
	    nodePos[ZC] = node->z;

	    if(param->xBoundType == Free &&  fabs(nodePos[XC]) < param->xBoundMax + param->rann/2 &&
	       param->yBoundType == Free &&  fabs(nodePos[YC]) < param->yBoundMax + param->rann/2 &&
	       param->zBoundType == Free &&  fabs(nodePos[ZC]) < param->zBoundMax + param->rann/2)continue;

	    if(InWhichGrain(home, nodePos[XC], nodePos[YC], nodePos[ZC]) != node->inGrain &&
	       node->inGrain - floor(node->inGrain) == 0){
	        oldPos[XC] = node->oldx;
	        oldPos[YC] = node->oldy;
	        oldPos[ZC] = node->oldz;
	        RepositionNode(home,oldPos,&node->myTag,1);
	        continue;
	    }

	    for(j=0; j<node->numNbrs; j++){
	        nbr = GetNeiborNode(home, node, j);
	        if(nbr == (Node_t *)NULL)continue;
	        if(nbr->constraint >= 1 &&  nbr->constraint <= 3)continue;
	        nbrPos[XC] = nbr->x;
	        nbrPos[YC] = nbr->y;
	        nbrPos[ZC] = nbr->z;
	        if(param->xBoundType == Free &&  fabs(nbrPos[XC]) > param->xBoundMax - param->rann/2 ||
	           param->yBoundType == Free &&  fabs(nbrPos[YC]) > param->yBoundMax - param->rann/2 ||
	           param->zBoundType == Free &&  fabs(nbrPos[ZC]) > param->zBoundMax - param->rann/2)continue;
		
	        if(InWhichGrain(home, nbrPos[XC], nbrPos[YC], nbrPos[ZC]) != nbr->inGrain &&
	           nbr->inGrain - floor(nbr->inGrain) == 0){
	            oldPos[XC] = nbr->oldx;
	            oldPos[YC] = nbr->oldy;
	            oldPos[ZC] = nbr->oldz;
	            RepositionNode(home,oldPos,&nbr->myTag,1);
	            continue;
	        }

		VECTOR_COPY(line.point, nbrPos);
	        line.direction[XC] = nodePos[XC] - nbrPos[XC];
	        line.direction[YC] = nodePos[YC] - nbrPos[YC];
	        line.direction[ZC] = nodePos[ZC] - nbrPos[ZC];
	        if(Normal(line.direction) < 1.0e-8){
	            Fatal("AdjustSurfaceNode: wrong direction.");
		}

	        NormalizeVec(line.direction);		

	        plane.normal[XC] = node->constraint==1;
	        plane.normal[YC] = node->constraint==2;
	        plane.normal[ZC] = node->constraint==3;

	        VECTOR_ZERO(plane.point);

	        if(node->constraint == 1){
	            plane.point[XC] = Sign(nodePos[XC]) * param->xBoundMax;
	        }
	        else if(node->constraint == 2){
	            plane.point[YC] = Sign(nodePos[YC]) * param->yBoundMax;
	        }
	        else if(node->constraint == 3){
	            plane.point[ZC] = Sign(nodePos[ZC]) * param->zBoundMax;
	        }

	        if(!IntersectionOfLinePlane(&line, &plane, interPos)){
                    Fatal("AdjustSurfaceNode: intesectionfaild");
                }

	        if(InWhichGrain(home, interPos[XC], interPos[YC], interPos[ZC]) != node->inGrain &&
	           node->inGrain - floor(node->inGrain) == 0){
	            Fatal("AdjustSurfaceBode: at %s, Line %d",__FILE__,__LINE__);
	            continue;
	        }

	        if((newNode = InsertNode(home,node,nbr,interPos,2,1)) == (Node_t *)NULL) {
	           Fatal("AdjustSurfaceNode: Insert node failed");
	        }

	        AssignConstraint(home, newNode, node->constraint);

	        MarkNodeForceObsolete(home, newNode);
	        MarkNodeForceObsolete(home, node);
	        MarkNodeForceObsolete(home, nbr);
	
	        newNode->flags |= NO_COLLISIONS;

	        newNode->vX = 0.0;
	        newNode->vY = 0.0;
	        newNode->vZ = 0.0;

	    } /*for(j=0; j< ...)*/
	}/*for(i=0; i< ...)*/	
	return;
}

#endif
void PrintTag(Node_t *node){

        if(node == (Node_t *)NULL){
            printf("(NULL)");
        }else{
            printf("[%d %d %d %3.1f]-{%e,%e,%e} ",node->myTag.domainID,
                   node->myTag.index,node->numNbrs,node->inGrain,
	           node->x,node->y,node->z);
	}
		

}

void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]){
    NormalizeVec(a);
    real8 u,v,w;
    
    u = a[XC];
    v = a[YC];
    w = a[ZC];
    R[XC][XC] = (1.0 - cos(theta))*u*u + cos(theta);
    R[XC][YC] = u*v*(1.0 - cos(theta)) - w*sin(theta);
    R[XC][ZC] = u*w*(1.0 - cos(theta)) + v*sin(theta);
    R[YC][XC] = u*v*(1.0 - cos(theta)) + w*sin(theta);
    R[YC][YC] = (1.0 - cos(theta))*v*v + cos(theta);
    R[YC][ZC] = v*w*(1.0 - cos(theta)) - u*sin(theta);
    R[ZC][XC] = u*w*(1.0 - cos(theta)) - v*sin(theta);
    R[ZC][YC] = v*w*(1.0 - cos(theta)) + u*sin(theta);
    R[ZC][ZC] = (1.0 - cos(theta))*w*w + cos(theta);
    return;
}

int SurfaceNode(Node_t *node){
	if(node->constraint >= 1 && node->constraint <= 3)return(1);
	return(0);
}

int ExceedBox(Param_t *param, real8 *pos){
	int flags = 0x0000;
        if (param->xBoundMax + 1.0E-6 < fabs(pos[XC]))flags |= EXCEEDX;
        if (param->yBoundMax + 1.0E-6 < fabs(pos[YC]))flags |= EXCEEDY;
        if (param->zBoundMax + 1.0E-6 < fabs(pos[ZC]))flags |= EXCEEDZ;
	return(flags);
}

int IntersectionOfLineBox(Param_t *param, Line_t *line, real8 interPos[2][3]){
	int	i,j,k;
	real8	cof;
	Plane_t plane;
	real8	box[3];

	box[XC] = param->xBoundMax;
	box[YC] = param->yBoundMax;
	box[ZC] = param->zBoundMax;

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
 *          f1        Array in which to return osmotic force at pos1
 *          f2        Array in which to return osmotic force at pos2
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
