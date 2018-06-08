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
	int		    nGrain, typei, GBID;
	real8		twinRatio;

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
	Matrix33Invert(home->rotMatrix, home->rotMatrixInverse);

	if(fabs(home->rotMatrix[0][0]) > 0.9999999999 &&  
	   fabs(home->rotMatrix[1][1]) > 0.9999999999 &&
	   fabs(home->rotMatrix[2][2]) > 0.9999999999){
	    param->useLabFrame = 0;
	}else{
	    param->useLabFrame = 1;
	}

    if(param->polycrystal == 0){
        if(param->xBounType == Periodic && 
           param->yBounType == Periodic && 
           param->zBounType == Periodic ){
            type = 1;
        }else type = 2;
    }else{
        if(param->nGrain == 1){
            if(param->xBounType == Periodic && 
               param->yBounType == Periodic && 
               param->zBounType == Periodic ){
                type = 3;
            }else{
                Fatal("at file %s line %d", __FILE__, __LINE__);
            }
        }else{
            if(param->xBounType == Periodic && 
               param->yBounType == Periodic && 
               param->zBounType == Periodic ){
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
        }

	    for(i=0; i<nGrain; i++){
	        home->grainKeys[i] = (Grain_t *)calloc(1, sizeof(Grain_t));
	        grain = home->grainKeys[i];
	        InitGrainStructure(grain, 2);
    
	        grain->type = (i==0) ? MATRIX_GRAIN : TWIN_GRAIN;
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

            if(grain->type == TWIN_GRAIN){
	            CopyMatrix33(TwinOrientation, grain->orientation);
            }else{
                CopyMatrix33(param->matrixOrientation, grain->orientation);
            }

	        invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	        if(!invFlag){
	            Fatal("Invert twin orientation faild!\n");
	        }
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

	    home->grainKeys[0] = (Grain_t *)calloc(1, sizeof(Grain_t));
        grain = home->grainKeys[0];
	    InitGrainStructure(grain, 1);
	    grain->type = TWIN_GRAIN;
        grain->grainID = 0;
        grain->sign[0] = 1.0;
        grain->GB[0] = home->GBKeys[0];
	    CopyMatrix33(param->matrixOrientation, grain->orientation);
	    invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	    if(!invFlag){
	        Fatal("Invert twin orientation faild!\n");
	    }

	    home->grainKeys[1] = (Grain_t *)calloc(1, sizeof(Grain_t));
	    grain = home->grainKeys[1];
	    InitGrainStructure(grain, 1);
	    grain->type = TWIN_GRAIN;
        grain->grainID = 1;
        grain->sign[0] = -1.0;
        grain->GB[0] = home->GBKeys[0];
    
	    CopyMatrix33(TwinOrientation, grain->orientation);
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
	numGrains = home->nGrains;

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

	    for (j = 0; j < grain->nBounds; j++){
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
            bLab[0] = node->burgX[segID];
            bLab[1] = node->burgY[segID];
            bLab[2] = node->burgZ[segID];
            Matrix33Vector3Multiply(home->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[0] = node->burgX[segID];
            bCryst[1] = node->burgY[segID];
            bCryst[2] = node->burgZ[segID];
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
	int		i, j, numGrains, nBounds;
	int		grainFlag;

	numGrains = home->nGrains;
	
	if(home->param->polycrystal == 0)return(0);

	ID = 0;
	
//	FoldBox(home->param, &x, &y, &z);
	for (i = 1; i < numGrains; i++){
	    if ((grain = home->grainKeys[i])==(Grain_t *)NULL){
	        continue;
	    }

	    grainFlag = 0;
	    nBounds = grain->nBounds;

	    for (j = 0; j < nBounds; j++){
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

int AdjustNodeDirection(Home_t *home, Node_t *node, real8 *vec){

    int 	NM;
    real8 	GC[3][3];

	NM = 0;

	GetGlideConstraintList(home, node, &NM, GC);
	ApplyConstraintsToVelocity(NM, GC, vec);

	return(NM);
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
		int			i, nBounds;
		
		grain = home->grainKeys[(int)inGrainPos];
		
		nBounds = grain->nBounds;
		for(i=0; i<nBounds; i++){
			vec[0] = node->x - grain->x0[i];
			vec[1] = node->y - grain->y0[i];
			vec[2] = node->z - grain->z0[i];
			plane[0] = grain->nx[i];
			plane[1] = grain->ny[i];
			plane[2] = grain->nz[i];
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
	real8	tempDistance, vec[3], vec1[3];
	int	existTB = 0;
	real8 	minDistance = 1e14;
	
	if(home->param->nGrain > 1)existTB = 1;
	if(Normal(GBPlane->normal) < 1.0e-3){	// the GB was not defined, so should define the nearest GB. 
	    int			i,j;
	    Grain_t 	*grain;
	    Plane_t 	plane;
	    	
	    for(i=existTB; i<home->nGrains; i++){
	        grain = home->grainKeys[i];
	        for(j=0; j<grain->nBounds; j++){
	            plane.normal[0] = grain->nx[j];
	            plane.normal[1] = grain->ny[j];
	            plane.normal[2] = grain->nz[j];
	            plane.point[0] = grain->x0[j];
	            plane.point[1] = grain->y0[j];
	            plane.point[2] = grain->z0[j];
	            	
	            PointPlaneIntersection(pos, plane.normal, plane.point, vec, &tempDistance);

	            if(fabs(tempDistance) < fabs(minDistance)){
	                minDistance = tempDistance;
	                VECTOR_COPY(GBPlane->normal,plane.normal);
	                VECTOR_COPY(GBPlane->point,plane.point);
	            }
	        }
	    }
		
	}else{
        PointPlaneIntersection(pos, GBPlane->normal, GBPlane->point, vec, &minDistance);	
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
    real8       vec1[3];
	int		    numGlideConstraints;
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
	
	pos[0] = node->x; pos[1] = node->y; pos[2] = node->z;
	FoldBox(home->param, pos, pos+1, pos+2);
	currInGrain = InWhichGrain(home, pos[0], pos[1], pos[2]);

	minDistance = ClosedGB(home, pos, &GBPlane);
    
	if(oldInGrain - floor(oldInGrain) == 0.5){
        if(fabs(minDistance) > 2.0*eps){
            printf("warning: the GB node (%d,%d) diviate far from GB, minDistance is %f\n"
                   ,node->myTag.domainID,node->myTag.index,
                   minDistance);
        }
        newPos[0] = pos[0] - GBPlane.normal[0] * minDistance;
        newPos[1] = pos[1] - GBPlane.normal[1] * minDistance;
        newPos[2] = pos[2] - GBPlane.normal[2] * minDistance;
	
	    if(!(NodeOnGB(home, node) & TB_NODE)){
	        newPos[0] =  node->oldx;
	        newPos[1] =  node->oldy;
	        newPos[2] =  node->oldz;
	    }

	}else{
        if(fabs(minDistance) > 1.0e1*eps){
            FormatVector(GBPlane.normal,"GB noromal");
            FormatVector(GBPlane.point,"GB noromal");
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
		
        GBPlane.point[0] -= GBPlane.normal[0] * eps * sb;
        GBPlane.point[1] -= GBPlane.normal[1] * eps * sb;
        GBPlane.point[2] -= GBPlane.normal[2] * eps * sb;
		
        VECTOR_COPY(vec, GBPlane.normal);
	    numGlideConstraints = AdjustNodeDirection(home, node, vec);
        if(Normal(vec) < 1.0E-4){
	        newPos[0] =  node->oldx;
	        newPos[1] =  node->oldy;
	        newPos[2] =  node->oldz;
            PrintNode(node);
            printf("PullToGB: node constraints is %d, and vec is zero\n", 
                   numGlideConstraints);
        }else{
            NormalizeVec(vec);
            switch(LinePlaneIntersection(vec, pos, GBPlane.normal, 
               GBPlane.point, vec1, newPos)){
                case POINT_INTERSECTION:
                    break;
                default:
                    printf("PullToGB: warning at Line %d", __LINE__);
                    break;
            }

            if(oldInGrain != InWhichGrain(home, newPos[0],newPos[1],newPos[2]) ){
                GBPlane.point[0] += (2*GBPlane.normal[0] * eps * sb);
                GBPlane.point[1] += (2*GBPlane.normal[1] * eps * sb);
                GBPlane.point[2] += (2*GBPlane.normal[2] * eps * sb);
                switch(LinePlaneIntersection(vec, pos, GBPlane.normal, 
                       GBPlane.point, vec1, newPos)){
                    case POINT_INTERSECTION:
                        break;
                    default:
                        printf("PullToGB: warning at Line %d", __LINE__);
                        break;
                }
            }
        }

        if(oldInGrain != InWhichGrain(home, newPos[0],newPos[1],newPos[2]) ){
            printf("PullToGB: at line %d\n", __LINE__);    
        }
        
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
	
	inGrainPos = GetInGrainFromPos(home,newPos[0],newPos[1],newPos[2]);

	if(inGrainPos != oldInGrain){
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

int ClosestPointOfPointLine(real8 *p, Line_t *line, real8 *outP){
	if(NodeOnLine(p, line)){
		VECTOR_COPY(outP, p);
		return(0);
	}
	
	real8 t,dotPLD, dotLPLD;
	
	dotPLD = DotProduct(p, line->direction);
	dotLPLD = DotProduct(line->point, line->direction);
	
	t = dotPLD - dotLPLD;
	
	outP[0] = line->point[0] + line->direction[0] * t;
	outP[1] = line->point[1] + line->direction[1] * t;
	outP[2] = line->point[2] + line->direction[2] * t;
	
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
	
	p1x = plane1->point[0]; p1y = plane1->point[1]; p1z = plane1->point[2]; 
	p2x = plane2->point[0]; p2y = plane2->point[1]; p2z = plane2->point[2]; 
	n1x = plane1->normal[0]; n1y = plane1->normal[1]; n1z = plane1->normal[2]; 
	n2x = plane2->normal[0]; n2y = plane2->normal[1]; n2z = plane2->normal[2]; 
	dx = line->direction[0]; dy = line->direction[1]; dz = line->direction[2];
	
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
	
	line->point[0] = -(dz*n1x*n2y*p1x - dx*n1z*n2y*p1x - dy*n1x*n2z*p1x + dx*n1y*n2z*p1x + 
					   dz*n1y*n2y*p1y - dy*n1z*n2y*p1y + dz*n1y*n2z*p1z - 
					   dy*n1z*n2z*p1z - dz*n1y*n2x*p2x + dy*n1z*n2x*p2x - 
					   dz*n1y*n2y*p2y + dy*n1z*n2y*p2y - dz*n1y*n2z*p2z + dy*n1z*n2z*p2z)/denominator;
	
	line->point[1] = -((-dz)*n1x*n2x*p1x + dx*n1z*n2x*p1x - dz*n1y*n2x*p1y + 
					   dy*n1z*n2x*p1y - dy*n1x*n2z*p1y + dx*n1y*n2z*p1y - dz*n1x*n2z*p1z + 
					   dx*n1z*n2z*p1z + dz*n1x*n2x*p2x - dx*n1z*n2x*p2x + 
					   dz*n1x*n2y*p2y - dx*n1z*n2y*p2y + dz*n1x*n2z*p2z - dx*n1z*n2z*p2z)/denominator;
	
	line->point[2] = -(dy*n1x*n2x*p1x - dx*n1y*n2x*p1x + dy*n1x*n2y*p1y - dx*n1y*n2y*p1y - 
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
	        nDir[0] = node->nx[i];
	        nDir[1] = node->ny[i];
	        nDir[2] = node->nz[i];
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


void PrintTag(Node_t *node){

        if(node == (Node_t *)NULL){
            printf("(NULL)");
        }else{
            printf("[%d %d %d %3.1f %d]-{%e,%e,%e} ",node->myTag.domainID,
                   node->myTag.index,node->numNbrs,node->inGrain,
	               node->constraint, node->x,node->y,node->z);
	}
		

}

void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]){
    NormalizeVec(a);
    real8 u,v,w;
    
    u = a[0];
    v = a[1];
    w = a[2];
    R[0][0] = (1.0 - cos(theta))*u*u + cos(theta);
    R[0][1] = u*v*(1.0 - cos(theta)) - w*sin(theta);
    R[0][2] = u*w*(1.0 - cos(theta)) + v*sin(theta);
    R[1][0] = u*v*(1.0 - cos(theta)) + w*sin(theta);
    R[1][1] = (1.0 - cos(theta))*v*v + cos(theta);
    R[1][2] = v*w*(1.0 - cos(theta)) - u*sin(theta);
    R[2][0] = u*w*(1.0 - cos(theta)) - v*sin(theta);
    R[2][1] = v*w*(1.0 - cos(theta)) + u*sin(theta);
    R[2][2] = (1.0 - cos(theta))*w*w + cos(theta);
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
