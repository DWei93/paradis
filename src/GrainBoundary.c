/****************************************************************************
* 
*	Module:		GrainBoundary.c
*
*	Description:	This module contains subroutines needed in 
*                  	grain/twin boundary model.
*
*	08/07/2016:	created by HD Fan, WLY 
*
*****************************************************************************/
#include <stdio.h>
#include <math.h>
#include "Home.h" 
#include "Grain.h"
#include "Util.h" 

#define ENERGY_RULER 1

Node_t *InsertNode(Home_t *home, Node_t *node1, Node_t *node2, real8 *p,int oldOrCurr, int globalOp){
	
	int	i,j,k,splitStatus,armCount = 1;
	real8	nodePos[3],newPos[3],nodeVel[3],newVel[3],vec[3];
	Node_t	*splitNode1,*splitNode2,*newNode;
	Param_t	*param;
	Tag_t	*tag1,*tag2,*tag3;	
	param = home->param;
	
	if(node1 == (Node_t *)NULL || node2 == (Node_t *)NULL){
		printf("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		return((Node_t *)NULL);
	}


/*          
 *          In case that node1 and node2 both have penetrated GB, the newPos 
 *          is in the other grain. So we use the old coordinates, this ensures
 *          that newPos is in the grain before penetration.
 */
	if(oldOrCurr == 0){
		
		if(InWhichGrain(home, node2->oldx,node2->oldy,node2->oldz) != node2->inGrain){
			vec[X] = node2->x - node1->oldx;
			vec[Y] = node2->y - node1->oldy;
			vec[Z] = node2->z - node1->oldz;
		}else{
			vec[X] = node2->oldx - node1->oldx;
			vec[Y] = node2->oldy - node1->oldy;
			vec[Z] = node2->oldz - node1->oldz;
		}
	}else if(oldOrCurr == 1){
		vec[X] = node2->x - node1->x;
		vec[Y] = node2->y - node1->y;
		vec[Z] = node2->z - node1->z;
	}

	ZImage(param, &vec[X], &vec[Y], &vec[Z]);

	if(oldOrCurr == 0){ 
		newPos[X] = node1->oldx + vec[X] * 0.5;
		newPos[Y] = node1->oldy + vec[Y] * 0.5;
		newPos[Z] = node1->oldz + vec[Z] * 0.5;

		if(InWhichGrain(home, newPos[X],newPos[Y],newPos[Z]) != node1->inGrain){
			printf("Warning: the position of insert node is error\n");
			FormatVector(newPos,"newPos");
		}
		
	}else if(oldOrCurr == 1){
		newPos[X] = node1->x + vec[X] * 0.5;
		newPos[Y] = node1->y + vec[Y] * 0.5;
		newPos[Z] = node1->z + vec[Z] * 0.5;
		
	}else{
		
		VECTOR_COPY(newPos, p);
	}

	newVel[X] = (node1->vX + node2->vX) * 0.5;
	newVel[Y] = (node1->vY + node2->vY) * 0.5;
	newVel[Z] = (node1->vZ + node2->vZ) * 0.5;
	
	/*
	 *          Get an estimate of forces that would change due to refining
	 *          the segment.
	 */
	/*EstRefinementForces(home, node1, node2, newPos, vec1,
		 f0Seg1, f1Seg1, f0Seg2, f1Seg2);*/

	FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

	/*
	 *          This should be a global operation distributed out to remote
	 *          domains only if the segment spans domains.
	 */
	
	
	nodePos[X] = node1->x;
	nodePos[Y] = node1->y;
	nodePos[Z] = node1->z;
	
	nodeVel[X] = node1->vX;
	nodeVel[Y] = node1->vY;
	nodeVel[Z] = node1->vZ;
	
	tag1 = &node1->myTag; 
	tag2 = &node2->myTag; 
/*
 *          Need a critical section now because we'll be inserting a
 *          node and changing connections between nodes...
 */

	int armIndex;

/*
 *      When refining a segment, we always move exactly one arm
 *      from the original node to the node being created...
 */
	armIndex = GetArmID(home, node1, node2);

	if (armIndex >= 0) {
			
	     splitStatus = SplitNode(home, OPCLASS_SEPARATION,
				     node1, nodePos, newPos,
				     nodeVel, newVel,
				     armCount, &armIndex,
				     globalOp, &splitNode1,
				     &splitNode2, 0);

	     if(splitStatus != SPLIT_SUCCESS)return((Node_t *)NULL);
	}

	
	newNode =(splitNode2 == node1)? splitNode1: splitNode2; 

	
	if(globalOp){
		MarkNodeForceObsolete(home, node1);
		MarkNodeForceObsolete(home, node2);
		MarkNodeForceObsolete(home, newNode);
	}

	

	return(newNode);
}



/****************************************************************************  
*
*	calculate the glide force towards GB of the node	
*
*****************************************************************************/
real8 FGlidePerLength(Home_t *home,Node_t *node,Plane_t *GBPlane){

	int	i;
	real8	normal[3],disLength;
	real8	intersection[3];
	real8	glideDir[3];
	real8	fGlide;
	Node_t  *nbr;
	
	if(node == (Node_t *)NULL){
		printf("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		return(0);
	}

	normal[0] = node->nx[0];	
	normal[1] = node->ny[0];	
	normal[2] = node->nz[0];	

	disLength = 0;

	for(i=0;i<node->numNbrs;i++){

	   nbr = GetNeighborNode(home,node,i);

	   if(nbr==(Node_t *)NULL) continue;
	   
	   disLength += NodeDistance(home,node,nbr)/2;
	}

	cross(normal,GBPlane->normal,intersection);
	cross(normal,intersection,glideDir);

	NormalizeVec(glideDir);

	fGlide = node->fX*glideDir[0]+
		 node->fY*glideDir[1]+
		 node->fZ*glideDir[2];

	return(fabs(fGlide)/disLength);
}

/****************************************************************************  
*
*	get the GB plane begween oldGrain and currGrain
*
*****************************************************************************/
void GetGBPlane(Home_t *home,real8 oldGrain, real8 currGrain, 
                Node_t *node, Plane_t *GBPlane){
	
	real8	matrix,twin;
	real8	x,y,z,nx,ny,nz,x0,y0,z0,yy,minDis = 1e14;
	int		i;
	Grain_t	*grain;
	
	if(node == (Node_t *)NULL){
            Fatal("Task %d: %s, line %d -- Node = NULL",
                  home->myDomain, __FILE__, __LINE__);
	}
	
	if(oldGrain == currGrain && oldGrain - floor(oldGrain) == 0){
            PrintNode(node);
            Fatal("GetGBPlane: the node was't changed into inGrain.");
	}
	
	if(oldGrain > 0 && currGrain > 0){
            PrintNode(node);
            Fatal("GetGBPlane: the node penetrated though 2 layer at least.");
	}
	
	if(oldGrain == 0 || oldGrain == 0.5){
            matrix = oldGrain;
            twin   = currGrain;
	}else{
            matrix = currGrain;
            twin   = oldGrain;
	}
	
	x = node->x; y = node->y; z = node->z;
	grain = home->grainKeys[(int)twin];
	for(i=0; i<grain->numBoundaries; i++){
		nx = grain->nx[i]; ny = grain->ny[i]; nz = grain->nz[i];
		x0 = grain->x0[i]; y0 = grain->y0[i]; z0 = grain->z0[i];
		yy = nx*(x - x0) + ny*(y - y0) + nz*(z - z0);
		if(fabs(yy) < minDis){
			GBPlane->normal[X] = nx; 	GBPlane->normal[Y] = ny; 	GBPlane->normal[Z] = nz; 
			GBPlane->point[X] = x0; 	GBPlane->point[Y] = y0; 	GBPlane->point[Z] = z0; 
			minDis = fabs(yy);
		}
	}
	return;
}


void GetIntersection(Plane_t *plane1, Plane_t *plane2, Line_t *inLine, Line_t *outLine){	
	
	cross(plane1->normal, plane2->normal,outLine->direction);
	NormalizeVec(outLine->direction);
	
	if(LineOnPlane(inLine, plane1)){
		if(!IntersectionOfLinePlane(inLine, plane2, outLine->point)){
			
			FormatVector(inLine->direction,"inLine.direction");
			FormatVector(plane2->normal,"plane2.normal");
			Fatal("GetIntersection: can't find the intersection");
		}
	}else{
		Fatal("the inLine is't on plane1");
	}
}

int InterPointOfLines(real8 *a_ori,real8 *b, real8 *c,real8 *d,real8 *r){

	real8	a[3];

	/*	
	 *	Sometimes, we call this function to update the coordinate of node, 
	 *	i.e. a and r all point to the node coordinate.
	 *	so we use a temporary array to store the original coordinate 	
	 */

	VECTOR_COPY(a,a_ori); 

	if(fabs(b[1]*d[0] - b[0]*d[1]) > 1e-3){

		r[0] = -((a[1]*b[0]*d[0] - a[0]*b[1]*d[0] - 
				  b[0]*c[1]*d[0] + b[0]*c[0]*d[1])/
				 (b[1]*d[0] - b[0]*d[1]));
		r[1] = -((b[1]*c[1]*d[0] - a[1]*b[0]*d[1] + 
				  a[0]*b[1]*d[1] - b[1]*c[0]*d[1])/
				 (-(b[1]*d[0]) + b[0]*d[1]));
		r[2] = -((a[2]*b[1]*d[0] - a[1]*b[2]*d[0] + 
				  b[2]*c[1]*d[0] - a[2]*b[0]*d[1] + 
				  a[0]*b[2]*d[1] - b[2]*c[0]*d[1])/
				 (-(b[1]*d[0]) + b[0]*d[1]));
		return(1);
	}

	if(fabs(b[2]*d[0] - b[0]*d[2]) > 1e-3){

		r[0] = -((a[2]*b[0]*d[0] - a[0]*b[2]*d[0] - 
          		  b[0]*c[2]*d[0] + b[0]*c[0]*d[2])/
				 (b[2]*d[0] - b[0]*d[2]));
		r[1] = -((-(a[2]*b[1]*d[0]) + a[1]*b[2]*d[0] + 
				  b[1]*c[2]*d[0] - a[1]*b[0]*d[2] + 
				  a[0]*b[1]*d[2] - b[1]*c[0]*d[2])/
				 (-(b[2]*d[0]) + b[0]*d[2]));
		r[2] = -((b[2]*c[2]*d[0] - a[2]*b[0]*d[2] + 
				  a[0]*b[2]*d[2] - b[2]*c[0]*d[2])/
				 (-(b[2]*d[0]) + b[0]*d[2]));
		return(1);
	}

	if(fabs(-(b[2]*d[1]) + b[1]*d[2]) > 1e-3){

		r[0] =  -((-(a[2]*b[0]*d[1]) + a[0]*b[2]*d[1] + 
				   b[0]*c[2]*d[1] + a[1]*b[0]*d[2] - 
				   a[0]*b[1]*d[2] - b[0]*c[1]*d[2])/
				  (-(b[2]*d[1]) + b[1]*d[2]));
		r[1] = -((a[2]*b[1]*d[1] - a[1]*b[2]*d[1] - 
				  b[1]*c[2]*d[1] + b[1]*c[1]*d[2])/
				 (b[2]*d[1] - b[1]*d[2]));
		r[2] = -((b[2]*c[2]*d[1] - a[2]*b[1]*d[2] + 
         		  a[1]*b[2]*d[2] - b[2]*c[1]*d[2])/
				 (-(b[2]*d[1]) + b[1]*d[2]));
		return(1);
	}

	printf("InterPointOfLines: cannot find intersection point. %f %f %f %f %f %f\n",
		   b[0],b[1],b[2],d[0],d[1],d[2]);
	
	r[0] = 0.0;
	r[1] = 0.0;
	r[2] = 0.0;

	return(0);
}

int MoveNodeToIntersection(Home_t *home, Node_t *node, real8 direction[3],
						   Line_t *interLine, real8 *nodePos){
	
	if(node == (Node_t *)NULL){
		Fatal("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
		
	}

	nodePos[0] = node->oldx;	
	nodePos[1] = node->oldy;	
	nodePos[2] = node->oldz;	
	
	if(InterPointOfLines(nodePos,direction,interLine->point,
						 interLine->direction,nodePos) == 0){
		PrintNode(node);
		Fatal("MoveNodeToIntersection: the node is not on the line");
	}
	
	FoldBox(home->param,nodePos,nodePos+1,nodePos+2);	
	
	if(!NodeOnLine(nodePos, interLine)){
		
		PullToLine(nodePos, interLine);
	}
	
	FoldBox(home->param,nodePos,nodePos+1,nodePos+2);

	return(1);
}

real8 SegementEnergy(Param_t *param, real8 cosTheta, real8 bMag){
	static real8 pi = 3.1415926535897932385;
	
	return(param->shearModulus * bMag * bMag * log(param->outRadius) *
		   (1.0 - param->pois * cosTheta * cosTheta)/
		   (4.0 * pi * (1.0 - param->pois)) );
		  
}

int EnergyRule(Home_t *home, Node_t *node, Plane_t *GBPlane, PeneRes_t *peneRes){
	
	int		i, isScrew = 0;
	real8	inBurg[3], inEnergy=0.0, cosTheta, bMag;
	real8	segEnergy, outEnergy=0.0;
	real8	vec[3], normal[3], temVec[3],burgMag;
	Param_t *param;
	Node_t 	*nbr;
	
	param = home->param;
	
#if 0
	for(i=0; i<node->numNbrs; i++){
		nbr = GetNeighborNode(home, node, i);
		if(nbr == (Node_t * )NULL){
			printf("WARNING: Neighbor not found at %s line %d\n",
				   __FILE__, __LINE__);
			continue;
		}
		inBurg[X] = node->burgX[i];
		inBurg[Y] = node->burgY[i];
		inBurg[Z] = node->burgZ[i];
		bMag = Normal(inBurg); 
		
		NormalizeVec(inBurg);
		vec[X] = nbr->x - node->x;
		vec[Y] = nbr->y - node->y;
		vec[Z] = nbr->z - node->z;
		ZImage(param, vec, vec+1, vec+2);
		NormalizeVec(vec);
		
		cosTheta = DotProduct(vec, inBurg);
		segEnergy = SegementEnergy(param, cosTheta, bMag);
		inEnergy += (segEnergy/2.0);
	}

	normal[X] = node->nx[1];
	normal[Y] = node->ny[1];
	normal[Z] = node->nz[1];
	temVec[X] = node->burgX[1];
	temVec[Y] = node->burgY[1];
	temVec[Z] = node->burgZ[1];
	burgMag = Normal(temVec);
	
	cross(normal, GBPlane->normal, vec);
	NormalizeVec(vec);
	
	for(i=0; i<peneRes->nSets; i++){
		outEnergy = 0.0;
		if(peneRes->effective[i] == 0)continue;
		cosTheta = DotProduct(peneRes->outBurg[i],vec);
		bMag = Normal(peneRes->outBurg[i]);
		segEnergy = SegementEnergy(param, cosTheta, bMag);
		outEnergy += segEnergy;
		
		cosTheta = DotProduct(peneRes->GBBurg[i],vec);
		bMag = Normal(peneRes->GBBurg[i]);
		segEnergy = SegementEnergy(param, cosTheta, bMag);
		outEnergy += segEnergy;
		
		outEnergy += burgMag * param->TBEnergy/param->burgMag;
		
		if(FGlidePerLength(home,node,GBPlane) * burgMag < outEnergy - inEnergy){
			peneRes->effective[i] = 0;
		
		}
	}
#endif 
	
	temVec[X] = node->burgX[0];
	temVec[Y] = node->burgY[0];
	temVec[Z] = node->burgZ[0];
	
	if(fabs(DotProduct(temVec, GBPlane->normal)) < 1.0e-2)isScrew = 1;
	

	if(isScrew){
	    if(fabs(FGlidePerLength(home,node,GBPlane))  < home->param->GBStrength ){
	        return(1);
	    }else{
	        return(1);
	    }
	}else{
	    if(fabs(FGlidePerLength(home,node,GBPlane))  < home->param->GBStrength ){
	        return(0);
	    }else{
	        return(1);
	    }
	}

	return(1);
}

int ScrewEnergyRule(Home_t *home, Node_t *node, real8 oldInGrain, 
					real8 corInGrain, Plane_t *GBPlane){

	return(1);
}

void InitPeneRes(PeneRes_t *peneRes){
	int	i;
	real8 zeroVec[3];
	
	zeroVec[X] = 0.0;
	zeroVec[Y] = 0.0;
	zeroVec[Z] = 0.0;
	
	for(i=0; i<peneRes->nSets; i++){
		peneRes->effective[i] = 0;
		VECTOR_COPY(peneRes->outNormal[i],zeroVec);
		VECTOR_COPY(peneRes->outBurg[i],zeroVec);
		VECTOR_COPY(peneRes->GBBurg[i],zeroVec);
		peneRes->powerDis[i] = -1.0e50;
		
	}
	return;
}

int GetPeneResults(Home_t *home, Node_t *node, Node_t *nbr, PeneRes_t *peneRes,
		   Grain_t *oldGrain, Grain_t *currGrain, Plane_t *GBPlane){
	
	real8		pos[3], burg[3], normal[3], crystalNormal[3], burg1[3], burg2[3];
	real8		crystalBurg[3], dot, vec[3], outNormal[3], zeroVec[3];
        real8           sb, TBBurgCryst[3];
	int 		isScrew = 0, armID, i, normalID;
	
	InitPeneRes(peneRes);

	pos[X] = node->x;
	pos[Y] = node->y;
	pos[Z] = node->z;
	
	zeroVec[X] = 0.0;
	zeroVec[Y] = 0.0;
	zeroVec[Z] = 0.0;
	
	armID = GetArmID(home, node, nbr);

	burg[X] = node->burgX[armID];
	burg[Y] = node->burgY[armID];
	burg[Z] = node->burgZ[armID];
	
	normal[X] = node->nx[armID];
	normal[Y] = node->ny[armID];
	normal[Z] = node->nz[armID];
	
	if(fabs(DotProduct(normal, GBPlane->normal)) > 0.95){
            FormatVector(normal,"GetPeneRes: arm noram parallel to GB");
            return(0);
        }

	dot = DotProduct(burg, GBPlane->normal);
 	if(fabs(dot) <1.0e-2)isScrew = 1;
        	

#if 0
	Matrix33Vector3Multiply(oldGrain->orientationInv, normal, crystalNormal);
	Matrix33Vector3Multiply(currGrain->orientation, crystalNormal, outNormal);
#endif
	
	Matrix33Vector3Multiply(home->param->rotateNormal, normal, outNormal);
        outNormal[X] = -outNormal[X];
        outNormal[Y] = -outNormal[Y];
        outNormal[Z] = -outNormal[Z];

	Matrix33Vector3Multiply(oldGrain->orientationInv, burg, crystalBurg);

	if(isScrew){
		peneRes->nSets = 2;
		
		for(i=0; i<peneRes->nSets; i++){
			peneRes->effective[i] = 1;
		}
	
		VECTOR_COPY(peneRes->outNormal[0],outNormal);
		VECTOR_COPY(peneRes->outBurg[0],burg);
		VECTOR_COPY(peneRes->GBBurg[0],zeroVec);

		peneRes->effective[1] = 1;

		VECTOR_COPY(peneRes->outNormal[1],GBPlane->normal);
		switch(BurgID_FCC(crystalBurg)){
		    case 1 :
		    burg1[X] = -1.0;
		    burg1[Y] = 2.0;
		    burg1[Z] = -1.0;
		        break;
		    case 3 :
		    burg1[X] = 2.0;
		    burg1[Y] = -1.0;
		    burg1[Z] = -1.0;
		        break;
		    case 5 :
		    burg1[X] = -1.0;
		    burg1[Y] = 2.0;
		    burg1[Z] = -1.0;
		        break;
		    default :
			Fatal("GetPeneRes: Wrong screw ID");
		}

		NormalizeVec(burg1);
		burg1[X] /= sqrt(3.0);	
		burg1[Y] /= sqrt(3.0);	
		burg1[Z] /= sqrt(3.0);	
		
	        Matrix33Vector3Multiply(oldGrain->orientation, burg1, burg2);

		if(DotProduct(burg2, burg) < 0.0){
		    burg2[X] = -burg2[X];
		    burg2[Y] = -burg2[Y];
		    burg2[Z] = -burg2[Z];
		}

		burg1[X] = burg[X] - burg2[X];
		burg1[Y] = burg[Y] - burg2[Y];
		burg1[Z] = burg[Z] - burg2[Z];

		VECTOR_COPY(peneRes->outBurg[1],burg1);
		VECTOR_COPY(peneRes->GBBurg[1],burg2);
		return(1);
	}


        if(BurgID_FCC(crystalBurg) % 2 != 0){
            FormatVector(crystalBurg,"crystalBurg");
            FormatVector(burg,"burg");
            FormatVector(normal,"normal");
            FormatVector(outNormal,"outNormal");
            FormatVector(GBPlane->normal,"GBPlane->normal");
            printf("the burg is not a non-screw burg, dot %.15e\n",dot);
	    return(0);
        }
        
        if(crystalBurg[X] < 1.0e-6 &&
           crystalBurg[Y] < 1.0e-6 &&
           crystalBurg[Z] < 1.0e-6){
           sb = -1.0;
       }else sb = 1.0;

	Matrix33Vector3Multiply(oldGrain->orientationInv, normal, crystalNormal);
        normalID = GlidePlaneID_FCC(crystalNormal);
        GetTBBurg_FCC(normalID, TBBurgCryst);
        TBBurgCryst[X] *= sb;
        TBBurgCryst[Y] *= sb;
        TBBurgCryst[Z] *= sb;
	Matrix33Vector3Multiply(oldGrain->orientation, TBBurgCryst, burg2);

	burg1[X] = burg[X] - burg2[X];
	burg1[Y] = burg[Y] - burg2[Y];
	burg1[Z] = burg[Z] - burg2[Z];
	
        real8 burg1Mag;
        burg1Mag = Normal(burg1);
        if(fabs(burg1Mag - 1.0) > 1.0e-3){
#if 0
            FormatVector(burg,"burg");
            FormatVector(burg1,"burg1");
            FormatVector(burg2,"burg2");
            Fatal("GetPeneResults: burg1Mag is %f",burg1Mag);
#else
	    return(0);
#endif
	}

	peneRes->nSets = 1;
	VECTOR_COPY(peneRes->outNormal[0],outNormal);
	VECTOR_COPY(peneRes->outBurg[0],burg1);
	VECTOR_COPY(peneRes->GBBurg[0],burg2);
	peneRes->effective[0] = 1;

	return(1);
}


void PullToLine(real8 *p, Line_t *line){
	real8 t;
	
	t = line->direction[X]*(p[X] - line->point[X])+ 
		line->direction[Y]*(p[Y] - line->point[Y])+ 
		line->direction[Z]*(p[Z] - line->point[Z]);
	
	p[X] = line->point[X] + t*line->direction[X];
	p[Y] = line->point[Y] + t*line->direction[Y];
	p[Z] = line->point[Z] + t*line->direction[Z];
}

int PenetrationFirst(Home_t *home, Node_t *node, Grain_t *oldGrain,
					 Grain_t *newGrain, Plane_t *GBPlane, PeneRes_t *peneRes){

	real8		nDir1[3], nodeVel[3], currPos[3], oldPos[3], burg[3], bMag;
	real8 		burgStore[3], normalStore[3],eps, sb, newVelocity, CurrNbrInGrain;
	real8		velNormal, glideDir[3], burg1[3],burg2[3], newVel[3], oldVel[3],crystalNormal[3];
	real8		newPos[3], vec[3], vecNormal, crossDistance, TBInGrain,currInGrain, oldInGrain;
	real8		peneVel[3],node1Pos[3],node2Pos[3],node3Pos[3],centorNode12[3];
	real8		maxPowerDis = - 1.0e50;
	int			armID, switchNode=0,i, set, peneThough, maxPowerSet=-1;
	int			iScrew = 0, globalOp, mobError, forward, myDomian;
	int		flags;
	Node_t 		*nbrNode1,*nbrNode2,*newNode,*globalNode, *nbr;
	Node_t		*node1,*node2,*node3, *node1Bak,*bkupNodeList;
	Node_t		*node4,*node5;
	Param_t		*param;
	Plane_t		oldGlidePlane, outPlane, plane;							//Basic plane
	Line_t		pathLine, interLine, lineNodeOldNbr2, lineNbr1Node, estimatedPath;
	
	myDomian = home->myDomain;
	
	param = home->param;
	eps = param->rann;
	
	nDir1[0] = node->nx[0];
	nDir1[1] = node->ny[0];
	nDir1[2] = node->nz[0];
	
	nodeVel[0] = node->vX;
	nodeVel[1] = node->vY;
	nodeVel[2] = node->vZ;
	
	currPos[0] = node->x;
	currPos[1] = node->y;
	currPos[2] = node->z;
	
	oldPos[0] = node->oldx;
	oldPos[1] = node->oldy;
	oldPos[2] = node->oldz;
	
	burg[0] = node->burgX[0];
	burg[1] = node->burgY[0];
	burg[2] = node->burgZ[0];
	
	burgStore[X] = node->burgX[1];
	burgStore[Y] = node->burgY[1];
	burgStore[Z] = node->burgZ[1];
	
	normalStore[X] = node->nx[1];
	normalStore[Y] = node->ny[1];
	normalStore[Z] = node->nz[1];
	
	flags = node->flags;
	
	bMag = Normal(burg);
	
	if(fabs(DotProduct(burg, GBPlane->normal))<0.001){
		iScrew=1;
	}else{
		iScrew = 0;
	}
	
	if(fabs(1.0-fabs(DotProduct(nDir1,GBPlane->normal))) < 1.0e-2)return(0);	// the glide plane parallel to th GB.
	
	cross(nDir1,GBPlane->normal,interLine.direction);
	NormalizeVec(interLine.direction);
	
	VECTOR_COPY(oldGlidePlane.normal, nDir1);
	VECTOR_COPY(oldGlidePlane.point, oldPos);
	VECTOR_COPY(pathLine.point, oldPos);
	VECTOR_COPY(plane.normal, nDir1);
	
	oldInGrain = node->inGrain;
	currInGrain = InWhichGrain(home, currPos[0], currPos[1], currPos[2]);
	
/*	
 *	find the dislocation glide direction	
 */
	if(Normal(nodeVel) < 1.0)return(0);
	
	velNormal = DotProduct(nodeVel, nDir1);

	glideDir[0] = nodeVel[0] - velNormal*nDir1[0];
	glideDir[1] = nodeVel[1] - velNormal*nDir1[1];
	glideDir[2] = nodeVel[2] - velNormal*nDir1[2];
	
	NormalizeVec(glideDir);	
	
	if(fabs(DotProduct(glideDir,interLine.direction))>0.995) {
		cross(nDir1,interLine.direction,glideDir);
		NormalizeVec(glideDir);
	}
	VECTOR_COPY(pathLine.direction, glideDir);
	
	GetIntersection(&oldGlidePlane, GBPlane, &pathLine, &interLine);

	
	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	
	if(nbrNode1 == (Node_t *)NULL || nbrNode2 == (Node_t *)NULL)return(0);
	
	bkupNodeList = (Node_t *)malloc(5 * sizeof(Node_t));
	
	FoldBox(home->param, &nbrNode1->x, &nbrNode1->y, &nbrNode1->z);
	FoldBox(home->param, &nbrNode2->x, &nbrNode2->y, &nbrNode2->z);

	if(nbrNode1==(Node_t *)NULL || nbrNode1==(Node_t *)NULL) return(0);
	
	BackupNode(home, nbrNode1, &bkupNodeList[0]);
	BackupNode(home, node, &bkupNodeList[1]);
	BackupNode(home, nbrNode2, &bkupNodeList[2]);
	
#if 0
	if(nbrNode1->inGrain - floor(nbrNode1->inGrain) == 0.5 ||
	   nbrNode2->inGrain - floor(nbrNode2->inGrain) == 0.5){
		PrintNode(nbrNode1);
		PrintNode(nbrNode2);
		Fatal("PenetrationFirst: neighbor is a GB node");
	}
#endif
	
	node1 = (Node_t *)NULL; 
	node2 = (Node_t *)NULL; 
	node3 = (Node_t *)NULL; 
	node4 = (Node_t *)NULL;
	
//	MarkNodeForceObsolete(home, nbrNode1);
//	MarkNodeForceObsolete(home, nbrNode2);
	
	
// node2 pos
	FoldBox(param,&nbrNode2->x,&nbrNode2->y,&nbrNode2->z);
	CurrNbrInGrain = InWhichGrain(home,nbrNode2->x,nbrNode2->y,nbrNode2->z);
	
	if( CurrNbrInGrain != nbrNode2->inGrain || CurrNbrInGrain == -1.0){
		
            vec[0] = nbrNode2->oldx - oldPos[0];
	    vec[1] = nbrNode2->oldy - oldPos[1];
	    vec[2] = nbrNode2->oldz - oldPos[2];

	}else {
		
            vec[0] = nbrNode2->x - oldPos[0];
	    vec[1] = nbrNode2->y - oldPos[1];
	    vec[2] = nbrNode2->z - oldPos[2];
	}
	
	ZImage(param, &vec[0], &vec[1], &vec[2]);
	
	newPos[X] = oldPos[X] + 0.5 * vec[X];
	newPos[Y] = oldPos[Y] + 0.5 * vec[Y];
	newPos[Z] = oldPos[Z] + 0.5 * vec[Z];
	
	FoldBox(param, newPos, newPos+1, newPos+2);
	VECTOR_COPY(plane.point,newPos);
	
	ClosetPointToGB(param, &plane, GBPlane, newPos);
	
	VECTOR_COPY(node2Pos,newPos);

//node1pos
	CurrNbrInGrain = InWhichGrain(home,nbrNode1->x,nbrNode1->y,nbrNode1->z);
	
	if( CurrNbrInGrain != nbrNode1->inGrain || CurrNbrInGrain == -1.0){
		
            vec[0] = nbrNode1->oldx - oldPos[0];
	    vec[1] = nbrNode1->oldy - oldPos[1];
	    vec[2] = nbrNode1->oldz - oldPos[2];
		
	}else {

            vec[0] = nbrNode1->x - oldPos[0];
	    vec[1] = nbrNode1->y - oldPos[1];
	    vec[2] = nbrNode1->z - oldPos[2];
	}
	
	ZImage(param, &vec[0], &vec[1], &vec[2]);
	
	newPos[X] = oldPos[X] + 0.5 * vec[X];
	newPos[Y] = oldPos[Y] + 0.5 * vec[Y];
	newPos[Z] = oldPos[Z] + 0.5 * vec[Z];
	
	FoldBox(param, newPos, newPos+1, newPos+2);
	VECTOR_COPY(plane.point,newPos);
	
	ClosetPointToGB(param, &plane, GBPlane, newPos);
	
	TBInGrain = GetInGrainFromPos(home, newPos[0], newPos[1], newPos[2]);
	if(TBInGrain - floor(TBInGrain) == 0){
            Fatal("the node is not pulled to the  GB");
	}
	
	VECTOR_COPY(node1Pos,newPos);

// centorNode12
	vec[X] = node1Pos[X] - node2Pos[X];
	vec[Y] = node1Pos[Y] - node2Pos[Y];
	vec[Z] = node1Pos[Z] - node2Pos[Z];
	ZImage(param, &vec[0], &vec[1], &vec[2]);
	centorNode12[X] = node2Pos[X] + 0.5 * vec[X];
	centorNode12[Y] = node2Pos[Y] + 0.5 * vec[Y];
	centorNode12[Z] = node2Pos[Z] + 0.5 * vec[Z];
	FoldBox(param, centorNode12, centorNode12+1, centorNode12+2);
	
//node3Pos

//	Matrix33Vector3Multiply(oldGrain->orientationInv, nDir1, crystalNormal);
//	Matrix33Vector3Multiply(newGrain->orientation, crystalNormal, outPlane.normal);

        
	Matrix33Vector3Multiply(param->rotateNormal, nDir1, outPlane.normal);
	VECTOR_COPY(outPlane.point,centorNode12);
	
	if(oldGrain->type == MATRIX_GRAIN && newGrain->type == TWIN_GRAIN){
            sb = -1.0;
	}else if(oldGrain->type == TWIN_GRAIN && newGrain->type == MATRIX_GRAIN){
            sb =  1.0;
	}else{
            Fatal("wrong grain type");
	}
	GBPlane->point[X] += sb * eps * GBPlane->normal[X];
	GBPlane->point[Y] += sb * eps * GBPlane->normal[Y];
	GBPlane->point[Z] += sb * eps * GBPlane->normal[Z];
	
	ClosetPointToGB(param, &outPlane, GBPlane, newPos);
	
	if(InWhichGrain(home,newPos[0],newPos[1],newPos[2]) != currInGrain){
#if 0
	    PrintNode(nbrNode1);
	    PrintNode(node);
	    PrintNode(nbrNode2);
            printf("Reposite new node faild\n");
	    VECTOR_COPY(newPos, centorNode12);
#else
	    return(0);
#endif
	}
	VECTOR_COPY(node3Pos,newPos);
	
	
/*	
 *	Insert a new node as a node on GB	
 */
	globalOp = 0;
	
	if((globalNode = InsertNode(home,node,nbrNode2,vec,0,globalOp)) == (Node_t *)NULL) {
	    return(0);
	}

	switchNode = 0;
	if(GetArmID(home, node, nbrNode1) == -1)switchNode=1;
	
	if(switchNode){
            node1 = globalNode;
            node2 = node;
		
	}else{
	    node1 = node;
	    node2 = globalNode;
	}
	
	RepositionNode(home, node2Pos, &node2->myTag, globalOp);
	RepositionNode(home, node1Pos, &node1->myTag, globalOp);
	BackupNode(home, node1, &bkupNodeList[3]);
	BackupNode(home, node2, &bkupNodeList[4]);
	
	for(set=0; set < peneRes->nSets; set++){
		
		if(peneRes->effective[set] == 0)continue;
		
		peneRes->powerDis[set] = 0.0;
		
		node1->inGrain = oldInGrain;
		node2->inGrain = oldInGrain;
		
		VECTOR_COPY(burg1, peneRes->outBurg[set]);
		VECTOR_COPY(burg2, peneRes->GBBurg[set]);
		
		peneThough = (fabs(DotProduct(peneRes->outNormal[set], GBPlane->normal)) < 0.995) ? 1 : 0;
		if(peneThough){
			node3 = InsertNode(home,node1,node2,node3Pos,2, globalOp);
		}else{
			node3 = InsertNode(home,node1,node2,centorNode12,2, globalOp);
		}
	  
		if(node3==(Node_t *)NULL){
			Fatal("PenetrationFirst: failed to creat node.\n");
		}
		
		if((GetArmID(home, node1, node2) != -1))Fatal("switchNode");
		

		ChangeArmBurg(home, node1, &node3->myTag, burg1[0], burg1[1], burg1[2],
					  peneRes->outNormal[set][X],peneRes->outNormal[set][Y],peneRes->outNormal[set][Z],
					  globalOp,DEL_SEG_NONE);
		ChangeArmBurg(home, node3, &node1->myTag, -burg1[0], -burg1[1], -burg1[2],
					 peneRes->outNormal[set][X],peneRes->outNormal[set][Y],peneRes->outNormal[set][Z],
					  globalOp,DEL_SEG_NONE);
		ChangeArmBurg(home, node3, &node2->myTag, burg1[0], burg1[1], burg1[2],
					 peneRes->outNormal[set][X],peneRes->outNormal[set][Y],peneRes->outNormal[set][Z],
					  globalOp,DEL_SEG_NONE);
		ChangeArmBurg(home, node2, &node3->myTag, -burg1[0], -burg1[1], -burg1[2],
					 peneRes->outNormal[set][X],peneRes->outNormal[set][Y],peneRes->outNormal[set][Z],
					  globalOp,DEL_SEG_NONE);

		node1->inGrain = TBInGrain;
		node2->inGrain = TBInGrain;
		if(peneThough){
			node3->inGrain = currInGrain;
		}else{
			node3->inGrain = TBInGrain;
		}

		node1->constraint = UNCONSTRAINED;
		node2->constraint = UNCONSTRAINED;

		
		if(!IsZeroVector(burg2)){
			InsertArm(home, node1, &node2->myTag, burg2[0], burg2[1],  burg2[2],
					  GBPlane->normal[0],GBPlane->normal[1],GBPlane->normal[2],globalOp);
			InsertArm(home, node2, &node1->myTag, -burg2[0], -burg2[1], -burg2[2],
					  GBPlane->normal[0],GBPlane->normal[1],GBPlane->normal[2],globalOp);
			
			if((node4 = InsertNode(home,node1,node2, vec, 1,globalOp)) == (Node_t *)NULL) {
				Fatal("PenetrationFirst: Insert node4 faild.");
			}
		}

		SetOneNodeForce(home,node3);
		mobError = EvaluateMobility(home,node3);
		
		peneRes->powerDis[set] += node3->fX * node3->vX + 
                                          node3->fY * node3->vY + 
                                          node3->fZ * node3->vZ;

		if(!IsZeroVector(burg2)){
			SetOneNodeForce(home,node4);
			mobError = EvaluateMobility(home,node4);
			
			peneRes->powerDis[set] += node4->fX * node4->vX + 
                                                  node4->fY * node4->vY + 
                                                  node4->fZ * node4->vZ;
		}

		newVel[X] = node3->vX; 
		newVel[Y] = node3->vY; 
		newVel[Z] = node3->vZ;
		
		forward = 1;
		if(peneThough){
			VECTOR_COPY(peneVel,newVel)
			newVelocity = Normal(newVel);
			NormalizeVec(newVel);
			
			vec[X] = newPos[X] - outPlane.point[X];
			vec[Y] = newPos[Y] - outPlane.point[Y];
			vec[Z] = newPos[Z] - outPlane.point[Z];
			ZImage(param, vec, vec+1, vec+2);
			NormalizeVec(vec);
			
			forward = (DotProduct(newVel,vec) > 0) ? 1 : 0;
			forward &= newVelocity > 1.0;
		}
		if(!forward)peneRes->effective[set] = 0;
		
		node3->numNbrs = 0;
		RemoveNode(home, node3, globalOp);
		
		if(!IsZeroVector(burg2)){
			node4->numNbrs = 0;
			RemoveNode(home, node4, globalOp);
		}
		
		RestoreNode(home, node1, &bkupNodeList[3]);
		RestoreNode(home, node2, &bkupNodeList[4]);
	
	}//end of for set

	globalNode->numNbrs = 0;
	RemoveNode(home, globalNode, globalOp);
	
	if(globalOp){
		AssignGrain(home, node, oldInGrain);
	}else{
		node->inGrain = oldInGrain;
	}

	RestoreNode(home, nbrNode1, &bkupNodeList[0]);
	RestoreNode(home, node, &bkupNodeList[1]);
	RestoreNode(home, nbrNode2, &bkupNodeList[2]);

	for (i = 0; i < 5; i++) {
	    FreeNodeArrays(&bkupNodeList[i]);
	}
	
	free(bkupNodeList);
	bkupNodeList = (Node_t *)NULL;
	
	node->x = currPos[0];
	node->y = currPos[1];
	node->z = currPos[2];
	node->oldx = oldPos[0];
	node->oldy = oldPos[1];
	node->oldz = oldPos[2];
	node->flags = flags;
	
	for(i=0; i<peneRes->nSets; i++){
		if(peneRes->effective[i] == 0)continue;
		
		if(i==0){
			maxPowerDis = peneRes->powerDis[i];
			maxPowerSet = i;
		}else{
			maxPowerDis = MAX(maxPowerDis,peneRes->powerDis[i]);
			if(maxPowerDis == peneRes->powerDis[i]){
				maxPowerSet = i;
			}
		}
	}
	if(maxPowerDis < 0){
//		printf("penetrationfirst: maxPowerDis is %e\n",maxPowerDis);
		return(0);
	}
	

	globalOp = 1;
	
	if((globalNode = InsertNode(home,node,nbrNode2,vec, 0,globalOp)) == (Node_t *)NULL) {
		Fatal("insert globalNode failed");
		
	}
	
	switchNode = 0;
	if(GetArmID(home, node, nbrNode1) == -1)switchNode=1;
	
	if(switchNode){
		node1 = globalNode;
		node2 = node;
		
	}else{
		node1 = node;
		node2 = globalNode;
	}
	
	RepositionNode(home, node1Pos, &node1->myTag, globalOp);
	RepositionNode(home, node2Pos, &node2->myTag, globalOp);
	
	peneThough = (fabs(DotProduct(peneRes->outNormal[maxPowerSet], GBPlane->normal)) < 0.995) ? 1 : 0;
	
	if(peneThough){
		node1->inGrain = currInGrain;		// ensure node3 is in current grain.
		node2->inGrain = currInGrain;
		node3 = InsertNode(home,node1,node2,node3Pos, 2, globalOp);		//node3 is in current domain.
	}else{
		node1->inGrain = TBInGrain;		// ensure node3 is in TB.
		node2->inGrain = TBInGrain;
		node3 = InsertNode(home,node1,node2,centorNode12, 1, globalOp);		//node3 is in TB.
	}

	AssignGrain(home, node1, TBInGrain);
	AssignGrain(home, node2, TBInGrain);

	if(node3 == (Node_t *)NULL){
	    Fatal("PenetrationFirst: failed to creat node.\n");
	}

	VECTOR_COPY(burg1,peneRes->outBurg[maxPowerSet]);
	VECTOR_COPY(burg2,peneRes->GBBurg[maxPowerSet]);
	ChangeArmBurg(home, node1, &node3->myTag, burg1[0], burg1[1], burg1[2],
				  peneRes->outNormal[maxPowerSet][X],peneRes->outNormal[maxPowerSet][Y],
				  peneRes->outNormal[maxPowerSet][Z],
				  globalOp,DEL_SEG_NONE);
	ChangeArmBurg(home, node3, &node1->myTag, -burg1[0], -burg1[1], -burg1[2],
				  peneRes->outNormal[maxPowerSet][X],peneRes->outNormal[maxPowerSet][Y],
				  peneRes->outNormal[maxPowerSet][Z],
				  globalOp,DEL_SEG_NONE);
	ChangeArmBurg(home, node3, &node2->myTag, burg1[0], burg1[1], burg1[2],
				  peneRes->outNormal[maxPowerSet][X],peneRes->outNormal[maxPowerSet][Y],
				  peneRes->outNormal[maxPowerSet][Z],
				  globalOp,DEL_SEG_NONE);
	ChangeArmBurg(home, node2, &node3->myTag, -burg1[0], -burg1[1], -burg1[2],
				  peneRes->outNormal[maxPowerSet][X],peneRes->outNormal[maxPowerSet][Y],
				  peneRes->outNormal[maxPowerSet][Z],
				  globalOp,DEL_SEG_NONE);

	if(!IsZeroVector(burg2)){
		InsertArm(home, node1, &node2->myTag, burg2[0], burg2[1],  burg2[2],
				  GBPlane->normal[0],GBPlane->normal[1],GBPlane->normal[2],globalOp);
		InsertArm(home, node2, &node1->myTag, -burg2[0], -burg2[1], -burg2[2],
				  GBPlane->normal[0],GBPlane->normal[1],GBPlane->normal[2],globalOp);
		if((node4 = InsertNode(home,node1,node2, vec, 1,globalOp)) == (Node_t *)NULL) {
			Fatal("PenetrationFirst: Insert node4 faild.");
		}
                node4->flags != NO_MESH_COARSEN;
		
	}
	
	MarkNodeForceObsolete(home, node1);
	MarkNodeForceObsolete(home, node2);
	
	node3->flags |= (NO_COLLISIONS | NODE_RESET_FORCES | NO_MESH_COARSEN);
	
	node1->oldvX = 0.0;
	node1->oldvY = 0.0;
	node1->oldvZ = 0.0;
	
	node2->oldvX = 0.0;
	node2->oldvY = 0.0;
	node2->oldvZ = 0.0;
	
	node1->vX = 0.0;
	node1->vY = 0.0;
	node1->vZ = 0.0;
	
	node2->vX = 0.0;
	node2->vY = 0.0;
	node2->vZ = 0.0;
	
	return(1);
}



int PenetrationSide(Home_t *home, Node_t *node, Grain_t *oldGrain, 
					Grain_t *currGrain, Plane_t *GBPlane){
	
	int		i, globalOp=1;
	real8 		nDir1[3], burg[3], nodeVel[3], oldPos[3], currPos[3];
	real8		glideDir[3], velNormal, currInGrain, newPos[3], vec[3], vector[3];
	real8		currNbrInGrain;
	Node_t		*nbrNode1, *nbrNode2, *nbr, *mergedNode;
	Plane_t		oldGlidePlane, plane;
	Line_t		interLine, pathLine;
        Param_t         *param;
        param = home->param;
	
//	printf("penetration side: ");
		
	nDir1[0] = node->nx[0];
	nDir1[1] = node->ny[0];
	nDir1[2] = node->nz[0];
	
	burg[0] = node->burgX[0];
	burg[1] = node->burgY[0];
	burg[2] = node->burgZ[0];
		
	
	nodeVel[0] = node->vX;
	nodeVel[1] = node->vY;
	nodeVel[2] = node->vZ;
	
	currPos[0] = node->x;
	currPos[1] = node->y;
	currPos[2] = node->z;
	
	oldPos[0] = node->oldx;
	oldPos[1] = node->oldy;
	oldPos[2] = node->oldz;
	
	if(fabs(1.0 - fabs(DotProduct(nDir1,GBPlane->normal))) < 1.0e-2)return(0);
	
	cross(nDir1,GBPlane->normal,interLine.direction);

	NormalizeVec(interLine.direction);

	VECTOR_COPY(oldGlidePlane.normal, nDir1);
	VECTOR_COPY(plane.normal, nDir1);
	VECTOR_COPY(oldGlidePlane.point, oldPos);
	VECTOR_COPY(pathLine.point, oldPos);
	
	currInGrain = InWhichGrain(home, currPos[0], currPos[1], currPos[2]);
/*	
 *	find the dislocation glide direction	
 */
	
	if(Normal(nodeVel) < 1.0)return(0);
	velNormal = DotProduct(nodeVel, nDir1);
	
	glideDir[0] = nodeVel[0] - velNormal*nDir1[0];
	glideDir[1] = nodeVel[1] - velNormal*nDir1[1];
	glideDir[2] = nodeVel[2] - velNormal*nDir1[2];
	
	NormalizeVec(glideDir);	
	
	if(fabs(DotProduct(glideDir,interLine.direction))>0.995) {
		cross(nDir1,interLine.direction,glideDir);
		NormalizeVec(glideDir);
	}
	
	VECTOR_COPY(pathLine.direction, glideDir);
	

	GetIntersection(&oldGlidePlane, GBPlane, &pathLine, &interLine);
	
	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	
	FoldBox(home->param, &nbrNode1->x, &nbrNode1->y, &nbrNode1->z);
	FoldBox(home->param, &nbrNode2->x, &nbrNode2->y, &nbrNode2->z);
	
	if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL) return(0);
	if(SurfaceNode(nbrNode1) && (nbrNode1->inGrain - floor(nbrNode1->inGrain) == 0.5) ||
	   SurfaceNode(nbrNode2) && (nbrNode2->inGrain - floor(nbrNode2->inGrain) == 0.5))return(0);
	
/*	nbrNode2 is always intersection node.	*/
	if(NodeOnGB(home,nbrNode1)){
		Node_t *nodeTemp;
		nodeTemp = nbrNode1;
		nbrNode1 = nbrNode2;
		nbrNode2 = nodeTemp;
	}
	
	if(nbrNode2->numNbrs > 3)return(0);
//	PrintTag(node); 
//	PrintTag(nbrNode1); 
//	PrintTag(nbrNode2); 
//	printf("\n");
	
	currNbrInGrain = InWhichGrain(home,nbrNode1->x,nbrNode1->y,nbrNode1->z);
	
	if(currNbrInGrain == nbrNode1->inGrain || currNbrInGrain == -1){
			
	    vec[0] = nbrNode1->x - currPos[0];
	    vec[1] = nbrNode1->y - currPos[1];
	    vec[2] = nbrNode1->z - currPos[2];
	
            ZImage(home->param, &vec[0], &vec[1], &vec[2]);
	
            VECTOR_COPY(vector,vec);
            NormalizeVec(vec);
            if(fabs(DotProduct(vec,interLine.direction)) < 0.995){
	        InterPointOfLines(currPos,vec,interLine.point,
                                  interLine.direction,newPos);
            }else{
                ZImage(param, vector, vector+1, vector+2);

                plane.point[X] = currPos[X] + 0.5 * vector[X];
                plane.point[Y] = currPos[Y] + 0.5 * vector[Y];
                plane.point[Z] = currPos[Z] + 0.5 * vector[Z];
                FoldBox(param, plane.point,plane.point+1,plane.point+2);

                ClosetPointToGB(param, &plane, GBPlane, newPos);
//              printf("PenetrationSide: approximate intersection 1\n");

            }

	}else{
	    vec[0] = nbrNode1->oldx - node->oldx;
	    vec[1] = nbrNode1->oldy - node->oldy;
	    vec[2] = nbrNode1->oldz - node->oldz;
	
            ZImage(home->param, &vec[0], &vec[1], &vec[2]);
	
	    newPos[0] = node->oldx + vec[0] * 0.5;
	    newPos[1] = node->oldy + vec[1] * 0.5;
	    newPos[2] = node->oldz + vec[2] * 0.5;
	
	    FoldBox(home->param,&newPos[0],&newPos[1],&newPos[2]);
	
            if(fabs(DotProduct(glideDir,interLine.direction)) < 0.995){
	        InterPointOfLines(newPos,glideDir,interLine.point,
                                  interLine.direction,newPos);
            }else{
                VECTOR_COPY(plane.point,newPos);
                ClosetPointToGB(param, &plane, GBPlane, newPos);
//              printf("PenetrationSide: approximate intersection 2\n");
               
            }
	}
	
	FoldBox(home->param,&newPos[0],&newPos[1],&newPos[2]);	
	RepositionNode(home, newPos, &nbrNode2->myTag, 1);
	
	if(InWhichGrain(home,nbrNode1->x,nbrNode1->y,nbrNode1->z) 
	!= nbrNode1->inGrain){
		
	    vec[0] = nbrNode1->oldx - nbrNode2->x;
	    vec[1] = nbrNode1->oldy - nbrNode2->y;
	    vec[2] = nbrNode1->oldz - nbrNode2->z;
	
	}else{
            vec[0] = nbrNode1->x - nbrNode2->x;
            vec[1] = nbrNode1->y - nbrNode2->y;
            vec[2] = nbrNode1->z - nbrNode2->z;
	}
	
	 ZImage(home->param, &vec[0], &vec[1], &vec[2]);
	 
	 newPos[0] = nbrNode2->x + vec[0] * 0.5;
	 newPos[1] = nbrNode2->y + vec[1] * 0.5;
	 newPos[2] = nbrNode2->z + vec[2] * 0.5;
	 
	 FoldBox(home->param,&newPos[0],&newPos[1],&newPos[2]);
	 
         if(InWhichGrain(home,newPos[X],newPos[Y],newPos[Z]) != node->inGrain){
             PullBack(home,node,node->inGrain);
             printf("PenetrationSide: wrong node position.\n");
         }else{
	     RepositionNode(home, newPos, &node->myTag, 1);
         }

	 for(i=0;i<node->numNbrs;i++){
	 
		 if((nbr = GetNeighborNode(home,node,i))!= (Node_t *)NULL)
			 nbr->flags |= NODE_RESET_FORCES;
	 }
	 
	 for(i=0;i<nbrNode2->numNbrs;i++){
	 
		 if((nbr = GetNeighborNode(home,nbrNode2,i))!= (Node_t *)NULL)
			 nbr->flags |= NODE_RESET_FORCES;
	 }
	 
	 return(1);
	
}

int PenetrationLast(Home_t *home, Node_t *node, Grain_t *oldGrain, 
					Grain_t *newGrain, Plane_t *GBPlane){
	real8 		burg[3], vec[3], mergePos[3],TBInGrain;
	Node_t 		*nbrNode1, *nbrNode2, *mergedNode;
	int			mergeStatus,i,j, globalOp,q,thisDomain;
	int crossSlipNode = 0;

	thisDomain = home->myDomain;
	
	burg[0] = node->burgX[0];
	burg[1] = node->burgY[0];
	burg[2] = node->burgZ[0];
	
	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL)
		return(0);
	
	if(SurfaceNode(nbrNode1) && (nbrNode1->inGrain - floor(nbrNode1->inGrain) == 0.5) ||
	   SurfaceNode(nbrNode2) && (nbrNode2->inGrain - floor(nbrNode2->inGrain) == 0.5))return(0);
	globalOp = 1;
	
	FoldBox(home->param, &nbrNode1->x, &nbrNode1->y, &nbrNode1->z);
	FoldBox(home->param, &nbrNode2->x, &nbrNode2->y, &nbrNode2->z);
	
	TBInGrain = nbrNode1->inGrain;

	node->flags |= NODE_RESET_FORCES;
	MarkNodeForceObsolete(home, nbrNode1);
	MarkNodeForceObsolete(home, nbrNode2);
	
	if(nbrNode1->numNbrs > 3 && nbrNode2->numNbrs > 3){
	    printf("Warning PenetrationLast: two neighbor TB node can't mege\n");
	    return(0);
	}

	if(node->numNbrs == 2){       
	    if(fabs(node->nx[0] * node->nx[1] +
	       node->ny[0] * node->ny[1] +
	       node->nz[0] * node->nz[1]) < 0.995){
	        crossSlipNode = 1;
	    }
	}

	if(crossSlipNode){
	    vec[X] = node->x - node->oldx;
	    vec[Y] = node->y - node->oldy;
	    vec[Z] = node->z - node->oldz;
	    ZImage(home->param, &vec[0], &vec[1], &vec[2]);
	    if(Normal(vec) < 1e-6)return(0);

	    NormalizeVec(vec);
	    Line_t line;
	    VECTOR_COPY(line.direction, vec);
	    line.point[X] = node->oldx;
	    line.point[Y] = node->oldy;
	    line.point[Z] = node->oldz;
	    
	    if(!IntersectionOfLinePlane(&line, GBPlane, mergePos)){
	        Fatal("PenetrationLast: intersect failed");
	    }
	}else{
	    vec[0] = nbrNode1->x - nbrNode2->x;
	    vec[1] = nbrNode1->y - nbrNode2->y;
	    vec[2] = nbrNode1->z - nbrNode2->z;
    	
	    ZImage(home->param, &vec[0], &vec[1], &vec[2]);
    	
	    mergePos[0] = nbrNode2->x + vec[0] * 0.5;
	    mergePos[1] = nbrNode2->y + vec[1] * 0.5;
	    mergePos[2] = nbrNode2->z + vec[2] * 0.5;
	}

	if(nbrNode1->numNbrs > 3){
	    mergePos[X] = nbrNode1->x;
	    mergePos[Y] = nbrNode1->y;
	    mergePos[Z] = nbrNode1->z;
	}
	if(nbrNode2->numNbrs > 3){
	    mergePos[X] = nbrNode2->x;
	    mergePos[Y] = nbrNode2->y;
	    mergePos[Z] = nbrNode2->z;
	}

	PullBack(home,node, node->inGrain);
	FoldBox(home->param,&mergePos[0],&mergePos[1],&mergePos[2]);

	MergeNode(home, OPCLASS_SEPARATION, nbrNode1,nbrNode2,
			  mergePos, &mergedNode, &mergeStatus, globalOp);
	
	if((mergeStatus & MERGE_SUCCESS) == 0){
	    printf("PenetrationLast: (%d,%d)--(%d,%d)--(%d,%d) mergeStatus is %d\n",
	           nbrNode1->myTag.domainID,nbrNode1->myTag.index,
	           node->myTag.domainID,node->myTag.index,
	           nbrNode2->myTag.domainID,nbrNode2->myTag.index, mergeStatus);
	    return(0);
	}
	return(1);
}

int PenetrationSurface(Home_t *home, Node_t *node, Grain_t *oldGrain, 
						Grain_t *newGrain, Plane_t *GBPlane, PeneRes_t *peneRes)
{
	
	if(node->constraint <1 || node->constraint >3)return(0);
	
	Node_t 		*nbrNode, *mergedNode, *newNode,*secNbrNode, *node1, *node2, *node3;
	real8		newPos[3], vec[3], CurrNbrInGrain, currInGrain,oldInGrain,TBInGrain;
	real8		oldPos[3], burg[3], sb, pos1[3], pos3[3],currPos[3],eps;
	real8		nDir[3], nSur[3];
	int		mergeStatus, storeConstraint,i,penetration = 0;
	Plane_t		plane;
	Param_t 	*param;
	Line_t		line;
	
	param = home->param;
	
	VECTOR_ZERO(vec);
	
	eps = param->rann;
	
	nbrNode = GetNeighborNode(home, node, 0);
	
	if(nbrNode == (Node_t *)NULL){
	    PrintNode(node);
	    printf("PenetrationSurface: can not found neighbor node\n");
	    return(0);
	}
	
	oldPos[0] = node->oldx;
	oldPos[1] = node->oldy;
	oldPos[2] = node->oldz;

	currPos[X] = node->x;
	currPos[Y] = node->y;
	currPos[Z] = node->z;

	nDir[X] = node->nx[0];
	nDir[Y] = node->ny[0];
	nDir[Z] = node->nz[0];
	
	nSur[X] = (node->constraint == 1);
	nSur[Y] = (node->constraint == 2);
	nSur[Z] = (node->constraint == 3);
	
	cross(nDir,nSur,line.direction);
	if(Normal(line.direction) < 1e-8){
	    Fatal("PenetrationSurface: line direction = 0");
	}
	NormalizeVec(line.direction);

	VECTOR_COPY(line.point, currPos);
	
	if(!IntersectionOfLinePlane(&line, GBPlane, newPos)){
	    Fatal("PenetrationSurface: intersect failed");
	}

	RepositionNode(home, newPos, &nbrNode->myTag, 1);
	AssignConstraint(home, nbrNode, node->constraint);

#if 0
	PullBack(home, node, node->inGrain);
#endif

	ChangeArmBurg(home, node, &nbrNode->myTag, 0.0, 0.0, 0.0,
				  0.0, 0.0, 0.0, 1, DEL_SEG_NONE);
	ChangeArmBurg(home, nbrNode, &node->myTag, 0.0, 0.0, 0.0,
				  0.0, 0.0, 0.0, 1, DEL_SEG_NONE);

	if(node->numNbrs == 0){
	    RemoveNode(home,node,1);
	}else{
	    PrintNode(node);
	    Fatal("the surface node can not removed");
	}

#if 0
	plane.normal[X] = node->nx[0];
	plane.normal[Y] = node->ny[0];
	plane.normal[Z] = node->nz[0];
	burg[X] = node->burgX[0];
	burg[Y] = node->burgY[0];
	burg[Z] = node->burgZ[0];
	
	storeConstraint = node->constraint;
	node->constraint = UNCONSTRAINED;
	
	CurrNbrInGrain = InWhichGrain(home,nbrNode->x,nbrNode->y,nbrNode->z);
	currInGrain = InWhichGrain(home,node->x,node->y,node->z);
	oldInGrain = InWhichGrain(home,node->oldx,node->oldy,node->oldz);
	
	if(CurrNbrInGrain != nbrNode->inGrain || CurrNbrInGrain == -1){
	    vec[X] = nbrNode->oldx - oldPos[0];
	    vec[Y] = nbrNode->oldy - oldPos[1];
	    vec[Z] = nbrNode->oldz - oldPos[2];
	}else{
		vec[X] = nbrNode->x - oldPos[0];
	    vec[Y] = nbrNode->y - oldPos[1];
	    vec[Z] = nbrNode->z - oldPos[2];
	}
	ZImage(param, &vec[0], &vec[1], &vec[2]);
	
	newPos[X] = oldPos[X] + 0.5 * vec[X];
	newPos[Y] = oldPos[Y] + 0.5 * vec[Y];
	newPos[Z] = oldPos[Z] + 0.5 * vec[Z];
	FoldBox(param, newPos, newPos+1, newPos+2);
	VECTOR_COPY(plane.point, newPos);
	ClosetPointToGB(param, &plane, GBPlane, newPos);
	
	newNode = InsertNode(home,node, nbrNode, newPos, 2, 1);
	FormatVector(newPos,"new node position");
	
	if(newNode == (Node_t *)NULL){
		Fatal("PenetrationSurface: insert node failed.");
	}
	
	TBInGrain = GetInGrainFromPos(home,newPos[X],newPos[Y],newPos[Z]);
	if(TBInGrain - floor(TBInGrain) == 0.0){
		Fatal("PenetrationSurface: TBInGrain is %3.1f", TBInGrain);
	}
		
	if(GetArmID(home,node, nbrNode) == -1){
		node1 = node;
		node2 = newNode;
		
	}else{
		node1 = newNode;
		node2 = node;
		printf("switch the arm\n");
	}
	
	
//	RepositionNode(home, newPos, &node2->myTag, 1);
	AssignGrain(home,node2, TBInGrain);
	AssignConstraint(home, node2, UNCONSTRAINED);
	
// node3 pos	
	VECTOR_COPY(plane.point, currPos);
	ClosetPointToGB(param, &plane, GBPlane, pos3);
	
	if(oldGrain->type == MATRIX_GRAIN && newGrain->type == TWIN_GRAIN){
		sb = -1.0;
	}else if(oldGrain->type == TWIN_GRAIN && newGrain->type == MATRIX_GRAIN){
		sb =  1.0;
	}else{
		Fatal("PenetrationSurface: wrong grain type");
	}
	GBPlane->point[X] += sb * eps * GBPlane->normal[X];
	GBPlane->point[Y] += sb * eps * GBPlane->normal[Y];
	GBPlane->point[Z] += sb * eps * GBPlane->normal[Z];
	
	VECTOR_COPY(plane.point, pos3);
	
	for(i=0; i<peneRes->nSets; i++){
		if(i == 1){
			Fatal("PenetrationSurface: can not support muti penetrationresults.");
		}
		VECTOR_COPY(plane.normal,peneRes->outNormal[i]);
		ClosetPointToGB(param, &plane, GBPlane, pos1);
		FormatVector(pos1,"node1 position");
		
		if(fabs(peneRes->outBurg[i][X]) < 1.0e-3 &&
				fabs(peneRes->outBurg[i][Y]) < 1.0e-3 &&
				fabs(peneRes->outBurg[i][Z]) < 1.0e-3){
			Fatal("PenetrationSurface: can not support this set.");
		}
		
		RepositionNode(home, pos1, &node1->myTag, 1);
		ChangeArmBurg(home, node1, &node2->myTag, peneRes->outBurg[i][X],
					  peneRes->outBurg[i][Y], peneRes->outBurg[i][Z],
					  peneRes->outNormal[i][X],peneRes->outNormal[i][Y],peneRes->outNormal[i][Z],
					  1,DEL_SEG_NONE);
		
		ChangeArmBurg(home, node2, &node1->myTag, -peneRes->outBurg[i][X],
					  -peneRes->outBurg[i][Y], -peneRes->outBurg[i][Z],
					  peneRes->outNormal[i][X],peneRes->outNormal[i][Y],peneRes->outNormal[i][Z],
					  1,DEL_SEG_NONE);
		
		printf("node1<-->node2\n");
		PrintNode(node1);
		PrintNode(node2);
		
		node1->flags |= (NO_COLLISIONS | NODE_RESET_FORCES | NO_MESH_COARSEN);
		node2->flags |= (NO_COLLISIONS | NODE_RESET_FORCES | NO_MESH_COARSEN);
		
		node1->oldvX = 0.0;
		node1->oldvY = 0.0;
		node1->oldvZ = 0.0;
		
		node2->oldvX = 0.0;
		node2->oldvY = 0.0;
		node2->oldvZ = 0.0;
		
		node1->vX = 0.0;
		node1->vY = 0.0;
		node1->vZ = 0.0;
		
		node2->vX = 0.0;
		node2->vY = 0.0;
		node2->vZ = 0.0;
		
	
		if(fabs(peneRes->GBBurg[i][X]) < 1.0e-3 &&
			fabs(peneRes->GBBurg[i][Y]) < 1.0e-3 &&
			fabs(peneRes->GBBurg[i][Z]) < 1.0e-3){
			AssignGrain(home,node1, currInGrain);
			AssignConstraint(home, node1, storeConstraint);
		
		}else{
			node3 = InsertNode(home,node1,node2,pos3,2,1);
			if(node3 == (Node_t *)NULL){
				Fatal("PenetrationSurface: insert node failed.");
			}
			
			if(GetArmID(home, node1, node2 )!= -1){
				Fatal("PenetrationSurface: arm switched");
			}
			
			printf("node1<-->node3<-->node2\n");
			PrintNode(node1);
			PrintNode(node3);
			PrintNode(node2);
			
			
			AssignGrain(home,node1, currInGrain);
			AssignConstraint(home, node1, storeConstraint);
			
			ChangeConnection(home, node1, &node3->myTag, &node2->myTag, 1);
			InsertArm(home, node2, &node1->myTag, -peneRes->outBurg[i][X],
					  -peneRes->outBurg[i][Y], -peneRes->outBurg[i][Z],
					  peneRes->outNormal[i][X],peneRes->outNormal[i][Y],peneRes->outNormal[i][Z],1);
			ChangeArmBurg(home, node3, &node1->myTag, 0,0,0,0,0,0,1,DEL_SEG_NONE);
			
			ChangeArmBurg(home, node3, &node2->myTag, peneRes->GBBurg[i][X],
					  peneRes->GBBurg[i][Y], peneRes->GBBurg[i][Z],
					  GBPlane->normal[X],GBPlane->normal[Y],GBPlane->normal[Z],1,DEL_SEG_NONE);
			ChangeArmBurg(home, node2, &node3->myTag, -peneRes->GBBurg[i][X],
					  -peneRes->GBBurg[i][Y], -peneRes->GBBurg[i][Z],
					  GBPlane->normal[X],GBPlane->normal[Y],GBPlane->normal[Z],1,DEL_SEG_NONE);
			AssignGrain(home,node3, TBInGrain);
			AssignConstraint(home, node3, storeConstraint);
			
			node3->flags |= (NO_COLLISIONS | NODE_RESET_FORCES | NO_MESH_COARSEN);
			node3->oldvX = 0.0;
			node3->oldvY = 0.0;
			node3->oldvZ = 0.0;
			
			node3->vX = 0.0;
			node3->vY = 0.0;
			node3->vZ = 0.0;
			
			printf("penetration surfaace success\n");
			PrintNode(node1);
			PrintNode(node2);
			PrintNode(node3);
			
			
		}
		
	}
#endif
	return(1);
	
}

int FindNeiborGBNodes(Home_t *home, Node_t *node, Node_t **nbr1, Node_t **nbr2){
	int 	i,firstGBNode = 1;
	Node_t 	*nbr;
	
	*nbr1 = (Node_t *)NULL;
	*nbr1 = (Node_t *)NULL;
	
	for(i=0; i<node->numNbrs; i++){
		if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL){
			continue;
		}
		
		if(nbr->inGrain - floor(nbr->inGrain) == 0)continue;
		
		if(firstGBNode){
			*nbr1 = nbr;
			firstGBNode = 0;
		}else{
			*nbr2 = nbr;
			
			printf("FindNeiborGBNodes:\n");
			PrintNode(*nbr1);
			PrintNode(*nbr2);
			return(1);
		}
	}
	
	return(0);
}

int SidePeneNode(Home_t *home, Node_t *node){
	int 	i;
	Node_t 	*nbr;
	real8	firstInGrain = -1;
	
	if(node->inGrain - floor(node->inGrain) == 0)return(0);
	if(node->numNbrs > 3 || node->numNbrs < 2)return(0);
	
	for(i=0; i<node->numNbrs; i++){
		if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL){
			continue;
		}
		
		if(nbr->inGrain - floor(nbr->inGrain) == 0.5)continue;
		
		if(firstInGrain == -1){
			firstInGrain = nbr->inGrain;
		}else{
			if(nbr->inGrain != firstInGrain){
				return(1);
			}
		}
	}
	return(0);
}

int FirstPeneNode(Home_t *home, Node_t *node){
	Node_t	*nbrNode1, *nbrNode2;
	int		isScrew, armID;
	
	if(node->numNbrs != 2)return(0);
	if(node->inGrain - floor(node->inGrain) == 0.5)return(0);
	
	nbrNode1 = GetNeighborNode(home, node, 0);
	nbrNode2 = GetNeighborNode(home, node, 1);
	
	if(SidePeneNode(home,nbrNode1) && SidePeneNode(home,nbrNode2)){
		return(1);
	}
	return(0);
}


void HandleGrainBoundary(Home_t *home){
	
	int 		i,j,k,m,numNbrOnGB, globalOp, mergeStatus;
	int	 	grainNode = 0, endNode, newConstraint;
	real8		x,y,z,burgj[3], burgk[3], vec[3], mergPos[3];
	real8		oldInGrain, inGrainPos, corInGrain;	
	real8		burg[3], normal[3], nDir1[3], nDir2[3];
	real8		mergePos[3], distance, newInGrain, newPos[3];
	Param_t		*param;
	Node_t		*node, *nbrNode1, *nbrNode2, *nbr1, *nbr2, *mergedNode;
	Plane_t		GBPlane, norPlane, surPlane;
	Line_t		line;
	Grain_t		*oldGrain, *currGrain;
	PeneRes_t	peneRes;

	
	param = home->param;
	
	peneRes.nSets = MAX_PENE_SETS;
	
	peneRes.effective = (int *) calloc(1, sizeof(int) * MAX_PENE_SETS);
	peneRes.powerDis = (real8 *) calloc(1, sizeof(real8) * MAX_PENE_SETS);
	peneRes.outNormal = (real8 (*)[3]) calloc(1, sizeof(real8) * MAX_PENE_SETS * 3);
	peneRes.outBurg = (real8 (*)[3]) calloc(1, sizeof(real8) * MAX_PENE_SETS * 3);
	peneRes.GBBurg = (real8 (*)[3]) calloc(1, sizeof(real8) * MAX_PENE_SETS * 3);
	
	
	for(i=0;i<home->newNodeKeyPtr;i++){

            if((node = home->nodeKeys[i]) == (Node_t *) NULL)continue;

            numNbrOnGB = 0;
            oldInGrain = node->inGrain;
            FoldBox(param, &node->x, &node->y, &node->z);
            x = node->x; y = node->y; z=node->z;
            	
            inGrainPos = GetInGrainFromPos(home,x,y,z);
            	
            if(oldInGrain == BOUNDARY_NODE)continue;
            	
            if(oldInGrain - floor(oldInGrain) == BOUNDARY_NODE){
            	
                if(inGrainPos != oldInGrain){
                    PullToGB(home, node, oldInGrain);
                }
                continue;
            }
		
            if(node->numNbrs < 1)continue;

            corInGrain = InWhichGrain(home, x, y, z);
		
            if(corInGrain == oldInGrain ){
                continue;
            }


            if(corInGrain == -1){
//              PrintNode(node);
//              printf("Just On TB\n");
	        continue;
            }	
//          param->peneTB = 0;
		
            if(!param->peneTB || node->numNbrs > 2){
                PullBack(home, node, oldInGrain);
                continue;
            }
		
/*	
 *	    Due to different slip planes of two segments of crossSlip node,
 *	    we only allow such node penetrateLast.  	
 */		
            int crossSlipNode = 0;
	    if(node->numNbrs == 2){       
	        if(fabs(node->nx[0] * node->nx[1] +
	                node->ny[0] * node->ny[1] +
                        node->nz[0] * node->nz[1]) < 0.995){
                    
                    crossSlipNode = 1;
                }
	    }
		
            GetGBPlane(home, oldInGrain, corInGrain, node, &GBPlane);

/*	
 *	    There are 4 scenarios: 
 *	  	PenetrationFirst: the first penetrating node.
 *		PenetrationSide: one side penetrating nodes.
 *		PenetrationLast: the last node.
 *		PenetrationSurface: for the surface node.
 */
            oldGrain  = home->grainKeys[(int)oldInGrain];
            currGrain = home->grainKeys[(int)corInGrain];
		
            if(node->constraint < 4 && node->constraint >0){

                if(node->numNbrs != 1){
                    PullBack(home, node, oldInGrain);
                    continue;
                }
			
                nbrNode1 = GetNeighborNode(home,node,0);

                if(nbrNode1 == (Node_t *)NULL){
                    Fatal("HandleGrainBoundary: nbrNode1 not found");
                }
	        if(nbrNode1->inGrain - floor(nbrNode1->inGrain) == 0){
	            PullBack(home, node, oldInGrain);
	            continue;
	        }

	        if(!PenetrationSurface(home, node, oldGrain, currGrain, &GBPlane, &peneRes)){
	            PullBack(home, node, oldInGrain);
	        }

	        continue;
	    }
		    
	    nbrNode1 = GetNeighborNode(home,node,0);
	    nbrNode2 = GetNeighborNode(home,node,1);
			
	    if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL){
	        printf("warning: nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL\n");
	        PullBack(home,node,oldInGrain);
	        continue;
	    }		
	    if(nbrNode1 == nbrNode2) {
	        printf("warning: nbrNode1 == nbrNode2\n");
	        PullBack(home,node,oldInGrain); 
	        continue;
	    }		
	    numNbrOnGB = 0;
	
/*
 *          If the neighbor is not only intersect GB node, it is also junction node,
 *          pull back the penetrated node.
 */
	    if(NodeOnGB(home,nbrNode1)) {
	        if(!(NodeOnGB(home, nbrNode1) & TB_NODE) &&
	           nbrNode1->numNbrs == 3){
	            PullBack(home, node, node->inGrain);
	            continue;
	        } 
	        numNbrOnGB++;	
	    }

	    if(NodeOnGB(home,nbrNode2)){
	        if(!(NodeOnGB(home, nbrNode2) & TB_NODE) &&
	           nbrNode2->numNbrs == 3){
	            PullBack(home, node, node->inGrain);
	            continue;
	        }	        
	        numNbrOnGB++;	
	    }
			
	    if(numNbrOnGB != 2 && crossSlipNode == 1){
/*
 *		CrossSlip node is only allowed to penetrationLast. 
 */
//              printf("warning: numNbrOnGB != 2 && crossSlipNode == 1\n");
                PullBack(home,node, oldInGrain);
                continue;
            }
            	
	    if(numNbrOnGB == 0){
                if(!GetPeneResults(home, node, nbrNode2, &peneRes,
                   oldGrain, currGrain, &GBPlane)){
		    PullBack(home,node, oldInGrain);
		    continue;
		}
            }
/*		
 *	    Node: this energy rule is used for nonscrew dislocation, 
 *	    since we donot know the resultant slip plane for nonscrew 
 *	    dislocation at this time. 
 */

#if ENERGY_RULER
	    if(EnergyRule(home, node, &GBPlane, &peneRes) == 0){
                PullBack(home,node, oldInGrain);
                continue;
	    }
#endif

            burg[0] = node->burgX[0];
            burg[1] = node->burgY[0];
            burg[2] = node->burgZ[0];
            normal[0] = GBPlane.normal[X];
            normal[1] = GBPlane.normal[Y];
            normal[2] = GBPlane.normal[Z];
		
            if((fabs(DotProduct(burg,normal))<0.001) &&
				(numNbrOnGB > 0)){
/*		
 *	        Energy rule for screw dislocation with at least a nbr on GB, so we can
 *		determine the resultant slip plane.	
 */
                if(ScrewEnergyRule(home,node,oldInGrain,
			   corInGrain,&GBPlane)==0){
	
                    PullBack(home,node, oldInGrain);
                    continue;
                }
            }

            if(numNbrOnGB == 0){

                if(!PenetrationFirst(home, node, oldGrain, currGrain, &GBPlane, &peneRes)){
                    PullBack(home,node, oldInGrain);
                    continue;
                }

            }else if(numNbrOnGB == 1){
                if(!PenetrationSide(home, node, oldGrain, currGrain, &GBPlane)){
                    PullBack(home,node, oldInGrain);
                    continue;
                }
			
            }else if(numNbrOnGB == 2){
                if(!PenetrationLast(home, node, oldGrain, currGrain, &GBPlane)){
                    PullBack(home,node, oldInGrain);
                    continue;
                }
            }


	}/*	end of for(i=0;i<home->newNodeKeyPtr ...*/
	
	
	if(param->peneTB == 0) return;

        DebugGB(home,3);

/*	
 *	Below we should handle the GB node that has many segments.
 *	If any GB node has two opposite burgers segments on the same slip
 *	plane and in the same grain, we disconnect such segments.
 */
	for(i=0;i<home->newNodeKeyPtr;i++){
            if((node = home->nodeKeys[i]) == (Node_t *) NULL)continue;
		
/*
 * 	   We'd better check double links here.
 */	
            char name[50];
            if(RemoveDoubleLinks(home,node,0) == 2){
                snprintf(name,50,"RDL%d.2",home->cycle);
                Tecplot(home, name, 0, 0, 0 ,0, 0, 1);
                printf("Unusual RemoveDoubleLinks , 2\n");
            }
            node->flags &= ~NODE_CHK_DBL_LINK; 
	    
	    grainNode = 0;
	    endNode = 0;
	   
            if(node->inGrain -floor(node->inGrain) == 0)grainNode = 1;
	    if(node->constraint == 7 && node->numNbrs > 1){
	        endNode = 1; 
	    }
	    
	    vec[X] = node->x;
	    vec[Y] = node->y;
	    vec[Z] = node->z;

	    VECTOR_ZERO(GBPlane.normal);
	    distance =  ClosedGB(home, vec, &GBPlane);

            if(node->numNbrs < 2){
		if(node->numNbrs == 1 && SurfaceNode(node) && 
		   node->inGrain -floor(node->inGrain) == 0.5){
		    nbr1 = GetNeighborNode(home, node, 0);
		    if(nbr1 != (Node_t *)NULL){
			if(nbr1->inGrain - floor(nbr1->inGrain) == 0){
			    AssignGrain(home, node, nbr1->inGrain);
			}
		    }
		}
		continue;
	    }
	    if(grainNode && node->numNbrs < 4)continue;
            	
            if(node->numNbrs == 2){
                nbr1 = GetNeighborNode(home,node,0);
                nbr2 = GetNeighborNode(home,node,1);
                if(nbr1==(Node_t *)NULL||nbr2==(Node_t *)NULL) continue;
                if(!(nbr1->inGrain - floor(nbr1->inGrain) == 0 &&
                   nbr2->inGrain == nbr1->inGrain)) continue;
            }
			
	    if(grainNode && fabs(distance) > 2*param->rann)continue;
            for(j=0;j<node->numNbrs;j++){
                nbr1 = GetNeighborNode(home,node,j);
                if(nbr1 == (Node_t *)NULL) continue;
                if(nbr1->myTag.domainID!=node->myTag.domainID) continue;
                	
                for(k=j+1;k<node->numNbrs;k++){
                    nbr2 = GetNeighborNode(home,node,k);
                    if(nbr2 == (Node_t *)NULL) continue;
					
		    if(nbr1->inGrain != nbr2->inGrain)continue;
                    if(nbr1 == nbr2) continue;
                    if(nbr2->myTag.domainID != node->myTag.domainID) continue;
					
                        nDir1[0] = node->nx[j];
                        nDir1[1] = node->ny[j];
                        nDir1[2] = node->nz[j];
					
                        nDir2[0] = node->nx[k];
                        nDir2[1] = node->ny[k];
                        nDir2[2] = node->nz[k];

                        if(fabs(DotProduct(nDir1,nDir2))<0.995) continue;
					
                        burg[0] = node->burgX[j] + node->burgX[k];
                        burg[1] = node->burgY[j] + node->burgY[k];
                        burg[2] = node->burgZ[j] + node->burgZ[k];
					
                        if(Normal(burg) > 0.001) continue;
					
                        globalOp = node->myTag.domainID!=nbr1->myTag.domainID||
                                   node->myTag.domainID!=nbr2->myTag.domainID;
					
	                globalOp = 1;
                        ChangeConnection(home,nbr1,&node->myTag,
                                         &nbr2->myTag,globalOp);
                        ChangeConnection(home,nbr2,&node->myTag,
                                         &nbr1->myTag,globalOp);
					
                        ChangeArmBurg(home, node , &nbr1->myTag,
                                      0,0,0,0,0,0,globalOp, DEL_SEG_NONE);
                        ChangeArmBurg(home, node , &nbr2->myTag,
                                     0,0,0,0,0,0,globalOp, DEL_SEG_NONE);
					
//                      printf("HandleGrainBoundary: sepreate GB success!\n");
                       
                        globalOp = nbr1->myTag.domainID != node->myTag.domainID;
	                globalOp = 1;
                        if(RemoveDoubleLinks(home,nbr1,globalOp) == 2){
                            snprintf(name,50,"RDL%d.3",home->cycle);
                            Tecplot(home, name, 0, 0, 0 ,0, 0, 1);
                            printf("Unusual RemoveDoubleLinks , 3\n");
                        }
			
                        globalOp = nbr2->myTag.domainID != node->myTag.domainID;
	                globalOp = 1;
                        if(RemoveDoubleLinks(home,nbr2,globalOp) == 2){
                            snprintf(name,50,"RDL%d.4",home->cycle);
                            Tecplot(home, name, 0, 0, 0 ,0, 0, 1);
                            printf("Unusual RemoveDoubleLinks , 4\n");
                        }
					
                        nbr1->flags &= ~NODE_CHK_DBL_LINK; 
                        nbr2->flags &= ~NODE_CHK_DBL_LINK;
					
                        j--;
                        break;
                    }
		}
		
		if(node->numNbrs == 0) RemoveNode(home,node,globalOp);
	}


	for(i=0;i<home->newNodeKeyPtr;i++){
            if((node = home->nodeKeys[i]) == (Node_t *) NULL)continue;
	    if(node->numNbrs != 2)continue;
	    if(node->inGrain - floor(node->inGrain) == 0)continue;
	    if(SurfaceNode(node))continue;

	    if(fabs(node->nx[0] * node->nx[1] +
	            node->ny[0] * node->ny[1] +
                    node->nz[0] * node->nz[1]) > 0.995){
                continue;
	    }
            nbr1 = GetNeighborNode(home,node,0);
            nbr2 = GetNeighborNode(home,node,1);
            if(nbr1==(Node_t *)NULL||nbr2==(Node_t *)NULL) continue;
	    if(!SurfaceNode(nbr1) && !SurfaceNode(nbr2))continue;
	    if(SurfaceNode(nbr1) && SurfaceNode(nbr2)){
	        AssignConstraint(home, node, nbr1->constraint);
	    }

	    if(SurfaceNode(nbr1)){
		if(NodeDistance(home, node, nbr1) > param->minSeg)continue;
	        vec[X] = nbr2->x; 
	        vec[Y] = nbr2->y; 
	        vec[Z] = nbr2->z; 
	        newInGrain = nbr2->inGrain; 
		newConstraint = nbr1->constraint;
	        norPlane.normal[X] = node->nx[1];
	        norPlane.normal[Y] = node->ny[1];
	        norPlane.normal[Z] = node->nz[1];
		
	        surPlane.point[X] = nbr1->x;
	        surPlane.point[Y] = nbr1->y;
	        surPlane.point[Z] = nbr1->z;
	    }else{
		if(NodeDistance(home, node, nbr2) > param->minSeg)continue;
	        vec[X] = nbr1->x; 
	        vec[Y] = nbr1->y; 
	        vec[Z] = nbr1->z; 
	        newInGrain = nbr1->inGrain; 
		newConstraint = nbr2->constraint;
	        norPlane.normal[X] = node->nx[0];
	        norPlane.normal[Y] = node->ny[0];
	        norPlane.normal[Z] = node->nz[0];
	        surPlane.point[X] = nbr2->x;
	        surPlane.point[Y] = nbr2->y;
	        surPlane.point[Z] = nbr2->z;
	    }
	    VECTOR_ZERO(GBPlane.normal);
	    distance =  ClosedGB(home, vec, &GBPlane);

	    GBPlane.point[X] += Sign(distance) * param->rann * GBPlane.normal[X];
	    GBPlane.point[Y] += Sign(distance) * param->rann * GBPlane.normal[Y];
	    GBPlane.point[Z] += Sign(distance) * param->rann * GBPlane.normal[Z];

	    norPlane.point[X] = node->x;
	    norPlane.point[Y] = node->y;
	    norPlane.point[Z] = node->z;

            if(!IntersectionOfPlanePlane(&norPlane,&GBPlane,&line)){
		PrintNode(node);
		PrintNode(nbr1);
		PrintNode(nbr2);
	        Fatal("HandleGrainBoundary: intersection of planes failed!");
	    }

	    VECTOR_ZERO(surPlane.normal);
	    surPlane.normal[newConstraint - 1] = 1.0;
	    if(!IntersectionOfLinePlane(&line, &surPlane, newPos)){
	        Fatal("HandleGrainBoundary: intersection of line and plane failed");
	    }

	    RepositionNode(home,newPos,&node->myTag,1);
	    AssignGrain(home, node, newInGrain);
	    AssignConstraint(home, node, newConstraint);
	    printf("node (%d,%d) new\n",node->myTag.domainID, node->myTag.index);
	}

	for(i=0;i<home->newNodeKeyPtr;i++){
            if((node = home->nodeKeys[i]) == (Node_t *)NULL)continue;
            if(node->inGrain - floor(node->inGrain) == 0.5)continue;
	    if(SurfaceNode(node))continue;
            if(node->inGrain != InWhichGrain(home, node->x, node->y, node->z)){
	        PullBack(home,node, node->inGrain);
                printf("HandleGrainBoundary: a penetrated not coopered\n");
            }
        }
	return;

}


void DoResistance(Home_t *home, Node_t *node, real8 *nForce){
	real8 	distance, pos[3],dir[3], factor, burgMag;
	real8	GBStress, lDir[3], nDir[3], gDir[3], length;
	real8 	burg[3], fricMag, dot, armGBForce[3];
	real8	rann, armf[3], forceGB[3], segCentor[3];
	real8	oriNForce[3], gForce[3];
	int	i, GBForce, fricForce, armGrain;
	Plane_t GBPlane;
	Node_t	*nbr;
	Param_t *param;
	
	param = home->param;
	pos[X] = node->x;
	pos[Y] = node->y;
	pos[Z] = node->z;
	VECTOR_ZERO(GBPlane.normal);
	VECTOR_ZERO(forceGB);

	nForce[X] = node->fX;
	nForce[Y] = node->fY;
	nForce[Z] = node->fZ;
	VECTOR_COPY(oriNForce, nForce);

	rann = 3.0 * param->rann;

	for(i=0; i<node->numNbrs; i++){
	    nbr = GetNeighborNode(home,node,i);
	    if(nbr == (Node_t *)NULL)continue;

	    armGrain = ChangeOrientation(home, node, nbr,0);
	    if(armGrain - floor(armGrain) == 0.5)continue;

	    lDir[X] = nbr->x - node->x;
	    lDir[Y] = nbr->y - node->y;
	    lDir[Z] = nbr->z - node->z;
	    ZImage(param, lDir, lDir+1, lDir+2);

	    segCentor[X] = node->x + 0.5 * lDir[X];
	    segCentor[Y] = node->y + 0.5 * lDir[Y];
	    segCentor[Z] = node->z + 0.5 * lDir[Z];

	    length = Normal(lDir);

	    if(length < 1e-8)continue;

	    NormalizeVec(lDir);
	    nDir[X] = node->nx[i];
	    nDir[Y] = node->ny[i];
	    nDir[Z] = node->nz[i];
	    cross(lDir, nDir, gDir);

	    if(Normal(gDir) < 1e-5) continue;
	    NormalizeVec(gDir);

	    burg[X] = node->burgX[i];
	    burg[Y] = node->burgY[i];
	    burg[Z] = node->burgZ[i];

	    armf[X] = node->armfx[i];
	    armf[Y] = node->armfy[i];
	    armf[Z] = node->armfz[i];

	    dot = DotProduct(armf, gDir);
	    gForce[X] = dot * gDir[X];
	    gForce[Y] = dot * gDir[Y];
	    gForce[Z] = dot * gDir[Z];

	    burgMag = sqrt(node->burgX[i] * node->burgX[i] +
			   node->burgY[i] * node->burgY[i] +
			   node->burgZ[i] * node->burgZ[i]);

	    fricMag = 0.5e0 * length * burgMag * param->fricStress;

	    if(Normal(gForce) <= fricMag){
		nForce[X] -= gForce[X];
		nForce[Y] -= gForce[Y];
		nForce[Z] -= gForce[Z];
		continue;
	    }

	    NormalizeVec(gForce);
	    nForce[X] -= (fricMag * gForce[X]);
	    nForce[Y] -= (fricMag * gForce[Y]);
	    nForce[Z] -= (fricMag * gForce[Z]);

	    armf[X] -= (fricMag * gForce[X]);
	    armf[Y] -= (fricMag * gForce[Y]);
	    armf[Z] -= (fricMag * gForce[Z]);

	    if(param->polycrystal == 0)continue;

	    GBForce = 1;
	    distance =  ClosedGB(home, segCentor, &GBPlane);
	    GBForce &= (fabs(distance) < 2.0e0 * rann);
	    GBForce &= (fabs(distance) > 1.0e-3);

/*
 *	    screw dislocation can easily cut though GB.
 */
	    dot = fabs(DotProduct(GBPlane.normal, burg));
	    GBForce &= (dot > 1.0e-1);

	    if(GBForce == 0)continue;

	    dot = DotProduct(armf, gDir);
	    gForce[X] = dot * gDir[X];
	    gForce[Y] = dot * gDir[Y];
	    gForce[Z] = dot * gDir[Z];

	    dir[X] = Sign(distance) * GBPlane.normal[X];
	    dir[Y] = Sign(distance) * GBPlane.normal[Y];
	    dir[Z] = Sign(distance) * GBPlane.normal[Z];
	    if(DotProduct(gForce, dir) >= 0.0)continue;

	    GBStress = param->GBStrength * pow(2.7182,-(distance*distance/rann/rann));
	    GBStress *= (length * burgMag * 0.5e0);
	    if(Normal(gForce) <= GBStress){
		nForce[X] -= gForce[X];
		nForce[Y] -= gForce[Y];
		nForce[Z] -= gForce[Z];
		continue;
	    }

	    NormalizeVec(gForce);
	    nForce[X] -= (GBStress * gForce[X]);
	    nForce[Y] -= (GBStress * gForce[Y]);
	    nForce[Z] -= (GBStress * gForce[Z]);

	    armf[X] -= (GBStress * gForce[X]);
	    armf[Y] -= (GBStress * gForce[Y]);
	    armf[Z] -= (GBStress * gForce[Z]);
	}

	if(DotProduct(nForce, oriNForce) < -1.0){
	    VECTOR_ZERO(nForce);
	}

	return;

}

void   ApplyPosConstraint(Param_t *param, Node_t *node){
	real8 disVec[3], maxDis= -1E50;
	int   i, constraintSurface;
	real8 pos[3];

	pos[X] = node->x;
	pos[Y] = node->y;
	pos[Z] = node->z;

	VECTOR_ZERO(disVec);

	if(param->xBoundType == Free &&  fabs(pos[X])>param->xBoundMax) { 
	    disVec[X] = fabs(pos[X]) - param->xBoundMax;
	}
	if(param->yBoundType == Free &&  fabs(pos[Y])>param->yBoundMax) { 
	    disVec[Y] = fabs(pos[Y]) - param->yBoundMax;
	}
	if(param->zBoundType == Free &&  fabs(pos[Z])>param->zBoundMax) { 
	    disVec[Z] = fabs(pos[Z]) - param->zBoundMax;
	}

	for(i=0; i<3; i++){
	    if(disVec[i] > maxDis){
	        maxDis = disVec[i];
	        constraintSurface = i + 1;
	    }
	}

	if(maxDis > 0){
	    node->constraint = constraintSurface;
	}
	if(maxDis < -5.0 && node->constraint >= 1 && 
	   node->constraint <=3 && node->numNbrs > 1){
	    node->constraint = 0;
	}
	return;
}


void   	AdjustSurfaceNode(Home_t *home)       //adjust surface node volecity
	{

	int 	i,j,switchArm, numConstraints;
	real8	nodePos[3], nbrPos[3];
	real8	interPos[2][3], oldPos[3], bakInter[3];
	real8	newPos[3], burg0[3], burgj[3];
 	Param_t *param;
	Node_t  *node, *nbr, *newNode, *temNode, *node1, *node2;
	int 	nodecst, surfaceNbr, constraint;
	Line_t	line;
	Plane_t	plane;
	param = home->param;
  
	for (i = 0; i < home->newNodeKeyPtr; i++){

        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

            node->flags &= ~NOT_HANDLE_SEG;
	    nodePos[X] = node->x;
	    nodePos[Y] = node->y;
	    nodePos[Z] = node->z;
	    ApplyPosConstraint(param, node);

	    if(SurfaceNode(node) && node->numNbrs > 1){
	        if(InWhichGrain(home, nodePos[X], nodePos[Y], nodePos[Z]) != node->inGrain &&
	           node->inGrain - floor(node->inGrain) == 0){
		    node->x = node->oldx;
		    node->y = node->oldy;
		    node->z = node->oldz;

	            nodePos[X] = node->x;
	            nodePos[Y] = node->y;
	            nodePos[Z] = node->z;
	            ApplyPosConstraint(param, node);
	        }
	    }

	    if(SurfaceNode(node) || node->numNbrs < 2) continue;
#if 0
	    surfaceNbr = 0;
	    if(node->numNbrs == 2){
	        for(j=0;j<node->numNbrs;j++){
		    if((nbr=GetNeighborNode(home,node,j))==(Node_t *)NULL) continue;
		    if(nbr->constraint >= 1 && nbr->constraint <= 3) surfaceNbr++;
		    if(surfaceNbr > 1 ) {
		        node->constraint = nbr->constraint;
		        break;
		    }
	        }
	    }
#endif
	}

	node = home->ghostNodeQ;

        while (node) {
	    nodePos[X] = node->x;
	    nodePos[Y] = node->y;
	    nodePos[Z] = node->z;
	    ApplyPosConstraint(param, node);

	    if(SurfaceNode(node) && node->numNbrs > 1){
	        if(InWhichGrain(home, nodePos[X], nodePos[Y], nodePos[Z]) != node->inGrain &&
	           node->inGrain - floor(node->inGrain) == 0){
		    node->x = node->oldx;
		    node->y = node->oldy;
		    node->z = node->oldz;

	            nodePos[X] = node->x;
	            nodePos[Y] = node->y;
	            nodePos[Z] = node->z;
	            ApplyPosConstraint(param, node);
	        }
	    }

            node = node->next;
	}

	for (i = 0;i < home->newNodeKeyPtr; i++){

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
	    if(SurfaceNode(node))continue;
	    nodePos[X] = node->x;
	    nodePos[Y] = node->y;
	    nodePos[Z] = node->z;

	    if(param->xBoundType == Free &&  fabs(nodePos[X]) >= param->xBoundMax  &&
	       param->yBoundType == Free &&  fabs(nodePos[Y]) >= param->yBoundMax  &&
	       param->zBoundType == Free &&  fabs(nodePos[Z]) >= param->zBoundMax)continue;

	    for(j=0; j < node->numNbrs; j++){

	        nbr = GetNeighborNode(home, node, j);
	        if(nbr == (Node_t *)NULL)continue;
		if(nbr->numNbrs < 2)continue;
	        if(!SurfaceNode(nbr))continue;
	        nbrPos[X] = nbr->x;
	        nbrPos[Y] = nbr->y;
	        nbrPos[Z] = nbr->z;

	        if(param->xBoundType == Free &&  fabs(nbrPos[X]) <= param->xBoundMax + 1.0E-10 &&
	           param->yBoundType == Free &&  fabs(nbrPos[Y]) <= param->yBoundMax + 1.0E-10 &&
	           param->zBoundType == Free &&  fabs(nbrPos[Z]) <= param->zBoundMax + 1.0E-10)continue;

		nodecst = node->constraint;
	        constraint = nbr->constraint;
		VECTOR_COPY(line.point, nodePos);
	        line.direction[X] = nbrPos[X] - nodePos[X];
	        line.direction[Y] = nbrPos[Y] - nodePos[Y];
	        line.direction[Z] = nbrPos[Z] - nodePos[Z];

	        if(Normal(line.direction) < 1.0e-3){
	           continue;
		}

	        NormalizeVec(line.direction);		
		
		VECTOR_ZERO(plane.normal);
		plane.normal[nbr->constraint - 1] = 1.0;

	        VECTOR_ZERO(plane.point);

	        if(nbr->constraint == 1){
	            plane.point[X] = Sign(nbrPos[X]) * param->xBoundMax;
	        }
	        else if(nbr->constraint == 2){
	            plane.point[Y] = Sign(nbrPos[Y]) * param->yBoundMax;
	        }
	        else if(nbr->constraint == 3){
	            plane.point[Z] = Sign(nbrPos[Z]) * param->zBoundMax;
	        }

	        VECTOR_ZERO(interPos[0]);
	        VECTOR_ZERO(interPos[1]);
	        if(!IntersectionOfLinePlane(&line, &plane, interPos[0])){
                    VECTOR_COPY(interPos[0], nbrPos);
                }

		if(fabs(interPos[0][X] - nbrPos[X]) < 1e-10 &&
		   fabs(interPos[0][Y] - nbrPos[Y]) < 1e-10 &&
		   fabs(interPos[0][Z] - nbrPos[Z]) < 1e-10){
		    interPos[0][X] = 0.99999 * nbrPos[X] + 0.00001 * nodePos[X]; 
		    interPos[0][Y] = 0.99999 * nbrPos[Y] + 0.00001 * nodePos[Y]; 
		    interPos[0][Z] = 0.99999 * nbrPos[Z] + 0.00001 * nodePos[Z]; 
		}

	        if((InWhichGrain(home, interPos[0][X], interPos[0][Y], interPos[0][Z]) != node->inGrain &&
	           node->inGrain - floor(node->inGrain) == 0) || 
		   (GetInGrainFromPos(home,interPos[0][X], interPos[0][Y], interPos[0][Z]) != node->inGrain &&
		   node->inGrain - floor(node->inGrain) == 0.5)){
		    PrintNode(node);
		    PrintNode(nbr);
		    FormatVector(interPos[0],"interPos");
	            printf("AdjustSurfaceNode: intersetion penetrated grain\n");
		    VECTOR_COPY(interPos[0], nbrPos);
	        }

		numConstraints = CheckReposition(home, node, interPos[0]);
	        if(numConstraints == 2){
	            real8	vec1[3], vec2[3];

	            line.direction[X] = interPos[0][X] - node->x;
	            line.direction[Y] = interPos[0][Y] - node->y;
	            line.direction[Z] = interPos[0][Z] - node->z;
	            VECTOR_COPY(bakInter,interPos[0]);
	            
	            if(Normal(line.direction) > 1.0e-3){
	                NormalizeVec(line.direction);		
	                if(!IntersectionOfLineBox(home->param, &line, interPos)){
			    PrintNode(node);
			    PrintNode(nbr);
	                    FormatVector(line.direction,"line.direction");
	                    FormatVector(line.point,"line.point");
	                    FormatVector(plane.normal,"plane.normal");
	                    FormatVector(plane.point,"point");
			    FormatVector(interPos[0],"interPos0");
			    FormatVector(interPos[1],"interPos0");
                            printf("AdjustSurfaceNode: numConstrins = 2 and can not find the inertsection.");
	                    continue;
                        }
	                vec1[X] = interPos[0][X] - nodePos[X];
	                vec1[Y] = interPos[0][Y] - nodePos[Y];
	                vec1[Z] = interPos[0][Z] - nodePos[Z];
	                vec2[X] = interPos[1][X] - nodePos[X];
	                vec2[Y] = interPos[1][Y] - nodePos[Y];
	                vec2[Z] = interPos[1][Z] - nodePos[Z];
	                if(Normal(vec1) > Normal(vec2))VECTOR_COPY(interPos[0], interPos[1]);
	                if(param->maxSeg < Normal(vec1) && param->maxSeg < Normal(vec2))continue;
	            }
		    printf("CR: %d, {%f,%f,%f}-->{%f,%f,%f},\n nbr {%f,%f,%f}, bakinter {%f,%f,%f}\n lineDir {%f,%f,%f}\n",
	                   numConstraints, node->x,node->y,node->z,
		           interPos[0][X], interPos[0][Y], interPos[0][Z],
	                   nbr->x, nbr->y, nbr->z,
	                   bakInter[X], bakInter[Y], bakInter[Z],
	                   line.direction[X], line.direction[Y], line.direction[Z]);
	        }
	        if( numConstraints < 3){
		    RepositionNode(home, interPos[0], &node->myTag, 1);
	            AssignConstraint(home, node, constraint);
	            MarkNodeForceObsolete(home, node);
		}

#if 0
	        printf("(%d,%d)<-(%d,%d)->(%d,%d)\n",node->myTag.domainID,node->myTag.index,
	               newNode->myTag.domainID, newNode->myTag.index,nbr->myTag.domainID, nbr->myTag.index);
	        FormatVector(nodePos,"nodePos");
	        FormatVector(interPos,"interPos");
	        FormatVector(nbrPos,"nbrPos");
#endif
	    } /*for(j=0; j< ...)*/
	}/*for(i=0; i< ...)*/	

	return;
}

int CheckReposition(Home_t *home, Node_t *node, real8 *newPos){

        int 	*numGlideConstraints = NULL, NM;
        real8 	(*glideConstraints)[3], GC[3][3], vec[3];

	NM = 0;
        numGlideConstraints = &NM;
        *numGlideConstraints = 0;
	
        glideConstraints = GC;
	vec[X] = newPos[X] - node->x;
	vec[Y] = newPos[Y] - node->y;
	vec[Z] = newPos[Z] - node->z;
	ZImage(home->param, vec, vec+1, vec+2);

	GetGlideConstraintList(home, node, numGlideConstraints, glideConstraints);
	ApplyConstraintsToVelocity(*numGlideConstraints,glideConstraints, vec);

	newPos[X] = node->x + vec[X];
	newPos[Y] = node->y + vec[Y];
	newPos[Z] = node->z + vec[Z];
	return(*numGlideConstraints);
}
