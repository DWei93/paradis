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


Node_t *InsertNode(Home_t *home, Node_t *node1, Node_t *node2, real8 *p, int globalOp, int opClass)
{

    int     armID, splitStatus;
    real8   nodePos[3], nodeVel[3];
    Node_t	*splitNode1,*splitNode2,*newNode;

    armID = GetArmID(home, node1, node2);
    if(armID < 0){
        Fatal("in %s at %d", __FILE__, __LINE__);
    }

    nodePos[0] = node1->x;
    nodePos[1] = node1->y;
    nodePos[2] = node1->z;
    
    nodeVel[0] = 0.5*(node1->vX + node2->vX);
    nodeVel[1] = 0.5*(node1->vY + node2->vY);
    nodeVel[2] = 0.5*(node1->vZ + node2->vZ);

    splitStatus = SplitNode(home, opClass,
                     node1, nodePos, p,
                     nodeVel, newVel,
                     1, &armIndex,
                     globalOp, &splitNode1,
                     &splitNode2, 0);

    if(splitStatus != SPLIT_SUCCESS)return((Node_t *)NULL);

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
real8 FGlidePerLength(Home_t *home, Node_t *node, Plane_t *GBPlane){

    int	    i;
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


real8 SegementEnergy(Param_t *param, real8 cosTheta, real8 bMag){
	static real8 pi = 3.1415926535897932385;
	
	return(param->shearModulus * bMag * bMag * log(param->outRadius) *
		   (1.0 - param->pois * cosTheta * cosTheta)/
		   (4.0 * pi * (1.0 - param->pois)) );
		  
}

int EnergyRule(Home_t *home, Node_t *node, Plane_t *GBPlane, Pene_t *pene){
	
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
		inBurg[0] = node->burgX[i];
		inBurg[1] = node->burgY[i];
		inBurg[2] = node->burgZ[i];
		bMag = Normal(inBurg); 
		
		NormalizeVec(inBurg);
		vec[0] = nbr->x - node->x;
		vec[1] = nbr->y - node->y;
		vec[2] = nbr->z - node->z;
		ZImage(param, vec, vec+1, vec+2);
		NormalizeVec(vec);
		
		cosTheta = DotProduct(vec, inBurg);
		segEnergy = SegementEnergy(param, cosTheta, bMag);
		inEnergy += (segEnergy/2.0);
	}

	normal[0] = node->nx[1];
	normal[1] = node->ny[1];
	normal[2] = node->nz[1];
	temVec[0] = node->burgX[1];
	temVec[1] = node->burgY[1];
	temVec[2] = node->burgZ[1];
	burgMag = Normal(temVec);
	
	cross(normal, GBPlane->normal, vec);
	NormalizeVec(vec);
	
	for(i=0; i<pene->nSets; i++){
		outEnergy = 0.0;
		if(pene->effective[i] == 0)continue;
		cosTheta = DotProduct(pene->outBurg[i],vec);
		bMag = Normal(pene->outBurg[i]);
		segEnergy = SegementEnergy(param, cosTheta, bMag);
		outEnergy += segEnergy;
		
		cosTheta = DotProduct(pene->GBBurg[i],vec);
		bMag = Normal(pene->GBBurg[i]);
		segEnergy = SegementEnergy(param, cosTheta, bMag);
		outEnergy += segEnergy;
		
		outEnergy += burgMag * param->TBEnergy/param->burgMag;
		
		if(FGlidePerLength(home,node,GBPlane) * burgMag < outEnergy - inEnergy){
			pene->effective[i] = 0;
		
		}
	}
#endif 
	
	temVec[0] = node->burgX[0];
	temVec[1] = node->burgY[0];
	temVec[2] = node->burgZ[0];
	
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


void FreePenStr(Pene_t *pene){
    if(pene->nSets > 0){
        int i, j;
        for(i=0; i<pene->nSets; i++){
            if(pene->res[i].nPlanes > 0){
                free(pene->res[i].nDissOnPlane);
                free(pene->res[i].inGrain);
                free(pene->res[i].passGB);
                free(pene->res[i].normal);

                pene->res[i].nDissOnPlane = NULL;
                pene->res[i].inGrain = NULL;
                pene->res[i].passGB = NULL;
                pene->res[i].normal = NULL;
            }

            if(pene->res[i].nDiss > 0){
                free(pene->res[i].burg);
                free(pene->res[i].pos);

                pene->res[i].burg = NULL;
                pene->res[i].pos = NULL;
            }

            pene->res[i].nPlanes = 0;
            pene->res[i].nDiss = 0;
            pene->res[i].effective = 0;
            pene->res[i].powerDis = 0;
        }
    }

    free(pene->res);
    pene->res = NULL;
    pene->nSets = 0;

    pene->gb = NULL;
    pene->node = NULL;

    VECTOR_ZERO(pene->line1.direction);
    VECTOR_ZERO(pene->line2.direction);
    VECTOR_ZERO(pene->line1.point);
    VECTOR_ZERO(pene->line2.point);
    return;
}

void InitPenStr(Pene_t *pene)
{
    int i, j; 
    pene->gb = NULL;
    pene->node = NULL;

    pene->nSets = MAX_PENE_SETS;

    pene->res = (Res_t *)calloc(MAX_PENE_SETS, sizeof(Res_t));

    for(i=0; i<pene->nSets; i++){
        pene->res[i].nPlanes = MAX_PENE_PLANES;
        pene->res[i].nDiss = MAX_PENE_DISS;
        pene->res[i].passGB = (bool *)calloc(pene->res[i].nPlanes, sizeof(bool));
        pene->res[i].nDissOnPlane = (int *)calloc(pene->res[i].nPlanes, sizeof(int));
        pene->res[i].inGrain = (int *)calloc(pene->res[i].nPlanes, sizeof(int));
        pene->res[i].normal = (real8 (*)[3])calloc(1,3*pene->res[i].nPlanes*sizeof(real8));

        pene->res[i].pos  = (real8 (*)[3])calloc(1, 3*pene->res[i].nDiss*sizeof(real8));
        pene->res[i].burg = (real8 (*)[3])calloc(1, 3*pene->res[i].nDiss*sizeof(real8));
    }
    
	return;
}


void InitPeneRes(Res_t *res, int nPlanes, int nDiss){
    if(res->nPlanes > 0){
        free(res->passGB);        
        free(res->inGrain);        
        free(res->normal);        
        free(res->nDissOnPlane);        
        
        res->passGB = NULL;
        res->inGrain = NULL;
        res->normal = NULL;
        res->nDissOnPlane = NULL;
    }
    res->nPlanes = nPlanes;

    if(res->nDiss > 0){
        free(res->pos);
        free(res->burg);
        
        res->pos = NULL;
        res->burg = NULL;
    }
    res->nDiss = nDiss;
    
    res->effective = 1;
    res->powerDis = -1.0E50;
    res->passGB = (bool *)calloc(1, nPlanes*sizeof(bool));
    res->nDissOnPlane = (int *)calloc(1, nPlanes*sizeof(int));
    res->inGrain = (int *)calloc(1, nPlanes*sizeof(int));
    res->normal = (real8 (*)[3])calloc(1, nPlanes*sizeof(real8)*3);

    res->pos = (real8 (*)[3])calloc(1, nDiss*sizeof(real8)*3);
    res->burg = (real8 (*)[3])calloc(1, nDiss*sizeof(real8)*3);
    return;
}

/********************************************************************************
 *   node1:     penetrated node
 *   node2:     nbr node
 * ******************************************************************************/
void FindIntersectionWithGB(Home_t *home, Node_t *node1, Node_t *node2, 
                            Plane_t *GBPlane, Line_t *line)
{
    int         armID, type;
    Param_t     *param;
    real8       n1Pos[3], n2Pos[3], n1OldPos[3], n2OldPos[3];
    real8       normal[3], centor[3];
    real8       vec1[3], vec2[3], vec3[3], vec4[3];
    real8       sign1, sign2, dis1, dis2, dis;

    param = home->param;

    armID = GetArmID(home, node1, node2);
    if(armID < 0) Fatal("In file %s at line %d", __FILE__, __LINE__);    

    normal[0] = node1->nx[armID];
    normal[1] = node1->ny[armID];
    normal[2] = node1->nz[armID];

    n1Pos[0] = node1->x;    
    n1Pos[1] = node1->y;    
    n1Pos[2] = node1->z;    

    n2Pos[0] = node2->x;    
    n2Pos[1] = node2->y;    
    n2Pos[2] = node2->z;    

    n1OldPos[0] = node1->oldx;    
    n1OldPos[1] = node1->oldy;    
    n1OldPos[2] = node1->oldz;    

    n2OldPos[0] = node2->oldx;    
    n2OldPos[1] = node2->oldy;    
    n2OldPos[2] = node2->oldz;    

    PBCPOSITION(param, n1Pos[0], n1Pos[1], n1Pos[2], 
                       n2Pos, n2Pos+1, n2Pos+2);
    PBCPOSITION(param, n1Pos[0], n1Pos[1], n1Pos[2], 
                       n1OldPos, n1OldPos+1, n1OldPos+2);
    PBCPOSITION(param, n1OldPos[0], n1OldPos[1], n1OldPos[2], 
                       n2OldPos, n2OldPos+1, n2OldPos+2);
    centor[0] = (n1Pos[0] +  n2Pos[0] + n1OldPos[0] + n2OldPos[0])*0.25;
    centor[1] = (n1Pos[1] +  n2Pos[1] + n1OldPos[1] + n2OldPos[1])*0.25;
    centor[2] = (n1Pos[2] +  n2Pos[2] + n1OldPos[2] + n2OldPos[2])*0.25;

    type = PlanePlaneIntersection(normal, centor, GBPlane->normal, 
           GBPlane->point, line->direction, line->point);
    
    if(type != LINE_INTERSECTION){
        Fatal("In file %s at line %d", __FILE__, __LINE__);        
    }

    PointPlaneIntersection(n1Pos, GBPlane->normal, GBPlane->point, vec1, &dis1);
    PointPlaneIntersection(n2Pos, GBPlane->normal, GBPlane->point, vec1, &dis2);
    sign1 = Sign(dis1);
    sign2 = Sign(dis2);
    if(sign1*sign2 > 0){
        vec1[0] = (n1Pos[0] +  n2Pos[0] - n1OldPos[0] - n2OldPos[0]);
        vec1[1] = (n1Pos[1] +  n2Pos[1] - n1OldPos[1] - n2OldPos[1]);
        vec1[2] = (n1Pos[2] +  n2Pos[2] - n1OldPos[2] - n2OldPos[2]);
        VECTOR_COPY(vec2, centor);
    }else{
        vec1[0] = n1Pos[0] - n2Pos[0];
        vec1[1] = n1Pos[1] - n2Pos[1];
        vec1[2] = n1Pos[2] - n2Pos[2];

        vec2[0] = (n1Pos[0] + n2Pos[0])*0.5;
        vec2[1] = (n1Pos[1] + n2Pos[1])*0.5;
        vec2[2] = (n1Pos[2] + n2Pos[2])*0.5;
    }
    if(Normal(vec1) < 1.0E-4)Fatal("In file %s at line %d", __FILE__, __LINE__);
    NormalizeVec(vec1);

    type = LineLineIntersection(vec1, vec2, line->direction, line->point, vec3, vec4, &dis);
    
    if(fabs(dis) > 10){
        printf("type %d, dis %f\n", type, dis);
        FormatVector(vec1, "vec1");
        printf("Warning in file %s at line %d\n", __FILE__, __LINE__);
    }
    
    VECTOR_COPY(line->point, vec4);
    return;
}

int GetPenResults(Home_t *home, Node_t *node, Pene_t *pene, GB_t *gb){
	
	real8		pos[3], burg[3], normal[3], crystalNormal[3], burg1[3], burg2[3];
	real8		crystalBurg[3], dot, vec[3], outNormal[3];
    real8       sb, TBBurgCryst[3], centorPos[3];
    real8       vec1[3], vec2[3], eps = 0.1, dis;
    int         currInGrain, currOnGB, TBID;
    Line_t      line;
	int 		isScrew = 0, disID, i, j, k, normalID, type;
    Node_t      *node1, *node2;
    Grain_t     *oldGrain;
    Plane_t     GBPlane;
	
    FreePenStr(pene);

    pene->node = node;
    pene->gb = gb;

    oldGrain  = home->grainKeys[node->inGrain];
	node1 = GetNeighborNode(home,node,0);
	node2 = GetNeighborNode(home,node,1);
    if(node1 == (Node_t *)NULL || node2 == (Node_t *)NULL){
        Fatal("In file %s at line %d", __FILE__, __LINE__);
    }

    VECTOR_COPY(GBPlane.normal, gb->n);
    VECTOR_COPY(GBPlane.point, gb->p);
    FindIntersectionWithGB(home, node, node1, GBPlane, &line);
    VECTOR_COPY(pene->line1.direction, line.direction);
    VECTOR_COPY(pene->line1.point, line.point);

    FindIntersectionWithGB(home, node, node2, GBPlane, &line);
    VECTOR_COPY(pene->line2.direction, line.direction);
    VECTOR_COPY(pene->line2.point, line.point);

    centorPos[0] = (pene->line1.point[0] + pene->line2.point[0])*0.5;
    centorPos[1] = (pene->line1.point[1] + pene->line2.point[1])*0.5;
    centorPos[2] = (pene->line1.point[2] + pene->line2.point[2])*0.5;

    type = LineLineIntersection(pene->line1.direction, pene->line1.point,
           pene->line2.direction, pene->line2.point, vec1, vec2, &dis);
    if(type != LINE_INTERSECTION && fabs(dis) > 5){
        printf("Warning: type %d in file %s at line %d, dis %f \n", 
                type, __FILE__, __LINE__, dis);
    }

	pos[0] = node->x;
	pos[1] = node->y;
	pos[2] = node->z;

    NodeLocation(home, pos[0], pos[1], pos[2],currInGrain, currOnGB);
    TBID = gb->GBID;

    vec[0] = node->burgX[1] + node->burgX[0];
    vec[1] = node->burgY[1] + node->burgY[0];
    vec[2] = node->burgZ[1] + node->burgZ[0];
    if(Normal(vec) > 0.01){
        PrintNode(node);
        Fatal("In file %s at line %d", __FILE__, __LINE__);
    }

	burg[0] = node->burgX[1];
	burg[1] = node->burgY[1];
	burg[2] = node->burgZ[1];
	
	normal[0] = node->nx[1];
	normal[1] = node->ny[1];
	normal[2] = node->nz[1];
	
	dot = DotProduct(burg, GBPlane->normal);
 	if(fabs(dot) <1.0e-2)isScrew = 1;
	
	Matrix33Vector3Multiply(home->param->rotateNormal, normal, outNormal);
    outNormal[0] = -outNormal[0];
    outNormal[1] = -outNormal[1];
    outNormal[2] = -outNormal[2];

	Matrix33Vector3Multiply(oldGrain->orientationInv, burg, crystalBurg);

	if(isScrew){
		pene->nSets = 3;
        pene->res = (Res_t *)calloc(pene->nSets, sizeof(Res_t));
/*
 *      set 1, secrew dislocation cut through the CTB like cross-slip
 */ 
        i = 0;
        j = 0;
        k = 0;
        InitPeneRes(&(pene->res[i]), 1, 1);
        pene->res[i].nDissOnPlane[j] = 1;

        pene->res[i].passGB[j] = 1;
        pene->res[i].inGrain[j] = currInGrain;
		VECTOR_COPY(pene->res[i].normal[j],outNormal);
		VECTOR_COPY(pene->res[i].burg[k],burg);

        cross(pene->line2.direction, pene->res[i].normal[j], vec1);
        if(Normal(vec1) < 1.0E-4){
            Fatal("In file %s at line %d", __FILE__, __LINE__);
        }
        NormalizeVec(vec1);
       
        if(pene->res[i].passGB[j]){ 
            type = PointPlaneIntersection(pos, GBPlane->normal, 
                   GBPlane->point, vec2, &dis);
            vec2[0] -= pos[0];
            vec2[1] -= pos[1];
            vec2[2] -= pos[2];
            if(Normal(vec2) < 1.0E-4){
                printf("Warning in file %s at line %d", __FILE__, __LINE__);
            }
            NormalizeVec(vec2);
            if(DotProduct(vec1, vec2) < 0){
                vec1[0] = -vec1[0];
                vec1[1] = -vec1[1];
                vec1[2] = -vec1[2];
            }
        }

        pene->res[i].pos[k][0] = centorPos[0] + vec1[0]*eps*((real8)k+1.0);
        pene->res[i].pos[k][1] = centorPos[1] + vec1[1]*eps*((real8)k+1.0);
        pene->res[i].pos[k][2] = centorPos[2] + vec1[2]*eps*((real8)k+1.0);
        
/*
 *      set 2 and 3, dislocation dissociate at the CTB
 */
        i = 1;
        j = 0;
        for(i=1; i<3; i++){
            disID = 0;
            InitPeneRes(&(pene->res[i]), 1, 2);
            for(j=0; j<pene->res[i].nPlanes; j++){
                pene->res[i].nDissOnPlane[j] = 2;
                pene->res[i].passGB[j] = 0;
                pene->res[i].inGrain[j] = TBInGrain;
		        VECTOR_COPY(pene->res[i].normal[j], GBPlane->normal);
     
                cross(pene->line2.direction, pene->res[i].normal[j], vec1);
                if(Normal(vec1) < 1.0E-4){
                    Fatal("In file %s at line %d", __FILE__, __LINE__);
                }
                NormalizeVec(vec1);
             
                if(pene->res[i].passGB[j]){ 
                    type = PointPlaneIntersection(pos, GBPlane->normal, 
                           GBPlane->point, vec2, &dis);
                    vec2[0] -= pos[0];
                    vec2[1] -= pos[1];
                    vec2[2] -= pos[2];
                    if(Normal(vec2) < 1.0E-4){
                        printf("Warning in file %s at line %d", __FILE__, __LINE__);
                    }
                    NormalizeVec(vec2);
                    if(DotProduct(vec1, vec2) < 0){
                        vec1[0] = -vec1[0];
                        vec1[1] = -vec1[1];
                        vec1[2] = -vec1[2];
                    }
                }
                for(k = disID; k<disID+pene->res[i].nDissOnPlane[j]; k++){
                    pene->res[i].pos[k][0] = centorPos[0] 
                                             + vec1[0]*eps*((real8)(k-disID)+1.0);
                    pene->res[i].pos[k][1] = centorPos[1] 
                                             + vec1[1]*eps*((real8)(k-disID)+1.0);
                    pene->res[i].pos[k][2] = centorPos[2] 
                                             + vec1[2]*eps*((real8)(k-disID)+1.0);
                }
                disID += pene->res[i].nDissOnPlane[j];
            }
        }


/*
 *      Find the out burg
 * */
		switch(BurgID_FCC(gb)){
		    case 1 :
		    burg1[0] = -1.0;
		    burg1[1] =  2.0;
		    burg1[2] = -1.0;
		        break;
		    case 3 :
		    burg1[0] =  2.0;
		    burg1[1] = -1.0;
		    burg1[2] = -1.0;
		        break;
		    case 5 :
		    burg1[0] = -1.0;
		    burg1[1] =  2.0;
		    burg1[2] = -1.0;
		        break;
		    default :
			Fatal("GetPenRes: Wrong screw ID");
		}

		NormalizeVec(burg1);
		burg1[0] /= sqrt(3.0);	
		burg1[1] /= sqrt(3.0);	
		burg1[2] /= sqrt(3.0);	
		
	    Matrix33Vector3Multiply(oldGrain->orientation, burg1, burg2);

		if(DotProduct(burg2, burg) < 0.0){
		    burg2[0] = -burg2[0];
		    burg2[1] = -burg2[1];
		    burg2[2] = -burg2[2];
		}

		burg1[0] = burg[0] - burg2[0];
		burg1[1] = burg[1] - burg2[1];
		burg1[2] = burg[2] - burg2[2];

        VECTOR_COPY(pene->res[1].burg[0], burg1);
        VECTOR_COPY(pene->res[1].burg[1], burg2);

        VECTOR_COPY(pene->res[2].burg[0], burg2);
        VECTOR_COPY(pene->res[2].burg[1], burg1);
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
        
        if(crystalBurg[0] < 1.0e-6 &&
           crystalBurg[1] < 1.0e-6 &&
           crystalBurg[2] < 1.0e-6){
           sb = -1.0;
       }else sb = 1.0;

	    Matrix33Vector3Multiply(oldGrain->orientationInv, normal, crystalNormal);
        normalID = GlidePlaneID_FCC(crystalNormal);
        GetTBBurg_FCC(normalID, TBBurgCryst);
        TBBurgCryst[0] *= sb;
        TBBurgCryst[1] *= sb;
        TBBurgCryst[2] *= sb;
	    Matrix33Vector3Multiply(oldGrain->orientation, TBBurgCryst, burg2);

	    burg1[0] = burg[0] - burg2[0];
	    burg1[1] = burg[1] - burg2[1];
	    burg1[2] = burg[2] - burg2[2];
	
        real8 burg1Mag;
        burg1Mag = Normal(burg1);
        if(fabs(burg1Mag - 1.0) > 1.0e-3){
#if 0
            FormatVector(burg,"burg");
            FormatVector(burg1,"burg1");
            FormatVector(burg2,"burg2");
            Fatal("GetPenResults: burg1Mag is %f",burg1Mag);
#else
	    return(0);
#endif
	    }

        i = 0;
	    pene->nSets = 1;
        pene->res = (Res_t *)calloc(pene->nSets, sizeof(Res_t));

        InitPeneRes(&(pene->res[i]), 2, 2);

        pene->res[i].nDissOnPlane[0] = 1;
        pene->res[i].nDissOnPlane[1] = 1;

        pene->res[i].passGB[0] = 1;
        pene->res[i].passGB[1] = 0;

        pene->res[i].inGrain[0] = currInGrain;      
        pene->res[i].inGrain[1] = TBInGrain;      

        VECTOR_COPY(pene->res[i].normal[0],outNormal);
        VECTOR_COPY(pene->res[i].normal[1],GBPlane->normal);

	    VECTOR_COPY(pene->res[i].burg[0],burg1);
	    VECTOR_COPY(pene->res[i].burg[1],burg2);
    
        disID = 0;
        for(j=0; j<pene->res[i].nPlanes; j++){
            cross(pene->line2.direction, pene->res[i].normal[j], vec1);
            if(Normal(vec1) < 1.0E-4){
                Fatal("In file %s at line %d", __FILE__, __LINE__);
            }
            NormalizeVec(vec1);
             
            if(pene->res[i].passGB[j]){ 
                type = PointPlaneIntersection(pos, GBPlane->normal, 
                       GBPlane->point, vec2, &dis);
                vec2[0] -= pos[0];
                vec2[1] -= pos[1];
                vec2[2] -= pos[2];
                if(Normal(vec2) < 1.0E-4){
                    printf("Warning in file %s at line %d", __FILE__, __LINE__);
                }
                NormalizeVec(vec2);
                if(DotProduct(vec1, vec2) < 0){
                    vec1[0] = -vec1[0];
                    vec1[1] = -vec1[1];
                    vec1[2] = -vec1[2];
                }
            }
            for(k = disID; k<disID+pene->res[i].nDissOnPlane[j]; k++){
                pene->res[i].pos[k][0] = centorPos[0] 
                                         + vec1[0]*eps*((real8)(k-disID)+1.0);
                pene->res[i].pos[k][1] = centorPos[1] 
                                         + vec1[1]*eps*((real8)(k-disID)+1.0);
                pene->res[i].pos[k][2] = centorPos[2] 
                                         + vec1[2]*eps*((real8)(k-disID)+1.0);
            }
            disID += pene->res[i].nDissOnPlane[j];
        }

	    return(1);
}


int PenetrationFirst(Home_t *home, Node_t *node, Grain_t *oldGrain,
					 Grain_t *newGrain, Plane_t *GBPlane, Pene_t *pene){

	real8		nDir1[3], nodeVel[3], currPos[3], oldPos[3], burg[3], bMag;
	real8 		burgStore[3], normalStore[3],eps, sb, newVelocity, CurrNbrInGrain;
	real8		velNormal, glideDir[3], burg1[3],burg2[3], newVel[3], oldVel[3],crystalNormal[3];
	real8		newPos[3], vec[3], vecNormal, crossDistance, TBInGrain,currInGrain, oldInGrain;
	real8		node1Pos[3],node2Pos[3], posClosedGB[3], dis;
	real8		maxPowerDis = - 1.0e50, totVel[3], totForce[3], dirGBNode[3];
	int			armID, switchNode=0, set, maxPowerSet=-1;
	int			iScrew = 0, globalOp, mobError, forward, myDomian;
	int		    numBkup, disID;
    int         i, j, k, type, breakLoop;
	Node_t 		*nbrNode1,*nbrNode2,*newNode,*globalNode, *nbr;
	Node_t		*node1,*node2,*bkupNodeList, **newNodePtrList;
	Node_t		*node4,*node5;
	Param_t		*param;
	Plane_t		oldGlidePlane, outPlane, plane;							//Basic plane
	Line_t		pathLine, interLine, lineNodeOldNbr2, lineNbr1Node, estimatedPath;
	
	myDomian = home->myDomain;
	
	param = home->param;

	nodeVel[0] = node->vX;
	nodeVel[1] = node->vY;
	nodeVel[2] = node->vZ;
	
	currPos[0] = node->x;
	currPos[1] = node->y;
	currPos[2] = node->z;

    type = PointPlaneIntersection(currPos, GBPlane->normal,
           GBPlane->point, posClosedGB, &dis);
    if(fabs(dis) < 1E-10)return(0);
    
    dirGBNode[0] = currPos[0] - posClosedGB[0];
    dirGBNode[1] = currPos[1] - posClosedGB[1];
    dirGBNode[2] = currPos[2] - posClosedGB[2];
    NormalizeVec(dirGBNode);

    oldPos[0] = node->oldx;
    oldPos[1] = node->oldy;
    oldPos[2] = node->oldz;

    currInGrain = InWhichGrain(home, currPos[0], currPos[1], currPos[2]); 
    oldInGrain = node->inGrain;
/*	
 *	find the dislocation glide direction	
 */
	if(Normal(nodeVel) < 1.0)return(0);

	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	
	if(nbrNode1 == (Node_t *)NULL || nbrNode2 == (Node_t *)NULL)return(0);
	
    numBkup = 5 + nbrNode1->numNbrs +nbrNode2->numNbrs;
	bkupNodeList = (Node_t *)calloc(numBkup, sizeof(Node_t));
	BackupNode(home, nbrNode1, &bkupNodeList[0]);
	BackupNode(home, node, &bkupNodeList[1]);
	BackupNode(home, nbrNode2, &bkupNodeList[2]);

    for(i=0; i<nbrNode1->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode1,i);
	    BackupNode(home, nbr, &bkupNodeList[i+3]);
    }

    for(i=0; i<nbrNode2->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode2,i);
	    BackupNode(home, nbr, &bkupNodeList[i+3+nbrNode1->numNbrs]);
    }

	node1 = (Node_t *)NULL; 
	node2 = (Node_t *)NULL; 
	
// node2 pos
    VECTOR_COPY(node2Pos, pene->line2.point);
    FoldBox(param, node2Pos, node2Pos+1, node2Pos+2);

//node1pos
    VECTOR_COPY(node1Pos, pene->line1.point);
    FoldBox(param, node1Pos, node1Pos+1, node1Pos+2);

    TBInGrain = NodeLocation(home, node1Pos[0], node1Pos[1], node1Pos[2]);

/*	
 *	Insert a new native node and select two segements to deal with
 */
	globalOp = 0;
    
	vec[0] = nbrNode2->oldx - node->oldx;
	vec[1] = nbrNode2->oldy - node->oldy;
	vec[2] = nbrNode2->oldz - node->oldz;
    ZImage(param, vec, vec+1, vec+2);

    newPos[0] = node->oldx + 0.5*vec[0];
    newPos[1] = node->oldy + 0.5*vec[1];
    newPos[2] = node->oldz + 0.5*vec[2];
    FoldBox(param, newPos, newPos+1, newPos+2);
    
	if((globalNode = InsertNode(home,node,nbrNode2,newPos,globalOp,OPCLASS_SEPARATION)) == (Node_t *)NULL) {
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
	BackupNode(home, node1, &bkupNodeList[numBkup-2]);
	BackupNode(home, node2, &bkupNodeList[numBkup-1]);
	
/*
 *  Test all possible resultant configurations.
 */
    newNodePtrList = (Node_t **)NULL;

    breakLoop = 0;
	for(i=0; i < pene->nSets; i++){
        if(pene->res[i].effective == 0)continue;
		pene->res[i].powerDis = 0.0;
        newNodePtrList = (Node_t **)calloc(pene->res[i].nDiss, sizeof(Node_t*));

	    ChangeArmBurg(home, node1, &node2->myTag, 0.0, 0.0, 0.0,
				      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
	    ChangeArmBurg(home, node2, &node1->myTag, 0.0, 0.0, 0.0,
				      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);

        node1->inGrain = TBInGrain;
        node2->inGrain = TBInGrain;

        disID = 0;
        for(j=0; j<pene->res[i].nPlanes; j++){
            for(k=disID; k<disID+pene->res[i].nDissOnPlane[j]; k++){
                newNodePtrList[k] = (Node_t *)NULL;
                InsertConnection(home, node1, node2, globalOp, pene->res[i].inGrain[j],
                                pene->res[i].pos[k], pene->res[i].burg[k],
                                pene->res[i].normal[j], &(newNodePtrList[k]));
            } // for k
            disID += pene->res[i].nDissOnPlane[j];
        } // for j
        if(disID != pene->res[i].nDiss){
            Fatal("in file %s at line %d", __FILE__, __LINE__);
        }
/*
 *      caculate the power dis
 */
        breakLoop = 0;
        disID = 0;
        for(j=0; j<pene->res[i].nPlanes; j++){
            if(breakLoop)break;
            VECTOR_ZERO(totVel);
            VECTOR_ZERO(totForce);
            for(k=disID; k<disID+pene->res[i].nDissOnPlane[j]; k++){
		        SetOneNodeForce(home,newNodePtrList[k]);
		        mobError = EvaluateMobility(home,newNodePtrList[k]);
		        totVel[0] += newNodePtrList[k]->vX; 
		        totVel[1] += newNodePtrList[k]->vY; 
		        totVel[2] += newNodePtrList[k]->vZ;

                totForce[0] += newNodePtrList[k]->fX;
                totForce[1] += newNodePtrList[k]->fY;
                totForce[2] += newNodePtrList[k]->fZ;

                if(pene->res[i].passGB[j]){
                    VECTOR_ZERO(newVel);
                    newVel[0] = newNodePtrList[k]->vX;
                    newVel[1] = newNodePtrList[k]->vY;
                    newVel[2] = newNodePtrList[k]->vZ;
                    NormalizeVec(newVel);

                    if(DotProduct(newVel, dirGBNode) <= 0){
                        pene->res[i].effective = 0;
                        breakLoop = 1;
                        break;
                    }
                } 
            } // for k
            disID += pene->res[i].nDissOnPlane[j];
            pene->res[i].powerDis += DotProduct(totVel, totForce);
        } // for j
        
/*
 *      Remove the local new nodes
 */
        for(k=0; k<pene->res[i].nDiss; k++){
            ForceRemoveNode(home, newNodePtrList[k], globalOp);
        }
        free(newNodePtrList);
        newNodePtrList = (Node_t **)NULL;

/*
 *      Restore node1-node2 configuration
 */
        RestoreNode(home, node1, &bkupNodeList[numBkup-2]);
	    RestoreNode(home, node2, &bkupNodeList[numBkup-1]);

    } // for i


/*
 *  Go back to original configuration
 */
	ForceRemoveNode(home, globalNode, globalOp);
    globalNode = (Node_t *)NULL;
	node->inGrain = oldInGrain;

	RestoreNode(home, nbrNode1, &bkupNodeList[0]);
	RestoreNode(home, node, &bkupNodeList[1]);
	RestoreNode(home, nbrNode2, &bkupNodeList[2]);

    for(i=0; i<nbrNode1->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode1,i);
	    RestoreNode(home, nbr, &bkupNodeList[i+3]);
    }
    for(i=0; i<nbrNode2->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode2,i);
	    RestoreNode(home, nbr, &bkupNodeList[i+3+nbrNode1->numNbrs]);
    }

	for (i = 0; i < numBkup; i++) {
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

/*
 *  find maximal power dissipation
 */
	maxPowerSet = -1;
    maxPowerDis = -1.0E50;
	for(set=0; set < pene->nSets; set++){
		if(pene->res[set].effective == 0)continue;
        maxPowerDis = MAX(maxPowerDis,pene->res[set].powerDis);
		if(maxPowerDis == pene->res[set].powerDis){
	        maxPowerSet = set;
		}
	}
	if(maxPowerDis < 0 || maxPowerSet < 0){
		printf("Penetrationfirst(%d,%d) failed. numSet %d: maxPowerDis is %e\n",node->myTag.domainID, 
               node->myTag.index, pene->nSets, maxPowerDis);
		return(0);
	}
	
/*
 *  create the two segements again with globalOp = 1
 */
	globalOp = 1;
	
	vec[0] = nbrNode2->oldx - node->oldx;
	vec[1] = nbrNode2->oldy - node->oldy;
	vec[2] = nbrNode2->oldz - node->oldz;
    ZImage(param, vec, vec+1, vec+2);

    newPos[0] = node->oldx + 0.5*vec[0];
    newPos[1] = node->oldy + 0.5*vec[1];
    newPos[2] = node->oldz + 0.5*vec[2];
    FoldBox(param, newPos, newPos+1, newPos+2);

	if((globalNode = InsertNode(home,node,nbrNode2,newPos,globalOp,OPCLASS_SEPARATION)) == (Node_t *)NULL) {
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
	
/*
 *  use the topological change of maximal power dissipation
 */

	ChangeArmBurg(home, node1, &node2->myTag, 0.0, 0.0, 0.0,
	              0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
	ChangeArmBurg(home, node2, &node1->myTag, 0.0, 0.0, 0.0,
				  0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);

	AssignGrain(home, node1, TBInGrain);
	AssignGrain(home, node2, TBInGrain);

    disID = 0;
    i = maxPowerSet;
    for(j=0; j<pene->res[i].nPlanes; j++){
        for(k=disID; k<disID + pene->res[i].nDissOnPlane[j]; k++){
            InsertConnection(home, node1, node2, globalOp, pene->res[i].inGrain[j],
                            pene->res[i].pos[k], pene->res[i].burg[k],
                            pene->res[i].normal[j], &(newNode));

	        MarkNodeForceObsolete(home, newNode);        
	        MarkNodeNoCollisions(home, newNode);        
	        MarkNodeNoMeshCoarsen(home, newNode);        
        } // for k
        disID += pene->res[i].nDissOnPlane[j];
    } // for j

	MarkNodeForceObsolete(home, node1);        
	MarkNodeNoCollisions(home, node1);        
	MarkNodeNoMeshCoarsen(home, node1);        

	MarkNodeForceObsolete(home, node2);        
	MarkNodeNoCollisions(home, node2);        
	MarkNodeNoMeshCoarsen(home, node2);        

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
	
	int		    i;
	real8		currNbrInGrain;
	real8		currInGrain, newPos[3], currNbrPos[3], dis;
	real8 		nDir1[3], burg[3], nodeVel[3], oldPos[3], currPos[3];
    real8       vec1[3], vec2[3], dir[3], vec[3], oldNbrPos[3];
	Line_t		interLine;
    int         mergeStatus, mergeDone;

	Node_t		*nbrNode1, *nbrNode2, *mergedNode;
    Param_t     *param;

	Plane_t		oldGlidePlane, plane;

    param = home->param;
		
	currPos[0] = node->x;
	currPos[1] = node->y;
	currPos[2] = node->z;
	currInGrain = InWhichGrain(home, currPos[0], currPos[1], currPos[2]);
	
	nDir1[0] = node->nx[0];
	nDir1[1] = node->ny[0];
	nDir1[2] = node->nz[0];
	
    if(PlanePlaneIntersection(nDir1, currPos, GBPlane->normal, 
       GBPlane->point, vec1, vec2) != LINE_INTERSECTION)return(0);
    VECTOR_COPY(interLine.direction, vec1);
    VECTOR_COPY(interLine.point, vec2);

	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	
	if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL) return(0);
	
/*	nbrNode2 is always intersection node.	*/
	if(NodeOnGB(home,nbrNode1)){
		Node_t *nodeTemp;
		nodeTemp = nbrNode1;
		nbrNode1 = nbrNode2;
		nbrNode2 = nodeTemp;
	}
    if(SurfaceNode(nbrNode1))return(0);
	
    currNbrPos[0] = nbrNode1->x;
    currNbrPos[1] = nbrNode1->y;
    currNbrPos[2] = nbrNode1->z;
    FoldBox(param, currNbrPos, currNbrPos+1, currNbrPos+2);
	currNbrInGrain = InWhichGrain(home,currNbrPos[0],currNbrPos[1],currNbrPos[2]);

    oldNbrPos[0] = nbrNode1->oldx;
    oldNbrPos[1] = nbrNode1->oldy;
    oldNbrPos[2] = nbrNode1->oldz;
	
	if(currNbrInGrain == nbrNode1->inGrain || currNbrInGrain == -1){
	    vec[0] = currNbrPos[0] - currPos[0];
	    vec[1] = currNbrPos[1] - currPos[1];
	    vec[2] = currNbrPos[2] - currPos[2];
	}else{
	    vec[0] = oldNbrPos[0] - currPos[0];
	    vec[1] = oldNbrPos[1] - currPos[1];
	    vec[2] = oldNbrPos[2] - currPos[2];
	}
	
    ZImage(param, &vec[0], &vec[1], &vec[2]);
    VECTOR_COPY(dir, vec);
    NormalizeVec(dir);
    switch (LineLineIntersection(dir, currPos, 
       interLine.direction, interLine.point, vec1, newPos, &dis)){
        case NO_INTERSECTION:
            printf("PenetrationSide: warning at %d\n",__LINE__);
            return(0);
            break;
        case LINE_INTERSECTION:
            newPos[0] = currPos[0] + 0.5*vec[0];
            newPos[1] = currPos[1] + 0.5*vec[1];
            newPos[2] = currPos[2] + 0.5*vec[2];
            if(PointLineIntersection(newPos, interLine.direction, 
               interLine.point, newPos, &dis) == NO_INTERSECTION){
                printf("PenetrationSide: warning at %d\n",__LINE__);
            }
            break;
        case POINT_INTERSECTION:
            break;
    }
	FoldBox(param,&newPos[0],&newPos[1],&newPos[2]);	

    currNbrPos[0] = nbrNode2->x;
    currNbrPos[1] = nbrNode2->y;
    currNbrPos[2] = nbrNode2->z;
    FoldBox(param, currNbrPos, currNbrPos+1, currNbrPos+2);

    vec[0] = newPos[0] - currNbrPos[0]; 
    vec[1] = newPos[1] - currNbrPos[1]; 
    vec[2] = newPos[2] - currNbrPos[2]; 
    ZImage(param, vec, vec+1, vec+2);

    switch(AdjustNodeDirection(home, nbrNode2, vec)){
        case 3:
            VECTOR_COPY(newPos, currNbrPos);
            break;
        case 2:
            newPos[0] = currNbrPos[0] + vec[0]; 
            newPos[1] = currNbrPos[1] + vec[1]; 
            newPos[2] = currNbrPos[2] + vec[2];
            break;
        default:
            printf("PenetrationSide: warning at %d\n",__LINE__);
    }
    FoldBox(param,newPos, newPos+1, newPos+2);
    MergeNode(home, OPCLASS_SEPARATION, node, nbrNode2, newPos,
              &mergedNode, &mergeStatus, 1);
    mergeDone = mergeStatus & MERGE_SUCCESS;
	
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
	    vec[0] = node->x - node->oldx;
	    vec[1] = node->y - node->oldy;
	    vec[2] = node->z - node->oldz;
	    ZImage(home->param, &vec[0], &vec[1], &vec[2]);
	    if(Normal(vec) < 1e-6)return(0);

	    NormalizeVec(vec);
	    Line_t line;
	    VECTOR_COPY(line.direction, vec);
	    line.point[0] = node->oldx;
	    line.point[1] = node->oldy;
	    line.point[2] = node->oldz;
	    
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
	    mergePos[0] = nbrNode1->x;
	    mergePos[1] = nbrNode1->y;
	    mergePos[2] = nbrNode1->z;
	}
	if(nbrNode2->numNbrs > 3){
	    mergePos[0] = nbrNode2->x;
	    mergePos[1] = nbrNode2->y;
	    mergePos[2] = nbrNode2->z;
	}

	PullBack(home,node, node->inGrain);
	FoldBox(home->param,&mergePos[0],&mergePos[1],&mergePos[2]);

	MergeNode(home, OPCLASS_SEPARATION, nbrNode1,nbrNode2,
			  mergePos, &mergedNode, &mergeStatus, globalOp);
	
	if((mergeStatus & MERGE_SUCCESS) == 0){
	    return(0);
	}
    
    if(mergedNode != (Node_t *)NULL){
        if(mergedNode->numNbrs == 2){
	        nbrNode1 = GetNeighborNode(home,mergedNode,0);
	        nbrNode2 = GetNeighborNode(home,mergedNode,1);
            if(nbrNode1 == (Node_t *)NULL || nbrNode2 == (Node_t *)NULL){
                Fatal("in file %s at line %d",__FILE__, __LINE__);
            }

            if(nbrNode1->inGrain == nbrNode2->inGrain && 
               nbrNode1->inGrain - floor(nbrNode1->inGrain) == 0){
	            AssignGrain(home, mergedNode, nbrNode1->inGrain);
            }
        }
    }

	return(1);
}

int PenetrationSurface(Home_t *home, Node_t *node, GB_t *gb)
{
    int         type;
    real8       n1[3], p1[3], n2[3], p2[3], n3[3], p[3];
    real8       lDir[3], lPos[3], vec[3], newPos[3];
    Node_t      *nbrNode;
	Param_t 	*param;

	param = home->param;

    if(!SurfaceNode(node))Fatal("in %s at %d", __FILE__, __LINE__);

	nbrNode = GetNeighborNode(home, node, 0);
    if(nbrNode == (Node_t *)NULL)Fatal("in %s at %d", __FILE__, __LINE__);
    if(nbrNode->myTag.domainID != home->myDomain || !IsOnGB(nbrNode))return(0);
    if(!Surface(node))Fatal("in %s at %d", __FILE__, __LINE__);

    VECTOR_ZERO(n1);
    n1[node->constraint - 1] = 1.0;
    p1[0] = node->x;
    p1[1] = node->y;
    p1[2] = node->z;

    VECTOR_COPY(n2, gb->n);
    p2[0] = nbrNode->x;
    p2[1] = nbrNode->y;
    p2[2] = nbrNode->z;

    type = PlanePlaneIntersection(n2, p2, n1, p1, lDir, lPos);
    if(type != LINE_INTERSECTION){
        printf("Warning in PenetrationSurface: in %s at %d", __FILE__, __LINE__);
        return(0);
    }

    n3[0] = node->nx[0];
    n3[1] = node->ny[0];
    n3[2] = node->nz[0];
    VECTOR_COPY(p3, p2);
    type = LinePlaneIntersection(lDir, lPos, n3, p3, vec, newPos);
    if(type != POINT_INTERSECTION){
        printf("Warning in PenetrationSurface: in %s at %d", __FILE__, __LINE__);
        return(0);
    }
    

    FoldBox(home->param, newPos, newPos+1, newPos+2);
	RepositionNode(home, newPos, &nbrNode->myTag, 1);
	AssignConstraint(home, nbrNode, node->constraint);

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

    return(1);
}


void HandleGrainBoundary(Home_t *home)
{
	
	int 		i,j,k,m,numNbrOnGB, globalOp = 1, mergeStatus;
	int	 	    grainNode = 0, endNode, newConstraint;
    int         newNodeKeyPtr, type;
	int 		oldInGrain, oldOnGB, currInGrain, currOnGB;	
    int         crossSlipNode, nbr1OnGB, nbr2OnGB;
	real8		x,y,z,burgj[3], burgk[3], vec[3], mergPos[3];
	real8		burg[3], normal[3], nDir1[3], nDir2[3];
	real8		mergePos[3], distance, newInGrain, newPos[3];
    real8       vec1[3], vec2[3];
	Param_t		*param;
	Node_t		*node, *nbrNode1, *nbrNode2, *nbr1, *nbr2, *mergedNode;
	Plane_t		GBPlane, norPlane, surPlane;
	Line_t		line;
	Grain_t		*oldGrain, *currGrain;
    GB_t        *gb;
	Pene_t	    pene;

	param = home->param;
	
    InitPenStr(&pene);
	newNodeKeyPtr = home->newNodeKeyPtr;
	for(i=0;i<newNodeKeyPtr;i++){
        
        if((node = home->nodeKeys[i]) == (Node_t *) NULL)continue;

        oldInGrain = node->inGrain;

        x = node->x; y = node->y; z=node->z;
        FoldBox(param, &x, &y, &z);
        NodeLocation(home, x, y, z, currInGrain, currOnGB);
            	
        if(IsOnGB(node)){
            if(currOnGB < 0){
                PullBack(home, node);
            }
            continue;
        }
		
        if(node->numNbrs < 1 node->constraint == PINNED_NODE)continue;
        if(currInGrain == oldInGrain || currInGrain < 0)continue;

        if(!param->peneTB || node->numNbrs > 2){
            PullBack(home, node, oldInGrain);
            continue;
        }
		
/*	
 *	    Due to different slip planes of two segments of crossSlip node,
 *	    we only allow such node penetrateLast.  	
 */		
        crossSlipNode = 0;
        if(node->numNbrs == 2){
	        if(fabs(node->nx[0] * node->nx[1] +
	                node->ny[0] * node->ny[1] +
                    node->nz[0] * node->nz[1]) < 0.995){
                crossSlipNode = 1;
            }
        }		

        oldGrain  = home->grainKeys[oldInGrain];
        currGrain = home->grainKeys[corInGrain];
        if(!GetPublicGB(oldGrain, currGrain, &gb)){
            Fatal("in %s at %d", __FILE__, __LINE__);
        }

/*	
 *	    There are 4 scenarios: 
 *	  	PenetrationFirst: the first penetrating node.
 *		PenetrationSide: one side penetrating nodes.
 *		PenetrationLast: the last node.
 *		PenetrationSurface: for the surface node.
 */
		
        if(SurfaceNode(node)){
            if(node->numNbrs != 1){
                PullBack(home, node, oldInGrain);
                continue;
            }
			
            nbrNode1 = GetNeighborNode(home,node,0);
            if(nbrNode1 == (Node_t *)NULL){
                Fatal("HandleGrainBoundary: nbrNode1 not found");
            }

	        if(!IsOnGB(nbrNode1) || nbrNode1->myTag.domainID != home->myDomain){
	            PullBack(home, node);
	            continue;
	        }

	        if(!PenetrationSurface(home, node, gb)){
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

/*
 *      If the neighbor is not only intersect GB node, it is also junction node,
 *      pull back the penetrated node.
 */
        numNbrOnGB = 0;
        nbr1OnGB = NodeOnGB(home, nbrNode1);
        nbr2OnGB = NodeOnGB(home, nbrNode2);

	    if(nbr1OnGB > 0) {
	        if(!(nbr1OnGB & TB_NODE) && nbrNode1->numNbrs == 3){
	           PullBack(home, node, node->inGrain);
	           continue;
	        } 
	        numNbrOnGB++;	
	    }

	    if(nbr2OnGB > 0) {
	        if(!(nbr2OnGB & TB_NODE) && nbrNode2->numNbrs == 3){
	           PullBack(home, node, node->inGrain);
	           continue;
	        } 
	        numNbrOnGB++;	
	    }
			
	    if(numNbrOnGB != 2 && crossSlipNode == 1){
/*
 *		CrossSlip node is only allowed to penetrationLast. 
 */
            PullBack(home,node, oldInGrain);
            continue;
        }
            	
	    if(numNbrOnGB == 0 && crossSlipNode != 1){
            if(!GetPenResults(home, node, &pene, gb)){
		        PullBack(home,node, oldInGrain);
		        continue;
		    }
        }
/*		
 *	        Node: this energy rule is used for nonscrew dislocation, 
 *	        since we donot know the resultant slip plane for nonscrew 
 *	        dislocation at this time. 
 */

#if ENERGY_RULER
	        if(EnergyRule(home, node, &GBPlane, &pene) == 0){
                PullBack(home,node, oldInGrain);
                continue;
	        }
#endif

            burg[0] = node->burgX[0];
            burg[1] = node->burgY[0];
            burg[2] = node->burgZ[0];
            normal[0] = GBPlane.normal[0];
            normal[1] = GBPlane.normal[1];
            normal[2] = GBPlane.normal[2];
		
            if((fabs(DotProduct(burg,normal))<0.001) &&
				(numNbrOnGB > 0)){
/*		
 *	        Energy rule for screw dislocation with at least a nbr on GB, so we can
 *		    determine the resultant slip plane.	
 */
                if(ScrewEnergyRule(home,node,oldInGrain,
			   corInGrain,&GBPlane)==0){
	
                    PullBack(home,node, oldInGrain);
                    continue;
                }
            }

            if(numNbrOnGB == 0){

                if(!PenetrationFirst(home, node, oldGrain, currGrain, &GBPlane, &pene)){
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
	
    FreePenStr(&pene);
	if(param->peneTB == 0) return;

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
	   
        if((NodeOnGB(home, node) & INTERSECT_NODE) == 0)continue;

	    vec[0] = node->x;
	    vec[1] = node->y;
	    vec[2] = node->z;

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

        if(node->numNbrs == 2){
            nbr1 = GetNeighborNode(home,node,0);
            nbr2 = GetNeighborNode(home,node,1);
            if(nbr1==(Node_t *)NULL||nbr2==(Node_t *)NULL) continue;
            if(!(nbr1->inGrain - floor(nbr1->inGrain) == 0 &&
               nbr2->inGrain == nbr1->inGrain)) continue;
        }
			
        for(j=0;j<node->numNbrs;j++){
            nbr1 = GetNeighborNode(home,node,j);
            if(nbr1 == (Node_t *)NULL) continue;
            if(nbr1->inGrain - floor(nbr1->inGrain) == 0.5)continue;
            if(nbr1->myTag.domainID != node->myTag.domainID) continue;
                	
            for(k=j+1;k<node->numNbrs;k++){
                nbr2 = GetNeighborNode(home,node,k);
                if(nbr2 == (Node_t *)NULL) continue;
	
		        if(nbr1->inGrain != nbr2->inGrain)continue;
                if(nbr1 == nbr2) continue;
                if(nbr2->inGrain - floor(nbr2->inGrain) == 0.5)continue;
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
					
	            globalOp = 1;
                ChangeConnection(home,nbr1,&node->myTag,
                                 &nbr2->myTag,globalOp);
                ChangeConnection(home,nbr2,&node->myTag,
                                 &nbr1->myTag,globalOp);
	
                ChangeArmBurg(home, node , &nbr1->myTag,
                              0,0,0,0,0,0,globalOp, DEL_SEG_NONE);
                ChangeArmBurg(home, node , &nbr2->myTag,
                             0,0,0,0,0,0,globalOp, DEL_SEG_NONE);
					
                if(RemoveDoubleLinks(home,nbr1,globalOp) == 2){
                    snprintf(name,50,"RDL%d.3",home->cycle);
                    Tecplot(home, name, 0, 0, 0 ,0, 0, 1);
                    printf("Unusual RemoveDoubleLinks , 3\n");
                }
			
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

#if 0
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
	        vec[0] = nbr2->x; 
	        vec[1] = nbr2->y; 
	        vec[2] = nbr2->z; 
	        newInGrain = nbr2->inGrain; 
		    newConstraint = nbr1->constraint;
	        norPlane.normal[0] = node->nx[1];
	        norPlane.normal[1] = node->ny[1];
	        norPlane.normal[2] = node->nz[1];
		
	        surPlane.point[0] = nbr1->x;
	        surPlane.point[1] = nbr1->y;
	        surPlane.point[2] = nbr1->z;
	    }else{
		    if(NodeDistance(home, node, nbr2) > param->minSeg)continue;
	        vec[0] = nbr1->x; 
	        vec[1] = nbr1->y; 
	        vec[2] = nbr1->z; 
	        newInGrain = nbr1->inGrain; 
		    newConstraint = nbr2->constraint;
	        norPlane.normal[0] = node->nx[0];
	        norPlane.normal[1] = node->ny[0];
	        norPlane.normal[2] = node->nz[0];
	        surPlane.point[0] = nbr2->x;
	        surPlane.point[1] = nbr2->y;
	        surPlane.point[2] = nbr2->z;
	    }
	    VECTOR_ZERO(GBPlane.normal);
	    distance =  ClosedGB(home, vec, &GBPlane);

	    GBPlane.point[0] += (Sign(distance) * param->rann * GBPlane.normal[0]);
	    GBPlane.point[1] += (Sign(distance) * param->rann * GBPlane.normal[1]);
	    GBPlane.point[2] += (Sign(distance) * param->rann * GBPlane.normal[2]);

	    norPlane.point[0] = node->x;
	    norPlane.point[1] = node->y;
	    norPlane.point[2] = node->z;

        type = PlanePlaneIntersection(norPlane.normal, norPlane.point,
               GBPlane.normal, GBPlane.point, vec1, vec2);

	    VECTOR_ZERO(surPlane.normal);
	    surPlane.normal[newConstraint - 1] = 1.0;

        switch (type){
            case LINE_INTERSECTION:
                if(LinePlaneIntersection(vec1, vec2, surPlane.normal,
                   surPlane.point, vec, newPos) != LINE_INTERSECTION){
                    printf("HandleGrounBoundary: warning at Line %d \n", __LINE__);
                    continue;
                }
                break;
            case PLANE_INTERSECTION:
                if(PlanePlaneIntersection(vec1, vec2, surPlane.normal,
                   surPlane.point, vec, newPos) == NO_INTERSECTION){
                    printf("HandleGrounBoundary: warning at Line %d \n", __LINE__);
                    continue;
                }
                FormatVector(vec1, "vec1");
                FormatVector(vec2, "vec2");
                FormatVector(surPlane.normal, "surface normal");
                FormatVector(surPlane.point, "surface point");
                
                printf("HandleGrounBoundary: warning at Line %d, type %d\n", __LINE__,type);
                break;
            defualt:
                printf("HandleGrounBoundary: warning at Line %d, type %d\n", __LINE__,type);
                break;
        }

        FoldBox(param, newPos, newPos+1, newPos+2);
        
	    RepositionNode(home,newPos,&node->myTag,1);
	    AssignGrain(home, node, newInGrain);
	    AssignConstraint(home, node, newConstraint);
	}
#endif
	for(i=0;i<home->newNodeKeyPtr;i++){
            if((node = home->nodeKeys[i]) == (Node_t *)NULL)continue;
            if(node->inGrain - floor(node->inGrain) == 0.5)continue;
	        if(SurfaceNode(node))continue;
            if(node->inGrain != InWhichGrain(home, node->x, node->y, node->z)){
                PrintTag(node);
	            PullBack(home,node, node->inGrain);
                printf("HandleGrainBoundary: a penetrated not coopered\n");
            }
        }
	return;

}

void DoResistance(Home_t *home, Node_t *node, real8 *nForce){
	real8 	distance, pos[3],dir[3], factor, burgMag;
	real8	GBForceMag, lDir[3], nDir[3], gDir[3], length;
	real8 	burg[3], fricMag, dot, armGBForce[3];
	real8	rann, armf[3], forceGB[3], segCentor[3];
	real8	oriNForce[3], gForce[3], burgNor[3], norForce[3];
	int	    i, cacuGBForce, fricForce, armInGrain, armOnGB;
    GB_t    *gb;
	Node_t	*nbr;
	Param_t *param;
	
	param = home->param;
	pos[0] = node->x;
	pos[1] = node->y;
	pos[2] = node->z;
	VECTOR_ZERO(forceGB);

	nForce[0] = node->fX;
	nForce[1] = node->fY;
	nForce[2] = node->fZ;

    if(param->polycrystal){
	    rann = 3.0 * param->rann;

	    for(i=0; i<node->numNbrs; i++){
	        nbr = GetNeighborNode(home,node,i);
	        if(nbr == (Node_t *)NULL)continue;
    
            GetArmGrain(home, node, nbr, armInGrain, armOnGB);
	        lDir[0] = nbr->x - node->x;
	        lDir[1] = nbr->y - node->y;
	        lDir[2] = nbr->z - node->z;
	        ZImage(param, lDir, lDir+1, lDir+2);
    
	        segCentor[0] = node->x + 0.5 * lDir[0];
	        segCentor[1] = node->y + 0.5 * lDir[1];
	        segCentor[2] = node->z + 0.5 * lDir[2];
    
	        length = Normal(lDir);
    
	        if(length < 1e-8)continue;
    
	        NormalizeVec(lDir);
	        nDir[0] = node->nx[i];
	        nDir[1] = node->ny[i];
	        nDir[2] = node->nz[i];
	        cross(lDir, nDir, gDir);
    
	        if(Normal(gDir) < 1e-5) continue;
	        NormalizeVec(gDir);
    
	        burg[0] = node->burgX[i];
	        burg[1] = node->burgY[i];
	        burg[2] = node->burgZ[i];
    
	        armf[0] = node->armfx[i];
	        armf[1] = node->armfy[i];
	        armf[2] = node->armfz[i];
    
	        dot = DotProduct(armf, gDir);
	        gForce[0] = dot * gDir[0];
	        gForce[1] = dot * gDir[1];
	        gForce[2] = dot * gDir[2];
    
            norForce[0] = armf[0] - gForce[0];
            norForce[1] = armf[1] - gForce[1];
            norForce[2] = armf[2] - gForce[2];
    
	        burgMag = sqrt(node->burgX[i] * node->burgX[i] +
			            node->burgY[i] * node->burgY[i] +
			            node->burgZ[i] * node->burgZ[i]);
    
	        distance =  ClosedGB(home, segCentor, &gb);
	        cacuGBForce &= (fabs(distance) < 2.0e0 * rann);
	        cacuGBForce &= (fabs(distance) > 1.0e-3);
	        cacuGBForce &= (armOnGB == -1);
    
/*
 *	        screw dislocation can easily cut though GB.
 */
            VECTOR_COPY(burgNor, burg);
            NormalizeVec(burgNor);
	        dot = fabs(DotProduct(gb->n, burgNor));
	        cacuGBForce &= (dot > 1.0e-1); // mon-screw dislocation
    
	        if(cacuGBForce == 1){
    
	            dir[0] = Sign(distance) * gb->n[0];
	            dir[1] = Sign(distance) * gb->n[1];
	            dir[2] = Sign(distance) * gb->n[2];
                if(DotProduct(gDir,dir) < 0){
                    gDir[0] = -gDir[0];
                    gDir[1] = -gDir[2];
                    gDir[2] = -gDir[3];
                }
    
	            GBForceMag = param->GBStrength * pow(2.7182,-(distance*distance/rann/rann));
	            GBForceMag *= (length * burgMag * 0.5e0);
                armGBForce[0] = GBForceMag * gDir[0];
                armGBForce[1] = GBForceMag * gDir[1];
                armGBForce[2] = GBForceMag * gDir[2];
    
                nForce[0] += armGBForce[0]; 
                nForce[1] += armGBForce[1]; 
                nForce[2] += armGBForce[2]; 
    
                armf[0] += armGBForce[0]; 
                armf[1] += armGBForce[1]; 
                armf[2] += armGBForce[2]; 
    
                gForce[0] += armGBForce[0]; 
                gForce[1] += armGBForce[1]; 
                gForce[2] += armGBForce[2]; 
	        }
        }
    }
    
	VECTOR_COPY(oriNForce, nForce);

	for(i=0; i<node->numNbrs; i++){
	    nbr = GetNeighborNode(home,node,i);
	    if(nbr == (Node_t *)NULL)continue;

	    lDir[0] = nbr->x - node->x;
	    lDir[1] = nbr->y - node->y;
	    lDir[2] = nbr->z - node->z;
	    ZImage(param, lDir, lDir+1, lDir+2);
	    length = Normal(lDir);

	    if(length < 1e-8)continue;

	    NormalizeVec(lDir);
	    nDir[0] = node->nx[i];
	    nDir[1] = node->ny[i];
	    nDir[2] = node->nz[i];
	    cross(lDir, nDir, gDir);

	    if(Normal(gDir) < 1e-5) continue;
	    NormalizeVec(gDir);

	    burg[0] = node->burgX[i];
	    burg[1] = node->burgY[i];
	    burg[2] = node->burgZ[i];

	    armf[0] = node->armfx[i];
	    armf[1] = node->armfy[i];
	    armf[2] = node->armfz[i];

	    dot = DotProduct(armf, gDir);
	    gForce[0] = dot * gDir[0];
	    gForce[1] = dot * gDir[1];
	    gForce[2] = dot * gDir[2];

	    burgMag = sqrt(node->burgX[i] * node->burgX[i] +
			           node->burgY[i] * node->burgY[i] +
			           node->burgZ[i] * node->burgZ[i]);

/*
 *      Friction Foreces
 */
	    fricMag = 0.5e0 * length * burgMag * param->fricStress;

	    if(Normal(gForce) <= fricMag){
		    nForce[0] -= gForce[0];
		    nForce[1] -= gForce[1];
		    nForce[2] -= gForce[2];

            armf[0] -= gForce[0];
            armf[1] -= gForce[1];
            armf[2] -= gForce[2];
	    }else{
	        NormalizeVec(gForce);
	        nForce[0] -= (fricMag * gForce[0]);
	        nForce[1] -= (fricMag * gForce[1]);
	        nForce[2] -= (fricMag * gForce[2]);

	        armf[0] -= (fricMag * gForce[0]);
	        armf[1] -= (fricMag * gForce[1]);
	        armf[2] -= (fricMag * gForce[2]);
        }

    }

    if(DotProduct(nForce,oriNForce) <= 0.0){
        VECTOR_ZERO(nForce);
    }

	return;

}

void ApplyPosConstraint(Param_t *param, Node_t *node){

	real8 disVec[3], maxDis= -1E50;
	int   i, constraintSurface;
	real8 pos[3];

	pos[0] = node->x;
	pos[1] = node->y;
	pos[2] = node->z;

	VECTOR_ZERO(disVec);

	if(param->xBoundType == Free &&  fabs(pos[0])>param->xBoundMax) { 
	    disVec[0] = fabs(pos[0]) - param->xBoundMax;
	}
	if(param->yBoundType == Free &&  fabs(pos[1])>param->yBoundMax) { 
	    disVec[1] = fabs(pos[1]) - param->yBoundMax;
	}
	if(param->zBoundType == Free &&  fabs(pos[2])>param->zBoundMax) { 
	    disVec[2] = fabs(pos[2]) - param->zBoundMax;
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


void AdjustSurfaceNode(Home_t *home)
{

	int 	i, j, numNodePtr, thisDomain, splitOK, switchNode;
    int     inGrain, onGB;
	real8	nodePos[3], nbrPos[3], nbrConstraint;
	real8	interPos[3], oldPos[3];
	real8	newPos[3], burg0[3], burgj[3], vec[3];
 	Param_t *param;
	Node_t  *node, *nbr, *newNode, *node1, *node2;
	int 	nodecst, surfaceNbr, constraint;
	Line_t	line;
	Plane_t	plane;
	param = home->param;
  
	for (i = 0; i < home->newNodeKeyPtr; i++){

        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
	    nodePos[0] = node->x;
	    nodePos[1] = node->y;
	    nodePos[2] = node->z;
	    ApplyPosConstraint(param, node);

	    if(SurfaceNode(node) && node->numNbrs > 1){
            NodeLocation(node, nodePos[0], nodePos[1], nodePos[2], inGrain, onGB);
	        if(inGrain != node->inGrain && !IsOnGB(node)){
		        node->x = node->oldx;
		        node->y = node->oldy;
		        node->z = node->oldz;

	            nodePos[0] = node->x;
	            nodePos[1] = node->y;
	            nodePos[2] = node->z;
	            ApplyPosConstraint(param, node);
	        }
	    }

	    if(SurfaceNode(node) || node->numNbrs < 2) continue;
#if 1
	    surfaceNbr = 0;
	    for(j=0;j<node->numNbrs;j++){
		    if((nbr=GetNeighborNode(home,node,j))==(Node_t *)NULL) continue;
		    if(SurfaceNode(nbr)){
                surfaceNbr++;
                nbrConstraint = nbr->constraint;
            }
		}
        if(surfaceNbr == node->numNbrs){
            AssignConstraint(home, node, nbrConstraint);
        }
#endif
	}

	node = home->ghostNodeQ;

    while (node) {
	    nodePos[0] = node->x;
	    nodePos[1] = node->y;
	    nodePos[2] = node->z;
	    ApplyPosConstraint(param, node);

	    if(SurfaceNode(node) && node->numNbrs > 1){
            NodeLocation(node, nodePos[0], nodePos[1], nodePos[2], inGrain, onGB);
	        if(inGrain != node->inGrain && !IsOnGB(node)){
		        node->x = node->oldx;
		        node->y = node->oldy;
		        node->z = node->oldz;

	            nodePos[0] = node->x;
	            nodePos[1] = node->y;
	            nodePos[2] = node->z;
	            ApplyPosConstraint(param, node);
	        }
	    }

        node = node->next;
	}

    numNodePtr = home->newNodeKeyPtr;
	for (i = 0; i<numNodePtr; i++){
        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
	    if(SurfaceNode(node))continue;
	    nodePos[0] = node->x;
	    nodePos[1] = node->y;
	    nodePos[2] = node->z;

	    if(param->xBoundType == Free &&  fabs(nodePos[0]) >= param->xBoundMax  ||
	       param->yBoundType == Free &&  fabs(nodePos[1]) >= param->yBoundMax  ||
	       param->zBoundType == Free &&  fabs(nodePos[2]) >= param->zBoundMax)continue;

	    for(j=0; j < node->numNbrs; j++){
	        nbr = GetNeighborNode(home, node, j);
	        if(nbr == (Node_t *)NULL)continue;
		    if(nbr->numNbrs < 2)continue;
	        if(!SurfaceNode(nbr))continue;
	        nbrPos[0] = nbr->x;
	        nbrPos[1] = nbr->y;
	        nbrPos[2] = nbr->z;

	        if(param->xBoundType == Free &&  fabs(nbrPos[0]) <= param->xBoundMax + 1.0E-10 &&
	           param->yBoundType == Free &&  fabs(nbrPos[1]) <= param->yBoundMax + 1.0E-10 &&
	           param->zBoundType == Free &&  fabs(nbrPos[2]) <= param->zBoundMax + 1.0E-10)continue;

		    nodecst = node->constraint;
	        constraint = nbr->constraint;
		    VECTOR_COPY(line.point, nodePos);
	        line.direction[0] = nbrPos[0] - nodePos[0];
	        line.direction[1] = nbrPos[1] - nodePos[1];
	        line.direction[2] = nbrPos[2] - nodePos[2];

	        if(Normal(line.direction) < 1.0e-3){
	            continue;
		    }

	        NormalizeVec(line.direction);		
		
		    VECTOR_ZERO(plane.normal);
		    plane.normal[nbr->constraint - 1] = 1.0;

	        VECTOR_ZERO(plane.point);

	        if(nbr->constraint == 1){
	            plane.point[0] = Sign(nbrPos[0]) * param->xBoundMax;
	        }else if(nbr->constraint == 2){
	            plane.point[1] = Sign(nbrPos[1]) * param->yBoundMax;
	        }else if(nbr->constraint == 3){
	            plane.point[2] = Sign(nbrPos[2]) * param->zBoundMax;
	        }

	        VECTOR_ZERO(interPos);

            if(LinePlaneIntersection(line.direction, line.point, plane.normal,
               plane.point, vec, interPos) != POINT_INTERSECTION){
                Fatal("AdjustSurfaceNode %d:", __LINE__);
            }

		    if(fabs(interPos[0] - nbrPos[0]) < 1e-10 &&
		       fabs(interPos[1] - nbrPos[1]) < 1e-10 &&
		       fabs(interPos[2] - nbrPos[2]) < 1e-10){
                printf("Warning: file %s at %d, too closed\n", __FILE__, __LINE__);
		        interPos[0] = 0.99999 * nbrPos[0] + 0.00001 * nodePos[0]; 
		        interPos[1] = 0.99999 * nbrPos[1] + 0.00001 * nodePos[1]; 
		        interPos[2] = 0.99999 * nbrPos[2] + 0.00001 * nodePos[2]; 
		    }

            NodeLocation(node, interPos[0], interPos[1], interPos[2], inGrain, onGB);
            if(inGrain != node->inGrain && !IsOnGB(node) ||
               onGB != node->inGrain && IsOnGB(node)){
		        VECTOR_COPY(interPos, nbrPos);
	        }

            splitOK = 0;
            if(DomainOwnsSeg(home, OPCLASS_REMESH,
               thisDomain, &nbr->myTag))splitOK = 1;     
            if(splitOK){
                newNode = InsertNode(home, node, nbr, interPos, 1, OPCLASS_REMESH);
                if(newNode == (Node_t *)NULL){
#if 1
                    PrintNode(node);
                    PrintNode(nbr);
                    printf("AdjustSurfaceNode: Insert failed\n");
#endif
                }
	            switchNode = 0;
	            if(GetArmID(home, node, nbr) != -1)switchNode=1;
	
	            if(switchNode){
                    node1 = newNode;
                    node2 = node;
	            }else{
	                node1 = node;
	                node2 = newNode;
	            }
                AssignConstraint(home, node1, nodecst);
                AssignConstraint(home, node2, constraint);
            }

	    } /*for(j=0; j< ...)*/
	}/*for(i=0; i< ...)*/	

	return;
}

