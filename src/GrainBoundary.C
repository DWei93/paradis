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

#include <vector>
#include <algorithm>
#include <numeric>

#define ENERGY_RULER 1

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
            int     type;
            GB_t    *gb;
            real8   nodePos[3], dis;
            
            nodePos[0] = node->x;
            nodePos[1] = node->y;
            nodePos[2] = node->z;
            gb = home->GBKeys[node->onGB];

            type = PointPlaneIntersection(nodePos, gb->n, gb->p, newPos, &dis);
            if(fabs(dis) > 5){
                PrintNode(node);
                FormatVector(newPos, "newPos");
                printf("Warning in PullBack, error is %f", dis);
            }
        }
    }else{
        if(currOnGB == -1){
            if(node->inGrain != currInGrain){
                printf("newPos--{%f,%f,%f}--(%d,%d)\n", newPos[0],
                       newPos[1],newPos[2], currInGrain, currOnGB);
                PrintNode(node);
                Tecplot(home, (char*)"NewPosToGB", 0, 0, 0 ,0, 0, 1);
                printf("Warning in %s at %d", __FILE__, __LINE__);

                if(!NewPosToGB(home, node, node->inGrain, currOnGB, newPos)){
                    PrintNode(node);
                    Fatal("in %s at %d", __FILE__, __LINE__);
                }
            }
        }else{
            if(!NewPosToGB(home, node, node->inGrain, currOnGB, newPos)){
                PrintNode(node);
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
        }
    }

	FoldBox(home->param, newPos, newPos+1, newPos+2);
    RepositionNode(home,newPos,&node->myTag,1);

	return;
}

Node_t *InsertNode(Home_t *home, Node_t *node1, Node_t *node2, real8 *p, int globalOp, int opClass)
{

    int     armID, splitStatus, armIndex;
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

    armIndex = GetArmID(home, node1, node2);

    splitStatus = SplitNode(home, opClass,
                     node1, nodePos, p,
                     nodeVel, nodeVel,
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

real8 FGlidePerLength(Home_t *home, Node_t *node){

    int	    i;
    real8	disLength;
    real8   vec[3];
    real8	fGlide;
    Node_t  *nbr;
    
    if(node == (Node_t *)NULL){
        printf("Task %d: %s, line %d -- Node = NULL",
               home->myDomain, __FILE__, __LINE__);
        return(0);
    }

    disLength = 0;

    for(i=0;i<node->numNbrs;i++){
       nbr = GetNeighborNode(home,node,i);
       if(nbr==(Node_t *)NULL) continue;
       disLength += NodeDistance(home,node,nbr)/2;
    }

    vec[0] = node->vX;
    vec[1] = node->vY;
    vec[2] = node->vZ;

    if(Normal(vec)<1E-2)return(0.0);
    NormalizeVec(vec);

    fGlide = node->fX*vec[0]+
         node->fY*vec[1]+
         node->fZ*vec[2];

    return(fabs(fGlide)/disLength);
}

real8 SegementEnergy(Param_t *param, real8 cosTheta, real8 bMag){
	static real8 pi = 3.1415926535897932385;
	
	return(param->shearModulus * bMag * bMag * log(param->outRadius) *
		   (1.0 - param->pois * cosTheta * cosTheta)/
		   (4.0 * pi * (1.0 - param->pois)) );
		  
}

int EnergyRule(Home_t *home, Node_t *node, Plane_t *GBPlane){
	
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


int ScrewEnergyRule(Home_t *home, Node_t *node, GB_t *gb){

	return(1);
}


void FreePenStr(Pene_t *pene){
    if(pene->nSets > 0){
        int i, j;
        for(i=0; i<pene->nSets; i++){
            if(pene->res[i].nDiss > 0){
                free(pene->res[i].burg);
                free(pene->res[i].pos);
                free(pene->res[i].normal);
                free(pene->res[i].inGrain);


                pene->res[i].burg = NULL;
                pene->res[i].pos = NULL;
                pene->res[i].normal = NULL;
                pene->res[i].inGrain = NULL;
            }

            pene->res[i].nDiss = 0;
            pene->res[i].effective = 0;
            pene->res[i].powerDis = 0;
        }
        free(pene->res);
    }

    pene->nSets = 0;

    pene->res = (Res_t *)NULL;
    pene->gb = (GB_t *)NULL;
    pene->node = (Node_t *)NULL;

    VECTOR_ZERO(pene->line1.direction);
    VECTOR_ZERO(pene->line2.direction);
    VECTOR_ZERO(pene->line1.point);
    VECTOR_ZERO(pene->line2.point);
    return;
}



void InitPeneRes(Res_t *res, int nDiss){

    if(res->nDiss > 0){
        free(res->pos);
        free(res->burg);
        free(res->inGrain);        
        free(res->normal);        
        
        res->pos = NULL;
        res->burg = NULL;
        res->inGrain = NULL;
        res->normal = NULL;
    }
    res->nDiss = nDiss;
    
    res->effective = 1;
    res->powerDis = -1.0E50;

    res->inGrain = (int *)malloc(nDiss*sizeof(int));
    res->normal = (real8 (*)[3])malloc(nDiss*sizeof(real8)*3);

    res->pos = (real8 (*)[3])malloc(nDiss*sizeof(real8)*3);
    res->burg = (real8 (*)[3])malloc(nDiss*sizeof(real8)*3);
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
        vec1[0] = (n1Pos[0] + n2Pos[0] - n1OldPos[0] - n2OldPos[0]);
        vec1[1] = (n1Pos[1] + n2Pos[1] - n1OldPos[1] - n2OldPos[1]);
        vec1[2] = (n1Pos[2] + n2Pos[2] - n1OldPos[2] - n2OldPos[2]);
        VECTOR_COPY(vec2, centor);
    }else{
        vec1[0] = n1Pos[0] - n2Pos[0];
        vec1[1] = n1Pos[1] - n2Pos[1];
        vec1[2] = n1Pos[2] - n2Pos[2];

        vec2[0] = (n1Pos[0] + n2Pos[0])*0.5;
        vec2[1] = (n1Pos[1] + n2Pos[1])*0.5;
        vec2[2] = (n1Pos[2] + n2Pos[2])*0.5;
    }
    if(Normal(vec1) < 1.0E-4){
        PrintNode(node1);
        PrintNode(node2);
        Fatal("In file %s at line %d", __FILE__, __LINE__);
    }
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

int GetPenResults(Home_t *home, Node_t *node, Pene_t *pene, GB_t *gb)
{
	
	real8		pos[3], burg[3], normal[3], resiBurg[3];
	real8		vec[3], crossProd[3], jth;
    real8       centorPos[3];
    real8       vec1[3], vec2[3], eps, dis;
    real8       barrierFactor;
    int         i, j, k, newPlaneID, localGBID, nDiss;
    int         currInGrain, currOnGB, nSets, type;
    bool        residualDis;
    Line_t      line;
    Node_t      *node1, *node2;
    Plane_t     GBPlane;
    Grain_t     *grain;

    BurgInfo_t  *gbd;
    
    std::vector<std::vector<double> >                       burgList;
    std::vector<std::vector< Dislocation_t> >               disArray; 
    std::vector< std::vector<Dislocation_t> >::iterator     itDisList;

    pene->node = node;
    pene->gb = gb;
    gbd = &(gb->burgData);
    eps = home->param->rann;

	node1 = GetNeighborNode(home,node,0);
	node2 = GetNeighborNode(home,node,1);
    if(node1 == (Node_t *)NULL || node2 == (Node_t *)NULL){
        Fatal("In file %s at line %d", __FILE__, __LINE__);
    }

    VECTOR_COPY(GBPlane.normal, gb->n);
    VECTOR_COPY(GBPlane.point, gb->p);
    FindIntersectionWithGB(home, node, node1, &GBPlane, &line);
    VECTOR_COPY(pene->line1.direction, line.direction);
    VECTOR_COPY(pene->line1.point, line.point);

    FindIntersectionWithGB(home, node, node2, &GBPlane, &line);
    VECTOR_COPY(pene->line2.direction, line.direction);
    VECTOR_COPY(pene->line2.point, line.point);

    centorPos[0] = (pene->line1.point[0] + pene->line2.point[0])*0.5;
    centorPos[1] = (pene->line1.point[1] + pene->line2.point[1])*0.5;
    centorPos[2] = (pene->line1.point[2] + pene->line2.point[2])*0.5;

    type = LineLineIntersection(pene->line1.direction, pene->line1.point,
           pene->line2.direction, pene->line2.point, vec1, vec2, &dis);
    if(type != LINE_INTERSECTION && fabs(dis) > 5){
        printf("Warning: type %d in file %s at line %d, dis is %f \n", 
                type, __FILE__, __LINE__, dis);
        return(0);
    }
#if CROSS_SLIP_GB_MODEL
    return(0);
#endif  

	pos[0] = node->x;
	pos[1] = node->y;
	pos[2] = node->z;

    NodeLocation(home, pos[0], pos[1], pos[2],currInGrain, currOnGB);

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

    
    if(!GeometricCriterion(gbd, burg, pene->line2.direction, normal, 1,
                       burgList, disArray))return(0);
	
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
    if(disArray.size() == 0)return(0);

    pene->nSets = disArray.size();
    pene->res = (Res_t *)calloc(pene->nSets, sizeof(Res_t));

    for(i=0; i<disArray.size(); i++){

        VECTOR_COPY(resiBurg, burg);
        for(j=0; j<disArray[i].size(); j++){
            resiBurg[0] -= burgList[disArray[i][j].burgID][0];
            resiBurg[1] -= burgList[disArray[i][j].burgID][1];
            resiBurg[2] -= burgList[disArray[i][j].burgID][2];
        }

        nDiss = disArray[i].size();
        if(Normal(resiBurg) > 1.0E-2)nDiss++;
        residualDis = (nDiss == disArray[i].size()) ? 0 : 1;

        InitPeneRes(&(pene->res[i]), nDiss);

        newPlaneID = -1;
        for(j=0; j<disArray[i].size(); j++){
            VECTOR_COPY(pene->res[i].burg[j], burgList[disArray[i][j].burgID]);
            VECTOR_COPY(pene->res[i].normal[j], 
                        gbd->planeList[disArray[i][j].planeID]);
            pene->res[i].inGrain[j] = gbd->planeInGrain[disArray[i][j].planeID];

            if(disArray[i][j].planeID != newPlaneID){
                newPlaneID = disArray[i][j].planeID;
                jth = 1.0;
            }else{
                jth += 1.0; 
            }
            cross(pene->line2.direction, pene->res[i].normal[j], vec);
            if(Normal(vec) < 1.0E-4){
                Fatal("In file %s at line %d", __FILE__, __LINE__);
            }
            NormalizeVec(vec);

            if(pene->res[i].inGrain[j] > -1){
                grain = home->grainKeys[pene->res[i].inGrain[j]];
                localGBID = GetGBIndexOfGrain(grain, gb->GBID);
                if(localGBID == -1){
                    Fatal("In file %s at line %d", __FILE__, __LINE__);
                }
                vec1[0] = -grain->sign[localGBID]*gb->n[0];
                vec1[1] = -grain->sign[localGBID]*gb->n[1];
                vec1[2] = -grain->sign[localGBID]*gb->n[2];
                
                if(DotProduct(vec, vec1) < 0){
                    vec[0] = -vec[0];
                    vec[1] = -vec[1];
                    vec[2] = -vec[2];
                }
            }else{
                cross(gb->n, pene->res[i].normal[j], vec1);
                if(Normal(vec1) > 1.0E-4){
                    Fatal("In file %s at line %d", __FILE__, __LINE__);
                }
            }

            pene->res[i].pos[j][0] = centorPos[0] + vec[0]*eps*jth;
            pene->res[i].pos[j][1] = centorPos[1] + vec[1]*eps*jth;
            pene->res[i].pos[j][2] = centorPos[2] + vec[2]*eps*jth;
            
        } /* for(j=0; j<disArray[i].size(); j++) */

        if(residualDis){
            VECTOR_COPY(pene->res[i].pos[j], centorPos);
            VECTOR_COPY(pene->res[i].burg[j], resiBurg);

            cross(resiBurg, pene->line2.direction, pene->res[i].normal[j]);
            NormalizeVec(pene->res[i].normal[j]);
        
            pene->res[i].inGrain[j] = -1;

        }

    }/* for(i=0; i<disArray.size(); i++) */

	return(1);
}

int PenetrationFirst(Home_t *home, Pene_t *pene)
#if CROSS_SLIP_GB_MODEL
{

    int         myDomain, type, currInGrain, currOnGB, oldInGrain;
    int         numBkup, i, j, globalOp, mobError, newPlaneID;
    int         planeID, localPlaneID, maxPowerSet, onGB;
    int         newInGrain, newOnGB;
    bool        switchedArm, passGB;
    real8       nodeVel[3], currPos[3], oldPos[3], posClosedGB[3];
    real8       dis, dirGBNode[3], node1Pos[3], node2Pos[3];
    real8       vec1[3], vec2[3], newVel[3], nForce[3], norm;
    real8       centorPos[3], lDir[3], peneLength, planeBurg[3];
    real8       maxPowerDis;
    Node_t      *node, *nbrNode1, *nbrNode2, *nbr, *node1, *node2;
    Node_t      *localNode, *newNode;
    Node_t      *bkupNodeList, **newNodePtrList; 
    GB_t        *gb;
    BurgInfo_t  *bd;
    
    real8       gammaS, gammaT;

    Param_t     *param;
	
	myDomain = home->myDomain;
	param = home->param;
    
    gammaS = 4.47E-2/param->burgMag;  /* units: Pa.b*/
    gammaT = 0.5*gammaS;

    node = pene->node;
    gb = pene->gb;
    bd = &(gb->burgData);

    centorPos[0] = (pene->line1.point[0] + pene->line2.point[0])*0.5;
    centorPos[1] = (pene->line1.point[1] + pene->line2.point[1])*0.5;
    centorPos[2] = (pene->line1.point[2] + pene->line2.point[2])*0.5;

    vec1[0] = (pene->line1.point[0] - pene->line2.point[0]);
    vec1[1] = (pene->line1.point[1] - pene->line2.point[1]);
    vec1[2] = (pene->line1.point[2] - pene->line2.point[2]);
    peneLength = Normal(vec1);

	nodeVel[0] = node->vX;
	nodeVel[1] = node->vY;
	nodeVel[2] = node->vZ;
	
	currPos[0] = node->x;
	currPos[1] = node->y;
	currPos[2] = node->z;

    oldPos[0] = node->oldx;
    oldPos[1] = node->oldy;
    oldPos[2] = node->oldz;

    type = PointPlaneIntersection(currPos, gb->n,
                                  gb->p, posClosedGB, &dis);
    if(fabs(dis) < 1E-10)return(0);

    vec1[0] = currPos[0] - oldPos[0];
    vec1[0] = currPos[0] - oldPos[0];
    vec1[0] = currPos[0] - oldPos[0];
    ZImage(param, vec1, vec1+1, vec1+2);
    if(DotProduct(vec1, nodeVel) <= 0)return(0);

    NodeLocation(home, currPos[0], currPos[1], currPos[2], currInGrain, currOnGB);
    oldInGrain = node->inGrain;

	if(Normal(nodeVel) < 1.0)return(0);

	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	
	if(nbrNode1 == (Node_t *)NULL || nbrNode2 == (Node_t *)NULL)return(0);
	
	node1 = (Node_t *)NULL; 
	node2 = (Node_t *)NULL; 
	
//node1pos
    VECTOR_COPY(node1Pos, pene->line1.point);
    FoldBox(param, node1Pos, node1Pos+1, node1Pos+2);

// node2 pos
    VECTOR_COPY(node2Pos, pene->line2.point);
    FoldBox(param, node2Pos, node2Pos+1, node2Pos+2);

/*	
 *	Insert a new native node and select two segements to deal with
 */
	globalOp = 1;

	if((localNode = InsertNode(home,node,nbrNode2,oldPos,
        globalOp,OPCLASS_SEPARATION)) == (Node_t *)NULL) {
	    return(0);
	}

    switchedArm = (GetArmID(home, node, nbrNode1) == -1);
	
	if(switchedArm){
        node1 = localNode;
        node2 = node;
	}else{
	    node1 = node;
	    node2 = localNode;
	}
	
	RepositionNode(home, node2Pos, &node2->myTag, globalOp);
	RepositionNode(home, node1Pos, &node1->myTag, globalOp);

    AssignGrain(home, node1, -1, gb->GBID);
    AssignGrain(home, node2, -1, gb->GBID);
    return(1);
}
#else
{

    int         myDomain, type, currInGrain, currOnGB, oldInGrain;
    int         numBkup, i, j, globalOp, mobError, newPlaneID;
    int         planeID, localPlaneID, maxPowerSet, onGB;
    int         newInGrain, newOnGB;
    bool        switchedArm, passGB;
    real8       nodeVel[3], currPos[3], oldPos[3], posClosedGB[3];
    real8       dis, dirGBNode[3], node1Pos[3], node2Pos[3];
    real8       vec1[3], vec2[3], newVel[3], nForce[3], norm;
    real8       centorPos[3], lDir[3], peneLength, planeBurg[3];
    real8       maxPowerDis;
    Node_t      *node, *nbrNode1, *nbrNode2, *nbr, *node1, *node2;
    Node_t      *localNode, *newNode;
    Node_t      *bkupNodeList, **newNodePtrList; 
    GB_t        *gb;
    BurgInfo_t  *bd;
    
    real8       gammaS, gammaT;

    Param_t     *param;
	
	myDomain = home->myDomain;
	param = home->param;
    
    gammaS = 4.47E-2/param->burgMag;  /* units: Pa.b*/
    gammaT = 0.5*gammaS;

    node = pene->node;
    gb = pene->gb;
    bd = &(gb->burgData);

    centorPos[0] = (pene->line1.point[0] + pene->line2.point[0])*0.5;
    centorPos[1] = (pene->line1.point[1] + pene->line2.point[1])*0.5;
    centorPos[2] = (pene->line1.point[2] + pene->line2.point[2])*0.5;

    vec1[0] = (pene->line1.point[0] - pene->line2.point[0]);
    vec1[1] = (pene->line1.point[1] - pene->line2.point[1]);
    vec1[2] = (pene->line1.point[2] - pene->line2.point[2]);
    peneLength = Normal(vec1);

	nodeVel[0] = node->vX;
	nodeVel[1] = node->vY;
	nodeVel[2] = node->vZ;
	
	currPos[0] = node->x;
	currPos[1] = node->y;
	currPos[2] = node->z;

    oldPos[0] = node->oldx;
    oldPos[1] = node->oldy;
    oldPos[2] = node->oldz;

    type = PointPlaneIntersection(currPos, gb->n,
                                  gb->p, posClosedGB, &dis);
    if(fabs(dis) < 1E-10)return(0);

    vec1[0] = currPos[0] - oldPos[0];
    vec1[0] = currPos[0] - oldPos[0];
    vec1[0] = currPos[0] - oldPos[0];
    ZImage(param, vec1, vec1+1, vec1+2);
    if(DotProduct(vec1, nodeVel) <= 0)return(0);
    
    dirGBNode[0] = currPos[0] - posClosedGB[0];
    dirGBNode[1] = currPos[1] - posClosedGB[1];
    dirGBNode[2] = currPos[2] - posClosedGB[2];
    NormalizeVec(dirGBNode);

    NodeLocation(home, currPos[0], currPos[1], currPos[2], currInGrain, currOnGB);
    oldInGrain = node->inGrain;

	if(Normal(nodeVel) < 1.0)return(0);

	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	
	if(nbrNode1 == (Node_t *)NULL || nbrNode2 == (Node_t *)NULL)return(0);
	
    numBkup = 2 + nbrNode1->numNbrs +nbrNode2->numNbrs;
	bkupNodeList = (Node_t *)calloc(numBkup, sizeof(Node_t));
	BackupNode(home, nbrNode1, &bkupNodeList[0]);
	BackupNode(home, node, &bkupNodeList[1]);
	BackupNode(home, nbrNode2, &bkupNodeList[2]);

    j = 3;
    for(i=0; i<nbrNode1->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode1,i);
        if(nbr == node)continue;
	    BackupNode(home, nbr, &bkupNodeList[j]);
        j++;
    }

    if(j != 3 + nbrNode1->numNbrs - 1){
        Fatal("in %s at %d", __FILE__, __LINE__);
    }

    for(i=0; i<nbrNode2->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode2,i);
        if(nbr == node)continue;
	    BackupNode(home, nbr, &bkupNodeList[j]);
        j++;
    }

	node1 = (Node_t *)NULL; 
	node2 = (Node_t *)NULL; 
	
//node1pos
    VECTOR_COPY(node1Pos, pene->line1.point);
    FoldBox(param, node1Pos, node1Pos+1, node1Pos+2);

// node2 pos
    VECTOR_COPY(node2Pos, pene->line2.point);
    FoldBox(param, node2Pos, node2Pos+1, node2Pos+2);

/*	
 *	Insert a new native node and select two segements to deal with
 */
	globalOp = 0;

	if((localNode = InsertNode(home,node,nbrNode2,oldPos,
        globalOp,OPCLASS_SEPARATION)) == (Node_t *)NULL) {
	    return(0);
	}

    switchedArm = (GetArmID(home, node, nbrNode1) == -1);
	
	if(switchedArm){
        node1 = localNode;
        node2 = node;
	}else{
	    node1 = node;
	    node2 = localNode;
	}
	
	RepositionNode(home, node2Pos, &node2->myTag, globalOp);
	RepositionNode(home, node1Pos, &node1->myTag, globalOp);
	BackupNode(home, localNode, &bkupNodeList[numBkup-1]);
	
/*
 *  Test all possible resultant configurations.
 */
    newNodePtrList = (Node_t **)NULL;

	for(i=0; i < pene->nSets; i++){
        if(pene->res[i].effective == 0)continue;
		pene->res[i].powerDis = 0.0;
        newNodePtrList = (Node_t **)malloc(pene->res[i].nDiss * sizeof(Node_t*));

	    ChangeArmBurg(home, node1, &node2->myTag, 0.0, 0.0, 0.0,
				      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
	    ChangeArmBurg(home, node2, &node1->myTag, 0.0, 0.0, 0.0,
				      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);

        node1->inGrain = -1;
        node2->inGrain = -1;
        node1->onGB = gb->GBID;
        node2->onGB = gb->GBID;

/*
 *      caculate the power dis
 */
        newPlaneID = -1;
        for(j=0; j<pene->res[i].nDiss; j++){

            planeID = PlaneID_FCC(bd, pene->res[i].normal[j]);
            if(planeID != newPlaneID){
                newPlaneID = planeID;
                localPlaneID = 0;
            }else{
                localPlaneID++;
                if(localPlaneID == 2)localPlaneID = 0;
            }

            onGB = (pene->res[i].inGrain[j] == -1) ? gb->GBID : -1;
            InsertConnection(home, node1, node2, globalOp, pene->res[i].inGrain[j], onGB,
                             pene->res[i].pos[j], pene->res[i].burg[j],
                             pene->res[i].normal[j], &(newNodePtrList[j]));
            
		    SetOneNodeForce(home, newNodePtrList[j]);
            
#if 0
            if(pene->res[i].inGrain[j] > -1){
                vec1[0] = centorPos[0] - pene->res[i].pos[j][0];
                vec1[1] = centorPos[1] - pene->res[i].pos[j][1];
                vec1[2] = centorPos[2] - pene->res[i].pos[j][2];
                if(localPlaneID == 1){
                    vec1[0] = -vec1[0];
                    vec1[1] = -vec1[1];
                    vec1[2] = -vec1[2];
                }
                newNodePtrList[j]->fX += (vec1[0]*peneLength*gammaS);
                newNodePtrList[j]->fY += (vec1[1]*peneLength*gammaS);
                newNodePtrList[j]->fZ += (vec1[2]*peneLength*gammaS);

                newNodePtrList[j]->armfx[0] += 0.5*(vec1[0]*peneLength*gammaS);
                newNodePtrList[j]->armfy[0] += 0.5*(vec1[1]*peneLength*gammaS);
                newNodePtrList[j]->armfz[0] += 0.5*(vec1[2]*peneLength*gammaS);

                newNodePtrList[j]->armfx[1] += 0.5*(vec1[0]*peneLength*gammaS);
                newNodePtrList[j]->armfy[1] += 0.5*(vec1[1]*peneLength*gammaS);
                newNodePtrList[j]->armfz[1] += 0.5*(vec1[2]*peneLength*gammaS);
            }
#endif

		    mobError = EvaluateMobility(home, newNodePtrList[j]);
            newVel[0] = newNodePtrList[j]->vX;
            newVel[1] = newNodePtrList[j]->vY;
            newVel[2] = newNodePtrList[j]->vZ;

            if(pene->res[i].inGrain[j] > -1){
                VECTOR_COPY(vec2, newVel);
                norm = Normal(vec2);

                passGB = (norm <= 1.0E-2) ? 0 : 1; 
                if(norm > 1.0E-2){
                    vec2[0] /= norm;
                    vec2[1] /= norm;
                    vec2[2] /= norm;
                
                    vec1[0] = centorPos[0] + vec2[0];
                    vec1[1] = centorPos[1] + vec2[1];
                    vec1[2] = centorPos[2] + vec2[2];
                    NodeLocation(home, vec1[0], vec1[1], vec1[2], newInGrain, newOnGB); 
                
                    if(newOnGB > -1 || newInGrain != pene->res[i].inGrain[j])passGB = 0;
                }

                if(!passGB){
                    pene->res[i].effective = 0;
                    pene->res[i].powerDis = -1.0E50;
                    ForceRemoveNode(home, newNodePtrList[j], globalOp);
                    newNodePtrList[j] = (Node_t *)NULL;
                    break;
                }
            }

            DoResistance(home, newNodePtrList[j], nForce);

            pene->res[i].powerDis += DotProduct(nForce, newVel);

            ForceRemoveNode(home, newNodePtrList[j], globalOp);
            newNodePtrList[j] = (Node_t *)NULL;
            
/*
 *          Restore node1-node2 configuration
 */
	        RestoreNode(home, localNode, &bkupNodeList[numBkup-1]);
	        RestoreNode(home, node, &bkupNodeList[1]);

	        if(switchedArm){
                node1 = localNode;
                node2 = node;
	        }else{
	            node1 = node;
	            node2 = localNode;
	        }

            node1->inGrain = -1;
            node2->inGrain = -1;
            node1->onGB = gb->GBID;
            node2->onGB = gb->GBID;

	        RepositionNode(home, node2Pos, &node2->myTag, globalOp);
	        RepositionNode(home, node1Pos, &node1->myTag, globalOp);

        } // for j
        
/*
 *      Remove the local new nodes
 */
        free(newNodePtrList);
        newNodePtrList = (Node_t **)NULL;

    } // for i


/*
 *  Go back to original configuration
 */
	ForceRemoveNode(home, localNode, globalOp);
    localNode = (Node_t *)NULL;

	RestoreNode(home, nbrNode1, &bkupNodeList[0]);
	RestoreNode(home, node, &bkupNodeList[1]);
	RestoreNode(home, nbrNode2, &bkupNodeList[2]);

    j = 3;
    for(i=0; i<nbrNode1->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode1,i);
        if(nbr == node)continue;
	    RestoreNode(home, nbr, &bkupNodeList[j]);
        j++;
    }

    for(i=0; i<nbrNode2->numNbrs; i++){
        nbr = GetNeighborNode(home,nbrNode2,i);
        if(nbr == node)continue;
	    RestoreNode(home, nbr, &bkupNodeList[j]);
        j++;
    }

	for (i = 0; i < numBkup; i++) {
	    FreeNodeArrays(&bkupNodeList[i]);
	}
    
	free(bkupNodeList);
	bkupNodeList = (Node_t *)NULL;
	
/*
 *  find maximal power dissipation
 */
	maxPowerSet = -1;
    maxPowerDis = -1.0E50;
	for(i=0; i < pene->nSets; i++){
		if(pene->res[i].effective == 0)continue;
        maxPowerDis = MAX(maxPowerDis,pene->res[i].powerDis);
		if(maxPowerDis == pene->res[i].powerDis){
	        maxPowerSet = i;
		}
	}
	if(maxPowerDis < 0 || maxPowerSet < 0){
//		printf("Penetrationfirst(%d,%d) failed. numSet %d: maxPowerDis is %e\n",node->myTag.domainID, 
//               node->myTag.index, pene->nSets, maxPowerDis);
		return(0);
	}
	
/*
 *  create the two segements again with globalOp = 1
 */
    
    globalOp = 1;
	
	if((localNode = InsertNode(home,node,nbrNode2,oldPos,
        globalOp,OPCLASS_SEPARATION)) == (Node_t *)NULL) {
	    return(0);
	}

    switchedArm = (GetArmID(home, node, nbrNode1) == -1);
	
	if(switchedArm){
        node1 = localNode;
        node2 = node;
	}else{
	    node1 = node;
	    node2 = localNode;
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

	AssignGrain(home, node1, -1, gb->GBID);
	AssignGrain(home, node2, -1, gb->GBID);

    i = maxPowerSet;

    newPlaneID = -1;
    for(j=0; j<pene->res[i].nDiss; j++){

        planeID = PlaneID_FCC(bd, pene->res[i].normal[j]);
        if(planeID != newPlaneID){
            newPlaneID = planeID;
            localPlaneID = 0;
        }else{
            localPlaneID++;
            if(localPlaneID == 2)localPlaneID = 0;
        }

        if(pene->res[i].inGrain[j] > -1){
            if(localPlaneID == 0){
                VECTOR_COPY(planeBurg, pene->res[i].burg[j]);
                if(j != pene->res[i].nDiss-1){
                    if(planeID == PlaneID_FCC(bd, pene->res[i].normal[j+1]))continue;
                }
            }else{
                planeBurg[0] += pene->res[i].burg[j][0];
                planeBurg[1] += pene->res[i].burg[j][1];
                planeBurg[2] += pene->res[i].burg[j][2];
            }
        }else{
            VECTOR_COPY(planeBurg, pene->res[i].burg[j]);
        }

        onGB = (pene->res[i].inGrain[j] == -1) ? gb->GBID : -1;
        InsertConnection(home, node1, node2, globalOp, pene->res[i].inGrain[j], onGB,
                         pene->res[i].pos[j], planeBurg,
                         pene->res[i].normal[j], &newNode);

        MarkNodeForceObsolete(home, newNode);        
        MarkNodeNoCollisions(home, newNode);        
        MarkNodeNoMeshCoarsen(home, newNode);        
    }

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
    
#if 0
    printf("PenetrationFirst: \n");
    PrintNode(node);
    PrintNode(node1);
    PrintNode(node2);
#endif

	return(1);
}
#endif



int PenetrationSide(Home_t *home, Node_t *node, Grain_t *oldGrain, 
					Grain_t *currGrain, Plane_t *GBPlane)
#if CROSS_SLIP_GB_MODEL

{
	
	int		    i, nodeDeletable;
	int         currInGrain, currOnGB;
	int         currNbrInGrain, currNbrOnGB;
	real8		newPos[3], currNbrPos[3], dis;
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
    NodeLocation(home, currPos[0], currPos[1], currPos[2], currInGrain, currOnGB);
	
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
	if(NodeOnGB(home,nbrNode1) & INTERSECT_NODE){
		Node_t *nodeTemp;
		nodeTemp = nbrNode1;
		nbrNode1 = nbrNode2;
		nbrNode2 = nodeTemp;
	}

    FindIntersectionWithGB(home, node, nbrNode1, GBPlane, &interLine);
    VECTOR_COPY(newPos, interLine.point);

    FoldBox(param, newPos, newPos+1, newPos+2);
	RepositionNode(home, newPos, &node->myTag, 1);

    AssignGrain(home, node, -1, nbrNode2->onGB);
    return(1);
}
#else
{
	
	int		    i, nodeDeletable;
	int         currInGrain, currOnGB;
	int         currNbrInGrain, currNbrOnGB;
	real8		newPos[3], currNbrPos[3], dis;
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
    NodeLocation(home, currPos[0], currPos[1], currPos[2], currInGrain, currOnGB);
	
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
	if(NodeOnGB(home,nbrNode1) & INTERSECT_NODE){
		Node_t *nodeTemp;
		nodeTemp = nbrNode1;
		nbrNode1 = nbrNode2;
		nbrNode2 = nodeTemp;
	}
    if(SurfaceNode(nbrNode1))return(0);

	nodeDeletable  = (node->myTag.domainID == home->myDomain); 
	nodeDeletable &= (node->constraint == UNCONSTRAINED);	

    nodeDeletable &= ((node->flags & NO_COLLISIONS) == 0);
	nodeDeletable &= ((node->flags & NO_MESH_COARSEN) == 0);

	for (i = 0; i < node->numNbrs; i++) {
	    nodeDeletable &= DomainOwnsSeg(home, OPCLASS_SEPARATION, home->myDomain,
			                    &node->nbrTag[i]);
	}
    if(nodeDeletable == 0 && nbrNode2->myTag.domainID != home->myDomain)return(0);

    FindIntersectionWithGB(home, node, nbrNode1, GBPlane, &interLine);
    VECTOR_COPY(newPos, interLine.point);
	FoldBox(param,&newPos[0],&newPos[1],&newPos[2]);	
    
    if(BeforeMergeNodes(home, node, nbrNode2, newPos, &dis)){
        if(fabs(dis) > 5.0){
            printf("Warning in PenetrationSide, dis %f\n", dis);
        }
    }else{
        return(0);
    } 

    FoldBox(param,newPos, newPos+1, newPos+2);
    MergeNode(home, OPCLASS_SEPARATION, node, nbrNode2, newPos,
              &mergedNode, &mergeStatus, 1);
    mergeDone = mergeStatus & MERGE_SUCCESS;
    if(mergeDone == 0)return(0);

	return(1);
}

#endif

int PenetrationLast(Home_t *home, Node_t *node, Grain_t *oldGrain, 
					Grain_t *newGrain, Plane_t *GBPlane)
#if CROSS_SLIP_GB_MODEL

{

	real8 		burg[3], vec[3], mergePos[3],TBInGrain;
    real8       nodePos[3], oldPos[3], vec1[3], vec2[3];
    real8       dis;
	Node_t 		*nbrNode1, *nbrNode2, *mergedNode;
	int			mergeStatus,i,j, globalOp,q,thisDomain;
	int         crossSlipNode = 0, type;
    Param_t     *param;

    param = home->param;
	thisDomain = home->myDomain;
	
	burg[0] = node->burgX[0];
	burg[1] = node->burgY[0];
	burg[2] = node->burgZ[0];
	
	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL)
		return(0);
    if(nbrNode1->myTag.domainID != thisDomain || 
       nbrNode2->myTag.domainID != thisDomain)return(0);
	
	if(SurfaceNode(nbrNode1) && IsOnGB(nbrNode1) ||
	   SurfaceNode(nbrNode2) && IsOnGB(nbrNode2))return(0);
	globalOp = 1;

    nodePos[0] = node->x;
    nodePos[1] = node->y;
    nodePos[2] = node->z;
    FoldBox(param, nodePos, nodePos+1, nodePos+2);

    oldPos[0] = node->oldx;
    oldPos[1] = node->oldy;
    oldPos[2] = node->oldz;
    PBCPOSITION(param, nodePos[0], nodePos[1], nodePos[2], 
                       oldPos, oldPos+1, oldPos+2);
    
    vec[0] = nodePos[0] - oldPos[0];
    vec[1] = nodePos[1] - oldPos[1];
    vec[2] = nodePos[2] - oldPos[2];
    if(Normal(vec) < 1E-2)return(0);

    NormalizeVec(vec);
    type = LinePlaneIntersection(vec, oldPos, GBPlane->normal, GBPlane->point, 
                                 vec1, vec2);
    
    switch(type){
        case POINT_INTERSECTION:
            VECTOR_COPY(mergePos, vec2);
            break;
        case LINE_INTERSECTION:
            vec[0] = 0.5*(nodePos[0] + oldPos[0]);
            vec[1] = 0.5*(nodePos[1] + oldPos[1]);
            vec[2] = 0.5*(nodePos[2] + oldPos[2]);
            PointPlaneIntersection(vec, GBPlane->normal, GBPlane->point,
                                   mergePos, &dis);
            break;
        default:
            Fatal("PenetrationLast: wrong intersection type");

    }
    FoldBox(param, mergePos, mergePos+1, mergePos+2);
	RepositionNode(home, mergePos, &node->myTag, 1);

    AssignGrain(home, node, -1, nbrNode2->onGB);

    return(1);
}
#else
{

	real8 		burg[3], vec[3], mergePos[3],TBInGrain;
    real8       nodePos[3], oldPos[3], vec1[3], vec2[3];
    real8       dis;
	Node_t 		*nbrNode1, *nbrNode2, *mergedNode;
	int			mergeStatus,i,j, globalOp,q,thisDomain;
	int         crossSlipNode = 0, type;
    Param_t     *param;

    param = home->param;
	thisDomain = home->myDomain;
	
	burg[0] = node->burgX[0];
	burg[1] = node->burgY[0];
	burg[2] = node->burgZ[0];
	
	nbrNode1 = GetNeighborNode(home,node,0);
	nbrNode2 = GetNeighborNode(home,node,1);
	if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL)
		return(0);
    if(nbrNode1->myTag.domainID != thisDomain || 
       nbrNode2->myTag.domainID != thisDomain)return(0);
	
	if(SurfaceNode(nbrNode1) && IsOnGB(nbrNode1) ||
	   SurfaceNode(nbrNode2) && IsOnGB(nbrNode2))return(0);
	globalOp = 1;

    nodePos[0] = node->x;
    nodePos[1] = node->y;
    nodePos[2] = node->z;
    FoldBox(param, nodePos, nodePos+1, nodePos+2);

    oldPos[0] = node->oldx;
    oldPos[1] = node->oldy;
    oldPos[2] = node->oldz;
    PBCPOSITION(param, nodePos[0], nodePos[1], nodePos[2], 
                       oldPos, oldPos+1, oldPos+2);
    
    vec[0] = nodePos[0] - oldPos[0];
    vec[1] = nodePos[1] - oldPos[1];
    vec[2] = nodePos[2] - oldPos[2];
    if(Normal(vec) < 1E-2)return(0);

    NormalizeVec(vec);
    type = LinePlaneIntersection(vec, oldPos, GBPlane->normal, GBPlane->point, 
                                 vec1, vec2);
    
    switch(type){
        case POINT_INTERSECTION:
            VECTOR_COPY(mergePos, vec2);
            break;
        case LINE_INTERSECTION:
            vec[0] = 0.5*(nodePos[0] + oldPos[0]);
            vec[1] = 0.5*(nodePos[1] + oldPos[1]);
            vec[2] = 0.5*(nodePos[2] + oldPos[2]);
            PointPlaneIntersection(vec, GBPlane->normal, GBPlane->point,
                                   mergePos, &dis);
            break;
        default:
            Fatal("PenetrationLast: wrong intersection type");

    }
    
    if(BeforeMergeNodes(home, nbrNode1, nbrNode2, mergePos, &dis)){
        if(fabs(dis) > 5.0){
            printf("Warning in PenetrationLast, dis %f\n", dis);
        }
    }else{
        return(0);
    } 

	PullBack(home,node);
	FoldBox(home->param,&mergePos[0],&mergePos[1],&mergePos[2]);

	MergeNode(home, OPCLASS_SEPARATION, nbrNode1,nbrNode2,
			  mergePos, &mergedNode, &mergeStatus, globalOp);
	
	if((mergeStatus & MERGE_SUCCESS) == 0){
	    return(0);
	}
    
    if(mergedNode != (Node_t *)NULL){
        SplitIntersectNode(home, mergedNode, OPCLASS_SEPARATION, 1);
    }
	return(1);
}
#endif

int PenetrationSurface(Home_t *home, Node_t *node, GB_t *gb)
#if CROSS_SLIP_GB_MODEL
{
    int         type;
    real8       n1[3], p1[3], n2[3], p2[3], n3[3], p3[3];
    real8       lDir[3], lPos[3], vec[3], newPos[3];
    Node_t      *nbrNode;
	Param_t 	*param;

	param = home->param;

    if(!SurfaceNode(node))Fatal("in %s at %d", __FILE__, __LINE__);

	nbrNode = GetNeighborNode(home, node, 0);
    if(nbrNode == (Node_t *)NULL)Fatal("in %s at %d", __FILE__, __LINE__);
    if(nbrNode->myTag.domainID != home->myDomain || !IsOnGB(nbrNode))return(0);
    if(!SurfaceNode(node))Fatal("in %s at %d", __FILE__, __LINE__);

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
	RepositionNode(home, newPos, &node->myTag, 1);
    AssignGrain(home, node, -1, nbrNode->onGB);

    return(1);
}

#else
{
    int         type;
    real8       n1[3], p1[3], n2[3], p2[3], n3[3], p3[3];
    real8       lDir[3], lPos[3], vec[3], newPos[3];
    Node_t      *nbrNode;
	Param_t 	*param;

	param = home->param;

    if(!SurfaceNode(node))Fatal("in %s at %d", __FILE__, __LINE__);

	nbrNode = GetNeighborNode(home, node, 0);
    if(nbrNode == (Node_t *)NULL)Fatal("in %s at %d", __FILE__, __LINE__);
    if(nbrNode->myTag.domainID != home->myDomain || !IsOnGB(nbrNode))return(0);
    if(!SurfaceNode(node))Fatal("in %s at %d", __FILE__, __LINE__);

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
#endif


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
	
    pene.nSets = 0;
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

        if(node->numNbrs == 1 && node->constraint == PINNED_NODE)continue;
        if(currInGrain == oldInGrain || currInGrain < 0)continue;

        if(!param->peneTB){
            PullBack(home, node);
            continue;
        }

        if(node->numNbrs > 2){
            if(FGlidePerLength(home, node) > param->GBStrength){
                vec1[0] = node->vX;
                vec1[1] = node->vY;
                vec1[2] = node->vZ;
                vec2[0] = node->oldx;
                vec2[1] = node->oldy;
                vec2[2] = node->oldz;
                NormalizeVec(vec1);
                oldGrain  = home->grainKeys[oldInGrain];
                currGrain = home->grainKeys[currInGrain];
                if(!GetPublicGB(oldGrain, currGrain, &gb)){
                    Fatal("in %s at %d", __FILE__, __LINE__);
                }
                if(LinePlaneIntersection(vec1, vec2, gb->n, gb->p, vec, newPos)
                   == POINT_INTERSECTION){
                    FoldBox(param, newPos, newPos+1, newPos+2);
	                RepositionNode(home, newPos, &node->myTag, 1);
                    AssignGrain(home, node, -1, gb->GBID);
                }else{
                    PullBack(home, node);
                }
    
            }else{
                PullBack(home, node);
            }
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
        currGrain = home->grainKeys[currInGrain];
        if(!GetPublicGB(oldGrain, currGrain, &gb)){
            Fatal("in %s at %d", __FILE__, __LINE__);
        }
        VECTOR_COPY(GBPlane.normal, gb->n);
        VECTOR_COPY(GBPlane.point, gb->p);

/*	
 *	    There are 4 scenarios: 
 *	  	PenetrationFirst: the first penetrating node.
 *		PenetrationSide: one side penetrating nodes.
 *		PenetrationLast: the last node.
 *		PenetrationSurface: for the surface node.
 */
		
        if(SurfaceNode(node)){
            if(node->numNbrs != 1){
                PullBack(home, node);
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
	            PullBack(home, node);
	        }
	        continue;
	    }
		    
	    nbrNode1 = GetNeighborNode(home,node,0);
	    nbrNode2 = GetNeighborNode(home,node,1);
			
	    if(nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL){
	        printf("warning: nbrNode1==(Node_t *)NULL || nbrNode2==(Node_t *)NULL\n");
	        PullBack(home,node);
	        continue;
	    }		
	    if(nbrNode1 == nbrNode2) {
	        printf("warning: nbrNode1 == nbrNode2\n");
	        PullBack(home,node); 
	        continue;
	    }		

/*
 *      If the neighbor is not only intersect GB node, it is also junction node,
 *      pull back the penetrated node.
 */
        numNbrOnGB = 0;
        nbr1OnGB = NodeOnGB(home, nbrNode1);
        nbr2OnGB = NodeOnGB(home, nbrNode2);

	    if(IsOnGB(nbrNode1)) {

#if CROSS_SLIP_GB_MODEL
#else
	        if((nbr1OnGB & TB_NODE) == 0 && nbrNode1->numNbrs == 3){
	           PullBack(home, node);
	           continue;
	        } 
#endif
	        numNbrOnGB++;	
	    }

	    if(IsOnGB(nbrNode2)) {
#if CROSS_SLIP_GB_MODEL
#else
	        if((nbr2OnGB & TB_NODE) == 0 && nbrNode2->numNbrs == 3){
	           PullBack(home, node);
	           continue;
	        } 
#endif
	        numNbrOnGB++;	
	    }
			
	    if(numNbrOnGB != 2 && crossSlipNode == 1){
/*
 *		CrossSlip node is only allowed to penetrationLast. 
 */
            PullBack(home,node);
            continue;
        }
            	
	    if(numNbrOnGB == 0 && crossSlipNode != 1){
            if(!GetPenResults(home, node, &pene, gb)){
		        PullBack(home,node);
		        continue;
		    }
        }
/*		
 *	    Node: this energy rule is used for nonscrew dislocation, 
 *	    since we donot know the resultant slip plane for nonscrew 
 *	    dislocation at this time. 
 */

#if ENERGY_RULER
	    if(EnergyRule(home, node, &GBPlane) == 0){
            PullBack(home,node);
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
 *	    Energy rule for screw dislocation with at least a nbr on GB, so we can
 *		determine the resultant slip plane.	
 */
            if(ScrewEnergyRule(home, node, gb)==0){
	
                PullBack(home,node);
                continue;
            }
        }

        if(numNbrOnGB == 0){
            if(!PenetrationFirst(home, &pene)){
                PullBack(home,node);
                continue;
            }
        }else if(numNbrOnGB == 1){
            if(!PenetrationSide(home, node, oldGrain, currGrain, &GBPlane)){
                PullBack(home,node);
                continue;
		    }
        }else if(numNbrOnGB == 2){
            if(!PenetrationLast(home, node, oldGrain, currGrain, &GBPlane)){
                PullBack(home,node);
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
		
#if 1

        SplitIntersectNode(home, node, OPCLASS_SEPARATION, 1);
#else
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

        if(node->numNbrs < 2){
		    if(node->numNbrs == 1 && SurfaceNode(node) && 
		       IsOnGB(node)){
		        nbr1 = GetNeighborNode(home, node, 0);
		        if(nbr1 != (Node_t *)NULL){
			        if(!IsOnGB(nbr1)){
			            AssignGrain(home, node, nbr1->inGrain, -1);
			        }
		        }
		    }
		    continue;
	    }

        if(node->numNbrs == 2){
            nbr1 = GetNeighborNode(home,node,0);
            nbr2 = GetNeighborNode(home,node,1);
            if(nbr1==(Node_t *)NULL||nbr2==(Node_t *)NULL) continue;
            if(!(!IsOnGB(nbr1) && !IsOnGB(nbr2))) continue;
        }
			
        for(j=0;j<node->numNbrs;j++){
            nbr1 = GetNeighborNode(home,node,j);
            if(nbr1 == (Node_t *)NULL) continue;
            if(IsOnGB(nbr1))continue;
            if(nbr1->myTag.domainID != node->myTag.domainID) continue;
                	
            for(k=j+1;k<node->numNbrs;k++){
                nbr2 = GetNeighborNode(home,node,k);
                if(nbr2 == (Node_t *)NULL) continue;
	
		        if(nbr1->inGrain != nbr2->inGrain)continue;
                if(nbr1 == nbr2) continue;
                if(IsOnGB(nbr2))continue;
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
#endif
	}

	for(i=0;i<home->newNodeKeyPtr;i++){
        if((node = home->nodeKeys[i]) == (Node_t *)NULL)continue;
        if(IsOnGB(node))continue;
	    if(SurfaceNode(node))continue;
        NodeLocation(home, node->x, node->y, node->z, currInGrain, currOnGB);
        if(node->inGrain != currInGrain){
            PrintTag(node);
	        PullBack(home,node);
            printf("HandleGrainBoundary: a penetrated not coopered\n");
        }
    }
	return;

}

void DoResistance(Home_t *home, Node_t *node, real8 *nForce)
{
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
            NodeLocation(home, nodePos[0], nodePos[1], nodePos[2], inGrain, onGB);
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
            NodeLocation(home, nodePos[0], nodePos[1], nodePos[2], inGrain, onGB);
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

            NodeLocation(home, interPos[0], interPos[1], interPos[2], inGrain, onGB);
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


void SplitIntersectNode(Home_t *home, Node_t *node, int opClass, int globalOp){

    int         i, j, BakOnGB;    
    Node_t      *nbr;
    real8       normal[3], crossProd[3];
    real8       newPos[3];
    GB_t        *gb;

    std::vector<int>    vec(2);
    std::vector<std::vector<int> >   nbrInGrain;

    if(home->myDomain != node->myTag.domainID)return;
    if(!IsOnGB(node))return;

    gb = home->GBKeys[node->onGB];

    for(i=0; i<node->numNbrs; i++){
	    nbr = GetNeighborNode(home, node, i);
        if(nbr == (Node_t *)NULL)continue;
        if(IsOnGB(nbr))continue;

        normal[0] = node->nx[i];
        normal[1] = node->ny[i];
        normal[2] = node->nz[i];
        cross(normal, gb->n, crossProd);
        
        if(Normal(crossProd) < 1.0E-4){
            printf("Warning: in SplitIntersectNode at %d", __LINE__);
            return;
        }
        
        for(j=0; j<nbrInGrain.size(); j++){
            if(nbr->inGrain == nbrInGrain[j][0])break;
        }
        if(j<nbrInGrain.size()){
            nbrInGrain[j].push_back(i);
        }else{
            vec[0] = nbr->inGrain;
            vec[1] = i;
            nbrInGrain.push_back(vec);
        }
    }
    if(nbrInGrain.size() == 0)return;

    for(i=0; i<nbrInGrain.size(); i++){
        if(nbrInGrain[i].size()-1 == node->numNbrs){
            if(!NewPosToGB(home, node, nbrInGrain[i][0], node->onGB, newPos)){
                PrintNode(node);
                printf("Warning in SplitIntersectNode at %d\n", __LINE__);
            }

            if(globalOp){
                AssignGrain(home, node, nbrInGrain[i][0], -1);
            }else{
                node->inGrain = nbrInGrain[i][0];
                node->onGB = -1;
            }

	        RepositionNode(home, newPos, &node->myTag, globalOp);
            return;
        }

        if(nbrInGrain[i].size() <= 2)continue;

        if(nbrInGrain[i].size()-1 < node->numNbrs){
            real8   nodePos[3], nodeVel[3], newVel[3];
            Node_t  *splitNode1, *splitNode2;
            int     splitStatus, nArms, *armList;
            real8   resBurg[3];

            nArms = nbrInGrain[i].size()-1;
            armList = (int *)malloc(nArms*sizeof(int));
            
            VECTOR_ZERO(resBurg);
            for(j=0; j<nArms; j++){
                armList[j] = nbrInGrain[i][1+j];
                resBurg[0] += node->burgX[armList[j]];
                resBurg[1] += node->burgY[armList[j]];
                resBurg[2] += node->burgZ[armList[j]];
            }
            if(Normal(resBurg) > 1.0E-2)continue;
    
            nodePos[0] = node->x;
            nodePos[1] = node->y;
            nodePos[2] = node->z;

            nodeVel[0] = node->vX;
            nodeVel[1] = node->vY;
            nodeVel[2] = node->vZ;
            VECTOR_ZERO(newVel);
            
            splitStatus = SplitNode(home, opClass,
                                    node, nodePos, nodePos,
                                    nodeVel, newVel,
                                    nArms,
                                    armList, globalOp, &splitNode1,
                                    &splitNode2, 0);

            free(armList);
            if (splitStatus != SPLIT_SUCCESS)continue;

            if(!NewPosToGB(home, splitNode2, nbrInGrain[i][0], node->onGB, newPos)){
                PrintNode(splitNode2);
                Tecplot(home, (char*)"NewPosToGB", 0, 0, 0 ,0, 0, 1);
                Fatal("Warning in SplitIntersectNode at %d\n", __LINE__);
                return;
            }
	        RepositionNode(home, newPos, &splitNode2->myTag, globalOp);
            
            if(globalOp){
                AssignGrain(home, splitNode2, nbrInGrain[i][0], -1);
            }else{
                splitNode2->inGrain = nbrInGrain[i][0];
                splitNode2->onGB = -1;
            }

            if(GetArmID(home, splitNode1, splitNode2) != -1){
                PrintNode(splitNode1);
                PrintNode(splitNode2);
                printf("Warning in SplitIntersectNode at %d", __LINE__);
            }

            if(node != splitNode1)return;
        }
    }
    
    return;
}



void GBCrossSlip(Home_t *home)
{
    int     i, j, k, nBurgsSplit, nNbrsOnGB;
    int     opClass, thisDomain;
    int     inGrain, onGB, globalOp;
    int     type, armID, normalID;
    int     maxPowerSet;
    int     mergeStatus;
    bool    firstCrossSlipNode, sideCrossSlipNode;
    bool    breakLoop;
    real8   test1, test2, test3, test1Max, test2Max, test3Max;
    real8   nodePos[3], nbr1Pos[3], nbr2Pos[3], gbNbrPos[3];
    real8   vec1[3], vec2[3], vec3[3], vec[3], crossProd[3];
    real8   noiseFactor, weightFactor, fnodeThreshold;
    real8   eps, thetacrit, sthetacrit, s2thetacrit, areamin;
    real8   dis, burg[3], normal[3];
    real8   newNodePos[3], newNbr1Pos[3], newNbr2Pos[3], newGBNbrPos[3];
    real8   barrierFactor, oriDis, newVel[3];
    real8   newPos[3],maxPowerDis ;
    real8  segForceOrig[4][3];
    Line_t  interLine;
    Node_t  *node, *nbr1, *nbr2, *nbr, *gbNbr, *newNode;
    Node_t  *mergedNode;
    Param_t *param;

    std::vector<real8>                                      powerDis;
    std::vector<std::vector<double> >                       burgList;
    std::vector<std::vector< Dislocation_t> >               disArray; 
    std::vector<std::vector<Dislocation_t> >::iterator      itDisList;
    std::vector<std::vector<std::vector<real8> > >          disPosArray;

    BurgInfo_t  *gbd;
    GB_t        *gb;

    param = home->param;
    thisDomain = home->myDomain;

    eps = 1.0e-06;
    thetacrit = 2.0 / 180.0 * M_PI;
    sthetacrit = sin(thetacrit);
    s2thetacrit = sthetacrit * sthetacrit;
    areamin = param->remeshAreaMin;

/*
 *  In order to cross-slip a node, the force on the cross-slip
 *  plane must exceed a value based on the force on the current
 *  plane plus an amount that should guarantee the force on
 *  the cross-slip plane is not just in the noise.  The next
 *  two values are used in calculating that needed force.
 */
    noiseFactor=1e-05;
    weightFactor=1.0;

/*
 *  For purposes of determining segment ownership, we need to treat
 *  cross slip operations the same as we do during collision handling.
 */
    opClass = OPCLASS_COLLISION;

/*
 *  Loop through all 2-nodes native to this domain
 */
    for (i = 0; i < home->newNodeKeyPtr; i++) {

        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
            continue;
        }

        if (node->constraint != UNCONSTRAINED) {
            continue;
        }

        if(node->numNbrs < 2)continue;
        if(!IsOnGB(node))continue;

        nBurgsSplit = 0;
        nNbrsOnGB = 0;
        for(j=0; j<node->numNbrs; j++){
            nbr = GetNeighborNode(home, node, j);
            if(nbr == (Node_t *)NULL){
                printf("WARNING: Neighbor not found at %s line %d\n",
                        __FILE__, __LINE__);
                break;
            }
            if(!DomainOwnsSeg(home, opClass, thisDomain, &node->nbrTag[j])){
                break;
            }

            if(IsOnGB(nbr)){
                nNbrsOnGB++;
                gbNbr = nbr;
            }

            if(sqrt(node->burgX[j]*node->burgX[j]+
               node->burgY[j]*node->burgY[j]+
               node->burgZ[j]*node->burgZ[j]) > 
               SHOCKLEY_BURG_MAG_FCC+1.0E-2) nBurgsSplit++;

        }
        if(j != node->numNbrs)continue;

        if(nBurgsSplit == 0 || nBurgsSplit > 2 || nNbrsOnGB == 0)continue;

        firstCrossSlipNode = (node->numNbrs == 2 && nNbrsOnGB == 2 && nBurgsSplit == 2);
        sideCrossSlipNode = (nNbrsOnGB == 1 && nBurgsSplit == 1);
        if(!(firstCrossSlipNode || sideCrossSlipNode))continue;
        
        nodePos[0] = node->x;   nodePos[1] = node->y;   nodePos[2] = node->z;
        
        gb = home->GBKeys[node->onGB];
        gbd = &(gb->burgData);

        if(firstCrossSlipNode){
            nbr1 = GetNeighborNode(home, node, 0);
            nbr2 = GetNeighborNode(home, node, 1);
            if(nbr1->myTag.domainID != thisDomain ||
               nbr2->myTag.domainID != thisDomain)continue;
            
            nbr1Pos[0] = nbr1->x;   nbr1Pos[1] = nbr1->y;   nbr1Pos[2] = nbr1->z;
            nbr2Pos[0] = nbr2->x;   nbr2Pos[1] = nbr2->y;   nbr2Pos[2] = nbr2->z;

            PBCPOSITION(param, nodePos[0], nodePos[1], nodePos[2],
                        nbr1Pos, nbr1Pos+1, nbr1Pos+2);
            PBCPOSITION(param, nodePos[0], nodePos[1], nodePos[2],
                        nbr2Pos, nbr2Pos+1, nbr2Pos+2);

//            if(GetPlaneID_FCC(home, node, nbr1) != GetPlaneID_FCC(home, node, nbr2))continue;

            vec1[0] = nbr1Pos[0] - nodePos[0];
            vec1[1] = nbr1Pos[1] - nodePos[1];
            vec1[2] = nbr1Pos[2] - nodePos[2];

            vec2[0] = nbr2Pos[0] - nodePos[0];
            vec2[1] = nbr2Pos[1] - nodePos[1];
            vec2[2] = nbr2Pos[2] - nodePos[2];

            vec3[0] = nbr2Pos[0] - nbr1Pos[0];
            vec3[1] = nbr2Pos[1] - nbr1Pos[1];
            vec3[2] = nbr2Pos[2] - nbr1Pos[2];
            
            test1Max = DotProduct(vec1,vec1);
            test2Max = DotProduct(vec2,vec2);
            test3Max = DotProduct(vec3,vec3);
        
            for(j=0; j<gbd->numKeyLines; j++){
                test3 = DotProduct(vec3, gbd->keyLineList[j]);
                test3 = test3 * test3;
                if(test3Max - test3 > test3Max * s2thetacrit)continue;

                test1 = DotProduct(vec1, gbd->keyLineList[j]);    
                test1 = test1 * test1;
                if(test1Max - test1 > 4.0 * test1Max * s2thetacrit)continue;

                test2 = DotProduct(vec2, gbd->keyLineList[j]);    
                test2 = test2 * test2;
                if(test2Max - test2 > 4.0 * test2Max * s2thetacrit)continue;

                break;
            }
            if(j == gbd->numKeyLines)continue;

        }else if(sideCrossSlipNode){
            gbNbrPos[0] = gbNbr->x;     gbNbrPos[1] = gbNbr->y;     gbNbrPos[2] = gbNbr->z;
            PBCPOSITION(param, nodePos[0], nodePos[1], nodePos[2],
                        gbNbrPos, gbNbrPos+1, gbNbrPos+2);
            
            vec1[0] = gbNbrPos[0] - nodePos[0];
            vec1[1] = gbNbrPos[1] - nodePos[1];
            vec1[2] = gbNbrPos[2] - nodePos[2];
            test1Max = DotProduct(vec1, vec1);
    
            for(j=0; j<gbd->numKeyLines; j++){
                test1 = DotProduct(vec1, gbd->keyLineList[j]);
                test1 = test1 * test1;
                if(test1Max - test1 > test1Max * s2thetacrit)continue;
                break;
            }
            if(j == gbd->numKeyLines)continue;
        }else{
            continue;
        }

        
        armID = GetArmID(home, node, gbNbr);
        burg[0] = node->burgX[armID];
        burg[1] = node->burgY[armID];
        burg[2] = node->burgZ[armID];
        normal[0] = node->nx[armID];
        normal[1] = node->ny[armID];
        normal[2] = node->nz[armID];

        NormalizeVec(vec1);
        VECTOR_COPY(interLine.direction,gbd->keyLineList[j]);
        type = PointPlaneIntersection(nodePos, gb->n, gb->p, interLine.point, &dis);
        if(fabs(dis) > 2)printf("Warning in %s at %d\n", __FILE__, __LINE__);

        if(firstCrossSlipNode){
            std::vector<std::vector< Dislocation_t> >().swap(disArray);
            std::vector<std::vector<double> >().swap(burgList);

            if(!GeometricCriterion(gbd, burg, interLine.direction, normal, 1,
                   burgList, disArray))continue;

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
            if(disArray.size() == 0)continue;
                      
            VECTOR_COPY(newNodePos, interLine.point);
            PointLineIntersection(nbr1Pos, interLine.direction, interLine.point,
                                 newNbr1Pos, &dis);
            PointLineIntersection(nbr2Pos, interLine.direction, interLine.point,
                                 newNbr2Pos, &dis);

            SaveCrossSlipInfo(node, nbr1, nbr2, 0, 1,
                              segForceOrig, nodePos, nbr1Pos,
                              nbr2Pos);

            oriDis = node->fX * node->vX + 
                     node->fY * node->vY + 
                     node->fZ * node->vZ;  

            oriDis *= (1.0+noiseFactor);

	        RepositionNode(home, newNodePos, &node->myTag, 0);
	        RepositionNode(home, newNbr1Pos, &nbr1->myTag, 0);
	        RepositionNode(home, newNbr2Pos, &nbr2->myTag, 0);

            powerDis.resize(disArray.size());
            disPosArray.resize(disArray.size());
            for(j=0; j<disArray.size(); j++){
                powerDis[j] = -oriDis;
                disPosArray[j].resize(disArray[j].size());
                for(k=0; k<disArray[j].size(); k++){
                    disPosArray[j][k].resize(3);
	                ChangeArmBurg(home, nbr1, &node->myTag, 
                                  burgList[disArray[j][k].burgID][0],
                                  burgList[disArray[j][k].burgID][1],
                                  burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  0, DEL_SEG_NONE);
                    
	                ChangeArmBurg(home, node, &nbr2->myTag, 
                                  burgList[disArray[j][k].burgID][0],
                                  burgList[disArray[j][k].burgID][1],
                                  burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  0, DEL_SEG_NONE);

	                ChangeArmBurg(home, nbr2, &node->myTag, 
                                  -burgList[disArray[j][k].burgID][0],
                                  -burgList[disArray[j][k].burgID][1],
                                  -burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  0, DEL_SEG_NONE);
                    
	                ChangeArmBurg(home, node, &nbr1->myTag, 
                                  -burgList[disArray[j][k].burgID][0],
                                  -burgList[disArray[j][k].burgID][1],
                                  -burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  0, DEL_SEG_NONE);

                    node->inGrain = gbd->planeInGrain[disArray[j][k].planeID];
                    node->onGB = (node->inGrain == -1) ? gb->GBID : -1; 

		            SetOneNodeForce(home, node);
                    EvaluateMobility(home, node);
                    vec[0] = node->vX;
                    vec[1] = node->vY;
                    vec[2] = node->vZ;

                    vec[0] *= (param->nextDT/10.0);
                    vec[1] *= (param->nextDT/10.0);
                    vec[2] *= (param->nextDT/10.0);

                    vec[0] += newNodePos[0];
                    vec[1] += newNodePos[1];
                    vec[2] += newNodePos[2];

                    VECTOR_COPY(disPosArray[j][k], vec);

                    if(!IsOnGB(node)){
                        NodeLocation(home, vec[0], vec[1], vec[2], inGrain, onGB);
                        if(onGB > -1 || inGrain != node->inGrain)break;
                    }

                    powerDis[j] += (node->fX * node->vX +
                                    node->fY * node->vY +
                                    node->fZ * node->vZ);
    
                }
                if(k<disArray[j].size())powerDis[j] = -oriDis;
            }

            maxPowerDis = -oriDis;
            maxPowerSet = -1;
            for(j=0; j<disArray.size(); j++){
               if(powerDis[j] > maxPowerDis){
                    maxPowerSet = j;
                    maxPowerDis = powerDis[j];
                }
            }

            if(maxPowerSet == -1 || maxPowerDis < 0){
                RestoreCrossSlipForce(node, nbr1, nbr2, 0,
                                      1, segForceOrig);
	            RepositionNode(home, nodePos, &node->myTag, 0);
	            RepositionNode(home, nbr1Pos, &nbr1->myTag, 0);
	            RepositionNode(home, nbr2Pos, &nbr2->myTag, 0);

	            ChangeArmBurg(home, nbr1, &node->myTag, 
                              burg[0],burg[1],burg[2],normal[0],normal[1],normal[2],
			                  0, DEL_SEG_NONE);
	            ChangeArmBurg(home, node, &nbr2->myTag, 
                              burg[0],burg[1],burg[2],normal[0],normal[1],normal[2],
			                  0, DEL_SEG_NONE);

	            ChangeArmBurg(home, nbr2, &node->myTag, 
                              -burg[0],-burg[1],-burg[2],normal[0],normal[1],normal[2],
			                  0, DEL_SEG_NONE);
	            ChangeArmBurg(home, node, &nbr1->myTag, 
                              -burg[0],-burg[1],-burg[2],normal[0],normal[1],normal[2],
			                  0, DEL_SEG_NONE);

                node->inGrain = -1;
                node->onGB = gb->GBID;
                continue;
            }else{
	            FoldBox(param,newNbr1Pos,newNbr1Pos+1,newNbr1Pos+2);
	            FoldBox(param,newNbr2Pos,newNbr2Pos+1,newNbr2Pos+2);
	            RepositionNode(home, newNbr1Pos, &nbr1->myTag, 1);
	            RepositionNode(home, newNbr2Pos, &nbr2->myTag, 1);
            }

            j = maxPowerSet;
            globalOp = 1;

            for(k=0; k<disArray[j].size(); k++){
                if(k==0){
	                ChangeArmBurg(home, nbr1, &node->myTag, 
                                  burgList[disArray[j][k].burgID][0],
                                  burgList[disArray[j][k].burgID][1],
                                  burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  globalOp, DEL_SEG_NONE);
                    
	                ChangeArmBurg(home, node, &nbr2->myTag, 
                                  burgList[disArray[j][k].burgID][0],
                                  burgList[disArray[j][k].burgID][1],
                                  burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  globalOp, DEL_SEG_NONE);

	                ChangeArmBurg(home, nbr2, &node->myTag, 
                                  -burgList[disArray[j][k].burgID][0],
                                  -burgList[disArray[j][k].burgID][1],
                                  -burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  globalOp, DEL_SEG_NONE);
                    
	                ChangeArmBurg(home, node, &nbr1->myTag, 
                                  -burgList[disArray[j][k].burgID][0],
                                  -burgList[disArray[j][k].burgID][1],
                                  -burgList[disArray[j][k].burgID][2],
                                  gbd->planeList[disArray[j][k].planeID][0],
                                  gbd->planeList[disArray[j][k].planeID][1],
                                  gbd->planeList[disArray[j][k].planeID][2],
				                  globalOp, DEL_SEG_NONE);

                    node->inGrain = gbd->planeInGrain[disArray[j][k].planeID];
                    node->onGB = (node->inGrain == -1) ? gb->GBID : -1; 

                    VECTOR_COPY(newNodePos, disPosArray[j][k]);
	                FoldBox(param,newNodePos,newNodePos+1,newNodePos+2);
	                RepositionNode(home, newNodePos, &node->myTag, globalOp);
                    AssignGrain(home, node, node->inGrain, node->onGB);
                    MarkNodeForceObsolete(home, node);        
                    MarkNodeNoCollisions(home, node);        
                    MarkNodeNoMeshCoarsen(home, node);        
                }else{
                    VECTOR_COPY(vec1, burgList[disArray[j][k].burgID]);
                    VECTOR_COPY(vec2, gbd->planeList[disArray[j][k].planeID]);
                   
                    inGrain = gbd->planeInGrain[disArray[j][k].planeID];
                    onGB = (inGrain == -1) ? gb->GBID : -1; 

                    VECTOR_COPY(newPos, disPosArray[j][k]);
	                FoldBox(param,newPos,newPos+1,newPos+2);
                    InsertConnection(home, nbr1, nbr2, globalOp, inGrain, onGB,
                         newPos, vec1, vec2, &newNode);
                    MarkNodeForceObsolete(home, newNode);        
                    MarkNodeNoCollisions(home, newNode);        
                    MarkNodeNoMeshCoarsen(home, newNode);        

                } 
            }
            
	        MarkNodeForceObsolete(home, nbr1);        
	        MarkNodeNoCollisions(home, nbr1);        
	        MarkNodeNoMeshCoarsen(home, nbr1);        
        
	        MarkNodeForceObsolete(home, nbr2);        
	        MarkNodeNoCollisions(home, nbr2);        
	        MarkNodeNoMeshCoarsen(home, nbr2);        
        
	        nbr1->oldvX = 0.0;
	        nbr1->oldvY = 0.0;
	        nbr1->oldvZ = 0.0;
	    
	        nbr2->oldvX = 0.0;
	        nbr2->oldvY = 0.0;
	        nbr2->oldvZ = 0.0;
	    
	        nbr1->vX = 0.0;
	        nbr1->vY = 0.0;
	        nbr1->vZ = 0.0;
	
	        nbr2->vX = 0.0;
	        nbr2->vY = 0.0;
	        nbr2->vZ = 0.0;
            
        }else if(sideCrossSlipNode){
            PointLineIntersection(gbNbrPos, interLine.direction, interLine.point,
                                  newPos, &dis);

            FoldBox(param,newPos, newPos+1, newPos+2);
            if(BeforeMergeNodes(home, node, gbNbr, newPos, &dis)){
                if(fabs(dis) > 5.0){
                printf("Warning in PenetrationSide, dis %f\n", dis);
                }
            }else{
                continue;
            } 
            FoldBox(param,newPos, newPos+1, newPos+2);
            MergeNode(home, OPCLASS_SEPARATION, node, gbNbr, newPos,
                      &mergedNode, &mergeStatus, globalOp);

	        if((mergeStatus & MERGE_SUCCESS) ==0)continue;
    
            if(mergedNode != (Node_t *)NULL){
                SplitIntersectNode(home, mergedNode, opClass, globalOp);
            }
        }

        
    }  /* end loop over all nodes */

    return;
}


bool AdjustSplitedNodePosition(Home_t *home, Node_t *node, real8 *currPos, real8 *pos)
{
    
    VECTOR_ZERO(pos);
    if(home->param->polycrystal == 0)return(0);

    int     currInGrain, currOnGB;
    GB_t    *GB;

    if(IsOnGB(node))return(0);
    NodeLocation(home, currPos[0], currPos[1],
                 currPos[2], currInGrain, currOnGB);
                                    
    if(currInGrain != node->inGrain || currOnGB > -1){
        if(currOnGB == -1){
            if(GetPublicGB(home->grainKeys[node->inGrain],
                           home->grainKeys[currInGrain], &GB)){
                currOnGB = GB->GBID;
            }else{
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
        }
        if(!NewPosToGB(home, node, node->inGrain, currOnGB, pos)){
            Fatal("in %s at %d", __FILE__, __LINE__);
        }
        return(1);
    }
    return(0);
}
