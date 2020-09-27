/**************************************************************************
 *
 *  Function    : Mobility_Wei
 *  Author      : DeAn Wei
 *  Date        : 12.10.2019
 *  Description : Generic Mobility Law for polycrystalline
 *  Returns:    0 on success
 *              1 if velcoity could not be determined
 *
 ***************************************************************************/

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include "Math.h"
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif

// unit: (Pa.s)^-1
inline real8 Mobility(real8 b, real8 B0, real8 sig0, real8 sig1, real8 cs, real8 sig)
{
    real8 sig_hat=(sig-sig0)/sig1;
    real8 B1 = b/cs;
    real8 B=B0+B1*sig1*(sig_hat+log(1+pow(2.71828,-sig_hat)));
    return 1.0/B;
}

int Mobility_Wei(Home_t *home, Node_t *node)
{
    int numNonZeroLenSegs = 0;
    int numGlideConstraints;
    int grainNode;
    Param_t *param;
    real8 nVel[3], Veldotlcr, forceGB[3];
    int i, j, nc, nconstraint, nlc, TBBurg, armID21, nodeOnGB;
    real8 normx[MAX_NBRS], normy[MAX_NBRS], normz[MAX_NBRS];
    real8 lineX[MAX_NBRS], lineY[MAX_NBRS], lineZ[MAX_NBRS];
    real8 a, b, velocity;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB, lxb;
    real8 lcx, lcy, lcz, normdotlc;
    Node_t *nbr;
    real8 MobScrew, MobEdge, MobTwinning, Mob;
    real8 bx, by, bz, br, dangle;
    real8 glideConstraintsCrystal[3][3], imageForce[3];
    real8 nForce[3], glideForce[3], glideConstraints[3][3], n[3], sig;
    real8 factor = 1.0, pinForce;

    static real8 params[4][3]={
                               {2.56938E-5, 1.03497E8 , 1.4982E7  },
                               {1.68337E-5, 6.50097E7 , 8.16172E7 },
                               {1.18461E-5, 1.05257E8 , 8.9693E7  },
                               {1.84995E-5, 2.12353E8 , 3.20935E7 } 
                              };
    static real8 cs = 1630.0;

    param = home->param;
	
    bool zeroLenJunc;
	
    nc = node->numNbrs;
	
/*
 *  If node is 'pinned' in place by constraints, the node may not be moved so just zero the velocity
 *  and return
 */
    if ((node->constraint == PINNED_NODE) ||
        NodeHasSessileBurg(home, node)){
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
    for (i = 0; i < nc; i++) {
        if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
        dx = node->x - nbr->x;
        dy = node->y - nbr->y;
        dz = node->z - nbr->z;
        if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
            numNonZeroLenSegs++;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  Find total dislocation length times drag coefficent (LtimesB)
 */
    LtimesB=0; lxb=0;
    for(j=0;j<nc;j++){
        if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
        dx=nbr->x - node->x;
        dy=nbr->y - node->y;
        dz=nbr->z - node->z;
        ZImage (param, &dx, &dy, &dz) ;

/*
 *      If needed, rotate the line sense from the laboratory frame to
 *      the crystal frame.
 */
        if (param->useLabFrame) {
            real8 dTmp[3] = {dx, dy, dz};
            real8 dRot[3];
            VectorLabToCrystal(home, node, nbr, dTmp, dRot);
            dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
        }

        lr=sqrt(dx*dx+dy*dy+dz*dz);
        
        if (lr==0){ 
		  /* zero arm segment can happen after node split 
           * it is OK to have plane normal vector == 0
           * Skip (do nothing)
           */
        }else{
            if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0)){
              printf("Mobility_FCC_0: (%d,%d) glide plane norm = 0\n"
                     "for segment with nonzero length lr = %e!\n",
                     node->myTag.domainID, node->myTag.index, lr);
            }

            lx=dx/lr; ly=dy/lr; lz=dz/lr;
        
            bx = node->burgX[j];
            by = node->burgY[j];
            bz = node->burgZ[j];
/*
 *          If needed, rotate the burgers vector from the laboratory frame to
 *          the crystal frame.
 */
	        TBBurg = 0;
            if (param->useLabFrame) {
                real8 bTmp[3] = {bx, by, bz};
                real8 bRot[3];

		        if(BurgType_FCC(bTmp) == SHOCKLEY_BURG)TBBurg = 1;
                VectorLabToCrystal(home, node, nbr, bTmp, bRot);
                bx = bRot[0]; by = bRot[1]; bz = bRot[2];
            }       

	        br = sqrt(bx*bx+by*by+bz*bz);
            lxb += (lr*br);
	        bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

            armID21 = GetArmID(home, nbr, node);
            if(armID21<0){
                PrintNode(home, node);
                PrintNode(home, nbr);
                Fatal("%s: armID21 <0", __func__);
            }
            glideForce[0] = node->armfx[j] + nbr->armfx[armID21];
            glideForce[1] = node->armfy[j] + nbr->armfy[armID21];
            glideForce[2] = node->armfz[j] + nbr->armfz[armID21];

            n[0] = node->nx[j]; n[1] = node->ny[j]; n[2] = node->nz[j];
            NormalComponent(glideForce, n, glideForce);
            sig = Normal(glideForce)/br/lr;
 
	        dangle = fabs(bx*lx+by*ly+bz*lz);
            if(TBBurg){ 
                MobScrew =Mobility(param->burgMag*br, params[2][0], params[2][1], params[2][2], cs, sig);
                MobEdge = Mobility(param->burgMag*br, params[3][0], params[3][1], params[3][2], cs, sig);
	        }else{
                MobScrew =Mobility(param->burgMag*br, params[0][0], params[0][1], params[0][2], cs, sig);
                MobEdge = Mobility(param->burgMag*br, params[1][0], params[1][1], params[1][2], cs, sig);
            }
	        Mob=MobEdge+(MobScrew-MobEdge)*dangle;
	        LtimesB+=(lr/Mob);
	    }
    }
    LtimesB/=2; lxb/=2;

/*
 *  Add resitance of friction force of GB force
 */
#if 0
	DoResistance(home, node, nForce);
    NodeImageForce(home, node, nForce);
    nForce[0] += node->fX;
    nForce[1] += node->fY;
    nForce[2] += node->fZ;
#else
    nForce[0] = node->fX;
    nForce[1] = node->fY;
    nForce[2] = node->fZ;
#endif
	/* Velocity is simply proportional to total force per unit length */
	GetGlideConstraintList(home, node, &numGlideConstraints, glideConstraints, zeroLenJunc);
    ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, nForce);
#if 1
    NodeImageForce(home, node, imageForce);
    ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, imageForce);
    if(DotProduct(nForce,imageForce)<0){
        if(Normal(nForce)>Normal(imageForce)){
            nForce[0] -= imageForce[0];
            nForce[1] -= imageForce[1];
            nForce[2] -= imageForce[2];
        }else{
            node->vX =0.0;
            node->vY =0.0;
            node->vZ =0.0;
            return (0);
        }
    }
#endif

#ifdef TB_PINNING
    nodeOnGB = NodeOnGB(home, node);
    if((nodeOnGB &INTERSECT_NODE)>0 && (numGlideConstraints ==2) && (nodeOnGB & GB_NODE)==0){
//      pinForce = ((param->tbStrength[0] < 1.0E5)?param->tbStrength[0]*lxb:(4.81851E8/param->grainSize + 5.9E8)*lxb);
        pinForce = 0.15*param->shearModulus;
        ParallelComponent(nForce, glideConstraints[2], glideForce);
        if(Normal(glideForce)<=pinForce){
            VECTOR_ZERO(nVel);
        }else{
            glideForce[0] -= (pinForce*glideConstraints[2][0]);
            glideForce[1] -= (pinForce*glideConstraints[2][1]);
            glideForce[2] -= (pinForce*glideConstraints[2][2]);
	        nVel[0] = glideForce[0]/LtimesB;
	        nVel[1] = glideForce[1]/LtimesB;
	        nVel[2] = glideForce[2]/LtimesB;
        }
    }else{
	    nVel[0] = nForce[0]/LtimesB;
	    nVel[1] = nForce[1]/LtimesB;
	    nVel[2] = nForce[2]/LtimesB;
    }
#else
	 nVel[0] = nForce[0]/LtimesB;
	 nVel[1] = nForce[1]/LtimesB;
	 nVel[2] = nForce[2]/LtimesB;
#endif

    if(Normal(nVel) > 1E-10 && SurfaceNode(node)){
        real8   vec[3], nor[3];
        VECTOR_ZERO(nor);
        nor[node->constraint-1] = 1.0;
        VECTOR_COPY(vec, nVel);
        NormalizeVec(vec);
        if(fabs(DotProduct(vec, nor)) > 0.01){
            FormatMatrix(glideConstraints, "glideConstraints");
            Fatal("in Mobility %d", numGlideConstraints);
        }   
    }

#if 0
	if (param->uselabframe) {
	    int   gindex;
	 
	    for (gindex = 0; gindex < *numglideconstraints; gindex++) {
	 
	        matrix33vector3multiply(home->rotmatrixinverse,
	                                glideconstraints[gindex],
	                                glideconstraintscrystal[gindex]);
	    }
	 
	 } else {
		 matrix33copy(glideconstraints, glideconstraintscrystal);
	 }
	
	 applyconstraintstovelocity(*numglideconstraints,
				    glideconstraintscrystal, nvel);
#endif

/* Any dislocation with glide plane not {111} type can only move along its length
 * This rule includes LC junction which is on {100} plane
 */
	 
#if 0
	if(nlc==1){ /* only one line constraint */
            for(i=0;i<nc;i++){
                if((lineX[i]!=0)||(lineY[i]!=0)||(lineZ[i]!=0)){ 
   	            lcx = lineX[i];
	            lcy = lineY[i];
	            lcz = lineZ[i];
                    break;
                }
            }
    
            /* project velocity along line */
            nVel[0] = Veldotlcr*lcx;
	    nVel[1] = Veldotlcr*lcy;
	    nVel[2] = Veldotlcr*lcz;
    	
	    if (nconstraint<=0){	
	        Fatal("MobilityLaw_FCC_0: nconstraint <= 0, nlc = 1 is impossible!");
	    }else if(nconstraint >= 1){ /* a few plane constraints and one line constraint */
			
		for(i=0;i<nc;i++){
		    normdotlc = normx[i]*lcx + normy[i]*lcy + normz[i]*lcz;
		    if(fabs(normdotlc)>FFACTOR_ORTH){
			/* set velocity to zero if line is not on every plane */
			nVel[0] = nVel[1] = nVel[2] = 0;
			break;
		    }else{
		        /* do nothing. Skip */
		    }        
		}
	    }
        
    }else if (nlc>=2){
        /* Velocity is zero when # of independnet lc constratins is equal to or more than 2 */
        nVel[0] = nVel[1] = nVel[2] = 0;
    }else{ 
        /* nlc == 0, do nothing */
    }
#endif

	node->vX = nVel[0];
	node->vY = nVel[1];
	node->vZ = nVel[2];

	if(SurfaceNode(node))factor = 20.0;
	    
	velocity = sqrt(node->vX * node->vX+
			   node->vY * node->vY+
			   node->vZ * node->vZ);
	if(velocity > param->vshr/factor && param->vshr > 0){
	    node->vX = node->vX * param->vshr /velocity/factor ;
	    node->vY = node->vY * param->vshr /velocity/factor ;
	    node->vZ = node->vZ * param->vshr /velocity/factor ;
	}

    return(0);
}
