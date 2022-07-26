/**************************************************************************
 *
 *  Function    : Mobility_FCC_Nolinear
 *  Author      : DeAn Wei
 *  Description : This is alternate version of a generic mobility function
 *                for FCC materials.  It is very similar in function and
 *                structure to the BCC_glide mobility although glide planes
 *                rather than burgers vectors are used to identify junctions.
 *                Additionally, this function attempts to identify nodes
 *                that both have short segments and have reversed the 
 *                their direction of velocity.  When found, these nodes
 *                are artifically slowed for a single timestep in an effort
 *                to dampen motion of flickering/oscillating nodes.
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
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

inline real8 B(real8 k, real8 vmax, real8 vc, real8 vMag){
    real8 tem=(vMag<vc)?k*(vmax-vMag):k*(vmax-vc);
    return 1.0/tem;
}

inline real8 FBack(real8 k, real8 vmax, real8 vc, real8 b, real8 vMag){
    real8 invK=1.0/k;
    real8 frac;
    if(vMag<vc){
        frac = vmax/(vmax-vMag);
        return (invK*log(frac) - b*vMag);
    }else{
        frac = vmax/(vmax-vc);
        return (invK*log(frac) - b*vc);
    }
    return 0.0;
}

int Mobility_FCC_Nolinear(Home_t *home, Node_t *node)
{
	int     i, j, k, numNbrs; 
	int     numNonZeroLenSegs = 0, hasShortSeg = 0;
	int	    numGlideConstraints;
	real8   dx, dy, dz;
	real8   nx, ny, nz, nx2, ny2, nz2; 
	real8   temp;
	real8   mag, mag2, halfMag, invMag;
	real8   invbMag, bMag2, costheta, costheta2, sintheta2;
	real8   Bscrew2, Beclimb2;
	real8   invBscrew2, invBedge2;
	real8   shortSegCutoff;
	real8   eps = 1.0e-12;
	real8   tor=1.e-5;
	real8   shortSegRatio = 1.0;
	real8   nForce[3]={node->fX, node->fY, node->fZ}, nVel[3];
	real8   Btotal[3][3];
	real8   glideConstraints[3][3];
	real8   invBtotal[3][3];

	
	Node_t  *nbrNode;
    bool    zeroLenJunc;

	Param_t *param;
	param  = home->param;
	numNbrs = node->numNbrs;
	shortSegCutoff = 0.5 * param->minSeg;
    real8   v_max=0.0, nEffNbrs=0.0, vnbr;
/*
 *  If node is 'pinned' in place by constraints, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
	if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
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
	for (i = 0; i < numNbrs; i++) {
	    if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
	        continue;
	    }

	    dx = node->x - nbrNode->x;
	    dy = node->y - nbrNode->y;
	    dz = node->z - nbrNode->z;
	    
	    ZImage(param, &dx, &dy, &dz);

/*
 *      This assumes that lengths are given in terms of burgers vector units
 */
	    if ((dx*dx + dy*dy + dz*dz) > eps) {
	        numNonZeroLenSegs++;
	    }
        vnbr= nbrNode->vX*nbrNode->vX + nbrNode->vY*nbrNode->vY + nbrNode->vZ*nbrNode->vZ; 
        if(vnbr>0){
            v_max += vnbr;
            nEffNbrs += 1.0;
        }
	}
    if(home->cycle<2){
        v_max = param->vshr;
    }else{
        if(v_max>0 && nEffNbrs > 0){
            v_max = 3.0*sqrt(v_max)/nEffNbrs;
        }
        v_max = MAX(v_max, param->vshr/10.0);
        v_max = MIN(v_max, param->vshr);
    }

	if (numNonZeroLenSegs == 0) {
	    node->vX = 0.0;
	    node->vY = 0.0;
	    node->vZ = 0.0;
	    return(0);
	}

/*
 *  Get the normal for the first segment in both the lab and crystal frames
 */
    GetGlideConstraintList(home, node, &numGlideConstraints, glideConstraints, zeroLenJunc);

    if(numGlideConstraints==3){
	    node->vX = 0.0e0;
	    node->vY = 0.0e0;
	    node->vZ = 0.0e0;
        return 0;
    }

#ifndef IMAGE_FORCE
    real8 imageForce[3], imag, nmag, iDir[3], ndot;
    NodeImageForce(home, node, imageForce); // nodal image force parallel to sement glide dirs
    imag = Normal(imageForce);
    if(imag > 1.0){ // if has image force
        iDir[0] = imageForce[0]/imag;
        iDir[1] = imageForce[1]/imag;
        iDir[2] = imageForce[2]/imag;

        ndot = DotProduct(nForce, iDir); // nodal force magitude along image force dir

        if(ndot < -imag){ // if the dirs are oppsite and abs(ndot) > imag, add the image force
            CoordAdd(nForce, imageForce, nForce);
        }else if(ndot < 0){ // if the dirs are oppsite and image > abs(ndot), animate the force conponent along the image force dir
            nForce[0] -= ndot*iDir[0];
            nForce[1] -= ndot*iDir[1];
            nForce[2] -= ndot*iDir[2];
        }
    }
#endif

#ifdef JOG_PINNING
    real8 pinForce = param->jogCoeff*param->shearModulus;
    int     nbrInGrain;
    real8   pinDir[3], pdot;   

    if(numGlideConstraints == 2 && node->onGB>-1 && numNbrs>1){
        nbrInGrain = -1;
        for(i=0; i<numNbrs; i++){
            if((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL)continue;
            if(nbrNode->inGrain>-1){
                if(nbrInGrain == -1){
                    nbrInGrain = nbrNode->inGrain; 
                }else if(nbrNode->inGrain != nbrInGrain){
                    break;
                }
            }
        }       

        if(i != numNbrs){
            pdot = DotProduct(nForce, glideConstraints[2]);

            if(pdot<0){
                pdot       = -pdot;
                pinDir[0]  = -glideConstraints[2][0];
                pinDir[1]  = -glideConstraints[2][1];
                pinDir[2]  = -glideConstraints[2][2];
            }else{
                VECTOR_COPY(pinDir, glideConstraints[2]);
            }

            if(pdot > pinForce){
                nForce[0] -= (pinForce*pinDir[0]);                
                nForce[1] -= (pinForce*pinDir[1]);                
                nForce[2] -= (pinForce*pinDir[2]);                
            }else{
	            node->vX = 0.0;
	            node->vY = 0.0;
	            node->vZ = 0.0;
	            return(0);
            }
        }
    }
#endif
/*
 *  v = v_max*(1 - exp(-k*tau))
 *                    drag coefficent                             back          force 
 *  fdrag = tau*b = b/(k*(v_max - v)) * v + b/k * log(v_max/(v_max - v)) -  b/(k*(v_max - v)) * v
 */
    static real8 vs_max  = 1454.62292/param->burgMag, ks  = 9.24E-9;
    static real8 ve_max  = 1529.63902/param->burgMag, ke  = 8.58E-9;
    static real8 vts_max = 1585.57223/param->burgMag, kts = 13.7005E-9;  // ks=7.91E-9/(1/sqrt[3]) (Pa*b)^-1
    static real8 vte_max = 1877.26443/param->burgMag, kte = 10.808E-9;   // ks=6.24E-9/(1/sqrt[3]) (Pa*b)^-1

    static real8 cutoff_coe = 0.99;
    static real8 vs_max_c   = vs_max*cutoff_coe, ve_max_c = ve_max*cutoff_coe, vts_max_c = vts_max*cutoff_coe, vte_max_c = vte_max*cutoff_coe;
    static real8 Beclimb    = 15.6E-4;
    real8  Bline,  BlmBecl; 

/*
 *  Begin construction of the node drag matrix
 *  Loop over all arms of the node, adding each arm's contribution
 *  to the drag matrix.
 */
	real8   Bglide;
    real8   b[3], d[3], m[3], n[3];
    real8   bMag, vTest[3]={node->vX, node->vY, node->vZ};
/*
 *  As far as possible, the test velocity should be equal to the original velocity, the old velocity, 
 *  or the velocities of neighbor nodes, with decreasing precedence.
 */
    if(DotProduct(vTest, vTest)<1.0E-10){
        VectorCopy(vTest, node->oldvX, node->oldvY, node->oldvZ);
        if(DotProduct(vTest, vTest)<1.0E-10){
            for(i=0; i<numNbrs; i++){
                if((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL)continue;
                VectorCopy(vTest, nbrNode->vX, nbrNode->vY, nbrNode->vZ);
                if(DotProduct(vTest, vTest)>1.0E-10)break;
            }       
        }
    }


/*
 *  Limit the direction of the vTest to be within the constraint and less than the maximum
 */
    real8 vMag;
    ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, vTest);
    vMag = Sign(DotProduct(nForce,vTest))*Normal(vTest);
    if(vMag>v_max){VectorTimesSelf(vTest, (v_max/vMag));}

    bool  tdKink;
    real8 fdrag[3], bForce[3], vDelta[3];
    real8 Bi[3][3], vgTest[3], vgMag, fTol=1, fErr=fTol+1;
    real8 lineTensor[3][3], glideTensor[3][3], climbTensor[3][3];
    real8 bs,be,inv_bs,inv_be,sbForceMag, ebForceMag, vRatio_s, vRatio_e, lterm, tau, bForceMag;
    bool partialBurg;
    int iter, maxStep=30;
    real8 v1g, v1l, v2g, v2l, va1gs, va1ge, va2ge, va2gs;
    vMag = Sign(DotProduct(nForce,vTest))*Normal(vTest);
    for(iter=0; iter<maxStep &&  fErr>fTol; iter++){
/*
 *      Initialze total drag tensor and drag force and recaculate vMag
 */
        NegativeVector(nForce, fdrag);
	    Matrix33Init(Btotal);

		for (i = 0; i < numNbrs; i++) {
            VectorCopy(b, node->burgX[i], node->burgY[i], node->burgZ[i]);
			bMag2 = (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
            bMag = sqrt(bMag2);
            partialBurg = ((BurgType_FCC(b) == SHOCKLEY_BURG)?true:false);
            NormalizeVec(b);
    
/*
 *          Calculate the length of the arm and its tangent line direction
 */         
		    nbrNode = GetNeighborNode(home, node, i);
            
		    if (nbrNode == (Node_t *)NULL) {
		        printf("WARNING: Neighbor not found at %s line %d\n",
		    	   __FILE__, __LINE__);
		        continue;
		    }
            
		    d[0] = nbrNode->x - node->x;
		    d[1] = nbrNode->y - node->y;
		    d[2] = nbrNode->z - node->z;
            
		    ZImage(param, &d[0], &d[1], &d[2]);
		    mag2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            
/*          
 *          If the segment is zero length (which can happen when
 *          the mobility function is being called from SplitMultiNodes())
 *          just skip the segment.
 */         
		    if (mag2 < eps) continue;
            
		    mag     = sqrt(mag2);
/*
 *          Fi = Sum_j l_ij*(sigma . b_ij) x tij = -Fdrag = Sum_j l_ij*b_ij*(Bij(v).v + fb) 
 */
		    halfMag = mag/2.0;
		    invMag  = 1.0 / mag;
            
		    d[0] *= invMag;
		    d[1] *= invMag;
		    d[2] *= invMag;
/*          
 *          If the node has a very short segment, we need to flag it
 *          for later
 */         
            if (mag < shortSegCutoff) {
                hasShortSeg = 1;
                shortSegRatio = MIN(shortSegRatio, mag/param->minSeg);
            }
            
/*          
 *          Calculate how close to screw the arm is
 */         
            costheta = (d[0]*b[0] + d[1]*b[1] + d[2]*b[2]);
            costheta2 = (costheta*costheta);
            sintheta2 = 1.0 - costheta2;
            if(sintheta2<0)sintheta2=0.0;
            
/*          
 *          arms not on [1 1 1] planes don't move as readily as
 *          other arms, so must be handled specially.
 *          
 *          If needed, rotate a copy of the glide plane vector from the
 *          laboratory frame to the crystal frame.
 */         
            VectorCopy(n, node->nx[i], node->ny[i], node->nz[i]);
            cross(n, d, m);
#ifdef DEBUG
            if(!Numerical(n) || !Numerical(d) || !Numerical(m)){
                PrintNode(home, node);
                Fatal("%s: at line %d, i=%d",__func__,__LINE__,i);
            }
#endif
            tdKink = (node->onGB > -1 && nbrNode->onGB > -1 && fabs(DotProduct(b, param->TBNormal))>5E-2);

            if((DotProduct(m,m) < 1.0E-2) || ((GetPlaneType_FCC(home, node, nbrNode)& GLIDE_PLANE)==0) || tdKink==true){
	            Matrix33Init(Bi);
                if (numNbrs == 2) {
/*
 *                  junction nonphysical node: B = Beclimb * I2
 */
                    Bi[0][0] = halfMag * Beclimb;
                    Bi[1][1] = halfMag * Beclimb;
                    Bi[2][2] = halfMag * Beclimb;
                    VECTOR_ZERO(bForce);
                } else {

/*
 *                  line velocity along line sense
 */                 
                    if(DotProduct(nForce,d)<0){d[0]=-d[0];d[1]=-d[1];d[2]=-d[2];}
                    ParallelComponent(vTest, d, vgTest); vgMag=Sign(DotProduct(vgTest,d))*Normal(vgTest);
                    if(!partialBurg){
                        bs = B(ks/bMag,vs_max,vs_max_c,vgMag);
                        sbForceMag = FBack(ks/bMag,vs_max,vs_max_c,bs,vgMag);
                        
                        be = B(ke/bMag,ve_max,ve_max_c,vgMag);
                        ebForceMag = FBack(ke/bMag,ve_max,ve_max_c,be,vgMag);
                    }else{
                        bs = B(kts,vts_max,vts_max_c,vgMag);
                        sbForceMag = FBack(kts,vts_max,vts_max_c,bs,vgMag);
                        
                        be = B(kte,vte_max,vte_max_c,vgMag);
                        ebForceMag = FBack(kte,vte_max,vte_max_c,be,vgMag);
                    }
                    Bline = bs*costheta2+be*sintheta2; bForceMag = sbForceMag*costheta2 + ebForceMag*sintheta2;
                    bForce[0] = halfMag*bForceMag*d[0];
                    bForce[1] = halfMag*bForceMag*d[1];
                    bForce[2] = halfMag*bForceMag*d[2];
                    
                    BlmBecl = Bline - Beclimb;

/*
 *                  junction physical node: B = Beclimb*(I2-lineTensor) + Bline*lineTensor
 */
                    Bi[0][0] = halfMag * (d[0]*d[0]*BlmBecl + Beclimb);
                    Bi[0][1] = halfMag * (d[0]*d[1]*BlmBecl);
                    Bi[0][2] = halfMag * (d[0]*d[2]*BlmBecl);
                    Bi[1][1] = halfMag * (d[1]*d[1]*BlmBecl + Beclimb);
                    Bi[1][2] = halfMag * (d[1]*d[2]*BlmBecl);
                    Bi[2][2] = halfMag * (d[2]*d[2]*BlmBecl + Beclimb);
                }
                Bi[1][0] = Bi[0][1];
                Bi[2][0] = Bi[0][2];
                Bi[2][1] = Bi[1][2];

                Btotal[0][0] += Bi[0][0];Btotal[0][1] += Bi[0][1];Btotal[0][2] += Bi[0][2];
                Btotal[1][0] += Bi[1][0];Btotal[1][1] += Bi[1][1];Btotal[1][2] += Bi[1][2];
                Btotal[2][0] += Bi[2][0];Btotal[2][1] += Bi[2][1];Btotal[2][2] += Bi[2][2];
            } else {
/*
 *              Assemby lineTensor, climbTensor and verify glide direction
 */
                Vec3TransposeAndMult(d, lineTensor);
                Vec3TransposeAndMult(n, climbTensor);
                NormalizeVec(m);

                if(DotProduct(nForce,m)<0){m[0]=-m[0];m[1]=-m[1];m[2]=-m[2];}
                if(DotProduct(nForce,d)<0){d[0]=-d[0];d[1]=-d[1];d[2]=-d[2];}
                Vec3TransposeAndMult(m, glideTensor);
                
/*
 *              glide velocity along glide direction
 */
                ParallelComponent(vTest, m, vgTest); vgMag=Sign(DotProduct(vgTest,m))*Normal(vgTest);
                if (i==0) v1g=vgMag; else v2g=vgMag;

                if(!partialBurg){
                    bs = B(ks/bMag,vs_max,vs_max_c,vgMag);
                    sbForceMag = FBack(ks/bMag,vs_max,vs_max_c,bs,vgMag);
                    
                    be = B(ke/bMag,ve_max,ve_max_c,vgMag);
                    ebForceMag = FBack(ke/bMag,ve_max,ve_max_c,be,vgMag);
                }else{
                    bs = B(kts,vts_max,vts_max_c,vgMag);
                    sbForceMag = FBack(kts,vts_max,vts_max_c,bs,vgMag);
                    
                    be = B(kte,vte_max,vte_max_c,vgMag);
                    ebForceMag = FBack(kte,vte_max,vte_max_c,be,vgMag);
                }
                Bglide = bs*costheta2+be*sintheta2; bForceMag = sbForceMag*costheta2 + ebForceMag*sintheta2;
                bForce[0] = halfMag*bForceMag*m[0];
                bForce[1] = halfMag*bForceMag*m[1];
                bForce[2] = halfMag*bForceMag*m[2];
#ifndef NO_LINE_MOBILITY
/*
 *              line velocity along line sense
 */
                ParallelComponent(vTest, d, vgTest); vgMag=Sign(DotProduct(vgTest,d))*Normal(vgTest);
                if (i==0) v1l=vgMag; else v2l=vgMag;
                if(!partialBurg){
                    bs = B(ks/bMag,vs_max,vs_max_c,vgMag);
                    sbForceMag = FBack(ks/bMag,vs_max,vs_max_c,bs,vgMag);
                    
                    be = B(ke/bMag,ve_max,ve_max_c,vgMag);
                    ebForceMag = FBack(ke/bMag,ve_max,ve_max_c,be,vgMag);
                }else{
                    bs = B(kts,vts_max,vts_max_c,vgMag);
                    sbForceMag = FBack(kts,vts_max,vts_max_c,bs,vgMag);
                    
                    be = B(kte,vte_max,vte_max_c,vgMag);
                    ebForceMag = FBack(kte,vte_max,vte_max_c,be,vgMag);
                }
                Bline = bs*costheta2+be*sintheta2; bForceMag = sbForceMag*costheta2 + ebForceMag*sintheta2;
                bForce[0] += halfMag*bForceMag*d[0];
                bForce[1] += halfMag*bForceMag*d[1];
                bForce[2] += halfMag*bForceMag*d[2];
#else
                Bline = Beclimb;
#endif
                for(j=0;j<3;j++){
                    for(k=0;k<3;k++){
                        Bi[j][k] = halfMag*(Bglide*glideTensor[j][k] + Beclimb*climbTensor[j][k] + Bline*lineTensor[j][k]);
                        Btotal[j][k]+=Bi[j][k];
                    }
                }
                
#ifdef DEBUG
                if(!Numerical(fdrag)){
                    FormatMatrix(Bi,"Bi");
                    FormatVector(vTest,"vTest");
                    FormatVector(bForce,"bForce");
                    Fatal("Normal(v)=%e",Normal(vTest)*param->burgMag);
                }
#endif
           }  /* End non-[0 0 1] arm */
           fdrag[0] += Bi[0][0]*vTest[0] + Bi[0][1]*vTest[1] + Bi[0][2]*vTest[2] + bForce[0];
           fdrag[1] += Bi[1][0]*vTest[0] + Bi[1][1]*vTest[1] + Bi[1][2]*vTest[2] + bForce[1];
           fdrag[2] += Bi[2][0]*vTest[0] + Bi[2][1]*vTest[1] + Bi[2][2]*vTest[2] + bForce[2];
       }  /* End loop over arms */
/*      
 *      At this point we should check if the matrix is invertable and
 *      if it isn't, find the eigen values and eigen vectors of the drag
 *      matrix, then invert the drag matrix keeping zero eigen values
 *      as zero.
 *      FIX ME!  For now, we're assuming the matrix is invertable.
 */

        if (Matrix33Invert(Btotal, invBtotal) == 0) {
           Fatal("Mobility_FCC_Nolinear: Cannot invert 3X3 matrix!");
        }

        Matrix33Vector3Multiply(invBtotal, fdrag, vDelta);
        fErr=Normal(fdrag);

        VectorMinusSelf(vTest, vDelta);
        vMag = Sign(DotProduct(nForce,vTest))*Normal(vTest);
        if(fabs(vMag)>2.0*v_max){
            VectorTimesSelf(vTest, fabs(v_max/vMag));
            break;
        }
    } 
#ifdef DEBUG
    if(iter >= maxStep){
        PrintNode(home,node);
        Fatal("%s: the velocity if the node can be converge in %d step\n", __func__, maxStep);
    }
#endif
    vMag = Normal(vTest);
    if(fabs(vMag)>v_max){
        VectorTimesSelf(vTest, fabs(v_max/vMag));
    }
    ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, vTest);
    node->vX = vTest[0];
    node->vY = vTest[1];
    node->vZ = vTest[2];

/*
 *  Sometimes nodes oscillate at high frequency causing timestep to drop.
 *  Usually, these are 2-nodes that have short segments attached.  So, if
 *  we find a 2-node that has a least 1 short segment, and the velocity
 *  of the node has changed direction since the previous cycle, slow
 *  the node down for this cycle.  Hopefully this will help dampen the
 *  oscillation and improve the timestep.
 */
    if ((numNbrs == 2) && (hasShortSeg)) {
        if (((node->vX * node->oldvX) < 0.0) ||
            ((node->vY * node->oldvY) < 0.0) ||
            ((node->vZ * node->oldvZ) < 0.0)) {
            node->vX *= shortSegRatio;
            node->vY *= shortSegRatio;
            node->vZ *= shortSegRatio;
        }
    }
	

#ifdef _FEM
/*
 *      The velocity of any surface node along the negative surface
 *      normal direction should be zero to prevent the node moving into
 *      the box.  Make a call to adjsut the velocity accordingly.
 *
 *      Note: If the node moves outside the box, it's allowed here, then
 *      position is adjusted in AdjustNodePosition().
 */
//        AdjustSurfaceNodeVel(home, node);
#endif

#ifdef DEBUG
    if(!Numerical(node->vX, node->vY, node->vZ)){
        PrintNode(home, node);
        Fatal("%s at line %d",__func__,__LINE__);
    }
#endif
    return(0);
}
