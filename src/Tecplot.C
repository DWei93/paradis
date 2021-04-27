#include "Home.h"
#include "Util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include "Math.h"
#include "Force.h"

static void NodeofEle(Param_t *param, int eleNum, int node[8]){

	int	nx,ny,nz;

	nx = eleNum%param->meshNumX;
	ny = eleNum%(param->meshNumX*param->meshNumY)/param->meshNumX; 
	nz = eleNum/(param->meshNumX*param->meshNumY); 

	node[0] = nx + ny*(param->meshNumX+1) + nz*(param->meshNumX+1)*(param->meshNumY+1);
	node[1] = node[0] + 1;
	node[2] = node[0] + param->meshNumX+1;
	node[3] = node[2] + 1;
	node[4] = node[0] + (param->meshNumX+1)*(param->meshNumY+1);
	node[5] = node[4] + 1;
	node[6] = node[4] + param->meshNumX+1;
	node[7] = node[6] + 1;
	
	return;
	
}	

/*---------------------------------------------------------------------------
 *
 *      Function:    Tecplot
 *      Description: Write segment data in tecplot format... Masato
 *
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          numArms          count of all arms in the problem space
 *
 *-------------------------------------------------------------------------*/
void Tecplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
             int writePrologue, int writeEpilogue, int numArms, int debug)
{
        int      i, j, thisDomain, arm;
        int      newNodeKeyPtr, narm;
        int      burgID=0, normalID; 
        int	     armInGrain, armOnGB;
        real8    x,y,z, sfForce[3];
        real8    x2,y2,z2;
        real8    x2o,y2o,z2o;
        real8    prx,pry,prz;
        real8    vx,vy,vz;
        real8    bX, bY, bZ;
        real8    Lx, Ly, Lz;
        real8    normal[3], normalCryst[3];
        real8	 minDistance, pos[3];
        real8	 velocity, displacement, vecOldNew[3];
        real8    burg[3], burgMag;
        char     fileName[256];
        Node_t   *node, *nbrNode;
        Param_t  *param;
        FILE     *fp;
        Plane_t	 plane;
        struct stat statbuf;
        real8   al, am, an, amag, sigijk, pstnijk, dpstnijk, ptstnijk, bl, bm, bn;
        real8   vec[3], f[3], f1[3], f2[3], length, fMag,ifMag,n[3],val, nodePos[3], nbrPos[3];
        int     armID21;
        real8   tau_g, tau_gg, tau_ext;

        int         index, constraint, nodeOnGB;
        param      = home->param;
        thisDomain = home->myDomain;

        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;
		
/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s.plt",
                     DIR_TECPLOT, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_TECPLOT, baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(fileName, &statbuf);
#endif
#endif

	if(numArms == 0 && debug == 1){
            for (i = 0; i < home->newNodeKeyPtr; i++) {
                if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                    continue;
                }
		for(j=0; j<node->numNbrs; j++){
                    nbrNode = GetNeighborNode(home, node, j);
                    if (nbrNode == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n",
                               __FILE__, __LINE__);
                        continue;
                    }
		    numArms++;
		}
	    }

	    node = home->ghostNodeQ;
	    while (node) {
	        if(node == (Node_t *)NULL)continue;

	        for (arm = 0; arm < node->numNbrs; arm++) {
	            nbrNode = GetNeighborNode(home, node, arm);
	            if(nbrNode == (Node_t *)NULL){
	                continue;
	            }
	            numArms++;			
	        }
            node = node->next;
	    }

	    if(numArms > 0){
		firstInGroup = 1;
		writePrologue = 1;
	    }
	}

/*
 *      If this process is the first member of its I/O group
 *      it needs to create the output file and do any intializations
 *      needed.
 */
        if (firstInGroup) {
            if(param->loadType!=6){
                al = param->edotdir[0];
                am = param->edotdir[1];
                an = param->edotdir[2];
                
                amag = sqrt(al*al+am*am+an*an);
                
                al /= amag;
                am /= amag;
                an /= amag;
                bl=al; bm=am; bn=an;
            }else{
                al = param->loadDir[0];
                am = param->loadDir[1];
                an = param->loadDir[2];
                bl = param->loadDir[3];
                bm = param->loadDir[4];
                bn = param->loadDir[5];
            
            }
            sigijk = param->appliedStress[0]*al*bl        +
                     param->appliedStress[1]*am*bm        +
                     param->appliedStress[2]*an*bn        +
                     param->appliedStress[3]*(am*bn+an*bm)+
                     param->appliedStress[4]*(an*bl+al*bn)+ 
                     param->appliedStress[5]*(al*bm+bl*am);
            
            pstnijk =  param->totpStn[0]*al*bl        +
                       param->totpStn[1]*am*bm        +
                       param->totpStn[2]*an*bn        +
                       param->totpStn[3]*(am*bn+an*bm)+
                       param->totpStn[4]*(an*bl+al*bn)+
                       param->totpStn[5]*(al*bm+bl*am);
            
            ptstnijk = param->totTwinpStn[0]*al*bl        +
                       param->totTwinpStn[1]*am*bm        +
                       param->totTwinpStn[2]*an*bn        +
                       param->totTwinpStn[3]*(am*bn+an*bm)+
                       param->totTwinpStn[4]*(an*bl+al*bn)+
                       param->totTwinpStn[5]*(al*bm+bl*am);
/*
 *          First task in the I/O group must open the data file for writing
 *          to overwrite any existing file of the same name.
 */
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
            }

            if (writePrologue) {
                static int print = 0;
                print++;
                if(print % 100 == 0){
                    printf(" +++ Writing Tecplot file(s) %s\n", baseFileName);
                }
#if CHECK_DOMAIN
		        fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,diplacement,numNbr"
                        ",inGrain,onGB, burgID,normalID,constraint,planeErr,domain,native\n");	
		        VECTOR_ZERO(plane.point);
#else
	        if(debug){
		        fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,numNbr,inGrain,onGB,burgID,normalID,constraint,index, dom, id, native\n");
		        fprintf(fp, "zone T = \"%d\" i = %d  F=POINT\n", home->cycle, numArms); 
	        }else{

#if CHECK_GB_TOPO
            static int num = 0, cycle = 0;
            if(cycle == home->cycle) {
                num++;
            }else{
                cycle = home->cycle;
                num = 0;
            }
		    fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,numNbr,inGrain,onGB,burgID,normalID,constraint,index, dom, id, fx, fy, fz, fMag, tau_g, ifMag\n");
		    fprintf(fp, "zone T = \"%d\" i = %d  F=POINT\n", home->cycle, numArms); 
            val = param->timeNow*1.0E9;
            val = fabs(param->eRate*param->timeNow*1.0E2);
            fprintf(fp, "AUXDATA strain = \"%10.7f\"\n", val);
            val = fabs(sigijk/1.0E6);
            fprintf(fp, "AUXDATA stress = \"%10.7f\"\n", val);
            if(param->maxSchmid >= 0.0){
                fprintf(fp, "AUXDATA shear_stress = \"%10.7f\"\n", val*param->maxSchmid);
            }
            val = fabs(pstnijk*1.0E2);
            fprintf(fp, "AUXDATA plastic_strain = \"%10.7f\"\n", val);
            fprintf(fp, "AUXDATA fileName = \"%s\"\n", baseFileName);
#else
		    fprintf(fp, "variables = X,Y,Z,V1,V2,V3,numNbr,inGrain,onGB,burgID,normalID,constraint,index, fx, fy, fz, fMag, tau_g, tau_sf, tau_ext,ifMag\n");
            val = param->timeNow*1.0E9;
		    fprintf(fp, "zone T = \"%d\", i = %d, SOLUTIONTIME = %d,  F=POINT\n", home->cycle, numArms, home->cycle); 
            val = fabs((param->loadType!=4)?param->eRate*param->timeNow*1.0E2:param->netCyclicStrain*1.0E2);
            fprintf(fp, "AUXDATA strain = \"%.2f\"\n", val);
            val = (param->loadType!=4)?fabs(sigijk/1.0E6):sigijk/1.0E6;
            fprintf(fp, "AUXDATA stress = \"%.2f\"\n", val);
            if(param->maxSchmid >= 0.0){
                fprintf(fp, "AUXDATA shear_stress = \"%.2f\"\n", val*param->maxSchmid);
            }
            val = (param->loadType!=4)?fabs(pstnijk*1.0E2):pstnijk*1.0E2;
            fprintf(fp, "AUXDATA plastic_strain = \"%7.3f\"\n", val);
            fprintf(fp, "AUXDATA twinning_deform = \"%7.3f\"\n", fabs(ptstnijk*100.0));
            fprintf(fp, "AUXDATA fileName = \"%s\"\n", baseFileName);
#endif
	    }
#endif	
				
       	    }
        } else {
/*
 *          Any task NOT first in its I/O group must open the data file
 *          in an append mode so everything gets added to the end of
 *          the file.
 */
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
            }
        }
 

/*
 *      Generate the plot data for the segments associated with this
 *      domain's data
 */
        newNodeKeyPtr = home->newNodeKeyPtr;
        real8 glideConstraints[3][3], nodeLen;
        int numGlideConstraints;
        bool zeroLenJunc;

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {;
                continue;
            }
        
            x = node->x;
            y = node->y;
            z = node->z;
			
            vecOldNew[0] = x - node->oldx;
            vecOldNew[1] = y - node->oldy;
            vecOldNew[2] = z - node->oldz;
            	
            ZImage(param, vecOldNew, vecOldNew+1, vecOldNew+2);
            displacement = Normal(vecOldNew);

	        GetGlideConstraintList(home, node, &numGlideConstraints, glideConstraints, zeroLenJunc);
            f[0] = node->fX; f[1]=node->fY; f[2]=node->fZ;
            ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, f);
            NodeImageForce(home, node,  f1);
            ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, f1);
            nodeLen =  NodeLengthXMag(home, node);
            f[0]/=(nodeLen*1E6); f[1]/=(nodeLen*1E6); f[2]/=(nodeLen*1E6);
            f1[0]/=(nodeLen*1E6); f1[1]/=(nodeLen*1E6); f1[2]/=(nodeLen*1E6);
            fMag=Normal(f); ifMag=Normal(f1);

            NormalizeVec(f1); NormalizeVec(f);
            if(!Numerical(f)){VECTOR_ZERO(f);}
            if(!Numerical(f1)){VECTOR_ZERO(f1);}
            if(!Numerical(ifMag))ifMag=0.0;
            if(!Numerical(fMag))fMag=0.0;
          
/*
 *          Don't assume PBC is enabled...
 */
            if (param->xBoundType == Periodic) {
                x = x - Lx*rint(x/Lx);
            }

            if (param->yBoundType == Periodic) {
                y = y - Ly*rint(y/Ly);
            }

            if (param->zBoundType == Periodic) {
                z = z - Lz*rint(z/Lz);
            }
        
            constraint = node->constraint;

            if(numGlideConstraints==3){
                constraint = ((node->onGB>-1)?7:4);
            }else if(numGlideConstraints==2){
                if(node->numNbrs==1)constraint=1; // surface node
                else if (node->numNbrs==2){
                    if(fabs(dot_product(node->nx[0],node->ny[0],node->nz[0], node->nx[1], node->ny[1], node->nz[1]))>0.98)constraint = 0; 
                    else constraint = (node->onGB>-1?5:2);
                }else{
                    if(node->onGB>-1){
                        if((NodeOnGB(home, node) & GB_NODE)>0)constraint=6; else constraint=5;
                    }else{
                        constraint = 3;
                    }
                }
            }else constraint=0;

            for (j = 0; j < node->numNbrs; j++) {
#if 0
/*
 *              This check should be obsolete; should never
 *              have index < zero at this point in execution
 */
                if ((index = node->nbrTag[j].index) < 0) {
                    continue; 
                }
#endif
        
/*
 *  FIX ME? This will result in segments crossing domain boundaries
 *  to be added to the file by both domains!  Is this what is wanted?
 */
				
/*
                if ((node->nbrTag[j].domainID == thisDomain) && 
                    (node->nbrTag[j].index < i)) {
                    continue;
                }
 */
                nbrNode = GetNeighborNode(home, node, j);
				
                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if(home->param->materialType == MAT_TYPE_FCC){
                    burgID = GetBurgID_FCC(home, node, nbrNode);
                }else{
                    Fatal("The burgers ID of non-FCC material has not defined\n");
                }
				
                normalID = GetPlaneID_FCC(home, node, nbrNode);

                x2o = nbrNode->x;
                y2o = nbrNode->y;
                z2o = nbrNode->z;
        
                prx = x2o-x;
                pry = y2o-y;
                prz = z2o-z;
        
                if (param->xBoundType == Periodic) {
                    prx = prx - Lx*rint(prx/Lx);
                }

                if (param->yBoundType == Periodic) {
                    pry = pry - Ly*rint(pry/Ly);
                }

                if (param->zBoundType == Periodic) {
                    prz = prz - Lz*rint(prz/Lz);
                }
        
                
                x2 = x + prx;
                y2 = y + pry;
                z2 = z + prz;

                vx =  node->vX * param->burgMag;
                vy =  node->vY * param->burgMag;
                vz =  node->vZ * param->burgMag;
	            velocity = sqrt(vx*vx + vy*vy + vz*vz);
                length = sqrt(prx*prx+pry*pry+prz*prz); 
                narm = node->numNbrs;
				
                GetArmGrain(home, node, nbrNode, armInGrain, armOnGB);
                vec[0]=prx; vec[1]=pry; vec[2]=prz;
                NormalizeVec(vec);
                n[0] = node->nx[j]; n[1]=node->ny[j]; n[2]=node->nz[j];

                armID21=GetArmID(home, nbrNode, node);
                if(armID21<0)Fatal("%s: armID 21 < 0",__func__);
                f1[0] =(node->armfx[j]+nbrNode->armfx[armID21])/1E6/length;
                f1[1] =(node->armfy[j]+nbrNode->armfy[armID21])/1E6/length;
                f1[2] =(node->armfz[j]+nbrNode->armfz[armID21])/1E6/length;
                
                if((!Numerical(f1))){VECTOR_ZERO(f1);}
                NormalComponent(f1,n,f1);
                tau_g=Normal(f1);
                if(!Numerical(tau_g))tau_g=0;

                NormalComponent(f1,vec,f1);
                tau_gg=Normal(f1);
                if(!Numerical(tau_gg))tau_gg=0;

                burg[0] = node->burgX[j];
                burg[1] = node->burgY[j];
                burg[2] = node->burgZ[j];
                nodePos[0] = x; nodePos[1]=y; nodePos[2]=z;
                nbrPos[0] = x2; nbrPos[1]=y2; nbrPos[2]=z2;
                tau_ext = ExtShearStress(param, burg, nodePos, nbrPos, n)/1E6;
                if(!Numerical(tau_ext))tau_ext=0.0;
#if 0
                SelfForce(0, 1, 0.324,burg[0],burg[1],burg[2],x,y,z,x2,y2,z2,param->rc,param->Ecore,f,f1);
                fMag=Normal(f);
#endif

                if(param->polycrystal==0){
                    index = normalID;
                }else{
                    if(AlignedVector(param->TBNormal, n)){
                        if(burgID>6 && burgID<22)index=4; else index=3;
                    }else{
                        vec[0]=node->burgX[j]; vec[1]=node->burgY[j]; vec[2]=node->burgZ[j];
                        if(PerpendicularDirection(param->TBNormal, vec)) index =2;
                        else index = 1;
                    }
                }
                if((GetPlaneType_FCC(home, node, nbrNode) & GLIDE_PLANE) == 0){
                    burgMag=Normal(burg);
                    if(burgID < 6) index=5;
                    else if (fabs(1.41421356-burgMag)<0.05)index=6;
                    else if(fabs(0.816497-burgMag)<0.05) index = 7; //stair rod
                    else if(fabs(1.63299-burgMag)<0.05) index = 8; //stair rod
                    else index = 9;

                }
/*
 *              This part is for output of junction investigation
 */
                StackingFaultForce(home, prx, pry, prz, param->sfEnergy/param->burgMag, node->signSF[j], n, f1, f2);
                tau_gg=Normal(f1)*2.0E-6/length; NormalizeVec(f1); 
                if(!Numerical(tau_gg))tau_gg=0.0;
				
#if CHECK_DOMAIN
                fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %7.4f %d %d %d %d %d %d %d %e %d %d\n", 
                        x,y,z,prx/2,pry/2,prz/2,velocity,displacement,narm,armInGrain,armOnGB,burgID,normalID,
                        node->constraint,minDistance,node->myTag.domainID,1);
#else
	            if(debug){
                    fprintf(fp, "%15.10f %15.10f %15.10f %6.1f %6.1f %6.1f %7.4f %d %d %d %d %d %d %d %d %d 1\n", 
                            x,y,z,prx/2,pry/2,prz/2,velocity,narm,node->inGrain,node->onGB,burgID,normalID,constraint,
                            index, node->myTag.domainID, node->myTag.index);
	            }else{
#if CHECK_GB_TOPO
                    fprintf(fp, "%15.10f %15.10f %15.10f %6.1f %6.1f %6.1f %7.4f %d %d %d %d %d %d %d %d %d %g %g %g %g %g %g\n", 
                            x,y,z,prx/2,pry/2,prz/2,velocity,narm,node->inGrain,node->onGB,burgID,normalID,constraint,
                            index, node->myTag.domainID, node->myTag.index, f[0], f[1], f[2], fMag, tau_g, ifMag);
#else
                    fprintf(fp, "%g %g %g %g %g %g %d %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n", 
                            x,y,z,prx/2,pry/2,prz/2,
                            narm,node->inGrain,node->onGB,burgID,normalID,constraint,
                            index, f[0],f[1],f[2],fMag,tau_g,tau_gg,tau_ext,ifMag);
#endif
	        }
#endif	
			
                
            }
        }
		
		
	if(debug==1){
		
	    node = home->ghostNodeQ;
	    while (node) {
	        if(node == (Node_t *)NULL)continue;

	        x = node->x;
	        y = node->y;
	        z = node->z;
			
	        vecOldNew[0] = x - node->oldx;
	        vecOldNew[1] = y - node->oldy;
	        vecOldNew[2] = z - node->oldz;
	        	
	        ZImage(param, vecOldNew, vecOldNew+1, vecOldNew+2);
	        	
	        displacement = Normal(vecOldNew);
	        	
	        if (param->xBoundType == Periodic) {
                    x = x - Lx*rint(x/Lx);
                }
			
                if (param->yBoundType == Periodic) {
                    y = y - Ly*rint(y/Ly);
                }
			
	        if (param->zBoundType == Periodic) {
                    z = z - Lz*rint(z/Lz);
	        }

	        for (arm = 0; arm < node->numNbrs; arm++) {
	            nbrNode = GetNeighborNode(home, node, arm);
	            if(nbrNode == (Node_t *)NULL){
	                continue;
	            }
				
                normalID = GetPlaneID_FCC(home, node, nbrNode);

	            if(home->param->materialType == MAT_TYPE_FCC){
	                burgID = GetBurgID_FCC(home, node, nbrNode);
	            }else{
	                Fatal("The burgers ID of non-FCC material has not defined\n");
	            }
		
	            x2o = nbrNode->x;
	            y2o = nbrNode->y;
	            z2o = nbrNode->z;
	
	            prx = x2o-x;
	            pry = y2o-y;
	            prz = z2o-z;
				
	            if (param->xBoundType == Periodic) {
	                prx = prx - Lx*rint(prx/Lx);
	            }
	            	
	            if (param->yBoundType == Periodic) {
	                pry = pry - Ly*rint(pry/Ly);
	            }
	            	
	            if (param->zBoundType == Periodic) {
	                prz = prz - Lz*rint(prz/Lz);
	            }
				
                x2 = x + prx;
                y2 = y + pry;
                z2 = z + prz;
				
                vx =  node->vX * param->burgMag;
                vy =  node->vY * param->burgMag;
                vz =  node->vZ * param->burgMag;

	            velocity = sqrt(vx*vx + vy*vy + vz*vz);
				
	            narm = node->numNbrs;

                GetArmGrain(home, node, nbrNode, armInGrain, armOnGB);
                GetArmGrain(home, node, nbrNode, armInGrain, armOnGB);
        
                if(armOnGB<0){
                    if(armInGrain == 0){
                        index = 3;
                    }else{
                        index = 6;
                    }
                }else{
                    if(burgID == 7){
                        index = 0;
                    }else if(burgID == 12){
                        index = 0;
                    }else if(burgID == 17){
                        index =0;
                    }else{
                        burg[0] = node->burgX[j];
                        burg[1] = node->burgY[j];
                        burg[2] = node->burgZ[j];
                        if(BurgType_FCC(burg) == BASIC_BURG){
                            index = 1;
                        }else{
                            index = 2;
                        }
                    }
                }

                constraint = node->constraint;

	            GetGlideConstraintList(home, node, &numGlideConstraints, glideConstraints, zeroLenJunc);
                f[0] = node->fX; f[1]=node->fY; f[2]=node->fZ;
                ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, f);
                NodeImageForce(home, node,  f1);
                ApplyConstraintsToVelocity(numGlideConstraints, glideConstraints, f1);
                nodeLen =  NodeLengthXMag(home, node);
                f[0]/=(nodeLen*1E6); f[1]/=(nodeLen*1E6); f[2]/=(nodeLen*1E6);
                f1[0]/=(nodeLen*1E6); f1[1]/=(nodeLen*1E6); f1[2]/=(nodeLen*1E6);
                fMag=Normal(f); ifMag=Normal(f1);
                
                NormalizeVec(f1); NormalizeVec(f);
                if(!Numerical(f)){VECTOR_ZERO(f);}
                if(!Numerical(f1)){VECTOR_ZERO(f1);}
                if(!Numerical(ifMag))ifMag=0.0;
                if(!Numerical(fMag))fMag=0.0;

                if((NodeOnGB(home, node) & INTERSECT_NODE)>0){
                    if(narm==2)constraint = 4; else if(narm==3) constraint = 5; else constraint =6;
                }

                fprintf(fp, "%15.10f %15.10f %15.10f %6.1f %6.1f %6.1f %7.4f %d %d %d %d %d %d %d %d %d 0 %g %g %g %g %g %g\n", 
                        x,y,z,prx/2,pry/2,prz/2,velocity,narm,node->inGrain,node->onGB,burgID,normalID,constraint,
                        index, node->myTag.domainID, node->myTag.index, f[0], f[1], f[2], fMag, tau_g, ifMag);
	        }
			
            node = node->next;
	    }
	}
/*
 *      Handle anything that needs to be done after all data has been
 *      written to the file
 */
        if (writeEpilogue) {
            fprintf(fp, "\n");
        }
      
        fclose(fp);
       
	return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    TecplotMesh
 *      Description: Write segment data in tecplot format with mesh. This is 
 *      	     used to output the plastic strain and dislocation density
 *      	     contour.
 *
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *          numSegs          count of all segments in the problem space
 *
 *-------------------------------------------------------------------------*/
void TecplotMesh(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
             int writePrologue, int writeEpilogue, int numSegs)
{
	int      i, j;
	int      newNodeKeyPtr;
	int	     eleNum,iBurg;
	real8    prx,pry,prz,burg[3];
	real8    bX, bY, bZ;
	real8    Lx, Ly, Lz;
	real8	 midPos[3],deltDen;
	real8	 burgMag;
	char     fileNameEle[256];
	Node_t   *node, *nbrNode;
	Param_t  *param;
	FILE     *fpEle;

	int		nodeEle[8],nx,ny,nz;
	int		imageNode[8],imageNum;
    int     meshNodeNumX, meshNodeNumY, meshNodeNumZ, totNodeNum;
    int     totEleNum;
	real8 	**nodeContour;
	real8	eleLx,eleLy,eleLz;
	real8	effPlaStn,currDisDen;

	real8	**totPlaStn;
	real8	**totDisDen;

    real8     *locVals, *gloVals;

	param      = home->param;

	Lx = param->Lx;
	Ly = param->Ly;
	Lz = param->Lz;

    meshNodeNumX = param->meshNumX + 1;
    meshNodeNumY = param->meshNumY + 1;
    meshNodeNumZ = param->meshNumZ + 1;
    totNodeNum = meshNodeNumX*meshNodeNumY*meshNodeNumZ;

    totEleNum = param->meshNumX*param->meshNumY*param->meshNumZ;
    nodeContour = (real8 **)malloc(totNodeNum*sizeof(real8*));
    for(i=0; i<totNodeNum; i++){
        nodeContour[i] = (real8*)malloc(3*sizeof(real8));
    }

    totPlaStn = (real8 **)malloc(totEleNum*sizeof(real8*));        
    totDisDen = (real8 **)malloc(totEleNum*sizeof(real8*));        

    locVals = (real8 *)calloc(totEleNum*6, sizeof(real8));
    gloVals = (real8 *)calloc(totEleNum*6, sizeof(real8));

    for(i=0; i<totEleNum*6; i++){
        locVals[i] = 0.0;
        gloVals[i] = 0.0;
    }
    
    for(i=0; i<totEleNum; i++){
        totPlaStn[i] = (real8 *)calloc(6,sizeof(real8));
        totDisDen[i] = (real8 *)calloc(6,sizeof(real8));
        for(j=0; j<6; j++){
            totPlaStn[i][j] = 0.0;
            totDisDen[i][j] = 0.0;
        }
    }

	eleLx = param->meshSizeX;
	eleLy = param->meshSizeY;
	eleLz = param->meshSizeZ;

	/*	Zero the arrays of dislocation density in each element.	*/

	for(i=0; i<totEleNum; i++){
	    home->disDenEle[i][0] = 0.0;
	    home->disDenEle[i][1] = 0.0;
	    home->disDenEle[i][2] = 0.0;
	    home->disDenEle[i][3] = 0.0;
	    home->disDenEle[i][4] = 0.0;
	    home->disDenEle[i][5] = 0.0;
	}
	
	/*
	 *      Set data file name.  Only append a sequence number to
	 *      the data file name if the data is to be spread across
	 *      multiple files.
	 */
	if (firstInGroup) {
		snprintf(fileNameEle, sizeof(fileNameEle), "%s/%s.plt",
				 DIR_TECMESH, baseFileName);
		/*
		 *          First task in the I/O group must open the data file for writing
		 *          to overwrite any existing file of the same name.
		 */

		if ((fpEle = fopen(fileNameEle, "w")) == (FILE *)NULL) {
			Fatal("tec_plot_mesh: Open error %d on %s\n", errno, fileNameEle);
		}
	}

	/*
	 *      Generate the plot data for the segments associated with this
	 *      domain's data
	 */
	newNodeKeyPtr = home->newNodeKeyPtr;

	for (i = 0; i < newNodeKeyPtr; i++) {

		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
		
		for (j = 0; j < node->numNbrs; j++) {

			nbrNode = GetNeighborNode(home, node, j);

			if (nbrNode == (Node_t *)NULL) continue;
			
			if (OrderNodes(node, nbrNode) != -1) continue;

			prx = nbrNode->x - node->x;
			pry = nbrNode->y - node->y;
			prz = nbrNode->z - node->z;
			
			if (param->xBoundType == Periodic) prx = prx - Lx*rint(prx/Lx);
			if (param->yBoundType == Periodic) pry = pry - Ly*rint(pry/Ly);
			if (param->zBoundType == Periodic) prz = prz - Lz*rint(prz/Lz);
			
			midPos[0] = node->x + prx*0.5;
			midPos[1] = node->y + pry*0.5;
			midPos[2] = node->z + prz*0.5;

			FoldBox(param,&midPos[0],&midPos[1],&midPos[2]);
			
			eleNum = Coor2EleNum(home, midPos);

			deltDen = sqrt(prx*prx + pry*pry + prz*prz)*param->burgVolFactor;

			bX = node->burgX[j];
			bY = node->burgY[j];
			bZ = node->burgZ[j];

			burgMag = sqrt(bX*bX + bY*bY +bZ*bZ);

			burg[0] = bX;
			burg[1] = bY;
			burg[2] = bZ;
			
			iBurg = BurgType_FCC(burg);

			home->disDenEle[eleNum][iBurg] += deltDen;
			
		}
	}

	/*	
	 *	Following is the output of dislocation density and plastic strain
	 *	in each element.
	 */
#ifdef PARALLEL

    for (i=0; i<totEleNum; i++){
        for(j = 0; j<6; j++){
            locVals[i*6+j] = home->disDenEle[i][j];
        }
    }

	MPI_Reduce(locVals, gloVals, 6*totEleNum, MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);

    for (i=0; i<totEleNum; i++){
        for(j = 0; j<6; j++){
            home->disDenEle[i][j] = gloVals[i*6+j];
        }
    }
#endif

	if (firstInGroup) {
		for(i=0;i<totEleNum;i++){
		    totPlaStn[i][0] = home->pstnEle[i][0];
		    totPlaStn[i][1] = home->pstnEle[i][1];
		    totPlaStn[i][2] = home->pstnEle[i][2];
		    totPlaStn[i][3] = home->pstnEle[i][3];
		    totPlaStn[i][4] = home->pstnEle[i][4];
		    totPlaStn[i][5] = home->pstnEle[i][5];
    
		    totDisDen[i][0] = home->disDenEle[i][0];
		    totDisDen[i][1] = home->disDenEle[i][1];
		    totDisDen[i][2] = home->disDenEle[i][2];
		    totDisDen[i][3] = home->disDenEle[i][3];
		    totDisDen[i][4] = home->disDenEle[i][4];
		    totDisDen[i][5] = home->disDenEle[i][5];
		}

	    for(i=0;i<totNodeNum;i++){
		    VECTOR_ZERO(nodeContour[i]);
	    }
	
//	    fprintf(fpEle,"solutiontime = %e\n", param->timeNow);
	    fprintf(fpEle, "Variables = \"X\", \"Y\", \"Z\",\"plaStn"
				"\",\"disDen\"\n"); 
	    fprintf(fpEle, "Zone T = \"%d\", SOLUTIONTIME = %d, N= %d  E= %d F=FEPOINT ET=Brick\n",
				home->cycle, home->cycle,totNodeNum,totEleNum); 

	    for(i=0;i<totEleNum;i++){

			effPlaStn=sqrt((totPlaStn[i][0]*totPlaStn[i][0] +
							totPlaStn[i][1]*totPlaStn[i][1] +
							totPlaStn[i][2]*totPlaStn[i][2] +
							2*totPlaStn[i][3]*totPlaStn[i][3] +
							2*totPlaStn[i][4]*totPlaStn[i][4] +
							2*totPlaStn[i][5]*totPlaStn[i][5])*2/3);

			currDisDen = totDisDen[i][0] + totDisDen[i][1] + 
						 totDisDen[i][2] + totDisDen[i][3] +
						 totDisDen[i][4] + totDisDen[i][5]; 

			NodeofEle(param, i,nodeEle);	
			
			for(j=0;j<8;j++){

				nodeContour[nodeEle[j]][0] += effPlaStn;
				nodeContour[nodeEle[j]][1] += currDisDen;
				nodeContour[nodeEle[j]][2] ++;

			}
			
	    }

	    for(i=0;i<totNodeNum;i++){

			imageNode[0] = i;		
			
			imageNum = 1;

			nx = i%meshNodeNumX;
			ny = i%(meshNodeNumX*meshNodeNumY)/meshNodeNumX; 
			nz = i/(meshNodeNumX*meshNodeNumY); 
			
			if(nodeContour[i][2] > 8){

				fprintf(fpEle,"%10.3e %10.2e %10.2e %10.2e %10.2e\n",
						nx*eleLx-Lx/2, ny*eleLy-Ly/2,nz*eleLz-Lz/2,
						nodeContour[i][0],nodeContour[i][1]);
				continue;
			}
			/*
			 * 		For surface with periodic boundary condition.
			 */
			if((nx==0||nx==meshNodeNumX)&&param->xBoundType == Periodic){
				
				imageNode[1] = i+ (nx==0?1:-1)*(meshNodeNumX-1);
				imageNum *= 2;
			}

			if((ny==0||ny==meshNodeNumY)&&param->yBoundType == Periodic){
				
				for(j=0;j<imageNum;j++){
					imageNode[imageNum+j] = imageNode[j] + (ny==0?1:-1)*
											meshNodeNumX*(meshNodeNumZ-1);
				}
				imageNum *= 2;
			}

			if((nz==0||nz==meshNodeNumZ)&&param->zBoundType == Periodic){
				
				for(j=0;j<imageNum;j++){
					imageNode[imageNum+j] = imageNode[j] + (nz==0?1:-1)*
											meshNodeNumX*meshNodeNumY*(meshNodeNumZ-1);
				}
				imageNum *= 2;
			}
			
			for(j=1;j<imageNum;j++){

				nodeContour[i][0] += nodeContour[imageNode[j]][0];
				nodeContour[i][1] += nodeContour[imageNode[j]][1];
			}
			
			nodeContour[i][2] *= imageNum;
			
			nodeContour[i][0] /= nodeContour[i][2];
			nodeContour[i][1] /= nodeContour[i][2];

			if(imageNum != 1) nodeContour[i][2] = 9;

			for(j=1;j<imageNum;j++){

				VECTOR_COPY(nodeContour[imageNode[j]], nodeContour[i]);
			}
			
			fprintf(fpEle,"%10.2e %10.2e %10.2e %10.2e %10.2e\n",
					nx*eleLx-Lx/2,ny*eleLy-Ly/2,nz*eleLz-Lz/2,
					nodeContour[i][0],nodeContour[i][1]);

	    }

	    for(i=0;i<totEleNum;i++){

			NodeofEle(param, i,nodeEle);	
			
			fprintf(fpEle, "%8d %8d %8d %8d %8d %8d %8d %8d\n", nodeEle[0]+1,
					nodeEle[1]+1,nodeEle[3]+1,nodeEle[2]+1,nodeEle[4]+1,
					nodeEle[5]+1,nodeEle[7]+1,nodeEle[6]+1);
	    }
        

	    fclose(fpEle);
	}

    for(i=0; i<totNodeNum; i++){
        free(nodeContour[i]);
    }
    free(nodeContour);
    nodeContour = NULL;

    for(i=0; i<totEleNum; i++){
        free(totPlaStn[i]);
        free(totDisDen[i]);
    }

    free(totPlaStn);
    free(totDisDen);
    totPlaStn = NULL;
    totDisDen = NULL;

    free(locVals);
    free(gloVals);

	return;
}
