#include "Home.h"
#include "Util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>

static void NodeofEle(int eleNum, int node[8]){

	int	nx,ny,nz;

	nx = eleNum%MESH_NUM;
	ny = eleNum%(MESH_NUM*MESH_NUM)/MESH_NUM; 
	nz = eleNum/(MESH_NUM*MESH_NUM); 

	node[0] = nx + ny*(MESH_NUM1) + nz*(MESH_NUM1)*(MESH_NUM1);
	node[1] = node[0] + 1;
	node[2] = node[0] + MESH_NUM1;
	node[3] = node[2] + 1;
	node[4] = node[0] + (MESH_NUM1)*(MESH_NUM1);
	node[5] = node[4] + 1;
	node[6] = node[4] + MESH_NUM1;
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
        int      burgID, normalID; 
        real8    x,y,z;
        real8    x2,y2,z2;
        real8    x2o,y2o,z2o;
        real8    prx,pry,prz;
        real8    vx,vy,vz;
        real8    bX, bY, bZ;
        real8    Lx, Ly, Lz;
        real8     normal[3], normalCryst[3];
        real8	 armInGrain;
        real8	 minDistance, pos[3];
        real8	 velocity, displacement, vecOldNew[3];
        char     fileName[256];
        Node_t   *node, *nbrNode;
        Param_t  *param;
        FILE     *fp;
        Plane_t	 plane;
        struct stat statbuf;


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
		fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,diplacement,numNbr,inGrain,burgID,normalID,constraint,planeErr,domain,native\n");	
		VECTOR_ZERO(plane.point);
#else
	        if(debug){
		    fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,numNbr,inGrain,burgID,normalID,constraint,displace,domain,index\n");
		    fprintf(fp, "zone T = \"%06d\" i = %d  F=POINT\n", home->cycle,numArms); 
	        }else{
		    fprintf(fp, "variables = X,Y,Z,V1,V2,V3,velocity,numNbr,inGrain,burgID,normalID,constraint\n");
		    fprintf(fp, "zone T = \"%10.5f\" i = %d  F=POINT\n", param->timeNow*1.0E9,numArms); 
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

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {;
                continue;
            }
        
            x = node->x;
            y = node->y;
            z = node->z;
			
            vecOldNew[X] = x - node->oldx;
            vecOldNew[Y] = y - node->oldy;
            vecOldNew[Z] = z - node->oldz;
            	
            ZImage(param, vecOldNew, vecOldNew+1, vecOldNew+2);
            	
            displacement = Normal(vecOldNew);
          
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
				
                normal[X] = node->nx[j];
                normal[Y] = node->ny[j];
                normal[Z] = node->nz[j];

                if(home->param->useLabFrame){
	            ChangeOrientation(home, node, nbrNode,3);
                    Matrix33Vector3Multiply(home->rotMatrixInverse, normal, normalCryst);
                }else{
                    VECTOR_COPY(normalCryst, normal);
                }

                normalID = GlidePlaneID_FCC(normalCryst);

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
				
                armInGrain = ChangeOrientation(home, node, nbrNode, 0);
        
/*
 *              This part is for output of junction investigation
 */
				
#if CHECK_DOMAIN
                fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %7.4f %d %3.1f %d %d %d %d %e %d %d\n", 
                        x,y,z,prx/2,pry/2,prz/2,velocity,displacement,narm,armInGrain,burgID,normalID,
                        node->constraint,minDistance,node->myTag.domainID,1);
#else
	        if(debug){
                    fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %d %3.1f %d %d %d %7.4f %d %d\n", 
                            x,y,z,prx/2,pry/2,prz/2,velocity,narm,node->inGrain,burgID,normalID,node->constraint,
	                    displacement,node->myTag.domainID,node->myTag.index);
	        }else{
                    fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %d %3.1f %d %d %d\n", 
                            x,y,z,prx/2,pry/2,prz/2,velocity,narm,node->inGrain,burgID,normalID,node->constraint);
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
			
	        vecOldNew[X] = x - node->oldx;
	        vecOldNew[Y] = y - node->oldy;
	        vecOldNew[Z] = z - node->oldz;
	        	
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
				
                    normal[X] = node->nx[arm];
                    normal[Y] = node->ny[arm];
                    normal[Z] = node->nz[arm];
    
                    if(home->param->useLabFrame){
	                ChangeOrientation(home, node, nbrNode,3);
                        Matrix33Vector3Multiply(home->rotMatrixInverse, normal, normalCryst);
                    }else{
                        VECTOR_COPY(normalCryst, normal);
                    }

                    normalID = GlidePlaneID_FCC(normalCryst);

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

	            armInGrain = ChangeOrientation(home, node, nbrNode, 0);

                    fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %d %3.1f %d %d %d %7.4f %d %d\n", 
                            x,y,z,prx/2,pry/2,prz/2,velocity,narm,node->inGrain,burgID,normalID,node->constraint,
	                    displacement,node->myTag.domainID,node->myTag.index);
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
	int	 eleNum,iBurg;
	real8    prx,pry,prz,burg[3];
	real8    bX, bY, bZ;
	real8    Lx, Ly, Lz;
	real8	 midPos[3],deltDen;
	real8	 burgMag;
	char     fileNameEle[256];
	Node_t   *node, *nbrNode;
	Param_t  *param;
	FILE     *fpEle;

	real8	totPlaStn[TOT_ELE_NUM][6];
	real8	totDisDen[TOT_ELE_NUM][4];

	param      = home->param;

	Lx = param->Lx;
	Ly = param->Ly;
	Lz = param->Lz;

	/*	Zero the arrays of dislocation density in each element.	*/

	for(i=0; i<TOT_ELE_NUM; i++){
	    param->disDenEle[i][0] = 0.0;
	    param->disDenEle[i][1] = 0.0;
	    param->disDenEle[i][2] = 0.0;
	    param->disDenEle[i][3] = 0.0;
	}
	
	/*
	 *      Set data file name.  Only append a sequence number to
	 *      the data file name if the data is to be spread across
	 *      multiple files.
	 */
	if (firstInGroup) {
		snprintf(fileNameEle, sizeof(fileNameEle), "%s/mesh%s.plt",
				 DIR_TECPLOT, baseFileName);
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

			burg[X] = bX;
			burg[Y] = bY;
			burg[Z] = bZ;
			
			iBurg = BurgType_FCC(home,burg);

			param->disDenEle[eleNum][iBurg] += deltDen;
			
		}
	}

	/*	
	 *	Following is the output of dislocation density and plastic strain
	 *	in each element.
	 */
#ifdef PARALLEL

	MPI_Reduce(param->pstnEle, totPlaStn, 6*TOT_ELE_NUM, MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(param->disDenEle, totDisDen, 4*TOT_ELE_NUM, MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for(i=0;i<TOT_ELE_NUM;i++){
	    totPlaStn[i][0] = param->pstnEle[i][0];
	    totPlaStn[i][1] = param->pstnEle[i][1];
	    totPlaStn[i][2] = param->pstnEle[i][2];
	    totPlaStn[i][3] = param->pstnEle[i][3];
	    totPlaStn[i][4] = param->pstnEle[i][4];
	    totPlaStn[i][5] = param->pstnEle[i][5];

	    totDisDen[i][0] = param->disDenEle[i][0];
	    totDisDen[i][1] = param->disDenEle[i][1];
	    totDisDen[i][2] = param->disDenEle[i][2];
	    totDisDen[i][3] = param->disDenEle[i][3];
	}

#endif

	if (firstInGroup) {
		
	    int		nodeEle[8],nx,ny,nz;
	    int		imageNode[8],imageNum;
	    real8	nodeContour[TOT_NODE_NUM][3];
	    real8	eleLx,eleLy,eleLz;
	    real8	effPlaStn,currDisDen;

	    eleLx = Lx/MESH_NUM;
	    eleLy = Ly/MESH_NUM;
	    eleLz = Lz/MESH_NUM/ASPECT_RATIO;

	    for(i=0;i<TOT_NODE_NUM;i++){
			
		VECTOR_ZERO(nodeContour[i]);
	    }
	
//	    fprintf(fpEle,"solutiontime = %e\n", param->timeNow);
	    fprintf(fpEle, "Variables = \"X\", \"Y\", \"Z\",\"plaStn"
				"\",\"disDen\"\n"); 
	    fprintf(fpEle, "Zone T = \"%10.5f\" N= %d  E= %d F=FEPOINT ET=Brick\n",
				param->timeNow*1.0E9, TOT_NODE_NUM,TOT_ELE_NUM); 

	    for(i=0;i<TOT_ELE_NUM;i++){

			effPlaStn=((real8)TOT_ELE_NUM) * sqrt((totPlaStn[i][0]*totPlaStn[i][0] +
							totPlaStn[i][1]*totPlaStn[i][1] +
							totPlaStn[i][2]*totPlaStn[i][2] +
							2*totPlaStn[i][3]*totPlaStn[i][3] +
							2*totPlaStn[i][4]*totPlaStn[i][4] +
							2*totPlaStn[i][5]*totPlaStn[i][5])*2/3);

			currDisDen = totDisDen[i][0] + totDisDen[i][1] + 
						 totDisDen[i][2] + totDisDen[i][3]; 

			NodeofEle(i,nodeEle);	
			
			for(j=0;j<8;j++){

				nodeContour[nodeEle[j]][0] += effPlaStn;
				nodeContour[nodeEle[j]][1] += currDisDen;
				nodeContour[nodeEle[j]][2] ++;

			}
			
	    }

	    for(i=0;i<TOT_NODE_NUM;i++){

			imageNode[0] = i;		
			
			imageNum = 1;

			nx = i%MESH_NUM1;
			ny = i%(MESH_NUM1*MESH_NUM1)/MESH_NUM1; 
			nz = i/(MESH_NUM1*MESH_NUM1); 
			
			if(nodeContour[i][2] >8){

				fprintf(fpEle,"%10.3e %10.2e %10.2e %10.2e %10.2e\n",
						nx*eleLx-Lx/2, ny*eleLy-Ly/2,nz*eleLz-Lz/2,
						nodeContour[i][0],nodeContour[i][1]);
				continue;
			}
			/*
			 * 		For surface with periodic boundary condition.
			 */
			if((nx==0||nx==MESH_NUM1)&&param->xBoundType == Periodic){
				
				imageNode[1] = i+ (nx==0?1:-1)*MESH_NUM;
				imageNum *= 2;
			}

			if((ny==0||ny==MESH_NUM1)&&param->yBoundType == Periodic){
				
				for(j=0;j<imageNum;j++){
					imageNode[imageNum+j] = imageNode[j] + (ny==0?1:-1)*
											MESH_NUM1*MESH_NUM;
				}
				imageNum *= 2;
			}

			if((nz==0||nz==MESH_NUM1)&&param->zBoundType == Periodic){
				
				for(j=0;j<imageNum;j++){
					imageNode[imageNum+j] = imageNode[j] + (nz==0?1:-1)*
											MESH_NUM1*MESH_NUM1*MESH_NUM;
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

	    for(i=0;i<TOT_ELE_NUM;i++){

			NodeofEle(i,nodeEle);	
			
			fprintf(fpEle, "%8d %8d %8d %8d %8d %8d %8d %8d\n", nodeEle[0]+1,
					nodeEle[1]+1,nodeEle[3]+1,nodeEle[2]+1,nodeEle[4]+1,
					nodeEle[5]+1,nodeEle[7]+1,nodeEle[6]+1);
	    }

	    fclose(fpEle);
	}


	return;
}
