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
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Home.h" 
#include "Grain.h"
#include "Util.h" 
#include <vector>
#include <algorithm>
#include <numeric>

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>


bool IsOnGB(Node_t *node)
{
    return(node->onGB > -1);
}

real8 Normal(std::vector<real8>& a)
{
    return( sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) );
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
	    PrintNode(home, nbr);
	}
	return;
}

void FormatMatrix(real8 M[3][3], const char *msg){
	int i;
	printf("%s:\n{",msg);
	for (i = 0; i < 3; i++){
	    printf("{%20.15f,%20.15f,%20.15f}", M[i][0], M[i][1], M[i][2]);
        if(i<2){
            printf(",\n");
        }else{
            printf("}\n");
        }
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


/** Draw polygon with tecplot formart
 *  \param[in] fp file pointer for output
 *  \param[in] f_vert vertices index list
 *  \param[in] v coordinates of vertices
 *  \param[in] startID the start index in f_vert for output polygon
 *  \param[in] zoneName the zone name of tecplot
 *  \param[in] vec the variales of the polygon */
void DrawPolygon(FILE *fp, std::vector<int> &f_vert,std::vector<double> &v,int startID, char *zoneName, std::vector<int> &vec, bool startFromCentor=false){
    real8 centor[3];
    int n;
    if(startFromCentor){
        CentroidOfPolygon(f_vert,v,startID, centor);
    }
    n=f_vert[startID++];
    if(startFromCentor)fprintf(fp,"Zone T = \"%s\", N = %d, E = %d, F = FEPOINT, ET = TRIANGLE\n", zoneName, n+1, n);
    else fprintf(fp,"Zone T = \"%s\", N = %d, E = %d, F = FEPOINT, ET = TRIANGLE\n", zoneName, n, n-2);

    std::vector<std::vector<int> > record(2);
    record[0].resize(n);    record[1].resize(n);

    for(int i=startID; i<startID+n; i++){
        fprintf(fp, "%g %g %g", v[f_vert[i]*3], v[f_vert[i]*3+1], v[f_vert[i]*3+2]);
        for(int j=0; j<vec.size(); j++)fprintf(fp," %d", vec[j]); 
        fprintf(fp, "\n");
        record[0][i-startID] = f_vert[i];   record[1][i-startID] = i-startID+1;        
    }

    if(startFromCentor){
        fprintf(fp, "%g %g %g", centor[0], centor[1], centor[2]);
        for(int j=0; j<vec.size(); j++)fprintf(fp," %d", vec[j]); 
        fprintf(fp, "\n");

        for(int i=startID; i<startID+n; i++){
            int j;
            for(j=0; j<n; j++)if(record[0][j] == f_vert[i])break;
            fprintf(fp, "%d %d",n+1, record[1][j]);
            for(j=0; j<n; j++)if(record[0][j] == f_vert[(i<startID+n-1?i+1:startID)])break;
            fprintf(fp, " %d\n", record[1][j]);
        }

        
    }else{
        for(int i=startID+1; i<startID+n-1; i++){
            int j;
            for(j=1; j<n-1; j++)if(record[0][j] == f_vert[i])break;
            fprintf(fp, "1 %d", record[1][j]);
            for(j=2; j<n; j++)if(record[0][j] == f_vert[i+1])break;
            fprintf(fp, " %d\n", record[1][j]);
        }
    }
    return;
}
void DrawTetrahedron(FILE *fp, real8 rotMatrix[3][3], real8 r, real8 pos[3], std::vector<int> &vec, char *baseName)
{
    const real8 a=0.57735026918962576451;

//  real8 d[3]={-a,-a,-a};
//  real8 a[3]={a,-a,a};
//  real8 b[3]={-a,a,a};
//  real8 c[3]={a,a,-a};

    real8 vCryst[12]={-a,-a,-a,  a,-a,a,  -a,a,a,  a,a,-a}, vLab[12]; 
    std::vector<real8> vertices(12);
    for(int i=0;i<4;i++){Matrix33Vector3Multiply(rotMatrix,vCryst+3*i,vLab+3*i); normalize_vec(vLab+i*3);}
    for(int i=0;i<12;i++)vertices[i]=r*vLab[i];
    for(int i=0;i<4;i++){vertices[3*i]+=pos[0]; vertices[3*i+1]+=pos[1];vertices[3*i+2]+=pos[2];}

    std::vector<int> f_vert(16);
    f_vert[0]=3; f_vert[1]=2; f_vert[2]=3; f_vert[3]=0;     // a plane  ID: 0
    f_vert[4]=3; f_vert[5]=1; f_vert[6]=3; f_vert[7]=0;     // b plane  ID: 1
    f_vert[8]=3; f_vert[9]=1; f_vert[10]=2; f_vert[11]=0;   // g plane  ID: 2
    f_vert[12]=3; f_vert[13]=1; f_vert[14]=2; f_vert[15]=3; // d plane  ID: 3

    char zoneName[256];
    snprintf(zoneName,sizeof(zoneName),"%s-bcd",baseName);
    vec[1]=0; DrawPolygon(fp, f_vert, vertices, 0, zoneName, vec);

    snprintf(zoneName,sizeof(zoneName),"%s-acd",baseName);
    vec[1]=1; DrawPolygon(fp, f_vert, vertices, 4, zoneName, vec);

    snprintf(zoneName,sizeof(zoneName),"%s-abd",baseName);
    vec[1]=2; DrawPolygon(fp, f_vert, vertices, 8, zoneName, vec);

    snprintf(zoneName,sizeof(zoneName),"%s-abc",baseName);
    vec[1]=3; DrawPolygon(fp, f_vert, vertices, 12, zoneName, vec, true);
    return;
}
void DrawGrainBoundary(FILE *fp, Param_t *param, real8 z, char *zoneName, std::vector<int> &vec)
{
    real8 x=param->xBoundMax, y=param->yBoundMax;
    real8 v[12]={x,y,z, x,-y,z, -x,-y,z, -x,y,z};
    std::vector<real8> vertices(12);
    for(int i=0; i<12; i++)vertices[i]=v[i];
    std::vector<int> f_vert(5); f_vert[0]=4; for(int i=1;i<5;i++)f_vert[i]=i-1;    
    DrawPolygon(fp, f_vert, vertices, 0, zoneName, vec);
    return;
}

bool cmp3(SlipSys_t p, SlipSys_t q){
    return(p.schmid > q.schmid);
}
void InitGrain(Home_t *home, Param_t *param){
	int		    i,j,k;
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

    real8       zaxis[3] = {0.0,0.0,1.0}, maxSchmid, dis, rad = 800, sqi3 = 0.5773502691896258;
    std::vector<SlipSys_t>  slipSys;
    SlipSys_t               slip;
    BurgInfo_t  *bd;
    int         slipIndex, inGrain;
    GB_t        *gb;

	real8   hMatrix, p[4][3];
	real8   bound, rho;
	int	    inner = 0, interType;
	Plane_t plane;
	Line_t  edge;

	/*the following variables are placed in param.h*/
	int		    numTGrains;
	real8		halfTwinThick;
	int		    invFlag;
	real8		lx=0, ly=0, lz=0;	/*half cube size in param.h*/

    int                     m, n, nmesh;
    real8                   vec1[3], vec2[3], vec3[3], interPos[3];
    real8                   poss[2][3], vecMag, length, vec4[3], msize = 10, tsize;
    std::vector<Plane_t>    surList(6);
    std::vector<std::vector<Plane_t> > surArray(2);
    std::vector<int>        records, records1;
    FILE     *fp, *fp1;
    char     fileName[256];

    int         th, color, l, burgID, maxBurg, numMaxSchmids;

    real8       cryst[3], lab[3];
    real8       ps[4][3];
    real8       cena[3] = {1,-1, 1}, cenb[3] = {-1,1,1}, cenc[3] = {1,1,-1}, cend[3] = {-1,-1,-1};

	lx = param->xBoundMax;
	ly = param->yBoundMax;
	lz = param->zBoundMax;
	maxSeg = param->maxSeg;

	for(i=0; i<3; i++){
	    VECTOR_COPY(home->rotMatrix[i], param->matrixOrientation[i]);
	}
	Matrix33Invert(home->rotMatrix, home->rotMatrixInverse);

    GenerateFCCBurgInfo(&(home->burgData), param->matrixOrientation, -1);
    GenerateFCCBurgInfo(&(home->crystBurgData), I, -1);

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

	if(param->polycrystal == 0)param->nGrain = 1;
	nGrain = param->nGrain;

    if(type == 3)Fatal("Not supposed this type now");
/*
 *	Type I: polycrystal = 0, nGrain = 0, PBC
 * 	origin PBC condition.
 */	

/*
 *	Type II: polycrystal = 0, nGrain = 0, Free surface.
 *	single crystal, free surface.
 */
	if(param->polycrystal ==0){
        home->nGrains = 1;
        home->nGBs = 0;

	    home->grainKeys	= (Grain_t **)malloc(home->nGrains * sizeof(Grain_t *));
        home->grainKeys[0] = (Grain_t *)calloc(1, sizeof(Grain_t));
        grain = home->grainKeys[0];
        grain->nBounds = 0;
        GenerateFCCBurgInfo(&(grain->burgData), param->matrixOrientation, -1);

        memcpy(grain->orientation, home->rotMatrix, sizeof(real8 [3][3]));
        memcpy(grain->orientationInv, home->rotMatrixInverse, sizeof(real8 [3][3]));

        bd = &(home->burgData);
        for(j=0; j<6; j++){
            VECTOR_COPY(slip.burg, bd->burgList[j]);
            slip.burgID = j;
            for(k=0; k<2; k++){
                slip.planeID = bd->burgOnPlaneID[j][k];
                VECTOR_COPY(slip.normal, bd->planeList[bd->burgOnPlaneID[j][k]]);
                slip.schmid = fabs(DotProduct(slip.normal, zaxis)*DotProduct(slip.burg,zaxis));
                grain->slipSys.push_back(slip);
            }
        }

//        VECTOR_COPY(zaxis, param->loadDir);

	    if(home->myDomain == 0){
            std::vector<SlipSys_t>().swap(slipSys);
            bd = &(home->burgData);
            for(j=0; j<6; j++){
                VECTOR_COPY(slip.burg, bd->burgList[j]);
                slip.burgID = j;
                for(k=0; k<2; k++){
                    slip.planeID = bd->burgOnPlaneID[j][k];
                    VECTOR_COPY(slip.normal, bd->planeList[bd->burgOnPlaneID[j][k]]);
                    slip.schmid = fabs(DotProduct(slip.normal, zaxis)*DotProduct(slip.burg,zaxis));
                    slipSys.push_back(slip);
                }
            }
            sort(slipSys.begin(), slipSys.end(),cmp3); 
    
            maxSchmid = slipSys[0].schmid;

            for(j=1; j<slipSys.size(); j++){
                if(fabs(slipSys[j].schmid-maxSchmid)>1.0E-3)break;
            }
            slipSys.resize(j);

            param->maxSchmid = maxSchmid;
	        FormatMatrix(param->matrixOrientation,"Orientation");

            for(i=0; i<6; i++){
                VECTOR_ZERO(surList[i].normal);
                VECTOR_ZERO(surList[i].point);
            }
            surList[0].normal[0] =-1.0;
            surList[1].normal[0] = 1.0;
            surList[2].normal[1] =-1.0;
            surList[3].normal[1] = 1.0;
            surList[4].normal[2] =-1.0;
            surList[5].normal[2] = 1.0;

            surList[0].point[0] = param->xBoundMin;
            surList[1].point[0] = param->xBoundMax;
            surList[2].point[1] = param->yBoundMin;
            surList[3].point[1] = param->yBoundMax;
            surList[4].point[2] = param->zBoundMin;
            surList[5].point[2] = param->zBoundMax;

		    snprintf(fileName, sizeof(fileName), "primarry_planes.plt");

		    if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
			    Fatal("InitGrain: Open error %d on %s\n", errno, fileName);
		    }
            
	        fprintf(fp, "Variables = x, y, z, v1, v2, v3, planeID, cosf\n");

		    snprintf(fileName, sizeof(fileName), "total.plt");
		    if ((fp1 = fopen(fileName, "w")) == (FILE *)NULL) {
			    Fatal("InitGrain: Open error %d on %s\n", errno, fileName);
		    }
            
	        fprintf(fp1, "Variables = x, y, z, v1, v2, v3, index, planID, burg, schmid\n");

            for(j=0; j<slipSys.size(); j++){
                VECTOR_COPY(plane.normal, slipSys[j].normal);
                VECTOR_ZERO(plane.point);
                inGrain = slipSys[j].inGrain;
        
                for(k=0; k<6; k++){
                    type = PlanePlaneIntersection(
                              plane.normal, plane.point,
                              surList[k].normal, surList[k].point, vec1, vec2); 
                    
                    if(type != LINE_INTERSECTION)continue;

                    th = 0;

                    for(m=0; m<6; m++){
                        if(m==k)continue;
                        type = LinePlaneIntersection(vec1, vec2, surList[m].normal,
                                     surList[m].point, vec3, interPos);
                        if(type != POINT_INTERSECTION)continue;
                        if(fabs(interPos[0]) > param->xBoundMax+1.0){
                            continue;
                        }
                        if(fabs(interPos[1]) > param->yBoundMax+1.0){
                            continue;
                        }
                        if(fabs(interPos[2]) > param->zBoundMax+1.0){
                            continue;
                        }
                        VECTOR_COPY(poss[th], interPos);
                        th++;

                        if(th == 2){

                            vec4[0] = poss[1][0] - poss[0][0];
                            vec4[1] = poss[1][1] - poss[0][1];
                            vec4[2] = poss[1][2] - poss[0][2];

                            nmesh = ((int)(Normal(vec4)/msize));
                            tsize = Normal(vec4)/((real8)nmesh);
                            NormalizeVec(vec4);
            
                            for(l=0; l<nmesh; l++){
	                            fprintf(fp, "%e %e %e %e %e %e %d %f\n", poss[0][0],
                                poss[0][1],poss[0][2], tsize*vec4[0],
                                tsize*vec4[1],tsize*vec4[2], slipSys[j].planeID,
                                    fabs(DotProduct(slipSys[j].normal, zaxis)));

	                            fprintf(fp1, "%e %e %e %e %e %e 1 %d -1 %f\n", poss[0][0],
                                            poss[0][1],poss[0][2], tsize*vec4[0],
                                            tsize*vec4[1],tsize*vec4[2],  slipSys[j].planeID,
                                            fabs(DotProduct(slipSys[j].normal, zaxis)));

                                poss[0][0] += tsize*vec4[0];
                                poss[0][1] += tsize*vec4[1];
                                poss[0][2] += tsize*vec4[2];

                            }

                            th = 0;
                            break;
                        }

                    }

                }

            }

            fclose(fp);

		    snprintf(fileName, sizeof(fileName), "Thompson.plt");

		    if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
			    Fatal("InitGrain: Open error %d on %s\n", errno, fileName);
		    }
            
	        fprintf(fp, "Variables = x, y, z, v1, v2, v3, burgID, schmid\n");

            NormalizeVec(cena);
            NormalizeVec(cenb);
            NormalizeVec(cenc);
            NormalizeVec(cend);

            VECTOR_ZERO(plane.normal);
            VECTOR_ZERO(plane.point);

            Matrix33Vector3Multiply(param->matrixOrientation, cena, ps[0]);
            Matrix33Vector3Multiply(param->matrixOrientation, cenb, ps[1]);
            Matrix33Vector3Multiply(param->matrixOrientation, cenc, ps[2]);
            Matrix33Vector3Multiply(param->matrixOrientation, cend, ps[3]);

            for(j=0; j<4; j++){
                ps[j][0] *= rad;
                ps[j][1] *= rad;
                ps[j][2] *= rad;

                ps[j][0] += (plane.point[0]);
                ps[j][1] += (plane.point[1]);
                ps[j][2] += (plane.point[2]);
            }


            for(j=0; j<3; j++){
                for(k=j+1; k<4; k++){
                    vec4[0] = ps[k][0] - ps[j][0];
                    vec4[1] = ps[k][1] - ps[j][1];
                    vec4[2] = ps[k][2] - ps[j][2];

                    VECTOR_COPY(poss[0], ps[j]);

                    nmesh = ((int)(Normal(vec4)/msize));
                    tsize = Normal(vec4)/((real8)nmesh);

                    NormalizeVec(vec4);

                    burgID = BurgID_FCC(bd, vec4);
                    for(l=0; l<slipSys.size(); l++){
                        if(burgID == slipSys[l].burgID)break;
                    }
                    if(l==slipSys.size()){
                        maxBurg = 0;
                    }else{
                        maxBurg = 1;
                    }

                    for(l=0; l<nmesh; l++){
	                    fprintf(fp, "%e %e %e %e %e %e %d %d\n", poss[0][0],
                        poss[0][1],poss[0][2], tsize*vec4[0],
                        tsize*vec4[1],tsize*vec4[2],
                        burgID, maxBurg);

	                    fprintf(fp1, "%e %e %e %e %e %e %d -1 %d %f\n", poss[0][0],
                        poss[0][1],poss[0][2], tsize*vec4[0],
                        tsize*vec4[1],tsize*vec4[2],
                        maxBurg, burgID, fabs(DotProduct(bd->burgList[burgID], zaxis)));
                        poss[0][0] += tsize*vec4[0];
                        poss[0][1] += tsize*vec4[1];
                        poss[0][2] += tsize*vec4[2];

                    }

                    
                }
            }

            fclose(fp);
            fclose(fp1);

	    }

#if PARALLEL
       MPI_Bcast(&(param->maxSchmid), sizeof(double), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
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
	twinRatio = param->twinRatio;
	
/*
 *	Initialize twin grain.
 */
	twinRatio = param->twinRatio;
	numTGrains = nGrain - 1;

#if 0
	Matrix33Vector3Multiply(param->matrixOrientation, param->TBNormalCryst, param->TBNormal);
    RotateWithAxis(param->TBNormal, PI, temR);
	Matrix33Mult33(temR, param->matrixOrientation, TwinOrientation);
    RotateWithAxis(param->TBNormal, PI, param->rotateNormal);

#else
        real8   CDeltaLab[3] = {-1,-1,2}, CDelta[3];
        real8   planeABC[3] = {1,1,1};
        NormalizeVec(CDeltaLab);
        NormalizeVec(planeABC); 
        Matrix33Vector3Multiply(param->matrixOrientation, CDeltaLab, CDelta);

        RotateWithAxis(CDelta, PI, temR);
        Matrix33Mult33(temR, param->matrixOrientation, TwinOrientation);

        Matrix33Vector3Multiply(param->matrixOrientation, planeABC, param->TBNormal);
        RotateWithAxis(param->TBNormal, PI, param->rotateNormal);

	    if(home->myDomain == 0){
	        FormatMatrix(param->matrixOrientation,"basis orientation");
	        printf("\n");
	        FormatMatrix(TwinOrientation,"twin orientation");
	        printf("\n");
	        FormatVector(param->TBNormal, "twin boundary normal");
	        printf("\n");

            real8 loadCrystDir[3], matrixOrientationInv[3][3];
	        invFlag = Matrix33Invert(param->matrixOrientation, matrixOrientationInv);
	        if(!invFlag){
	            Fatal("Invert matrix orientation faild!\n");
	        }
            Matrix33Vector3Multiply(matrixOrientationInv, param->edotdir, loadCrystDir);
            FormatVector(loadCrystDir, "loading crystal direction");
	    }

#endif
/*
 *	Type IV: polycrystal = 1, nGrain >= 2, PBC.
 *	polycrystal, twin, RVE, XOZ plane is twin plane
 */
    if(type == 4){
#if 0
	    halfTwinThick = twinRatio * ly / (real8)numTGrains;
	    HalfMatrixThick = ly / (real8)numTGrains - halfTwinThick;
	    if (HalfMatrixThick <= 0){
	        printf("the thickness of matrix is %f\n ", HalfMatrixThick);
	        Fatal("the thickness of matrix is negative value at %s line %d\n", __FILE__, __LINE__);
	    }
#else
        real8 thick = 2*lz/((real8)nGrain);       
#endif	    	

	    snprintf(fileName, sizeof(fileName), "planes.plt");
		if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
		    Fatal("InitGrain: Open error %d on %s\n", errno, fileName);
		}
        fprintf(fp,"variables = x, y, z, id1, id2\n");
        home->nGBs = home->nGrains;
        home->GBKeys = (GB_t **)malloc(home->nGrains*sizeof(GB_t *));   
        for(i=0; i<home->nGrains; i++){
            home->GBKeys[i] = (GB_t *)calloc(1, sizeof(GB_t));
            home->GBKeys[i]->GBID = i;
            VECTOR_ZERO(home->GBKeys[i]->n);
            VECTOR_ZERO(home->GBKeys[i]->p);
            home->GBKeys[i]->n[2] = 1.0;
#if 0
            home->GBKeys[i]->p[1] = -ly + HalfMatrixThick + 2*((real8)i)*halfTwinThick;
#else
            home->GBKeys[i]->p[2] = -lz + (((real8)i)+0.5)*thick;
#endif
            home->GBKeys[i]->nbrGrain[0] = i;
            home->GBKeys[i]->nbrGrain[1] = ((i == home->nGrains -1) ? 0 : i+1);
        }

	    for(i=0; i<nGrain; i++){
	        home->grainKeys[i] = (Grain_t *)calloc(1, sizeof(Grain_t));
	        grain = home->grainKeys[i];
	        InitGrainStructure(grain, 2);
    
	        grain->type = i==0?FOLDED_GRAIN:ODINARY_GRAIN;
#if 0
	        grain->cy = -ly + (2 * ((real8)i) - 1) * (HalfMatrixThick + halfTwinThick);
#else
	        grain->cz = -lz + double(i)*thick;
#endif
            grain->cx = 0.0;
            grain->cy = 0.0;
            grain->grainID = i;

            GBID = i -1;
            if(GBID < 0) GBID = nGrain - 1;
            grain->GB[0] = home->GBKeys[GBID];
            grain->sign[0] = -1.0;
            grain->GB[1] = home->GBKeys[i];
            grain->sign[1] =  1.0;

            if(i%2==1){
                memcpy(grain->orientation, TwinOrientation, sizeof(real8 [3][3]));
            }else{
                memcpy(grain->orientation, param->matrixOrientation, sizeof(real8 [3][3]));
            }

            GenerateFCCBurgInfo(&(grain->burgData), grain->orientation, grain->grainID);
	        invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	        if(!invFlag){
	            Fatal("Invert twin orientation faild!\n");
	        }
        }

        char zoneName[256]; std::vector<int> indexs(2); real8 r=thick/3.0; real8 tPos[3];
        indexs[0]=i; indexs[1]=4;
        for(i=0; i<home->nGBs; i++){
            GenerateGBBurgInfo(home, home->GBKeys[i]);    
  
            gb = home->GBKeys[i];

            snprintf(zoneName,sizeof(zoneName),"GB-%d",i);
            DrawGrainBoundary(fp, param, gb->p[2], zoneName, indexs);
        }

        grain = home->grainKeys[0];
        tPos[0]=grain->cx;  tPos[1]=grain->cy;
        tPos[2]=-lz;
        snprintf(zoneName,sizeof(zoneName),"tetrahedron-0-1");
        DrawTetrahedron(fp, grain->orientation, r, tPos, indexs, zoneName);
        
        tPos[2]=lz;
        snprintf(zoneName,sizeof(zoneName),"tetrahedron-0-2");
        DrawTetrahedron(fp, grain->orientation, r, tPos, indexs, zoneName);

        for(i=1; i<home->nGrains; i++){
            grain = home->grainKeys[i];
            tPos[0]=grain->cx;  tPos[1]=grain->cy;  tPos[2]=grain->cz-0.15*thick;
            snprintf(zoneName,sizeof(zoneName),"tetrahedron-%d",i);
            DrawTetrahedron(fp, grain->orientation, r, tPos, indexs, zoneName);
        }

        fclose(fp);
        return;
    }

	if(nGrain == 2){

/*
 *	    Type VI: polycrystal = 1, nGrain =2
 *	    simulat twin sample.
 */

	    VECTOR_COPY(plane.normal, param->TBNormal);
	    VECTOR_ZERO(plane.point);

        if(fabs(0.5-twinRatio) > 1.0E-3){

            if(param->xBoundType != Periodic){
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
	 	    	  		type = LinePlaneIntersection(edge.direction, edge.point, 
	 	    	  				plane.normal, plane.point, vec1, interPos);
	 	    	  		if(type != POINT_INTERSECTION){
	 	    	  			inner = 0;
	 	    	  			break;
	 	    	  		} 
          
	                    if(fabs(interPos[i]) > bound){
	                        inner = 0; break; }
	                }/* for (j=0; j ...)*/
	                if(inner == 1)break;
          
	            }/* for(i=0; i ...)*/
	            
	            if(inner == 0)Fatal("InitGrain: can not find TB plane");
            }else{
                plane.point[1] = param->yBoundMin + (1-twinRatio)*(param->yBoundMax - param->yBoundMin);
            }
        }


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
        GenerateFCCBurgInfo(&(grain->burgData), grain->orientation, grain->grainID);
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
        GenerateFCCBurgInfo(&(grain->burgData), grain->orientation, grain->grainID);
	    invFlag = Matrix33Invert(grain->orientation, grain->orientationInv);
	    if(!invFlag){
	        Fatal("Invert twin orientation faild!\n");
	    }
        for(i=0; i<home->nGBs; i++){
            GenerateGBBurgInfo(home, home->GBKeys[i]);    
        }

	    if(home->myDomain == 0){
	        FormatMatrix(param->matrixOrientation,"basis orientation");
	        printf("\n");
	        FormatMatrix(TwinOrientation,"twin orientation");
	        printf("\n");
	        FormatVector(param->TBNormal, "twin boundary normal");
	        printf("\n");
	    }

        if(home->myDomain == 0){

            for(i=0; i<6; i++){
                VECTOR_ZERO(surList[i].normal);
                VECTOR_ZERO(surList[i].point);
            }
            surList[0].normal[0] =-1.0;
            surList[1].normal[0] = 1.0;
            surList[2].normal[1] =-1.0;
            surList[3].normal[1] = 1.0;
            surList[4].normal[2] =-1.0;
            surList[5].normal[2] = 1.0;

            surList[0].point[0] = param->xBoundMin;
            surList[1].point[0] = param->xBoundMax;
            surList[2].point[1] = param->yBoundMin;
            surList[3].point[1] = param->yBoundMax;
            surList[4].point[2] = param->zBoundMin;
            surList[5].point[2] = param->zBoundMax;

		    snprintf(fileName, sizeof(fileName), "TB.plt");

		    if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
			    Fatal("InitGrain: Open error %d on %s\n", errno, fileName);
		    }
            
	        fprintf(fp, "Variables = x, y, z, v1, v2, v3, burgID, normalID, schmid, inGrain\n");

            for(i=0; i<6; i++){
                type = PlanePlaneIntersection(surList[i].normal, surList[i].point,
                                              plane.normal, plane.point, vec1, vec2); 
                if(type != LINE_INTERSECTION){
                    records1.push_back(i);
                    continue;
                }
                records.push_back(i);

                k = 0;
                for(j=0; j<6; j++){
                    if(j==i)continue;
                    type = LinePlaneIntersection(vec1, vec2, surList[j].normal, 
                                                 surList[j].point, vec3, interPos);
                    if(type != POINT_INTERSECTION){
                        continue;
                    }
                    if(fabs(interPos[0]) > param->xBoundMax+1.0){
                        continue;
                    }
                    if(fabs(interPos[1]) > param->yBoundMax+1.0){
                        continue;
                    }
                    if(fabs(interPos[2]) > param->zBoundMax+1.0){
                        continue;
                    }
                    VECTOR_COPY(poss[k], interPos);
                    k++;
                    if(k==2){
                        vec4[0] = poss[1][0] - poss[0][0];
                        vec4[1] = poss[1][1] - poss[0][1];
                        vec4[2] = poss[1][2] - poss[0][2];
                        length = Normal(vec4);
                        nmesh = int(length/msize);
                        tsize = length/((real8)nmesh);
                        NormalizeVec(vec4);
                        for(m=0; m<nmesh; m++){
#if 1
	                        fprintf(fp, "%e %e %e %e %e %e -1 3 -1.00 -1\n", poss[0][0],
                                    poss[0][1],poss[0][2], tsize*vec4[0],
                                    tsize*vec4[1], tsize*vec4[2]);
#endif
                            poss[0][0] += tsize*vec4[0];
                            poss[0][1] += tsize*vec4[1];
                            poss[0][2] += tsize*vec4[2];
                        }
                        k = 0;
                        break;
                    }
                }
            }

            surArray[0].push_back(plane);
            home->grainKeys[0]->surface.push_back(plane);
            plane.normal[0] = -plane.normal[0];
            plane.normal[1] = -plane.normal[1];
            plane.normal[2] = -plane.normal[2];

            surArray[1].push_back(plane);
            home->grainKeys[1]->surface.push_back(plane);
            plane.normal[0] = -plane.normal[0];
            plane.normal[1] = -plane.normal[1];
            plane.normal[2] = -plane.normal[2];

            for(i=0; i<records.size(); i++){
                surArray[0].push_back(surList[records[i]]);
                surArray[1].push_back(surList[records[i]]);
                home->grainKeys[0]->surface.push_back(surList[records[i]]);
                home->grainKeys[1]->surface.push_back(surList[records[i]]);
            }

            for(i=0; i<records1.size(); i++){
                if(DotProduct(surList[records1[i]].normal,plane.normal)<0){
                    surArray[0].push_back(surList[records1[i]]);
                    home->grainKeys[0]->surface.push_back(surList[records1[i]]);
                }else{
                    surArray[1].push_back(surList[records1[i]]);
                    home->grainKeys[1]->surface.push_back(surList[records1[i]]);
                }
            }
#if 0
            for(i=0; i<surArray[0].size(); i++){
                printf("0 %d: {%f,%f,%f}--{%f,%f,%f}\n", i,
                    surArray[0][i].normal[0],surArray[0][i].normal[1],surArray[0][i].normal[2],
                    surArray[0][i].point[0],surArray[0][i].point[1],surArray[0][i].point[2]);
            }
            for(i=0; i<surArray[1].size(); i++){
                printf("1 %d: {%f,%f,%f}--{%f,%f,%f}\n",i, 
                    surArray[1][i].normal[0],surArray[1][i].normal[1],surArray[1][i].normal[2],
                    surArray[1][i].point[0],surArray[1][i].point[1],surArray[1][i].point[2]);
            }
#endif

            NormalizeVec(cena);
            NormalizeVec(cenb);
            NormalizeVec(cenc);
            NormalizeVec(cend);

            gb = home->GBKeys[0];
            for(i=0; i<home->nGrains; i++){
                grain = home->grainKeys[i];
#if 1
                VECTOR_COPY(plane.normal, gb->n);
                VECTOR_COPY(plane.point, gb->p);

                plane.point[0] -= (grain->sign[0]*plane.normal[0]*500);
                plane.point[1] -= (grain->sign[0]*plane.normal[1]*500);
                plane.point[2] -= (grain->sign[0]*plane.normal[2]*500);

                Matrix33Vector3Multiply(grain->orientation, cena, ps[0]);
                Matrix33Vector3Multiply(grain->orientation, cenb, ps[1]);
                Matrix33Vector3Multiply(grain->orientation, cenc, ps[2]);
                Matrix33Vector3Multiply(grain->orientation, cend, ps[3]);

                for(j=0; j<4; j++){
                    ps[j][0] *= rad;
                    ps[j][1] *= rad;
                    ps[j][2] *= rad;

                    ps[j][0] += (plane.point[0]);
                    ps[j][1] += (plane.point[1]);
                    ps[j][2] += (plane.point[2]);
                }

                bd = &(grain->burgData);

                for(j=0; j<3; j++){
                    for(k=j+1; k<4; k++){
                        vec4[0] = ps[k][0] - ps[j][0];
                        vec4[1] = ps[k][1] - ps[j][1];
                        vec4[2] = ps[k][2] - ps[j][2];

                        VECTOR_COPY(poss[0], ps[j]);

                        nmesh = ((int)(Normal(vec4)/msize));
                        tsize = Normal(vec4)/((real8)nmesh);

                        NormalizeVec(vec4);

                        for(l=0; l<nmesh; l++){
	                        fprintf(fp, "%e %e %e %e %e %e %d -1 -1.0 %d\n", poss[0][0],
                            poss[0][1],poss[0][2], tsize*vec4[0],
                            tsize*vec4[1],tsize*vec4[2], BurgID_FCC(bd, vec4), i);

                            poss[0][0] += tsize*vec4[0];
                            poss[0][1] += tsize*vec4[1];
                            poss[0][2] += tsize*vec4[2];

                        }
                    
                    }
                }
                
#endif
                bd = &(grain->burgData);
                slip.inGrain = i;
                for(j=0; j<6; j++){
                    VECTOR_COPY(slip.burg, bd->burgList[j]);
                    slip.burgID = j;
                    for(k=0; k<2; k++){
                        slip.planeID = bd->burgOnPlaneID[j][k];
                        VECTOR_COPY(slip.normal, bd->planeList[bd->burgOnPlaneID[j][k]]);
                        slip.schmid = fabs(DotProduct(slip.normal, zaxis)*DotProduct(slip.burg,zaxis));
                        slipSys.push_back(slip);
                    }
                }
                sort(slipSys.begin(), slipSys.end(),cmp3); 
    
                maxSchmid = slipSys[0].schmid;
                printf("max %f\n",  maxSchmid);

                grain->maxSchmid = maxSchmid;

                for(j=1; j<slipSys.size(); j++){
                    if(fabs(slipSys[j].schmid-maxSchmid)>1.0E-3)break;
                }

                numMaxSchmids = j;

                for(j=0; j<numMaxSchmids; j++){
                    VECTOR_COPY(plane.normal, slipSys[j].normal);
                    VECTOR_ZERO(plane.point);
                    inGrain = slipSys[j].inGrain;
                    plane.point[0] += -10.0*surArray[inGrain][0].normal[0];
                    plane.point[1] += -10.0*surArray[inGrain][0].normal[1];
                    plane.point[2] += -10.0*surArray[inGrain][0].normal[2];


                    for(k=0; k<surArray[inGrain].size(); k++){
                        type = PlanePlaneIntersection(
                                  plane.normal, plane.point,
                                  surArray[inGrain][k].normal, surArray[inGrain][k].point, vec1, vec2); 
                        
                        if(type != LINE_INTERSECTION)continue;
                        th = 0;
                        for(m=0; m<surArray[inGrain].size(); m++){
                            if(m==k)continue;
                            type = LinePlaneIntersection(vec1, vec2, surArray[inGrain][m].normal,
                                         surArray[inGrain][m].point, vec3, interPos);

                            if(type != POINT_INTERSECTION)continue;

                            type = PointPlaneIntersection(interPos, surArray[inGrain][0].normal,
                                       surArray[inGrain][0].point, vec4, &dis);
                            if(dis > 0.1)continue;
                            if(fabs(interPos[0]) > param->xBoundMax+1.0){
                                continue;
                            }
                            if(fabs(interPos[1]) > param->yBoundMax+1.0){
                                continue;
                            }
                            if(fabs(interPos[2]) > param->zBoundMax+1.0){
                                continue;
                            }

                            VECTOR_COPY(poss[th], interPos);
                            th++;
#if 0
                            if(th == 2){
                                vec4[0] = poss[1][0] - poss[0][0];
                                vec4[1] = poss[1][1] - poss[0][1];
                                vec4[2] = poss[1][2] - poss[0][2];
                                length = Normal(vec4);
                                nmesh = (int)(length/msize);
                                tsize = length/((real8)nmesh);
                                NormalizeVec(vec4);
                                color = slipSys[j].burgID;
            
                                for(n=0; n<nmesh; n++){
	                                fprintf(fp, "%e %e %e %e %e %e %d %d %f %d\n", poss[0][0],
                                    poss[0][1],poss[0][2], tsize*vec4[0],
                                    tsize*vec4[1],tsize*vec4[2], slipSys[j].burgID, 
                                    slipSys[j].planeID, slipSys[j].schmid, inGrain);
                                    poss[0][0] += tsize*vec4[0];
                                    poss[0][1] += tsize*vec4[1];
                                    poss[0][2] += tsize*vec4[2];
                                }
                                th = -1;
                                break;
                            }
#endif
                        }
                    
                    }
                }

                swap(grain->slipSys, slipSys);

            }

            fclose(fp);
            
        }
	}else{
	    Fatal("Not suppose free surface polycrystal problem");
	}

	return;
}


void NodeLocation(Home_t *home, real8 pos[3], int& inGrain, int& onGB, real8 err){
    int         i, j, nGrains, type, breakLoop = 0;
    real8       vec[3], dis, sign;
	Grain_t		*grain;
	Param_t		*param;

    inGrain = -1;
    onGB = -1;
	param = home->param;

	if(param->polycrystal == 0){
        return;
    }

	nGrains = home->nGrains;

    FoldBox(param, pos, pos+1, pos+2);

    for(i=0; i<nGrains; i++){
        if(breakLoop)break;
        grain = home->grainKeys[i];
        for(j=0; j<grain->nBounds; j++){
            type = PointPlaneIntersection(pos, grain->GB[j]->n, 
                   grain->GB[j]->p, vec, &dis);
            dis *= grain->sign[j];

            if(fabs(dis)<err){
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
    if(inGrain == -1 && onGB ==  -1){
        Fatal("in %s at %d", __FILE__, __LINE__);
    }

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


BurgInfo_t* GetBurgData(Home_t *home, Node_t *node1, Node_t *node2)
{
        if(home->param->polycrystal){
            int     inGrain, onGB;
            Grain_t *grain;
            GB_t    *gb;

            GetArmGrain(home, node1, node2, inGrain, onGB);
            if(onGB > -1){
                gb = home->GBKeys[onGB];
                return(&(gb->burgData));
            }else{
                grain = home->grainKeys[inGrain];    
                return(&(grain->burgData));
            }
        }else{
            if(!home->param->useLabFrame){
                return(&(home->crystBurgData));               
            }else{
                return(&(home->burgData));               
            }
        }
        return((BurgInfo_t *)NULL);
}

BurgInfo_t* GetBurgData(Home_t *home, int inGrain, int onGB)
{
        if(home->param->polycrystal){
            Grain_t *grain;
            GB_t    *gb;

            if(inGrain == -1 && onGB ==  -1){
                Fatal("in %s at %d", __FILE__, __LINE__);
            }
            if(onGB > -1){
                gb = home->GBKeys[onGB];
                return(&(gb->burgData));
            }else{
                grain = home->grainKeys[inGrain];    
                return(&(grain->burgData));
            }
        }else{
            if(!home->param->useLabFrame){
                return(&(home->crystBurgData));               
            }else{
                return(&(home->burgData));               
            }
        }
        return((BurgInfo_t *)NULL);
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
#if 0
        if(i == burgData->numBurgVectors){
            printf("Warning: in %s at %d\n", __FILE__, __LINE__);
        }
#endif        
        return(i);
}

int PlaneID_FCC(BurgInfo_t *burgData, real8 plane[3])
{
        int     i;
        real8   crossProd[3];
        
        for(i=0; i<burgData->numPlanes; i++){
            cross(burgData->planeList[i], plane, crossProd);
            if(Normal(crossProd) < 1.0E-3)break;
        }
#if 0
        if(i == burgData->numPlanes){
            printf("Warning: in %s at %d\n", __FILE__, __LINE__);
        }
#endif        
        return(i);
}

int KeyLineID_FCC(BurgInfo_t *burgData, real8 line[3])
{
        int     i;
        real8   crossProd[3];
        
        for(i=0; i<burgData->numKeyLines; i++){
            cross(burgData->keyLineList[i], line, crossProd);
            if(Normal(crossProd) < 1.0E-3)break;
        }
        
        if(i == burgData->numKeyLines){
            printf("Warning: in %s at %d\n", __FILE__, __LINE__);
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
	if(fabs(burgMag - 1.0E0) < eps){
		return(BASIC_BURG);
	}
/*	
 *	2: Shockley burgers
 */

	
	if(fabs(burgMag - SHOCKLEY_BURG_MAG_FCC) < eps) {
		return(SHOCKLEY_BURG);
	}

/*	
 *	3: junction 1 burgers, like [Sqrt(2),0,0]
 */
	if(fabs(burgMag - JUNCTION1_BURG_MAG_FCC) < eps) {
		return(JUNCTION1_BURG);
	}

/*	
 *	3: junction 2 burgers, like DA+DB
 */
	if(fabs(burgMag - JUNCTION2_BURG_MAG_FCC) < eps) {
		return(JUNCTION2_BURG);
	}
	return(UNKNOWN_BURG);

}


int PlaneType_FCC(BurgInfo_t *burgData, real8 plane[3])
{
    int planeID;

    planeID = PlaneID_FCC(burgData, plane);

    if(planeID >= NUM_PLANES_IN_FCC) return(UNKOWN_PLANE);
 
    return(burgData->planeType[planeID]);
}

int GetBurgID_FCC(Home_t *home, Node_t *node, Node_t *nbr){
       
    int         armID;
    real8       burg[3];
    BurgInfo_t  *burgData;
    burgData = GetBurgData(home, node, nbr);

    armID = GetArmID(home, node, nbr);
    burg[0] = node->burgX[armID];
    burg[1] = node->burgY[armID];
    burg[2] = node->burgZ[armID];
        
    return(BurgID_FCC(burgData, burg)); 
}

int GetPlaneID_FCC(Home_t *home, Node_t *node, Node_t *nbr){
       
    int         armID;
    real8       plane[3];
    BurgInfo_t  *burgData;
    burgData = GetBurgData(home, node, nbr);

    armID = GetArmID(home, node, nbr);
    if(armID == -1 )Fatal("in %s at %d", __FILE__, __LINE__);
    plane[0] = node->nx[armID];
    plane[1] = node->ny[armID];
    plane[2] = node->nz[armID];
        
    return(PlaneID_FCC(burgData, plane)); 
}

int GetPlaneType_FCC(Home_t *home, Node_t *node, Node_t *nbr){
       
    int         armID;
    real8       plane[3];
    BurgInfo_t  *burgData;
    burgData = GetBurgData(home, node, nbr);

    armID = GetArmID(home, node, nbr);
    if(armID == -1 )Fatal("in %s at %d", __FILE__, __LINE__);
    plane[0] = node->nx[armID];
    plane[1] = node->ny[armID];
    plane[2] = node->nz[armID];
        
    return(PlaneType_FCC(burgData, plane)); 
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
	int 	i, numNbrs, type = 0;
	Node_t 	*nbr;
	
	if(node == (Node_t *)NULL){
		Fatal("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
	}

	if(!IsOnGB(node))return(INGRAIN_NODE);

	numNbrs = node->numNbrs;
	
	for(i=0; i<numNbrs; i++){
        if((nbr = GetNeighborNode(home, node, i)) 
            == (Node_t *)NULL){
            PrintNode(home, node);
#ifdef DEBUG
            Fatal("Task %d:  NodeOnGB, unkonwn neighbor, flag=%d\n",home->myDomain, home->flag);
#else
            Fatal("Task %d:  NodeOnGB, unkonwn neighbor\n",home->myDomain);
#endif
            continue;
        }
		
        if(IsOnGB(nbr)){
            type |= GB_NODE;
        }else{
            type |= INTERSECT_NODE;
        }
	}
	return(type);
}

int AdjustNodeDirection(Home_t *home, Node_t *node, real8 *vec){

    int 	NM;
    real8 	GC[3][3];
    bool    zeroLenJunc;

	NM = 0;

	GetGlideConstraintList(home, node, &NM, GC, zeroLenJunc);
	ApplyConstraintsToVelocity(NM, GC, vec);

	return(NM);
}


real8 ClosedGB(Home_t *home, real8 *pos, GB_t **GB){
	real8	    tempDistance, vec[3]={pos[0],pos[1],pos[2]},vec1[3];
	real8 	    minDistance = 1e14;
	GB_t 	    *temGB;
	    
    if(!home->param->polycrystal)return(0.0);
   
    int     i;
    FoldBox(home->param, vec, vec+1, vec+2);
    for(i=0; i<home->nGBs; i++){
        temGB = home->GBKeys[i];
        PointPlaneIntersection(vec, temGB->n, temGB->p, vec1, &tempDistance);
        if(fabs(tempDistance) < fabs(minDistance)){
            *GB = temGB;
            minDistance = tempDistance;
        }
    }
	return(minDistance);
}

int GetGBIndexOfGrain(Grain_t *grain, int onGB){
    int i;
    for(i=0; i<grain->nBounds; i++){
         if(grain->GB[i]->GBID == onGB)return(i);
    }
    return(-1);
}

bool GetPublicGB(Home_t *home, Grain_t *grain1, Grain_t *grain2, GB_t **GB, real8 *pos)
{
    *GB = (GB_t *)NULL;

    int id, i;
    if(home->nGrains!=2){
        for(i=0; i<grain1->nBounds; i++){
            id = GetGBIndexOfGrain(grain1, grain2->GB[i]->GBID);
            if(id > -1){
                *GB = grain1->GB[id];
                return(1);
            } 
        }
    }else{
        real8 dis=ClosedGB(home,pos,GB); return true;
    }
    return(0);
}

bool NewPosToGB(Home_t *home, Node_t *node, int inGrain, int onGB, real8 pos[3], real8 *newPos){
    
    real8       dis, vec[3], vec1[3], vec2[3], eps;
    real8       interPos[3], vec3[3], GC[3][3]; 
    int         type, NM, index, curInGrain, curOnGB;
    bool        rval;
    Param_t     *param;
    GB_t        *gb;
    Grain_t     *grain;
    Plane_t     plane;
    bool        zeroLenJunc;    

    NM = 0;
   
    param = home->param;

    if(param->polycrystal == 0){
        return(0);
    }

    if(node->constraint ==7){
        printf("Task %s %d: ",__func__,home->myDomain);
        PrintTag(home, node);
        printf("move from {%g,%g,%g} to grain %d\n",pos[0],pos[1],pos[2],inGrain);
    }
    if(onGB < 0 || onGB >= home->nGBs)Fatal("in %s at %d", __FILE__, __LINE__);
	if(node == (Node_t *)NULL){
        Fatal("Task %d: %s, line %d -- Node = NULL",
               home->myDomain, __FILE__, __LINE__);
	}

    eps = param->rTol;

	FoldBox(param, pos, pos+1, pos+2);

    NodeLocation(home, pos[0], pos[1], pos[2], curInGrain, curOnGB);
    if(inGrain == -1){
        rval = (curOnGB == onGB) ? 0 : 1;
    }else{
        rval = (inGrain == curInGrain) ? 0 : 1;
    }

    gb = home->GBKeys[onGB];
    VECTOR_COPY(plane.normal, gb->n);
    VECTOR_COPY(plane.point, gb->p);

    if(inGrain > -1){
        if(inGrain >= home->nGrains)Fatal("in %s at %d", __FILE__, __LINE__);
        grain = home->grainKeys[inGrain];
        index = GetGBIndexOfGrain(grain, onGB);
        if(index < 0)Fatal("in %s at %d", __FILE__, __LINE__);
        plane.point[0] -= (grain->sign[index] * eps * gb->n[0]);
        plane.point[1] -= (grain->sign[index] * eps * gb->n[1]);
        plane.point[2] -= (grain->sign[index] * eps * gb->n[2]);
    }
    
    PointPlaneIntersection(pos, plane.normal, plane.point, newPos, &dis);
    if(fabs(dis) < 1E-4){
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(rval);
    }

    vec[0] = newPos[0] - pos[0];
    vec[1] = newPos[1] - pos[1];
    vec[2] = newPos[2] - pos[2];

    if(Normal(vec) < 1E-10){
        PrintNode(home, node);
        printf("Warning in NewPosToGB: {%f,%f,%f} too short of the vector\n",
                vec[0], vec[1], vec[2]);
        VECTOR_COPY(newPos, pos);
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(rval);
    }

    NormalizeVec(vec);
    VECTOR_COPY(vec1, vec);
	GetGlideConstraintList(home, node, &NM, GC, zeroLenJunc, 0, 
                           plane.normal[0], plane.normal[1],
                           plane.normal[2]);
	ApplyConstraintsToVelocity(NM, GC, vec);

    if(Normal(vec) < 1.0E-6 && node->constraint != 7){
        while(NM >= 0){
            NM--;
            VECTOR_COPY(vec, vec1);
            ApplyConstraintsToVelocity(NM, GC, vec);
            if(Normal(vec) >= 1.0E-6)break;
        }
    }

    if(Normal(vec) < 1.0E-6){
        VECTOR_COPY(newPos, pos);
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(rval);
    }

    NormalizeVec(vec);
    type = LinePlaneIntersection(vec, pos, plane.normal, plane.point, vec1, vec2);
    if(type == NO_INTERSECTION){
#if 0
        printf("Warning: NewPosToGB, Line {%f,%f,%f}--{%f,%f,%f}; ", 
               vec[0], vec[1], vec[2], pos[0], pos[1], pos[2]);
        printf("Plane {%f,%f,%f}--{%f,%f,%f} == {%f,%f,%f}\n", 
               plane.normal[0],plane.normal[1],plane.normal[2],
               plane.point[0],plane.point[1],plane.point[2],
               vec2[0], vec2[1], vec2[2]);
#endif
        VECTOR_COPY(newPos, vec2);
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(0);
    }else{
        VECTOR_COPY(newPos, vec2);
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(0);
    }

    return(1);
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
            GB_t    *gb; real8 node1Pos[3]={node1->x, node1->y, node1->z};
            int     node1SegID;

            if(Connected(node1, node2, &node1SegID)){
                PrintNode(node1);
                PrintNode(node2);
#ifdef DEBUG
                Fatal("%d,%d: GetArmGrain: node1 in %d grain, but node2 in %d grain\n",
                      home->myDomain,home->flag,node1->inGrain, node2->inGrain);
#else
                Fatal("%d: GetArmGrain: node1 in %d grain, but node2 in %d grain\n",
                      home->myDomain,node1->inGrain, node2->inGrain);
#endif
            }

            if(GetPublicGB(home, home->grainKeys[node1->inGrain], 
                        home->grainKeys[node2->inGrain], &gb, node1Pos)){
                inGrain = -1;
                onGB = gb->GBID;
            }else{
                PrintNode(node1);
                PrintNode(node2);
#ifdef DEBUG
                Fatal("%d,%d: GetArmGrain: node1 in %d grain, but node2 in %d grain\n",
                      home->myDomain,home->flag,node1->inGrain, node2->inGrain);
#else
                Fatal("%d: GetArmGrain: node1 in %d grain, but node2 in %d grain\n",
                      home->myDomain,node1->inGrain, node2->inGrain);
#endif
            }
        }
    }

    if(inGrain == -1 && onGB ==  -1){
        Fatal("in %s at %d", __FILE__, __LINE__);
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
        printf("[%d %d %d %d %d %d]-{%f,%f,%f} ",node->myTag.domainID,
               node->myTag.index,node->numNbrs,node->inGrain,
               node->onGB, node->constraint, node->x,node->y,node->z);
	}
    return;
}
void PrintTag(Home_t *home, Node_t *node){

    if(node == (Node_t *)NULL){
        printf("(NULL)"); return;
    }else{
        printf("[%d %d %d %d %d %d]-{%f,%f,%f} ",node->myTag.domainID,
               node->myTag.index,node->numNbrs,node->inGrain,
               node->onGB, node->constraint, node->x,node->y,node->z);
	}
    int numGlideConstraints; bool  zeroLenJunc;
    real8 glideConstraints[3][3];
    GetGlideConstraintList(home, node, 
                       &numGlideConstraints, glideConstraints, zeroLenJunc,
                       1, 0.0, 0.0, 0.0, -1);
    if(numGlideConstraints== 1){
        printf("plane constraint: {%g,%g,%g} ", glideConstraints[0][0],glideConstraints[0][1],glideConstraints[0][2]);
    }else if(numGlideConstraints== 2){
        printf("line constraint: {%g,%g,%g} ", glideConstraints[2][0],glideConstraints[2][1],glideConstraints[2][2]);
    }else{
        printf("completely constraint");
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
void ImageForce(Home_t *home, Node_t *node1, Node_t *node2, real8 f1[3], real8 f2[3])
{
    Param_t *param = home->param;

    real8       dis, burgMag, burg[3], segCentor[3], tbStrength, xi[3];
    real8       ldir[3],pdir[3], a = 0.5*param->rTol, amax=10*param->rTol, err=1.0E-4; 
    real8       gausFac[2] = {0.211324865405187, // (1-sqrt(3)/3)/2
                              0.7886751345948128}; // (1+sqrt(3)/3)/2
    real8       coff=0.0, length, nor[3], gDir[3], gbNormal[3], pos[3], pBound[3], nBound[3], vel[3];
    GB_t        *gb;
    Grain_t     *grain;
    int         armID, i, j, type, armInGrain,armOnGB;

    VECTOR_ZERO(f1);
    VECTOR_ZERO(f2);

    armID = GetArmID(home, node1, node2);
    if(armID == -1)return;

    burg[0] = node1->burgX[armID];
    burg[1] = node1->burgY[armID];
    burg[2] = node1->burgZ[armID];

    nor[0] = node1->nx[armID];
    nor[1] = node1->ny[armID];
    nor[2] = node1->nz[armID];

    burgMag = Normal(burg);
    NormalizeVec(burg);

    xi[0] = node2->x - node1->x;
    xi[1] = node2->y - node1->y;
    xi[2] = node2->z - node1->z;
    ZImage(param, xi, xi+1, xi+2);

    length = Normal(xi);
    if(length < 1.0E-2)return;
    NormalizeVec(xi);
    GetArmGrain(home, node1, node2, armInGrain, armOnGB);

    vel[0] = 0.5*(node2->vX+node1->vX);
    vel[1] = 0.5*(node2->vY+node1->vY);
    vel[2] = 0.5*(node2->vZ+node1->vZ);

    for(i=0; i<2; i++){
        pos[0] = node1->x + length*gausFac[i]*xi[0];
        pos[1] = node1->y + length*gausFac[i]*xi[1];
        pos[2] = node1->z + length*gausFac[i]*xi[2];
        FoldBox(param, pos, pos+1, pos+2);
/*
 *      GB image force
 */
        if(param->gbStrength > 0){
            if(param->xBoundType==Periodic){
                nBound[0]=1.0; nBound[1]=0; nBound[2]=0;
                pBound[0]=param->xBoundMin; pBound[1]=0; pBound[2]=0;
                type = PlanePlaneIntersection(nor,pos,nBound,pBound,ldir,pdir);
                if(type==LINE_INTERSECTION){
                    CoordMinus(pos,pdir,ldir); 
                    if((dis=Normal(ldir))>err){
                        NormalizeVec(ldir);
                        coff = param->gbStrength * a/(fabs(dis)+a);
                        NormalComponent(ldir,xi,gDir);
                        if(DotProduct(vel,gDir)>0){
                            f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                            f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                            f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                        }
                    }
                }
                pBound[0]=param->xBoundMax; pBound[1]=0; pBound[2]=0;
                type = PlanePlaneIntersection(nor,pos,nBound,pBound,ldir,pdir);
                if(type==LINE_INTERSECTION){
                    CoordMinus(pos,pdir,ldir); 
                    if((dis=Normal(ldir))>err){
                        NormalizeVec(ldir);
                        coff = param->gbStrength * a/(fabs(dis)+a);
                        NormalComponent(ldir,xi,gDir);
                        if(DotProduct(vel,gDir)>0){
                            f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                            f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                            f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                        }
                    }
                }
            }

            if(param->yBoundType==Periodic){
                nBound[0]=0.0; nBound[1]=1.0; nBound[2]=0;
                pBound[0]=0; pBound[1]=param->yBoundMin; pBound[2]=0;
                type = PlanePlaneIntersection(nor,pos,nBound,pBound,ldir,pdir);
                if(type==LINE_INTERSECTION){
                    CoordMinus(pos,pdir,ldir); 
                    if((dis=Normal(ldir))>err){
                        NormalizeVec(ldir);
                        coff = param->gbStrength * a/(fabs(dis)+a);
                        NormalComponent(ldir,xi,gDir);
                        if(DotProduct(vel,gDir)>0){
                            f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                            f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                            f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                        }
                    }
                }
                pBound[0]=0; pBound[1]=param->yBoundMax; pBound[2]=0;
                type = PlanePlaneIntersection(nor,pos,nBound,pBound,ldir,pdir);
                if(type==LINE_INTERSECTION){
                    CoordMinus(pos,pdir,ldir); 
                    if((dis=Normal(ldir))>err){
                        NormalizeVec(ldir);
                        coff = param->gbStrength * a/(fabs(dis)+a);
                        NormalComponent(ldir,xi,gDir);
                        if(DotProduct(vel,gDir)>0){
                            f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                            f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                            f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                        }
                    }
                }
            }

            if(param->zBoundType==Periodic){
                nBound[0]=0.0; nBound[1]=0.0; nBound[2]=1.0;
                pBound[0]=0; pBound[1]=0; pBound[2]=param->zBoundMin;
                type = PlanePlaneIntersection(nor,pos,nBound,pBound,ldir,pdir);
                if(type==LINE_INTERSECTION){
                    CoordMinus(pos,pdir,ldir); 
                    if((dis=Normal(ldir))>err){
                        NormalizeVec(ldir);
                        coff = param->gbStrength * a/(fabs(dis)+a);
                        NormalComponent(ldir,xi,gDir);
                        if(DotProduct(vel,gDir)>0){
                            f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                            f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                            f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                        }
                    }
                }
                pBound[0]=0; pBound[1]=0; pBound[2]=param->zBoundMax;
                type = PlanePlaneIntersection(nor,pos,nBound,pBound,ldir,pdir);
                if(type==LINE_INTERSECTION){
                    CoordMinus(pos,pdir,ldir); 
                    if((dis=Normal(ldir))>err){
                        NormalizeVec(ldir);
                        coff = param->gbStrength * a/(fabs(dis)+a);
                        NormalComponent(ldir,xi,gDir);
                        if(DotProduct(vel,gDir)>0){
                            f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                            f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                            f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                        }
                    }
                }
            }
        }
/*
 *      TB image force
 */
        if(armInGrain<0 || param->polycrystal==0)continue; 
        tbStrength = (PerpendicularDirection(burg,param->TBNormal,0.05)==true?param->tbStrength[0]:param->tbStrength[1]);
        
        grain=home->grainKeys[armInGrain];
        for(j=0; j<grain->nBounds; j++){
            gb = grain->GB[j];
            type = PlanePlaneIntersection(nor,pos,gb->n,gb->p,ldir,pdir);
            if(type==LINE_INTERSECTION){
                CoordMinus(pos,pdir,ldir); 
                if((dis=Normal(ldir))>err){
                    NormalizeVec(ldir);
                    coff = tbStrength * a/(fabs(dis)+a);
                    NormalComponent(ldir,xi,gDir);
                    if(DotProduct(vel,gDir)>0){
                        f1[0] += 0.25*length*coff*gDir[0]*burgMag;
                        f1[1] += 0.25*length*coff*gDir[1]*burgMag;
                        f1[2] += 0.25*length*coff*gDir[2]*burgMag;
                    }
                }
            }
        } 
    }

    VECTOR_COPY(f2, f1);

    return;
}

void NodeImageForce(Home_t *home, Node_t *node, real8 *f)
{
    real8 armf1[3],armf2[3];
    VECTOR_ZERO(f); Node_t *nbr;
    for(int i=0; i<node->numNbrs; i++){
        nbr=GetNeighborNode(home, node, i); 
        if(nbr==(Node_t *)NULL)continue;
        ImageForce(home, node, nbr, armf1, armf2);
        f[0] += armf1[0];
        f[1] += armf1[1];
        f[2] += armf1[2];
    }
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
    int     type, armID12, armID21;
    int     numGlideConstraintsNode1, numGlideConstraintsNode2;
    real8   glideConstraintsNode1[3][3], glideConstraintsNode2[3][3]; 
    bool    zeroLenJunc1, zeroLenJunc2;
    real8   p1[3], p2[3], dis, err=home->param->rTol;
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

    Connected(node1, node2, &armID12);
    Connected(node2, node1, &armID21);

    GetGlideConstraintList(home, node1, 
                       &numGlideConstraintsNode1, glideConstraintsNode1, zeroLenJunc1,
                       1, 0.0, 0.0, 0.0, armID12);
    GetGlideConstraintList(home, node2, 
                       &numGlideConstraintsNode2, glideConstraintsNode2, zeroLenJunc2,
                       1, 0.0, 0.0, 0.0, armID21);
    
    switch (numGlideConstraintsNode1){
        case 1:
            switch (numGlideConstraintsNode2){
                case 1:
                    type = PlanePlaneIntersection(glideConstraintsNode1[0], p1,
                           glideConstraintsNode2[0], p2, vec1, vec2);
                    if(type == NO_INTERSECTION){
                        type = fabs(PointPlaneDistance(p1, glideConstraintsNode2[0], p2))<err?PLANE_INTERSECTION:NO_INTERSECTION;
                    }
                    break;
                case 2:
                    type = LinePlaneIntersection(glideConstraintsNode2[2], p2, glideConstraintsNode1[0], p1,
                           vec1, vec2, 0.05, err);
                    if(type == NO_INTERSECTION)
                        type=fabs(PointPlaneDistance(p2, glideConstraintsNode1[0], p1))<err?LINE_INTERSECTION:NO_INTERSECTION;
                    break;
                case 3:
                    PointPlaneIntersection(p2, glideConstraintsNode1[0], p1, vec2, &dis); 
                    type = fabs(dis)<err?POINT_INTERSECTION:NO_INTERSECTION;
                    VECTOR_COPY(vec2, p2);
                    break;
                default:
                    PrintNode(home, node1);
                    PrintNode(home, node2);
                    Fatal("NodesIntersectio: node2 %d",numGlideConstraintsNode2);
                    break;
            }
            break;
        case 2:
            switch (numGlideConstraintsNode2){
                case 1:
                    type = LinePlaneIntersection(glideConstraintsNode1[2], p1, glideConstraintsNode2[0], p2,
                           vec1, vec2, 0.05, err);
                    if(type == NO_INTERSECTION)
                        type=fabs(PointPlaneDistance(p1, glideConstraintsNode2[0], p2))<err?LINE_INTERSECTION:NO_INTERSECTION;
                    break;
                case 2:
                    LineLineIntersection(glideConstraintsNode1[2], p1, glideConstraintsNode2[2], p2,
                           vec1, vec2, &dis, 0.05, err);
                    if(fabs(DotProduct(glideConstraintsNode1[2],glideConstraintsNode2[2]))>0.95){
                        type = fabs(dis)<err?LINE_INTERSECTION:NO_INTERSECTION;
                        VECTOR_COPY(vec1, glideConstraintsNode1[2]);
                    }else{
                        type = fabs(dis)<err?POINT_INTERSECTION:NO_INTERSECTION;
                    }
                    CoordAdd(0.5,p1,0.5,p2,vec2);
                    break;
                case 3:
                    PointLineIntersection(p2, glideConstraintsNode1[2], p1, vec2, &dis); 
                    VECTOR_COPY(vec2, p2);
                    type = fabs(dis)<err?POINT_INTERSECTION:NO_INTERSECTION;
                    break;
                default:
                    PrintNode(home, node1);
                    PrintNode(home, node2);
                    Fatal("NodesIntersectio: node2 %d",numGlideConstraintsNode2);
                    break;
            }

            break;
        case 3:
            switch (numGlideConstraintsNode2){
                case 1:
                     PointPlaneIntersection(p1, glideConstraintsNode2[0], p2,
                           vec2, &dis); 
                    VECTOR_COPY(vec2, p1);
                    type = fabs(dis)<err?POINT_INTERSECTION:NO_INTERSECTION;
                    break;
                case 2:
                    PointLineIntersection(p1, glideConstraintsNode2[2], p2, vec2, &dis); 
                    VECTOR_COPY(vec2, p1);
                    type = fabs(dis)<err?POINT_INTERSECTION:NO_INTERSECTION;
                    break;
                case 3:
                    PointPointIntersection(p2, p1, vec2, &dis); 
                    VECTOR_COPY(vec2, p2);
                    type = fabs(dis)<err?POINT_INTERSECTION:NO_INTERSECTION;
                    break;
                default:
                    PrintNode(home, node1);
                    PrintNode(home, node2);
                    Fatal("NodesIntersectio: node2 %d",numGlideConstraintsNode2);
                    break;
            }
            break;
        default:
            PrintNode(home, node1);
            PrintNode(home, node2);
            Fatal("NodesIntersection: node1 %d",numGlideConstraintsNode1);
            break;
    }

    return(type);
}

bool BeforeMergeNodes(Home_t *home, Node_t *node1, Node_t *node2, real8 *pos, real8* dis, int &type)
{
    bool    canMerge = 1;
    real8   vec1[3]={0.0,0.0,0.0}, vec2[3]={0.0,0.0,0.0};
    real8   vec3[3]={0.0,0.0,0.0}, bakPos[3]={0.0,0.0,0.0};

    *dis = 0.0;
    type = NO_INTERSECTION;

    type = NodesIntersection(home, node1, node2, vec1, vec2); 
    PBCPOSITION(home->param, node1->x, node1->y, node1->z,
                pos, pos+1, pos+2);
    VECTOR_COPY(bakPos, pos);

    switch (type) {
        case NO_INTERSECTION:
            canMerge = 0;
            break;
        case POINT_INTERSECTION:
            if(POINT_INTERSECTION != PointPointIntersection(pos, vec2, vec3, dis)){
                VECTOR_COPY(pos, vec2);
            }
            break;
        case LINE_INTERSECTION:
            if(POINT_INTERSECTION != PointLineIntersection(pos, vec1, vec2, vec3, dis)){
                VECTOR_COPY(pos, vec3);
            }
            break;
        case PLANE_INTERSECTION:
            if(POINT_INTERSECTION != PointPlaneIntersection(pos, vec1, vec2, vec3, dis)){
                VECTOR_COPY(pos, vec3);
            }
            break;
        default:
            Fatal("BeforeMergeArm: wrong type %d", type);
            break;
    }

    CoordMinus(bakPos, pos, bakPos); ZImage(home->param, bakPos, bakPos+1, bakPos+2);
    *dis = Normal(bakPos);

    if(fabs(*dis) > 1000.0 && type != POINT_INTERSECTION){
        int     numGlideConstraintsNode1, numGlideConstraintsNode2;
        real8   glideConstraintsNode1[3][3], glideConstraintsNode2[3][3]; 
        bool    zeroLenJunc1, zeroLenJunc2;
        
        GetGlideConstraintList(home, node1, 
                       &numGlideConstraintsNode1, glideConstraintsNode1, zeroLenJunc1);
        GetGlideConstraintList(home, node2, 
                       &numGlideConstraintsNode2, glideConstraintsNode2, zeroLenJunc2);

        FormatVector(vec1, "vec1");
        FormatVector(vec2, "vec2");
        FormatVector(vec3, "vec3");

        FormatMatrix(glideConstraintsNode1, "numGlideConstraintsNode1");
        FormatMatrix(glideConstraintsNode2, "numGlideConstraintsNode2");

        PrintNode(home, node1);
        PrintNode(home, node2);

        Fatal("%d: Merge dis %f %d (%g,%g,%g)(%f,%f,%f) - %d,%d\n",
                home->myDomain, *dis, bakPos[0],bakPos[1], bakPos[2],type, pos[0],
                pos[1], pos[2],
                numGlideConstraintsNode1, numGlideConstraintsNode2);
    }

    return(canMerge);
}

void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3], int grainID)
{
        int         i;

        burgData->numBurgVectors = NUM_BURGS_IN_FCC;
        burgData->numPlanes   = NUM_PLANES_IN_FCC;
        burgData->numKeyLines = NUM_KEY_LINES_IN_FCC;

        burgData->burgList = (real8 (*)[3])calloc(3*burgData->numBurgVectors, sizeof(real8));
        burgData->burgType = (int *)calloc(burgData->numBurgVectors, sizeof(int));
        burgData->planeList = (real8 (*)[3])calloc(3*burgData->numPlanes, sizeof(real8));
        burgData->planeType = (int *)calloc(burgData->numPlanes, sizeof(int));
        burgData->planeInGrain = (int *)calloc(burgData->numPlanes, sizeof(int));
        burgData->keyLineList = (real8 (*)[3])calloc(3*burgData->numKeyLines, sizeof(real8));

        burgData->numPlanesPerBurg = (int *)calloc(burgData->numBurgVectors, sizeof(int));        
        burgData->numBurgsPerPlane = (int *)calloc(burgData->numPlanes, sizeof(int));        
        burgData->numPlanesPerKeyLine = (int *)calloc(burgData->numKeyLines, sizeof(int));        
        

        real8 burgList[NUM_BURGS_IN_FCC][3] = {
/* 
 *          basic Roman-Roman burgers vectors       (6)
 *          magnitude: 1.0
 */
           { 0.7071067811865475,       0.7071067811865475,        0.0e+00            },      // DC -- 0
           { 0.7071067811865475,      -0.7071067811865475,        0.0e+00            },      // BA -- 1
           { 0.7071067811865475,       0.0e+00,                   0.7071067811865475 },      // DA -- 2
           { 0.7071067811865475,       0.0e+00,                  -0.7071067811865475 },      // BC -- 3
           { 0.0e+00,                  0.7071067811865475,        0.7071067811865475 },      // DB -- 4
           { 0.0e+00,                  0.7071067811865475,       -0.7071067811865475 },      // AC -- 5
/* 
 *          Shocklay Roman-Greek burgers vectors    (12)
 *          magnitude: 1/Sqrt(3)
 */
           { 0.4714045207910317,       0.2357022603955158,        0.2357022603955158 },      // Db -- 6
           {-0.4714045207910317,       0.2357022603955158,        0.2357022603955158 },      // dB -- 7
           { 0.4714045207910317,      -0.2357022603955158,        0.2357022603955158 },      // gA -- 8
           { 0.4714045207910317,       0.2357022603955158,       -0.2357022603955158 },      // aC -- 9
           { 0.2357022603955158,       0.4714045207910317,        0.2357022603955158 },      // Da -- 10
           {-0.2357022603955158,       0.4714045207910317,        0.2357022603955158 },      // gB -- 11
           { 0.2357022603955158,      -0.4714045207910317,        0.2357022603955158 },      // dA -- 12
           { 0.2357022603955158,       0.4714045207910317,       -0.2357022603955158 },      // bC -- 13
           { 0.2357022603955158,       0.2357022603955158,        0.4714045207910317 },      // Dg -- 14
           {-0.2357022603955158,       0.2357022603955158,        0.4714045207910317 },      // aB -- 15
           { 0.2357022603955158,      -0.2357022603955158,        0.4714045207910317 },      // bA -- 16
           { 0.2357022603955158,       0.2357022603955158,       -0.4714045207910317 },      // dC -- 17
/* 
 *          compound Roman-Roman burgers vectors    (3) 
 *          magnitude: Sqrt(2)
 */
           { 1.414213562373095,        0.0e+00,                   0.0e+00            },      // DC+BA -- 18
           { 0.0e+00,                  1.414213562373095,         0.0e+00            },      // DC+AB -- 19
           { 0.0e+00,                  0.0e+00,                   1.414213562373095  }       // CA+DB -- 20
        };

        for(i=0; i<NUM_BURGS_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, burgList[i], burgData->burgList[i]);
        }

        int burgType[NUM_BURGS_IN_FCC] = {EDGE_BURG, EDGE_BURG, EDGE_BURG, EDGE_BURG, EDGE_BURG, EDGE_BURG,
                PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, 
                PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG, PARTIAL_BURG,
                JUNCTION_BURG, JUNCTION_BURG, JUNCTION_BURG};
        memcpy(burgData->burgType, burgType, sizeof(int [NUM_BURGS_IN_FCC]));


        int numPlanesPerBurg[NUM_BURGS_IN_FCC] = { 3, 3, 3, 3, 3, 3, 
                                                   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                   4, 4, 4 }; 
        memcpy(burgData->numPlanesPerBurg,numPlanesPerBurg,sizeof(int [NUM_BURGS_IN_FCC]));

        int burgFirstPlaneIndex[NUM_BURGS_IN_FCC] = { 0, 3, 6, 9, 12, 15, 
                                                     18,20,22,24,26,28,30,32,34,36,38,40,
                                                     42,46,50};

        real8 planeList[NUM_PLANES_IN_FCC][3] = {
            {-0.5773502691896258,       0.5773502691896258,       -0.5773502691896258 },    // a plane  ID: 0
            { 0.5773502691896258,      -0.5773502691896258,       -0.5773502691896258 },    // b plane  ID: 1
            {-0.5773502691896258,      -0.5773502691896258,        0.5773502691896258 },    // g plane  ID: 2
            { 0.5773502691896258,       0.5773502691896258,        0.5773502691896258 },    // d plane  ID: 3
                                                                                      
            { 1.0e+00,                  0.0e+00,                   0.0e+00            },    // (100) junction plane     ID: 4
            { 0.0e+00,                  1.0e+00,                   0.0e+00            },    // (010) junction plane     ID: 5
            { 0.0e+00,                  0.0e+00,                   1.0e+00            },    // (001) junction plane     ID: 6
                                                                                      
            { 0.0e+00,                  0.7071067811865475,        0.7071067811865475 },    // (011) junction plane     ID: 7
            { 0.0e+00,                 -0.7071067811865475,        0.7071067811865475 },    // (0-11) junction plane    ID: 8
                                                                                      
            { 0.7071067811865475,       0.0e+00,                   0.7071067811865475 },    // (101) junction plane     ID: 9
            {-0.7071067811865475,       0.0e+00,                   0.7071067811865475 },    // (-101) junction plane    ID: 10
                                                                                      
            { 0.7071067811865475,       0.7071067811865475,        0.0e+00            },    // (110) junction plane     ID: 11
            {-0.7071067811865475,       0.7071067811865475,        0.0e+00            }     // (-110) junction plane    ID: 12
        };

        for(i=0; i<NUM_PLANES_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, planeList[i], burgData->planeList[i]);
        }
        
        int planeType[NUM_PLANES_IN_FCC] = { GLIDE_PLANE, GLIDE_PLANE, GLIDE_PLANE, GLIDE_PLANE,
            JUNCTION_PLANE1, JUNCTION_PLANE1, JUNCTION_PLANE1, JUNCTION_PLANE2, 
            JUNCTION_PLANE2, JUNCTION_PLANE2, JUNCTION_PLANE2, JUNCTION_PLANE2, JUNCTION_PLANE2 }; 

        memcpy(burgData->planeType,planeType,NUM_PLANES_IN_FCC*sizeof(int));

        int planeInGrain[NUM_PLANES_IN_FCC] = { grainID, grainID, grainID, grainID,                          
                                                grainID, grainID, grainID, grainID, 
                                                grainID, grainID, grainID, grainID, grainID};     
        memcpy(burgData->planeInGrain,planeInGrain,NUM_PLANES_IN_FCC*sizeof(int));

        int numBurgsPerPlane[NUM_PLANES_IN_FCC] = { 6, 6, 6, 6,                    // basic glide planes
                                                    4, 4, 4,                       // Junction plane 1 
                                                    3, 3, 3, 3, 3, 3               // Junction plane 2
                                                  }; 

        memcpy(burgData->numBurgsPerPlane, numBurgsPerPlane,sizeof(int [NUM_PLANES_IN_FCC]));

        int planeFirstBurgIndex[NUM_PLANES_IN_FCC] = { 0, 6, 12, 18,
                                                       24,28,32,
                                                       36,39,42,45,48,51 };  

        real8 keyLineList[NUM_KEY_LINES_IN_FCC][3] = {
           { 0.7071067811865475,       0.7071067811865475,        0.0e+00            },      // DC -- 0
           { 0.7071067811865475,      -0.7071067811865475,        0.0e+00            },      // BA -- 1
           { 0.7071067811865475,       0.0e+00,                   0.7071067811865475 },      // DA -- 2
           { 0.7071067811865475,       0.0e+00,                  -0.7071067811865475 },      // BC -- 3
           { 0.0e+00,                  0.7071067811865475,        0.7071067811865475 },      // DB -- 4
           { 0.0e+00,                  0.7071067811865475,       -0.7071067811865475 }       // AC -- 5
        };

        for(i=0; i<NUM_KEY_LINES_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, keyLineList[i], burgData->keyLineList[i]);
        }

        int numPlanesPerKeyLine[NUM_KEY_LINES_IN_FCC] = { 3, 3, 3, 3, 3, 3};
 
        memcpy(burgData->numPlanesPerKeyLine,numPlanesPerKeyLine,sizeof(int [NUM_KEY_LINES_IN_FCC]));

        int keyLineFirstPlaneIndex[NUM_BURGS_IN_FCC] = { 0, 3, 6, 9, 12, 15 }; 

        int burgOnPlaneID[54] = {
                                  0,  1,  6,       /* DC (0) on plane  */
                                  2,  3,  6,       /* BA (1) on plane  */
                                  1,  2,  5,       /* DA (2) on plane  */
                                  0,  3,  5,       /* BC (3) on plane  */
                                  0,  2,  4,       /* DB (4) on plane  */
                                  1,  3,  4,       /* AC (5) on plane  */
                                  1,  8,           /* Db (6) */
                                  3,  8,           /* dB (7) */
                                  2,  7,           /* gA (8) */
                                  0,  7,           /* aC (9) */
                                  0, 10,           /* Da (10) */
                                  2,  9,           /* gB (11) */
                                  3, 10,           /* dA (12) */
                                  1,  9,           /* bC (13) */
                                  2, 12,           /* Dg (14) */
                                  0, 11,           /* aB (15) */
                                  1, 11,           /* bA (16) */
                                  3, 12,           /* dC (17) */
                                  5,  6,  7,  8,   /* [100] (18) */
                                  4,  6,  9, 10,   /* [010] (19) */
                                  4,  5, 11, 12    /* [001] (20) */
                                };
        burgData->burgOnPlaneID = (int **)calloc(NUM_BURGS_IN_FCC, sizeof(int *));       
        for(i=0; i<NUM_BURGS_IN_FCC; i++){
            burgData->burgOnPlaneID[i] = (int *)calloc(numPlanesPerBurg[i], sizeof(int));
            memcpy(burgData->burgOnPlaneID[i], burgOnPlaneID+burgFirstPlaneIndex[i], 
                   sizeof(int [numPlanesPerBurg[i]]));
        }
        
        int planeOnBurgID[54] = {  0,  3,  4,  9, 10, 15,     // plane (0)
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

        burgData->planeOnBurgID = (int **)calloc(NUM_PLANES_IN_FCC, sizeof(int *));       
        for(i=0; i<NUM_PLANES_IN_FCC; i++){
            burgData->planeOnBurgID[i] = (int *)calloc(numBurgsPerPlane[i], sizeof(int));
            memcpy(burgData->planeOnBurgID[i], planeOnBurgID+planeFirstBurgIndex[i], 
                   sizeof(int [numBurgsPerPlane[i]]));
        }

        int keyLineOnPlaneID[18] = {
                                     0, 1, 6,
                                     2, 3, 6,
                                     1, 2, 5,
                                     0, 3, 5,
                                     0, 2, 4,
                                     1, 3, 4
                                   };
        burgData->keyLineOnPlaneID = (int **)calloc(NUM_KEY_LINES_IN_FCC, sizeof(int *));       
        for(i=0; i<NUM_KEY_LINES_IN_FCC; i++){
            burgData->keyLineOnPlaneID[i] = (int *)calloc(numPlanesPerKeyLine[i], sizeof(int));
            memcpy(burgData->keyLineOnPlaneID[i], keyLineOnPlaneID+keyLineFirstPlaneIndex[i], 
                   sizeof(int [numPlanesPerKeyLine[i]]));
        }
        return;
}

void FreeBurgInfoStr(BurgInfo_t *burgData)
{
    
    int     i;

    free(burgData->burgList); 
    free(burgData->planeList); 
    free(burgData->planeType); 
    free(burgData->burgType); 
    free(burgData->planeInGrain); 
    free(burgData->keyLineList);

    if(burgData->burgOnPlaneID != (int **)NULL){
        for(i=0; i<burgData->numBurgVectors; i++){
            free(*(burgData->burgOnPlaneID+i)); 
        }
    }

    if(burgData->planeOnBurgID != (int **)NULL){
        for(i=0; i<burgData->numPlanes; i++){
            free(*(burgData->planeOnBurgID+i));
        }
    }

    if(burgData->keyLineOnPlaneID != (int **)NULL){
        for(i=0; i<burgData->numKeyLines; i++){
            free(*(burgData->keyLineOnPlaneID+i));
        }
    }

    if(burgData->numPlanesPerBurg != (int *)NULL){
        free(burgData->numPlanesPerBurg); 
    }

    if(burgData->numBurgsPerPlane != (int *)NULL){
        free(burgData->numBurgsPerPlane); 
    }

    if(burgData->numPlanesPerKeyLine != (int *)NULL){
        free(burgData->numPlanesPerKeyLine); 
    }

    burgData->burgList = NULL;
    burgData->planeList = NULL;
    burgData->planeType = NULL;
    burgData->burgType = NULL;
    burgData->planeInGrain = NULL;
    burgData->keyLineList = NULL;

    burgData->numPlanesPerBurg = NULL;
    burgData->numPlanesPerKeyLine = NULL;
    burgData->numBurgsPerPlane = NULL;

    burgData->burgOnPlaneID = NULL;
    burgData->planeOnBurgID = NULL;
    burgData->keyLineOnPlaneID = NULL;

    return;
}

void CompressDirectionList(std::vector< std::vector<real8> >& dirList)
{
    int     i, j, nRows;
    real8   crossProdt[3], vec1[3], vec2[3];
    int     newNRows = 0;
    std::vector< std::vector<real8> >::iterator  it;

    nRows = dirList.size();
    for(i=0; i<nRows; i++){
        for(j=0; j<newNRows; j++){
            VECTOR_COPY(vec1, dirList[i]);
            VECTOR_COPY(vec2, dirList[j]);
            cross(vec1, vec2, crossProdt);
            if(Normal(crossProdt) < 1.0E-2)break;

        }
        if(j == newNRows){
            VECTOR_COPY(dirList[newNRows], dirList[i])
            newNRows++;
        }
    }
    
    for(it = dirList.end(); it != dirList.begin()+newNRows; it--){
        dirList.erase(it);
    }
    return;
}

void CompressBurgListWithType(std::vector<std::vector<real8> >& burgList,
                              std::vector<int>& burgType)
{
    int     i, j, nRows;
    real8   crossProdt[3], vec1[3], vec2[3];
    int     newNRows = 0;
    std::vector< std::vector<real8> >::iterator  it;
    std::vector<int>::iterator  it2;

    nRows = burgList.size();
    if(nRows != burgType.size())Fatal("in %s at %d", __FILE__, __LINE__);
    for(i=0; i<nRows; i++){
        for(j=0; j<newNRows; j++){
            VECTOR_COPY(vec1, burgList[i]);
            VECTOR_COPY(vec2, burgList[j]);
            cross(vec1, vec2, crossProdt);
            if(Normal(crossProdt) < 1.0E-2 && fabs(Normal(burgList[i]) - 
               Normal(burgList[j])) < 1.0E-4 &&  burgType[i] == burgType[j])break;

        }
        if(j == newNRows){
            VECTOR_COPY(burgList[newNRows], burgList[i])
            burgType[newNRows] = burgType[i];
            newNRows++;
        }
    }
    
    for(it = burgList.end(); it != burgList.begin()+newNRows; it--){
        burgList.erase(it);
    }
//    for(it2 = burgType.end(); it2 != burgType.begin()+newNRows; it2--){
//        burgType.erase(it2);
//    }
        burgType.resize(newNRows);

    return;
}
void GenerateGBBurgInfo(Home_t *home, GB_t *GB)
{

    int         i, j, numRows;
    real8       crossProdt[3], vec[3];
    Grain_t     *grain1, *grain2;
    BurgInfo_t  *bd1, *bd2, *gbd;
    std::vector<double>                 vecSTL(3);
   
    grain1 = home->grainKeys[GB->nbrGrain[0]];
    grain2 = home->grainKeys[GB->nbrGrain[1]];

    bd1 = &(grain1->burgData);
    bd2 = &(grain2->burgData);
    gbd = &(GB->burgData);

    numRows = bd1->numPlanes + bd2->numPlanes;
    gbd->numPlanes = numRows;

    gbd->planeList = (real8 (*)[3])malloc(3*numRows*sizeof(real8));
    gbd->planeType = (int *)malloc(numRows*sizeof(int));
    gbd->planeInGrain = (int *)malloc(numRows*sizeof(int));

    memcpy(gbd->planeList, bd1->planeList, bd1->numPlanes*sizeof(real8 [3]));
    memcpy(gbd->planeList+bd1->numPlanes, bd2->planeList, 
           bd2->numPlanes*sizeof(real8 [3]));

    memcpy(gbd->planeType, bd1->planeType, bd1->numPlanes*sizeof(int));
    memcpy(gbd->planeType+bd1->numPlanes, bd2->planeType, 
           bd2->numPlanes*sizeof(int));

    memcpy(gbd->planeInGrain, bd1->planeInGrain, bd1->numPlanes*sizeof(int));
    memcpy(gbd->planeInGrain+bd1->numPlanes, bd2->planeInGrain, 
           bd2->numPlanes*sizeof(int));

    for(i=0; i<numRows; i++){
        cross(gbd->planeList[i], GB->n, crossProdt);
        if(Normal(crossProdt) < 1.0E-2){
            gbd->planeType[i] |= GB_PLANE;
            gbd->planeInGrain[i] = -1;
        }
    }

    std::vector<std::vector<real8> >                burgList;
    std::vector<int>                                burgType;
    for(i=0; i<bd1->numBurgVectors; i++){
        VECTOR_COPY(vecSTL, bd1->burgList[i]);
        burgList.push_back(vecSTL);
        burgType.push_back(bd1->burgType[i]);
    }
    for(i=0; i<bd2->numBurgVectors; i++){
        VECTOR_COPY(vecSTL, bd2->burgList[i]);
        burgList.push_back(vecSTL);
        burgType.push_back(bd2->burgType[i]);
    }
    CompressBurgListWithType(burgList, burgType);

    numRows = burgList.size();
    gbd->numBurgVectors = numRows;
    gbd->burgList = (real8 (*)[3])malloc(3*numRows*sizeof(real8));
    gbd->burgType = (int *)malloc(numRows*sizeof(int));
    gbd->numPlanesPerBurg = (int *)malloc(numRows*sizeof(int));
    std::vector<std::vector<int> >        burgOnPlaneID(numRows);

    for(i=0; i<numRows; i++){
        VECTOR_COPY(gbd->burgList[i], burgList[i]);
        gbd->burgType[i] = burgType[i];
        gbd->numPlanesPerBurg[i] = 0;
        VECTOR_COPY(vec, burgList[i]);
        NormalizeVec(vec);
        for(j=0; j<gbd->numPlanes; j++){
            if(fabs(DotProduct(vec, gbd->planeList[j])) < 1.0E-2){
                gbd->numPlanesPerBurg[i]++;
                burgOnPlaneID[i].push_back(j);
            }
        }
    }
    
    std::vector<std::vector<int> >      planeOnBurgID(gbd->numPlanes);
    gbd->burgOnPlaneID = (int **)malloc(gbd->numBurgVectors*sizeof(int *));
    for(i=0; i<gbd->numBurgVectors; i++){
        gbd->burgOnPlaneID[i] = (int *)malloc(gbd->numPlanesPerBurg[i]*sizeof(int));
        for(j=0; j<gbd->numPlanesPerBurg[i]; j++){
            gbd->burgOnPlaneID[i][j] = burgOnPlaneID[i][j];
            planeOnBurgID[burgOnPlaneID[i][j]].push_back(i);
        }
    }

    gbd->numBurgsPerPlane = (int *)malloc(gbd->numPlanes*sizeof(int));
    gbd->planeOnBurgID = (int **)malloc(gbd->numPlanes*sizeof(int *));
    for(i=0; i<gbd->numPlanes; i++){
        gbd->numBurgsPerPlane[i] = planeOnBurgID[i].size();
        gbd->planeOnBurgID[i] = (int *)malloc(gbd->numBurgsPerPlane[i]*sizeof(int));
        for(j=0; j<gbd->numBurgsPerPlane[i]; j++){
            gbd->planeOnBurgID[i][j] = planeOnBurgID[i][j];
        }
    }
    
    std::vector<std::vector<double> >   keyLineList;

    for(i=0; i<gbd->numPlanes; i++){
        cross(gbd->planeList[i], GB->n, crossProdt);
        if(Normal(crossProdt) < 1E-2)continue;
        if((gbd->planeType[i] & GLIDE_PLANE) == 0)continue;
        NormalizeVec(crossProdt);
        VECTOR_COPY(vecSTL, crossProdt);
        keyLineList.push_back(vecSTL);
    }

    CompressDirectionList(keyLineList);

    gbd->numKeyLines = keyLineList.size();
    gbd->keyLineList = (real8 (*)[3])malloc(3*gbd->numKeyLines*sizeof(real8));

    for(i=0; i<keyLineList.size(); i++){
        VECTOR_COPY(gbd->keyLineList[i], keyLineList[i]);
    }
    
#if 0
    printf("Burgs info:\n");
    for(i=0; i<gbd->numBurgVectors; i++){
        printf("%d(%d): {%f,%f,%f} numPlanes %d, ", i,gbd->burgType[i], 
               gbd->burgList[i][0], gbd->burgList[i][1],gbd->burgList[i][2],
               gbd->numPlanesPerBurg[i]);
        for(j=0; j<gbd->numPlanesPerBurg[i]; j++){
            printf("%d ", gbd->burgOnPlaneID[i][j]);
        }
        printf("\n");
    }

    printf("\nPlanes info:\n");
    for(i=0; i<gbd->numPlanes; i++){
        printf("%d(%d,%d): {%f,%f,%f} numBurgs %d, ", i ,gbd->planeInGrain[i], gbd->planeType[i],
               gbd->planeList[i][0], gbd->planeList[i][1],gbd->planeList[i][2],
               gbd->numBurgsPerPlane[i]);
        for(j=0; j<gbd->numBurgsPerPlane[i]; j++){
            printf("%d ", gbd->planeOnBurgID[i][j]);
        }
        printf("\n");
        
    }
#endif
    return;
}


void SumOfValuesIsConstant(int n, int sum, std::vector<int>& valueList, int tot, 
                            std::vector<std::vector<int> >& valueLists, int maxVar)
{
    int i, j, nextSum;

    if(n == 0)return;
    for(i=0; i<= maxVar && i<= sum; i++){
        valueList[n-1] = i;
        nextSum = sum - i;
        if(nextSum > (n-1)*maxVar)continue;
        SumOfValuesIsConstant(n-1, nextSum, valueList, tot
                                 ,valueLists, maxVar);
    }

    if(n-1 == 0 && nextSum <= maxVar){
        valueLists.push_back(valueList);
    }
    return;
}


bool cmp(std::vector<int > p, std::vector<int > q)
{
    if(p.back() != q.back()){
        return(p.back()<q.back());
    }else{
        int sum1 = accumulate(p.begin(), p.end()-1, 0);
        int sum2 = accumulate(q.begin(), q.end()-1, 0);
        return(sum1 < sum2);
    }
}

bool cmp2(Dislocation_t p, Dislocation_t q)
{
    return(p.planeID < q.planeID);
}

void Iter(std::vector<int>& disLocalBurgID, 
          std::vector< std::vector<int> >& burgOnPlaneID,
          std::vector<int>& planeIDList, 
          std::vector<int>& lastPlaneIDList,
          std::vector< std::vector<int> >& planeIDLists,
          int i, int j)
{
    int m,n;
    
    planeIDList[i] = burgOnPlaneID[disLocalBurgID[i]][j];

    if(i == disLocalBurgID.size() - 1){
        for(m=0; m<planeIDList.size(); m++){
            lastPlaneIDList[m] = planeIDList[m];
        }
        planeIDLists.push_back(planeIDList);
        if(j == burgOnPlaneID[disLocalBurgID[i]].size() -1){
#if DEBUG_DisassembleDislocation
            printf("Iter:\n");
            for(m=0; m<planeIDLists.size(); m++){
                for(n=0; n<planeIDLists[m].size(); n++){
                    printf("(%d)[%d] ", planeIDLists[m][n],
                            disLocalBurgID[n]);
                }
                printf("\n");
            }
#endif
            return;
        }else{
            i=0;
        }
    }else{
        i++;
    }
    
    for(m=0; m<burgOnPlaneID[disLocalBurgID[i]].size(); m++){
        if(lastPlaneIDList[i] == burgOnPlaneID[disLocalBurgID[i]][m]){
            j = m+1;
            break;
        }
    }
    if(m >= burgOnPlaneID[disLocalBurgID[i]].size() - 1)j=0;

    Iter(disLocalBurgID, burgOnPlaneID, planeIDList,
          lastPlaneIDList, planeIDLists, i, j);
    return;
}


void CompressVectorList(std::vector< std::vector<int> >& list)
{
    int i, j, k, nRows = 0;
    for(i=0; i<list.size(); i++){
        for(j=0; j<nRows; j++){
            if(list[i].size() == list[j].size()){
                for(k=0; k<list[i].size(); k++){
                    if(list[i][k] != list[j][k])break;
                }
                if(k == list[i].size())break;
            }
            
        }
        if(j == nRows){
            for(k=0; k<list[nRows].size(); k++){
                list[nRows][k] = list[i][k];
            }
            nRows++;
        }
    }
    list.resize(nRows);
}

void CompressDislocationArray(std::vector< std::vector<Dislocation_t> >& list)
{
    int i, j, k, nRows = 0;
    for(i=0; i<list.size(); i++){
        for(j=0; j<nRows; j++){
            if(list[i].size() == list[j].size()){
                for(k=0; k<list[i].size(); k++){
                    if(list[i][k].planeID != list[j][k].planeID ||
                    list[i][k].burgID != list[j][k].burgID)break;
                }
                if(k == list[i].size())break;
            }
            
        }
        if(j == nRows){
            for(k=0; k<list[nRows].size(); k++){
                list[nRows][k].planeID = list[i][k].planeID;
                list[nRows][k].burgID = list[i][k].burgID;
            }
            nRows++;
        }
    }
    list.resize(nRows);
}
void IterFun1(int n, int m, std::vector<int>& list,
       std::vector<std::vector<int> >& lists)
{
    if(m < 1 || n < 0)return;
    int i;

    
    for(i=0; i< n+1; i++){
        list[m-1] = i;
        IterFun1(n-i,m-1,list,lists);
    }
    if(m==1){
        list[0] = n; 
        std::vector<int>  IDList(list.size());
        std::vector<int>  sortList(list.size());
        for(i=0; i<list.size(); i++){
            IDList[i] = i;
        }
        while(next_permutation(IDList.begin(), IDList.end())){
            for(i=0; i<list.size(); i++){
                sortList[i] = list[IDList[i]];
            }
            lists.push_back(sortList);
        }

        if(list.size() == 1){
            std::vector<int> a;
            a.push_back(n);
            lists.push_back(a);
        }
    }

    return;
}


void IterFun2(std::vector<int>& oriVec, int iLayer, 
             std::vector<int>& list, 
             std::vector< std::vector<int> >& lists)
{
    int i;
    if(iLayer == oriVec.size()){
        lists.push_back(list);
        return;
    }

    for(i=0; i<oriVec[iLayer]; i++){
        list[iLayer] = i;
        IterFun2(oriVec, iLayer+1, list, lists);
    }

    return;
}


int GeometricCriterion(BurgInfo_t *bd, real8 *burg, real8 *lDir,
                       real8 *plane, int flag, 
                       std::vector< std::vector<double> >&       burgList, 
                       std::vector<std::vector<Dislocation_t> >& disArray,
                       int maxDissPerBasicBurg, int maxDiss) 
{
    int     i, j, k, m, nPlanes = 0, nBurgs = 0;
    int     sum, onGBPlane;
    int     nBurgSets, burgType;
    real8   vec[3];
    real8   crossProd[3], resiBurg[3];
    real8   resiBurgMag;
    std::vector< std::vector<double> >              planeList;
    std::vector<int>                                planeIDList;
    std::vector< std::vector<int> >                 burgOnPlaneID;
    std::vector< std::vector<int> >                 setLists;
    std::vector<int>                                valueList;
    std::vector< std::vector<int> >::iterator       iter;
    std::vector<int>                                planeID;
    std::vector<double>                             vector(3);
    std::vector<int>                                disOnBurgID;
    std::vector<int>                                disOnPlaneID;
    std::vector<int>                                oriVec;
    std::vector<int>                                orderList;
    std::vector<std::vector<int> >                  disOnPlaneIDList;
    std::vector<std::vector<int> >                  disOnBurgIDList;
    std::vector<std::vector<int> >                  orderArray;

    std::vector<std::vector<int> >                  recordList(3);
    std::vector<std::vector<int> >                  recordArray;
    
    std::vector<int>                                sortList;
    std::vector<Dislocation_t>                      disList;
    std::vector<Dislocation_t>                      tempDisList;
    std::vector<std::vector<Dislocation_t> >        tempDisArray;
    
    std::vector<std::vector<std::vector<int> > >    disOnPlaneIDArray;

#if DEBUG_DisassembleDislocation 
    printf("Burg Vector: {%7.4f,%7.4f,%7.4f}\n", burg[0],burg[1],burg[2]);
    printf("Line Vector: {%7.4f,%7.4f,%7.4f}\n", lDir[0],lDir[1],lDir[2]);
    if(flag == 1){
        printf("origin Plane Normal cross-slip: {%7.4f,%7.4f,%7.4f}\n\n", plane[0],plane[1],plane[2]);
    }else if(flag == 2){
        printf("origin Plane Normal decompose: {%7.4f,%7.4f,%7.4f}\n\n", plane[0],plane[1],plane[2]);
    }
    printf("Basic glide planes : \n");
#endif

    std::vector<std::vector<double> >().swap(burgList);
    std::vector<std::vector<Dislocation_t> >().swap(disArray);

    for(i=0; i<bd->numPlanes; i++){
        planeID.push_back(nPlanes);
        if((bd->planeType[i] & GLIDE_PLANE) == 0)continue;
        if(fabs(DotProduct(lDir,bd->planeList[i])) < 1.0E-2){
            cross(bd->planeList[i], plane, crossProd);
            switch(flag){
                case 0:
                    VECTOR_COPY(vector, bd->planeList[i]);
                    planeList.push_back(vector);
                    planeIDList.push_back(i);
                    nPlanes++;
                    break;
                case 1:
                    if(Normal(crossProd) > 1.0E-2){
                        VECTOR_COPY(vector, bd->planeList[i]);
                        planeList.push_back(vector);
                        planeIDList.push_back(i);
                        nPlanes++;
                    }
                    break;
                case 2:
                    if(Normal(crossProd) < 1.0E-2){
                        VECTOR_COPY(vector, bd->planeList[i]);
                        planeList.push_back(vector);
                        planeIDList.push_back(i);
                        nPlanes++;
                    }
                    break;
                default:
                    break;
            }
        }
    }

#if DEBUG_DisassembleDislocation
    for(i=0; i<planeIDList.size(); i++){
        printf("%d(%d): {%7.4f,%7.4f,%7.4f}\n",
               i, planeIDList[i], planeList[i][0],
               planeList[i][1], planeList[i][2]);
    }
#endif

    j = 0;
    nBurgs = 0;
    for(i=0; i<bd->numPlanes; i++){
        if(j == nPlanes)continue;
        if(i == planeIDList[j]){
            j++;
            for(k=0; k<bd->numBurgsPerPlane[i]; k++){
                burgType = BurgType_FCC(bd->burgList[bd->planeOnBurgID[i][k]]);
                if(flag == 2 && burgType == BASIC_BURG)continue;
                if(!(burgType != SHOCKLEY_BURG || burgType != BASIC_BURG))continue;
                VECTOR_COPY(vector, bd->burgList[bd->planeOnBurgID[i][k]]);
                if(DotProduct(vector, burg) < 0){
                    vector[0] = -vector[0];
                    vector[1] = -vector[1];
                    vector[2] = -vector[2];
                }
                burgList.push_back(vector);
                nBurgs++;
            }
        }
    }
    
    CompressDirectionList(burgList);

    if(burgList.size() == 0)return(0); 

#if DEBUG_DisassembleDislocation
    printf("\nBasic Shockley burgs: \n");
    for(i=0; i<burgList.size(); i++){
        printf("%2d: {%7.4f,%7.4f,%7.4f}\n",
               i, burgList[i][0], burgList[i][1],
               burgList[i][2]);
    }
    printf("\n");
#endif
    
    valueList.resize(burgList.size());
    for(sum=0; sum < maxDiss+1; sum++){
        SumOfValuesIsConstant(burgList.size(), sum, valueList, 
                              burgList.size(), setLists, maxDissPerBasicBurg);
    }

    if(setLists.size() == 0)return(0);
    
    for(i=0; i<setLists.size(); i++){
        resiBurg[0] = burg[0];
        resiBurg[1] = burg[1];
        resiBurg[2] = burg[2];

        for(j=0; j<setLists[i].size(); j++){
            resiBurg[0] -= ((real8)setLists[i][j] * burgList[j][0]); 
            resiBurg[1] -= ((real8)setLists[i][j] * burgList[j][1]); 
            resiBurg[2] -= ((real8)setLists[i][j] * burgList[j][2]); 
        }
        
        resiBurgMag = Normal(resiBurg);
        for(j=0; j<setLists[i].size(); j++){
            VECTOR_COPY(vec, burgList[j]);
            cross(resiBurg, vec, crossProd);
            if(Normal(crossProd) < 1.0E-2 && 
               BurgType_FCC(resiBurg) == SHOCKLEY_BURG){
                break;    
            }
        }

        if(j != setLists[i].size() && resiBurgMag > 1.0E-2){
            setLists[i].push_back(100);
        }else{
            setLists[i].push_back((int)(resiBurgMag*100));
        }
    }

    sort(setLists.begin(), setLists.end(),cmp); 
    
    nBurgSets = 0;
    for(i=0; i<setLists.size(); i++){
        if(setLists[i][setLists[i].size()-1] > 
           (int)(SHOCKLEY_BURG_MAG_FCC * 100 + 1.0)){
            break;
        }
#if 1
        if(setLists[0][setLists[i].size()-1] == 0){
            if(setLists[i][setLists[i].size()-1] > 0){
                break;
            }
        } 
#endif

        nBurgSets++;
    }
    for(iter = setLists.end(); iter != setLists.begin()+nBurgSets; iter--){
        setLists.erase(iter);
    }

    if(nBurgSets == 0)return(0);
#if DEBUG_DisassembleDislocation
    printf("SetLists: \n");
    for(i=0; i<setLists.size(); i++){
        printf("Set %2d: ", i);
        for(j=0; j<setLists[i].size(); j++){
            printf("%4d ", setLists[i][j]);
        }
        printf("\n");
    }
    printf("\n");
#endif

    burgOnPlaneID.resize(burgList.size());
    for(i=0; i<burgOnPlaneID.size(); i++){
        onGBPlane = 0;
        for(j=0; j<nPlanes; j++){
            if(fabs(DotProduct(burgList[i], 
                planeList[j])) < 1.0E-2){
                if((bd->planeType[planeIDList[j]] & GB_PLANE) > 0){
                    onGBPlane ++;
                    if(onGBPlane == 2)continue;
                }
                if((bd->planeType[planeIDList[j]] & GLIDE_PLANE) == 0)continue;
                burgOnPlaneID[i].push_back(planeIDList[j]);
            }
        }
    }

#if DEBUG_DisassembleDislocation
//    burgOnPlaneID[3].push_back(15);
//    burgOnPlaneID[7].push_back(3);
//    printf("burgOnPlaneID[2].size() %d\n", burgOnPlaneID[2].size());
    for(i=0; i<burgOnPlaneID.size(); i++){
        printf("Shockley burg %2d on plane:", i);
        for(j=0; j<burgOnPlaneID[i].size(); j++){
            printf(" %2d(%d)--{%7.4f,%7.4f,%7.4f} ", planeID[burgOnPlaneID[i][j]],
                   burgOnPlaneID[i][j],
                   bd->planeList[burgOnPlaneID[i][j]][0],
                   bd->planeList[burgOnPlaneID[i][j]][1],
                   bd->planeList[burgOnPlaneID[i][j]][2]);
        }
        printf("\n");
    }
    printf("\n");
#endif
    for(i=0; i<setLists.size(); i++){
        std::vector<int>().swap(disOnBurgID);
        for(j=0; j<burgList.size(); j++){
            for(k=0; k<setLists[i][j]; k++){
                disOnBurgID.push_back(j);
            }
        }
        disOnBurgIDList.push_back(disOnBurgID);
    }

    for(i=0; i<burgOnPlaneID.size(); i++){
        for(j=0; j<burgOnPlaneID[i].size(); j++){
            recordList[0].push_back(i);
            recordList[1].push_back(j);
            recordList[2].push_back(0);
        }
    }

    for(i=0; i<disOnBurgIDList.size(); i++){
        disOnPlaneID.resize(disOnBurgIDList[i].size());
        oriVec.resize(disOnBurgIDList[i].size());
        orderList.resize(disOnBurgIDList[i].size());
        std::vector<std::vector<int> > ().swap(orderArray);
        for(j=0; j<disOnBurgIDList[i].size(); j++){
            oriVec[j] = burgOnPlaneID[disOnBurgIDList[i][j]].size();
        }

        IterFun2(oriVec, 0, orderList, orderArray);        

        std::vector<std::vector<int> >().swap(disOnPlaneIDList);
        disOnPlaneID.resize(disOnBurgIDList[i].size());
        std::vector<std::vector<int> >().swap(recordArray);
        for(k=0; k<orderArray.size(); k++){
            for(m=0; m<recordList[0].size(); m++){
                recordList[2][m] = 0;           
            }

            for(j=0; j<disOnBurgIDList[i].size(); j++){
                for(m=0; m<recordList[0].size(); m++){
                    if(disOnBurgIDList[i][j] == recordList[0][m] &&
                       orderArray[k][j] == recordList[1][m]){
                        recordList[2][m]++;
                    }
                }
                disOnPlaneID[j] = burgOnPlaneID[disOnBurgIDList[i][j]][orderArray[k][j]];
            }


            if(recordArray.size() == 0){
                recordArray.push_back(recordList[2]);
                disOnPlaneIDList.push_back(disOnPlaneID);
            }else{
                for(j=0; j<recordArray.size(); j++){
                    for(m=0; m<recordList[0].size(); m++){
                        if(recordList[2][m] != recordArray[j][m])break;
                    }
                    if(m == recordList[0].size())break;
                }
                if(j == recordArray.size()){
                    recordArray.push_back(recordList[2]);
                    disOnPlaneIDList.push_back(disOnPlaneID);
                }
            }
        }
        disOnPlaneIDArray.push_back(disOnPlaneIDList);
    }

    for(i=0; i<disOnPlaneIDArray.size(); i++){
        for(j=0; j<disOnPlaneIDArray[i].size(); j++){
            disList.resize(disOnPlaneIDArray[i][j].size());
            for(k=0; k<disOnPlaneIDArray[i][j].size(); k++){
                disList[k].burgID = disOnBurgIDList[i][k];
                disList[k].planeID = disOnPlaneIDArray[i][j][k];
                VECTOR_COPY(vec, burgList[disList[k].burgID]);
                disList[k].energyFactor = Normal(vec)/Normal(burg);
                disList[k].energyFactor = disList[k].energyFactor*
                                            disList[k].energyFactor;
            }

            sortList.resize(disList.size());
            tempDisList.resize(disList.size());

            std::vector<std::vector<Dislocation_t> >().swap(tempDisArray);

            for(k=0; k<sortList.size(); k++)sortList[k] = k;
            
            tempDisArray.push_back(disList);
#if 0
            while(next_permutation(sortList.begin(), sortList.end())){
                for(k=0; k<sortList.size(); k++){
                    tempDisList[k] = disList[sortList[k]];
                }
                
                sort(tempDisList.begin(), tempDisList.end(), cmp2);
                tempDisArray.push_back(tempDisList);
            }


            CompressDislocationArray(tempDisArray);
#endif
            for(k=0; k<tempDisArray.size(); k++){
                disArray.push_back(tempDisArray[k]);
            }
        }
    }

    return(1);
}

void FindGlidePlane(Home_t *home, real8 burg[3], real8 lDir[3],
                    real8 glidePlane[3], Node_t *node1, Node_t *node2)
{
    int         armID, burgID, normalID=-1, lineID;
    int         inGrain, onGB;
    int         burgType, interType;
    int         times;
    real8       normBurg[3], crossProdt[3];
    real8       vec1[3], vec2[3];
    Param_t     *param = home->param;
    BurgInfo_t  *bd = GetBurgData(home, node1, node2);

#if 1
#else
    FindPreciseGlidePlane(home, burg, lDir, glidePlane, node1, node2);
    return;
#endif

    srand(time(NULL));

    burgID = BurgID_FCC(bd, burg);
    GetArmGrain(home, node1, node2, inGrain, onGB);
    burgType = BurgType_FCC(burg);

    if(Normal(lDir) < 1.0E-4){
        interType = NodesIntersection(home, node1, node2, vec1, vec2);   
        switch(interType){
            case LINE_INTERSECTION:
                VECTOR_COPY(lDir, vec1);
                break;
            case PLANE_INTERSECTION:
                VECTOR_COPY(glidePlane, vec1);
                return;
                break;
            default:
                VECTOR_COPY(glidePlane, burg);
                NormalizeVec(glidePlane);
                return;
                break;
        }
    }

    if(Normal(lDir) < 1.0E-4){
        if(burgID<bd->numBurgVectors){
            printf("Warning: define a palne with zero length, arm (%d,%d)\n", inGrain, onGB);
            times = 0;
            while(normalID == -1 && times < 100){
                times++;
                normalID = bd->burgOnPlaneID[burgID][rand()%bd->numPlanesPerBurg[burgID]];
                if(burgType == BASIC_BURG){
                    if((bd->planeType[normalID] & GLIDE_PLANE) == 0){
                        normalID = -1;
                        continue;
                    }
                }
                if(bd->planeInGrain[normalID] != inGrain){
                    normalID = -1;
                }else{
                    VECTOR_COPY(glidePlane, bd->planeList[normalID]);
                    break;
                }

                if(times >= 100){
                    VECTOR_COPY(glidePlane, burg);
                    NormalizeVec(glidePlane);
                    printf("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
                    return;
                }
            }
            return;
        }else{
#if 1
            PrintNode(home, node1);
            PrintNeighbor(home, node1);
            PrintNode(home, node2);
            PrintNeighbor(home, node2);
#endif
            VECTOR_COPY(glidePlane, burg);
            NormalizeVec(glidePlane);
            printf("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
            return;
        }
    }

    NormalizeVec(lDir);
    VECTOR_COPY(normBurg, burg);
    NormalizeVec(normBurg);

    cross(normBurg, lDir, crossProdt);
    if(Normal(crossProdt) < 1.0E-4){
        if(burgID >= bd->numBurgVectors){
            FormatVector(burg, "burg");
            FormatVector(lDir, "lDir");
            VECTOR_COPY(glidePlane, burg);
            NormalizeVec(glidePlane);
            printf("Warning: can not find a glide plane without burgID\n");
            return;
        }else{
            times = 0;
            while(normalID == -1){
                times++;
                normalID = bd->burgOnPlaneID[burgID][rand()%bd->numPlanesPerBurg[burgID]];
                if(burgType == BASIC_BURG){
                    if((bd->planeType[normalID] & GLIDE_PLANE) == 0){
                        normalID = -1;
                        continue;
                    }
                }
                if(bd->planeInGrain[normalID] != inGrain){
                    normalID = -1;
                }else{
                    VECTOR_COPY(glidePlane, bd->planeList[normalID]);
                    break;
                }
                if(times >= 100){
                    VECTOR_COPY(glidePlane, burg);
                    NormalizeVec(glidePlane);
#if 0
                    FormatVector(burg, "burg");
                    PrintNode(home, node1);
                    Fatal("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
#endif
                    return;
                }
            }
        }
    }else{
        NormalizeVec(crossProdt);
        normalID = PlaneID_FCC(bd, crossProdt);
        if(normalID < bd->numPlanes){
            VECTOR_COPY(glidePlane, bd->planeList[normalID]);
        }else{
            VECTOR_COPY(glidePlane, crossProdt);
        }
    }
    
    return;
}


void StackingFaultForce(Home_t *home, real8 dx, real8 dy, real8 dz, real8 sfEnergy, 
                        real8 signSF, real8 normal[3], real8 f1[3], real8 f2[3])
{
    real8       vec[3], length, dir[3];
    VECTOR_ZERO(f1);
    VECTOR_ZERO(f2);
    if(fabs(signSF) < 1.0E-4)return;

    vec[0] = dx;
    vec[1] = dy;
    vec[2] = dz;

    length = Normal(vec);
    if(length < 1.0E-20)return;
    NormalizeVec(vec);

    cross(normal, vec, dir);
    if(Normal(dir) < 1.0E-4)return;
    NormalizeVec(dir);

    f1[0] = length*sfEnergy*0.5*dir[0]*signSF;
    f1[1] = length*sfEnergy*0.5*dir[1]*signSF;
    f1[2] = length*sfEnergy*0.5*dir[2]*signSF;

    VECTOR_COPY(f2, f1);
    return;
}

       
void Coor2EleNum(Param_t *param, double x, double y, real8 *nx, real8 *ny){

	double	eleLx,eleLy,eleLz;

	eleLx = param->fMeshSizeX;
	eleLy = param->fMeshSizeY;

	*nx = floor((x - param->xBoundMin)/eleLx);
	*ny = floor((y - param->yBoundMin)/eleLy);

    return;
}

void RandomInternalForce(Home_t *home, real8 x1, real8 y1, real8 z1, real8 x2, real8 y2, real8 z2,
                         real8 dx, real8 dy, real8 dz, real8 bx, real8 by, real8 bz,
                         real8 f1[3], real8 f2[3]) 
{
    int     i, a, b, c, d;
    real8   imin, imax, jmin, jmax, numX1, numY1;
    real8   vec[3], x0, y0, z0, A, sigR;
    real8   x, y, z, length;
    real8   E = 2.718281828459045;
    real8   gausFac1 = 0.211324865405187; // (1-sqrt(3)/3)/2
    real8   gausFac2 = 0.7886751345948128; // (1+sqrt(3)/3)/2
    real8   l1, l2, nx, ny, f, totForce, norX, norY;

    real8   sigma[3][3], burg[3], strb[3], xi[3];
    real8   *globalArray;
    int     armID;

    Param_t *param;
    param = home->param;

    globalArray = home->randomForceTable;

    VECTOR_ZERO(sigma[0]);
    VECTOR_ZERO(sigma[1]);
    VECTOR_ZERO(sigma[2]);
    
    burg[0] = bx;
    burg[1] = by;
    burg[2] = bz;

    xi[0] = dx;
    xi[1] = dy;
    xi[2] = dz;
    length = Normal(vec);
    NormalizeVec(xi);

    VECTOR_ZERO(f1);
    VECTOR_ZERO(f2);

    l1 = param->fMeshSizeX;
    l2 = param->fMeshSizeY;

    nx = param->fNumX;
    ny = param->fNumY;

    totForce = 0.0;
    for(i=0; i<2; i++){
        if(i==0){
            x = x1 + gausFac1 * dx; 
            y = y1 + gausFac1 * dy; 
            z = z1 + gausFac1 * dz; 
        }else{
            x = x1 + gausFac2 * dx; 
            y = y1 + gausFac2 * dy; 
            z = z1 + gausFac2 * dz; 
        }

        FoldBox(param, &x, &y, &z);
        Coor2EleNum(param, x, y, &imin, &jmin);
        if(imin<0 || imin >= nx || jmin<0 || jmin>=ny){
            Fatal("(%f,%f,%f)%f,%f,%f,%f,%f,%f\n", x, y, z, imin, jmin, nx, ny, param->xBoundMin, param->yBoundMin);
        }

        imax = imin + 1.0;
        jmax = jmin + 1.0;
        if(imax > nx) imax -= nx;
        if(jmax > ny) jmax -= ny;
 
        numX1 = nx;
        numY1 = ny;
        norX = (x-(param->xBoundMin+(imin*l1)))/l1;
        norY = (y-(param->yBoundMin+(jmin*l2)))/l2;
 
        if(norX < -1.0E-10 || norX > 1+1.0E-10 ||
           norY < -1.0E-10 || norY > 1+1.0E-10){
            printf("normalized factors: %f %f, %f %f %f %f %f %f\n",
                   norX, norY, x, y, param->xBoundMin, param->yBoundMin, l1, l2);
            Fatal("interpolation point is out of the cube");
        }
 
        a = (int)(jmin*4+imin*4*numY1);
        b = (int)(jmax*4+imin*4*numY1);
        c = (int)(jmin*4+imax*4*numY1);
        d = (int)(jmax*4+imax*4*numY1);
        
        f = BicubicInterpolation(
                          globalArray[0+a],    globalArray[1+a]*l1,
                          globalArray[2+a]*l2, globalArray[3+a]*l1*l2,
                 
                          globalArray[0+b],    globalArray[1+b]*l1,
                          globalArray[2+b]*l2, globalArray[3+b]*l1*l2,
                 
                          globalArray[0+c],    globalArray[1+c]*l1,
                          globalArray[2+c]*l2, globalArray[3+c]*l1*l2,
                 
                          globalArray[0+d],    globalArray[1+d]*l1,
                          globalArray[2+d]*l2, globalArray[3+d]*l1*l2,
                          norX, norY); 
        totForce += f;
    }
#if 0
    for(i=0; i<param->nPotParams; i++){

        A  = home->potParams[i][0];
        x0 = home->potParams[i][1];
        y0 = home->potParams[i][2];
        z0 = home->potParams[i][3];
        
        vec[0] = x0 - x;
        vec[1] = y0 - y;
        vec[2] = z0 - z;
        ZImage(home->param, vec, vec+1, vec+2);
        x0 = x + vec[0];
        y0 = y + vec[1];
        z0 = z + vec[2];

//        sigma[0][2] += ((A*(x - x0))/(pow(E,(pow(x - x0, 2.0) + pow(y - y0,2.0) + pow(z - z0, 2.0))/
//                 (2.0*pow(sigR,2)))*pow(sigR,2)));
        sigma[0][2] = 0.0;
        sigma[1][2] += ((A*(y - y0))/(pow(E,-(pow(x - x0, 2.0) + pow(y - y0,2.0) + pow(z - z0, 2.0))/
                 (2.0*pow(sigR,2)))*pow(sigR,2)));
        sigma[0][1] = 0.0;
//        sigma[0][1] += ((A*(z - z0))/(pow(E,(pow(x - x0, 2.0) + pow(y - y0,2.0) + pow(z - z0, 2.0))/
//                 (2.0*pow(sigR,2)))*pow(sigR,2)));

    }
#endif

    sigma[0][2] = 0.0;
    sigma[1][2] = 0.5*totForce/param->burgMag;
    sigma[0][1] = 0.0;
    if(isnan(sigma[1][2]))sigma[1][2] = 0.0;

    sigma[2][0] = sigma[0][2];
    sigma[2][1] = sigma[1][2];
    sigma[1][0] = sigma[0][1];

    for (i = 0; i < 3; i++) {
        strb[i] = sigma[i][0]*burg[0]+sigma[i][1]*burg[1]+sigma[i][2]*burg[2];
    }

    cross(strb, xi, f1);
    
    f1[0] *= 0.5;
    f1[1] *= 0.5;
    f1[2] *= 0.5;

    VECTOR_COPY(f2, f1);
    return;
}

bool InCuboid(Param_t *param, real8 pos[3], real8 eps)
{
    if(param->xBoundType != Periodic || param->gbStrength>0){
        if(fabs(pos[0]) > param->xBoundMax + eps)return(0);
    }
    
    if(param->yBoundType != Periodic || param->gbStrength>0){
        if(fabs(pos[1]) > param->yBoundMax + eps)return(0);
    }

    if(param->zBoundType != Periodic || param->gbStrength>0){
        if(fabs(pos[2]) > param->zBoundMax + eps)return(0);
    }
    return(1);
}

void LineBoxIntersection(Param_t *param, real8 dir[3], real8 pos[3], real8 interPos[2][3])
{
    int                     i;
    std::vector<Plane_t>    surface;
    Plane_t                 plane;

    if(param->xBoundType != Periodic){
        VECTOR_ZERO(plane.normal);
        VECTOR_ZERO(plane.point);
        
        plane.normal[0] = 1;
        plane.point[0] = param->xBoundMin;
        surface.push_back(plane);

        plane.point[0] = param->xBoundMax;
        surface.push_back(plane);

    }
    
    if(param->yBoundType != Periodic){
        VECTOR_ZERO(plane.normal);
        VECTOR_ZERO(plane.point);
        
        plane.normal[1] = 1;
        plane.point[1] = param->yBoundMin;
        surface.push_back(plane);

        plane.point[1] = param->yBoundMax;
        surface.push_back(plane);

    }

    if(param->zBoundType != Periodic){
        VECTOR_ZERO(plane.normal);
        VECTOR_ZERO(plane.point);
        
        plane.normal[2] = 1;
        plane.point[2] = param->zBoundMin;
        surface.push_back(plane);

        plane.point[2] = param->zBoundMax;
        surface.push_back(plane);
    }


    for(i=0;i<surface.size();i++){
        
    }
    
    return;
}


