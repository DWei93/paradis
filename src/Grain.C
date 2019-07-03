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


#ifdef GEOS
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <geos/unload.h>
#include <geos_c.h>
#endif

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


void InitGrainStructure(Grain_t *grain, int nGBs){
	grain->nGBs = nGBs;
    grain->GB = (GB_t **)malloc(nGBs*sizeof(GB_t *));	
    grain->sign = (real8 *)malloc(nGBs*sizeof(real8));
    return;
}


bool cmp3(SlipSys_t p, SlipSys_t q){
    return(p.schmid > q.schmid);
}

void DrawLine(FILE *fp, real8 *pos1, real8 *pos2, real8 meshSize, std::vector<int> &indexs)
{
    real8 vec[3], length, nums;
    vec[0] = pos2[0]-pos1[0];vec[1] = pos2[1]-pos1[1];vec[2] = pos2[2]-pos1[2];
    if((length = Normal(vec))<1E-10)return;
    NormalizeVec(vec);
    nums = floor(length/meshSize);
    meshSize = length/nums;

    for(real8 i=0; i<nums; i++){
        fprintf(fp, "%g %g %g %g %g %g",pos1[0]+i*meshSize*vec[0], 
                pos1[1]+i*meshSize*vec[1],pos1[2]+i*meshSize*vec[2], meshSize*vec[0], meshSize*vec[1], meshSize*vec[2]);
        for(auto j=0; j<indexs.size(); j++)fprintf(fp," %d",indexs[j]);
        fprintf(fp,"\n");
    }

    return;
}

void DrawPolygon(FILE *fp, std::vector<int> &f_vert,std::vector<double> &v,int startID, char *zoneName, std::vector<int> &vec){

	int n=f_vert[startID++];

    fprintf(fp,"Zone T = \"%s\", N = %d, E = %d, F = FEPOINT, ET = TRIANGLE\n",
            zoneName, n, n-2);

    std::vector<std::vector<int> > record(2);
    record[0].resize(n);    record[1].resize(n);

    for(int i=startID; i<startID+n; i++){
        fprintf(fp, "%g %g %g", v[f_vert[i]*3], v[f_vert[i]*3+1], v[f_vert[i]*3+2]);
        for(unsigned int j=0; j<vec.size(); j++)fprintf(fp," %d", vec[j]); 
        fprintf(fp, "\n");
        record[0][i-startID] = f_vert[i];   record[1][i-startID] = i-startID+1;        
    }

    for(int i=startID+1; i<startID+n-1; i++){
        int j;
        for(j=1; j<n-1; j++)if(record[0][j] == f_vert[i])break;
        fprintf(fp, "1 %d", record[1][j]);
        for(j=2; j<n; j++)if(record[0][j] == f_vert[i+1])break;
        fprintf(fp, " %d\n", record[1][j]);
    }
    return;
}

void CentroidOfPolygon(std::vector<int> &f_vert,std::vector<double> &v,int startID, real8 *centroid){

	int n=f_vert[startID++];

    real8 pos0[3], pos1[3], pos2[3], cpd[3], areas=0, area;
    VECTOR_ZERO(centroid); 
    for(int i=startID+1; i<startID+n-1; i++){
        pos1[0]=v[f_vert[i]*3]-v[f_vert[startID]*3];
        pos1[1]=v[f_vert[i]*3+1]-v[f_vert[startID]*3+1];
        pos1[2]=v[f_vert[i]*3+2]-v[f_vert[startID]*3+2];
        pos2[0]=v[f_vert[i+1]*3]-v[f_vert[startID]*3];
        pos2[1]=v[f_vert[i+1]*3+1]-v[f_vert[startID]*3+1];
        pos2[2]=v[f_vert[i+1]*3+2]-v[f_vert[startID]*3+2];

        pos0[0] = (pos1[0]+pos2[0])/3.0;
        pos0[1] = (pos1[1]+pos2[1])/3.0;
        pos0[2] = (pos1[2]+pos2[2])/3.0;
        cross(pos1, pos2, cpd);area = Normal(cpd);
        areas += area;
        centroid[0] += area*pos0[0]; 
        centroid[1] += area*pos0[1]; 
        centroid[2] += area*pos0[2]; 
    }
    centroid[0] /= areas; centroid[1] /= areas; centroid[2] /= areas;
    centroid[0] += v[f_vert[startID]*3];
    centroid[1] += v[f_vert[startID]*3+1];
    centroid[2] += v[f_vert[startID]*3+2];
    return;
}

int AddAndSearchVertice(std::vector<real8> &vertices, real8 x, real8 y, real8 z)
{
    if(vertices.size() > 0){
        unsigned int i;
        for(i=0; i<vertices.size(); i+=3){
            if(fabs(vertices[i]-x)<1.0E-12 &&
               fabs(vertices[i+1]-y)<1.0E-12 &&
               fabs(vertices[i+2]-z)<1.0E-12)break;
        }
        if(i == vertices.size()){
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
        return (int)i/3;

    }else{
        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);
        return(0);
    }
}

void EdgeListToPolygon(std::vector<Edge_t> &edges, std::vector<int> &polygon){
    polygon.clear();
    if(edges.empty())Fatal("empty edges");
    polygon.push_back(edges[0].vertID1); 
    if(edges.size()==1){
       polygon.push_back(edges[0].vertID2); return;
    }

    int next=edges[0].vertID2;
    do{
        polygon.push_back(next); 
        for(auto i=1; i<edges.size(); i++){
            if(edges[i].vertID1 == next){
                next = edges[i].vertID2; 
                edges[i].vertID1 = -edges[i].vertID1-1; edges[i].vertID2 = -edges[i].vertID2-1;
                break;
            }else if(edges[i].vertID2 == next){
                next = edges[i].vertID1; 
                edges[i].vertID1 = -edges[i].vertID1-1; edges[i].vertID2 = -edges[i].vertID2-1;
                break;
            }
        }
    }while(next != edges[0].vertID1);
    for(auto i=1; i<edges.size(); i++){
        edges[i].vertID1 = -edges[i].vertID1-1; edges[i].vertID2 = -edges[i].vertID2-1;
    }
    return; 
}

int AddAndSearchEdge(std::vector<Edge_t> &edges, int vertID1, int vertID2, int ofFace, int type, int ofGB=-1, int ofIB=-1)
{
    if(vertID1 == vertID2)Fatal("AddAndSearchSegment: verticeID1 == verticeID2");
    if(vertID2 < vertID1){
        int tempID = vertID1; vertID1 = vertID2; vertID2=tempID;
    }

    if(edges.size()>0){
        unsigned int i;
        for(i=0; i<edges.size(); i++){
            if(vertID1 == edges[i].vertID1 && vertID2 == edges[i].vertID2){
                edges[i].type = type; edges[i].ofIB[0]=ofGB; edges[i].ofIB[1]=ofIB;
                if(edges[i].ofFaces.size()>0){
                    unsigned int j; for(j=0; j<edges[i].ofFaces.size(); j++)if(edges[i].ofFaces[j]==ofFace)break;
                    if(j==edges[i].ofFaces.size())edges[i].ofFaces.push_back(ofFace);
                }else{
                    edges[i].ofFaces.push_back(ofFace);
                }
                break;
            }
        }
        if(i==edges.size()){
            Edge_t edge;
            edge.vertID1 = vertID1; edge.vertID2 = vertID2; edge.type=type;
            edge.ofFaces.push_back(ofFace);
            edge.ofIB[0]=ofGB; edge.ofIB[1]=ofIB;
            edges.push_back(edge);
        }
        return (int)i;
    }else{
        Edge_t edge;
        edge.vertID1 = vertID1; edge.vertID2 = vertID2; edge.type=type;
        edge.ofFaces.push_back(ofFace);
        edge.ofIB[0]=ofGB; edge.ofIB[1]=ofIB;
        edges.push_back(edge);
        return 0;
    }
}

void EdgeCopy(Edge_t &a, Edge_t &b){
    a.type = b.type; a.vertID1=b.vertID1; a.vertID2=b.vertID2;
    a.ofIB[0]=b.ofIB[0]; a.ofIB[1]=b.ofIB[1];
    a.ofFaces.resize(b.ofFaces.size());
    for(auto i=0;i<a.ofFaces.size();i++)a.ofFaces[i]=b.ofFaces[i];
    return;
}

void InitGrain(Home_t *home, Param_t *param){

    /** allocate and construct the container
     */
    home->con = new container(param->xBoundMin,param->xBoundMax,param->yBoundMin,
                    param->yBoundMax, param->zBoundMin, param->zBoundMax, 1, 1, 1,
                    false, false, false, 8);
    container *con = home->con;

    /** import the voro++ input file in domain 0 and boardcast it.
     */
    const int MAXLINELENGTH = 1024;
    FILE *fp;
    char *strp, str[MAXLINELENGTH];
    if(home->myDomain == 0){
        if((fp = fopen(home->voroFile, "r+")) == NULL){
            Fatal("Can not open voro++ file %s\n", home->voroFile);
        }

        home->nGrains = 0;
        while(1){
            strp = fgets(str, MAXLINELENGTH, fp);
            if(feof(fp))break;
            if(str == "\n" || str == NULL)continue;
            home->nGrains++;
        }
        fclose(fp);
    }
#if PARALLEL
    MPI_Bcast(&(home->nGrains), sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
#endif

    real8   *dataList = new real8[home->nGrains*(con->ps+1)];
    if(home->myDomain == 0){
        if((fp = fopen(home->voroFile, "r+")) == NULL){
            Fatal("Can not open voro++ file %s\n", home->voroFile);
        }
        
        char *strok; int i=0;
        while(1){
            strp = fgets(str, MAXLINELENGTH, fp);
            if(feof(fp))break;
            if(str == "\n" || str == NULL)continue;

            dataList[i++] = atof(strtok(str, " "));
            for(int j=1; j<(con->ps+1); j++)dataList[i++] = atof(strtok(NULL, " "));
        }
        fclose(fp);
    }
#if PARALLEL
    MPI_Bcast((char*)dataList, home->nGrains*(con->ps+1)*sizeof(real8), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

#if DEBUG
    if(home->myDomain ==1){
        printf("nGrains %d\n", home->nGrains);
        for(int i=0; i<home->nGrains*(con->ps+1); i++)printf("%lg ",dataList[i]);
        printf("\n");
    }
#endif

    /** put all voronoicell seed
     */
    for(int i=0; i<home->nGrains; i++){
        real8 ibNormal[3]={dataList[i*(con->ps+1)+5], dataList[i*(con->ps+1)+6], dataList[i*(con->ps+1)+7]};
        real8 rotMatrix[3][3]={
            {dataList[i*(con->ps+1)+8], dataList[i*(con->ps+1)+9], dataList[i*(con->ps+1)+10]}, 
            {dataList[i*(con->ps+1)+11], dataList[i*(con->ps+1)+12], dataList[i*(con->ps+1)+13]}, 
            {dataList[i*(con->ps+1)+14], dataList[i*(con->ps+1)+15], dataList[i*(con->ps+1)+16]}};
        con->put((int)dataList[i*(con->ps+1)],
        dataList[i*(con->ps+1)+1], dataList[i*(con->ps+1)+2], dataList[i*(con->ps+1)+3],
        dataList[i*(con->ps+1)+4], ibNormal, rotMatrix);
    }

    delete [] dataList;

    /** allocate the loop class home->cl and voronicell_neighbor
     */
    home->cl = new c_loop_all(*con);
    home->cells = (voronoicell_neighbor **)malloc(home->nGrains*sizeof(voronoicell_neighbor *));
    for(int i=0; i<home->nGrains; i++)home->cells[i] = new voronoicell_neighbor;
    voronoicell_neighbor **cells = home->cells;
    c_loop_all &cl = *home->cl;

    /** Initialize home->grainKeys
     */
    home->grainKeys = new Grain_t *[home->nGrains];
    for(int i=0; i<home->nGrains; i++)home->grainKeys[i] = new Grain_t;

    /** open the file polygons.plt
     */
    if(home->myDomain ==0) {
        if ((fp = fopen("polygons.plt", "w")) == (FILE *)NULL) {
         Fatal("tec_plot: Open error %d on %s\n", errno, "polygons.plt");
        }
        fprintf(fp,"variables = x, y, z, type, grain1, grain2\n");
    }

    /** compute voronoicells, cout all faces and draw them out.
     */
    int     nFaces =0;
	std::vector<int> neigh,f_vert, convertList;
    std::vector<double> vertices, normals, gloVertices;
    real8   ibThick, ibNormal[3], testPoint[3], testRsq;
    std::vector<real8>  ibRsq;
    if(cl.start()) for(int i=0; i<home->nGrains; i++, ibRsq.clear()){
        int k =0; std::vector<int> info(3);
        Grain_t *grain = home->grainKeys[i];
        voronoicell_neighbor *cell = cells[i];
         
        /** compute the voronoicell
         */
        con->compute_cell(*cells[i], cl);

        /** Now, what we should construct in Grain_t are:
         *      int         id; 
	     *      real8	    pos[3]; 
         */
        grain->id = cl.pid();
        grain->pos[0] = cl.x(); grain->pos[1] = cl.y(); grain->pos[2] = cl.z();

		// Gather information about the computed Voronoi cell
		cell->neighbors(neigh);
		cell->face_vertices(f_vert);
		cell->vertices(grain->pos[0],grain->pos[1],grain->pos[2],vertices);

		for(auto j=0, k=0; j<neigh.size();j++) {
			// Skip if the neighbor information is smaller than
			// this grain ID, to avoid double counting.
            if(neigh[j] > grain->id || neigh[j] < 0){
                char zoneName[256];
                if(neigh[j] < 0){
                    snprintf(zoneName,sizeof(zoneName),"FS.%d",nFaces);
                    info[0] = 0;
                }else{
                    info[0] = 1;
                     snprintf(zoneName,sizeof(zoneName),"GB.%d",nFaces);
                }
                 info[1] = grain->id; info[2] = neigh[j];
                if(home->myDomain ==0)DrawPolygon(fp, f_vert, vertices, k, zoneName, info);

                nFaces++;
            }
			k+=f_vert[k]+1;
        }
        if(!cl.inc())break;
    }


    /** alloccate home=>faces
     */
    home->nFaces = nFaces;
    home->faces = new Face_t *[home->nFaces];
    std::vector<Edge_t> edges;
    for(auto i=0; i<home->nFaces; i++)home->faces[i] = new Face_t; 
    if(cl.start()) for(int i=0, nFaces =0; i<home->nGrains; i++, ibRsq.clear()){
        int k =0; std::vector<int> info(3);
        real8   dis; int segID;
        Grain_t *grain = home->grainKeys[i];
        voronoicell_neighbor *cell = cells[i];
        
		// Gather information about the computed Voronoi cell
		cell->neighbors(neigh);
		cell->face_vertices(f_vert);
		cell->vertices(grain->pos[0],grain->pos[1],grain->pos[2],vertices);
        cell->normals(normals);

        // Generate convertion list and add local vertices into global vertices.
        convertList.resize(vertices.size()/3);
        for(auto j=0; j<convertList.size(); j++)
            convertList[j] = AddAndSearchVertice(gloVertices, vertices[3*j], vertices[3*j+1], vertices[3*j+2]);
#if 0
        if(home->myDomain == 0){
            printf("Grain %d: ", grain->id); for(auto j=0; j<convertList.size(); j++){
                printf(" %d(%g,%g,%g)", convertList[j],gloVertices[convertList[j]*3],
                      gloVertices[convertList[j]*3+1],gloVertices[convertList[j]*3+2]);
            }
            printf("\n");
        }
#endif

        // Get the thickness and the normal of inner boundary of the grain and normalize the normal
        ibThick = cl.ib_thick();
        ibNormal[0] = cl.ib_nx(); ibNormal[1] = cl.ib_ny(); ibNormal[2] = cl.ib_nz();
        NormalizeVec(ibNormal); VECTOR_COPY(grain->ibN, ibNormal);

        // Get the first initial rotation matrix and normalize it
        cl.rotMatrix(grain->rotMatrixs[0]);
        NormalizeVec(grain->rotMatrixs[0][0]);NormalizeVec(grain->rotMatrixs[0][1]);NormalizeVec(grain->rotMatrixs[0][2]);
	    Matrix33Invert(grain->rotMatrixs[0], grain->rotMatrixsInv[0]);

        if(ibThick>0){
            real8 temR[3][3]; RotateWithAxis(grain->ibN, PI, temR);
            Matrix33Mult33(temR, grain->rotMatrixs[0], grain->rotMatrixs[1]);
	        Matrix33Invert(grain->rotMatrixs[1], grain->rotMatrixsInv[1]);

            testRsq = ibThick;  while(cell->plane_intersects_guess(ibNormal[0], ibNormal[1], ibNormal[2], testRsq)){
                ibRsq.push_back(testRsq);
                testRsq += 2*ibThick;
            }
            testRsq = ibThick;  while(cell->plane_intersects_guess(-ibNormal[0], -ibNormal[1], -ibNormal[2], testRsq)){
                ibRsq.insert(ibRsq.begin(),-testRsq);
                testRsq += 2*ibThick;
            }
            
            /** Assembly inner interfaces of grain
             */
            grain->nIBs = (int)ibRsq.size(); grain->ibPs = new real8 [grain->nIBs][3];
            for(auto j=0; j<ibRsq.size(); j++){
                grain->ibPs[j][0] = grain->pos[0] + 0.5*grain->ibN[0]*ibRsq[j];
                grain->ibPs[j][1] = grain->pos[1] + 0.5*grain->ibN[1]*ibRsq[j];
                grain->ibPs[j][2] = grain->pos[2] + 0.5*grain->ibN[2]*ibRsq[j];
            }
            
            /** allocate face pointers of the grain
             */
            grain->nFaces = (int)neigh.size();
            grain->faces = new Face_t *[grain->nFaces];
            grain->sign = new real8 [grain->nFaces];
        }

		for(auto j=0, k=0; j<neigh.size();j++) {
            /** convert to global andeliminate duplicate segments.
             */
            for(auto l=k+1; l<k+f_vert[k]+1; l++)f_vert[l]=convertList[f_vert[l]];

            if(neigh[j] > grain->id || neigh[j] < 0){
                for(auto l=k+1; l<k+f_vert[k]+1; l++){
                    if(l==k+f_vert[k])segID = AddAndSearchEdge(edges, f_vert[l],f_vert[k+1], nFaces, 0);
                    else segID = AddAndSearchEdge(edges, f_vert[l],f_vert[l+1], nFaces, 0);
                }
            
                grain->faces[j] = home->faces[nFaces]; grain->sign[j]=1.0;
                Face_t *face = home->faces[nFaces]; face->id = nFaces;
                face->n[0] = normals[j*3]; face->n[1] = normals[j*3+1]; face->n[2] = normals[j*3+2];
                face->p[0] = gloVertices[3*f_vert[k+1]];
                face->p[1] = gloVertices[1+3*f_vert[k+1]];
                face->p[2] = gloVertices[2+3*f_vert[k+1]];
                face->nbrGrain[0] = grain->id; face->nbrGrain[1] = neigh[j];
                face->type = (neigh[j]>=0)?GB_FACE:FREE_SURFACE;
                PointPlaneIntersection(grain->pos, face->n, face->p, face->p, &dis);

                face->nEdges = f_vert[k];
                face->edges = new real8 [face->nEdges*6];
                for(auto l=0; l<face->nEdges; l++){
                    int next = (l<face->nEdges-1)?l+1:0; 
                    real8 gDir[3], lDir[3], p[3], crpd[3];
                    lDir[0] = gloVertices[f_vert[k+1+next]*3] - gloVertices[f_vert[k+1+l]*3];
                    lDir[1] = gloVertices[f_vert[k+1+next]*3+1] - gloVertices[f_vert[k+1+l]*3+1];
                    lDir[2] = gloVertices[f_vert[k+1+next]*3+2] - gloVertices[f_vert[k+1+l]*3+2];
                    if(Normal(lDir) < 1E-10)Fatal("the two vertices are too close");
                    NormalizeVec(lDir);
                    p[0] = 0.5*(gloVertices[f_vert[k+1+next]*3] + gloVertices[f_vert[k+1+l]*3]);
                    p[1] = 0.5*(gloVertices[f_vert[k+1+next]*3+1] + gloVertices[f_vert[k+1+l]*3+1]);
                    p[2] = 0.5*(gloVertices[f_vert[k+1+next]*3+2] + gloVertices[f_vert[k+1+l]*3+2]);
                    
                    cross(face->n, lDir, crpd); NormalizeVec(crpd);
                    if(crpd[0]*(face->p[0]-p[0])+crpd[1]*(face->p[1]-p[1])+crpd[2]*(face->p[2]-p[2])>0){
                        crpd[0] = -crpd[0];crpd[1] = -crpd[1];crpd[2] = -crpd[2];
                    }
                    VECTOR_COPY(&(face->edges[l*6]), crpd); VECTOR_COPY(&(face->edges[l*6+3]), p);
#if 0
                    real8 vec1[3], vec2[3];
                    if(LINE_INTERSECTION != LinePlaneIntersection(&(face->edges[l*6]), &(face->edges[l*6+3]),
                       face->n, face->p, vec1, vec2)){
                        Fatal("Wrong!");
                    }else {
                        real8 dis1, dis2, dis3;
                        PointPlaneIntersection(&(face->edges[l*6+3]),  face->n, face->p, vec1, &dis1);
                        PointPlaneIntersection(&(gloVertices[f_vert[k+1+next]*3]),  face->n, face->p, vec1, &dis2);
                        PointPlaneIntersection(&(gloVertices[f_vert[k+1+l]*3]),  face->n, face->p, vec1, &dis2);
                        printf("dis error %e %e %e\n", dis1, dis2, dis3);
                    }
#endif
                }

                nFaces++;
            }
            
			k+=f_vert[k]+1;
        }

        grain->ibsID = new int [grain->nIBs];
        grain->ibsNumEdges = new int [grain->nIBs]; grain->ibsEdges = new real8 *[grain->nIBs];
        for(auto m=0; m<grain->nIBs; m++){
            grain->ibsID[m] = IB_START_ID+GRAIN_ID_FACTOR*grain->id+m;
            std::vector<Edge_t> ibEdges;
            real8 p1[3], p2[3], p[3], t; int v1,v2;
            for(auto j=0, k=0; j<neigh.size();j++) {
                Edge_t edge; edge.vertID1=-1; edge.vertID2=-1;
                edge.type=1; edge.ofIB[0]=grain->id; edge.ofIB[1]=m;
                for(auto l=k+1; l<k+f_vert[k]+1; l++){
                    v1=f_vert[l]; v2=(l<k+f_vert[k]?f_vert[l+1]:f_vert[k+1]);
                    p1[0] = gloVertices[v1*3];p1[1] = gloVertices[v1*3+1];p1[2] = gloVertices[v1*3+2];
                    p2[0] = gloVertices[v2*3];p2[1] = gloVertices[v2*3+1];p2[2] = gloVertices[v2*3+2];

                    if(POINT_INTERSECTION == SegmentPlaneIntersection(p1, p2, grain->ibN, grain->ibPs[m], t)){
                        p[0] = p1[0]+t*(p2[0]-p1[0]); p[1] = p1[1]+t*(p2[1]-p1[1]); p[2] = p1[2]+t*(p2[2]-p1[2]);
                        if(edge.vertID1 ==-1){
                            edge.vertID1 = AddAndSearchVertice(gloVertices,p[0],p[1],p[2]);
                        }else if(edge.vertID2 == -1){
                            edge.vertID2 = AddAndSearchVertice(gloVertices,p[0],p[1],p[2]);
                        }else Fatal("There is third intersection point");
                    }
                }
                if(edge.vertID2 != -1){
                    edges.push_back(edge);
                    ibEdges.push_back(edge);
                }

                k+=f_vert[k]+1;
            }
            std::vector<int> polygon; real8 centroid[3],dis;
            EdgeListToPolygon(ibEdges,polygon); polygon.insert(polygon.begin(), (int)polygon.size());
            CentroidOfPolygon(polygon,gloVertices,0, centroid);
            PointPlaneIntersection(centroid, grain->ibN, grain->ibPs[m], grain->ibPs[m], &dis);
        
            real8 gDir[3], lDir[3], crpd[3];
            grain->ibsNumEdges[m] = (int)ibEdges.size(); grain->ibsEdges[m] = new real8 [grain->ibsNumEdges[m]*6];
            for(auto n=0; n<grain->ibsNumEdges[m]; n++){
                lDir[0] = gloVertices[ibEdges[n].vertID2*3] - gloVertices[ibEdges[n].vertID1*3];
                lDir[1] = gloVertices[ibEdges[n].vertID2*3+1] - gloVertices[ibEdges[n].vertID1*3+1];
                lDir[2] = gloVertices[ibEdges[n].vertID2*3+2] - gloVertices[ibEdges[n].vertID1*3+2];
                if(Normal(lDir) < 1E-10)Fatal("the two vertices are too close");NormalizeVec(lDir);
                p[0] = 0.5*(gloVertices[ibEdges[n].vertID2*3] + gloVertices[ibEdges[n].vertID1*3]);
                p[1] = 0.5*(gloVertices[ibEdges[n].vertID2*3+1] + gloVertices[ibEdges[n].vertID1*3+1]);
                p[2] = 0.5*(gloVertices[ibEdges[n].vertID2*3+2] + gloVertices[ibEdges[n].vertID1*3+2]);

                cross(grain->ibN, lDir, crpd); NormalizeVec(crpd);
                if(crpd[0]*(grain->ibPs[m][0]-p[0])+crpd[1]*(grain->ibPs[m][1]-p[1])+crpd[2]*(grain->ibPs[m][2]-p[2])>0){
                    crpd[0] = -crpd[0];crpd[1] = -crpd[1];crpd[2] = -crpd[2];
                }
                VECTOR_COPY(&(grain->ibsEdges[m][n*6]), crpd); VECTOR_COPY(&(grain->ibsEdges[m][n*6+3]), p);
            }
            
            if(home->myDomain ==0){
                char zoneName[256];std::vector<int> info(3);
                snprintf(zoneName,sizeof(zoneName),"IB.%d.%d",grain->id, m);
        
                info[0] = 2; info[1] = grain->id; info[2] = m;
                DrawPolygon(fp, polygon, gloVertices, 0, zoneName, info);
            }
        }

        if(!cl.inc())break;
    }

    if(home->myDomain ==0)fclose(fp);

    if(cl.start()) for(int i=0, nFaces=0; i<home->nGrains; i++){
		cells[i]->neighbors(neigh);
        Grain_t *grain = home->grainKeys[i];
        Face_t  *face;
        for(auto j=0; j<neigh.size(); j++){
            if(neigh[j] < grain->id && neigh[j] >=0){
                int k =0; for(k=0; k<home->nFaces; k++){
                    face = home->faces[k];
                    if(face->nbrGrain[0] == neigh[j] && face->nbrGrain[1]==grain->id)break;
                }
                if(k==home->nFaces)Fatal("can not find the %dth face of the grain %d in Home: %d", j, grain->id, neigh[j]);
                grain->faces[j] = face; grain->sign[j]=-1.0;
            }
        }
        if(!cl.inc())break;
    }

    /*  copy vertices into home
     */
    home->nVertices = (int)gloVertices.size()/3; home->vertices = new real8 [3*home->nVertices];
    for(auto i=0; i<3*home->nVertices; i++)home->vertices[i]=gloVertices[i];
    home->nEdges = (int)edges.size(); home->edges = new Edge_t [home->nEdges];
    for(auto i=0;i<home->nEdges;i++)EdgeCopy(home->edges[i], edges[i]);
    
    if(home->myDomain == 0){
        if ((fp = fopen("lines.plt", "w")) == (FILE *)NULL) {
             Fatal("tec_plot: Open error %d on %s\n", errno, "lines.plt");
        }
        fprintf(fp,"variables = x, y, z, vx, vy, vz, segID, type, nFaces\n");

        for(auto i=0; i<edges.size(); i++){
            real8 pos1[3], pos2[3]; std::vector<int> indexs(3);
            indexs[0]= (int)i; indexs[2]=(int)edges[i].ofFaces.size(); indexs[1] = edges[i].type;
            pos1[0] = gloVertices[edges[i].vertID1*3];
            pos1[1] = gloVertices[edges[i].vertID1*3+1];
            pos1[2] = gloVertices[edges[i].vertID1*3+2];
            pos2[0] = gloVertices[edges[i].vertID2*3];
            pos2[1] = gloVertices[edges[i].vertID2*3+1];
            pos2[2] = gloVertices[edges[i].vertID2*3+2];
            DrawLine(fp, pos1, pos2, 10, indexs);
        }   
        fclose(fp);

        if ((fp = fopen("centers.plt", "w")) == (FILE *)NULL) {
            Fatal("tec_plot: Open error %d on %s\n", errno, "centers.plt");
        }
        fprintf(fp,"variables = x, y, z, center_type\n");
        for(auto i=0; i<home->nGrains; i++){
            Grain_t *grain = home->grainKeys[i];
            fprintf(fp, "%g %g %g 0\n", grain->pos[0], grain->pos[1], grain->pos[2]);
            for(auto j=0; j<grain->nFaces; j++)
                fprintf(fp,"%g %g %g 1\n",grain->faces[j]->p[0],grain->faces[j]->p[1],grain->faces[j]->p[2]);
            for(auto j=0; j<grain->nIBs; j++)
                fprintf(fp, "%g %g %g 2\n",grain->ibPs[j][0], grain->ibPs[j][1], grain->ibPs[j][2]);
        }
        fclose(fp);
#ifdef DEBUG
        if ((fp = fopen("vertices.plt", "w")) == (FILE *)NULL) {
            Fatal("tec_plot: Open error %d on %s\n", errno, "vertices.plt");
        }
        fprintf(fp,"variables = x, y, z, place[0], place[1], place[2], place[3]\n");
        int place[4]; real8 pos[3];
        for(auto i=0; i<home->nVertices; i++){
            pos[0]=home->vertices[i*3];pos[1]=home->vertices[i*3+1];pos[2]=home->vertices[i*3+2];
            NodeLocation(home, pos, place);
            fprintf(fp, "%g %g %g %d %d %d %d\n",pos[0],pos[1],pos[2],place[0],place[1],place[2],place[3]);
        }
#endif
    }
    
    return;
}

void NodeLocation(Home_t *home, real8 *pos, int *place)
{
    *place = -1; place[1]=0; place[2]=0; place[3]=0;
    
    int i, j, k, l;
    bool breakLoop=0; int type; real8 vec[3], dis, interPos[3];
    Grain_t *grain; Face_t *face;
    for(i=0; i<home->nGrains; i++){
        grain = home->grainKeys[i];
        for(j=0; j<grain->nFaces; j++){
            face=grain->faces[j];
            type = PointPlaneIntersection(pos, face->n, face->p, vec, &dis);
            dis *= grain->sign[j];
            if(type == POINT_INTERSECTION){
                breakLoop=true;
                if(POINT_INTERSECTION==SegmentPolygonIntersection2D(face->p, pos, face->edges, face->nEdges, interPos)){
                    interPos[0] -= pos[0];interPos[1] -= pos[1];interPos[2] -= pos[2];
                    if(Normal(interPos) > TWIN_LAYER_RADIUS)breakLoop=false;
                }
                if(breakLoop){
                    place[0] = face->type==FREE_SURFACE?ON_FREE_SURFACE:ON_GB; place[1]=face->id;
                }
            }
            if(dis > 0 || breakLoop)break;
        }
        if(breakLoop)break;
        if(j==grain->nFaces){
            place[0] = IN_GRAIN; place[1]=grain->id; break;
        }
    }
#ifdef DEBUG
    if(i==home->nGrains)Fatal("NodeLocation: can not find location of (%g,%g,%g)",*pos,pos[1], pos[2]);
#endif
    switch(place[0]){
        case IN_GRAIN:
            grain=home->grainKeys[place[1]];
            for(i=0; i<grain->nIBs; i++){
                if(PointPlaneIntersection(pos, grain->ibN, grain->ibPs[i], vec, &dis)==POINT_INTERSECTION){
                    place[0]=ON_IB; break;
                }
                if(dis<0)break;
            }
            place[2]=i;break;
        case ON_GB:
            face=home->faces[place[1]];
            for(i=0; i<2; i++){
#ifdef DEBUG
                if(face->nbrGrain[i] < 0)Fatal("nbrGrain[%d]=%d of face %d < 0",i, face->nbrGrain[i], face->id);
#endif
                grain=home->grainKeys[face->nbrGrain[i]];
                for(j=0; j<grain->nIBs; j++){
                    if(PointPlaneIntersection(pos, grain->ibN, grain->ibPs[j], vec, &dis)==POINT_INTERSECTION){
                        place[0]=ON_IB; place[1]=grain->id; place[2]=j; place[3]=0; break;
                    }
                    if(dis<0)break;
                }

                if(place[0]==ON_IB)break;
                else *(place+2+i)=j;
            }break;
        case ON_FREE_SURFACE:
            face=home->faces[place[1]]; grain=home->grainKeys[face->nbrGrain[0]];
#ifdef DEBUG
                if(face->nbrGrain[0] < 0)Fatal("nbrGrain[%d]=%d of face %d < 0",0, face->nbrGrain[0], face->id);
#endif
            for(i=0; i<grain->nIBs; i++){
                if(PointPlaneIntersection(pos, grain->ibN, grain->ibPs[i], vec, &dis)==POINT_INTERSECTION){
                    place[0]=ON_IB; place[1]=grain->id; place[2]=i; place[3]=0; break;
                }
                if(dis<0)break;
            }
            if(place[0]==ON_IB)break;
            else {
                place[2]=grain->id; place[3]=i;
            }break;
        default:
            Fatal("can not supose this location type %d", place[0]);break;
    }
    return;
}
void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int *place)
{
    real8 pos[3] = {x, y, z};
    NodeLocation(home, pos, place); 
    return;
}

void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int& inGrain, int& onGB){
    int         i, j, nGrains, type, breakLoop = 0;
    real8       vec[3], dis, pos[3], sign;
	Grain_t		*grain;
	Param_t		*param;

    inGrain = -1;
    onGB = -1;
	param = home->param;

	if(param->polycrystal == 0){
        return;
    }

	nGrains = home->nGrains;
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    FoldBox(param, pos, pos+1, pos+2);

    for(i=0; i<nGrains; i++){
        if(breakLoop)break;
        grain = home->grainKeys[i];
        for(j=0; j<grain->nGBs; j++){
            type = PointPlaneIntersection(pos, grain->GB[j]->n, 
                   grain->GB[j]->p, vec, &dis);
            dis *= grain->sign[j];

            if(type == POINT_INTERSECTION){
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
                    if(j == grain->nGBs - 1){
                        inGrain = grain->grainID;
                        onGB = -1;
                        breakLoop = 1;
                        break;
                    }
                }
            }
        }   /*for(j=0; j<grain->nGBs; j++) */
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
		printf("Task %d: %s, line %d -- Node = NULL",
			   home->myDomain, __FILE__, __LINE__);
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
            printf("Warning in NodeOnGB, unkonwn neighbor\n");
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
	real8	    tempDistance, vec[3];
	real8 	    minDistance = 1e14;
	GB_t 	    *temGB;
	    	
    if(!home->param->polycrystal)return(0.0);
    int     i;
    for(i=0; i<home->nGBs; i++){
        temGB = home->GBKeys[i];
        PointPlaneIntersection(pos, temGB->n, temGB->p, vec, &tempDistance);
        if(fabs(tempDistance) < fabs(minDistance)){
            *GB = temGB;
            minDistance = tempDistance;
        }
    }
	return(minDistance);
}

int GetGBIndexOfGrain(Grain_t *grain, int onGB){
    int i;
    for(i=0; i<grain->nGBs; i++){
         if(grain->GB[i]->GBID == onGB)return(i);
    }
    return(-1);
}

bool GetPublicGB(Grain_t *grain1, Grain_t *grain2, GB_t **GB)
{
    *GB = (GB_t *)NULL;

    int id, i;
    for(i=0; i<grain1->nGBs; i++){
        id = GetGBIndexOfGrain(grain1, grain2->GB[i]->GBID);
        if(id > -1){
            *GB = grain1->GB[id];
            return(1);
        } 
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

    if(onGB < 0 || onGB >= home->nGBs)Fatal("in %s at %d", __FILE__, __LINE__);
	if(node == (Node_t *)NULL){
        Fatal("Task %d: %s, line %d -- Node = NULL",
               home->myDomain, __FILE__, __LINE__);
	}

#if ENABLE_COMPLICATED_REACTION
    eps = param->rann/2.0;
#else
    eps = 1.0;
#endif

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
            GB_t    *gb;
            int     node1SegID;

            if(Connected(node1, node2, &node1SegID)){
                PrintNode(home, node1);
                PrintNode(home, node2);
                Fatal("GetArmGrain: node1 in %d grain, but node2 in %d grain\n",
                      node1->inGrain, node2->inGrain);
            }

            if(GetPublicGB(home->grainKeys[node1->inGrain], 
                        home->grainKeys[node2->inGrain], &gb)){
                inGrain = -1;
                onGB = gb->GBID;
            }else{
                PrintNode(home, node1);
                PrintNode(home, node2);
                Fatal("Warning in GetArmGrain: node1 in %d grain, but node2 in %d grain\n",
                      node1->inGrain, node2->inGrain);
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

    real8       dis, burgMag, burg[3], segCentor[3];
    real8       vec[3], a = 5.8685446009;
    real8       coff, length, nor[3], gDir[3], gbNormal[3];
    GB_t        *gb;
    Param_t     *param;
    int         armID;

    VECTOR_ZERO(f1);
    VECTOR_ZERO(f2);

    if(!home->param->polycrystal)return;
    if(IsOnGB(node1) || IsOnGB(node2))return;

    param = home->param;

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
    if(fabs(DotProduct(burg, param->TBNormal)) < 1.0E-2)return;

    vec[0] = node2->x - node1->x;
    vec[1] = node2->y - node1->y;
    vec[2] = node2->z - node1->z;
    ZImage(param, vec, vec+1, vec+2);

    length = Normal(vec);
    if(length < 1.0E-2)return;

    segCentor[0] = node1->x + 0.5*vec[0];
    segCentor[1] = node1->y + 0.5*vec[1];
    segCentor[2] = node1->z + 0.5*vec[2];
    FoldBox(param, segCentor, segCentor+1, segCentor+2);

	dis =  ClosedGB(home, segCentor, &gb);
    if(fabs(dis) > 100 || fabs(dis) < 1.0E-3)return;

    NormalizeVec(vec);
    cross(nor, vec, gDir);
    if(Normal(gDir) < 1E-3)return;
    NormalizeVec(gDir);

    gbNormal[0] = ((real8)Sign(dis)) * gb->n[0];
    gbNormal[1] = ((real8)Sign(dis)) * gb->n[1];
    gbNormal[2] = ((real8)Sign(dis)) * gb->n[2];

    if(DotProduct(gbNormal, gDir) < 0.0){
        gDir[0] = -gDir[0];
        gDir[1] = -gDir[1];
        gDir[2] = -gDir[2];
    }

    coff = param->GBStrength * a/(fabs(dis)+a);

    f1[0] = 0.5*length*coff*gDir[0];
    f1[1] = 0.5*length*coff*gDir[1];
    f1[2] = 0.5*length*coff*gDir[2];

    VECTOR_COPY(f2, f1);

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
                    break;
                case 2:
                    type = LinePlaneIntersection(glideConstraintsNode2[2], p2, glideConstraintsNode1[0], p1,
                           vec1, vec2);
                    break;
                case 3:
                    type = PointPlaneIntersection(p2, glideConstraintsNode1[0], p1, vec2, &dis); 
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
                           vec1, vec2); 
                    break;
                case 2:
                    type = LineLineIntersection(glideConstraintsNode1[2], p1, glideConstraintsNode2[2], p2,
                           vec1, vec2, &dis); 
                    break;
                case 3:
                    type = PointLineIntersection(p2, glideConstraintsNode1[2], p1, vec2, &dis); 
                    VECTOR_COPY(vec2, p2);
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
                    type = PointPlaneIntersection(p1, glideConstraintsNode2[0], p2,
                           vec2, &dis); 
                    VECTOR_COPY(vec2, p1);
                    break;
                case 2:
                    type = PointLineIntersection(p1, glideConstraintsNode2[2], p2, vec2, &dis); 
                    VECTOR_COPY(vec2, p1);
                    break;
                case 3:
                    type = PointPointIntersection(p2, p1, vec2, &dis); 
                    VECTOR_COPY(vec2, p2);
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
            Fatal("BeforeMergeArm: wrong type %d", __LINE__);
            break;
    }

    if(fabs(*dis) > 500.0 && type != POINT_INTERSECTION){
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

        Fatal("%d: Merge dis %f %d (%f,%f,%f) - %d,%d\n",
                home->myDomain, *dis, type, pos[0],
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
            case POINT_INTERSECTION:
                VECTOR_COPY(glidePlane, burg);
                NormalizeVec(glidePlane);
                return;
                break;
            default:
                Fatal("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
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
                    printf("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
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
    if(param->xBoundType != Periodic){
        if(fabs(pos[0]) > param->xBoundMax + eps)return(0);
    }
    
    if(param->yBoundType != Periodic){
        if(fabs(pos[1]) > param->yBoundMax + eps)return(0);
    }

    if(param->zBoundType != Periodic){
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

#ifdef VORO
void VoroTest()
{
// This function returns a random floating point number between 0 and 1

    const double x_min=-1,x_max=1;
    const double y_min=-1,y_max=1;
    const double z_min=-2,z_max=2;
    const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

// Set up the number of blocks that the container is divided into
    const int n_x=10,n_y=10,n_z=10;

// Set the number of particles that are going to be randomly introduced
    const int particles=2;

// This function returns a random double between 0 and 1

	int i;
	double x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Randomly add particles into the container
#if 0
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}
#else
    con.put(0, 0, y_min+0.25*(y_max-y_min), 0);
    con.put(1, 0, y_min+0.75*(y_max-y_min), 0);
#endif
    c_loop_all vl(con);
	voronoicell v, *vcs;
    vcs = new voronoicell[particles];
    i = 0;
	if(vl.start()) do if(con.compute_cell(v,vl)) {
        vcs[i++] = v;
	} while(vl.inc());

    for(i=0; i<particles; i++){
        voronoicell &vr = vcs[i];
        std::vector<real8>  vec;
	    printf("Total vertices      : %d\n",vr.p);
	    printf("Vertex positions    : ");vr.output_vertices();puts("");
	    printf("Vertex orders       : ");vr.output_vertex_orders();puts("");
        printf("Vertex positions (2): ");
        vr.vertices(vec);
        for(int j=0; j<vr.p; j++){
            printf("(%g,%g,%g) ", vec[3*j], vec[3*j+1], vec[3*j+2]);
        }puts("");

	    printf("Max rad. sq. vertex : %g\n",0.25*vr.max_radius_squared());

	    // Output edge-based statistics
	    printf("Total edges         : %d\n",vr.number_of_edges());
	    printf("Total edge distance : %g\n",vr.total_edge_distance());
	    printf("Face perimeters     : ");vr.output_face_perimeters();puts("");

	    // Output face-based statistics
	    printf("Total faces         : %d\n",vr.number_of_faces());
	    printf("Surface area        : %g\n",vr.surface_area());
	    printf("Face freq. table    : ");vr.output_face_freq_table();puts("");
	    printf("Face orders         : ");vr.output_face_orders();puts("");
	    printf("Face areas          : ");vr.output_face_areas();puts("");
	    printf("Face normals        : ");vr.output_normals();puts("");
	    printf("Face vertices       : ");vr.output_face_vertices();puts("\n");
    }

	// Sum up the volumes, and check that this matches the container volume
	double vvol=con.sum_cell_volumes();
	printf("Container volume : %g\n"
	       "Voronoi volume   : %g\n"
	       "Difference       : %g\n",cvol,vvol,vvol-cvol);

	// Output the particle positions in gnuplot format
	con.draw_particles("random_points_p.gnu");

	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("random_points_v.gnu");
    delete [] vcs;
    return;
} 
#endif

#ifdef GEOS
static void notice(const char *fmt, ...)
{
    std::fprintf( stdout, "NOTICE: ");

    va_list ap;
    va_start(ap, fmt);
    std::vfprintf(stdout, fmt, ap);
    va_end(ap);

    std::fprintf(stdout, "\n");
}

void GeosTest()
{

    GEOSGeometry* geom1_;
    GEOSGeometry* geom2_;
    GEOSWKTWriter* w_;
    GEOSGeom_destroy(geom1_);
    GEOSGeom_destroy(geom2_);
    GEOSWKTWriter_destroy(w_);
    geom1_ = 0;
    geom2_ = 0;
    initGEOS(notice, notice);
    w_ = GEOSWKTWriter_create();
    GEOSWKTWriter_setTrim(w_, 1);
    
    geom1_ = GEOSGeomFromWKT("LINESTRING(0 0, 10 10.2, 10 0, 0 10)");
    geom2_ = GEOSNode(geom1_);
    
    GEOSNormalize(geom2_);
    char* wkt_c = GEOSWKTWriter_write(w_, geom2_);
    std::string out(wkt_c); 
    free(wkt_c);
    
    printf("out: %s\n", out.c_str());
    GEOSGeom_destroy(geom1_);
    GEOSGeom_destroy(geom2_);
    GEOSWKTWriter_destroy(w_);
    geom1_ = 0;
    geom2_ = 0;
    finishGEOS();
    
    return;
}
#endif
