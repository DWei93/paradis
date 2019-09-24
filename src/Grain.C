/****************************************************************************
* 
*	Module:		    Grain.C
*
*	Description:	This module contains basic subroutines needed in 
*                  	voronoi polycrystal model.
*
*	12/06/2019:	    created by DeAn Wei
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
#include "Math.h"

#ifdef GEOS
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <geos/unload.h>
#include <geos_c.h>
#endif

/**************************************************************InitGrain*************************************************************/
//{
bool CmpKeyNodeT(EdgeJoint_t &p, EdgeJoint_t &q){ return p.dis<q.dis;}
void DrawLine(FILE *fp, real8 *pos1, real8 *pos2, real8 meshSize, std::vector<int> &indexs)
{
    real8 vec[3], length, nums;
    CoordMinus(pos2, pos1, vec);
    if((length = normal(vec))<1E-20)return;
    normalize_vec(vec);
    nums = floor(length/meshSize);
    meshSize = length/nums;

    for(real8 i=0; i<nums; i++){
        fprintf(fp, "%g %g %g %g %g %g",pos1[0]+i*meshSize*vec[0], 
                pos1[1]+i*meshSize*vec[1],pos1[2]+i*meshSize*vec[2], meshSize*vec[0], meshSize*vec[1], meshSize*vec[2]);
        for(auto j=0; j<int(indexs.size()); j++)fprintf(fp," %d",indexs[j]);
        fprintf(fp,"\n");
    }

    return;
}

/** Draw polygon with tecplot formart
 *  \param[in] fp file pointer for output
 *  \param[in] f_vert vertices index list
 *  \param[in] v coordinates of vertices
 *  \param[in] startID the start index in f_vert for output polygon
 *  \param[in] zoneName the zone name of tecplot
 *  \param[in] vec the variales of the polygon */
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

int AddAndSearchVertice(std::vector<real8> &vertices, std::vector<std::vector<int> > &verticesOnGrain, real8 x, real8 y, real8 z, int onGrain)
{
#ifdef DEBUG
    if(verticesOnGrain.size() != vertices.size()/3)Fatal("%s: wrong",__func__);
#endif
    if(vertices.size() > 0){
        unsigned int i;
        for(i=0; i<vertices.size(); i+=3){
            if(fabs(vertices[i]-x)<LENGTH_ERROR &&
               fabs(vertices[i+1]-y)<LENGTH_ERROR &&
               fabs(vertices[i+2]-z)<LENGTH_ERROR)break;
        }
        if(i == vertices.size()){
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
            std::vector<int> verticeOnGrain(1);
            verticeOnGrain[0]=onGrain;
            verticesOnGrain.push_back(verticeOnGrain);
            return int(i/3);
        }else{
            i/=3;
#ifdef DEBUG
            if(verticesOnGrain[i].size()==0)Fatal("%s: the vertice was not defined location initially", __func__);
#endif
            unsigned int j; for(j=0; j<verticesOnGrain[i].size(); j++){
                if(onGrain == verticesOnGrain[i][j])break;
            }
            if(j==verticesOnGrain[i].size())verticesOnGrain[i].push_back(onGrain);
            return int(i);
        }

    }else{
        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);
        verticesOnGrain.resize(1);
        verticesOnGrain[0].push_back(onGrain);
        return(0);
    }
}

/** Convert edge vector list to polygon vector list */
void EdgeListToPolygon(std::vector<EdgeVec_t> &edges, std::vector<int> &polygon){

#ifdef DEBUG
    if(edges.empty())Fatal("%s: empty edges",__func__);
#endif

    polygon.clear(); polygon.push_back(edges[0].vertID1); 
    if(edges.size()==1){
       polygon.push_back(edges[0].vertID2); return;
    }

    int next=edges[0].vertID2;
    do{
        polygon.push_back(next); 
        for(auto i=1; i<int(edges.size()); i++){
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
    for(auto i=1; i<int(edges.size()); i++){
        edges[i].vertID1 = -edges[i].vertID1-1; edges[i].vertID2 = -edges[i].vertID2-1;
    }

    for(auto i=0; i<int(edges.size())-1; i++){
        if((edges[i].vertID1==polygon[i] && edges[i].vertID2==polygon[i+1]) ||
           (edges[i].vertID2==polygon[i] && edges[i].vertID1==polygon[i+1]) )continue;
        unsigned int j; for(j=i+1; j<edges.size(); j++){
            if((edges[j].vertID1==polygon[i] && edges[j].vertID2==polygon[i+1]) ||
               (edges[j].vertID2==polygon[i] && edges[j].vertID1==polygon[i+1]) ){
                std::swap(edges[i], edges[j]);break;}
        }
        if(j==edges.size())Fatal("%s: can not find the edge",__func__);
    }

    polygon.insert(polygon.begin(), (int)polygon.size());
    return; 
}

int AddAndSearchEdge(std::vector<EdgeVec_t> &edges, int vertID1, int vertID2, int ofFace)
{
    if(vertID1 == vertID2)Fatal("AddAndSearchSegment: verticeID1 == verticeID2");
    if(vertID2 < vertID1){int tempID = vertID1; vertID1 = vertID2; vertID2=tempID;}

    if(edges.size()>0){
        unsigned int i;
        for(i=0; i<edges.size(); i++){
            if(vertID1 == edges[i].vertID1 && vertID2 == edges[i].vertID2){
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
            EdgeVec_t edge;
            edge.vertID1 = vertID1; edge.vertID2 = vertID2;
            edge.ofFaces.push_back(ofFace);
            edges.push_back(edge);
        }
        return (int)i;
    }else{
        EdgeVec_t edge;
        edge.vertID1 = vertID1; edge.vertID2 = vertID2;
        edge.ofFaces.push_back(ofFace);
        edges.push_back(edge);
        return 0;
    }
}

void EdgeCopy(Edge_t *a, EdgeVec_t &b){
    a->vertID1 = b.vertID1; a->vertID2 = b.vertID2; 
    a->nFaces = int(b.ofFaces.size());
    if(a->nFaces>0)a->ofFaces = (int *)malloc(a->nFaces*sizeof(int));
    for(auto i=0; i<a->nFaces; i++)a->ofFaces[i] = b.ofFaces[i];

    a->nIBs = (b.ibs.empty()?0:int(b.ibs.size()));
    if(a->nIBs>0){
        a->ibInGrain = (int *)malloc(a->nIBs*sizeof(int));
        a->ibID = (int *)malloc(a->nIBs*sizeof(int));
        a->dis = (real8 *)malloc(a->nIBs*sizeof(real8));
        a->ibSign = (real8 *)malloc(a->nIBs*sizeof(real8));
        for(auto i=0; i<a->nIBs; i++){
            a->ibInGrain[i] = b.ibs[i].inGrain;
            a->ibID[i] = b.ibs[i].ID;
            a->dis[i] = b.ibs[i].dis;
            a->ibSign[i] = b.ibs[i].sign; 
        }
    }
    return;
}

inline void InitBurgInfo(BurgInfo_t *bInf)
{bInf->numPlanes = 0; bInf->numBurgVectors=0; bInf->numKeyLines=0; return;}

/** Rotation matrix about axis in current coordinate system
 *  \param[in] a rotation axis
 *  \param[in] theta rotation angle
 *  \param[out] R rotation matrix     */
void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]){
    
    normalize_vec(a);
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

#ifdef DEBUG
inline void NodeLocationDebug(Home_t *home, real8 *pos, int *place, real8 error=FACE_ERROR)
{
    NodeLocation(home, pos, place, error);
    if(!AtLocation(home, place, pos, error)){
        FormatPlace(place);
        Fatal("%s: position (%g,%g,%g) at wrong place [%d,%d,%d,%d]",
             __func__, pos[0],pos[1],pos[2],place[0],place[1],place[2],place[3]);
    }
}
#endif

#ifdef DEBUG
    #define NodeLocation(a,b,c,d) NodeLocationDebug(a,b,c,d)
#endif

void InitGrain(Home_t *home, Param_t *param)
{

    /** allocate and construct the non-periodical container
     */
    home->con = new container(param->xBoundMin,param->xBoundMax,param->yBoundMin,
                    param->yBoundMax, param->zBoundMin, param->zBoundMax, 1, 1, 1,
                    false, false, false, 8);
    container *con = home->con;

    /* Only support FCC material at present*/
    param->materialType = MAT_TYPE_FCC;

    /** import the voro++ input file in domain 0 and boardcast it.
     */
    const int MAXLINELENGTH = 1024;
    FILE *fp=(FILE *)NULL;
    char str[MAXLINELENGTH];
    if(home->myDomain == 0){
        if((fp = fopen(home->voroFile, "r+")) == NULL){
            Fatal("Can not open voro++ file %s\n", home->voroFile);
        }

        home->nGrains = 0;
        while(1){
            fgets(str, MAXLINELENGTH, fp);
            if(feof(fp))break; 
            if(str == (char *)NULL){continue;} if(*str == '\n'){continue;}
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
        
        int i=0;
        while(1){
            fgets(str, MAXLINELENGTH, fp);
            if(feof(fp))break;
            if(str == (char *)NULL){continue;} if(*str == '\n'){continue;}

            dataList[i++] = atof(strtok(str, " "));
            for(int j=1; j<(con->ps+1); j++)dataList[i++] = atof(strtok(NULL, " "));
        }
        fclose(fp);
    }
#if PARALLEL
    MPI_Bcast((char*)dataList, home->nGrains*(con->ps+1)*sizeof(real8), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

#ifdef DEBUG
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
    param->polycrystal=1;
    home->cells = (voronoicell_neighbor **)malloc(home->nGrains*sizeof(voronoicell_neighbor *));
    for(int i=0; i<home->nGrains; i++)home->cells[i] = new voronoicell_neighbor;
    voronoicell_neighbor **cells = home->cells;
    c_loop_all &cl = *home->cl;

    /** Initialize home->grainKeys
     */
    home->grainKeys = (Grain_t **)malloc(home->nGrains*sizeof(Grain_t*));
    for(int i=0; i<home->nGrains; i++)home->grainKeys[i] = (Grain_t *)malloc(sizeof(Grain_t));

    /** open the file polygons.plt
     */
    if(home->myDomain ==0) {
        if ((fp = fopen("polygons.plt", "w")) == (FILE *)NULL) {
            Fatal("%s: Open error %d on %s\n", __func__, errno, "polygons.plt");
        }
        fprintf(fp,"variables = x, y, z, type, grain, neiOrIBID\n");
    }

    /** compute voronoicells, cout all faces and draw them out.
     */
    int     nFaces =0;
	std::vector<int> neigh,f_vert, convertList;
    std::vector<double> vertices, normals, gloVertices;
    std::vector<std::vector<int> > verticesOnGrain;
    real8   ibThick, ibNormal[3], testRsq;
    std::vector<real8>  ibRsq;
    param->volume = 0;
    if(cl.start()) for(int i=0; i<home->nGrains; i++, ibRsq.clear()){
        std::vector<int> info(3);
        Grain_t *grain = home->grainKeys[i];
        voronoicell_neighbor *cell = cells[i];
         
        /** compute the voronoicell
         */
        con->compute_cell(*cells[i], cl);

        /** Now, what we should construct in Grain_t are:
         *      int         id; 
	     *      real8	    pos[3], volume, rMax; 
         */
        grain->id = cl.pid(); grain->volume=cell->volume(); param->volume += grain->volume;
        grain->pos[0] = cl.x(); grain->pos[1] = cl.y(); grain->pos[2] = cl.z();

		// Gather information about the computed Voronoi cell
		cell->neighbors(neigh);
		cell->face_vertices(f_vert);
		cell->vertices(grain->pos[0],grain->pos[1],grain->pos[2],vertices);

        real8 pos[3]; grain->rMax=0;
        for(auto j=0; j<int(vertices.size())/3; j++){
            Coordinate(vertices,j,pos);
            grain->rMax=MAX(grain->rMax,PointPointDistance(grain->pos, pos));
        }

		for(auto j=0, k=0; j<int(neigh.size());j++) {
			// Skip if the neighbor information is smaller than
			// this grain ID, to avoid double counting.
            if(neigh[j] > grain->id || neigh[j] < 0){
                char zoneName[256];
                if(neigh[j] < 0){
                    snprintf(zoneName,sizeof(zoneName),"FS-%d",nFaces);
                    info[0] = 0;
                }else{
                    info[0] = 1;
                     snprintf(zoneName,sizeof(zoneName),"GB-%d",nFaces);
                }
                 info[1] = grain->id; info[2] = neigh[j];
                if(home->myDomain ==0)DrawPolygon(fp, f_vert, vertices, k, zoneName, info);

                nFaces++;
            }
			k+=f_vert[k]+1;
        }
        if(!cl.inc())break;
    }

    /** alloccate home->faces
     */
    home->nFaces = nFaces;
    home->faces = (Face_t **)malloc(home->nFaces*sizeof(Face_t*));
    std::vector<EdgeVec_t> edges;
    for(auto i=0; i<home->nFaces; i++)home->faces[i] = (Face_t *)malloc(sizeof(Face_t)); 
    if(cl.start()) for(int i=0, nFaces =0; i<home->nGrains; i++, ibRsq.clear()){
        std::vector<int> info(3);
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
        for(auto j=0; j<int(convertList.size()); j++)
            convertList[j] = AddAndSearchVertice(gloVertices, verticesOnGrain, vertices[3*j], vertices[3*j+1], vertices[3*j+2], i);
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
        normalize_vec(ibNormal); VECTOR_COPY(grain->ibN, ibNormal);

        // Get the first initial rotation matrix and normalize it
        cl.rotMatrix(grain->rotMatrixs[0]);
        normalize_vec(grain->rotMatrixs[0][0]);normalize_vec(grain->rotMatrixs[0][1]);normalize_vec(grain->rotMatrixs[0][2]);
	    Matrix33Invert(grain->rotMatrixs[0], grain->rotMatrixsInv[0]);

        // Generate Burgers Info
        GenerateBurgInfo(param->materialType, grain->rotMatrixs[0], &(grain->bInfs[0]));

        // generate positions and normals of inner boundaries
        if(ibThick>0){
            real8 temR[3][3]; RotateWithAxis(grain->ibN, PI, temR);
            Matrix33Mult33(temR, grain->rotMatrixs[0], grain->rotMatrixs[1]);
	        Matrix33Invert(grain->rotMatrixs[1], grain->rotMatrixsInv[1]);
            GenerateBurgInfo(param->materialType, grain->rotMatrixs[1], &(grain->bInfs[1]));

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
            grain->nIBs = (int)ibRsq.size();
            grain->ibPs = (real8 (*)[3])malloc(grain->nIBs*3*sizeof(real8));
            grain->ibAreas = (real8 *)malloc(grain->nIBs*sizeof(real8));
            for(auto j=0; j<int(ibRsq.size()); j++){
                grain->ibPs[j][0] = grain->pos[0] + 0.5*grain->ibN[0]*ibRsq[j];
                grain->ibPs[j][1] = grain->pos[1] + 0.5*grain->ibN[1]*ibRsq[j];
                grain->ibPs[j][2] = grain->pos[2] + 0.5*grain->ibN[2]*ibRsq[j];
            }
            
            /** allocate face pointers of the grain
             */
            grain->nFaces = (int)neigh.size();
            grain->faces = (Face_t **)malloc(grain->nFaces*sizeof(Face_t *));
            grain->sign = (real8*)malloc(grain->nFaces*sizeof(real8));
        }else grain->nIBs=0;

        /** Assemble faces of grain i
         */
		for(auto j=0, k=0; j<int(neigh.size());j++) {
            /** convert to global and eliminate duplicate segments.
             */
            for(auto l=k+1; l<k+int(f_vert[k])+1; l++)f_vert[l]=convertList[f_vert[l]];

            /** currentlt, all faces connnected neighbor grain lower than this grain have been constructed.
             *  so we just construt face with higher grain ID and free surface here */ 
            if(neigh[j] > grain->id || neigh[j] < 0){
                grain->faces[j] = home->faces[nFaces]; grain->sign[j]=1.0;
                Face_t *face = home->faces[nFaces]; face->id = nFaces;
                face->nEdges = f_vert[k];
                face->edgesID = (int *)malloc(face->nEdges*sizeof(int));

                for(auto l=k+1; l<k+f_vert[k]+1; l++){
                    if(l==k+f_vert[k])segID = AddAndSearchEdge(edges, f_vert[l],f_vert[k+1], nFaces);
                    else segID = AddAndSearchEdge(edges, f_vert[l],f_vert[l+1], nFaces);
                    face->edgesID[l-k-1] = segID;
                }
            
                face->n[0] = normals[j*3]; face->n[1] = normals[j*3+1]; face->n[2] = normals[j*3+2];
                face->p[0] = gloVertices[3*f_vert[k+1]];
                face->p[1] = gloVertices[1+3*f_vert[k+1]];
                face->p[2] = gloVertices[2+3*f_vert[k+1]];
                face->nbrGrain[0] = grain->id; face->nbrGrain[1] = neigh[j];
                face->type = (neigh[j]>=0)?GB_FACE:FREE_SURFACE;
                PointPlaneDistance(grain->pos, face->n, face->p, dis, face->p);

                face->edges = (real8*)malloc(face->nEdges*6*sizeof(real8));
                for(auto l=0; l<face->nEdges; l++){
                    int next = (l<face->nEdges-1)?l+1:0; 
                    real8 lDir[3], p[3], crpd[3];
#ifdef DEBUG
                    if(f_vert[k+1+next]==f_vert[k+1+l])Fatal("something wrong");
#endif
                    lDir[0] = gloVertices[f_vert[k+1+next]*3] - gloVertices[f_vert[k+1+l]*3];
                    lDir[1] = gloVertices[f_vert[k+1+next]*3+1] - gloVertices[f_vert[k+1+l]*3+1];
                    lDir[2] = gloVertices[f_vert[k+1+next]*3+2] - gloVertices[f_vert[k+1+l]*3+2];
                    if(normal(lDir) < 1E-10)Fatal("the two vertices are too close %e",normal(lDir));
                    normalize_vec(lDir);
                    p[0] = 0.5*(gloVertices[f_vert[k+1+next]*3] + gloVertices[f_vert[k+1+l]*3]);
                    p[1] = 0.5*(gloVertices[f_vert[k+1+next]*3+1] + gloVertices[f_vert[k+1+l]*3+1]);
                    p[2] = 0.5*(gloVertices[f_vert[k+1+next]*3+2] + gloVertices[f_vert[k+1+l]*3+2]);
                    
                    cross_product(face->n, lDir, crpd); normalize_vec(crpd);
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
        } // for(auto j=0, k=0; j<neigh.size();j++) 

        if(!cl.inc())break;
    }

    if(cl.start()) for(int i=0; i<home->nGrains; i++){
        Face_t  *face=(Face_t *)NULL;
        std::vector<real8> grainVertices;
        voronoicell_neighbor *cell = cells[i];
        Grain_t *grain = home->grainKeys[i];

		// Gather information about the computed Voronoi cell
		cell->neighbors(neigh);
		cell->face_vertices(f_vert);
        std::vector<int> glo_f_vert(f_vert.size());
        for(auto j=0; j<int(f_vert.size()); j++)glo_f_vert[j] = f_vert[j];
		cell->vertices(grain->pos[0],grain->pos[1],grain->pos[2],vertices);

        // Generate convertion list and add local vertices into global vertices.
        convertList.resize(vertices.size()/3);
        for(auto j=0; j<int(convertList.size()); j++)
            convertList[j] = AddAndSearchVertice(gloVertices, verticesOnGrain, vertices[3*j], vertices[3*j+1], vertices[3*j+2], i);

        /** grain claims the face which is back to it
         */
        for(auto j=0, k=0; j<int(neigh.size()); j++){
            /** convert to global and eliminate duplicate segments.
             */
            for(auto l=k+1; l<k+int(f_vert[k])+1; l++)glo_f_vert[l]=convertList[f_vert[l]];
            if(neigh[j] < grain->id && neigh[j] >=0){
                int l; for(l=0; l<home->nFaces; l++){
                    face = home->faces[l];
                    if(face->nbrGrain[0] == neigh[j] && face->nbrGrain[1]==grain->id)break;
                }
                if(l==home->nFaces)Fatal("%dth face of the grain %d can not be claimed by faces in home, negihbor grain is %d", j, grain->id, neigh[j]);
                grain->faces[j] = face; grain->sign[j]=-1.0;
            }
			k+=f_vert[k]+1;
        }

        /**  Now, we already have whole information of face and basic information of inner boundary.
         *   It is time to assemble edge information of inner boudaries */
        if(grain->nIBs!=0){
            grain->ibsNumEdges = (int*)malloc(grain->nIBs*sizeof(int));
            grain->ibsEdges = (real8**)malloc(grain->nIBs*sizeof(real8*));
            grain->ibsEdgesFaceID = (int **)malloc(grain->nIBs*sizeof(int*));
            grain->ibPolygons = (int **)malloc(grain->nIBs*sizeof(int*));
        }
        for(auto m=0; m<grain->nIBs; m++){
            real8 p1[3], p2[3], p[3]; int v1,v2;
    
            // assemble edges of inner boundary m
            std::vector<EdgeVec_t> ibEdges;
            for(auto j=0, k=0; j<int(neigh.size());j++){
                face = grain->faces[j];
                EdgeVec_t edge; edge.vertID1=-1; edge.vertID2=-1; edge.ofFaces.push_back(face->id);
                // loop through all edges of the face, find intersections between face and inner boundar m
                for(auto l=k+1; l<k+int(glo_f_vert[k])+1; l++){
                    v1=glo_f_vert[l]; v2=(l<k+glo_f_vert[k]?glo_f_vert[l+1]:glo_f_vert[k+1]);
                    Coordinate(gloVertices, v1, p1);  Coordinate(gloVertices, v2, p2);
                    
                    // if find the intersection, add intersection to IB edge as vertice and face edge as node
                    if(POINT_INTERSECTION == SegmentPlaneIntersection(p1, p2, grain->ibN, grain->ibPs[m], p)){
                        if(edge.vertID1 == -1){
                            edge.vertID1 = AddAndSearchVertice(grainVertices,p[0],p[1],p[2]);
                        }else if(edge.vertID2 == -1){
                            edge.vertID2 = AddAndSearchVertice(grainVertices,p[0],p[1],p[2]);
                        }else Fatal("There is third intersection point, too complicated, we don not like it.");
#ifdef DEBUG
                        int edgesSize = int(edges.size());
#endif               
                        int faceEdgeID = AddAndSearchEdge(edges, v1, v2, face->id); 
#ifdef DEBUG         
                        if(edgesSize == faceEdgeID)Fatal("%s: edges size changed during searching", __func__);
#endif               
                        // Avoid double counting key nodes on face edge
                        if(!edges[faceEdgeID].ibs.empty()){
                            int n; for(n=0; n<int(edges[faceEdgeID].ibs.size());n++){
                                if(i==edges[faceEdgeID].ibs[n].inGrain && m==edges[faceEdgeID].ibs[n].ID)break;
                            }
                            if(n!=int(edges[faceEdgeID].ibs.size()))continue;
                        }

                        // assemble key nodes on face edge
                        EdgeJoint_t ib; real8 dir[3];
                        ib.inGrain = i; ib.ID=m; 
                        if(v1<v2){
                            CoordMinus(p,p1,dir);
                            ib.dis=(PointPointDistance(p,p1)/PointPointDistance(p1,p2));
                        }else{ 
                            CoordMinus(p,p2,dir);
                            ib.dis=(PointPointDistance(p,p2)/PointPointDistance(p1,p2)); 
                        }
#ifdef DEBUG
                        if(!Numerical(ib.dis))Fatal("%s: a intersection node of inner boundary and face edge is not numerical",__func__);
#endif
                        ib.sign = DotProduct(dir, grain->ibN)<0?-1.0:1.0;
                        edges[faceEdgeID].ibs.push_back(ib);
                        
                    }
                }

                // if the edge is found, store it and draw edge line
                if(edge.vertID2 != -1)ibEdges.push_back(edge);

                k+=glo_f_vert[k]+1;
            }
            grain->ibsNumEdges[m]=(int)ibEdges.size();
            std::vector<int> polygon; real8 centroid[3],dis;
            
            // Generate polygon for outputing
            EdgeListToPolygon(ibEdges,polygon); 

            // caculate centroid and areaof polygon and copy it to grain->ibPs
            CentroidOfPolygon(polygon, grainVertices, 0, centroid);
            grain->ibAreas[m] = AreaOfPolygon(polygon, grainVertices, 0);
            PointPlaneDistance(centroid, grain->ibN, grain->ibPs[m], dis, grain->ibPs[m]);

#ifdef DEBUG
            if(fabs(dis)>1E-10)Fatal("%s at line %d, the distance to centorid (%e) is too large",__FILE__,__LINE__, dis);
#endif

            real8 lDir[3], crpd[3];
            grain->ibsEdges[m] = (real8*)malloc(grain->ibsNumEdges[m]*6*sizeof(real8));
            grain->ibsEdgesFaceID[m] = (int *)malloc(grain->ibsNumEdges[m]*sizeof(int));
            grain->ibPolygons[m] = (int *)malloc(grain->ibsNumEdges[m]*sizeof(int));

            for(auto n=0; n<grain->ibsNumEdges[m]; n++){
                lDir[0] = grainVertices[ibEdges[n].vertID2*3] - grainVertices[ibEdges[n].vertID1*3];
                lDir[1] = grainVertices[ibEdges[n].vertID2*3+1] - grainVertices[ibEdges[n].vertID1*3+1];
                lDir[2] = grainVertices[ibEdges[n].vertID2*3+2] - grainVertices[ibEdges[n].vertID1*3+2];
                if(normal(lDir) < 1E-11){Fatal("the two grainVertices are too close");} normalize_vec(lDir);
                p[0] = 0.5*(grainVertices[ibEdges[n].vertID2*3] +   grainVertices[ibEdges[n].vertID1*3]);
                p[1] = 0.5*(grainVertices[ibEdges[n].vertID2*3+1] + grainVertices[ibEdges[n].vertID1*3+1]);
                p[2] = 0.5*(grainVertices[ibEdges[n].vertID2*3+2] + grainVertices[ibEdges[n].vertID1*3+2]);
                
                // We comsider a plane which is both perpendicular to inner boundary and line sense of the edge,
                // and cacualte the normal of plane (point to outside of the inner boundary)
                cross_product(grain->ibN, lDir, crpd); normalize_vec(crpd);
                if(crpd[0]*(grain->ibPs[m][0]-p[0])+crpd[1]*(grain->ibPs[m][1]-p[1])+crpd[2]*(grain->ibPs[m][2]-p[2])>0){
                    crpd[0] = -crpd[0];crpd[1] = -crpd[1];crpd[2] = -crpd[2];
                }
                VECTOR_COPY(&(grain->ibsEdges[m][n*6]), crpd); VECTOR_COPY(&(grain->ibsEdges[m][n*6+3]), p);
                
                // store the face ID that the edge on for boundary topology
                grain->ibsEdgesFaceID[m][n] = ibEdges[n].ofFaces[0];
            }
            
            if(home->myDomain ==0){
                char zoneName[256];std::vector<int> info(3);
                snprintf(zoneName,sizeof(zoneName),"IB-%d-%d",grain->id, m);
                info[0] = 2; info[1] = grain->id; info[2] = m;
                DrawPolygon(fp, polygon, grainVertices, 0, zoneName, info);
            }
            
            // store inner boundary polygon
            for(auto n=0; n<grain->ibsNumEdges[m]; n++)grain->ibPolygons[m][n]=polygon[n+1];
        }

        if(!grainVertices.empty()){
            grain->nIBVertices=int(grainVertices.size()/3);
            grain->ibVertices = (real8 *)malloc(3*grain->nIBVertices*sizeof(real8));
            grain->ibVertOnEdge = (int *)malloc(grain->nIBVertices*sizeof(int));
            grain->ibVertOnSubedge = (int *)malloc(grain->nIBVertices*sizeof(int));
            for(auto j=0; j<grain->nIBVertices;j++){
                grain->ibVertices[3*j]=grainVertices[3*j];
                grain->ibVertices[3*j+1]=grainVertices[3*j+1];
                grain->ibVertices[3*j+2]=grainVertices[3*j+2];
                grain->ibVertOnEdge[j]=-1;
                grain->ibVertOnSubedge[j]=-1;
            }
        }else grain->nIBVertices=0;

        grainVertices.clear();
        if(!cl.inc())break;
    }
    if(home->myDomain ==0){fclose(fp);}

    /*  copy vertices into home
     */
    home->nVertices = (int)gloVertices.size()/3; 
    home->vertices = (real8*)malloc(3*home->nVertices*sizeof(real8));
    for(auto i=0; i<3*home->nVertices; i++)home->vertices[i]=gloVertices[i];
    home->nEdges = int(edges.size());
    home->edges = (Edge_t **)malloc(home->nEdges*sizeof(Edge_t *));
    for(auto i=0;i<home->nEdges;i++){
        home->edges[i] = (Edge_t *)malloc(sizeof(Edge_t));
        sort(edges[i].ibs.begin(), edges[i].ibs.end(), CmpKeyNodeT);
        EdgeCopy(home->edges[i], edges[i]);

        real8 pos[3]; Edge_t *edge=home->edges[i];
        for(auto j=0; j<edge->nIBs; j++){
            CoordAdd(1.0-edge->dis[j], Coordinate(home->vertices, edge->vertID1),
                      edge->dis[j], Coordinate(home->vertices, edge->vertID2), pos);
            Grain_t *grain=home->grainKeys[edge->ibInGrain[j]];
            int k; for(k=0; k<grain->nIBVertices; k++){
                if(PointPointDistance(Coordinate(grain->ibVertices,k), pos)<FACE_ERROR)break;
            }
            if(k==grain->nIBVertices)Fatal("%s: can not claim the sub edge node", __func__);
            grain->ibVertOnEdge[k]=i; grain->ibVertOnSubedge[k]=j;
        }
    }

#ifdef DEBUG
    for(auto i=0; i<home->nGrains; i++){
        Grain_t *grain=home->grainKeys[i];
        for(auto j=0; j<grain->nIBVertices; j++){
            if(grain->ibVertOnEdge[j]==-1 || grain->ibVertOnSubedge[j]==-1){
                Fatal("%s: ibVertOnEdge and ibVertOnSubedge not defineded",__func__);
            }
        }
    }
#endif
   
    /** Assembly burg informations of GBs
     */ 
    for(auto i=0; i<home->nFaces; i++){
        Face_t *face=home->faces[i]; 
        InitBurgInfo(&face->bInfs[0][0]);
        InitBurgInfo(&face->bInfs[0][1]);
        InitBurgInfo(&face->bInfs[1][0]);
        InitBurgInfo(&face->bInfs[1][1]);
        if(face->type == FREE_SURFACE)continue;
        Grain_t *grain1=home->grainKeys[face->nbrGrain[0]];
        Grain_t *grain2=home->grainKeys[face->nbrGrain[1]];
        GenerateGBBurgInfo(&face->bInfs[0][0], &grain1->bInfs[0], face->n, &grain2->bInfs[0]);
        if(grain1->nIBs>0)GenerateGBBurgInfo(&face->bInfs[1][0], &grain1->bInfs[1], face->n, &grain2->bInfs[0]);
        if(grain2->nIBs>0)GenerateGBBurgInfo(&face->bInfs[0][1], &grain1->bInfs[0], face->n, &grain2->bInfs[1]);
        if(grain1->nIBs>0&&grain2->nIBs>0)GenerateGBBurgInfo(&face->bInfs[1][1], &grain1->bInfs[1], face->n, &grain2->bInfs[1]);
    }

    /** Assembly burg informations of IBs
     */ 
    for(auto i=0; i<home->nGrains; i++){
        Grain_t *grain=home->grainKeys[i];
        InitBurgInfo(&grain->ibBurgInfs[0]); InitBurgInfo(&grain->ibBurgInfs[1]);
        if(grain->nIBs==0)continue;
        GenerateGBBurgInfo(&grain->ibBurgInfs[0], &grain->bInfs[0], grain->ibN, &grain->bInfs[1]);
        GenerateGBBurgInfo(&grain->ibBurgInfs[1], &grain->bInfs[1], grain->ibN, &grain->bInfs[0]);
    }

    if(home->myDomain == 0){
        if ((fp = fopen("face-edges.plt", "w")) == (FILE *)NULL) {
             Fatal("%s: Open error %d on %s\n",__func__, errno, "face-edges.plt");
        }
        fprintf(fp,"variables = x, y, z, vx, vy, vz, edge-type\n");
        std::vector<int> indexs(1);
        for(auto i=0; i<int(edges.size()); i++){
            real8 pos1[3], pos2[3];
            indexs[0]= 0;
            Coordinate(gloVertices, edges[i].vertID1, pos1);
            Coordinate(gloVertices, edges[i].vertID2, pos2);
            DrawLine(fp, pos1, pos2, 10, indexs);
        }   

        for(auto i=0; i<home->nGrains; i++){
            Grain_t *grain = home->grainKeys[i];
            real8 *p1, *p2; indexs[0]=1;
            for(auto j=0; j<grain->nIBs; j++){
                for(auto k=0; k<grain->ibsNumEdges[j]; k++){
                    p1=Coordinate(grain->ibVertices, grain->ibPolygons[j][k]);
                    p2=(k<grain->ibsNumEdges[j]-1)?
                        Coordinate(grain->ibVertices, grain->ibPolygons[j][k+1]):
                        Coordinate(grain->ibVertices, grain->ibPolygons[j][0]);
                    DrawLine(fp, p1, p2, 10, indexs);
                }
            }
        }

        fclose(fp);

#ifdef DEBUG
        int place[PLACE_LEN]; real8 pos[3];

        if ((fp = fopen("keypoints.plt", "w")) == (FILE *)NULL) {
            Fatal("%s: Open error %d on %s\n", __func__,errno, "keypoints.plt");
        }

        fprintf(fp,"variables = x, y, z, pointType, place[0], place[1], place[2], place[3]\n");

        /** output centroid of grain, face, inner boundary with pointType 0, 1, 2, respectivetly */
        for(auto i=0; i<home->nGrains; i++){
            Grain_t *grain = home->grainKeys[i];

            VECTOR_COPY(pos,  grain->pos);
            NodeLocation(home, pos, place, FACE_ERROR);
            fprintf(fp, "%g %g %g %d %d %d %d %d\n", pos[0], pos[1], pos[2],
                    0, place[0], place[1], place[2], place[3]);

            for(auto j=0; j<grain->nFaces; j++){
                VECTOR_COPY(pos,  grain->faces[j]->p);
                NodeLocation(home, pos, place, FACE_ERROR);
                fprintf(fp, "%g %g %g %d %d %d %d %d\n", pos[0], pos[1], pos[2],
                        1, place[0], place[1], place[2], place[3]);
            }
            for(auto j=0; j<grain->nIBs; j++){
                VECTOR_COPY(pos, grain->ibPs[j]);
                NodeLocation(home, pos, place, FACE_ERROR);
                fprintf(fp, "%g %g %g %d %d %d %d %d\n", pos[0], pos[1], pos[2],
                        2, place[0], place[1], place[2], place[3]);
            }
        }

        /** output all vertices with pointType 3 */
        for(auto i=0; i<home->nVertices; i++){
            Coordinate(home->vertices, i, pos);
            NodeLocation(home, pos, place, FACE_ERROR);
            fprintf(fp, "%g %g %g %d %d %d %d %d\n", pos[0], pos[1], pos[2],
                    3, place[0], place[1], place[2], place[3]);
        }

        for(auto i=0; i<home->nEdges; i++){
            Edge_t *edge = home->edges[i];
            for(auto j=0; j<edge->nIBs; j++){
                CoordAdd(1.0-edge->dis[j], Coordinate(home->vertices,edge->vertID1), 
                             edge->dis[j], Coordinate(home->vertices,edge->vertID2), pos);
                NodeLocation(home, pos, place, FACE_ERROR);
                fprintf(fp, "%g %g %g %d %d %d %d %d\n", pos[0], pos[1], pos[2],
                        4, place[0], place[1], place[2], place[3]);
            }
        }
        fclose(fp);
#endif
    }

#ifdef GSL
    GSLTest();
#endif
    return;
}
#ifdef DEBUG
    #undef NodeLocation
#endif
//}
/*******************************************************End of InitGrain*************************************************************/

/** -----------------------------------------------------------------------------------------------------
 *  Table: define the location or arm
 *  grain ID: range [0,home->nGrains)
 *  sub-grain (SG) ID: range [-grain->nIBs-1,-1] in grain, [-edge->nIBs-1,-1] on face edge
 *  inner boundary (IB) ID: range [0,grain->nIBs) in grain, [0,edge->nIBs) on face edge
 *  face ID: range [0,home->nFaces)
 *  free surface (FS) ID: range [-6,-1]
 *  edge ID: range [0, home->nEdges)
 *  ----------------------------------------------------------------------------------------------
 *  TYPE (place[0])   | place[1]   | place[2]       | place[3]              | note
 *  ------------------|------------|----------------|-----------------------|---------------------
 *  IN_GRAIN          | grain ID   | SG ID          |                       |
 *  ON_IB             | grain ID   | IB ID1         |                       |
 *  ON_FACE           | face ID    | SG (IB) ID1    | SG (IB) ID2 (if has)  | 
 *  ON_FACE_EDGE      | edge ID    | SG (IB) ID     |                       | immoveable if it has any arm in bulk
 *  AT_VERTICE        | vertice ID |                |                       | immoveable
 *  -------------------------------------------------------------------------------------------------------*/
void FormatPlace(int *place){
    switch(place[0]){
        case IN_GRAIN:
            printf("In grain %d, sub grain %d\n", place[1], place[2]); break;
        case ON_IB:
            printf("On inner boundary, in grain %d, IB index %d\n", place[1], place[2]); break;
        case ON_FACE:
            printf("On face %d", place[1]);
            for(auto i=2; i<4; i++){
                if(place[i]<0)printf(", in neighbor grain %d's subgrain %d",i-2,place[i]);
                else printf(", on neghbor grain %d's inner edge %d", i-2,place[i]);
            }printf("\n");break;
        case ON_FACE_EDGE:
            printf("on face edge %d, sub edge id %d\n", place[1], place[2]); break;
        case AT_VERTICE:
            printf("at vertice %d\n", place[1]); break;
        default:
            printf("Warning for format place\n");break;
    }
    return;
}

/** According to the current node coordinates, roughly judge the node position.
  * \param[in] pos the coordinates of node
  * \param[out] place the rough location of the node
  * \param[in] error the error, absoule margin of face and edge
  * Note that the point should not be outside of the geometry */
void NodeLocation(Home_t *home, real8 *pos, int *place, real8 error)
{
    int i, j, k, edgeID=0;
    bool breakLoop=false; real8 dis, interPos[3];
    Grain_t *grain; Face_t *face;

    real8 rx, ry, rz; int grainID; Grain_t **grainList; int grainListSize, nGrains;
    if(home->con->find_voronoi_cell(pos[0], pos[1], pos[2], rx, ry, rz, grainID)){
        grain = home->grainKeys[grainID];
        grainListSize=grain->nFaces+1; nGrains=grainListSize;
        grainList = (Grain_t **)malloc((grain->nFaces+1)*sizeof(Grain_t *));
        grainList[0] = grain;
        for(i=0; i<grain->nFaces; i++){
            face = grain->faces[i];
            if(face->type==FREE_SURFACE){grainList[i+1]=(Grain_t *)NULL; continue;}
            if(grain->id != face->nbrGrain[0])grainList[i+1]=home->grainKeys[face->nbrGrain[0]];
            else grainList[i+1]=home->grainKeys[face->nbrGrain[1]]; 
        }
    }else {
        grainListSize = 0;
        nGrains = home->nGrains; 
        grainList = home->grainKeys;

#ifdef DEBUG
#if    DEBUG==2
        printf("Warning: can not find the voronoicell of position (%g,%g,%g)\n", pos[0], pos[1], pos[2]);
#endif
#endif
    }

    // initialize the location type
    *place = -1; *(place+3)=-1;
    // Loop through all grains
    for(i=0; i<nGrains; i++){
        grain = grainList[i];
        // Loop through all faces of the grain
        for(j=0; j<grain->nFaces; j++){
            face=grain->faces[j];
            PointPlaneDistance(pos, face->n, face->p, dis);
            dis *= grain->sign[j];

            // if the displacement to face is lager than tolerance, check if the position is inside of the face
            if(fabs(dis) < error){
                breakLoop=true;
                if(POINT_INTERSECTION!=ClipSegmentInConcavePolyhedron(face->p, pos, face->edges, face->nEdges, -error, edgeID, interPos)){
                    // if there is no intersection on inner edge, the face is found.
                    if(breakLoop){place[0] = ON_FACE; place[1]=face->id;}
                }else{
                    int edgeID2;
                    // if there is an intersection on inner edge and no intersection on outer edge, the point is on the edge or at vertice.
                    if(POINT_INTERSECTION!=ClipSegmentInConcavePolyhedron(face->p, pos, face->edges, face->nEdges, error, edgeID2, interPos)){
                    // if there is an intersection on inner edge and no intersection on outer edge, the point is on the edge or at vertice.
                        // if the node is at one end of edge, the node should be at vertice and return func.
                        Edge_t *edge = home->edges[face->edgesID[edgeID]];
                        real8 len=PointPointDistance(Coordinate(home->vertices,edge->vertID2), Coordinate(home->vertices,edge->vertID1));
                        real8 localError=error/len;
                        real8 dis=PointPointDistance(pos, Coordinate(home->vertices,edge->vertID1))/len;
#ifdef  DEBUG
                        if(!Numerical(dis))Fatal("%s, variable dis is not numerical.",__func__);
#endif
                        if(dis<localError){place[0]=AT_VERTICE; place[1]=edge->vertID1;} // at vertice 1 of the edge
                        else if(dis>1.0-localError){place[0]=AT_VERTICE; place[1]=edge->vertID2;} // at vertice 2 of the edge
                        else {
                            place[0]=ON_FACE_EDGE; place[1]=face->edgesID[edgeID]; place[2]=InverseIndex(edge->nIBs); // at the face edge
                            for(k=0; k<edge->nIBs; k++){
                                if(fabs(dis-edge->dis[k])<localError){place[2]=k; break;} //at node k of the face edge
                                else if(dis<edge->dis[k]){place[2]=InverseIndex(k); break;} // at subedge -1-k of the face edge
                                else continue;
                            }
                        }

                        if(grainListSize>0)free(grainList);return;
                    }else{
                        // if the point is out side of outer edge, the loop continues to find next face
                        breakLoop = false;
                    }
                }
            }
            // if the distance between the node and plane is greater than error or the face which the node located on has been found, 
            // break out the inner loop.
            if(dis > error || breakLoop)break;
        }
        // if the face is found, break out the outer loop
        if(breakLoop)break;

        // if the node is in the polyhedron, place[0] and place[1] can be defined.
        if(j==grain->nFaces){place[0] = IN_GRAIN; place[1]=grain->id; break;}
    }

#ifdef DEBUG
    if(i==nGrains)Fatal("%s: can not find location of (%g,%g,%g)",__func__,*pos,pos[1], pos[2]);
#endif

    // Depending on the location type, continue to determine the subgrain or IB ID
    switch(place[0]){
        case IN_GRAIN:
            grain=home->grainKeys[place[1]];
            for(i=0; i<grain->nIBs; i++){
                PointPlaneDistance(pos, grain->ibN, grain->ibPs[i], dis);
                //if the distance is within the margin of error, determine IB index.
                if(fabs(dis)<error){place[0]=ON_IB; place[2]= i; break;}
                if(dis<0)break;
            }
            if(place[0]==IN_GRAIN)place[2]=InverseIndex(i);
            break;
        case ON_FACE:
            face=home->faces[place[1]];
            for(i=0; i<2; i++){
                // skip the second neighbor grain that the node on free surface
                if(face->nbrGrain[i]<0){place[2+i]=face->nbrGrain[i]; continue;}
                grain=home->grainKeys[face->nbrGrain[i]];
                
                place[2+i]=InverseIndex(grain->nIBs); // initialize place[2+i] to top subgrain ID
                for(j=0; j<grain->nIBs; j++){
                    PointPlaneDistance(pos, grain->ibN, grain->ibPs[j], dis);
                    // if the node is both on face and on IB, it is on inner boundary edge in ith nbrGrain.
                    if(fabs(dis)<error){place[2+i]=j; break;}
                    // else if the point is outside of cutting IB, it is in subgrain
                    else if(dis<0){place[2+i]=InverseIndex(j); break;}
                }
            }
            break;
        default:
            Fatal("can not supose this location type %d", place[0]);break;
    }
    if(grainListSize>0)free(grainList);return;
}

/** According to the current node coordinates, judge if the node is the grain and subgrain
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error inner error of grain
 * \return false if the position are not at the location, true otherwise. */
bool InGrain(Home_t *home, int *place, real8 *pos, real8 error)
{
    int i; real8 dis;
#ifdef DEBUG
    if(place[0]!=IN_GRAIN)Fatal("%s: place[0] = %d",__func__, place[0]);
    if(place[1]<0 || place[1]>home->nGrains)Fatal("%s: place[1] = %d",__func__, place[1]);
#endif    
    // make sure the node is in the grain
    Grain_t *grain = home->grainKeys[place[1]];
#ifdef DEBUG
    if(place[2]<-grain->nIBs-1||place[2]>-1)Fatal("%s: no sub-grain %d in grain %d",__func__,place[2], place[1]);
#endif
    /*Check if the node is outside of inner boundary firstly*/
    int ibID = InverseIndex(place[2]);
    if(ibID<grain->nIBs){
        PointPlaneDistance(pos, grain->ibN, grain->ibPs[ibID], dis);
        if(dis > -error)return false;
    }
    if(ibID>0){
        PointPlaneDistance(pos, grain->ibN, grain->ibPs[ibID-1], dis);
        if(dis < error)return false;
    }

    return InGrain(grain, pos, error);
}

/** According to the current position of the node, determine whether the node 
 *  is on the face, on the subface or on the IB edge.
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error inner error of the face
 * \return false if the position are not at the location, true otherwise. */
bool OnFace(Home_t *home, int *place, real8 *pos, real8 error)
{
#ifdef DEBUG
    if(place[1]<0 || place[1]>=home->nFaces)Fatal("%s: wrong face ID", __func__);
#endif
    Face_t *face = home->faces[place[1]];
    real8  dis;

    // First, if the point is not on face with error of FACE_ERROR, false
    PointPlaneDistance(pos, face->n, face->p, dis);
    if(fabs(dis)>FACE_ERROR)return false;

    // Second, check if the point matchs of inner edges of face
    int ibID;  Grain_t *grain;
    for(auto i=0; i<2; i++){
        if(face->nbrGrain[i]<0)continue;
        grain = home->grainKeys[face->nbrGrain[i]];
    
        // If on a IB, locate IB directly and check if it is on that IB.
        if(place[2+i]>=0){
            PointPlaneDistance(pos, grain->ibN, grain->ibPs[place[2+i]], dis);
            if(fabs(dis)>FACE_ERROR)return false;
        }else{
            ibID = InverseIndex(place[2+i]);
            if(ibID<grain->nIBs){  // if has upper boundary
                PointPlaneDistance(pos, grain->ibN, grain->ibPs[ibID], dis);
                if(dis>-error)return false;
            }
            if(ibID>0){ // if has lower boundary
                PointPlaneDistance(pos, grain->ibN, grain->ibPs[ibID-1], dis);
                if(dis<error)return false;
            }
        }
    }

    // Last, if the point is outside of inner edges of the face, false
    if(!InConcavePolyhedron(pos, face->edges, face->nEdges, -error))return false;

    return true;
}


/** According to the current position of the node, determine whether the position
 *  is on the inner boundary
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error error for determination of inner boundary
 * \return false if the position are not at the location, true otherwise. */
bool OnIB(Home_t *home, int *place, real8 *pos, real8 error)
{
#ifdef DEBUG
    if(place[0]!=ON_IB)Fatal("%s: wrong place",__func__);
    if(place[1]<0 || place[1]>=home->nGrains)Fatal("%s: the grain index %d is not correct", __func__, place[1]);
#endif
    Grain_t *grain = home->grainKeys[place[1]];
#ifdef DEBUG
    if(place[2]<0 || place[2]>=grain->nIBs)Fatal("%s: the IB index %d is not correct",__func__, place[2]);
#endif    
    
    int &ibID=place[2]; real8 dis;
    PointPlaneDistance(pos, grain->ibN, grain->ibPs[ibID], dis);
    if(fabs(dis) > error)return false;

    if(!InConcavePolyhedron(pos, grain->ibsEdges[ibID], grain->ibsNumEdges[ibID], -error)) return false;

    return true;
}

/** According to the current position of the node, determine whether the position
 *  is on the face edge
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error error for determination of edge and joint
 * \return false if the position are not at the location, true otherwise. */
bool OnFaceEdge(Home_t *home, int *place, real8 *pos, real8 error)
{
#ifdef DEBUG
    if(place[0] != ON_FACE_EDGE)Fatal("%s line %d: place[0]=%d",__func__,__LINE__,place[0]);
    if(place[1]<0 || place[1]>home->nEdges)Fatal("%s line %d: place[1]=%d", __func__,__LINE__,place[1]); 
#endif
    Edge_t *edge = home->edges[place[1]];
    real8  *p1, *p2, vec[3], lDir[3], dis;
    p1=Coordinate(home->vertices, edge->vertID1);
    p2=Coordinate(home->vertices, edge->vertID2); 
    CoordMinus(pos, p1, vec); CoordMinus(p2, p1, lDir);
  
#ifdef DEBUG
    if(normal(lDir)<FACE_ERROR)Fatal("%s: two vertice are too close",__func__);
#endif
    // if the point is not in range of the segment, it is definitely not on the edge.
    real8 len=normal(lDir), dotRatio=DotProduct(lDir, vec)/len;
    if(dotRatio < error || dotRatio > len-error)return false;

    // if the point is not at the line, it is definitely not on the edge. 
    normalize_vec(lDir);
    PointLineIntersection(pos, lDir, p1, vec, &dis); 
    if(fabs(dis)>fabs(error))return false;

#ifdef DEBUG
    if(place[2] >= edge->nIBs || place[2]<-edge->nIBs-1)Fatal("%s: the id %d is not in the range of edge subindexs",__func__, place[2]);
#endif
    if(place[2]>-1)return (fabs(dotRatio-edge->dis[place[2]]*len)<error);
    else {
        int ibID=InverseIndex(place[2]);
        if(ibID<edge->nIBs)if(edge->dis[ibID]*len-dotRatio > error)return false;
        if(ibID>0)if(edge->dis[ibID-1]*len-dotRatio < -error)return false;
    }
    return true;
}

/** According to the current position of the node, determine whether the position
 *  is at the vertice
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error error for determination of edge and joint
 * \return false if the position are not at the location, true otherwise. */
bool AtVertice(Home_t *home, int *place, real8 *pos, real8 error)
{
#ifdef DEBUG
    if(place[0] != AT_VERTICE)Fatal("%s: wrong place %d",__func__, place[0]);
    if(place[1]<0 || place[1]>= home->nVertices)Fatal("%s: wrong vertice ID %s",__func__, place[1]);
#endif
    return (PointPointDistance(Coordinate(home->vertices, place[1]), pos)<error);
}

/** According to the current position of the node, determine whether the position
 *  is at the location
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error error for determination of edge and joint
 * \return false if the position are not at the location, true otherwise. */
bool AtLocation(Home_t *home, int *place, real8 *pos, real8 error)
{
    switch (place[0]){
        case IN_GRAIN: return InGrain(home, place, pos, error);
        case ON_IB: return OnIB(home, place, pos, error); 
        case ON_FACE: return OnFace(home, place, pos, error); 
        case ON_FACE_EDGE: return OnFaceEdge(home, place, pos, error); 
        case AT_VERTICE: return AtVertice(home, place, pos, error);  
        default: Fatal("%s: there is no location type %d",__func__, place[0]);
    }
    return false;
}

/** Find the intersection between nodal path and grain (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of the node (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathInGrain(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace)
{
    Face_t *face;  bool status=false; real8 nFace[3];

#ifdef DEBUG
    if(place[0]!=IN_GRAIN)Fatal("%s: place[0] = %d",__func__, place[0]);
    if(place[1]<0 || place[1]>=home->nGrains)Fatal("%s: place[1] = %d",__func__, place[1]);
    if(!InGrain(home, place, pos1))Fatal("%s: the position 1 is not in corrected grain.",__func__);
#endif    
    Grain_t *grain = home->grainKeys[place[1]];
    int i, ibID, edgeID, tempEdgeID;
    real8 tempPos[3];
    for(i=0; i<grain->nFaces; i++){
        face=grain->faces[i];
        if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, face->n, face->p, pos2)){
            VECTOR_COPY(tempPos, pos2);
            // if the position 2 is within the inner boundary, position 2 is on the face
            if(POINT_INTERSECTION != ClipSegmentInConcavePolyhedron(face->p, tempPos, face->edges, face->nEdges, -FACE_ERROR, edgeID)){
                status = true;
                newPlace[0] = ON_FACE;
                newPlace[1] = face->id;
                // loop through two neighbor grains of the face
                for(auto j=0; j<2; j++){
                    // if the face is free surface (no second neighbor grain), save the surface ID intp newPlace and continue
                    if(face->nbrGrain[j]<0){newPlace[2+j]=face->nbrGrain[j]; continue;}
                
                    // if the neighbor grain is current grain, determine whether the segment crosses the two inner boundaries
                    if(face->nbrGrain[j]==grain->id){
#ifdef DEBUG    
                        if(place[2]<-grain->nIBs-1||place[2]>-1)Fatal("%s: no sub-grain %d in grain %d",__func__,place[2], place[1]);
#endif          
                        newPlace[2+j]=place[2]; //the default is old subgrain index
                        ibID = InverseIndex(place[2]);
                        
                        if(ibID<grain->nIBs){ // if has upper innner boundary
                            if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], pos2)){
                                // if the path is clipped on inner boundary now, there are two posibilities of the node: on IB or on IB edge
                                int localEdgeID;
                                if(ClipSegmentInConcavePolyhedron(grain->ibPs[ibID], pos2, grain->ibsEdges[ibID], 
                                   grain->ibsNumEdges[ibID], localEdgeID)==POINT_INTERSECTION){
                                    // if there is a intersection to IB edge, the node is on IB edge
                                    newPlace[2+j]=ibID; continue;
                                }else{
                                    newPlace[0]=ON_IB; newPlace[1]=grain->id; newPlace[2]=ibID; return status;
                                }
                            }
                        }
                
                        if(ibID>0){ // if has lower boundary
                            ibID--;
                            if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], pos2)){
                                // if the path is clipped on inner boundary now, there are two posibilities of node: on IB or on IB edge
                                int localEdgeID;
                                if(ClipSegmentInConcavePolyhedron(grain->ibPs[ibID], pos2, grain->ibsEdges[ibID], 
                                   grain->ibsNumEdges[ibID], localEdgeID)==POINT_INTERSECTION){
                                    // if there is a intersection to IB edge, the node is on IB edge
                                    newPlace[2+j]=ibID; continue;
                                }else{
                                    newPlace[0]=ON_IB; newPlace[1]=grain->id; newPlace[2]=ibID; return status;
                                }
                            }
                        }
                    } /*if the neighbor is the adjacent grain*/ else {
                        real8 dis;
                        Grain_t *nbrGrain = home->grainKeys[face->nbrGrain[j]];
                        newPlace[2+j] = InverseIndex(nbrGrain->nIBs);
                        for(auto k=0; k<nbrGrain->nIBs; k++){
                            PointPlaneDistance(pos2, nbrGrain->ibN, nbrGrain->ibPs[k], dis);
                            if(dis<0){newPlace[2+j]=InverseIndex(k); break;}
                        }
                    }
                }
                return status;
            }
            /* if there is an intersection on inner face edge*/ 
            else {
                VECTOR_COPY(tempPos, pos2);
                // if no intersection outer face edge, it means that position 2 is on the edge already
                if(POINT_INTERSECTION != ClipSegmentInConcavePolyhedron(face->p, tempPos, face->edges, face->nEdges, FACE_ERROR, tempEdgeID)){
                    Edge_t *edge = home->edges[face->edgesID[edgeID]];
                    real8 len=PointPointDistance(Coordinate(home->vertices,edge->vertID2), Coordinate(home->vertices,edge->vertID1));
                    real8 localError=FACE_ERROR/len;
                    real8 dis=PointPointDistance(pos2, Coordinate(home->vertices,edge->vertID1))/len;
#ifdef  DEBUG
                    if(!Numerical(dis))Fatal("%s, variable dis is not numerical.",__func__);
#endif
                    if(dis<localError){newPlace[0]=AT_VERTICE; newPlace[1]=edge->vertID1;} // at vertice 1 of the edge
                    else if(dis>1.0-localError){newPlace[0]=AT_VERTICE; newPlace[1]=edge->vertID2;} // at vertice 2 of the edge
                    else {
                        newPlace[0]=ON_FACE_EDGE; newPlace[1]=face->edgesID[edgeID]; newPlace[2]=InverseIndex(edge->nIBs); // at the face edge
                        for(auto k=0; k<edge->nIBs; k++){
                            if(fabs(dis-edge->dis[k])<localError){newPlace[2]=k; break;} //at node k of the face edge
                            else if(dis<edge->dis[k]){newPlace[2]=InverseIndex(k); break;} // at subedge -1-k of the face edge
                            else continue;
                        }
                    }
                    return true;
                }
                // else if there has intersection on outer edge, means that the position 2 is not in the face, just continue it to loop through another face
            }
        }else{
            // if there is no intersection on the face, keep looking for intersection on inner margin with -SAFE_FACE_ERROR
            // and then clip the path and make sure that the posistion of node is not too close to the face.
            nFace[0] = grain->sign[i]*face->n[0]; nFace[1] = grain->sign[i]*face->n[1]; nFace[2] = grain->sign[i]*face->n[2];
            SegmentPlaneIntersection(pos1, pos2, nFace, face->p, -SAFE_FACE_ERROR, pos2);
        }
    }

    ibID = InverseIndex(place[2]); 
#ifdef DEBUG
    if(ibID<0||ibID>grain->nIBs)Fatal("%s: no subgrain %d in grain %d",__func__,place[2], place[1]);
#endif
    if(ibID<grain->nIBs){ // if has upper boundary
        if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], pos2)){
            newPlace[0]=ON_IB; newPlace[1]=grain->id; newPlace[2]=ibID;
            return true; 
        }

        // if there is an intersection point on inner plane of subgrain with a margin of -SAFE_FACE_ERROR
        // clip the segment without changing place
        SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], -SAFE_FACE_ERROR, pos2);
    }

    if(ibID>0){ // if has lower boundary
        ibID--;
        if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], pos2)){
            newPlace[0]=ON_IB; newPlace[1]=grain->id; newPlace[2]=ibID;
            return true; 
        }

        // if there is an intersection point on inner plane of subgrain with a margin of SAFE_FACE_ERROR
        // clip the segment without changing place
        SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], SAFE_FACE_ERROR, pos2);
    }

    return status;
}


/** Find the intersection between nodal path and face (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position of the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the node, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of intersection (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathOnFace(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace)
{  
    int edgeID, ibID; real8 dis;Grain_t *grain;
#ifdef DEBUG
    if(place[1]<0 || place[1]>=home->nFaces)Fatal("%s: wrong face ID", __func__);
#endif
    Face_t *face = home->faces[place[1]];
#ifdef DEBUG
    if(!OnFace(home, place, pos1[0], pos1[1], pos1[2]))Fatal("%s: the first position (%f,%f,%f) is not on the face.",
        __func__, pos1[0], pos1[1], pos1[2]);
#endif

    // Make sure the current position is on the face, clip the path if point is not.
    PointPlaneDistance(pos2, face->n, face->p, dis, pos2);
#ifdef DEBUG
    if(fabs(dis) > SAFE_FACE_ERROR)Fatal("%s: the point 2 is far away to the face, distance %e", __func__, dis);
#endif
   
    // check if the path has crossed the face boundary 
    if(POINT_INTERSECTION==ClipSegmentInConcavePolyhedron(pos1, pos2, face->edges, face->nEdges, edgeID, pos2)){
#ifdef DEBUG
        if(face->edgesID[edgeID] < 0 || face->edgesID[edgeID]>=home->nEdges)Fatal("%s line %d: wrong global edge ID",__func__, __LINE__);
#endif
        // the clipping point may be on the face edge or at the vertice of the edge
        Edge_t *edge = home->edges[face->edgesID[edgeID]];
        real8  *p1, *p2, len; 
        p1=Coordinate(home->vertices, edge->vertID1); 
        p2=Coordinate(home->vertices, edge->vertID2);
        len=PointPointDistance(p1,p2);
        real8 localError=FACE_ERROR/len;
        real8 dis=PointPointDistance(pos2,p1)/len;

#ifdef  DEBUG
        if(!Numerical(dis))Fatal("%s, variable dis is not numerical.",__func__);
#endif
        if(dis<localError){newPlace[0]=AT_VERTICE; newPlace[1]=edge->vertID1;} // at vertice 1 of the edge
        else if(dis>1.0-localError){newPlace[0]=AT_VERTICE; newPlace[1]=edge->vertID2;} // at vertice 2 of the edge
        else {
            newPlace[0]=ON_FACE_EDGE; newPlace[1]=face->edgesID[edgeID]; newPlace[2]=InverseIndex(edge->nIBs); // at the face edge
            for(auto k=0; k<edge->nIBs; k++){
                if(fabs(dis-edge->dis[k])<localError){newPlace[2]=k; break;} //at node k of the face edge
                else if(dis<edge->dis[k]){newPlace[2]=InverseIndex(k); break;} // at subedge -1-k of the face edge
                else continue;
            }
        }
        return true;
    }
    // if there has no intersection on face edges, clip the path by inner margins to ensure there is 
    // no misunderstanding of the place.
    else ClipSegmentInConcavePolyhedron(pos1, pos2, face->edges, face->nEdges,-SAFE_FACE_ERROR, edgeID, pos2);
    
    // if the node is neither on face edge nor at vertice, initialize newPlace, it may not be used here, but we have to plan for it
    newPlace[0] = ON_FACE; newPlace[1]=face->id;

    bool status = false;
    // loop through grains on both side of the face, clip the path acoording to the nodal place
    for(auto i=0; i<2; i++){
        // if no neighbor grain on this side (face is free surface), copy the surface ID and continue
        if(face->nbrGrain[i]<0){newPlace[2+i]=face->nbrGrain[i]; continue;}

        grain = home->grainKeys[face->nbrGrain[i]];
        // if the index is IB index already, clip the sgement by the inner boundary
        if(place[2+i]>=0){
            SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[place[2+i]], pos2);
        }else /* if the segment is also in subgrain of the neighbor */{
            ibID = InverseIndex(place[2+i]);

            if(ibID<grain->nIBs){ // if has upper boundary
                if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], pos2)){
                    newPlace[2+i] = ibID; status=true; continue;
                }
                // if ther is an intersection point on inner plane of subgrain with a margin of -SAFE_FACE_ERROR
                // clip the place without changing place
                else SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], -SAFE_FACE_ERROR, pos2);
            }

            if(ibID>0){ // if has lower boundary
                ibID--; 
                if(POINT_INTERSECTION == SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], pos2)){
                    newPlace[2+i] = ibID; status=true; continue;
                }
                // if ther is an intersection point on inner plane of subgrain with a margin of SAFE_FACE_ERROR
                // clip the path without changing place
                else SegmentPlaneIntersection(pos1, pos2, grain->ibN, grain->ibPs[ibID], SAFE_FACE_ERROR, pos2);
            }
            newPlace[2+i] = place[2+i];
        }
    }/*end for(i=0;i<2;i++) */
    return status;
}

/** Find the intersection between nodal path and inner boundary (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position of the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the node, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of intersection (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathOnIB(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace)
{
#ifdef DEBUG
    if(place[0] != ON_IB)Fatal("%s: wrong place type %d",__func__, place[0]);
    if(place[1]<0 || place[1]>= home->nGrains)Fatal("%s: can not find grain %d",__func__, place[1]);
    if(!OnIB(home, place, pos1, FACE_ERROR))Fatal("%s: the position 1 is not on the IB",__func__);
#endif
    Grain_t *grain = home->grainKeys[place[1]];
#ifdef DEBUG
    if(place[2]<0 || place[2]>=grain->nIBs)Fatal("%s: can not find ib %d in grain %d",__func__, place[2], place[1]);
#endif    
    real8 dis;
    int edgeID, ibID=place[2];
    PointPlaneDistance(pos2, grain->ibN, grain->ibPs[ibID], dis, pos2);
#ifdef DEBUG
    if(fabs(dis) > SAFE_FACE_ERROR)Fatal("%s: the point 2 is far away to the inner boundary, distance %e", __func__, dis);
#endif

    // if there is an intersection on edge of inner boundary
    if(POINT_INTERSECTION==ClipSegmentInConcavePolyhedron(pos1, pos2, grain->ibsEdges[ibID], grain->ibsNumEdges[ibID], edgeID, pos2)){
        int vertID1=grain->ibPolygons[ibID][edgeID];
        int vertID2 = (edgeID==grain->ibsNumEdges[ibID]-1)?grain->ibPolygons[ibID][0]:grain->ibPolygons[ibID][edgeID+1];
        int vertID = -1;
        if(PointPointDistance(pos2, Coordinate(grain->ibVertices, vertID1))<FACE_ERROR)vertID = vertID1;
        if(PointPointDistance(pos2, Coordinate(grain->ibVertices, vertID2))<FACE_ERROR)vertID = vertID2;
        if(vertID != -1){
            newPlace[0] = ON_FACE_EDGE; newPlace[1] = grain->ibVertOnEdge[vertID]; 
            newPlace[2]=grain->ibVertOnSubedge[vertID]; newPlace[3]=-1;
        }else {
            newPlace[0]=ON_FACE; 
            newPlace[1]=grain->ibsEdgesFaceID[ibID][edgeID];
            Face_t *face=home->faces[newPlace[1]];
            for(auto i=0; i<2; i++){
                if(face->nbrGrain[i]<0){newPlace[2+i]=face->nbrGrain[i]; continue;}
                newPlace[i+2] = WhichSubgrain(home->grainKeys[face->nbrGrain[i]], pos2, FACE_ERROR);
            }
        }
        return true;
    }
    // else, clip the path to make sure that pos is on the IB
    else ClipSegmentInConcavePolyhedron(pos1, pos2, grain->ibsEdges[ibID], grain->ibsNumEdges[ibID],-SAFE_FACE_ERROR, edgeID, pos2);
    return false;
}

/** Find the intersection between nodal path and face edge (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position of the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the node, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of intersection (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathOnFaceEdge(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace)
{
#ifdef DEBUG
    if(place[0] != ON_FACE_EDGE)Fatal("%s: wrong location type %d", __func__, place[0]);
    if(place[1]<0 || place[1]>=home->nEdges)Fatal("%s: wrong edge ID %d", __func__, place[1]);
    if(!OnFaceEdge(home,place,pos1,FACE_ERROR))Fatal("%s: the position 1 is not on face edge",__func__);
#endif
    Edge_t *edge = home->edges[place[1]];
    real8 *p1 = Coordinate(home->vertices, edge->vertID1); 
    real8 *p2 = Coordinate(home->vertices, edge->vertID2); 
    real8 lDir[3], len; CoordMinus(p2,p1,lDir); len=normal(lDir);
    normalize_vec(lDir);

    real8 dis; PointLineIntersection(pos2, lDir, p1, pos2, &dis);
#ifdef DEBUG
    if(dis>SAFE_FACE_ERROR)Fatal("%s: the position 2 is far away from the edge (%e)",__func__, dis);
#endif
    dis = PointPointDistance(p1, pos2)/len;

    // if the place is at the joint of edge
    if(place[2]<0){
        int jointID = InverseIndex(place[2]);
        real8 localErr=FACE_ERROR/len, upperDis = 1.0, lowerDis=0.0;
        if(jointID < edge->nIBs)upperDis = edge->dis[jointID];
        else if(jointID > 0)lowerDis = edge->dis[jointID-1];

        if(dis-upperDis>-localErr){
            if(jointID==edge->nIBs){newPlace[0]=AT_VERTICE; newPlace[1]=edge->vertID2;}
            else {newPlace[0]=ON_FACE_EDGE; newPlace[1]=place[1]; newPlace[2]=jointID; newPlace[3]=-1;}
            CoordAdd(1.0-upperDis, p1, upperDis, p2, pos2);
            return true;
        }

        if(dis-lowerDis<localErr){
            if(jointID==0){newPlace[0]==AT_VERTICE; newPlace[1]=edge->vertID1;}
            else {newPlace[0]=ON_FACE_EDGE; newPlace[1]=place[1]; newPlace[2]=jointID-1;newPlace[3]=-1;}
            CoordAdd(1.0-lowerDis, p1, lowerDis, p2, pos2);
            return true;
        }
        
    }else{
        if(fabs(dis-edge->dis[place[2]])>FACE_ERROR){
            Fatal("%s: joint on edge can not is inmoveable, move distance %e",__func__, dis-edge->dis[place[2]]);
        }
    }

    return false;
}

/** Find the intersection between nodal path and vertice (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position of the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the node, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of intersection (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathAtVertice(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace)
{
#ifdef DEBUG
    if(place[0]!=AT_VERTICE)Fatal("%s: wrong location type %d",__func__,place[0]);
    if(place[1]<0 || place[1]>=home->nVertices)Fatal("%s: can not find vertice %d",__func__,place[1]);
    if(!AtVertice(home, place, pos1, FACE_ERROR))Fatal("%s: position 1 is not at vertice",__func__);
#endif

    real8 dis, *p=Coordinate(home->vertices, place[1]);
    if((dis=PointPointDistance(p, pos2))>FACE_ERROR){
#ifdef DEBUG
        if(dis>SAFE_FACE_ERROR)Fatal("%s: distance of the path is %e", __func__, dis);
#endif
        VECTOR_COPY(pos2,p);
    }
    return true;
}

/*-------------------------------------ArmLocation Functions ---------------------------------*/
/** Get adjucent grains of the face
 *  \param[in] faceID face index
 *  \param[in] id1,id2 subface indexs
 *  \param[out] grains index list of grains and subgrains */
void AdjacentGrainsOfFace(Home_t *home, int faceID, int id1, int id2, std::vector<int> &grains)
{   
#ifdef DEBUG
    if(faceID<0 || faceID>= home->nFaces)Fatal("%s: can not find face %d",__func__, faceID);
#endif
    grains.clear(); int place[2]={id1,id2};
    Face_t *face = home->faces[faceID];
    for(int i=0;i<2;i++){
        if(face->nbrGrain[i]<0)continue;
        if(place[i]<0){
            grains.push_back(face->nbrGrain[i]); grains.push_back(InverseIndex(place[i]));
        }else{
            grains.push_back(face->nbrGrain[i]); grains.push_back(place[i]);
            grains.push_back(face->nbrGrain[i]); grains.push_back(place[i]+1);
        }
    }
    return;
}

/** Get adjucent faces of the face edge
 *  \param[in] edgeID edge index
 *  \param[in] subEdgeID subedge index
 *  \param[out] faces index list of faces and subfaces */
void AdjacentFacesOfFaceEdges(Home_t *home, int edgeID, int subEdgeID, std::vector<int> &faces)
{
#ifdef DEBUG
    if(edgeID<0 || edgeID>=home->nEdges)Fatal("%s: can not find edge %d",__func__, edgeID);
#endif
    faces.clear();
    Edge_t *edge = home->edges[edgeID];
#ifdef DEBUG
    if(subEdgeID<edge->nIBs-1 || subEdgeID>=edge->nIBs)
        Fatal("%s: can not find sub edge %d in edge %d",__func__,subEdgeID, edgeID);
#endif


    if(subEdgeID>=0){
        int ibInGrain = edge->ibInGrain[subEdgeID];
        int ibID = edge->ibID[subEdgeID];
        real8 p[3]; int j,k,ids[2];
        CoordAdd(1.0-edge->dis[subEdgeID],Coordinate(home->vertices, edge->vertID1),
                     edge->dis[subEdgeID],Coordinate(home->vertices, edge->vertID2), p);                
        Face_t *face;
        for(int i=0; i<edge->nFaces; i++){
            face = home->faces[edge->ofFaces[i]];
            for(j=0; j<2; j++){
                if(face->nbrGrain[i]==ibInGrain){
                    ids[j]=ibID; k=(j==0)?1:0;
                    if(face->nbrGrain[k]<0)ids[k]=face->nbrGrain[k];
                    else ids[k]=WhichSubgrain(home->grainKeys[face->nbrGrain[k]], p, FACE_ERROR);
                    faces.push_back(face->id); faces.push_back(ids[0]); faces.push_back(ids[1]);
                    break;
                }
            }
        }
#if 0 // unfinish
    }else{
        int ibID = InverseIndex(edge->ibID[subEdgeID]);
        real8 upperDis = (ibID<edge->nIBs)?edge->dis[ibID]:1.0;
        real8 lowerDis = (ibID>0)?edge->dis[ibID-1]:0.0;
#endif   
    }
    return;
}

inline void ArmLocationGrainGrain(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    char msg[256]; 
    if(p1[1]!=p2[1]){snprintf(msg,sizeof(msg),"%s: the arm intersects grain boundary",__func__); throw msg;}
    if(p1[2]!=p2[2]){snprintf(msg,sizeof(msg),"%s: the arm intersects inner boundary",__func__); throw msg;}
#endif
    return;
}
inline void ArmLocationGrainIB(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    char msg[256]; 
    if(p1[1]!=p2[1]){snprintf(msg,sizeof(msg),"%s: the arm intersects grain boundary",__func__); throw msg;}
    int ibID=InverseIndex(p1[2]);
    if(ibID<p2[2] || ibID>p2[2]+1){snprintf(msg,sizeof(msg),"%s: the arm intersects inner boundary",__func__); throw msg;} 
#endif
    return;
}
inline void ArmLocationGrainFace(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    char msg[256]; std::vector<int> grains;
    AdjacentGrainsOfFace(home, p2[1], p2[2], p2[3], grains);
    unsigned int i=0; for(; i<grains.size(); i=i+2){
        if(grains[i]==p1[1] && grains[i+1]==p1[2])break;
    }
    if(i==grains.size()){
        snprintf(msg,sizeof(msg),"%s: can not find grain (%d,%d) in adjacent grain list of face (%d,%d,%d)",
                 __func__, p1[1],p1[2],p2[1],p2[2],p2[3]);
        throw msg;
    }
#endif
    return;
}
inline void ArmLocationGrainFaceEdge(Home_t *home, int *p1, int *p2, int *p)
{
    
    return;
}
inline void ArmLocationGrainVertice(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}

inline void ArmLocationIBIB(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    if(p1[1]!=p2[1] || p1[2]!=p2[2]){
        char msg[256];
        snprintf(msg,sizeof(msg),"%s: two inner boundaries (%d,%d,%d,%d) are not same",
                 __func__, p1[1],p1[2],p2[1],p2[2]);
        throw msg;
    }
#endif
    return;
}
inline void ArmLocationIBFace(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}
inline void ArmLocationIBFaceEdge(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}
inline void ArmLocationIBVertice(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}

inline void ArmLocationFaceFace(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    char msg[256];
    if(p1[1]!=p2[1]){
        snprintf(msg,sizeof(msg),"%s: two different face %d %d",__func__,p1[1],p2[1]);
        throw msg;
    }
    if(p1[1]<0 || p1[1]>=home->nFaces){
        snprintf(msg,sizeof(msg),"%s: no face %d",__func__,p1[1]);
        throw msg;
    }
#endif
    if(p1[2]==p2[2] && p1[3]==p2[3])return;
#ifdef DEBUG
    for(int i=2; i<4; i++){
        if(p2[i]<0 && p1[i]>-1)if(InverseIndex(p2[i])<p1[i]+1 && InverseIndex(p2[i])>=p1[i])continue;
        if(p1[i]<0 && p2[i]>-1)if(InverseIndex(p1[i])<p2[i]+1 && InverseIndex(p1[i])>=p2[i])continue;
        snprintf(msg,sizeof(msg),"%s: the arm interscts inner boundary %d %d",__func__,p1[i],p2[i]);
        throw msg;
    }
#endif
    p[2] = p1[2]<=p2[2]?p1[2]:p2[2];
    p[3] = p1[3]<=p2[3]?p1[3]:p2[3];

    return;
}
inline void ArmLocationFaceFaceEdge(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}
inline void ArmLocationFaceVertice(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}


inline void ArmLocationFaceEdgeFaceEdge(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    char msg[256];
    if(p1[1]!=p2[1]){
        snprintf(msg,sizeof(msg),"%s: the arms (%d,%d)--(%d,%d) intersects two edges",
            __func__,p1[1],p1[2],p1[2],p2[2]);
        throw msg;
    }
    if(p1[2]>-1 && p2[2]>-1){
        snprintf(msg,sizeof(msg),"%s: the ends of arm (%d--%d) is on subedges",
                 __func__,p1[2],p2[2]);
        throw msg;
    }
#endif
    p[2]=p1[2]<=p2[2]?p1[2]:p2[2];
    return;
}
inline void ArmLocationFaceEdgeVertice(Home_t *home, int *p1, int *p2, int *p)
{
    return;
}

inline void ArmLocationVerticeVertice(Home_t *home, int *p1, int *p2, int *p)
{
#ifdef DEBUG
    char msg[256];
    snprintf(msg,sizeof(msg),"%s: the arms (%d,%d)--(%d,%d) at vertice",
        __func__,p1[1],p1[2],p1[2],p2[2]);
    throw msg;
#endif
    return;
}

/** Get Arm location 
 *  \param[in] place1, place2 two locations of the arm ends
 *  \param[out] place arm location*/
void ArmLocation(Home_t *home, int *place1, int *place2, int *place)
{
    int *p1, *p2;
    if(place1[0]<=place2[0]){p1=place1; p2=place2;}else{p1=place2; p2=place1;}
    for(auto i=0; i<PLACE_LEN; i++)place[i]=p1[i];

    try{
        switch(p1[0]){
            case IN_GRAIN:
                switch(p2[0]){
                    case IN_GRAIN: ArmLocationGrainGrain(home,p1,p2,place); break;
                    case ON_IB: ArmLocationGrainIB(home,p1,p2,place); break;
                    case ON_FACE: ArmLocationGrainFace(home,p1,p2,place); break;
                    case ON_FACE_EDGE: ArmLocationGrainFaceEdge(home,p1,p2,place); break;
                    case AT_VERTICE: ArmLocationGrainVertice(home,p1,p2,place); break;
                    default: Fatal("%s at line %d: not support location type %d",__func__,__LINE__, p2[0]); break;
                }break;
            case ON_IB:
                switch(p2[0]){
                    case ON_IB: ArmLocationIBIB(home,p1,p2,place); break;
                    case ON_FACE: ArmLocationIBFace(home,p1,p2,place); break;
                    case ON_FACE_EDGE: ArmLocationIBFaceEdge(home,p1,p2,place); break;
                    case AT_VERTICE: ArmLocationIBVertice(home,p1,p2,place); break;
                    default: Fatal("%s at line %d: not support location type %d",__func__,__LINE__, p2[0]); break;
                }break;
        
            case ON_FACE:
                switch(p2[0]){
                    case ON_FACE: ArmLocationFaceFace(home,p1,p2,place); break;
                    case ON_FACE_EDGE: ArmLocationFaceFaceEdge(home,p1,p2,place); break;
                    case AT_VERTICE: ArmLocationFaceVertice(home,p1,p2,place); break;
                    default: Fatal("%s at line %d: not support location type %d",__func__,__LINE__, p2[0]); break;
                }break;
            case ON_FACE_EDGE:
                switch(p2[0]){
                    case ON_FACE_EDGE: ArmLocationFaceEdgeFaceEdge(home,p1,p2,place); break;
                    case AT_VERTICE: ArmLocationFaceEdgeVertice(home,p1,p2,place); break;
                    default: Fatal("%s at line %d: not support location type %d",__func__,__LINE__, p2[0]); break;
                }break;
            case AT_VERTICE:
                switch(p2[0]){
                    case AT_VERTICE: ArmLocationVerticeVertice(home,p1,p2,place); break;
                    default: Fatal("%s at line %d: not support location type %d",__func__,__LINE__, p2[0]); break;
                }break;
            default: Fatal("%s at line %d: not support location type %d",__func__,__LINE__, p1[0]); break;
                
        }
    }catch(const char *msg){ 

    }
    return;
}
/*---------------------------------End of ArmLocation Functions ---------------------------------*/

/** Get Burgers information of arm
 * \param[in] place arm location
 * \return burg information pointer */
BurgInfo_t *BurgInformation(Home_t *home, int *place)
{
    int ibID, subgrainID;
    Grain_t *grain; Face_t *face; Edge_t *edge;
    switch(place[0]){
        case IN_GRAIN: return &(home->grainKeys[place[1]]->bInfs[InverseIndex(place[2])%2]); break;
        case ON_IB: return &(home->grainKeys[place[1]]->ibBurgInfs[place[2]%2]); break;
        case ON_FACE:
            face = home->faces[place[1]]; 
            if(face->type==FREE_SURFACE){
                grain = home->grainKeys[face->nbrGrain[0]];
                if(place[2]>-1) return &(grain->ibBurgInfs[place[2]%2]);
                else return &(grain->bInfs[InverseIndex(place[2])%2]);
            }else{
                if(place[2]>-1 || place[3]>-1){
                    return (place[2]>-1?&(home->grainKeys[face->nbrGrain[0]]->ibBurgInfs[place[2]%2]):
                            &(home->grainKeys[face->nbrGrain[1]]->ibBurgInfs[place[3]%2]));
                } else return &(face->bInfs[InverseIndex(place[2])%2][InverseIndex(place[3])%2]);
            }
            break;
        case ON_FACE_EDGE:
            edge=home->edges[place[1]];
#ifdef DEBUG 
            if(place[2]>-1)Fatal("%s: arm location can not at edge joint", __func__);
#endif
            ibID = InverseIndex(place[3]);
            if(ibID == edge->nIBs){
                grain = home->grainKeys[edge->ibInGrain[ibID-1]];
                subgrainID = edge->ibSign[ibID-1]==1.0?edge->ibID[ibID-1]+1:edge->ibID[ibID-1];
            }else{
                grain = home->grainKeys[edge->ibInGrain[ibID]];
                subgrainID = edge->ibSign[ibID]==1.0?edge->ibID[ibID]:edge->ibID[ibID]+1;
            }
            return &grain->bInfs[subgrainID%2];
        default:
            Fatal("%s: can not support arm location type %d",__func__,place[0]); break;
    }
    return (BurgInfo_t *)NULL;
}

/** Match Burgers vector that have aligned direction 
 *  \return return matched Burgers vector index*/
int MatchBurgID(BurgInfo_t *bInf, real8 *burg)
{
    real8 mBurg[3], b[3]={burg[0],burg[1],burg[2]};
    normalize_vec(b);
    for(int i=0; i<bInf->numBurgVectors; i++){
        VECTOR_COPY(mBurg, bInf->burgList[i]);
        normalize_vec(mBurg);
        if(AlignedDirection(mBurg, b))return i;
    }
    return bInf->numBurgVectors;
}

/** Get glide constraint of the node
 *  \param[in] node the node to be caculated
 *  \param[out] vec constraint normal (plane constraint) or direction (line direction)
 *  \return constrain type 0: no constraint; 1: plane constraint; 2: line constraint; 3: fixed*/
inline int AddToPlaneConstraint(real8 *v, real8 nx, real8 ny, real8 nz){
    if(!AlignedDirection(v[0],v[1],v[2],nx,ny,nz)){
        cross_product(v[0],v[1],v[2],nx,ny,nz,v);
        normalize_vec(v); return 2;
    }else return 1;
}
inline int AddToLineConstraint(real8 *v, real8 nx, real8 ny, real8 nz){
    return PerpendicularDirection(v[0],v[1],v[2],nx,ny,nz)==true?2:3;
}
inline bool FindKeyLine(BurgInfo_t *bInf, real8 nx, real8 ny, real8 nz, real8 *kDir)
{
    for(auto i=0; i<bInf->numKeyLines; i++)
        if(PerpendicularDirection(bInf->keyLineList[i][0],bInf->keyLineList[i][1],bInf->keyLineList[i][2],
                                  nx, ny, nz)){VECTOR_COPY(kDir, bInf->keyLineList[i]); return true;}
    return false;
}
int GlideConstraint(Home_t *home, Node_t *node, real8 *vec)
{
    // Pinned node 
    if(node->constraint == PINNED_NODE)return 3;

    real8 lDir[3], kDir[3];
    int i, j=1, jArm = -1;
    Node_t *nbr; BurgInfo_t *bInf;
    bool jDefined = false;
    Param_t *param=home->param;

#ifdef DEBUG
    if(node->numNbrs==0)Fatal("%s: something wrong with the arm",__func__);
#endif

    // Loop through arms, find the first junction arm (if has)
    for(i=0; i<node->numNbrs; i++){
#ifdef DEBUG
        /* Make sure that normal is always perpendicular to Burgers*/
        if(!PerpendicularDirection(node->nx[i],node->ny[i],node->nz[i],
            node->burgX[i],node->burgY[i],node->burgZ[i])){
            PrintNode(home, node);
            Fatal("%s: arm %d, the Burgers vector is not perpendicular to the normal",__func__,i);
        }
#endif
        /* We need neighbor node to get the Burgers information, if we can not find the neighbor, 
         * what might come up something unexpected.*/
        if((nbr = GetNeighborNode(home, node, i))==(Node_t *)NULL){
            PrintNode(home,node); 
            if(node->myTag.domainID==home->myDomain){
                Fatal("%s: can not find the neighbor %d of the native node",__func__,i);
            }
            /* FIXME: If node is ghost, the unfounded neighbor should be ghost or secondary ghost node.
             * To avoid more serious problems, completlely constraint the node*/
            else{
                printf("Task %d, Warning in %s, can not find the secondary ghost neighbor node (%d)\n",
                        home->myDomain,__func__,i);
                return 3;
            }
        }
        
        // Get the Burgers information
        bInf = BurgInformation(home, node->place, nbr->place);

        // if the plane is not glide plane, the arm should b a junction, but we are not sure
        // that the line sense of this arm can be defined
        if(GLIDE_PLANE != PlaneType(bInf, node->nx[i], node->ny[i], node->nz[i])){
            
            // loop through all key lines in the data base, find the key line on the junction plane
            if(FindKeyLine(bInf, node->nx[i], node->ny[i], node->nz[i], kDir)){
#ifdef DEBUG
                if(param->allowFuzzyGlidePlanes==0){
                    lDir[0]=nbr->x-node->x; lDir[1]=nbr->y-node->y; lDir[2]=nbr->z-node->z;
                    if(DotProduct(lDir,lDir)>1E-15){
                        normalize_vec(lDir);
                        if(!AlignedDirection(lDir, kDir)){
                            Fatal("%s: the junction has climbed without allowing fuzzy glide planes",__func__);
                        }
                    }
                }
#endif
                // if the line constraint is found, the junction can be defined.
                if(jDefined==false){VECTOR_COPY(vec, kDir); jDefined=true; jArm=i;}

                // if there has more than 2 junctions with different line constraints, means that
                // the node is fixed.
                else if(!AlignedDirection(kDir,vec))return 3;

            } /* if the key line has been found */  else jArm = i;
        } /*if the junction plane has been found */
    }/* for(i=0; i<node->numNbrs; i++) */
    
#ifdef DEBUG
    if(jDefined==false && jArm >= 0){
        PrintNode(home, node);
        Fatal("%s: can not define junction plane");
    } 
#endif

    // if there is no junction 
    if(jArm<0){
        // assign first normal as initial plane constraint
        vec[0]=node->nx[0]; vec[1]=node->ny[0]; vec[2]=node->nz[0];
        for(i=1; i<node->numNbrs; i++){
            j=AddToPlaneConstraint(vec, node->nx[i], node->ny[i], node->nz[i]); 
            if(j==2){
                for(auto k=i+1; k<node->numNbrs; k++){
                    j=AddToLineConstraint(vec,node->nx[k], node->ny[k], node->nz[k]);
                    if(j==3)return 3;
                }
                break;
            }
        }
    }
    // if junction has been found line constraint has been defined, just add other glide constraints to line constraint
    else if(jDefined==true){
        for(i=0; i<node->numNbrs; i++)
            if(3==AddToLineConstraint(vec, node->nx[i], node->ny[i], node->nz[i])) return 3;
        j=2;
    }
    // if the junction arm has been found that it intersects with surface (only one arm),
    else if(jDefined==false && node->numNbrs==1){
        return 3;
    }
    else Fatal("%s: there is fourth possiblities!",__func__);

    // Add extra location constraint to glide constraint
    Grain_t *grain; Face_t *face; Edge_t *edge; int *place = node->place;
    switch(place[0]){
        case IN_GRAIN: break;
        case ON_IB: 
#ifdef DEBUG
            if(place[1]<0 || place[1]>=home->nGrains)Fatal("%s: can not find grain %d", __func__, place[1]);
#endif
            grain=home->grainKeys[place[1]];
            if(j==1) j=AddToPlaneConstraint(vec, grain->ibN[0], grain->ibN[1], grain->ibN[2]);
            else j = AddToLineConstraint(vec, grain->ibN[0], grain->ibN[1], grain->ibN[2]);
            break;
        case ON_FACE:
#ifdef DEBUG
            if(place[1]<0 || place[1]>=home->nFaces)Fatal("%s: can not find face %d", __func__, place[1]);
#endif
            face=home->faces[place[1]];
            if(j==1) j=AddToPlaneConstraint(vec, face->n[0], face->n[1], face->n[1]);
            else j = AddToLineConstraint(vec, face->n[0], face->n[1], face->n[1]);
            break;
        case ON_FACE_EDGE:
#ifdef DEBUG
            if(place[1]<0 || place[1]>=home->nEdges)Fatal("%s: can not find edge %d", __func__, place[1]);
#endif
            edge = home->edges[place[1]];
            CoordMinus(Coordinate(home->vertices, edge->vertID2), Coordinate(home->vertices, edge->vertID1), lDir);
            normalize_vec(lDir);
            if(j==1){if(!PerpendicularDirection(lDir, vec))return 3; else{j=2; VECTOR_COPY(vec,lDir);}}
            else return (AlignedDirection(lDir,vec)==true)?2:3;
            break;
        case AT_VERTICE: return 3; break;
        default: Fatal("%s: not support location type %d",__func__, place[0]); break;
    }
    return j;
}

/** Find plane of the arm, if there are several glide planes, we chose the glide plane with 
 *  maximal Schmid factor if FEM is disabled, or chose the glide plane with highest local dissipation
 *  if FEM is abled.
 *  \param[in] burg Burgers vector of the arm
 *  \param[in,out] lDir line direction of the arm. Note: no climb component in the line direction
 *  \param[in] node1, node the ends node of the arm
 *  \param[out] plane plane of the arm 
 *  \param[out] found true is found, otherwise false*/
void FindGlidePlane(Home_t *home, real8 *burg, real8 *lDir, Node_t *node1, Node_t *node2, real8 *plane, bool &found)
{
    found = false;
    real8 burgNorm[3]={burg[0],burg[1],burg[2]}; normalize_vec(burgNorm);
    BurgInfo_t  *bInf = BurgInformation(home, node1->place, node2->place); 

    // check glide constraints
    real8 v1[3], v2[3], *v; int t1, t2, t;
    t1 = GlideConstraint(home, node1, v1);
    t2 = GlideConstraint(home, node2, v2);
    
    char msg[256];
    if(t1==t2){
        if(t1==1){
            if(!AlignedDirection(v1,v2)){
                snprintf(msg, sizeof(msg),"%s at %d: two different glide planes",__func__,__LINE__); 
                throw msg;
            }
        }else if(t1==2){
            if(!AlignedDirection(v1,v2)){
                snprintf(msg, sizeof(msg),"%s at %d: two different line constraints",__func__,__LINE__);
                throw msg;
            }
        }else if(t1==3){
            snprintf(msg, sizeof(msg),"%s at %d: the two nodes are both pinned",__func__,__LINE__);
            throw msg;
        }
    }else if((t1==1&&t2==2) || (t2==1&&t1==2)){
        if(!PerpendicularDirection(v1,v2)){
            snprintf(msg, sizeof(msg),"%s at %d: line constraint is not perpendicular to glide plane",__func__,__LINE__);
            throw msg;
        }
    }
    
    // the arm constraint type should be the minimum of the nodal constraints
    if(t2<t1){ t=t2; v=v2; }else{t=t1; v=v1;}

    // if the plane constraint, just copy the plane
    if(t==1){
        VECTOR_COPY(plane, v);
#ifdef DEBUG
        /* If the line direction is defined */
        if(IsNormalized(lDir)){
            if(!PerpendicularDirection(v, lDir)){
                Fatal("%s: line direction does not match the plane constraint",__func__);
            }
        }
        if(!PerpendicularDirection(v,burgNorm)){
            Fatal("%s: Burgers direction does not match the plane constraint",__func__);
        }
#endif
        found=true; return;
    }

    // if the line constraint (junction), just cross product to get the glide plane
    else if (t==2){
        cross_product(v, burgNorm, plane);
#ifdef DEBUG
        if(IsNormalized(lDir)){
            if(!AlignedDirection(v,lDir)){
                Fatal("%s: line direction does not match the line constraint",__func__);
            }
        }
        if(normal(plane)<1E-2){
            Fatal("%s: a screw junction is occured",__func__);
        }
#endif
        normalize_vec(plane);
        int planeID = PlaneID(bInf, plane);
        if(planeID < bInf->numPlanes){VECTOR_COPY(plane, bInf->planeList[planeID]);}
#ifdef DEBUG
        else printf("warning in %s: there is not junction plane (%g,%g,%g) in Burger data base\n",__func__,plane[0],plane[1],plane[2]);
#endif
        found=true; return;
    }

    return;
}

/** Get the distance between point and polygon
 *  \param[in] p the coordinate of the point
 *  \param[in] v coordinates lis of vertices
 *  \param[in edges indexs and factor of polygon vertices
 *  \return distance                    */
double PointConcavePolygonDistance(double *p, std::vector<double> &v, std::vector<EdgeVec_t> &edges){
    real8 p1[3], p2[3], min=1E50, dis;
    for(auto i=0; i<int(edges.size()); i++){
        Coordinate(v, edges[i].vertID1, p1);
        Coordinate(v, edges[i].vertID2, p2);
        dis = (PointSegmentDistance(p, p1, p2)*edges[i].factor);
        min = (min>dis)?dis:min;
    }
    return min;
}

/** Get edges of intersection between plane and grain
 *  \param[in] n the normal vector of plane
 *  \param[in] p a point on the plane, this point has to be in the grain.
 *  \param[in] place the grin location
 *  \param[out] vertice coordinates of vertices
 *  \param[out] edges edge list     */
void PlaneGrainIntersection(Home_t *home, real8 *n, real8 *p, int *place, std::vector<real8> &vertices, std::vector<EdgeVec_t> &edges)
{
    vertices.clear(); edges.clear();

    int ibID;
    real8 *p1, *p2, interPos[3];
    real8 pos1[3], pos2[3], pos[3], dot1, dot2, t;
    Face_t *face; Edge_t *gloEdge;
    Grain_t *grain=home->grainKeys[place[1]];
    for(auto i=0; i<grain->nFaces; i++){
        face = grain->faces[i];
        EdgeVec_t itEdge; itEdge.vertID1=-1; itEdge.vertID2=-1; 
        itEdge.factor= (face->type==FREE_SURFACE)?1.0:0.5;
        itEdge.ofFaces.push_back(face->id);
        for(auto j=0; j<face->nEdges; j++){
            gloEdge = home->edges[face->edgesID[j]];
            p1 = Coordinate(home->vertices, gloEdge->vertID1);
            p2 = Coordinate(home->vertices, gloEdge->vertID2);

            if(POINT_INTERSECTION==SegmentPlaneIntersection(p1,p2,n,p,interPos)){
                if(itEdge.vertID1==-1){
                    itEdge.vertID1 = AddAndSearchVertice(vertices, interPos[0], interPos[1], interPos[2]);
                }else if (itEdge.vertID2==-1){
                    itEdge.vertID2 = AddAndSearchVertice(vertices, interPos[0], interPos[1], interPos[2]);
                }else Fatal("There is third intersection point, too complicated, we don not like it.");
            }
        }

        if(itEdge.vertID2!=-1){
            if(grain->nIBs>0){
                ibID=InverseIndex(place[2]);
                Coordinate(vertices, itEdge.vertID1, pos1);
                Coordinate(vertices, itEdge.vertID2, pos2);

                if(ibID<grain->nIBs){ // if has upper boundary
                    dot1 = grain->ibN[0]*(pos1[0]-grain->ibPs[ibID][0])+grain->ibN[1]*(pos1[1]-grain->ibPs[ibID][1])
                          +grain->ibN[2]*(pos1[2]-grain->ibPs[ibID][2]);
                    dot2 = grain->ibN[0]*(pos2[0]-grain->ibPs[ibID][0])+grain->ibN[1]*(pos2[1]-grain->ibPs[ibID][1])
                          +grain->ibN[2]*(pos2[2]-grain->ibPs[ibID][2]);
                    if(dot1>0 && dot2>0)continue;
                    if(dot1*dot2<0){
                        t=fabs(dot1)/(fabs(dot1)+fabs(dot2));
                        CoordAdd(1.0-t, pos1, t, pos2, pos);
                        if(dot1>0){VECTOR_COPY(pos1, pos);}
                        else {VECTOR_COPY(pos2, pos);}
                    }
                }
  
                if(ibID>0){ // if has lower boundary
                    ibID--;
                    dot1 = grain->ibN[0]*(pos1[0]-grain->ibPs[ibID][0])+grain->ibN[1]*(pos1[1]-grain->ibPs[ibID][1])
                          +grain->ibN[2]*(pos1[2]-grain->ibPs[ibID][2]);
                    dot2 = grain->ibN[0]*(pos2[0]-grain->ibPs[ibID][0])+grain->ibN[1]*(pos2[1]-grain->ibPs[ibID][1])
                          +grain->ibN[2]*(pos2[2]-grain->ibPs[ibID][2]);
                    if(dot1<0 && dot2<0)continue;
                    if(dot1*dot2<0){
                        t=fabs(dot1)/(fabs(dot1)+fabs(dot2));
                        CoordAdd(1.0-t, pos1, t, pos2, pos);
                        if(dot1<0){VECTOR_COPY(pos1, pos);}
                        else {VECTOR_COPY(pos2, pos);}
                    }
                }
                vertices[3*itEdge.vertID1]=pos1[0];
                vertices[3*itEdge.vertID1+1]=pos1[1];
                vertices[3*itEdge.vertID1+2]=pos1[2];
                vertices[3*itEdge.vertID2]=pos2[0];
                vertices[3*itEdge.vertID2+1]=pos2[1];
                vertices[3*itEdge.vertID2+2]=pos2[2];
            }
            edges.push_back(itEdge);
        }
    }

    if(grain->nIBs>0){
        real8 intersection[6];
        ibID=InverseIndex(place[2]);
        if(ibID<grain->nIBs){ // if has upper boundary
            if(PlaneConcavePolygonIntersection(n,p,grain->ibVertices,
                grain->ibPolygons[ibID],grain->ibsNumEdges[ibID],intersection)){
                EdgeVec_t itEdge;  itEdge.factor=0.5;
                itEdge.vertID1=AddAndSearchVertice(vertices, intersection[0], intersection[1], intersection[2]);
                itEdge.vertID2=AddAndSearchVertice(vertices, intersection[3], intersection[4], intersection[5]);
                edges.push_back(itEdge);
            }
        }
  
        if(ibID>0){ // if has lower boundary
            ibID--;
            if(PlaneConcavePolygonIntersection(n,p,grain->ibVertices,
                grain->ibPolygons[ibID],grain->ibsNumEdges[ibID],intersection)){
                EdgeVec_t itEdge; itEdge.factor=0.5; 
                itEdge.vertID1=AddAndSearchVertice(vertices, intersection[0], intersection[1], intersection[2]);
                itEdge.vertID2=AddAndSearchVertice(vertices, intersection[3], intersection[4], intersection[5]);
                edges.push_back(itEdge);
            }
        }
    }
    return;
}

void GenerateFCCBurgInfo(real8 rotMatrix[3][3], BurgInfo_t *burgData)
{
        int         i;

        burgData->numBurgVectors = NUM_BURGS_IN_FCC;
        burgData->numPlanes   = NUM_PLANES_IN_FCC;
        burgData->numKeyLines = NUM_KEY_LINES_IN_FCC;

        burgData->burgList = (real8 (*)[3])malloc(3*burgData->numBurgVectors*sizeof(real8));
        burgData->burgType = (int *)malloc(burgData->numBurgVectors*sizeof(int));
        burgData->planeList = (real8 (*)[3])malloc(3*burgData->numPlanes*sizeof(real8));
        burgData->planeType = (int *)malloc(burgData->numPlanes*sizeof(int));
        burgData->cutNormals = (real8 (*)[3])malloc(3*burgData->numPlanes*sizeof(real8));
        burgData->keyLineList = (real8 (*)[3])malloc(3*burgData->numKeyLines*sizeof(real8));

        burgData->numPlanesPerBurg = (int *)malloc(burgData->numBurgVectors*sizeof(int));        
        burgData->numBurgsPerPlane = (int *)malloc(burgData->numPlanes*sizeof(int));        
        burgData->numPlanesPerKeyLine = (int *)malloc(burgData->numKeyLines*sizeof(int));        
        

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
            VECTOR_ZERO(burgData->cutNormals[i]);
        }
        
        int planeType[NUM_PLANES_IN_FCC] = { GLIDE_PLANE, GLIDE_PLANE, GLIDE_PLANE, GLIDE_PLANE,
            JUNCTION_PLANE1, JUNCTION_PLANE1, JUNCTION_PLANE1, JUNCTION_PLANE2, 
            JUNCTION_PLANE2, JUNCTION_PLANE2, JUNCTION_PLANE2, JUNCTION_PLANE2, JUNCTION_PLANE2 }; 

        memcpy(burgData->planeType,planeType,NUM_PLANES_IN_FCC*sizeof(int));

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
        burgData->burgOnPlaneID = (int **)malloc(NUM_BURGS_IN_FCC*sizeof(int *));       
        for(i=0; i<NUM_BURGS_IN_FCC; i++){
            burgData->burgOnPlaneID[i] = (int *)malloc(numPlanesPerBurg[i]*sizeof(int));
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

        burgData->planeOnBurgID = (int **)malloc(NUM_PLANES_IN_FCC* sizeof(int *));       
        for(i=0; i<NUM_PLANES_IN_FCC; i++){
            burgData->planeOnBurgID[i] = (int *)malloc(numBurgsPerPlane[i]* sizeof(int));
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
        burgData->keyLineOnPlaneID = (int **)malloc(NUM_KEY_LINES_IN_FCC* sizeof(int *));       
        for(i=0; i<NUM_KEY_LINES_IN_FCC; i++){
            burgData->keyLineOnPlaneID[i] = (int *)malloc(numPlanesPerKeyLine[i]* sizeof(int));
            memcpy(burgData->keyLineOnPlaneID[i], keyLineOnPlaneID+keyLineFirstPlaneIndex[i], 
                   sizeof(int [numPlanesPerKeyLine[i]]));
        }
        return;
}

/** Generate Burgers data information
 *  \param[in] matType material type
 *  \param[in] rotMatrix rotation matrix 
 *  \param[out] burgData Burgers data */
void GenerateBurgInfo(int matType, real8 rotMatrix[3][3], BurgInfo_t *burgData){
    switch (matType){
        case MAT_TYPE_FCC:
            GenerateFCCBurgInfo(rotMatrix, burgData);            
            break;
        default:
            Fatal("%s: can not support this material type %d", matType);
    }
    return;
}

void CompressBurgListWithType(std::vector<std::vector<real8> >& burgList, std::vector<int>& burgType)
{
    int     i, j, nRows;
    real8   crossProdt[3], vec1[3], vec2[3];
    int     newNRows = 0;
    std::vector< std::vector<real8> >::iterator  it;
    std::vector<int>::iterator  it2;

    nRows = burgList.size();
    if(nRows != int(burgType.size()))Fatal("in %s at %d", __FILE__, __LINE__);
    for(i=0; i<nRows; i++){
        for(j=0; j<newNRows; j++){
            VECTOR_COPY(vec1, burgList[i]);
            VECTOR_COPY(vec2, burgList[j]);
            cross_product(vec1, vec2, crossProdt);
            if(AlignedVector(vec1,vec2)==true &&  burgType[i] == burgType[j])break;

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
/** Assemble GB Burgers information by bonding Burgers information of adjacent grains at grain boundary
 *  \param[out] bInf Burgers informatiob of the GB
 *  \param[in] bInf1 bInf2 two Burgers information of adjacent grains
 *  \param[in] GBNormal normal of grain boundary for cacualte cutting plane */
void GenerateGBBurgInfo(BurgInfo_t *bInf, BurgInfo_t *bInf1, real8 *GBNormal, BurgInfo_t *bInf2)
{
    /** Assign plane Information and key line Info */
    bInf->numPlanes = bInf1->numPlanes + bInf2->numPlanes;
    bInf->planeList = (real8 (*)[3])malloc(3*bInf->numPlanes*sizeof(real8));
    memcpy(bInf->planeList, bInf1->planeList, bInf1->numPlanes*sizeof(real8 [3]));
    memcpy(bInf->planeList + bInf1->numPlanes, bInf2->planeList, bInf2->numPlanes*sizeof(real8 [3]));

    bInf->planeType = (int *)malloc(bInf->numPlanes*sizeof(int));
    memcpy(bInf->planeType, bInf1->planeType, bInf1->numPlanes*sizeof(int));
    memcpy(bInf->planeType + bInf1->numPlanes, bInf2->planeType, bInf2->numPlanes*sizeof(int));

    real8 crossProdt[3], dot;
    std::vector<std::vector<double> >   keyLineList; 
    std::vector<real8> keyLineDir(3); 
    std::vector<std::vector<int> > keyLineOnPlaneID;
    bInf->cutNormals = (real8 (*)[3])malloc(3*bInf->numPlanes*sizeof(real8));
    for(auto i=0; i<bInf->numPlanes; i++){
        cross_product(bInf->planeList[i], GBNormal, crossProdt); normalize_vec(crossProdt);
        if(normal(crossProdt) < 1.0E-2){
            bInf->planeType[i] |= GB_PLANE; VECTOR_ZERO(bInf->cutNormals[i]);
        }else{
            // if the GB cut the plane, we record the info.
            cross_product(bInf->planeList[i], crossProdt, bInf->cutNormals[i]); normalize_vec(bInf->cutNormals[i]);
            dot = DotProduct(bInf->cutNormals[i], GBNormal);
            if((dot>0&&i<bInf1->numPlanes) || (dot<0&&i>=bInf1->numPlanes)){
                bInf->cutNormals[i][0] = -bInf->cutNormals[i][0];
                bInf->cutNormals[i][1] = -bInf->cutNormals[i][1];
                bInf->cutNormals[i][2] = -bInf->cutNormals[i][2];
            }
            // if the plane is glide plane, record the intersection line
            if((bInf->planeType[i] & GLIDE_PLANE)> 0){
                int j=0; for(;j<int(keyLineList.size());j++){
                    if(fabs(DotProduct(keyLineList[j],crossProdt))>0.99)break;
                }
                if(j==int(keyLineList.size())){
                    VECTOR_COPY(keyLineDir, crossProdt);
                    keyLineList.push_back(keyLineDir);
                    keyLineOnPlaneID.resize(j+1);
                }
                keyLineOnPlaneID[j].push_back(i);
            }
        }
    }

    bInf->numKeyLines = int(keyLineList.size());
    bInf->keyLineList = (real8 (*)[3])malloc(3*bInf->numKeyLines*sizeof(real8));
    bInf->numPlanesPerKeyLine = (int *)malloc(bInf->numKeyLines*sizeof(int));
    bInf->keyLineOnPlaneID = (int **)malloc(bInf->numKeyLines*sizeof(int*));
    for(auto i=0; i<int(keyLineList.size()); i++){
        VECTOR_COPY(bInf->keyLineList[i], keyLineList[i]);
        bInf->numPlanesPerKeyLine[i]=(int)keyLineOnPlaneID[i].size();
        bInf->keyLineOnPlaneID[i] = (int *)malloc(bInf->numPlanesPerKeyLine[i]*sizeof(int));
        for(auto j=0; j<bInf->numPlanesPerKeyLine[i]; j++){
            bInf->keyLineOnPlaneID[i][j]=keyLineOnPlaneID[i][j];
        }
    }
    keyLineList.clear(); keyLineOnPlaneID.clear();
    
    /** Assign burg Information and compress the burg data
     */
    std::vector<int>    burgType(bInf1->numBurgVectors+bInf2->numBurgVectors); 
    std::vector<std::vector<real8> > burgList(bInf1->numBurgVectors+bInf2->numBurgVectors);
    for(auto i=0;i<bInf1->numBurgVectors;i++){
        burgList[i].resize(3); 
        burgList[i][0]=bInf1->burgList[i][0];
        burgList[i][1]=bInf1->burgList[i][1];
        burgList[i][2]=bInf1->burgList[i][2];
        burgType[i]=bInf1->burgType[i];
    }
    for(auto i=0;i<bInf2->numBurgVectors;i++){
        int j=i+bInf1->numBurgVectors;
        burgList[j].resize(3); 
        burgList[j][0]=bInf2->burgList[i][0];
        burgList[j][1]=bInf2->burgList[i][1];
        burgList[j][2]=bInf2->burgList[i][2];
        burgType[j]=bInf2->burgType[i];
    }
    CompressBurgListWithType(burgList, burgType);

    bInf->numBurgVectors = (int)burgList.size();
    bInf->burgList = (real8 (*)[3])malloc(3*bInf->numBurgVectors*sizeof(real8));
    bInf->burgType = (int *)malloc(bInf->numBurgVectors*sizeof(int));
    bInf->numPlanesPerBurg = (int *)malloc(bInf->numBurgVectors*sizeof(int));
    std::vector<std::vector<int> >        burgOnPlaneID(bInf->numBurgVectors);

    // copy stl vectors to list count burgOnPlaneID
    real8 burgNormalized[3];
    for(auto i=0; i<bInf->numBurgVectors; i++){
        VECTOR_COPY(bInf->burgList[i], burgList[i]);
        bInf->burgType[i] = burgType[i];
        bInf->numPlanesPerBurg[i] = 0;
        VECTOR_COPY(burgNormalized, burgList[i]);
        normalize_vec(burgNormalized);
        for(auto j=0; j<bInf->numPlanes; j++){
            if(fabs(DotProduct(burgNormalized, bInf->planeList[j])) < 1.0E-2){
                bInf->numPlanesPerBurg[i]++;
                burgOnPlaneID[i].push_back(j);
            }
        }
    }
    burgList.clear(); burgType.clear();
    
    // Copy burgOnPlaneID and assembly planeOnBurgID
    std::vector<std::vector<int> >      planeOnBurgID(bInf->numPlanes);
    bInf->burgOnPlaneID = (int **)malloc(bInf->numBurgVectors*sizeof(int *));
    for(auto i=0; i<bInf->numBurgVectors; i++){
        bInf->burgOnPlaneID[i] = (int *)malloc(bInf->numPlanesPerBurg[i]*sizeof(int));
        for(auto j=0; j<bInf->numPlanesPerBurg[i]; j++){
            bInf->burgOnPlaneID[i][j] = burgOnPlaneID[i][j];
            planeOnBurgID[burgOnPlaneID[i][j]].push_back(i);
        }
    }
    burgOnPlaneID.clear();

    // Copy planeOnBurgID
    bInf->numBurgsPerPlane = (int *)malloc(bInf->numPlanes*sizeof(int));
    bInf->planeOnBurgID = (int **)malloc(bInf->numPlanes*sizeof(int *));
    for(auto i=0; i<bInf->numPlanes; i++){
        bInf->numBurgsPerPlane[i] = planeOnBurgID[i].size();
        bInf->planeOnBurgID[i] = (int *)malloc(bInf->numBurgsPerPlane[i]*sizeof(int));
        for(auto j=0; j<bInf->numBurgsPerPlane[i]; j++){
            bInf->planeOnBurgID[i][j] = planeOnBurgID[i][j];
        }
    }
    planeOnBurgID.clear();
    
    return;
}

/** Free Burgers information structure */
void FreeBurgInfoStr(BurgInfo_t *burgData)
{
    
    int     i;

    free(burgData->burgList); 
    free(burgData->planeList); 
    free(burgData->cutNormals); 
    free(burgData->planeType); 
    free(burgData->burgType); 
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
    burgData->cutNormals = NULL;
    burgData->keyLineList = NULL;

    burgData->numPlanesPerBurg = NULL;
    burgData->numPlanesPerKeyLine = NULL;
    burgData->numBurgsPerPlane = NULL;

    burgData->burgOnPlaneID = NULL;
    burgData->planeOnBurgID = NULL;
    burgData->keyLineOnPlaneID = NULL;

    return;
}

/* Return type of Burgers vector based on magnintude of the Burgers vector */
int	BurgType_FCC(real8 burg[3]){
	real8 eps = 1e-2;
	real8 burgMag;
	burgMag = normal(burg);

    //	0: zero burgers
    if(burgMag < eps) return(ZERO_BURG);

    // 	1: basic Burgers vector
	if(fabs(burgMag - 1.0) < eps) return(BASIC_BURG);
	
    //	2: Shockley burgers
	if(fabs(burgMag - SHOCKLEY_BURG_MAG_FCC) < eps) return(SHOCKLEY_BURG);
	
    //	3: junction 1 burgers, like [Sqrt(2),0,0]
	if(fabs(burgMag - JUNCTION1_BURG_MAG_FCC) < eps) return(JUNCTION1_BURG);
	
    //	3: junction 2 Burgers, like DA+DB
	if(fabs(burgMag - JUNCTION2_BURG_MAG_FCC) < eps) return(JUNCTION2_BURG);

	return(UNKNOWN_BURG);
}

/** Convert crystallographic (local) vector to labortorial (global) vector */
void VectorCrystalToLab(Home_t *home, int *place, real8 *cryst, real8 *lab)
{
    Face_t *face; Edge_t *edge; int jointID;
    switch(place[0]){
        case IN_GRAIN:
            Matrix33Vector3Multiply(home->grainKeys[place[1]]->rotMatrixs[InverseIndex(place[2])%2], cryst, lab); break;
        case ON_IB:
            Matrix33Vector3Multiply(home->grainKeys[place[1]]->rotMatrixs[place[2]%2],cryst,lab); break;
        case ON_FACE:
            face = home->faces[place[1]];
            if(place[2]<0) Matrix33Vector3Multiply(home->grainKeys[face->nbrGrain[0]]->rotMatrixs[InverseIndex(place[2])%2],cryst, lab);
            else Matrix33Vector3Multiply(home->grainKeys[face->nbrGrain[0]]->rotMatrixs[place[2]%2],cryst,lab);
            break;
        case ON_FACE_EDGE:
            edge = home->edges[place[1]];
#ifdef DEBUG
            if(place[2]>-1)Fatal("arm location can not at joint of edge");
#endif
            jointID = InverseIndex(place[2]);
            if(jointID==edge->nIBs)jointID--;
            Matrix33Vector3Multiply(home->grainKeys[edge->ibInGrain[jointID]]->rotMatrixs[edge->ibID[jointID]%2],cryst,lab);
            break;
        default: Fatal("%s: can not support arm location %d",__func__,place[0]);
    }
}

/** Convert labortorial (global) vector to crystallographic (local) vector */
void VectorLabToCrystal(Home_t *home, int *place, real8 *lab, real8 *cryst)
{
    Face_t *face; Edge_t *edge; int jointID;
    switch(place[0]){
        case IN_GRAIN:
            Matrix33Vector3Multiply(home->grainKeys[place[1]]->rotMatrixsInv[InverseIndex(place[2])%2], lab,cryst); break;
        case ON_IB:
            Matrix33Vector3Multiply(home->grainKeys[place[1]]->rotMatrixsInv[place[2]%2],lab,cryst); break;
        case ON_FACE:
            face = home->faces[place[1]];
            if(place[2]<0) Matrix33Vector3Multiply(home->grainKeys[face->nbrGrain[0]]->rotMatrixsInv[InverseIndex(place[2])%2],lab,cryst);
            else Matrix33Vector3Multiply(home->grainKeys[face->nbrGrain[0]]->rotMatrixsInv[place[2]%2],lab,cryst);
            break;
        case ON_FACE_EDGE:
            edge = home->edges[place[1]];
#ifdef DEBUG
            if(place[2]>-1)Fatal("arm location can not at joint of edge");
#endif
            jointID = InverseIndex(place[2]);
            if(jointID==edge->nIBs)jointID--;
            Matrix33Vector3Multiply(home->grainKeys[edge->ibInGrain[jointID]]->rotMatrixsInv[edge->ibID[jointID]%2],lab,cryst);
            break;
        default: Fatal("%s: can not support arm location %d",__func__,place[0]);
    }
}
//------------------------------------------------------------------------------

bool IsOnGB(Node_t *node)
{
    return(node->onGB > -1);
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


void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int& inGrain, int& onGB){
    int         i, j, nGrains, type, breakLoop = 0;
    real8       vec[3], dis, pos[3];
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

    if(normal(vec) < 1E-10){
        PrintNode(home, node);
        printf("Warning in NewPosToGB: {%f,%f,%f} too short of the vector\n",
                vec[0], vec[1], vec[2]);
        VECTOR_COPY(newPos, pos);
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(rval);
    }

    normalize_vec(vec);
    VECTOR_COPY(vec1, vec);
	GetGlideConstraintList(home, node, &NM, GC, zeroLenJunc, 0, 
                           plane.normal[0], plane.normal[1],
                           plane.normal[2]);
	ApplyConstraintsToVelocity(NM, GC, vec);

    if(normal(vec) < 1.0E-6 && node->constraint != 7){
        while(NM >= 0){
            NM--;
            VECTOR_COPY(vec, vec1);
            ApplyConstraintsToVelocity(NM, GC, vec);
            if(normal(vec) >= 1.0E-6)break;
        }
    }

    if(normal(vec) < 1.0E-6){
        VECTOR_COPY(newPos, pos);
        FoldBox(param, newPos, newPos+1, newPos+2);
        return(rval);
    }

    normalize_vec(vec);
    type = LinePlaneIntersection(vec, pos, plane.normal, plane.point, vec1, vec2);
    if((type&NO_INTERSECTION)>0){
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

	return(normal(vec));
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
	
	normalize_vec(vec);
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
	normalize_vec(vec1);
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
                            found = 1; inGrain = GB1->nbrGrain[i];
                            break;
                        }
                    }
                }
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
    GB_t        *gb= (GB_t *)NULL;
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

    burgMag = normal(burg);
    normalize_vec(burg);
    if(fabs(DotProduct(burg, param->TBNormal)) < 1.0E-2)return;

    vec[0] = node2->x - node1->x;
    vec[1] = node2->y - node1->y;
    vec[2] = node2->z - node1->z;
    ZImage(param, vec, vec+1, vec+2);

    length = normal(vec);
    if(length < 1.0E-2)return;

    segCentor[0] = node1->x + 0.5*vec[0];
    segCentor[1] = node1->y + 0.5*vec[1];
    segCentor[2] = node1->z + 0.5*vec[2];
    FoldBox(param, segCentor, segCentor+1, segCentor+2);

	dis =  ClosedGB(home, segCentor, &gb);
    if(fabs(dis) > 100 || fabs(dis) < 1.0E-3)return;

    normalize_vec(vec);
    cross_product(nor, vec, gDir);
    if(normal(gDir) < 1E-3)return;
    normalize_vec(gDir);

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
    int     type=-1, armID12, armID21;
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
        case NO_INTERSECTION|OUTSIDE_SPACE: canMerge=0; break;
        case NO_INTERSECTION|INSIDE_SPACE: canMerge=0; break;
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
            cross_product(vec1, vec2, crossProdt);
            if(normal(crossProdt) < 1.0E-2)break;

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


void SumOfValuesIsConstant(int n, int sum, std::vector<int>& valueList, int tot, 
                            std::vector<std::vector<int> >& valueLists, int maxVar)
{
    int i, j, nextSum=0;

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


bool cmp(std::vector<int > &p, std::vector<int > &q)
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
            cross_product(bd->planeList[i], plane, crossProd);
            switch(flag){
                case 0:
                    VECTOR_COPY(vector, bd->planeList[i]);
                    planeList.push_back(vector);
                    planeIDList.push_back(i);
                    nPlanes++;
                    break;
                case 1:
                    if(normal(crossProd) > 1.0E-2){
                        VECTOR_COPY(vector, bd->planeList[i]);
                        planeList.push_back(vector);
                        planeIDList.push_back(i);
                        nPlanes++;
                    }
                    break;
                case 2:
                    if(normal(crossProd) < 1.0E-2){
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
        
        resiBurgMag = normal(resiBurg);
        for(j=0; j<setLists[i].size(); j++){
            VECTOR_COPY(vec, burgList[j]);
            cross_product(resiBurg, vec, crossProd);
            if(normal(crossProd) < 1.0E-2 && 
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
                disList[k].energyFactor = normal(vec)/normal(burg);
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

void FindGlidePlane(Home_t *home, real8 burg[3], real8 lDir[3], real8 glidePlane[3], Node_t *node1, Node_t *node2)
{
    int         armID, burgID, normalID=-1, lineID;
    int         inGrain, onGB;
    int         burgType, interType;
    int         times;
    real8       normBurg[3], crossProdt[3];
    real8       vec1[3], vec2[3];
    Param_t     *param = home->param;
    BurgInfo_t  *bd = GetBurgData(home, node1, node2);

    srand(time(NULL));

    burgID = BurgID(bd, burg);
    GetArmGrain(home, node1, node2, inGrain, onGB);
    burgType = BurgType_FCC(burg);

    if(normal(lDir) < 1.0E-4){
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
                normalize_vec(glidePlane);
                return;
                break;
            default:
                Fatal("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
                break;
        }
    }

    if(normal(lDir) < 1.0E-4){
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
#if 0
                if(bd->planeInGrain[normalID] != inGrain){
                    normalID = -1;
                }else{
                    VECTOR_COPY(glidePlane, bd->planeList[normalID]);
                    break;
                }
#endif

                if(times >= 100){
                    VECTOR_COPY(glidePlane, burg);
                    normalize_vec(glidePlane);
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
            normalize_vec(glidePlane);
            printf("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
            return;
        }
    }

    normalize_vec(lDir);
    VECTOR_COPY(normBurg, burg);
    normalize_vec(normBurg);

    cross_product(normBurg, lDir, crossProdt);
    if(normal(crossProdt) < 1.0E-4){
        if(burgID >= bd->numBurgVectors){
            FormatVector(burg, "burg");
            FormatVector(lDir, "lDir");
            VECTOR_COPY(glidePlane, burg);
            normalize_vec(glidePlane);
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
#if 0
                if(bd->planeInGrain[normalID] != inGrain){
                    normalID = -1;
                }else{
                    VECTOR_COPY(glidePlane, bd->planeList[normalID]);
                    break;
                }
                if(times >= 100){
                    VECTOR_COPY(glidePlane, burg);
                    normalize_vec(glidePlane);
                    printf("Warning in %s at %d: chose a non-glide plane\n", __FILE__, __LINE__);
                    return;
                }
#endif
            }
        }
    }else{
        normalize_vec(crossProdt);
        normalID = PlaneID(bd, crossProdt);
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
    normalize_vec(vec);

    cross_product(normal, vec, dir);
    if(Normal(dir) < 1.0E-4)return;
    normalize_vec(dir);

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
    length = normal(xi);
    normalize_vec(xi);

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

    cross_product(strb, xi, f1);
    
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
