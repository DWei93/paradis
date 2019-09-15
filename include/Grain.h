#ifndef	_GRAIN_H
#define	_GRAIN_H

#include "Typedefs.h"
#include "Home.h"
#include "Grain.h"
#include "Util.h" 
#include "Tag.h" 
#include "ParadisProto.h"

#include <vector>

#ifdef VORO
#include "voro++.hh"
using namespace voro;
#endif

/*
 * 	Mesh number used for output of plastic strain field.
 */

#define FACE_ERROR			1e-11 
#define	SAFE_FACE_ERROR	    2e-11
#define IB_START_ID                 10000000
#define GRAIN_ID_FACTOR             1000

#define GLIDE_PLANE         0x0001
#define JUNCTION_PLANE1     0x0002
#define JUNCTION_PLANE2     0x0004
#define GB_PLANE            0x0008
#define UNKOWN_PLANE        0x0010

#define EDGE_BURG       0x0001
#define PARTIAL_BURG    0x0002
#define JUNCTION_BURG   0x0004
#define UNKOWN_BURG     0x0008

/** Convert inner boundary index to subgrain index or
    Convert subgrain index to inner boundary index */
#define InverseIndex(a) (-1-(a))

typedef struct{
    int     vertID1, vertID2;
    int     nFaces, *ofFaces;
    int     nIBs, *ibInGrain, *ibID;
    real8   *dis, *ibSign;
}Edge_t;

typedef struct{
    int     inGrain, ID;    
    real8   dis, sign;
}EdgeJoint_t;

typedef struct{
    int     vertID1, vertID2;
    double  factor;
    std::vector<int> ofFaces;
    std::vector<EdgeJoint_t> ibs;
}EdgeVec_t;

typedef	enum{
    ZERO_BURG,
    BASIC_BURG,
    SHOCKLEY_BURG,
    JUNCTION1_BURG,
    JUNCTION2_BURG, 
    UNKNOWN_BURG 
} BurgTypeFCC_t;

typedef enum{
    FREE_SURFACE,
    GB_FACE
} FaceType_t;

/*
 *     Define a structure to hold the burgers vectors and associated
 *     glide planes for material for which we dynamically calculate
 *     those values.
 */
struct _burgInfo {
        int     numBurgVectors;           /* Total number of burgers vectors */
                                          /* in <burgList> below             */

        int     numPlanes;                /* Total number of planes in */
                                          /* <planeList> below         */

        int     numKeyLines;              /* Total number of planes in */
                                          /* <keyLineList> below         */

        int     *numPlanesPerBurg;        /* Number of planes associated with */
                                          /* each burgers vector in the list  */

        int     *numPlanesPerKeyLine;     /* Number of planes associated with */
                                          /* each key line in the list  */

        int     *numBurgsPerPlane;        /* Number of burgers associated with */
                                          /* each plane in the list  */

        real8   (*burgList)[3];           /* Array of burgers vectors */
        int     *burgType;

        real8   (*planeList)[3];          /* Array of planes */
        real8   (*cutNormals)[3];          /* Array of cut normal of planes, available for incomplete glide plane*/
        int     *planeType;               /* Array of type of planes*/


        real8   (*keyLineList)[3];         /* Array of key line vectors */

        int     **burgOnPlaneID;

        int     **keyLineOnPlaneID;
        
        int     **planeOnBurgID;

};

typedef struct{
    int         id, type, nbrGrain[2];
    real8       n[3], p[3];
    BurgInfo_t  bInfs[2][2];
    int         nEdges;
    real8       *edges;
    int         *edgesID;
}Face_t;

typedef struct {
	real8 normal[3];
	real8 point[3];
} Plane_t; 

typedef struct {
	real8 direction[3];
	real8 point[3];
} Line_t;

typedef struct { 
    int         GBID;
    int         type;
    int         nbrGrain[2];
    BurgInfo_t  burgData;

    real8       n[3];
    real8       p[3];

	int		    nBounds;
    int         *boundgb;
    real8       **boundn, **boundp;
} GB_t;

struct _grain{
    int         id;
	real8	    pos[3], rMax, volume;
    BurgInfo_t  bInfs[2], ibBurgInfs[2];

    int         nIBs, nIBVertices, **ibPolygons;
    int         *ibVertOnEdge, *ibVertOnSubedge;
    real8       *ibVertices;
    real8       ibN[3], (*ibPs)[3], *ibAreas;
    int         *ibsNumEdges, **ibsEdgesFaceID;
    real8       **ibsEdges;

    int         nFaces;
    Face_t      **faces;
    real8       *sign;
    real8       rotMatrixs[2][3][3], rotMatrixsInv[2][3][3];

    /** Below is ready to be abandoned.
     */
	real8	    orientation[3][3];
	real8	    orientationInv[3][3];
    int         nGBs;
    GB_t        **GB;
    real8       maxSchmid;
    int         grainID;
	int		    type;	
    BurgInfo_t  burgData;

    std::vector<SlipSys_t>      slipSys;
    std::vector<Plane_t>        surface;
};
//------------------------------------------------
#define	FOLDED_GRAIN		0
#define ODINARY_GRAIN		1
#define MATRIX_BOUNDARY		0
#define TWIN_BOUNDARY		1

/*	define some node types	*/
#define INGRAIN_NODE	    0x0001	
#define INTERSECT_NODE	    0x0002	
#define GB_NODE		        0x0004	
#define	UNKNOWN_NODE	    0x0010

#define INBOX   0x0000
#define EXCEEDX 0x0001
#define EXCEEDY 0x0002
#define EXCEEDZ 0x0004

#define NUM_BURGS_IN_FCC        21
#define NUM_PLANES_IN_FCC       13
#define NUM_KEY_LINES_IN_FCC    6
#define NUM_SLIP_SYSTEMS_IN_FCC 12

#define SHOCKLEY_BURG_MAG_FCC   0.5773502691896258
#define JUNCTION1_BURG_MAG_FCC  1.414213562373095
#define JUNCTION2_BURG_MAG_FCC  1.732050807568877


struct _slipSys{
     real8  burg[3];
     real8  normal[3];
     real8  schmid;
     int    inGrain;
     int    burgID;
     int    planeID;
}; 




typedef struct {
    bool    effective;
    real8   powerDis;
    int     nDiss;          

    int     *inGrain;       
    real8   (*normal)[3];   
    real8   (*pos)[3];     
    real8   (*burg)[3];    
} Res_t;

struct _penetration{
    Node_t  *node;          // penetration node
    GB_t    *gb;

    Line_t  line1;
    Line_t  line2;
	int		nSets;          // num of posible penetration sets
    Res_t   *res;           // results list of penentration sets
};
typedef struct {
    int     planeID;
    int     burgID;
    
    real8   energyFactor;
}Dislocation_t;


/** Rotation matrix about axis in current coordinate system
 *  \param[in] a rotation axis
 *  \param[in] theta rotation angle
 *  \param[out] R rotation matrix     */
void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]);

/*Initialize geomety and lattice*/
void InitGrain(Home_t *home, Param_t *param);

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

/** Format place information*/
void FormatPlace(int *place);

/** Determine if p1 and p2 are the same place */
inline bool SamePlace(int *p1, int *p2){
    if(p1[0] == IN_GRAIN || p1[0] == ON_IB || p1[0]==ON_FACE_EDGE ||
       p2[0] == IN_GRAIN || p2[0] == ON_IB || p2[0]==ON_FACE_EDGE )return p1[0]==p2[0] && p1[1]==p2[1] && p1[2] == p2[2];
    else if(p1[0] == ON_FACE || p2[0] == ON_FACE) return p1[0]==p2[0] && p1[1]==p2[1] && p1[2] == p2[2] && p1[3] == p2[3];
    else if(p1[0] == AT_VERTICE || p1[0] == AT_VERTICE) return p1[0]==p2[0] && p1[1]==p2[1];
    else  Fatal("%s: can not support location type %d",__func__,__LINE__);
}

/** According to the current node coordinates, roughly judge the node position.
  * \param[in] pos: the coordinates of node
  * \param[out] place: the rough location of the node
  * \param[in] error: the error, absoule margin of face and edge
  * Note that the point should not be outside of the geometry */
void NodeLocation(Home_t *home, real8 *pos, int *place, real8 error=FACE_ERROR);
inline void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int *place, real8 error=FACE_ERROR);
void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int *place, real8 error)
{ real8 pos[3] = {x, y, z}; return NodeLocation(home, pos, place, error);}

/** Judge if the node is in the grain
 *  \param[in] grain the grain
 *  \param[in] pos position of the node
 *  \param[in] margin inner margin
 *  \return true if the node is in the grain, false otherwise*/
inline bool InGrain(Grain_t *grain, real8 *pos, real8 margin=FACE_ERROR){
    real8 dis;
    // make sure the node is in the grain
    for(auto i=0; i<grain->nFaces; i++){
        PointPlaneDistance(pos, grain->faces[i]->n, grain->faces[i]->p, dis);
        dis *= grain->sign[i];
        if(dis>-margin)return false;
    }
    return true;
}

/** Judge which subgrain the node is in, but there is no guarantee that the node is in the grain
 *  \param[in] grain the grain
 *  \param[in] pos position of the node
 *  \param[in] error inner boundary error
 *  \return subgrain ID (<0) or inner boundary ID (>=0)*/
inline int WhichSubgrain(Grain_t *grain, real8 *pos, real8 error=FACE_ERROR){
    real8 dis;
    for(auto i=0; i<grain->nIBs; i++){
        PointPlaneDistance(pos, grain->ibN, grain->ibPs[i], dis);
        if(fabs(dis)<error)return i;
        else if (dis<0) return InverseIndex(i);
        else continue;
    }    
    return InverseIndex(grain->nIBs);
}

/** According to the current position of the node, determine whether the position
 *  is at the location
 * \param[in] pos the coordinates of node
 * \param[in] place the location of the node
 * \param[in] error error for determination of edge and joint
 * \return false if the position are not at the location, true otherwise. */
bool InGrain(Home_t *home, int *place, real8 *pos, real8 error=FACE_ERROR);
inline bool InGrain(Home_t *home, int *place, real8 x, real8 y, real8 z, real8 error=FACE_ERROR)
{real8 pos[3]={x,y,z}; return InGrain(home, place, pos, error);}
bool OnIB(Home_t *home, int *place, real8 *pos, real8 error=FACE_ERROR);
inline bool OnIB(Home_t *home, int *place, real8 x, real8 y, real8 z, real8 error=FACE_ERROR)
{real8 pos[3]={x,y,z}; return OnIB(home, place, pos, error);}
bool OnFace(Home_t *home, int *place, real8 *pos, real8 error=FACE_ERROR);
inline bool OnFace(Home_t *home, int *place, real8 x, real8 y, real8 z, real8 error=FACE_ERROR)
{real8 pos[3]={x,y,z}; return OnFace(home, place, pos, error);}
bool OnFaceEdge(Home_t *home, int *place, real8 *pos, real8 error=FACE_ERROR);
inline bool OnFaceEdge(Home_t *home, int *place, real8 x, real8 y, real8 z, real8 error=FACE_ERROR)
{real8 pos[3]={x,y,z}; return OnFaceEdge(home, place, pos, error);}
bool AtVertice(Home_t *home, int *place, real8 *pos, real8 error=FACE_ERROR);
inline bool AtVertice(Home_t *home, int *place, real8 x, real8 y, real8 z, real8 error=FACE_ERROR)
{real8 pos[3]={x,y,z}; return AtVertice(home,place,pos,error); }
bool AtLocation(Home_t *home, int *place, real8 *pos, real8 error=FACE_ERROR);
inline bool AtLocation(Home_t *home, int *place, real8 x, real8 y, real8 z, real8 error=FACE_ERROR)
{real8 pos[3]={x,y,z}; return AtLocation(home,place,pos,error); }

/** Find the intersection between nodal path and place (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position of the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the node, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of intersection (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathInGrain(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);
bool ClipPathOnIB(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);
bool ClipPathOnFace(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);
bool ClipPathOnFaceEdge(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);
bool ClipPathAtVertice(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);
bool ClipPathAtLocation(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);

/** Generate Burgers data information
 *  \param[in] matType material type
 *  \param[in] rotMatrix rotation matrix 
 *  \param[out] burgData Burgers data */
void GenerateBurgInfo(int matType, real8 rotMatrix[3][3], BurgInfo_t *burgData);

/** Assemble GB Burgers information by bonding Burgers information of adjacent grains at grain boundary
 *  \param[out] bInf Burgers informatiob of the GB
 *  \param[in] bInf1 bInf2 two Burgers information of adjacent grains
 *  \param[in] GBNormal normal of grain boundary for cacualte cutting plane */
void GenerateGBBurgInfo(BurgInfo_t *bInf, BurgInfo_t *bInf1, real8 *GBNormal, BurgInfo_t *bInf2);

/** Free Burgers information structure */
void FreeBurgInfoStr(BurgInfo_t *burgData);

/** Get Arm location 
 *  \param[in] place1, place2 two locations of the arm ends
 *  \param[out] place arm location*/
void ArmLocation(Home_t *home, int *place1, int *place2, int *place);

/*------------------------------------------Burgers Information Functions----------------------------*/
/** Get Burgers information of arm
 * \param[in] place arm location
 * \return burg information pointer */
BurgInfo_t *BurgInformation(Home_t *home, int *place);
inline BurgInfo_t *BurgInformation(Home_t *home, int *place1, int *place2){
    int place[PLACE_LEN]; ArmLocation(home, place1, place2, place);
    return BurgInformation(home, place);
}

/** Get Burgers index from Burgers information */
inline int BurgID(BurgInfo_t *bInf, real8 burg[3]){
    for(int i=0; i<bInf->numBurgVectors; i++){
        if(AlignedVector(bInf->burgList[i], burg))return i;
    }
    return bInf->numBurgVectors;
}
inline int BurgID(Home_t *home, Node_t *node, Node_t *nbr){
    real8       burg[3];
    BurgInfo_t  *burgData = BurgInformation(home, node->place, nbr->place);
    int armID = GetArmID(home, node, nbr);
#ifdef DEBUG
    if(armID == -1 )Fatal("%s: wrong connection", __func__);
#endif
    burg[0] = node->burgX[armID];
    burg[1] = node->burgY[armID];
    burg[2] = node->burgZ[armID];
    return(BurgID(burgData, burg)); 
}

/** Match Burgers vector that have aligned direction 
 *  \return return matched Burgers vector index*/
int MatchBurgID(BurgInfo_t *bInf, real8 burg[3]);
inline int MatchBurgID(Home_t *home, Node_t *node, Node_t *nbr){
    real8       burg[3];
    BurgInfo_t  *burgData = BurgInformation(home, node->place, nbr->place);
    int armID = GetArmID(home, node, nbr);
#ifdef DEBUG
    if(armID == -1 )Fatal("%s: wrong connection", __func__);
#endif
    burg[0] = node->burgX[armID];
    burg[1] = node->burgY[armID];
    burg[2] = node->burgZ[armID];
    return(MatchBurgID(burgData, burg)); 
}

/** Get Burgers type from Burgers information */
inline int BurgType(BurgInfo_t *burgData, real8 burg[3]){
    int burgID = BurgID(burgData, burg);
    if(burgID >= burgData->numBurgVectors) return(UNKOWN_BURG);
    return(burgData->burgType[burgID]);
}
inline int BurgType(BurgInfo_t *burgData, real8 bx, real8 by, real8 bz)
{ real8 b[3]={bx, by, bz}; return BurgType(burgData, b);}

inline int BurgType(Home_t *home, Node_t *node, Node_t *nbr){
    real8       burg[3];
    BurgInfo_t *burgData = BurgInformation(home, node->place, nbr->place);
    int armID = GetArmID(home, node, nbr);
#ifdef DEBUG
    if(armID == -1 )Fatal("%s: wrong connection", __func__);
#endif
    burg[0] = node->burgX[armID];
    burg[1] = node->burgY[armID];
    burg[2] = node->burgZ[armID];
    return(BurgType(burgData, burg)); 
}

/** Get plane index from Burgers information */
inline int PlaneID(BurgInfo_t *burgData, real8 plane[3]){
    for(int i=0; i<burgData->numPlanes; i++){
        if(AlignedDirection(burgData->planeList[i], plane))return i;
    }
   return burgData->numPlanes;
}
inline int PlaneID(Home_t *home, Node_t *node, Node_t *nbr){
    real8       plane[3];
    BurgInfo_t  *burgData = BurgInformation(home, node->place, nbr->place);
    int armID = GetArmID(home, node, nbr);
#ifdef DEBUG
    if(armID == -1 )Fatal("%s: wrong connection", __func__);
#endif
    plane[0] = node->nx[armID];
    plane[1] = node->ny[armID];
    plane[2] = node->nz[armID];
    return(PlaneID(burgData, plane)); 
}

/** Get plane type from Burgers information */
inline int PlaneType(BurgInfo_t *burgData, real8 plane[3]){
    int planeID = PlaneID(burgData, plane);
    if(planeID >= burgData->numPlanes) return(UNKOWN_PLANE);
    return(burgData->planeType[planeID]);
}
inline int PlaneType(BurgInfo_t *burgData, real8 nx, real8 ny, real8 nz)
{ real8 n[3]={nx, ny, nz}; return PlaneType(burgData, n);}

inline int PlaneType(Home_t *home, Node_t *node, Node_t *nbr){
    real8       plane[3];
    BurgInfo_t *burgData = BurgInformation(home, node->place, nbr->place);
    int armID = GetArmID(home, node, nbr);
#ifdef DEBUG
    if(armID == -1 )Fatal("%s: wrong connection", __func__);
#endif
    plane[0] = node->nx[armID];
    plane[1] = node->ny[armID];
    plane[2] = node->nz[armID];
    return(PlaneType(burgData, plane)); 
}

inline int KeyLineID_FCC(BurgInfo_t *burgData, real8 line[3]){
    for(int i=0; i<burgData->numKeyLines; i++){
        if(AlignedDirection(burgData->keyLineList[i], line))return i;
    }
    return burgData->numKeyLines;
}

/* Return type of Burgers vector based on magnintude of the Burgers vector */
int	BurgType_FCC(real8 burg[3]);

/** Convert crystallographic (local) vector to labortorial (global) vector */
void VectorCrystalToLab(Home_t *home, int *place, real8 *cryst, real8 *lab);
inline void VectorCrystalToLab(Home_t *home, Node_t *node1, Node_t *node2, real8 *cryst, real8 *lab){
    int place[PLACE_LEN]; ArmLocation(home, node1->place, node2->place, place);
    return VectorCrystalToLab(home, place, cryst, lab);
}

/** Convert labortorial (global) vector to crystallographic (local) vector */
void VectorLabToCrystal(Home_t *home, int *place, real8 *lab, real8 *cryst);
inline void VectorLabToCrystal(Home_t *home, Node_t *node1, Node_t *node2, real8 *lab, real8 *cryst){
    int place[PLACE_LEN]; ArmLocation(home, node1->place, node2->place, place);
    return VectorLabToCrystal(home, place, lab, cryst);
}
/*-------------------------------------End of Burgers Information Functions-------------------------*/


/** get glide constraint of the node
 *  \param[in] node the node to be caculated
 *  \param[out] vec constraint normal (plane constraint) or direction (line direction)
 *  \return constrain type 0: no constraint; 1: plane constraint; 2: line constraint; 3: fixed*/
int GlideConstraint(Home_t *home, Node_t *node, real8 *vec);

/** Find plane of the arm, if there are several glide planes, we chose the glide plane with 
 *  maximal Schmid factor if FEM is disabled, or chose the glide plane with highest local dissipation
 *  if FEM is abled.
 *  \param[in] burg Burgers vector of the arm
 *  \param[in,out] lDir line direction of the arm. Note: no climb component in the line direction
 *  \param[in] node1, node the ends node of the arm
 *  \param[out] plane plane of the arm 
 *  \param[out] found true is found, otherwise false*/
void FindGlidePlane(Home_t *home, real8 *burg, real8 *lDir, Node_t *node1, Node_t *node2, real8 *plane, bool &found);

/** Get the vector of nodal movement under the constraint type 
 *  \param[in] node the node pointer to be caculated 
 *  \param[in,out] velocity the velocity to be corrected along the constraint
 *  \return constraint type   */
inline int ApplyConstraintsToVelocity(Home_t *home, Node_t *node, real8 *velocity){
    real8 vec[3];
    int type=GlideConstraint(home, node, vec);
    if(type == 1)NormalComponent(velocity, vec, velocity);
    else if(type == 2)ParallelComponent(velocity,vec, velocity);
    else if(type == 3){VECTOR_ZERO(velocity);}
    else Fatal("%s: wrong constraint type %d", __func__, type);
    return type;
}

/** Get edges of intersection between plane and grain
 *  \param[in] n the normal vector of plane
 *  \param[in] p a point on the plane, this point has to be in the grain.
 *  \param[in] place the grin location
 *  \param[out] vertice coordinates of vertices
 *  \param[out] edges edge list     */
void PlaneGrainIntersection(Home_t *home, real8 *n, real8 *p, int *place, std::vector<real8> &vertices, std::vector<EdgeVec_t> &edges);

/** Draw polygon with tecplot formart
 *  \param[in] fp file pointer for output
 *  \param[in] f_vert vertices index list
 *  \param[in] v coordinates of vertices
 *  \param[in] startID the start index in f_vert for output polygon
 *  \param[in] zoneName the zone name of tecplot
 *  \param[in] vec the variales of the polygon */
void DrawPolygon(FILE *fp, std::vector<int> &f_vert,std::vector<double> &v,int startID, char *zoneName, std::vector<int> &vec);
/** Draw line between two points*/
void DrawLine(FILE *fp, real8 *pos1, real8 *pos2, real8 meshSize, std::vector<int> &indexs);
/** Convert edge vector list to polygon vector list */
void EdgeListToPolygon(std::vector<EdgeVec_t> &edges, std::vector<int> &polygon);

/** Get the distance between point and polygon
 *  \param[in] p the coordinate of the point
 *  \param[in] v coordinates lis of vertices
 *  \param[in edges indexs and factor of polygon vertices
 *  \return distance                    */
double PointConcavePolygonDistance(double *p, std::vector<double> &v, std::vector<EdgeVec_t> &edges);

//-----------------------------------------------------------------------------------------------
Node_t *InsertNode(Home_t *home, Node_t *node1, Node_t *node2, real8 *p, int globalOp, int opClass);
void InitOrientation(real8 MatrixOrientation[3][3], real8 TwinOrientation[3][3]);
void FreeGrainStructure(Grain_t *grain);
void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int& inGrain, int& onGB);
Node_t *GetinDataNodeFromTag(InData_t *inData, Tag_t tag);
BurgInfo_t* GetBurgData(Home_t *home, Node_t *node1, Node_t *node2);
BurgInfo_t* GetBurgData(Home_t *home, int inGrain, int onGB);
int HasShortSeg(Home_t *home, Node_t *node, real8 shortSegLenSq);
void AssignGrain (Home_t *home, Node_t *node, int inGrain, int onGB);
void AssignConstraint(Home_t *home, Node_t *node, int constraint) ;
int GetGBIndexOfGrain(Grain_t *grain, int onGB);
bool GetPublicGB(Grain_t *grain1, Grain_t *grain2, GB_t **GB);
void AddOpGB (Home_t *home, OpType_t type, int dom1, int idx1, int inGrain, int onGB);
int NodeOnGB(Home_t *home, Node_t *node);
real8 ClosedGB(Home_t *home, real8 *pos, GB_t **GB);
bool NewPosToGB(Home_t *home, Node_t *node, int inGrain, int onGB, real8 pos[3], real8 *newPos);
void PullBack(Home_t *home, Node_t *node);
void PenetrationRule(Home_t *home, Node_t *node);
void HandleGrainBoundary(Home_t *home);

real8 FGlidePerLength(Home_t *home,Node_t *node,Plane_t *GBPlane);
real8 Distance(Home_t *home,real8 *p1,real8 *p2);
int NodeOnLine(real8 *p, Line_t *line);
void PrintTag(Node_t *node);
void FreePenStr(Pene_t *pene);
int  GetPenResults(Home_t *home, Node_t *node, Pene_t *pene, GB_t *gb);

/*
 * geometrical relationship
 */

real8 PointPlaneDistance(real8 *point, Plane_t *plane);
int IntersectionOfLinePlane(Line_t *line, Plane_t *plane, real8 *p);
int NodeOnLine(real8 *p, Line_t *line);
int LineOnPlane(Line_t *line, Plane_t *plane);
int ClosestPointOfPointLine(real8 *p, Line_t *line, real8 *outP);
int IntersectionOfPlanePlane(Plane_t *plane1, Plane_t *plane2, Line_t *line);
void GetArmGrain(Home_t *home, Node_t *node1, Node_t *node2, int &inGrain, int &onGB);

/*
 * 	others
 */

int	IsZeroVector(real8 *vec);

int  Coor2EleNum(Home_t *home,real8 position[3]);
void AdjustSurfaceNode(Home_t *home);
int SameTag(Tag_t *tag1, Tag_t *tag2);

void DoResistance(Home_t *home, Node_t *node, real8 *nForce);
int SurfaceNode(Node_t *node);
int AdjustNodeDirection(Home_t *home, Node_t *node, real8 *vec);
int ExceedBox(Param_t *param, real8 *pos);
int IntersectionOfLineBox(Param_t *param, Line_t *line, real8 interPos[2][3]);
void ImageForce(Home_t *home, Node_t *node1, Node_t *node2, real8 f1[3], real8 f2[3]);
bool BeforeMergeNodes(Home_t *home, Node_t *node1, Node_t *node2, real8 *pos, real8* dis, int &type);
bool IsOnGB(Node_t *node);
int GeometricCriterion(BurgInfo_t *bd, real8 *burg, real8 *lDir,
                       real8 *plane, int flag, 
                       std::vector< std::vector<double> >&       burgList, 
                       std::vector<std::vector<Dislocation_t> >& disArray,
                       int maxDissPerBasicBurg = 1, int maxDiss = 4);
void SplitIntersectNode(Home_t *home, Node_t *node, int opClass, int globalOp);
void GBCrossSlip(Home_t *home);
void FindGlidePlane(Home_t *home, real8 burg[3], real8 lDir[3],
                    real8 glidePlane[3], Node_t *node1, Node_t *node2);
void StackingFaultForce(Home_t *home, real8 dx, real8 dy, real8 dz, real8 sfEnergy, 
                        real8 signSF, real8 *normal, real8 *f1, real8 *f2);

void GeneratePotentialParamters(Home_t *home, Param_t *param); 
void RandomInternalForce(Home_t *home, real8 x1, real8 y1, real8 z1, real8 x2, real8 y2, real8 z2,
                         real8 dx, real8 dy, real8 dz, real8 bx, real8 by, real8 bz,
                         real8 f1[3], real8 f2[3]); 
bool AdjustSplitedNodePosition(Home_t *home, Node_t *node, real8 *currPos, real8 *pos);
bool InCuboid(Param_t *param, real8 pos[3], real8 eps);
void VoroTest();
void GeosTest();
#endif
