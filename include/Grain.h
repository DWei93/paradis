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

typedef struct{
    int type;   // 0: Grain boundary edge
                // 1: Inner boubdary edge
    int vertID1, vertID2;
    std::vector<int> ofFaces;
    int ofIB[2];
}Edge_t;

/*
 * 	Mesh number used for output of plastic strain field.
 */

#define FACE_ERROR			1e-3 
#define	SAFE_FACE_ERROR	    2e-3
#define IB_START_ID                 10000000
#define GRAIN_ID_FACTOR             1000

#define GLIDE_PLANE         0x0001
#define JUNCTION_PLANE1     0x0002
#define JUNCTION_PLANE2     0x0004
#define GB_PLANE            0x0008
#define UNKOWN_PLANE        0x0010

/** Convert inner boundary index to subgrain index or
    Convert subgrain index to inner boundary index */
#define InverseIndex(a) (-1-(a))

//------------------------------------------------
#define	FOLDED_GRAIN		0
#define ODINARY_GRAIN		1
#define MATRIX_BOUNDARY		0
#define TWIN_BOUNDARY		1

/*	define some node types	*/
#define INGRAIN_NODE	    0x0001	
#define INTERSECT_NODE	    0x0002	
#define GB_NODE		        0x0004	
#define	JUNCTION_NODE	    0x0008
#define	UNKNOWN_NODE	    0x0010

#define INBOX   0x0000
#define EXCEEDX 0x0001
#define EXCEEDY 0x0002
#define EXCEEDZ 0x0004

#define NUM_BURGS_IN_FCC        21
#define NUM_PLANES_IN_FCC       13
#define NUM_KEY_LINES_IN_FCC    6


#define SHOCKLEY_BURG_MAG_FCC   0.5773502691896258
#define JUNCTION1_BURG_MAG_FCC  1.414213562373095
#define JUNCTION2_BURG_MAG_FCC  1.732050807568877

#define EDGE_BURG       0x0001
#define PARTIAL_BURG    0x0002
#define JUNCTION_BURG   0x0004

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


typedef struct {
	real8 normal[3];
	real8 point[3];
} Plane_t; 

struct _slipSys{
     real8  burg[3];
     real8  normal[3];
     real8  schmid;
     int    inGrain;
     int    burgID;
     int    planeID;
}; 

struct _grain{
    int         id;

	real8	    pos[3];

    BurgInfo_t  bInfs[2], ibBurgInfs[2];

    int         nIBs;
    real8       ibN[3], (*ibPs)[3];
    int         *ibsNumEdges, *ibsID;
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

typedef struct {
	real8 direction[3];
	real8 point[3];
} Line_t;


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

/*Initialize geomety and lattice*/
void InitGrain(Home_t *home, Param_t *param);


/** According to the current node coordinates, roughly judge the node position.
  * \param[in] pos: the coordinates of node
  * \param[out] place: the rough location of the node
  * \param[in] error: the error, absoule margin of face and edge
  * Note that the point should not be outside of the geometry */
void NodeLocation(Home_t *home, real8 *pos, int *place, real8 error=FACE_ERROR);
inline void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int *place, real8 error=FACE_ERROR);
void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int *place, real8 error)
{ real8 pos[3] = {x, y, z}; return NodeLocation(home, pos, place, error);}

/** According to the current node coordinates, judge if the node is the grain and subgrain
 * \param[in] pos: the coordinates of node
 * \param[in] place: the location of the node
 * \return false if the original position are not at the location, true otherwise. */
bool InGrain(Home_t *home, int *place, real8 *pos);
inline bool InGrain(Home_t *home, int *place, real8 x, real8 y, real8 z)
{   real8 pos[3]={x,y,z};
    return InGrain(home, place, pos);
}

/** Find the intersection between segment and grain (if has) and return the new place of intersection point.
 * \param[in] place the location of the segment.
 * \param[in] pos1 the position of first end point of the segnent, it has to be in the place.
 * \param[in] (x1,y1,z1) the position of first end point of the segnent, it has to be in the place.
 * \param[in,out] pos2 the position of second end point, moved to intersection if it was outside of the place. 
 * \param[in] (x2,y2,z2) the position of second end point of the segnent
 * \param[out] newPlace the new place of intersection (if has) 
 * return true if the intersection is found, otherwise false */
bool ClipSegmentInGrain(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);
inline bool ClipSegmentInGrain(Home_t *home, int *place, real8 x1, real8 y1, real8 z1, real8 x2, real8 y2, real8 z2, int *newPlace)
{
    real8 pos1[3]={x1, y1, z1}; real8 pos2[3]={x2, y2, z2};
    return ClipSegmentInGrain(home, place, pos1, pos2, newPlace);
}

bool CorrectedLocation(Home_t *home, int *place, real8 *pos);
inline bool IBFaceID(const int &grainID, const int &ibID, int &faceID);
inline bool IBGrainID(const int faceID, int &grainID, int &ibID);
void GetArmLocation(Home_t *home, Node_t *node1, Node_t *node2, int *place);
void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3]);
void GenerateGBBurgInfo(BurgInfo_t *bInf, BurgInfo_t *bInf1, real8 *GBNormal, BurgInfo_t *bInf2);

//-----------------------------------------------------------------------------------------------
Node_t *InsertNode(Home_t *home, Node_t *node1, Node_t *node2, real8 *p, int globalOp, int opClass);
void InitOrientation(real8 MatrixOrientation[3][3], real8 TwinOrientation[3][3]);
void FreeGrainStructure(Grain_t *grain);
void NodeLocation(Home_t *home, real8 x, real8 y, real8 z, int& inGrain, int& onGB);
Node_t *GetinDataNodeFromTag(InData_t *inData, Tag_t tag);
BurgInfo_t* GetBurgData(Home_t *home, Node_t *node1, Node_t *node2);
BurgInfo_t* GetBurgData(Home_t *home, int inGrain, int onGB);
int BurgID_FCC(BurgInfo_t *burgData, real8 burg[3]);
int	BurgType_FCC(real8 burg[3]);
int PlaneID_FCC(BurgInfo_t *burgData, real8 plane[3]);
int PlaneType_FCC(BurgInfo_t *burgData, int planeID);
int GetBurgID_FCC(Home_t *home, Node_t *node, Node_t *nbr);
int GetPlaneID_FCC(Home_t *home, Node_t *node, Node_t *nbr);
int GetPlaneType_FCC(Home_t *home, Node_t *node, Node_t *nbr);
int HasShortSeg(Home_t *home, Node_t *node, real8 shortSegLenSq);
void AssignGrain (Home_t *home, Node_t *node, int inGrain, int onGB);
void AssignConstraint(Home_t *home, Node_t *node, int constraint) ;
int GetGBIndexOfGrain(Grain_t *grain, int onGB);
bool GetPublicGB(Grain_t *grain1, Grain_t *grain2, GB_t **GB);
void AddOpGB (Home_t *home, OpType_t type, int dom1, int idx1, int inGrain, int onGB);
int NodeOnGB(Home_t *home, Node_t *node);
void VectorLabToCrystal(Home_t *home, Node_t *node1, Node_t *node2, 
                        real8 lab[3], real8 *cryst);
void VectorCrystalToLab(Home_t *home, Node_t *node1, Node_t *node2, 
                        real8 cryst[3], real8 *lab);
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

void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]);
void DoResistance(Home_t *home, Node_t *node, real8 *nForce);
int SurfaceNode(Node_t *node);
int AdjustNodeDirection(Home_t *home, Node_t *node, real8 *vec);
int ExceedBox(Param_t *param, real8 *pos);
int IntersectionOfLineBox(Param_t *param, Line_t *line, real8 interPos[2][3]);
void ImageForce(Home_t *home, Node_t *node1, Node_t *node2, real8 f1[3], real8 f2[3]);
bool BeforeMergeNodes(Home_t *home, Node_t *node1, Node_t *node2, real8 *pos, real8* dis, int &type);
bool IsOnGB(Node_t *node);
void FreeBurgInfoStr(BurgInfo_t *burgData);
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
