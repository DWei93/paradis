#ifndef	_GRAIN_H
#define	_GRAIN_H

#include "Typedefs.h"
#include "Home.h"
#include "Grain.h"
#include "Util.h" 
#include "Tag.h" 
#include "ParadisProto.h"

#include <vector>

typedef struct {
    int     planeID;
    int     burgID;
    
    real8   energyFactor;
}Dislocation_t;
/*
 * 	Mesh number used for output of plastic strain field.
 */

#define GRAIN_DEBUG			1

#define TWIN_LAYER_RADIUS			1e-2 
#define	SAFE_LAYER_RADIUS			1
#define ERROR_TWIN_LAYER_RADIUS		500

#define	FOLDED_GRAIN		0
#define ODINARY_GRAIN			1
#define MATRIX_BOUNDARY		0
#define TWIN_BOUNDARY		1

#define MATRIX_NODE			0
#define BOUNDARY_NODE		0.5
#define TWIN_1_NODE			1

#define DELETE_REPEAT_FLAG  0x80000000

#define FIRST_PENETRATION	0
#define SIDE_PENETRATION	1
#define LAST_PENETRATION	3

#define MAX_PENE_SETS               3
#define MAX_PENE_DISS               3

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

#define GLIDE_PLANE         0x0001
#define JUNCTION_PLANE1     0x0002
#define JUNCTION_PLANE2     0x0004
#define GB_PLANE            0x0008
#define UNKOWN_PLANE        0x0010

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
        int     *planeType;               /* Array of type of planes*/

        int     *planeInGrain;

        real8   (*keyLineList)[3];         /* Array of key line vectors */

        int     **burgOnPlaneID;

        int     **keyLineOnPlaneID;
        
        int     **planeOnBurgID;

};

typedef struct {
    int         GBID;
    int         type;
    int         nbrGrain[2];
    BurgInfo_t  burgData;

    real8       n[3];
    real8       p[3];
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
	real8	    cx, cy, cz;
	real8	    orientation[3][3];
	real8	    orientationInv[3][3];
    real8       maxSchmid;
    int         grainID;
	int		    type;	/*0: folded grain*/
					    /*1: ordinary grain*/
    BurgInfo_t  burgData;

	int		    nBounds;
    GB_t        **GB;
    real8       *sign;

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

Node_t *InsertNode(Home_t *home, Node_t *node1, Node_t *node2, real8 *p, int globalOp, int opClass);
void PrintNeighbor(Home_t *home, Node_t *node);
void FormatMatrix(real8 M[3][3], const char *msg);
void FormatMatrix66(real8 M[6][6], const char *msg);
void FormatVector(real8 vec[3], const char *msg);
void InitOrientation(real8 MatrixOrientation[3][3], real8 TwinOrientation[3][3]);
void FreeGrainStructure(Grain_t *grain);
void InitGrain(Home_t *home, Param_t *param);
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
void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3], int grainID = -1);
void GenerateGBBurgInfo(Home_t *home, GB_t *GB);
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
void GeosTest();
#endif
