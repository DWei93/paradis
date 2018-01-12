#ifndef	_GRAIN_H
#define	_GRAIN_H

#include "Typedefs.h"
#include "Home.h"
#include "Grain.h"
#include "Util.h" 
#include "Tag.h" 
#include "ParadisProto.h"

/*
 * 	Mesh number used for output of plastic strain field.
 */
#define MESH_NUM	41
#define ASPECT_RATIO	2
#define MESH_NUM1	(MESH_NUM+1)
#define TOT_ELE_NUM	(MESH_NUM*MESH_NUM*MESH_NUM*ASPECT_RATIO)
#define TOT_NODE_NUM	(MESH_NUM1*MESH_NUM1*(ASPECT_RATIO*MESH_NUM+1))

#define GRAIN_DEBUG			1

#define TWIN_LAYER_RADIUS			1e-2 
#define	SAFE_LAYER_RADIUS			1
#define ERROR_TWIN_LAYER_RADIUS		500

#define	MATRIX_GRAIN		0
#define TWIN_GRAIN			1
#define MATRIX_BOUNDARY		0
#define TWIN_BOUNDARY		1

#define MATRIX_NODE			0
#define BOUNDARY_NODE		0.5
#define TWIN_1_NODE			1

#define DELETE_REPEAT_FLAG  0x80000000

#define FIRST_PENETRATION	0
#define SIDE_PENETRATION	1
#define LAST_PENETRATION	3

/*	define some node types	*/
#define INGRAIN_NODE	0x0000	// node within grains
#define INTERSECT_NODE	0x0001	// intersection node between TB node and node within grain
#define TB_NODE		0x0002	// node with a segment with TB burgers
#define	RESIDUAL_NODE	0x0004	// other node except TB node

/*	define some burg types	*/
#define ZERO_BURG		0x0000	// zero burgers
#define BASIC_BURG		0x0001	// basic burgers
#define TB_BURG			0x0002	// TB burgers
#define UNKNOWN_BURG	0x0004	

#define INBOX   0x0000
#define EXCEEDX 0x0001
#define EXCEEDY 0x0002
#define EXCEEDZ 0x0004
struct _grain{
	real8	cx, cy, cz;
	real8	orientation[3][3];
	real8	orientationInv[3][3];
	int		type;	/*0: Matrix*/
					/*1: Twin*/

	/*
	 *	Boundaries	
	 *	nx * (x - x0) + ny * (y - y0) + nz * (z - z0) =0
	 *	(nx, ny, nz) is outer normal orientation
	 */
	int		numBoundaries;
	int		*BoundaryIndex;
	int		*BoundaryType;	/*0: grain boundary*/
							/*1: twin boundary*/

	real8	*nx, *ny, *nz;
	real8	*x0, *y0, *z0;

};

typedef struct {
	real8 normal[3];
	real8 point[3];
} Plane_t; 


typedef struct {
	real8 direction[3];
	real8 point[3];
} Line_t;

struct _peneResults{
	int		nSets;
	real8	(*outNormal)[3];
	real8	(*outBurg)[3];
	real8	(*GBBurg)[3];
	real8	*powerDis;
	int		*effective;
	
};

void PrintNeighbor(Home_t *home, Node_t *node);
void FormatMatrix(real8 M[3][3], char *msg);
void FormatMatrix66(real8 M[6][6], char *msg);
void FormatVector(real8 vec[3], char *msg);
void CopyMatrix33(real8 M1[3][3], real8 M2[3][3]);
void InitOrientation(real8 MatrixOrientation[3][3], real8 TwinOrientation[3][3]);
void InitGrain(Home_t *home, Param_t *param);
real8 GetInGrainFromPos(Home_t *home, real8 x, real8 y, real8 z);
Node_t *GetinDataNodeFromTag(InData_t *inData, Tag_t tag);
int JudgeNodeInCuboid(real8 *x1, real8 *y1, real8 *z1, real8 x2, real8 y2, real8 z2);
int GlidePlaneID_FCC(real8 nCryst[3]);
void GetTBBurg_FCC(int nID, real8 burgTB[3]);
int	BurgType_FCC(Home_t *home, real8 *burg);
int BurgID_FCC(real8 bCryst[3]);
int GetBurgID_FCC(Home_t *home, Node_t *node, Node_t *nbr);
int HasShortSeg(Home_t *home, Node_t *node, real8 shortSegLenSq);
void AssignGrain (Home_t *home, Node_t *node, real8 inGrain);
void AssignConstraint(Home_t *home, Node_t *node, int constraint) ;
void AddOpGB (Home_t *home, OpType_t type, int dom1, int idx1, real8 inGrain);
real8 InWhichGrain(Home_t *home, real8 x, real8 y, real8 z);
int NodeOnGB(Home_t *home, Node_t *node);
real8 CorrectInGrain(Home_t *home, Node_t *node, real8 inGrainPos);
real8 ChangeOrientation(Home_t *home, Node_t *node, Node_t *nbr, int flag);
real8 ClosedGB(Home_t *home, real8 *pos, Plane_t *GBPlane);
void PullToGB(Home_t *home, Node_t *node, real8 oldInGrain);
void PullBack(Home_t *home, Node_t *node, real8 oldInGrain);
void PenetrationRule(Home_t *home, Node_t *node);
void HandleGrainBoundary(Home_t *home);
void GetGBPlane(Home_t *home,real8 oldGrain, real8 currGrain, 
				Node_t *node, Plane_t *GBPlane);

real8 FGlidePerLength(Home_t *home,Node_t *node,Plane_t *GBPlane);
real8 Distance(Home_t *home,real8 *p1,real8 *p2);
int NodeOnLine(real8 *p, Line_t *line);
void PullToLine(real8 *p, Line_t *line);
void PrintTag(Node_t *node);
void InitPeneRes(PeneRes_t *peneRes);
int GetPeneResults(Home_t *home, Node_t *node, Node_t *nbr, PeneRes_t *peneRes,
				   Grain_t *oldGrain, Grain_t *currGrain, Plane_t *GBPlane);

/*
 * geometrical relationship
 */

real8 PointPlaneDistance(real8 *point, Plane_t *plane);
int IntersectionOfLinePlane(Line_t *line, Plane_t *plane, real8 *p);
int NodeOnLine(real8 *p, Line_t *line);
int LineOnPlane(Line_t *line, Plane_t *plane);
int ClosestPointOfPointLine(real8 *p, Line_t *line, real8 *outP);
int IntersectionOfPlanePlane(Plane_t *plane1, Plane_t *plane2, Line_t *line);
real8 GetArmInGrain(real8 inGrain1, real8 inGrain2);
void ClosetPointToGB(Param_t *param, Plane_t *plane, Plane_t *GBPlane, real8 *p);

/*
 * 	others
 */

int	IsZeroVector(real8 *vec);

int  Coor2EleNum(Home_t *home,real8 position[3]);
void AdjustSurfaceNodeVel(Param_t *param, Node_t *node);
void AdjustSurfaceNode(Home_t *home);
int SameTag(Tag_t *tag1, Tag_t *tag2);

void RotateWithAxis(real8 a[3], real8 theta, real8 R[3][3]);
void DoResistance(Home_t *home, Node_t *node, real8 *nForce);
int SurfaceNode(Node_t *node);
int CheckReposition(Home_t *home, Node_t *node, real8 *newPos);
int ExceedBox(Param_t *param, real8 *pos);
int IntersectionOfLineBox(Param_t *param, Line_t *line, real8 interPos[2][3]);
#endif
