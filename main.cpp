#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>

#define VECTOR_ZERO(a) {a[0]=0; a[1]=0; a[2]=0;}
#define VECTOR_COPY(a,b) {(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];}
#define DotProduct(vec1,vec2)       \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))

#define real8   double

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

#define EDGE_BURG       0
#define PARTIAL_BURG    1
#define JUNCTION_BURG   2

using namespace std;
    
typedef struct {
    int     planeID;
    int     burgID;
    
    real8   energyFactor;
}Dislocation_t;


typedef	enum{
    ZERO_BURG,
    BASIC_BURG,
    SHOCKLEY_BURG,
    JUNCTION1_BURG,
    JUNCTION2_BURG, 
    UNKNOWN_BURG 
} BurgTypeFCC_t;

typedef struct {
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

} BurgInfo_t;


typedef struct {
    int         GBID;
    int         type;
    int         nbrGrain[2];
    BurgInfo_t  burgData;

    real8       n[3];
    real8       p[3];
} GB_t;

typedef struct  {
	real8	    cx, cy, cz;
	real8	    orientation[3][3];
	real8	    orientationInv[3][3];
    int         grainID;
	int		    type;	/*0: folded grain*/
					    /*1: ordinary grain*/
    BurgInfo_t  burgData;

	int		    nBounds;
    GB_t        **GB;
    real8       *sign;

} Grain_t;

real8 Normal(std::vector<real8>& a)
{
    return( sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) );
}
typedef struct {
    GB_t    **GBKeys;
    Grain_t **grainKeys;    

} Home_t; 

void cross(real8 a[3], real8 b[3], real8 c[3])
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
real8 Normal(real8 a[3])
{
        return( sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) );
}

int	BurgType_FCC(real8 burg[3]);
void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3], int grainID = -1);
void GenerateGBBurgInfo(Home_t *home, GB_t *GB);
void FreeBurgInfoStr(BurgInfo_t *burgData);
int GeometricCriterion(BurgInfo_t *bd, real8 *burg, real8 *lDir,
                       real8 *plane, int flag, 
                       vector< vector<double> >&  burgList, 
                       vector<vector<Dislocation_t> >&    disArray,
                       int maxDissPerBasicBurg = 1, int maxDiss = 4);

void InitGrain(Home_t *home);
int BurgID_FCC(BurgInfo_t *burgData, real8 burg[3]);

void Matrix33Vector3Multiply(real8 A[3][3], real8 x[3], real8 y[3])
{
        y[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
        y[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
        y[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];

        return;
}

void NormalizeVec(real8 vec[3])
{
        real8 a2, a;

        a2 = (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        if (a2 > 0.0) {
            a = sqrt(a2);
            vec[0] /= a;
            vec[1] /= a;
            vec[2] /= a;
        }

        return;
}

void InitGrainStructure(Grain_t *grain, int numBounds){
	grain->nBounds = numBounds;
    grain->GB = (GB_t **)malloc(numBounds*sizeof(GB_t *));	
    grain->sign = (real8 *)malloc(numBounds*sizeof(real8));
    return;
}

void InitHome(Home_t *home){
    int     i;
    
    home->GBKeys = (GB_t **)malloc( 1*sizeof(GB_t *));
    for(i=0; i<1; i++){
        home->GBKeys[i] = (GB_t *)calloc(1, sizeof(GB_t));
    }
    home->grainKeys = (Grain_t **)malloc(2*sizeof(Grain_t *));

    for(i=0; i<2; i++){
        home->grainKeys[i] = (Grain_t *)calloc(1, sizeof(Grain_t));
    }
    real8       I[3][3] = {
        {1.0e+00, 0.0e+00, 0.0e+00},
        {0.0e+00, 1.0e+00, 0.0e+00},
        {0.0e+00, 0.0e+00, 1.0e+00}
    };

    InitGrain(home);
    return;
}

void InitGrain(Home_t *home){
    int     i, j, k;
    GB_t    *gb;
    Grain_t *grain;

    i = 0;
    gb = home->GBKeys[0];
    gb->GBID = i;
    gb->nbrGrain[0] = 0;
    gb->nbrGrain[1] = 1;
    VECTOR_ZERO(gb->n);
    VECTOR_ZERO(gb->p);
    gb->n[1] = 1.0;

    real8 rotMatrix[3][3] = 

#if 0
    {{-0.788675, 0.57735, 
  0.211325}, {-0.57735, -0.57735, -0.57735}, {-0.211325, -0.57735, 
  0.788675}};
#else
    {
        {7.071067811865475e-01,-7.071067811865475e-01,0.000000000000000e+00},
        {5.773502691896258e-01,5.773502691896258e-01,5.773502691896258e-01},
        {-4.082482904638631e-01,-4.082482904638631e-01,8.164965809277261e-01}
    };
#endif
    real8 rotTwin[3][3] = 
#if 1        
        {{-0.707107, 0.707107, 0.}, {0.57735, 0.57735, 0.57735}, {0.408248, 
  0.408248, -0.816497}};
#else
    {{0., 0.707107, -0.707107}, {-0.57735, -0.57735, -0.57735}, 
{-0.816497, 0.408248, 0.408248}};
#endif
    for(i=0; i<2; i++){
        grain = home->grainKeys[i];
        grain->grainID = i;
        InitGrainStructure(grain, 1);
        grain->GB[0] = gb;
        if(i==0){
            memcpy(grain->orientation, rotMatrix, sizeof(real8 [3][3]));
        }else{
            memcpy(grain->orientation, rotTwin, sizeof(real8 [3][3]));
        }

        GenerateFCCBurgInfo(&(grain->burgData), grain->orientation, grain->grainID);
    }
    
    GenerateGBBurgInfo(home, gb); 

#if 1
    bool     stat;
    int      crossSlip = 1;
   real8   burgCry[3] = {0.707107, 0.0, 0.707107};    // AD
 //real8   burgCry[3] = {0.707107, -0.707107, 0.};     //AB
  //  real8   burgCry[3] = {0.3, -0.7, 0.6};     //randon
    real8   lineCry[3] = {0.707107, -0.707107, 0.};     //AB
  //  real8   planeCry[3]= {0.57735, 0.57735, 0.57735}; // ABC 
    real8   planeCry[3]= {0.57735, 0.57735, -0.57735}; // ABD 
 //   real8   planeCry[3]= {0.57735, 0.57735, 0.57735};   // ABC
    int     grainID = 0;
    real8   burg[3], line[3], plane[3], vec[3], crossProd[3];

    burgCry[0] *= 2.0;
    burgCry[1] *= 2.0;
    burgCry[2] *= 2.0;

    vector< vector<double> >                    burgList;
    std::vector<std::vector<Dislocation_t> >    disArray;
    std::vector< std::vector<Dislocation_t> >::iterator     itDisList;
    std::vector<Dislocation_t>::iterator                    itDis;
    BurgInfo_t  *gbd;

    gbd = &(gb->burgData);

    grain = home->grainKeys[grainID];
    Matrix33Vector3Multiply(grain->orientation, burgCry, burg);
    Matrix33Vector3Multiply(grain->orientation, lineCry, line);
    Matrix33Vector3Multiply(grain->orientation, planeCry,plane);
    stat = GeometricCriterion(&(gb->burgData), burg, line, plane, crossSlip,
                           burgList, disArray);

    printf("Testing: cross-slip %d\nA test Burg in grain %d:\n", crossSlip, grainID);
    printf("crystal burg {%7.4f,%7.4f,%7.4f} to global burg {%7.4f,%7.4f,%7.4f}\n",
           burgCry[0],burgCry[1],burgCry[2], burg[0],burg[1],burg[2]);
    printf("crystal line {%7.4f,%7.4f,%7.4f} to global line {%7.4f,%7.4f,%7.4f}\n",
           lineCry[0],lineCry[1],lineCry[2], line[0],line[1],line[2]);
    printf("crystal plane {%7.4f,%7.4f,%7.4f} to global plane {%7.4f,%7.4f,%7.4f}\n",
           planeCry[0],planeCry[1],planeCry[2], plane[0],plane[1],plane[2]);


    if(stat){
        printf("\n Basic shockley burg\n");
        for(i=0; i<burgList.size(); i++){
            printf("%d: {%7.4f,%7.4f,%7.4f}\n",i, 
                    burgList[i][0], burgList[i][1],burgList[i][2]);
        }


        printf("\nDislocation Arrray:\n");
        for(i=0; i<disArray.size(); i++){
            printf("%d: ",i);
            for(j=0; j<disArray[i].size(); j++){
#if 1
                printf("(%d<%d>)[%d] ", disArray[i][j].planeID, 
                        gb->burgData.planeInGrain[disArray[i][j].planeID], 
                        disArray[i][j].burgID);
#else
                printf("(%5.3f,%5.3f,%5.3f)[%5.3f,%5.3f,%5.3f] ",
                       gb->burgData.planeList[disArray[i][j].planeID][0],
                       gb->burgData.planeList[disArray[i][j].planeID][1],
                       gb->burgData.planeList[disArray[i][j].planeID][2],
                       burgList[disArray[i][j].burgID][0],
                       burgList[disArray[i][j].burgID][1],
                       burgList[disArray[i][j].burgID][2]);
#endif
            }
            printf("\n");
        }

        real8   barrierFactor;
        for(j=0; j<disArray.size(); ){
            barrierFactor = -1.0;
            for(k=0; k<disArray[j].size(); k++){
                barrierFactor += disArray[j][k].energyFactor;
                cross(gb->n, gbd->planeList[disArray[j][k].planeID], crossProd);
                VECTOR_COPY(vec, burgList[disArray[j][k].burgID]);
                if(Normal(crossProd) < 1.0E-2){
                    if(BurgType_FCC(vec) == BASIC_BURG)break;
                }else{
                    if(BurgType_FCC(vec) != BASIC_BURG)break;
                }
            }

            if(k != disArray[j].size() || j > 0 && 
               barrierFactor > 0){
                itDisList = disArray.begin() + j;
                disArray.erase(itDisList);
            }else{
                j++;
            }
        }
        printf("\nDislocation Arrray:\n");
        for(i=0; i<disArray.size(); i++){
            printf("%d: ",i);
            for(j=0; j<disArray[i].size(); j++){
#if 1
                printf("(%d<%d>)[%d] ", disArray[i][j].planeID, 
                        gb->burgData.planeInGrain[disArray[i][j].planeID], 
                        disArray[i][j].burgID);
#else
                printf("(%5.3f,%5.3f,%5.3f)[%5.3f,%5.3f,%5.3f] ",
                       gb->burgData.planeList[disArray[i][j].planeID][0],
                       gb->burgData.planeList[disArray[i][j].planeID][1],
                       gb->burgData.planeList[disArray[i][j].planeID][2],
                       burgList[disArray[i][j].burgID][0],
                       burgList[disArray[i][j].burgID][1],
                       burgList[disArray[i][j].burgID][2]);
#endif
            }
            printf("\n");
        }
    }else{
        printf("zero results\n");
    }
#endif


    return;

}

void FreeHome(Home_t *home){
    int i,j;
    for(i=0; i<1 ; i++){
        free(home->GBKeys[i]); home->GBKeys[i] = NULL;
    }

    for(i=0; i<2; i++){
        free(home->grainKeys[i]->GB); home->grainKeys[i]->GB = NULL;
        free(home->grainKeys[i]); home->grainKeys[i] = NULL;
    }
    
    free(home->GBKeys);
    free(home->grainKeys);

    home->GBKeys = NULL;
    home->grainKeys = NULL;

    return;
}

void F(Home_t *home){
    int     i, j;
    GB_t    *GB;
    Grain_t *grain;
    int     GBID;
    for(i=0; i<11; i++){
        GB = home->GBKeys[i];
        GB->GBID = i;
        VECTOR_ZERO(GB->n);
        VECTOR_ZERO(GB->p);
        GB->n[1] = (real8)i;
        if(i==5){
            GB->n[1] = 1;
        }
    }

    for(i=0; i<2; i++){
        grain = home->grainKeys[i];
        grain->nBounds = 6;
        grain->GB = (GB_t **)malloc(6*sizeof(GB_t *));

        for(j=0; j<6; j++){
            GBID = (i==0) ? 1 : 0;
            GBID += (j-1+i*6);
            grain->GB[j] = home->GBKeys[GBID];
        }
    }

    return;
}

void G(Home_t *home){
    printf("1: %f\n", home->GBKeys[5]->n[1]);
    printf("2: %f %f\n", home->grainKeys[0]->GB[5]->n[1], 
            home->grainKeys[1]->GB[0]->n[1]);
    return;
}

void F1(int& inGrain){
    inGrain = 2;
    return;
}

void F2(real8 M1[3][3], real8 M2[3][3]){
    memcpy(M1,M2,sizeof(real8 [3][3]));
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
    free(burgData->planeInGrain); 
    free(burgData->keyLineList);

    for(i=0; i<burgData->numBurgVectors; i++){
        free(*(burgData->burgOnPlaneID+i)); 
    }

    for(i=0; i<burgData->numPlanes; i++){
        free(*(burgData->planeOnBurgID+i));
    }

    for(i=0; i<burgData->numKeyLines; i++){
        free(*(burgData->keyLineOnPlaneID+i));
    }

    free(burgData->numPlanesPerBurg); 
    free(burgData->numBurgsPerPlane); 
    free(burgData->numPlanesPerKeyLine); 

    burgData->burgList = NULL;
    burgData->planeList = NULL;
    burgData->planeType = NULL;
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
    vector< vector<real8> >::iterator  it;

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
    if(nRows != burgType.size())printf("in %s at %d", __FILE__, __LINE__);
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

    return;
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


void SumOfValuesIsConstant(int n, int sum, vector<int>& valueList, int tot, 
                            vector<vector<int> >& valueLists, int maxVar)
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

int BurgID_FCC(BurgInfo_t *burgData, real8 burg[3])
{
        int     i;
        real8   crossProd[3];
        
        for(i=0; i<burgData->numBurgVectors; i++){
            cross(burgData->burgList[i], burg, crossProd);
            if(Normal(crossProd) < 1.0E-3)break;
        }
        
        return(i);
}

#define DEBUG_DisassembleDislocation 1
bool cmp(const vector<int >& p, const vector<int >& q)
{
    if(p.back() != q.back()){
        return(p.back()<q.back());
    }else{
        int sum1 = accumulate(p.begin(), p.end()-1, 0);
        int sum2 = accumulate(q.begin(), q.end()-1, 0);
        return(sum1 < sum2);
    }
}

bool cmp2(Dislocation_t& p, Dislocation_t& q)
{
    return(p.planeID < q.planeID);
}

void Iter(vector<int>& disLocalBurgID, 
          vector< vector<int> >& burgOnPlaneID,
          vector<int>& planeIDList, 
          vector<int>& lastPlaneIDList,
          vector< vector<int> >& planeIDLists,
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
#if 1
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


void CompressVectorList(vector< vector<int> >& list)
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

void CompressDislocationArray(vector< vector<Dislocation_t> >& list)
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
void IterFun1(int n, int m, vector<int>& list,
       vector<vector<int> >& lists)
{
    if(m < 1 || n < 0)return;
    int i;

    
    for(i=0; i< n+1; i++){
        list[m-1] = i;
        IterFun1(n-i,m-1,list,lists);
    }
    if(m==1){
        list[0] = n; 
        vector<int>  IDList(list.size());
        vector<int>  sortList(list.size());
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
            vector<int> a;
            a.push_back(n);
            lists.push_back(a);
        }
    }

    return;
}


void IterFun2(vector<int>& oriVec, int iLayer, 
             vector<int>& list, 
             vector< vector<int> >& lists)
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

int main(){
    Home_t  home;

    InitHome(&home);

    FreeHome(&home);
    return(0);
}
