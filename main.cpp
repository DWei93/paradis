#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define VECTOR_ZERO(a) {a[0]=0; a[1]=0; a[2]=0;}
#define VECTOR_COPY(a,b) {(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];}

#define real8   double

#define NUM_BURGS_IN_FCC        21
#define NUM_PLANES_IN_FCC       13
#define NUM_KEY_LINES_IN_FCC    6


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

void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3]);
void GenerateGBBurgInfo(Home_t *home, GB_t *GB);
void FreeBurgInfoStr(BurgInfo_t *burgData);

void InitGrain(Home_t *home);

void Matrix33Vector3Multiply(real8 A[3][3], real8 x[3], real8 y[3])
{
        y[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
        y[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
        y[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];

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
    int     i;
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

    real8 rotMatrix[3][3] = {
        {7.071067811865475e-01,-7.071067811865475e-01,0.000000000000000e+00},
        {5.773502691896258e-01,5.773502691896258e-01,5.773502691896258e-01},
        {-4.082482904638631e-01,-4.082482904638631e-01,8.164965809277261e-01}
    };

    real8 rotTwin[3][3] = {
        {-7.071067811865476e-01,7.071067811865474e-01,-7.070501591499380e-17},
        {-5.773502691896257e-01,-5.773502691896260e-01,-5.773502691896258e-01},
        {-4.082482904638631e-01,-4.082482904638631e-01,8.164965809277261e-01}
    };

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

        GenerateFCCBurgInfo(&(grain->burgData), grain->orientation);
    }
    
    GenerateGBBurgInfo(home, gb); 

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

void GenerateFCCBurgInfo(BurgInfo_t *burgData, real8 rotMatrix[3][3])
{
        int         i, j;

        
        burgData->numBurgVectors = NUM_BURGS_IN_FCC;
        burgData->numPlanes   = NUM_PLANES_IN_FCC;
        burgData->numKeyLines = NUM_KEY_LINES_IN_FCC;

        burgData->burgList = (real8 (*)[3])calloc(3*burgData->numBurgVectors, sizeof(real8));
        burgData->planeList = (real8 (*)[3])calloc(3*burgData->numPlanes, sizeof(real8));
        burgData->planeType = (int *)calloc(burgData->numPlanes, sizeof(int));
        burgData->keyLineList = (real8 (*)[3])calloc(3*burgData->numKeyLines, sizeof(real8));

        burgData->numPlanesPerBurg = (int *)calloc(burgData->numBurgVectors, sizeof(int));        
        burgData->numBurgsPerPlane = (int *)calloc(burgData->numPlanes, sizeof(int));        
        burgData->numPlanesPerKeyLine = (int *)calloc(burgData->numKeyLines, sizeof(int));        
        

        real8 burgList[NUM_BURGS_IN_FCC][3] = {
/* 
 *          basic Roman-Roman burgers vectors       (6)
 *          magnitude: 1.0
 */
           { 0.7071067811865475,    0.7071067811865475,     0.0e+00   },   // DC -- 0
           { 0.7071067811865475,    -0.7071067811865475,    0.0e+00   },   // BA -- 1
           { 0.7071067811865475,    0.0e+00   ,    0.7071067811865475 },   // DA -- 2
           { 0.7071067811865475,    0.0e+00   ,    -0.7071067811865475},   // BC -- 3
           { 0.0e+00   ,   0.7071067811865475,     0.7071067811865475 },   // DB -- 4
           { 0.0e+00   ,   0.7071067811865475,     -0.7071067811865475},   // AC -- 5
/* 
 *          Shocklay Roman-Greek burgers vectors    (12)
 *          magnitude: 1/Sqrt(3)
 */
           { 0.4714045207910317,       0.2357022603955158,        0.2357022603955158 },      // Db -- 6
           {-0.4714045207910317,       0.2357022603955158,        0.2357022603955158 },      // dB -- 7
           { 0.4714045207910317,       -0.2357022603955158,       0.2357022603955158 },      // gA -- 8
           { 0.4714045207910317,       0.2357022603955158,        -0.2357022603955158},      // aC -- 9
           { 0.2357022603955158,       0.4714045207910317,        0.2357022603955158 },      // Da -- 10
           {-0.2357022603955158,       0.4714045207910317,        0.2357022603955158 },      // gB -- 11
           { 0.2357022603955158,       -0.4714045207910317,       0.2357022603955158 },      // dA -- 12
           { 0.2357022603955158,       0.4714045207910317,        -0.2357022603955158},      // bC -- 13
           { 0.2357022603955158,       0.2357022603955158,        0.4714045207910317 },      // Dg -- 14
           {-0.2357022603955158,       0.2357022603955158,        0.4714045207910317 },      // aB -- 15
           { 0.2357022603955158,       -0.2357022603955158,       0.4714045207910317 },      // bA -- 16
           { 0.2357022603955158,       0.2357022603955158,        -0.4714045207910317},      // dC -- 17
/* 
 *          compound Roman-Roman burgers vectors    (3) 
 *          magnitude: Sqrt(2)
 */
           { 1.414213562373095,         0.0e+00,                    0.0e+00           },   // DC+BA -- 18
           { 0.0e+00,                   1.414213562373095,          0.0e+00           },   // DC+AB -- 19
           { 0.0e+00,                   0.0e+00,                    1.414213562373095 }    // CA+DB -- 20
        };

        for(i=0; i<NUM_BURGS_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, burgList[i], burgData->burgList[i]);
        }

        int numPlanesPerBurg[NUM_BURGS_IN_FCC] = { 3, 3, 3, 3, 3, 3, 
                                                   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                   4, 4, 4 }; 
        memcpy(burgData->numPlanesPerBurg,numPlanesPerBurg,sizeof(int [NUM_BURGS_IN_FCC]));

        int burgFirstPlaneIndex[NUM_BURGS_IN_FCC] = { 0, 3, 6, 9, 12, 15, 
                                                     18,20,22,24,26,28,30,32,34,36,38,40,
                                                     42,46,50};

        real8 planeList[NUM_PLANES_IN_FCC][3] = {
            {-0.5773502691896258,       0.5773502691896258,       -0.5773502691896258},    // a plane  ID: 0
            { 0.5773502691896258,      -0.5773502691896258,       -0.5773502691896258},    // b plane  ID: 1
            {-0.5773502691896258,      -0.5773502691896258,        0.5773502691896258},    // g plane  ID: 2
            { 0.5773502691896258,       0.5773502691896258,        0.5773502691896258},    // d plane  ID: 3

            { 1.0e+00,                  0.0e+00,                   0.0e+00           },    // (100) junction plane     ID: 4
            { 0.0e+00,                  1.0e+00,                   0.0e+00           },    // (010) junction plane     ID: 5
            { 0.0e+00,                  0.0e+00,                   1.0e+00           },    // (001) junction plane     ID: 6
              
            { 0.0e+00,                  0.7071067811865475,        0.7071067811865475},    // (011) junction plane     ID: 7
            { 0.0e+00,                 -0.7071067811865475,        0.7071067811865475},    // (0-11) junction plane    ID: 8

            { 0.7071067811865475,       0.0e+00,                   0.7071067811865475},    // (101) junction plane     ID: 9
            {-0.7071067811865475,       0.0e+00,                   0.7071067811865475},    // (-101) junction plane    ID: 10

            { 0.7071067811865475,       0.7071067811865475,        0.0e+00           },    // (110) junction plane     ID: 11
            {-0.7071067811865475,       0.7071067811865475,        0.0e+00           }     // (-110) junction plane    ID: 12
        };

        for(i=0; i<NUM_PLANES_IN_FCC; i++){
            Matrix33Vector3Multiply(rotMatrix, planeList[i], burgData->planeList[i]);
        }
        int planeType[NUM_PLANES_IN_FCC] = { 0, 0, 0, 0,                          // basic glide planes
                                             1, 1, 1, 1, 1, 1, 1, 1, 1 };         // junction planes
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
           { 0.7071067811865475,    0.7071067811865475,     0.0e+00   },   // DC -- 0
           { 0.7071067811865475,    -0.7071067811865475,    0.0e+00   },   // BA -- 1
           { 0.7071067811865475,    0.0e+00   ,    0.7071067811865475 },   // DA -- 2
           { 0.7071067811865475,    0.0e+00   ,    -0.7071067811865475},   // BC -- 3
           { 0.0e+00   ,   0.7071067811865475,     0.7071067811865475 },   // DB -- 4
           { 0.0e+00   ,   0.7071067811865475,     -0.7071067811865475}    // AC -- 5
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
    burgData->keyLineList = NULL;

    burgData->numPlanesPerBurg = NULL;
    burgData->numPlanesPerKeyLine = NULL;
    burgData->numBurgsPerPlane = NULL;

    burgData->burgOnPlaneID = NULL;
    burgData->planeOnBurgID = NULL;
    burgData->keyLineOnPlaneID = NULL;

    return;
}

int CompressDirectionList(real8 (*dirList)[3], int *IDList, int nRows)
{
    int     i, j;
    real8   crossProdt[3];
    int     newNRows = 0;

    for(i=0; i<nRows; i++){
        IDList[i] = i;
    }

    for(i=0; i<nRows; i++){
        for(j=0; j<newNRows; j++){
            cross(dirList[i], dirList[j], crossProdt);
            if(Normal(crossProdt) < 1.0E-2)break;

        }

        IDList[i] = j;

        if(j == newNRows){
            VECTOR_COPY(dirList[newNRows], dirList[i])
            newNRows++;
        }
    }
    return(newNRows);

}

void AddDirectionList(real8 (*list1)[3], int n1, real8 (*list2)[3], int n2, real8 (*list)[3])
{
    int     i, j;
    int     newRows;

    for(i=0; i<n1+n2; i++){
        if(i<n1){
            VECTOR_COPY(list[i], list1[i]);
        }else{
            VECTOR_COPY(list[i], list2[i-n1]);
        }
    }

    return;    
}

void GenerateGBBurgInfo(Home_t *home, GB_t *GB)
{

    int         i, newRows, numRows;
    Grain_t     *grain1, *grain2;
    BurgInfo_t  *bd1, *bd2;
    
    grain1 = home->grainKeys[GB->nbrGrain[0]];
    grain2 = home->grainKeys[GB->nbrGrain[1]];

    bd1 = &(grain1->burgData);
    bd2 = &(grain2->burgData);

    numRows = bd1->numPlanes + bd2->numPlanes + 1;
    real8 (*planeList)[3] = (real8 (*)[3])malloc(3*numRows*sizeof(real8));
    
    int *planeIDList  = (int *)malloc(numRows*sizeof(int));

    AddDirectionList(bd1->planeList, bd1->numPlanes,
                     bd2->planeList, bd2->numPlanes, planeList);
    VECTOR_COPY(planeList[numRows-1], GB->n);

    newRows = CompressDirectionList(planeList, planeIDList, numRows); 

    for(i = 0; i<numRows; i++){
        printf("{%-7.3f,%-7.3f,%-7.3f} -- %2d\n", planeList[i][0],
                planeList[i][1],planeList[i][2], planeIDList[i]);
    }

    return;
}

int main(){
    Home_t  home;

    InitHome(&home);
    printf("1\n");


    BurgInfo_t  burgData;
    real8       rotMatrix[3][3] = {
        {7.874747423221179e-01,-5.805866356860647e-01,-2.068881066360532e-01},
        {2.157549463427914e-01,5.740956585181632e-01,-7.898506048609546e-01},
        {5.773502691896258e-01,5.773502691896258e-01,5.773502691896258e-01}
    };

    GenerateFCCBurgInfo(&burgData, rotMatrix);

    FreeBurgInfoStr(&burgData);

    FreeHome(&home);
    return(0);
}
