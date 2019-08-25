#ifndef _ParadisGen_h
#define _ParadisGen_h

#include <vector>

/*
 *      Define a structure containing all items corresponding to
 *      all command line options that have associated values.  This
 *      gives us an easy way to pass the command line arg values around
 *      to various functions.
 */
typedef struct {
        int     cubeLength;
        int     frLenMin;
        int     frLenMax;
        int     hexl;
        int     interstitialLoops;
        int     loopType;
        real8   maxSegLen;
        int     numChains;
        int     numFRSrcs;
        int     numLoops;
        char    *outputFile;
        int     pbcVal;
        real8   radius;
        int     seed;
        int     type;
        real8   xSurf[2];
        real8   ySurf[2];
        real8   zSurf[2];
/*
 *      The following are not input variables, but it is useful to include
 *      them in this structure to pass them around with the values from
 *      the command line arguments.
 */		
	    int	    polycrystal;	
	    real8	twinRatio;
	    int	    nGrain;
	    real8	heightRatio;
	    int	    deltaNums;
	    int	    nLayers;
	    real8   xMatDir[3];
	    real8   yMatDir[3];
	    real8   zMatDir[3];
	    real8	TBAngle;
	    int	    srcLenType;
	    real8	maxDensity; // l_tot/volume, unit=1/b^2
        int     extended;
	    int     burgID, planeID;
        real8   TBCryst[3];
        char    outFile[256];
        real8   loadDir[3]; 
        char    *voroFile;

        std::vector<std::vector<int> >    disList;
}InArgs_t;

/*
 *      Prototype the various global functions needed during the initial
 *      problem generation.
 */
void  CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreateScrewConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreateFiniteMixedConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  CreatePrismaticLoop(Home_t *home, InData_t *inData, int cubeLength,
          int loopType, int useinterstitial, int numSegs, real8 radius,
          int seed, real8 *totDislocLen, int dislocType);
void CreateFRSource(Home_t *home, InData_t *inData, InArgs_t *inArgs, real8 &totDislocLen);
void CreateFRSource(Home_t *home, InData_t *inData, InArgs_t *inArgs,
					real8 *totDislocLen);
void CreateFCCConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
						 real8 *totDislocLen);
void  CreateFCCIrradConfig(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, int numLoops, real8 hexl,
          real8 *totDislocLen, int dislocType);
void  CreateFCCPerfectLoop(Home_t *home, InData_t *inData, int cubeLength,
          int numChains, int seed, real8 *totDislocLen, int dislocType);
void  InitRemesh(InData_t *inData, int domValue, int startIndex);
real8 randm(int *seed);
void  WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock);
int GetBurgID_FCCFromArm(Home_t *home, Node_t *node, int segID);


/*
 *      For each type of nodal configuration that can be created
 *      define a name and integer id to be associated with the type
 */
#define FTYPE_SCREW		0
#define FTYPE_FINITE_MIXED	1
#define FTYPE_PRISMATIC_LOOP	2
#define FTYPE_FRANK_READ	3
#define FTYPE_FCC		4
#define FTYPE_FCC_IRRAD		5
#define FTYPE_FCC_PERFECT_LOOP	6
#define FTYPE_EDGE		7
#define FTYPE_MAX		8

#define FNAME_SCREW		"screw"
#define FNAME_EDGE		"edge"
#define FNAME_FINITE_MIXED	"finite-mixed"
#define FNAME_PRISMATIC_LOOP	"prismatic-loop"
#define FNAME_FRANK_READ	"frank-read-src"
#define FNAME_FCC		"fcc"
#define FNAME_FCC_IRRAD		"fcc-irrad"
#define FNAME_FCC_PERFECT_LOOP	"fcc-perfect-loop"


/*
 *      Define a structure to hold a nodal configuration type, name, and
 *      a pointer to the function to invoke to create that type of nodal
 *      configuration.
 */
typedef struct {
        int   funcType;
        const char  *funcName;
/*
        void  (* func)();
*/
} FuncData_t;


/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.
 */
typedef	enum{
	OPT_CUBEL = 0,
	OPT_HEIGHTRATIO,
	OPT_DELTA_NUMS,
	OPT_FRLEN,
	OPT_HELP,
	OPT_HEXSIZE,
	OPT_LOOPTYPE,
	OPT_MAXSEG,
	OPT_NCHAINS,
	OPT_NFRSRCS,
	OPT_NLOOPS,
	OPT_OUTFILE,
	OPT_PBC,
	OPT_RADIUS,
	OPT_SEED,
	OPT_TYPE,
	OPT_VACANCY,
	OPT_XSURF,
	OPT_YSURF,
	OPT_ZSURF,
	OPT_POLY,
	OPT_NGRAIN,
	OPT_TWIN_RATIO,
	OPT_NLAYERS,
	OPT_XMATDIR,
	OPT_YMATDIR,
	OPT_ZMATDIR,
	OPT_TBANGLE,
	OPT_SRCLENTYPE,
	OPT_MAXDENSITY,
    OPT_EXTENDED,
    OPT_BURGID,
    OPT_PLANEID,
    OPT_TBCRYST,
    OPT_OUT_FILE,
    OPT_LOAD_DIR,
    OPT_DIS_LIST,
    OPT_VORO_FILE,
	OPT_MAX
	
}ParadisGenArgs_t;



/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and a flag indicating if the option is paired with
 *      a value or not.
 */
typedef struct {
        int   optType;
        const char  *optName;
        int   optMinAbbrev;
        int   optPaired;
} Option_t;




#endif /* _ParadisGen_h */
