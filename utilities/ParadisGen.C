/****************************************************************************
 *
 *      Module:         ParadisGen.c
 *      Description:    Contains the main routine for controlling
 *                      generation of nodal data for a ParaDiS run, as
 *                      well as functions for command line parsing,
 *                      sanity checking, and a small help package.
 *
 *      Included functions:
 *              main()
 *              CheckArgSanity()
 *              GetInArgs()
 *              InitDefaultValues()
 *              PrintHelp()
 *              Usage()
 *
 *      Usage:
 *
 *      This utility can generate initial dislocation configurations
 *      of various types.  Since the information needed to generate these
 *      different types of configurations varies, not all options are
 *      used for all types of configurations.  The command lines below
 *      show the command line options available for use when creating
 *      the specified types (as indicated by the <type> command line option)
 *      of dislocation configurations.
 *
 *      For <type> "screw" or "edge" or "fcc"
 *
 *          paradisgen [-cubel <cubelength>] [-help] [-maxseg <maxseglen]   \
 *                     [-nchains <numchains>] [-seed <seedval]              \
 *                     [-type <configtype>]
 *
 *      For <type> "finite-mixed"
 *
 *          paradisgen [-cubel <cubelength>] [-help] [-maxseg <maxseglen]   \
 *                     [-nchains <numchains>] [-pbc pbcVal]                 \
 *                     [-seed <seedval] [-type <configtype>]                \
 *                     [-xsurf min[,max]] [-ysurf min[,max]]                \
 *                     [-zsurf min[,max]]
 *
 *      For <type> "prismatic-loop"
 *
 *          paradisgen [-cubel <cubelength>] [-help] [-looptype <type>]     \
 *                     [-maxseg <maxseglen] [-nloops <numloops>]            \
 *                     [-radius <radius>] [-seed <seedval]                  \
 *                     [-type <configtype>] [-vacancy]
 *
 *      For <type> "frank-read-src"
 *
 *          paradisgen [-cubel <cubelength>] [-frlen minsize[,maxsize]]     \
 *                     [-help] [-maxseg <maxseglen] [-nfrsrcs <numsrcs>]    \
 *                     [-seed <seedval] [-type <configtype>]
 *
 *      For <type> "fcc-irrad"
 *
 *          paradisgen [-cubel <cubelength>] [-help] [-hexsize <size>]      \
 *                     [-maxseg <maxseglen] [-nchains <numchains>]          \
 *                     [-nloops <numloops>] [-seed <seedval]                \
 *                     [-type <configtype>]
 *
 *      Overview:
 *
 *      This utility will generate nodal data suitable for a ParaDiS run and
 *      write the data into a an output file.  The nodal data generated
 *      will be determined by the user supplied command line options (or
 *      default values if none are specified).  Not all options are applicable
 *      to each type of nodal configuration; unnecessary options will be
 *      ignored.  
 *      
 *      All options may be abbreviated to the shortest non-ambiguous
 *      abbreviation of the option.  For example, '-nc' is a valid abbreviation
 *      for '-nchains' however, '-n' could refer to either '-nchains' or
 *      '-nloops' and hence is invalid.  Most options take on default values
 *      if none are specified on the command line.  
 *      
 ***************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "Typedefs.h"
#include "InData.h"
#include "ParadisGen.h"
#include "Restart.h"
#include "Decomp.h"
#include "Mobility.h"
#include "Grain.h"


/*
 *      Define and initialize an array of structures for mapping
 *      the types of dislocation creation functions to a corresponging
 *      integer type
 */
FuncData_t funcData[FTYPE_MAX] = {
        {FTYPE_SCREW,            FNAME_SCREW},
        {FTYPE_FINITE_MIXED,     FNAME_FINITE_MIXED},
        {FTYPE_PRISMATIC_LOOP,   FNAME_PRISMATIC_LOOP},
        {FTYPE_FRANK_READ,       FNAME_FRANK_READ},
        {FTYPE_FCC,              FNAME_FCC},
        {FTYPE_FCC_IRRAD,        FNAME_FCC_IRRAD},
        {FTYPE_FCC_PERFECT_LOOP, FNAME_FCC_PERFECT_LOOP},
        {FTYPE_EDGE,             FNAME_EDGE}
};


/*
 *      Define and initialize an array of structures containing
 *      all the possible command line arguments, and some info
 *      determining how it is treated.
 *
 *      option    option  #characters     1 unless option
 *      type      name    in unique       has no associated
 *                        abbreviation    value
 */
Option_t optList[OPT_MAX] = {
        {OPT_CUBEL,     "cubel",    	3, 1},
	{OPT_HEIGHTRATIO,"heightRatio",	3, 1},
	{OPT_DELTA_NUMS,"deltaNums",	3, 1},
        {OPT_FRLEN,     "frlen",    	3, 1},
        {OPT_HELP,      "help",     	3, 0},
        {OPT_HEXSIZE,   "hexsize",  	3, 1},
        {OPT_LOOPTYPE,  "looptype", 	3, 1},
        {OPT_MAXSEG,    "maxseg",   	4, 1},
        {OPT_NCHAINS,   "nchains",  	3, 1},
        {OPT_NFRSRCS,   "nfrsrcs",  	3, 1},
        {OPT_NLOOPS,    "nloops",   	3, 1},
        {OPT_OUTFILE,   "outfile",  	4, 1},
        {OPT_PBC,       "pbc",      	3, 1},
        {OPT_RADIUS,    "radius",   	3, 1},
        {OPT_SEED,      "seed",     	3, 1},
        {OPT_TYPE,      "type",     	3, 1},
        {OPT_VACANCY,   "vacancy",  	3, 0},
        {OPT_XSURF,     "xsurf",    	3, 1},
        {OPT_YSURF,     "ysurf",    	3, 1},
        {OPT_ZSURF,     "zsurf",    	3, 1},
	{OPT_POLY,	"polycrystal",	3, 1},
	{OPT_NGRAIN,	"nGrain",	3, 1},
	{OPT_TWIN_RATIO,"twinRatio",	3, 1},
	{OPT_NLAYERS,   "nLayers",	3, 1},
	{OPT_XMATDIR,   "xMatDir",	3, 1},
	{OPT_YMATDIR,   "yMatDir",	3, 1},
	{OPT_ZMATDIR,   "zMatDir",	3, 1},
	{OPT_TBANGLE,   "TBAngle",	3, 1},
	{OPT_SRCLENTYPE,"srcLenType",	3, 1},
	{OPT_MAXDENSITY,"maxDensity",	4, 1}

};


/*---------------------------------------------------------------------------
 *
 *      Function:       Usage
 *      Description:    Print out a brief message indicating the possible
 *                      command line options and terminated.
 *
 *      Arguments:
 *          program    name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void Usage(char *program)
{
        printf("\n");
        printf("Usage:\n");
        printf("\n");
        printf("This tool can generate initial dislocation configurations\n");
        printf("of various types.  Since the information needed to generate\n");
        printf("these different types of configurations varies, not all \n");
        printf("options are used for all types of configurations.  The \n");
        printf("command lines below show the command line options available\n");
        printf("for use when creating the specified types (as indicated by\n");
        printf("the <type> command line option) of configurations.\n");
        printf("\n");
        printf("    For <type> 'screw' or 'edge' or 'fcc'\n");
        printf("\n");
        printf("    %10s [-cubel <cubelength>] [-help] \n", program);
        printf("               [-maxseg <maxseglen] [-nchains <numchains>]\n");
        printf("               [-seed <seedval] [-type <configtype>]\n");
        printf("\n");
        printf("    For <type> 'finite-mixed'\n");
        printf("\n");
        printf("    %10s [-cubel <cubelength>] [-help] \n", program);
        printf("               [-maxseg <maxseglen] [-nchains <numchains>] \n");
        printf("               [-pbc pbcVal] [-seed <seedval] \n");
        printf("               [-type <configtype>] [-xsurf min[,max]]\n");
        printf("               [-ysurf min[,max]] [-zsurf min[,max]]\n");
        printf("\n");
        printf("    For <type> 'prismatic-loop'\n");
        printf("\n");
        printf("    %10s [-cubel <cubelength>] [-help] \n", program);
        printf("               [-looptype <type>] [-maxseg <maxseglen]\n");
        printf("               [-nloops <numloops>] [-radius <radius>] \n");
        printf("               [-seed <seedval] [-type <configtype>] \n");
        printf("               [-vacancy]\n");
        printf("\n");
        printf("    For <type> 'frank-read-src (FCC materials)'\n");
        printf("\n");
        printf("    %10s 	   [-cubel <cubelength>] \n", program);
        printf("               [-srcLenType <FR souce length type>] \n");
        printf("               [-frlen minsize[,maxsize]] [-help] \n");
        printf("               [-maxseg <maxseglen] [-nfrsrcs <numsrcs>]\n");
        printf("               [-seed <seedval] [-type <configtype>]\n");
	printf("			   [-polycrystal <1 = ploy>]  [-nGrain <num grains>]\n");
	printf("			   [-twinRatio <twinRatio>] \n");
	printf("			   [-heightRatio <heightRatio>]  [-deltaNums <deltaNums of FR>]\n");
        printf("\n");
        printf("    For <type> 'fcc-irrad'\n");
        printf("\n");
        printf("    %10s [-cubel <cubelength>] [-help] \n", program);
        printf("               [-hexsize <size>] [-maxseg <maxseglen] \n");
        printf("               [-nchains <numchains>] [-nloops <numloops>] \n");
        printf("               [-seed <seedval] [-type <configtype>]\n");
        printf("\n");

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PrintHelp
 *      Description:    Print out a detailed description of the available
 *                      program options, what they represent, how they
 *                      relate to each other, which are interdependent, or
 *                      mutually exclusive, default values, etc.
 *
 *      Arguments:
 *          program    name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void PrintHelp(char *program)
{
    Usage(program);

    printf("    Options may be abbreviated to the shortest non-ambiguous\n");
    printf("    abbreviation of the option.  For example, '-nc' is a valid \n");
    printf("    abbreviation for '-nchains' however, '-n' could refer to \n");
    printf("    either '-nchains' or '-nloops' and hence is invalid.\n\n");
    printf("Options:\n\n");
    printf("  -cubel    Define the length (units of b) of a single side\n");
    printf("            of the cubic problem space to be generated.\n\n");
    printf("  -frlen    Length (units of b) of Frank-Read sources\n");
    printf("            to generate\n\n");
    printf("  -help     Prints this help information.\n\n");
    printf("  -hexsize  Define the size of the hexagonal interstitial loops\n");
    printf("            for the fcc-irrad node configuration.\n\n");
    printf("  -looptype Specifies the type of prismatic loops to generate.\n");
    printf("            0 == mixture of [1 1 1] and [1 0 0] types\n");
    printf("            1 == all [1 1 1] type loops\n");
    printf("            2 == all [1 0 0] type loops\n");
    printf("            This value only applies if the selected dislocation\n");
    printf("            type is prismatic-loops.\n\n");
    printf("  -maxseg   Define the maximum segment length (units of b).\n\n");
    printf("  -nchains  Define the number of chains to create.  This value\n");
    printf("            is ignored for the prismatic-loop and frank-read src\n");
    printf("            configurations and must be a multiple of 12 for\n");
    printf("            fcc-perfect-loop node configurations.\n\n");
    printf("  -nfrsrcs  Define the number of Frank-Read sources to create\n\n");
    printf("  -nloops   Define the number of loops to generate in the\n");
    printf("            fcc-irrad or prismatic loop configurations\n\n");
    printf("  -outfile  Name of the file into which to write the\n");
    printf("            nodal data. \n\n");
    printf("  -pbc      Integer value specifying dimensions in which \n");
    printf("            periodic boundaries are enabled.  Default is \n");
    printf("            no dimensions.  Values can be bitwise OR'ed\n");
    printf("            together to yield periodic boundaries in any\n");
    printf("            combination of dimensions.  Only applies if\n");
    printf("            configuration type is finitemixed\n");
    printf("              1 == periodic in X dimensions.\n");
    printf("              2 == periodic in Y dimensions.\n");
    printf("              4 == periodic in Z dimensions.\n\n");
    printf("  -radius   Radius (units of b) of loop created for the\n");
    printf("            prismatic-loop and rhombohedral-va node configur-\n");
    printf("            ations. This value is ignored for all other nodal\n");
    printf("            configurations.\n\n");
    printf("  -seed     Define value to seed the random number\n");
    printf("            generator.\n\n");
    printf("  -type     Indicate the type of dislocations to create\n");
    printf("            for the initial nodal configuration.  Type can be\n");
    printf("            specified via the integer type value or by type\n");
    printf("            name.  Valid types are:\n\n");
    printf("             0 or screw        Random screw dislocations\n");
    printf("             1 or finite-mixed Random screw and edge \n");
    printf("               dislocations terminating \n");
    printf("               at the edges of the problem\n");
    printf("               space rather than wrapping\n");
    printf("               around\n");
    printf("             2 or prismatic-loop Dislocation loops.  Unless\n");
    printf("               otherwise specified, loops\n");
    printf("               are interstitial loops\n");
    printf("             3 or frank-read-src Frank-Read sources\n");
    printf("             4 or fcc            Random line dislocations\n");
    printf("               suitable for use with \n");
    printf("               MobilityRule_FCC1.\n");
    printf("             5 or fcc-irrad      Mixture of random line\n");
    printf("               dislocations and hexagonal\n");
    printf("               interstitial loops\n");
    printf("             6 or fcc-perfect-loop Not fully implemented\n");
    printf("             7 or edge          Random dislocations in \n");
    printf("               [100] directions. Not\n");
    printf("               pure edge, but not screw.\n");
    printf("  -vacancy  Specifies prismatic loops are to be vacancy loops\n");
    printf("            rather than interstitial loops.  This option\n");
    printf("            applies only if the selected dislocation type\n");
    printf("            is prismatic-loops or rhombohedral-va.\n\n");
    printf("  -xsurf,   Used to specify the coordinates of the surfaces\n");
    printf("  -ysurf,   when free surfaces are used.  Applies only.\n");
    printf("  -zsurf    to the finite-mixed configuration type.  Surface.\n");
    printf("            coordinates default to -.5*cubel and .5*cubel if\n");
    printf("            not specified.  If only min surface coordinate\n");
    printf("            is specified, max defaults to -min value.\n");
    printf("\n");

    printf("\n");
 
    exit(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:      InitDefaultValues
 *      Description:   Set default values for all items that may be
 *                     specified via command line arguments.
 *
 *-------------------------------------------------------------------------*/
static void InitDefaultValues(InArgs_t *inArgs)
{
        inArgs->cubeLength = 35000;
        inArgs->heightRatio = 1;   
        inArgs->hexl       = 50;
        inArgs->loopType   = 1;
        inArgs->interstitialLoops = 1;
        inArgs->maxSegLen  = 200.0;
        inArgs->numChains  = 2;
        inArgs->numFRSrcs  = 16;
        inArgs->frLenMin   = 3500;
        inArgs->frLenMax   = 3500;
        inArgs->numLoops   = 1000;
        inArgs->outputFile = (char*)"paradis.data";
        inArgs->radius     = inArgs->maxSegLen / 2.0;
        inArgs->seed       = time(0) + getpid();
        inArgs->type       = FTYPE_SCREW;
        inArgs->pbcVal     = 0;  /* pbc in no dimensions */
        inArgs->polycrystal = 1;
        inArgs->nGrain		= 2;
        inArgs->twinRatio	= 0.5;
        	
        inArgs->deltaNums	= 0;
        inArgs->nLayers		= 1;
        inArgs->srcLenType	= 0;
		
        inArgs->xSurf[0]   = 0;
        inArgs->xSurf[1]   = 0;
        inArgs->ySurf[0]   = 0;
        inArgs->ySurf[1]   = 0;
        inArgs->zSurf[0]   = 0;
        inArgs->zSurf[1]   = 0;

	inArgs->xMatDir[0] =  1;
	inArgs->xMatDir[1] = -1;
	inArgs->xMatDir[2] =  0;
	    
	inArgs->yMatDir[0] =  1;
	inArgs->yMatDir[1] =  1;
	inArgs->yMatDir[2] =  1;

	inArgs->zMatDir[0] = -1;
	inArgs->zMatDir[1] = -1;
	inArgs->zMatDir[2] =  2;

	inArgs->TBAngle = 0;
	inArgs->maxDensity = 5.0E14;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       CheckArgSanity
 *      Description:    Check the current values for the input arguments
 *                      (possibly a combination of default and user supplied
 *                      values), and verify that the values are reasonable.
 *
 *                      If there are any problems or inconsistencies, just
 *                      exit with an error message.
 *
 *-------------------------------------------------------------------------*/
static void CheckArgSanity(InArgs_t *inArgs)
{
        if (inArgs->cubeLength < 1)
            Fatal("Cube length must be greater than zero.\n");
		

        if ((inArgs->type < 0) || (inArgs->type > FTYPE_MAX))
            Fatal("Invalid type specified.\n");

        if (inArgs->seed < 0)
            Fatal("Seed value must be greater than zero\n");

        if ((inArgs->type == FTYPE_FCC_PERFECT_LOOP) &&
            ((inArgs->numChains < 12) || ((inArgs->numChains % 12) != 0)))
            Fatal("For type %s configurations, nchains must be\n"
                  "a multiple of 12.\n", funcData[inArgs->type].funcName);

        if (inArgs->radius > (inArgs->cubeLength / 2.0))
            Fatal("Radius must be less than 1/2 the cube length.\n");

        if ((inArgs->hexl < 1) && (inArgs->type == FTYPE_FCC_PERFECT_LOOP))
            Fatal("For type %s configurations, the hexsize value\n"
                  "must be greater than zero\n",
                  funcData[inArgs->type].funcName);

        if ((inArgs->numLoops < 1)&&(inArgs->type == FTYPE_FCC_PERFECT_LOOP))
            Fatal("For type %s configurations, the numloops value\n"
                  "must be specified\n", funcData[inArgs->type].funcName);

        if (inArgs->type == FTYPE_FRANK_READ) {
                if (inArgs->numFRSrcs < 0) {
                        Fatal("Number of frank read sources must be > 0");
                }
                if (inArgs->frLenMin < 0.0) {
                    inArgs->frLenMin = (int) floor(0.5 * inArgs->cubeLength);
                }
                if (inArgs->frLenMax < inArgs->frLenMin) {
                    inArgs->frLenMax = (int) floor(1.1 * inArgs->frLenMin);
                }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PrintArgs
 *      Description:    Dump out the arguments (command line or default
 *                      values) that are being used to generate the problem.
 *                      Only those arguments appropriate to the type of
 *                      dislocation configuration being created will be
 *                      displayed.
 *
 *-------------------------------------------------------------------------*/
static void PrintArgs(InArgs_t *inArgs)
{
        printf("\n");
        printf("Dislocation Type:         %s\n",
               funcData[inArgs->type].funcName);
        printf("Cube Length:              %d\n",inArgs->cubeLength);
		
		
        printf("Maximum Segment Length    %lf\n",inArgs->maxSegLen);
        printf("Data File:                %s\n",inArgs->outputFile);

        switch (inArgs->type) {
            case FTYPE_EDGE:
            case FTYPE_SCREW:
            case FTYPE_FINITE_MIXED:
		printf("Number of Loops:          %d\n", inArgs->numLoops);
                printf("Loop Radius:              %lf\n", inArgs->radius);
                printf("Seed Value:               %d\n", inArgs->seed);
                break;

            case FTYPE_PRISMATIC_LOOP:
                printf("Loop type:                %s %s\n",
                (inArgs->interstitialLoops ? "interstitial loops of" :
                                             "vacancy loops of"),
                (inArgs->loopType == 2 ?  "[1 0 0] types" :
                 (inArgs->loopType == 1 ? "[1 1 1] types" :
                                          "mixed types")));
                printf("Number of Loops:          %d\n", inArgs->numLoops);
                printf("Loop Radius:              %lf\n", inArgs->radius);
                printf("Seed Value:               %d\n", inArgs->seed);
                break;

            case FTYPE_FRANK_READ:
                printf("Source length distribution type:        %d\n", inArgs->srcLenType);
                printf("Number of sources:        %d\n", inArgs->numFRSrcs);
                if (inArgs->frLenMax > inArgs->frLenMin) {
                    printf("Length of sources:        %d-%d\n",
                           inArgs->frLenMin, inArgs->frLenMax);
                } else {
                    printf("Length of sources:        %d\n", inArgs->frLenMin);
                }
                printf("Seed Value:               %d\n", inArgs->seed);
	            printf("HeightRatio:              %f\n",inArgs->heightRatio);
	            printf("deltaNums:                %d\n",inArgs->deltaNums);
	            printf("nLayers:                  %d\n",inArgs->nLayers);
                break;

            case FTYPE_FCC:
                printf("Number of Chains:         %d\n", inArgs->numChains);
                printf("Seed Value:               %d\n", inArgs->seed);
                break;

            case FTYPE_FCC_IRRAD:
                printf("Number of Chains:         %d\n", inArgs->numChains);
                printf("Number of Loops:          %d\n", inArgs->numLoops);
                printf("Size of Hexagonal Loops:  %d\n", inArgs->hexl);
                printf("Seed Value:               %d\n", inArgs->seed);
                break;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       GetInArgs
 *      Description:    Parse the command line arguments.  Verify all options
 *                      are valid keywords and that all options requiring
 *                      associated values have them, and some simple sanity
 *                      checks.  All user supplied values will be stored in
 *                      a the inArgs structure, and a more in-depth sanity
 *                      check on those values is done elsewhere.
 *      Arguments
 *              argc    count of command line args
 *              argv    command line args
 *              inArgs  structure to hold user supplied values.  This
 *                      structure should have been populated with any
 *                      appropriate default values before this function
 *                      is called.
 *
 *-------------------------------------------------------------------------*/
static void GetInArgs(int argc, char *argv[], InArgs_t *inArgs)
{
        int     setDefaultRadius = 1;
        int     i, j, k;
	int	existDir[3];
	int	orientationType, existCount;
        char    *argName;
        char    *argValue;
        char    *token;

	VECTOR_ZERO(existDir);
	orientationType = 0;

        for (i = 1; i < argc; i++) {
/*
 *          If the option doesn't begin with a '-' something
 *          is wrong, so notify the user and terminate.
 */
            if (argv[i][0] != '-') {
                Usage(argv[0]);
                exit(1);
            }

            argName = &argv[i][1];

/*
 *          Scan the array of valid options for the user supplied
 *          option name.  (This may be any unique abbreviation of
 *          one of any of the options.  If we don't find the option
 *          notify the user and terminate.
 */
            for (j = 0; j < OPT_MAX; j++) {
                if (!strncmp(argName, optList[j].optName,
                    optList[j].optMinAbbrev)) {
                    break;
                }
            }

            if (j == OPT_MAX) {
                Usage(argv[0]);
                exit(1);
            }

/*
 *          Verify that there is an associated value if the specified
 *          option is supposed to be paired with a value.
 */
            if (optList[j].optPaired) {
                if (i+1 >= argc) {
                    Usage(argv[0]);
                    exit(1);
                } else {
                    argValue = argv[++i];
                }
            }

/*
 *          Do any option-specific processing...
 */
            switch (j)  {
                case OPT_CUBEL:

/*
 *                  The default radius value is dependent on the cube
 *                  length and should be reset when the cube length
 *                  changes, UNLESS the user has explicitly set the
 *                  radius value himself in which case we don't reset
 *                  it.
 */
                    inArgs->cubeLength = atoi(argValue);
                    if (setDefaultRadius)
                        inArgs->radius = inArgs->cubeLength / 4;
                    break;
					
	        case OPT_HEIGHTRATIO:
	            inArgs->heightRatio = atof(argValue);
	            break;
	            	
	        case OPT_DELTA_NUMS:
	            inArgs->deltaNums = atoi(argValue);
	            break;
					
	        case OPT_FRLEN:
                    token = strtok(argValue, ",");
                    inArgs->frLenMin = atoi(token);
                    token = strtok(NULL, ",");
                    if (token == (char *)NULL) {
                        inArgs->frLenMax = inArgs->frLenMin;
                        break;
                    }
                    inArgs->frLenMax = atoi(token);
                    break;
                case OPT_HELP:
                    PrintHelp(argv[0]);
                    exit(0);
                case OPT_HEXSIZE:
                    inArgs->hexl = atoi(argValue);
                    break;
                case OPT_LOOPTYPE:
                    inArgs->loopType = atoi(argValue);
                    break;
                case OPT_MAXSEG:
                    inArgs->maxSegLen = atof(argValue);
                    break;
                case OPT_NCHAINS:
                    inArgs->numChains = atoi(argValue);
                    break;
                case OPT_NFRSRCS:
                    inArgs->numFRSrcs = atoi(argValue);
                    break;
                case OPT_NLOOPS:
                    inArgs->numLoops = atoi(argValue);
                    break;
                case OPT_OUTFILE:
                    inArgs->outputFile = argValue;
                    break;
                case OPT_PBC:
                    inArgs->pbcVal = atoi(argValue);
                    break;
                case OPT_RADIUS:
                    inArgs->radius = atof(argValue);
                    setDefaultRadius = 0;
                    break;
                case OPT_SEED:
                    inArgs->seed = atoi(argValue);
                    break;
                case OPT_TYPE:
/*
 *                  Just to be nice, we allow the user to specify the
 *                  type of nodal configuration by either the integer
 *                  id number, or the string name.  If the first character
 *                  of the specified type is not an integer, assume
 *                  the user named the configration type and look it
 *                  up in the table.  (See the help package or include
 *                  file for the names and ids of the configuration types)
 */
                    if ((argValue[0] < '0') || (argValue[0] > '9')) {
                        for (k = 0; k < FTYPE_MAX; k++) {
                            if (!strcmp(argValue,
                                funcData[k].funcName)) {
                                break;
                            }
                        }
                        inArgs->type = k;
                    } else {
                        inArgs->type = atoi(argValue);
                    }
                    break;
                case OPT_VACANCY:
                    inArgs->interstitialLoops = 0;
                    break;
                case OPT_XSURF:
                    token = strtok(argValue, ",");
                    inArgs->xSurf[0] = atof(token);
                    token = strtok(NULL, ",");
                    if (token == (char *)NULL) {
                        inArgs->xSurf[1] = -inArgs->xSurf[0];
                        break;
                    }
                    inArgs->xSurf[1] = atof(token);
                    break;
                case OPT_YSURF:
                    token = strtok(argValue, ",");
                    inArgs->ySurf[0] = atof(token);
                    token = strtok(NULL, ",");
                    if (token == (char *)NULL) {
                        inArgs->ySurf[1] = -inArgs->ySurf[0];
                        break;
                    }
                    inArgs->ySurf[1] = atof(token);
                    break;
                case OPT_ZSURF:
                    token = strtok(argValue, ",");
                    inArgs->zSurf[0] = atof(token);
                    token = strtok(NULL, ",");
                    if (token == (char *)NULL) {
                        inArgs->zSurf[1] = -inArgs->zSurf[0];
                        break;
                    }
                    inArgs->zSurf[1] = atof(token);
                        break;
	        case OPT_POLY:
	            inArgs->polycrystal	= atof(argValue);
	                break;
	        case OPT_NGRAIN:
	            inArgs->nGrain = atoi(argValue);
	                break;
	        case OPT_TWIN_RATIO:
	            inArgs->twinRatio = atof(argValue);
	                break;
	        case OPT_NLAYERS:
	            inArgs->nLayers = atoi(argValue);
	                break;
	        case OPT_XMATDIR:
	            token = strtok(argValue,",");
	            inArgs->xMatDir[0] = atof(token);
	            token = strtok(NULL,",");
                    if (token == (char *)NULL) Fatal("orientation wrong");
	            inArgs->xMatDir[1] = atof(token);
	            token = strtok(NULL,",");
                    if (token == (char *)NULL) Fatal("orientation wrong");
	            inArgs->xMatDir[2] = atof(token);
	            NormalizeVec(inArgs->xMatDir);
	            existDir[0] = 1;
	            orientationType = 1;
	                break;
	        case OPT_YMATDIR:
	            token = strtok(argValue,",");
	            inArgs->yMatDir[0] = atof(token);
	            token = strtok(NULL,",");
                    if (token == (char *)NULL) Fatal("orientation wrong");
	            inArgs->yMatDir[1] = atof(token);
	            token = strtok(NULL,",");
                    if (token == (char *)NULL) Fatal("orientation wrong");
	            inArgs->yMatDir[2] = atof(token);
	            NormalizeVec(inArgs->yMatDir);
	            existDir[1] = 1;
	            orientationType = 1;
	                break;
	        case OPT_ZMATDIR:
	            token = strtok(argValue,",");
	            inArgs->zMatDir[0] = atof(token);
	            token = strtok(NULL,",");
                    if (token == (char *)NULL) Fatal("orientation wrong");
	            inArgs->zMatDir[1] = atof(token);
	            token = strtok(NULL,",");
                    if (token == (char *)NULL) Fatal("orientation wrong");
	            inArgs->zMatDir[2] = atof(token);
	            NormalizeVec(inArgs->yMatDir);
	            existDir[1] = 1;
	            orientationType = 1;
	                break;
	        case OPT_TBANGLE:
	            inArgs->TBAngle = atof(argValue);
	            orientationType = 4;
	                break;
	        case OPT_SRCLENTYPE:
	            inArgs->srcLenType = atoi(argValue);
	                break;
	        case OPT_MAXDENSITY:
	            inArgs->maxDensity = atof(argValue);
	                break;
            }
        }

	NormalizeVec(inArgs->xMatDir);
	NormalizeVec(inArgs->yMatDir);
	NormalizeVec(inArgs->zMatDir);

	if(orientationType == 0)return;

	if(orientationType == 1){
	    existCount = 0;
	    for(i=0; i<3; i++){
	        if(existDir[i] == 0)orientationType = i+1;
		if(existDir[i] == 1)existCount++;
	    }
    
	    if(existCount == 2 && orientationType >= 1 && orientationType <= 3){
	        if(orientationType == 1){
	            cross(inArgs->yMatDir, inArgs->zMatDir, inArgs->xMatDir);
		    NormalizeVec(inArgs->xMatDir);
	        }else if (orientationType == 2){
	            cross(inArgs->zMatDir, inArgs->xMatDir, inArgs->yMatDir);
		    NormalizeVec(inArgs->yMatDir);
	        }else{
	            cross(inArgs->xMatDir, inArgs->yMatDir, inArgs->zMatDir);
		    NormalizeVec(inArgs->zMatDir);
	        }
	    }
	}else{
	    real8    a1[3] = {1,-1,0}, temR1[3][3];
	    real8    a2[3] = {1, 0,0}, temR2[3][3], temR3[3][3];
	    real8    a3[3] = {0, 0,1};
	    real8    matrixOrientation[3][3];
	    real8    angle, angleZTB = 9.553166181245092e-01;
	    int	     seed;

	    seed = inArgs->seed;
	    NormalizeVec(a1);
	    RotateWithAxis(a1, angleZTB, temR1);

	    NormalizeVec(a2);
	    RotateWithAxis(a2, inArgs->TBAngle, temR2);

	    Matrix33Mult33(temR2, temR1, temR3);

	    NormalizeVec(a3);
	    angle = 2 * PI * randm(&seed);
	    RotateWithAxis(a3, angle, temR1);

	    Matrix33Mult33(temR1, temR3, matrixOrientation);

	    VECTOR_COPY(inArgs->xMatDir,matrixOrientation[0])
	    VECTOR_COPY(inArgs->yMatDir,matrixOrientation[1])
	    VECTOR_COPY(inArgs->zMatDir,matrixOrientation[2])

	}

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       WriteInitialNodeData
 *      Description:    Writes the nodal data for all nodes contained
 *                      in the inData->node list to the specified data
 *                      file.  The first call to ths function will
 *                      open up the output file, write the version
 *                      and header data and so on before writing nodal
 *                      data.  The final call (i.e. lastBlock == 1) will
 *                      insert the final node count into the data file
 *                      and clean up.
 *      Arguments:
 *
 *          lastBlock       set to 1 if if the supplied nodal data
 *                          is the last block of dfata to be written
 *                          to the specified file, set to zero in
 *                          all other cases.
 *
 * 
 *-------------------------------------------------------------------------*/
void WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock)
{
        int             i, iArm, numSeg = 0;
        Node_t          *node;
        Param_t         *param;
        static int      totalNodeCount = 0;
        static FILE     *fp = (FILE *)NULL;
	    static FILE     *fpTecplot = (FILE *)NULL;
    	Node_t	        *nbr;
	    Tag_t	        nbrTag;
	    real8	        v1, v2, v3;
	    int	            burgID, planeID;
        BurgInfo_t      *burgData;


        param = inData->param;
/*
 *      If this is the first block of data being written to
 *      this file (i.e. fp == NULL), do some basic initialization
 */
        if (fp == (FILE *)NULL) {

            fp = fopen(param->node_data_file, "w");
            if (!fp) {
                Fatal("%s: error %d opening file %s\n",
                      "WriteInitialNodeData", errno,
                      param->node_data_file);
            }

/*
 *          Write data file parameters
 */
            WriteParam(home->dataParamList, -1, fp);
			
/*
 * 			Werite the grain orientations into the data file
 */	
			
/*
			
			if(param->polycrystal){
				fprintf(fp, "Grain = %d\n",home->newGrainKeyPtr);
				Grain_t	*grain;
				
				for(i =  0; i< home->newGrainKeyPtr; i++){
					fprintf(fp, "%d \n",i);
					grain = home->grainKeys[i];
					fprintf(fp, " %16.10e %16.10e %16.10e\n",grain->orientation[0][0],
						grain->orientation[0][1],grain->orientation[0][2]);
					fprintf(fp, " %16.10e %16.10e %16.10e\n",grain->orientation[1][0],
						grain->orientation[1][1],grain->orientation[1][2]);
					fprintf(fp, " %16.10e %16.10e %16.10e\n",grain->orientation[2][0],
						grain->orientation[2][1],grain->orientation[2][2]);
			
				}
				fprintf(fp, "# End of grain orientation\n\n");
			}
*/

/*
 *          Write the domain decomposition into the data file
 */
            fprintf(fp, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");

            fprintf(fp, "domainDecomposition = \n");
            WriteDecompBounds(home, fp);

			

			
            fprintf(fp, "nodalData = \n");
            fprintf(fp, "#  Primary lines: node_tag, x, y, z, "
                    "num_arms, constraint, inGrain, onGB\n");
            fprintf(fp, "#  Secondary lines: arm_tag, burgx, burgy, "
                    "burgz, nx, ny, nz\n");
        }
		
		if(fpTecplot == (FILE *)NULL){
			fpTecplot = fopen("DDInitial.plt","w");
			fprintf(fpTecplot,"variables = X,Y,Z,V1,V2,V3,inGrain,onGB,burgID,planeID,constraint,numNbrs\n");
			
		}

/*
 *      Dump all the nodal data in this block to the file
 */
        totalNodeCount += inData->nodeCount;
		
		for (i = 0; i < inData->nodeCount; i++){
			node = &inData->node[i];
			for (iArm = 0; iArm < node->numNbrs; iArm++){
				numSeg++;
			}
		}
		fprintf(fpTecplot,"zone i=%d  F=POINT\n",numSeg);

        for (i = 0; i < inData->nodeCount; i++) {

            node = &inData->node[i];

            fprintf(fp,
                    " %d,%d %.4f %.4f %.4f %d %d %d %d\n",
                    node->myTag.domainID, node->myTag.index, 
                    node->x, node->y, node->z, node->numNbrs,
                    node->constraint,node->inGrain, node->onGB);

            for (iArm = 0; iArm < node->numNbrs; iArm++) {
	            if(home->param->materialType == MAT_TYPE_FCC){
                    real8 burg[3], plane[3];
                
                    if(home->param->polycrystal){
                        burgData = &(home->grainKeys[node->inGrain]->burgData);
                    }else{
                        burgData = &(home->burgData);
                    }

                    burg[0] = node->burgX[iArm];          
                    burg[1] = node->burgY[iArm];          
                    burg[2] = node->burgZ[iArm];          

                    plane[0] = node->nx[iArm];
                    plane[1] = node->ny[iArm];
                    plane[2] = node->nz[iArm];

                    burgID = BurgID_FCC(burgData, burg);
                    planeID = PlaneID_FCC(burgData, plane);
	            }else{
	                Fatal("The burgers ID of non-FCC material has not defined\n");
	            }

                fprintf(fp, "%d,%d %.15e %.15e %.15e\n"
			                "%.15e %.15e %.15e\n",
                        node->nbrTag[iArm].domainID,
                        node->nbrTag[iArm].index,
                        node->burgX[iArm], node->burgY[iArm],
                        node->burgZ[iArm], node->nx[iArm],
                        node->ny[iArm], node->nz[iArm]);
				
				nbrTag.domainID = node->nbrTag[iArm].domainID;
				nbrTag.index	= node->nbrTag[iArm].index;
				nbr	= GetinDataNodeFromTag(inData, nbrTag);
				if(nbr == (Node_t *)NULL)continue;
					
				v1 = 0.5 * (nbr->x - node->x);
				v2 = 0.5 * (nbr->y - node->y);	
				v3 = 0.5 * (nbr->z - node->z);	
				fprintf(fpTecplot, "%16.10e %16.10e %16.10e "
								   "%16.10e %16.10e %16.10e "
						           "%d %d %d %d %d %d\n",
						node->x, node->y, node->z, v1, v2, v3,
						node->inGrain,node->onGB,burgID, planeID,
					    node->constraint, node->numNbrs); 		
            }
        }


/*
 *      If this is the last block of data being written to the
 *      file, seek back to the node count in the file, overwrite
 *      it with the final total, and cleanup.
 */
		
        if (lastBlock) {
            fclose(fp);
	    fclose(fpTecplot);
            fp = (FILE *)NULL;
            printf("\nTotal node count:         %d\n", totalNodeCount);
        }

        return;
}


main(int argc, char *argv[])
{
        int             maxSide, minSide,maxHigh,minHigh, memSize;
        real8           cellVolume, bMag2, burgVolFactor;
        real8           totDislocDensity;
        real8           totDislocLen = 0.0;
        InArgs_t        inArgs;
        InData_t        inData;
        Param_t         *param;
        Home_t          home;


        Meminfo(&memSize);

        memset(&home, 0, sizeof(Home_t));
        memset(&inData, 0, sizeof(InData_t));

        inData.param = (Param_t *)calloc(1, sizeof(Param_t));
        home.dataParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));
        home.param = inData.param;
        home.numDomains = 1;

        InitDefaultValues(&inArgs);
        GetInArgs(argc, argv, &inArgs);

        CheckArgSanity(&inArgs);

        DataParamInit(inData.param, home.dataParamList);

/*
 *      Need to initialize some param_t structure values that are
 *      needed by most of the node generation functions.
 */
        param = inData.param;

        maxSide = inArgs.cubeLength / 2;
        minSide = -maxSide;

	if(inArgs.pbcVal == 0){
	    minHigh = minSide;
	    maxHigh = maxSide;
	}else{
	    maxHigh = (int)(inArgs.heightRatio * (real8)maxSide);
	    minHigh = -maxHigh;
	}
	
        param->minSideX = minSide;
        param->minSideY = minSide;
        param->minSideZ = minHigh;

        param->maxSideX = maxSide;
        param->maxSideY = maxSide;
        param->maxSideZ = maxHigh;
		
	param->xBoundMin = minSide;
	param->yBoundMin = minSide;
	param->zBoundMin = minHigh;
		
	param->xBoundMax = maxSide;
	param->yBoundMax = maxSide;
	param->zBoundMax = maxHigh;
		
        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

        param->invLx = 1.0 / param->Lx;
        param->invLy = 1.0 / param->Ly;
        param->invLz = 1.0 / param->Lz;

        param->maxSeg = inArgs.maxSegLen;
	param->minSeg = inArgs.maxSegLen/2.0;
		
	param->polycrystal = inArgs.polycrystal;
	param->nGrain = inArgs.nGrain;
	param->twinRatio = inArgs.twinRatio;

	param->xBoundType = ((inArgs.pbcVal & 0x01) == 0 ? Periodic : Free);
	param->yBoundType = ((inArgs.pbcVal & 0x02) == 0 ? Periodic : Free);
	param->zBoundType = ((inArgs.pbcVal & 0x04) == 0 ? Periodic : Free);

	VECTOR_COPY(param->matrixOrientation[0],inArgs.xMatDir);
	VECTOR_COPY(param->matrixOrientation[1],inArgs.yMatDir);
	VECTOR_COPY(param->matrixOrientation[2],inArgs.zMatDir);

	InitGrain(&home, param);
		
    strncpy(param->node_data_file, inArgs.outputFile,
            sizeof(param->node_data_file)-1);

#if 1
/*
 * Temporary code until [1 0 0] types are supported again
 */
if ((inArgs.type == FTYPE_PRISMATIC_LOOP) && (inArgs.loopType != 1)) {
    printf("\n***\n*** Currently only [1 1 1] type prismatic loops can be\n");
    printf("*** generated.  Please use the '-loopType 1' option\n***\n\n");
    exit(1);
}
#endif
        PrintArgs(&inArgs); 
/*
 *      And some additional initializations for when we write the
 *      data portion of the restart file.
 */
        param->minCoordinates[0] = minHigh;
        param->minCoordinates[1] = minHigh;
        param->minCoordinates[2] = minHigh;

        param->maxCoordinates[0] = maxHigh;
        param->maxCoordinates[1] = maxHigh;
        param->maxCoordinates[2] = maxHigh;

        param->decompType = param->dataDecompType;

/*
 *      Generate an initial domain decomposition.  Uses default
 *      geometry from param->dataDecompGeometry.
 */
        param->nXdoms = param->dataDecompGeometry[0];
        param->nYdoms = param->dataDecompGeometry[1];
        param->nZdoms = param->dataDecompGeometry[2];

        for (home.xMaxLevel = 0; param->nXdoms >> home.xMaxLevel != 1;
             home.xMaxLevel++);

        for (home.yMaxLevel = 0; param->nYdoms >> home.yMaxLevel != 1;
             home.yMaxLevel++);

        for (home.zMaxLevel = 0; param->nZdoms >> home.zMaxLevel != 1;
             home.zMaxLevel++);

        UniformDecomp(&home, &home.decomp);
		
		
		

/*
 *      All that's left is to invoke the proper function to
 *      generate nodal data for dislocations of the specified
 *      type.
 */
        switch (inArgs.type) {
        case FTYPE_SCREW:
                CreateScrewConfig(&home, &inData, inArgs.cubeLength,
                                  inArgs.numChains, inArgs.seed,
                                  &totDislocLen, inArgs.type);
                break;
        case FTYPE_EDGE:
                CreateEdges(&home, &inData, inArgs.cubeLength,
                            inArgs.numChains, inArgs.seed,
                            &totDislocLen, inArgs.type);
                break;
        case FTYPE_FINITE_MIXED:
/*
 *              Currently, this is the only configuration type that supports 
 *              a mix of periodic and free surfaces so only use the user
 *              provided boundary conditions here.  If no surface coordinates
 *              have been explicitly provided, use the appropriate simulation
 *              boundaries as the free surfaces.
 */
                param->xBoundType = ((inArgs.pbcVal & 0x01) > 0 ? Periodic : Free);
                param->yBoundType = ((inArgs.pbcVal & 0x02) > 0 ? Periodic : Free);
                param->zBoundType = ((inArgs.pbcVal & 0x04) > 0 ? Periodic : Free);

                if (inArgs.xSurf[0] != inArgs.xSurf[1]) {
                    param->xBoundMin = inArgs.xSurf[0];
                    param->xBoundMax = inArgs.xSurf[1];
                } else {
                    param->xBoundMin = param->minCoordinates[0];
                    param->xBoundMax = param->maxCoordinates[0];
                }

                if (inArgs.ySurf[0] != inArgs.ySurf[1]) {
                    param->yBoundMin = inArgs.ySurf[0];
                    param->yBoundMax = inArgs.ySurf[1];
                } else {
                    param->yBoundMin = param->minCoordinates[1];
                    param->yBoundMax = param->maxCoordinates[1];
                }

                if (inArgs.zSurf[0] != inArgs.zSurf[1]) {
                    param->zBoundMin = inArgs.zSurf[0];
                    param->zBoundMax = inArgs.zSurf[1];
                } else {
                    param->zBoundMin = param->minCoordinates[2];
                    param->zBoundMax = param->maxCoordinates[2];
                }

                CreateFiniteMixedConfig(&home, &inData, inArgs.cubeLength,
                                        inArgs.numChains, inArgs.seed,
                                        &totDislocLen, inArgs.type);
                break;
        case FTYPE_PRISMATIC_LOOP:
                CreatePrismaticLoop(&home, &inData, inArgs.cubeLength,
                                    inArgs.loopType,
                                    inArgs.interstitialLoops,
                                    inArgs.numLoops, inArgs.radius,
                                    inArgs.seed, &totDislocLen, inArgs.type);
                break;
        case FTYPE_FRANK_READ:
                CreateFRSource(&home, &inData, &inArgs,
                               &totDislocLen);
                break;
        case FTYPE_FCC:
			CreateFCCConfig(&home, &inData, inArgs.cubeLength,
							inArgs.numChains, inArgs.seed, &totDislocLen,
							inArgs.type);
                break;
        case FTYPE_FCC_IRRAD:
                CreateFCCIrradConfig(&home, &inData, inArgs.cubeLength,
                                     inArgs.numChains, inArgs.seed,
                                     inArgs.numLoops, inArgs.hexl,
                                     &totDislocLen, inArgs.type);
                break;
        case FTYPE_FCC_PERFECT_LOOP:
                printf("\n");
                printf("Support for %s type dislocations is %s.\n",
                       funcData[FTYPE_FCC_PERFECT_LOOP].funcName,
                       "not yet fully implemented");
                printf("\n");
                exit(1);
/*
                CreateFCCPerfectLoop();
                break;
*/
        }

/*
 *      Calculate the total dislocation density and print it.
 */
        cellVolume       = param->Lx * param->Ly * param->Lz;
        bMag2            = 2.556e-10 * 2.556e-10;
        burgVolFactor    = 1.0 / (bMag2 * cellVolume);
        totDislocDensity = totDislocLen * burgVolFactor;

        printf("Total density:            %e\n", totDislocDensity);
        Meminfo(&memSize);
        printf("\nEstimated memory usage = %d Kbytes\n\n", memSize);

        exit(0);
}


