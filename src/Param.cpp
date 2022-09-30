/**************************************************************************
 *
 *      Module:       Param.c
 *      Description:  Contains functions for binding control and data
 *                    file parameter names to code variables.
 *
 *      Public functions:
 *          CtrlParamInit()
 *          DataParamInit()
 *          MarkParamDisabled()
 *          MarkParamEnabled()
 *
 *************************************************************************/
#include <stdarg.h>
#include <ctype.h>
#include "Home.h"
#include "Parse.h"


/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamDisabled
 *      Description:  Explicitly mark a control file parameter 
 *                    to be disabled.  Primarily used to
 *                    prevent writing to the restart file any
 *                    control file parameters that are not
 *                    applicable or appropriate to the current
 *                    execution of the code.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure 
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly disabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamDisabled(ParamList_t *CPList, const char *name)
{
        int  paramIndex;

        if ((paramIndex = LookupParam(CPList,(char*) name)) >= 0) {
            CPList->varList[paramIndex].flags |= VFLAG_DISABLED;
        } else {
            Fatal("MarkParamDisabled: unknown parameter %s\n", name);
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamEnabled
 *      Description:  Explicitly mark a control file parameter 
 *                    to be enabled.  Primarily used when
 *                    explicitly controlling whether a control
 *                    parameter is to be written to a restart
 *                    file.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure 
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly enabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamEnabled(ParamList_t *CPList, const char *name)
{
        int  paramIndex;

        if ((paramIndex = LookupParam(CPList,(char*) name)) >= 0) {
            CPList->varList[paramIndex].flags &= ~VFLAG_DISABLED;
        } else {
            Fatal("MarkParamEnabled: unknown parameter %s\n", name);
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     CtrlParamInit
 *      Description:  Bind all valid control file parameters to
 *                    the associated code variables.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure 
 *                    associated with the control file parameters.
 *
 *----------------------------------------------------------------------*/
void CtrlParamInit(Param_t *param, ParamList_t *CPList)
{
/*
 *      Note: Parameters need only be initialized if their
 *      default values are non-zero.
 */
        CPList->paramCnt = 0;
        CPList->varList = (VarData_t *)NULL;

/*
 *      Simulation cell and processor setup
 */
        BindVar(CPList, "Simulation cell and processor setup",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "numXdoms", &param->nXdoms, V_INT, 1, VFLAG_NULL);
        param->nXdoms = 1;

        BindVar(CPList, "numYdoms", &param->nYdoms, V_INT, 1, VFLAG_NULL);
        param->nYdoms = 1;

        BindVar(CPList, "numZdoms", &param->nZdoms, V_INT, 1, VFLAG_NULL);
        param->nZdoms = 1;

        BindVar(CPList, "numXcells", &param->nXcells, V_INT, 1, VFLAG_NULL);
        param->nXcells = 3;   /* Must be >= 3 */

        BindVar(CPList, "numYcells", &param->nYcells, V_INT, 1, VFLAG_NULL);
        param->nYcells = 3;   /* Must be >= 3 */

#ifdef _BGP
        BindVar(CPList, "taskMappingMode", &param->taskMappingMode, V_INT, 1,
                VFLAG_NULL);
        param->taskMappingMode = 1;  /* use user-supplied decomp by default */
#endif

        BindVar(CPList, "numZcells", &param->nZcells, V_INT, 1, VFLAG_NULL);
        param->nZcells = 3;   /* Must be >= 3 */

#if 0
        BindVar(CPList, "xBoundType", &param->xBoundType, V_INT, 1, VFLAG_NULL);
        param->xBoundType = Periodic;

        BindVar(CPList, "yBoundType", &param->yBoundType, V_INT, 1, VFLAG_NULL);
        param->yBoundType = Periodic;

        BindVar(CPList, "zBoundType", &param->zBoundType, V_INT, 1, VFLAG_NULL);
        param->zBoundType = Periodic;
#endif
        BindVar(CPList, "decompType", &param->decompType, V_INT, 1, VFLAG_NULL);
        param->decompType = 1;

        BindVar(CPList, "DLBfreq", &param->DLBfreq, V_INT, 1, VFLAG_NULL);
        param->DLBfreq = 3;

#if 0
        BindVar(CPList, "xBoundMin", &param->xBoundMin, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "xBoundMax", &param->xBoundMax, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "yBoundMin", &param->yBoundMin, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "yBoundMax", &param->yBoundMax, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "zBoundMin", &param->zBoundMin, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "zBoundMax", &param->zBoundMax, V_DBL, 1, VFLAG_NULL);
#endif
/*
 *      Simulation time and timestepping controls
 */
        BindVar(CPList, "Simulation time and timestepping controls",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "cycleStart", &param->cycleStart, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "wallTime", &param->wallTime, V_DBL, 1, VFLAG_NULL);
        param->wallTime = -1;
        BindVar(CPList, "nextSubmit", &param->nextSubmit, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "dtGBCrossSlip", &param->dtGBCrossSlip, V_DBL, 1, VFLAG_NULL);
        param->dtGBCrossSlip = 5.0E-6;

        BindVar(CPList, "timeNextGBCrossSlip", &param->timeNextGBCrossSlip, V_DBL, 1, VFLAG_NULL);
        param->timeNextGBCrossSlip = 0.0E0;

        BindVar(CPList, "maxstep", &param->maxstep, V_INT, 1, VFLAG_NULL);
        param->maxstep = 100;

        BindVar(CPList, "maxStn", &param->maxStn, V_DBL, 1, VFLAG_NULL);
        param->maxStn = 1.0E-2;

        BindVar(CPList, "maxStr", &param->maxStr, V_DBL, 1, VFLAG_NULL);
        param->maxStr = -2.0E8;

        BindVar(CPList, "maxPst", &param->maxPst, V_DBL, 1, VFLAG_NULL);
        param->maxPst = 1.0E-2;

        BindVar(CPList, "maxDisDen", &param->maxDisDen, V_DBL, 1, VFLAG_NULL);
        param->maxDisDen = 1.0E20;

        BindVar(CPList, "timeNow", &param->timeNow, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "timeStart", &param->timeStart, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "timestepIntegrator", param->timestepIntegrator,
                V_STRING, 1, VFLAG_NULL);
        strcpy(param->timestepIntegrator, "trapezoid");

        BindVar(CPList, "deltaTT", &param->deltaTT, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "maxDT", &param->maxDT, V_DBL, 1, VFLAG_NULL);
        param->maxDT = 1.0e-08;

        BindVar(CPList, "lastDT", &param->lastDT, V_DBL, 1, VFLAG_NULL);
        param->lastDT = 1.0e-013;

        BindVar(CPList, "nextDT", &param->nextDT, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "dtIncrementFact", &param->dtIncrementFact,
                V_DBL, 1, VFLAG_NULL);
        param->dtIncrementFact = 1.2;

        BindVar(CPList, "dtDecrementFact", &param->dtDecrementFact,
                V_DBL, 1, VFLAG_NULL);
        param->dtDecrementFact = 0.5;

        BindVar(CPList, "dtExponent", &param->dtExponent, V_DBL, 1, VFLAG_NULL);
        param->dtExponent = 4.0;

        BindVar(CPList, "dtVariableAdjustment",
                &param->dtVariableAdjustment, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "rTol", &param->rTol, V_DBL, 1, VFLAG_NULL);
        param->rTol = -1.0;

        BindVar(CPList, "rmax", &param->rmax, V_DBL, 1, VFLAG_NULL);
        param->rmax = 100; 

/*
 *      Discretization controls and controls for topological changes
 */
        BindVar(CPList, "Discretization and topological change controls",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);
#if 0
        BindVar(CPList, "maxSeg", &param->maxSeg, V_DBL, 1, VFLAG_NULL);
        param->maxSeg = -1.0;
#endif
        BindVar(CPList, "minSeg", &param->minSeg, V_DBL, 1, VFLAG_NULL);
        param->minSeg = -1.0;
        
        BindVar(CPList, "rann", &param->rann, V_DBL, 1, VFLAG_NULL);
        param->rann = -1.0;

        BindVar(CPList, "remeshRule", &param->remeshRule, V_INT, 1, VFLAG_NULL);
        param->remeshRule = 2;

        BindVar(CPList, "splitMultiNodeFreq", &param->splitMultiNodeFreq,
                V_INT, 1, VFLAG_NULL);
        param->splitMultiNodeFreq = 1;

        BindVar(CPList, "collisionMethod", &param->collisionMethod, V_INT, 1,
                VFLAG_NULL);
        param->collisionMethod = 2;

/*
 *      FMM controls
 */
        BindVar(CPList, "Fast Multipole Method controls", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "fmEnabled", &param->fmEnabled, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "fmMPOrder", &param->fmMPOrder, V_INT, 1, VFLAG_NULL);
        param->fmMPOrder = 2;

        BindVar(CPList, "fmTaylorOrder", &param->fmTaylorOrder, V_INT,
                1, VFLAG_NULL);
        param->fmTaylorOrder = 5;

        BindVar(CPList, "fmCorrectionTbl", param->fmCorrectionTbl, V_STRING,
                1, VFLAG_NULL);
        strcpy(param->fmCorrectionTbl,"inputs/fm-ctab.Ta.600K.0GPa.m2.t5.dat");

/*
 *      Identify tables needed for remote force calculations if
 *      FMM is not enabled
 */
        BindVar(CPList, "Tables for non-FMM far-field force calcs",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "Rijmfile", param->Rijmfile, V_STRING, 1, VFLAG_NULL);
        strcpy(param->Rijmfile,"inputs/Rijm.cube.out");

        BindVar(CPList, "RijmPBCfile", param->RijmPBCfile, V_STRING,
                1, VFLAG_NULL);
        strcpy(param->RijmPBCfile,"inputs/RijmPBC.cube.out");

/*
 *      Loading condition parameters
 */
        BindVar(CPList, "Loading conditions", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "TempK", &param->TempK, V_DBL, 1, VFLAG_NULL);
        param->TempK = 600;

        BindVar(CPList, "loadType", &param->loadType, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "appliedStress", param->appliedStress, V_DBL, 6,
                VFLAG_NULL);

        BindVar(CPList, "eRate", &param->eRate, V_DBL, 1, VFLAG_NULL);
        param->eRate = 1.0;

        BindVar(CPList, "indxErate", &param->indxErate, V_INT, 1, VFLAG_NULL);
        param->indxErate = 1;

        BindVar(CPList, "edotdir",&param->edotdir, V_DBL, 3, VFLAG_NULL);
        param->edotdir[0] = 0.0;
        param->edotdir[1] = 0.0;
        param->edotdir[2] = 0.0;
        
        BindVar(CPList, "edotPolar",&param->edotPolar, V_DBL, 2, VFLAG_NULL);
        param->edotPolar[0] = 0.0;
        param->edotPolar[1] = 0.0;

        BindVar(CPList, "strRate",&param->strRate, V_DBL, 1, VFLAG_NULL);
        param->strRate = -1.0E15;

        BindVar(CPList, "loadDir",&param->loadDir, V_DBL, 6, VFLAG_NULL);
        param->loadDir[0] = 0.0;
        param->loadDir[1] = 0.0;
        param->loadDir[2] = 1.0;
        param->loadDir[3] = 0.0;
        param->loadDir[4] = 0.0;
        param->loadDir[5] = 0.0;

        BindVar(CPList, "plaRateCtr",&param->plaRateCtr, V_DBL, 1, VFLAG_NULL);
        param->plaRateCtr = -1.0;

        BindVar(CPList, "cTimeOld", &param->cTimeOld, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "dCyclicStrain", &param->dCyclicStrain, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "netCyclicStrain", &param->netCyclicStrain,
                V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "numLoadCycle", &param->numLoadCycle, V_INT, 1,
                VFLAG_NULL);

        BindVar(CPList, "eAmp", &param->eAmp, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "maxCycle", &param->maxCycle, V_DBL, 1, VFLAG_NULL);
/*
 *      To specify input sample axes in laboratory frame
 */
        BindVar(CPList, "useLabFrame", &param->useLabFrame, V_INT, 1,
                VFLAG_NULL);

        BindVar(CPList, "labFrameXDir", param->labFrameXDir, V_DBL, 3,
                VFLAG_NULL);
        param->labFrameXDir[0] = 1.0;
        param->labFrameXDir[1] = 0.0;
        param->labFrameXDir[2] = 0.0;

        BindVar(CPList, "labFrameYDir", param->labFrameYDir, V_DBL, 3,
                VFLAG_NULL);
        param->labFrameYDir[0] = 0.0;
        param->labFrameYDir[1] = 1.0;
        param->labFrameYDir[2] = 0.0;

        BindVar(CPList, "labFrameZDir", param->labFrameZDir, V_DBL, 3,
                VFLAG_NULL);
        param->labFrameZDir[0] = 1.0;
        param->labFrameZDir[1] = 0.0;
        param->labFrameZDir[2] = 0.0;

/*
 *      Parameters for material specific constants and mobility values
 */
        BindVar(CPList, "Material and mobility parameters",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "mobilityLaw", param->mobilityLaw, V_STRING, 1,
                VFLAG_NULL);
        strcpy(param->mobilityLaw, "BCC_0");

        BindVar(CPList, "vacancyConc", &param->vacancyConc, V_DBL, 1,
                VFLAG_NULL);
        param->vacancyConc = -1.0; /* signifies not used */

        BindVar(CPList, "vacancyConcEquilibrium",
                &param->vacancyConcEquilibrium, V_DBL, 1, VFLAG_NULL);
        param->vacancyConcEquilibrium = -1.0; /* signifies not used */

        BindVar(CPList, "shearModulus", &param->shearModulus, V_DBL, 1,
                VFLAG_NULL);
        param->shearModulus = 6.488424e+10; /* Ta: units=Pa @t=600K, p=0Pa*/

        BindVar(CPList, "pois", &param->pois, V_DBL, 1, VFLAG_NULL);
        param->pois = 3.327533e-01;    /* Ta: temp=600K, pressure=0Pa */

        BindVar(CPList, "massDensity", &param->massDensity, V_DBL,
                1, VFLAG_NULL);
        param->massDensity = 8900;
    
        BindVar(CPList, "outRadius", &param->outRadius, V_DBL, 1,
                VFLAG_NULL);
        param->outRadius = 1.0e4; 
        
        BindVar(CPList, "TBSFEnergy", &param->TBSFEnergy, V_DBL, 1,
                VFLAG_NULL);
        param->TBSFEnergy = 0.0; 
        
        BindVar(CPList, "TBEnergy", &param->TBEnergy, V_DBL, 1,
                VFLAG_NULL);
        param->TBEnergy = 0.0;  

        BindVar(CPList, "sfEnergy", &param->sfEnergy, V_DBL, 1,
                VFLAG_NULL);
        param->sfEnergy = 29.498E-3;  /* Jm^-2*/ 

        BindVar(CPList, "burgMag", &param->burgMag, V_DBL, 1, VFLAG_NULL);
        param->burgMag = 2.875401e-10;  /* Ta: units=m @t=600K, pressure=0Pa */

        BindVar(CPList, "YoungsModulus", &param->YoungsModulus, V_DBL, 1,
                VFLAG_NULL);
        param->YoungsModulus=2.0 * param->shearModulus * (1.0 + param->pois); /* Ta: units=Pa @t=600K, p=0Pa */

        BindVar(CPList, "elasticConstantMatrix", param->elasticConstantMatrix,
                V_DBL, 36, VFLAG_NULL);
        
        BindVar(CPList, "rc", &param->rc, V_DBL, 1, VFLAG_NULL);
        param->rc = -1.0;

        BindVar(CPList, "Ecore", &param->Ecore, V_DBL, 1, VFLAG_NULL);
        param->Ecore = -1.0;

        BindVar(CPList, "MobScrew", &param->MobScrew, V_DBL, 1, VFLAG_NULL);
        param->MobScrew = 10.0;

        BindVar(CPList, "MobEdge", &param->MobEdge, V_DBL, 1, VFLAG_NULL);
        param->MobEdge  = 10.0;

        BindVar(CPList, "MobTwinning", &param->MobTwinning, V_DBL, 1, VFLAG_NULL);
        param->MobTwinning  = 10.0;

        BindVar(CPList, "MobClimb", &param->MobClimb, V_DBL, 1, VFLAG_NULL);
        param->MobClimb = 1.0e-02;

/*
 *      List of burgers vectors/line directions to be considered sessile
 */
        BindVar(CPList, "sessileburgspec", param->sessileburgspec,
                V_DBL, 30, VFLAG_NULL);

        BindVar(CPList, "sessilelinespec", param->sessilelinespec,
                V_DBL, 30, VFLAG_NULL);

        BindVar(CPList, "includeInertia", &param->includeInertia, V_INT,
                1, VFLAG_NULL);


/*
 *      Flux decomposition
 */
        BindVar(CPList, "Flux decomposition",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "totstraintensor", param->totstraintensor, V_DBL, 6,  VFLAG_NULL);
        param->totstraintensor[0] = 0.0;
        param->totstraintensor[1] = 0.0;
        param->totstraintensor[2] = 0.0;
        param->totstraintensor[3] = 0.0;
        param->totstraintensor[4] = 0.0;
        param->totstraintensor[5] = 0.0;

        BindVar(CPList, "strain", &param->strain, V_DBL, 1,  VFLAG_NULL);
        param->strain = 0.0;

        BindVar(CPList, "accElast", &param->accElast, V_INT, 1,  VFLAG_NULL);
        param->accElast = 0;

        BindVar(CPList, "accPst", &param->accPst, V_DBL, 1,  VFLAG_NULL);
        param->accPst = 2.0E-5;

        BindVar(CPList, "totpStn", param->totpStn, V_DBL, 6, VFLAG_NULL);

        BindVar(CPList, "totpSpn", param->totpSpn, V_DBL, 6, VFLAG_NULL);

        BindVar(CPList, "totTwinpStn",  param->totTwinpStn, V_DBL, 6, VFLAG_NULL);

        BindVar(CPList, "Ltot", param->Ltot, V_DBL, 16, VFLAG_NULL);
        BindVar(CPList, "FCC_Ltot", param->FCC_Ltot, V_DBL, 24, VFLAG_NULL);

        BindVar(CPList, "fluxtot", param->fluxtot, V_DBL, 28, VFLAG_NULL);
        BindVar(CPList, "FCC_fluxtot", param->FCC_fluxtot, V_DBL, 42,
                VFLAG_NULL);

/*
 *      Total system density: this is informational only and recalculated
 *      each time a restart file is written.
 */
        BindVar(CPList, "Total density. Informational only; ignored on input",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "disloDensity", &param->disloDensity, V_DBL, 1, VFLAG_NULL);
        BindVar(CPList, "disloDensity_TB", &param->disloDensity_TB, V_DBL, 1, VFLAG_NULL);

/*
 *      Velocity statistics
 */
        BindVar(CPList, "Velocity statistics",
                (void *)NULL, V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "vAverage", &param->vAverage, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "vStDev", &param->vStDev, V_DBL, 1, VFLAG_NULL);


/*
 *      I/O controls and options
 */
        BindVar(CPList, "I/O controls and parameters", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "dirname", param->dirname, V_STRING, 1, VFLAG_NULL);

        BindVar(CPList, "writeBinRestart", &param->writeBinRestart, V_INT,
                1, VFLAG_NULL);

        BindVar(CPList, "skipIO", &param->skipIO, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "numIOGroups", &param->numIOGroups, V_INT, 1,
                VFLAG_NULL);
        param->numIOGroups = 1;

/* 
 *      segment/arm files
 */
        BindVar(CPList, "armfile", &param->armfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "armfilefreq", &param->armfilefreq, V_INT, 1,
                VFLAG_NULL);
        param->armfilefreq = 100;

        BindVar(CPList, "armfiledt", &param->armfiledt, V_DBL, 1, VFLAG_NULL);
        param->armfiledt = -1.0;

        BindVar(CPList, "armfiletime", &param->armfiletime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "armfilecounter", &param->armfilecounter, V_INT, 1,
                VFLAG_NULL);


/* 
 *      flux decomposition files
 */
        BindVar(CPList, "fluxfile", &param->fluxfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "fluxfreq", &param->fluxfreq, V_INT, 1, VFLAG_NULL);
        param->fluxfreq = 100;

        BindVar(CPList, "fluxdt", &param->fluxdt, V_DBL, 1, VFLAG_NULL);
        param->fluxdt = -1.0;

        BindVar(CPList, "fluxtime", &param->fluxtime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "fluxcounter", &param->fluxcounter, V_INT, 1,
                VFLAG_NULL);


/* 
 *      chain fragment files -- to be post-processed into a data format
 *      suitable for use with the VisIt visualization software.
 */
        BindVar(CPList, "fragfile", &param->fragfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "fragfreq", &param->fragfreq, V_INT, 1, VFLAG_NULL);
        param->fragfreq = 100;

        BindVar(CPList, "fragdt", &param->fragdt, V_DBL, 1, VFLAG_NULL);
        param->fragdt = -1.0;

        BindVar(CPList, "fragtime", &param->fragtime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "fragcounter", &param->fragcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      gnuplot files
 */
        BindVar(CPList, "gnuplot", &param->gnuplot, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "gnuplotfreq", &param->gnuplotfreq, V_INT, 1,
                VFLAG_NULL);
        param->gnuplotfreq = 100;

        BindVar(CPList, "gnuplotdt", &param->gnuplotdt, V_DBL, 1, VFLAG_NULL);
        param->gnuplotdt = -1.0;

        BindVar(CPList, "gnuplottime", &param->gnuplottime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "gnuplotcounter", &param->gnuplotcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      pole figure files
 */
        BindVar(CPList, "polefigfile", &param->polefigfile, V_INT, 1,
                VFLAG_NULL);

        BindVar(CPList, "polefigfreq", &param->polefigfreq, V_INT, 1,
                VFLAG_NULL);
        param->polefigfreq = 100;

        BindVar(CPList, "polefigdt", &param->polefigdt, V_DBL, 1, VFLAG_NULL);
        param->polefigdt = -1.0;

        BindVar(CPList, "polefigtime", &param->polefigtime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "polefilecounter", &param->polefigcounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      povray files
 */
        BindVar(CPList, "povray", &param->povray, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "povrayfreq", &param->povrayfreq, V_INT, 1, VFLAG_NULL);
        param->povrayfreq = 100;

        BindVar(CPList, "povraydt", &param->povraydt, V_DBL, 1, VFLAG_NULL);
        param->povraydt = -1.0;

        BindVar(CPList, "povraytime", &param->povraytime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "povraycounter", &param->povraycounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      postscript files
 */
        BindVar(CPList, "psfile", &param->psfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "psfilefreq", &param->psfilefreq, V_INT, 1, VFLAG_NULL);
        param->psfilefreq = 100;

        BindVar(CPList, "psfiledt", &param->psfiledt, V_DBL, 1, VFLAG_NULL);
        param->psfiledt = -1.0;

        BindVar(CPList, "psfiletime", &param->psfiletime, V_DBL, 1, VFLAG_NULL);


/*
 *      restart files
 */
        BindVar(CPList, "savecn", &param->savecn, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "savecnfreq", &param->savecnfreq, V_INT, 1, VFLAG_NULL);
        param->savecnfreq = 100;

        BindVar(CPList, "savecndt", &param->savecndt, V_DBL, 1, VFLAG_NULL);
        param->savecndt = -1.0;

        BindVar(CPList, "savecntime", &param->savecntime, V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "savecncounter", &param->savecncounter, V_INT, 1,
                VFLAG_NULL);


/*
 *      properties files
 */
        BindVar(CPList, "saveprop", &param->saveprop, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "savepropfreq", &param->savepropfreq, V_INT, 1,
                VFLAG_NULL);
        param->savepropfreq = 100;

        BindVar(CPList, "savepropdt", &param->savepropdt, V_DBL, 1, VFLAG_NULL);
        param->savepropdt = -1.0;

        BindVar(CPList, "saveproptime", &param->saveproptime, V_DBL, 1,
                VFLAG_NULL);


/*
 *      timer files
 */
        BindVar(CPList, "savetimers", &param->savetimers, V_INT, 1, VFLAG_NULL);
        param->savetimers = 0;

        BindVar(CPList, "savetimersfreq", &param->savetimersfreq, V_INT, 1,
                VFLAG_NULL);
        param->savetimersfreq = 100;

        BindVar(CPList, "savetimersdt", &param->savetimersdt, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "savetimerstime", &param->savetimerstime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "savetimerscounter", &param->savetimerscounter,
                V_INT, 1, VFLAG_NULL);
        param->savetimerscounter = 0;


/*
 *      tecplot files
 */
        BindVar(CPList, "tecplot", &param->tecplot, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "tecplotfreq", &param->tecplotfreq, V_INT, 1, VFLAG_NULL);

        param->tecplotfreq = 4000;

        BindVar(CPList, "tecplotdstrain", &param->tecplotdstrain, V_DBL, 1, VFLAG_NULL);
        param->tecplotdstrain = -1.0;

        BindVar(CPList, "tecplotstrain", &param->tecplotstrain, V_DBL, 1, VFLAG_NULL);
        param->tecplotstrain = 0.0;

        BindVar(CPList, "tecmesh", &param->tecmesh, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "tecmeshfreq", &param->tecmeshfreq, V_INT, 1,
                VFLAG_NULL);
        param->tecmeshfreq = 40000;

        BindVar(CPList, "tecplotdt", &param->tecplotdt, V_DBL, 1, VFLAG_NULL);
        param->tecplotdt = -1.0;

        BindVar(CPList, "tecplottime", &param->tecplottime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "tecmeshdstrain", &param->tecmeshdstrain, V_DBL, 1, VFLAG_NULL);
        param->tecmeshdstrain = -1.0;

        BindVar(CPList, "tecmeshstrain", &param->tecmeshstrain, V_DBL, 1, VFLAG_NULL);
        param->tecmeshstrain = 0.0;

        BindVar(CPList, "tecplotcounter", &param->tecplotcounter, V_INT, 1,
                VFLAG_NULL);
        param->tecplotcounter = 0;

        BindVar(CPList, "tecmeshcounter", &param->tecmeshcounter, V_INT, 1,
                VFLAG_NULL);
        param->tecmeshcounter = 0;

/*
 *      nodal velocity data
 */
        BindVar(CPList, "velfile", &param->velfile, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "velfilefreq", &param->velfilefreq, V_INT, 1,
                VFLAG_NULL);
        param->velfilefreq = 100;

        BindVar(CPList, "velfiledt", &param->velfiledt, V_DBL, 1, VFLAG_NULL);
        param->velfiledt = -1.0;

        BindVar(CPList, "velfiletime", &param->velfiletime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "velfilecounter", &param->velfilecounter, V_INT, 1,
                VFLAG_NULL);

/*
 *      nodal force data
 */
        BindVar(CPList, "writeForce", &param->writeForce, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeForceFreq", &param->writeForceFreq, V_INT, 1,
                VFLAG_NULL);
        param->writeForceFreq = 100;

        BindVar(CPList, "writeForceDT", &param->writeForceDT, V_DBL, 1,
                VFLAG_NULL);
        param->writeForceDT = -1.0;

        BindVar(CPList, "writeForceTime", &param->writeForceTime, V_DBL, 1,
                VFLAG_NULL);

        BindVar(CPList, "writeForceCounter", &param->writeForceCounter, V_INT,
                1, VFLAG_NULL);


/*
 *      Parameters related to creation of VisIt output files
 */
        BindVar(CPList, "writeVisit", &param->writeVisit, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitFreq", &param->writeVisitFreq,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitCounter", &param->writeVisitCounter,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitSegments", &param->writeVisitSegments,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitSegmentsAsText",
                &param->writeVisitSegmentsAsText, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitNodes", &param->writeVisitNodes,
                V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitNodesAsText",
                &param->writeVisitNodesAsText, V_INT, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitDT", &param->writeVisitDT,
                V_DBL, 1, VFLAG_NULL);

        BindVar(CPList, "writeVisitTime", &param->writeVisitTime,
                V_DBL, 1, VFLAG_NULL);

        param->writeVisitFreq = 100;
        param->writeVisitDT = -1.0;

/*
 *      X-window display
 */
        BindVar(CPList, "winDefaultsFile", param->winDefaultsFile, V_STRING, 1,
                VFLAG_NULL);
        strcpy(param->winDefaultsFile, "inputs/paradis.xdefaults");


/*
 *      3D-dislocation density data
 */
        BindVar(CPList, "savedensityspec", param->savedensityspec, V_INT, 3,
                VFLAG_NULL);


/*
 *      Miscellaneous parameters
 */
        BindVar(CPList, "Miscellaneous parameters", (void *)NULL,
                V_COMMENT, 1, VFLAG_NULL);

        BindVar(CPList, "enforceGlidePlanes", &param->enforceGlidePlanes,
                V_INT, 1, VFLAG_NULL);
        param->enforceGlidePlanes = 1;

        BindVar(CPList, "enableCrossSlip", &param->enableCrossSlip,
                V_INT, 1, VFLAG_NULL);
        param->enableCrossSlip = -1;

        BindVar(CPList, "TensionFactor", &param->TensionFactor, V_DBL, 1,
                VFLAG_NULL);
        param->TensionFactor = 1;

        BindVar(CPList, "elasticinteraction", &param->elasticinteraction,
                V_INT, 1, VFLAG_NULL);
        param->elasticinteraction = 1;
        
        BindVar(CPList, "meshNumX", &param->meshNumX, V_INT, 1,
                VFLAG_NULL);
        param->meshNumX = 11;

        BindVar(CPList, "meshNumY", &param->meshNumY, V_INT, 1,
                VFLAG_NULL);
        param->meshNumY = 11;

        BindVar(CPList, "meshNumZ", &param->meshNumZ, V_INT, 1,
                VFLAG_NULL);
        param->meshNumZ = 11;

        BindVar(CPList, "position", &param->position, V_DBL, 1,
                VFLAG_NULL);
        param->position = 0.0;

        BindVar(CPList, "loadTime", &param->loadTime, V_DBL, 1,
                VFLAG_NULL);
        param->loadTime = -1;

        BindVar(CPList, "correLength", &param->correLength, V_DBL, 1,
                VFLAG_NULL);
        param->correLength = 1;

        BindVar(CPList, "disRule", &param->disRule, V_INT, 1,
                VFLAG_NULL);
        param->disRule = 0;

        BindVar(CPList, "targetStress", &param->targetStress, V_DBL, 1,
                VFLAG_NULL);
        param->targetStress = 300;

        BindVar(CPList, "criticalEff", &param->criticalEff, V_DBL, 1,
                VFLAG_NULL);
        param->criticalEff = 5;
/*
 *      Grain parameters
 */
        BindVar(CPList, "peneTB", &param->peneTB,
                V_INT, 1, VFLAG_NULL);
        param->peneTB = 1;
        BindVar(CPList, "gbStrength", &param->gbStrength,
                V_DBL, 1, VFLAG_NULL);
        param->gbStrength = -1;

        BindVar(CPList, "stressFieldPenetrable", &param->stressFieldPenetrable, V_INT, 1, VFLAG_NULL);
        param->stressFieldPenetrable = 1;
            
        BindVar(CPList, "tbStrength", &param->tbStrength,
                V_DBL, 2, VFLAG_NULL);
        param->tbStrength[0] = -1.0;
        param->tbStrength[1] = -1.0;

        BindVar(CPList, "EsciagParams", &param->EsciagParams, V_DBL, 3, VFLAG_NULL);
        param->EsciagParams[0] = -1;
        param->EsciagParams[1] = 0.0;
        param->EsciagParams[2] = 0.0;

        BindVar(CPList, "grainSize", &param->grainSize,
                V_DBL, 1, VFLAG_NULL);
        param->grainSize = -1;
        
        BindVar(CPList, "a0", &param->a0, V_DBL, 4, VFLAG_NULL);
        param->a0[0] = -0.1;
        param->a0[1] = -0.1;
        param->a0[2] = -0.1;
        param->a0[3] = -0.1;

        BindVar(CPList, "fricStress", &param->fricStress,
                V_DBL, 1, VFLAG_NULL);
        param->fricStress = 0.0e0;

        BindVar(CPList, "preRann", &param->preRann, V_DBL, 1, VFLAG_NULL);
        param->preRann = -1.0;

        BindVar(CPList, "minDT", &param->minDT, V_DBL, 1, VFLAG_NULL);
        param->minDT = 1.0E-16;
        
        BindVar(CPList, "limitDisplacement", &param->limitDisplacement, V_INT, 1, VFLAG_NULL);
        param->limitDisplacement = 0;

        BindVar(CPList, "randomInternalForce", &param->randomInternalForce, V_INT, 1, VFLAG_NULL);
        param->randomInternalForce = 0;
        
        BindVar(CPList, "scale", &param->scale, V_DBL, 1, VFLAG_NULL);
        param->scale = 1.0;

#ifdef JOG_PINNING
        BindVar(CPList, "jogCoeff", &param->jogCoeff, V_DBL, 1, VFLAG_NULL);
        param->jogCoeff = 0.5;
#endif

#ifdef VELOCITY_STATISTICS
        BindVar(CPList, "accumulation", &param->accumulation, V_DBL, 8, VFLAG_NULL);
        param->accumulation[0] = 0.0;
        param->accumulation[1] = 0.0;
        param->accumulation[2] = 0.0;
        param->accumulation[3] = 0.0;
        param->accumulation[4] = 0.0;
        param->accumulation[5] = 0.0;
        param->accumulation[6] = 0.0;
        param->accumulation[7] = 0.0;
        
        BindVar(CPList, "disLength", &param->disLength, V_DBL, 3, VFLAG_NULL);
        param->disLength[0] = 0.0;
        param->disLength[1] = 0.0;
        param->disLength[2] = 0.0;

        BindVar(CPList, "cumulativeTime", &param->cumulativeTime, V_DBL, 1, VFLAG_NULL);
        param->cumulativeTime = 0.0;

        BindVar(CPList, "aveLength", &param->aveLength, V_DBL, 3, VFLAG_NULL);
        param->aveLength[0] = 0.0;
        param->aveLength[1] = 0.0;
        param->aveLength[2] = 0.0;
#endif

#ifdef _FEM
/*
 *      Add in some FEM specific stuff for M. Tang
 */
        BindVar(CPList, "BC_type", &param->BC_type, V_INT, 1, VFLAG_NULL);
        param->BC_type = 5;

        BindVar(CPList, "mesh_type", &param->mesh_type, V_INT, 1, VFLAG_NULL);
        param->mesh_type = 1;

        BindVar(CPList, "dirmax", &param->dirmax, V_INT, 1, VFLAG_NULL);
        param->dirmax = 2000;

        BindVar(CPList, "fem_nx", &param->fem_nx, V_INT, 1, VFLAG_NULL);
        param->fem_nx = 10;

        BindVar(CPList, "fem_ny", &param->fem_ny, V_INT, 1, VFLAG_NULL);
        param->fem_ny = 10;

        BindVar(CPList, "fem_nz", &param->fem_nz, V_INT, 1, VFLAG_NULL);
        param->fem_nz = 10;

        BindVar(CPList, "fem_radius", &param->fem_radius, V_DBL, 1, VFLAG_NULL);
        param->fem_radius = 1750.0;

        BindVar(CPList, "fem_height", &param->fem_height, V_DBL, 1, VFLAG_NULL);
        param->fem_height = 17500.0;

        BindVar(CPList, "fem_nr", &param->fem_nr, V_INT, 1, VFLAG_NULL);
        param->fem_nr = 5;

        BindVar(CPList, "fem_nh", &param->fem_nh, V_INT, 1, VFLAG_NULL);
        param->fem_nh = 10;

        BindVar(CPList, "fem_void_radius", &param->fem_void_radius, V_DBL,
                1, VFLAG_NULL);
        param->fem_void_radius = 500.0;

        BindVar(CPList, "fem_cube_edge_length", &param->fem_cube_edge_length,
                V_DBL, 1, VFLAG_NULL);
        param->fem_cube_edge_length = 1500.0;

        BindVar(CPList, "fem_numelm_arc", &param->fem_numelm_arc, V_INT,
                1, VFLAG_NULL);
        param->fem_numelm_arc = 5;

        BindVar(CPList, "fem_base_diagonal_half",
                &param->fem_base_diagonal_half, V_DBL, 1, VFLAG_NULL);
        param->fem_base_diagonal_half = 500.0;

        BindVar(CPList, "fem_grid_base", &param->fem_grid_base, V_INT,
                1, VFLAG_NULL);
        param->fem_grid_base = 3;

        BindVar(CPList, "fem_grid_vertical", &param->fem_grid_vertical,
                V_INT, 1, VFLAG_NULL);
        param->fem_grid_vertical = 3;

        BindVar(CPList, "fem_ageom_x", param->fem_ageom_x, V_DBL, 3,
                VFLAG_NULL);
        param->fem_ageom_x[0] = 1.0;
        param->fem_ageom_x[1] = 0.0;
        param->fem_ageom_x[2] = 0.0;

        BindVar(CPList, "fem_ageom_z", param->fem_ageom_z, V_DBL, 3,
                VFLAG_NULL);
        param->fem_ageom_z[0] = 0.0;
        param->fem_ageom_z[1] = 0.0;
        param->fem_ageom_z[2] = 1.0;
#endif

#ifdef PETSC
        BindVar(CPList, "yamlFile", param->yamlFile, V_STRING,
                1, VFLAG_NULL);
        strcpy(param->yamlFile,"options.yaml");
        BindVar(CPList, "optsFile", param->optsFile, V_STRING,
                1, VFLAG_NULL);
        strcpy(param->optsFile,"options.inf");

        BindVar(CPList, "velocityKernel", &param->velocityKernel, V_INT, 1, VFLAG_NULL);
        param->velocityKernel = 0;
#else
        BindVar(CPList, "velocityKernel", &param->velocityKernel, V_INT, 1, VFLAG_NULL);
        param->velocityKernel = 0;
#endif
        return;
}


void DataParamInit(Param_t *param, ParamList_t *DPList)
{

/*
 *      Note: Parameters need only be initialized if their
 *      default values are non-zero.
 */
        DPList->paramCnt = 0;
        DPList->varList = (VarData_t *)NULL;

        BindVar(DPList, "dataFileVersion", &param->dataFileVersion, V_INT,
                1, VFLAG_NULL);
        param->dataFileVersion = NODEDATA_FILE_VERSION;

        BindVar(DPList, "numFileSegments", &param->numFileSegments, V_INT,
                1, VFLAG_NULL);
        param->numFileSegments = 1;

        BindVar(DPList, "minCoordinates", param->minCoordinates, V_DBL, 3,
                VFLAG_NULL);

        BindVar(DPList, "maxCoordinates", param->maxCoordinates, V_DBL, 3,
                VFLAG_NULL);

        BindVar(DPList, "nodeCount", &param->nodeCount, V_INT, 1, VFLAG_NULL);

        BindVar(DPList, "dataDecompType", &param->dataDecompType, V_INT, 1,
                VFLAG_NULL);
        param->dataDecompType = 2;

        BindVar(DPList, "dataDecompGeometry", &param->dataDecompGeometry,
                V_INT, 3, VFLAG_NULL);
        param->dataDecompGeometry[0] = 1;
        param->dataDecompGeometry[1] = 1;
        param->dataDecompGeometry[2] = 1;
        
        BindVar(DPList, "maxSeg", &param->maxSeg, V_DBL, 1, VFLAG_NULL);
        param->maxSeg = 200;
        
        BindVar(DPList, "maxSchmid", &param->maxSchmid, V_DBL, 1, VFLAG_NULL);
        param->maxSchmid = -1;
        
        BindVar(DPList, "minSeg", &param->minSeg, V_DBL, 1, VFLAG_NULL);
        param->minSeg = 100;
        BindVar(DPList, "polycrystal", &param->polycrystal, V_INT, 1, VFLAG_NULL);
        param->polycrystal = 1;
            
        BindVar(DPList, "nGrain", &param->nGrain, V_INT, 1, VFLAG_NULL);
        param->nGrain = 2;
            
        BindVar(DPList, "twinRatio", &param->twinRatio, V_DBL, 1, VFLAG_NULL);
        param->twinRatio = 0.5;
        
        BindVar(DPList, "gradientStructure", param->gradientStructure, V_STRING, 1, VFLAG_NULL);
        snprintf(param->gradientStructure,sizeof(param->gradientStructure),"0");

        BindVar(DPList, "xBoundType", &param->xBoundType, V_INT, 1, VFLAG_NULL);
        param->xBoundType = (BoundType_t)0;
    
        BindVar(DPList, "yBoundType", &param->yBoundType, V_INT, 1, VFLAG_NULL);
        param->yBoundType = (BoundType_t)0;
    
        BindVar(DPList, "zBoundType", &param->zBoundType, V_INT, 1, VFLAG_NULL);
        param->zBoundType = (BoundType_t)0;
    
        BindVar(DPList, "xBoundMin", &param->xBoundMin, V_DBL, 1, VFLAG_NULL);
        BindVar(DPList, "yBoundMin", &param->yBoundMin, V_DBL, 1, VFLAG_NULL);
        BindVar(DPList, "zBoundMin", &param->zBoundMin, V_DBL, 1, VFLAG_NULL);
    
        BindVar(DPList, "xBoundMax", &param->xBoundMax, V_DBL, 1, VFLAG_NULL);
        BindVar(DPList, "yBoundMax", &param->yBoundMax, V_DBL, 1, VFLAG_NULL);
        BindVar(DPList, "zBoundMax", &param->zBoundMax, V_DBL, 1, VFLAG_NULL);
    
        BindVar(DPList, "matrixOrientationX", &param->matrixOrientation[0], V_DBL, 3, VFLAG_NULL);
        BindVar(DPList, "matrixOrientationY", &param->matrixOrientation[1], V_DBL, 3, VFLAG_NULL);
        BindVar(DPList, "matrixOrientationZ", &param->matrixOrientation[2], V_DBL, 3, VFLAG_NULL);
    
        BindVar(DPList, "TBNormal", &param->TBNormal, V_DBL, 3,
                VFLAG_NULL);
        param->TBNormal[0] = 1.0;  
        param->TBNormal[1] = 1.0;  
        param->TBNormal[2] = 1.0;  
    
        BindVar(DPList, "TBNormalCryst", &param->TBNormalCryst, V_DBL, 3,
                VFLAG_NULL);
        param->TBNormalCryst[0] = 1.0;  
        param->TBNormalCryst[1] = 1.0;  
        param->TBNormalCryst[2] = 1.0;  
    
        param->preDecreaseDT = 0;
        param->markPreDTDisable = 0;
        param->finish = 0;

        return;
}
