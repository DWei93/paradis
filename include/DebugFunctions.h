#ifndef _DEBUGFUNCTIONS_H
#define _DEBUGFUNCTIONS_H
/****************************************************************************
 *
 *      DebugFunctions.h  Contains prototypes for various debug-only functions
 *
 ***************************************************************************/

#ifdef DEBUG_PLASTIC_STRAIN
struct dp_t{
    Node_t *node;
    BurgInfo_t *bd;
    int numNbrs;
    int burgIDList[4];
    int planeIDList[4];
    real8 len;
    real8 dis;
    real8 nodePst;                
};
real8 DumpDeltaPlasticStrain(Home_t *home, int pstType = 0, bool printFile = false, int flag = 0); // plastic strain default
real8 DeltaStress(Home_t *home, int pstType, real8 &dpstn);
#endif

void CheckSegLengths(Home_t *home, const char *msg);
void CheckForNANS(Home_t *home, real8 flag);
void CheckForEmptySimulation(Home_t *home);
void CheckForUndefinedPlanes(Home_t *home, const char *msg);
void CheckAll(Home_t *home, const char *format, ...);
void CheckLocal(Home_t *home, int flag);
void CheckNodeNANS(Node_t *node);
void DebugGB(Home_t *home, int flag);

#endif
