#ifndef _DEBUGFUNCTIONS_H
#define _DEBUGFUNCTIONS_H
/****************************************************************************
 *
 *      DebugFunctions.h  Contains prototypes for various debug-only functions
 *
 ***************************************************************************/

void CheckSegLengths(Home_t *home, const char *msg);
void CheckForNANS(Home_t *home, real8 flag);
void CheckForEmptySimulation(Home_t *home);
void CheckForUndefinedPlanes(Home_t *home, const char *msg);
void CheckAll(Home_t *home, const char *format, ...);
void CheckLocal(Home_t *home, int flag);
void CheckNodeNANS(Node_t *node);
void DebugGB(Home_t *home, int flag);
#endif
