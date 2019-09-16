#ifndef	_GRAIN_BOUNDARY_H
#define	_GRAIN_BOUNDARY_H

#include "Typedefs.h"
#include "Home.h"
#include "Grain.h"
#include "Util.h" 
#include "Tag.h" 
#include "ParadisProto.h"
#include <vector>

/** Assign nodal location
 *  \param[in] node the node
 *  \param[in] newPlace the new location
 *  \param[in] globalOp true if add global operation list */
inline void AssignLocation(Home_t *home, Node_t *node, int *newPlace, bool globalOp) 
{
    node->place[0]=newPlace[0]; node->place[1]=newPlace[1];
    node->place[2]=newPlace[2]; node->place[3]=newPlace[3];
    if(globalOp)AddOp(home, CHANGE_LOCATION, node->myTag.domainID, node->myTag.index,   
                      newPlace[0], newPlace[1], newPlace[2], newPlace[3],
                      0,0,0, 0,0,0, 0,0,0);
	return;
}

/** Assign nodal constraint
 *  \param[in] node the node
 *  \param[in] constraint the new constraint
 *  \param[in] globalOp true if add global operation list */
inline void AssignConstraint(Home_t *home, Node_t *node, int constraint, bool globalOp) 
{
    node->constraint = constraint;
    if(globalOp)AddOp(home, CHANGE_LOCATION, node->myTag.domainID, node->myTag.index,   
                      constraint, 0, 0, 0,
                      0,0,0, 0,0,0, 0,0,0);
    return;
}

/** The main clipping function for clipping all nodes in parallel */
void Clipping(Home_t *home);

#endif
