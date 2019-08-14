#ifndef	_GRAIN_BOUNDARY_H
#define	_GRAIN_BOUNDARY_H

#include "Typedefs.h"
#include "Home.h"
#include "Grain.h"
#include "Util.h" 
#include "Tag.h" 
#include "ParadisProto.h"
#include <vector>

/** Find the intersection between nodal path and grain (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of the node (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathInGrain(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);

/** Find the intersection between nodal path and face (if has) and return the new place of intersection point.
 * \param[in] place the location of the node.
 * \param[in] pos1 the old position of the node, it has to be in the place.
 * \param[in,out] pos2 the current position of the node, moved to intersection if it was outside of the place. 
 * \param[out] newPlace the new place of intersection (if has) 
 * \return true if the place is changed, otherwise false */
bool ClipPathOnFace(Home_t *home, int *place, real8 *pos1, real8 *pos2, int *newPlace);

#endif
