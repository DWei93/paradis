/***************************************************************************
 *
 *	OpList.h	Define the struct that holds all data needed
 *			to handle a cross-domain topology change.
 *
 **************************************************************************/


#ifndef _OPLIST_H
#define _OPLIST_H

#include "Typedefs.h"
#include "Node.h"
#include "Tag.h"

#define OpBlock_Count 500

struct _operate {
	OpType_t	type;
	int			dom1;
	int			idx1;
	int			dom2;
	int			idx2;
	int			dom3;
	int			idx3;
	real8		bx;
	real8		by;
	real8		bz;
	real8		x;
	real8		y;
	real8		z;
	real8		nx;
	real8		ny;
	real8		nz;
};

/*
 *      Prototype functions related to managing the remote operation list
 */
void AddOp(Home_t *home, OpType_t type, int dom1, int idx1,
        int dom2, int idx2, int dom3, int idx3,
        real8 bx, real8 by, real8 bz, real8 x, real8 y, real8 z,
        real8 nx, real8 ny, real8 nz);
inline void UpdateVelocity(Home_t *home, Node_t *node, real8 *vel){
	node->vX = vel[0];
	node->vY = vel[1];
	node->vZ = vel[2];
	AddOp(home, UPDATE_VELOCITY, node->myTag.domainID,node->myTag.index,
			-1,-1,-1,-1,0.0,0.0,0.0,0.0,0.0,0.0,
			vel[0],vel[1],vel[2]);
	return;
}
inline void UpdateGrain(Home_t *home, Node_t *node, int inGrain, int onGB){
	node->inGrain=inGrain; node->onGB =onGB;
	AddOp(home, UPDATE_GRAIN, node->myTag.domainID,node->myTag.index,
			inGrain,onGB,-1,-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	return;
}
inline void UpdateConstraint(Home_t *home, Node_t *node, int constraint){
	node->constraint = constraint;
	AddOp(home, UPDATE_CONSTRAINT, node->myTag.domainID,node->myTag.index,
			constraint,-1,-1,-1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
	return;
}

void ClearOpList(Home_t *home);
void ExtendOpList(Home_t *home);
void FreeOpList(Home_t *home);
void InitOpList(Home_t *home);
void PrintOpList(Home_t *home);

#endif /* _OPLIST_H */
