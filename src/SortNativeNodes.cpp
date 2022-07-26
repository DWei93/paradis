/**************************************************************************
 *
 *  Function    : SortNativeNodes 
 *  Description : Queue each native node onto the cell it currently falls into.
 *
 *    7/31/01  Replaced sideX with minSideX, maxSideX, etc     t.p.
 *
 **************************************************************************/

#include "Home.h"
#include "Node.h"
#include "Cell.h"
#include "Util.h"

void AssignNodeToCell(Home_t *home, Node_t *node)
{
        Param_t *param;
        int iCell, jCell, kCell, cellIdx;
        real8 probXmin, probYmin, probZmin;
        real8 cellXsize, cellYsize, cellZsize;
        Cell_t *cell;

        if (node == (Node_t *)NULL) return;

/*
 *      set the lower limit of the base area (excluding possible periodic cells)
 *      and the size of each cell
 */
        param = home->param;
        probXmin = param->minSideX;
        probYmin = param->minSideY;
        probZmin = param->minSideZ;

        cellXsize = (param->maxSideX - param->minSideX) / param->nXcells;
        cellYsize = (param->maxSideY - param->minSideY) / param->nYcells;
        cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells;

        real8 x = node->x;
        real8 y = node->y;
        real8 z = node->z;
/*
 *      put the node on its proper cell. If the index exceeds this domains
 *      range of native cells, put in the nearest native cell
 */
        if(param->xBoundType == Periodic)if(x>=param->maxSideX)x=probXmin;
        iCell = (int)((x - probXmin) / cellXsize);
        if (iCell < param->iCellNatMin) iCell = param->iCellNatMin;
        if (iCell > param->iCellNatMax) iCell = param->iCellNatMax;
        iCell++;   /* compensate for periodic cells */

        if(param->yBoundType == Periodic)if(y>=param->maxSideY)y=probYmin;
        jCell = (int)((y - probYmin) / cellYsize);
        if (jCell < param->jCellNatMin) jCell = param->jCellNatMin;
        if (jCell > param->jCellNatMax) jCell = param->jCellNatMax;
        jCell++;   /* compensate for periodic cells */

        if(param->zBoundType == Periodic)if(z>=param->maxSideZ)z=probZmin;
        kCell = (int)((z - probZmin) / cellZsize);
        if (kCell < param->kCellNatMin) kCell = param->kCellNatMin;
        if (kCell > param->kCellNatMax) kCell = param->kCellNatMax;
        kCell++;   /* compensate for periodic cells */

        cellIdx = EncodeCellIdx (home, iCell, jCell, kCell);
        cell = home->cellKeys[cellIdx];
        node->nextInCell = cell->nodeQ;
        cell->nodeQ = node;
        cell->nodeCount++;

        node->cellIdx = cellIdx;

        return;
}


void SortNativeNodes (Home_t *home)
{
   Param_t *param;
   int i, iCell, jCell, kCell, cellIdx ;
   real8 probXmin, probYmin, probZmin ;
   real8 cellXsize, cellYsize, cellZsize ;
   Cell_t *cell ;
   Node_t *node ;

   real8 x,y,z;

   TimerStart(home, SORT_NATIVE_NODES);

/* set the lower limit of the base area (excluding possible periodic cells)
 * and the size of each cell
 */

   param = home->param;
   probXmin = param->minSideX ;
   probYmin = param->minSideY ;
   probZmin = param->minSideZ ;

   cellXsize = (param->maxSideX - param->minSideX) / param->nXcells ;
   cellYsize = (param->maxSideY - param->minSideY) / param->nYcells ;
   cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells ;

/* Loop through cells and reinitialize their node 
   queues to empty.  */

   for (i = 0 ; i < home->cellCount ; i++) {

      cell = home->cellKeys[home->cellList[i]] ;
      cell->nodeQ = 0 ;
      cell->nodeCount = 0 ;
   }

/* Loop thru active nodes, putting them in their proper cell. If the
 * index exceeds this domains range of native cells, put in the nearest 
 * native cell
 */

   for (i = 0 ; i < home->newNodeKeyPtr ; i++) {

        node = home->nodeKeys[i] ;
        if (!node) continue ;
        
        x=node->x; y=node->y; z=node->z;
        
        if(param->xBoundType == Periodic)if(x>=param->maxSideX)x=probXmin;
        iCell = (int)((x - probXmin) / cellXsize);
        if (iCell < param->iCellNatMin) iCell = param->iCellNatMin;
        if (iCell > param->iCellNatMax) iCell = param->iCellNatMax;
        iCell++;   /* compensate for periodic cells */
        
        if(param->yBoundType == Periodic)if(y>=param->maxSideY)y=probYmin;
        jCell = (int)((y - probYmin) / cellYsize);
        if (jCell < param->jCellNatMin) jCell = param->jCellNatMin;
        if (jCell > param->jCellNatMax) jCell = param->jCellNatMax;
        jCell++;   /* compensate for periodic cells */
        
        if(param->zBoundType == Periodic)if(z>=param->maxSideZ)z=probZmin;
        kCell = (int)((z - probZmin) / cellZsize);
        if (kCell < param->kCellNatMin) kCell = param->kCellNatMin;
        if (kCell > param->kCellNatMax) kCell = param->kCellNatMax;
        kCell++;   /* compensate for periodic cells */
        
        cellIdx = EncodeCellIdx (home, iCell, jCell, kCell) ;
        cell = home->cellKeys[cellIdx] ;
        node->nextInCell = cell->nodeQ ;
        cell->nodeQ = node ;
        cell->nodeCount++ ;
        
        node->cellIdx = cellIdx ;
        node->native = 1 ;

   }

   TimerStop(home, SORT_NATIVE_NODES);

   return;
}
