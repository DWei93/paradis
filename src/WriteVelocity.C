#include <sys/types.h>
#include <sys/stat.h>
#include "Home.h"
#include "Util.h"
/*---------------------------------------------------------------------------
 *
 *      Function:    WriteVelocity
 *      Description: Write the velocity of all active nodes to a file.
 *
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *
 *-------------------------------------------------------------------------*/
void WriteVelocity(Home_t *home, char *baseFileName, int ioGroup,
                   int firstInGroup, int writePrologue, int writeEpilogue)
{
        int      i, newNodeKeyPtr, j;
        real8    vx, vy, vz, burgMag, velocity, length;
	real8    velDir[3], fNor, gForce[3], armf[3];
	real8	 lDir[3], nDir[3], gDir[3], dot;
        char     fileName[256];
        Node_t   *node, *nbr;
        Param_t  *param = home->param;
        FILE     *fp;
        struct stat statbuf;

/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s", DIR_VELOCITY,
                     baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d", DIR_VELOCITY,
                     baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS

/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (firstInGroup) {
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteVelocity: Open error %d on %s\n", errno, fileName);
            }

	    if(writePrologue){
		fprintf(fp, "variables = X,Y,Z,gF_X,gF_Y,gF_z,vel,length,sign\n");	
	    }
        } else {
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WriteVelocity: Open error %d on %s\n", errno, fileName);
            }
        }

/*
 *      Write velocity data for all nodes in this domain
 */
        burgMag = param->burgMag;
        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

	    length = 0;
	    VECTOR_ZERO(gForce);
	    for(j=0; j<node->numNbrs; j++){
                nbr = GetNeighborNode(home, node, j);
		if(nbr == (Node_t *)NULL)continue;
		length += (NodeDistance(home, node, nbr)/2.0);
		
	        lDir[X] = nbr->x - node->x;
	        lDir[Y] = nbr->y - node->y;
	        lDir[Z] = nbr->z - node->z;
	        ZImage(param, lDir, lDir+1, lDir+2);

	        if(Normal(lDir) < 1e-5) continue;
	        NormalizeVec(lDir);
	        nDir[X] = node->nx[j];
	        nDir[Y] = node->ny[j];
	        nDir[Z] = node->nz[j];
	        cross(lDir, nDir, gDir);

	        if(Normal(gDir) < 1e-5) continue;
	        NormalizeVec(gDir);

	        armf[X] = node->armfx[j];
	        armf[Y] = node->armfy[j];
	        armf[Z] = node->armfz[j];
    
	        dot = DotProduct(armf, gDir);
	        gForce[X] += (dot * gDir[X]);
	        gForce[Y] += (dot * gDir[Y]);
	        gForce[Z] += (dot * gDir[Z]);
	    }
	    
            vx =  node->vX;
            vy =  node->vY;
            vz =  node->vZ;
	    velocity = sqrt(vx*vx +vy*vy + vz*vz);

	    velDir[X]= node->vX;
	    velDir[Y]= node->vY;
	    velDir[Z]= node->vZ;
	    NormalizeVec(velDir);

	    fNor = sqrt(node->fX * node->fX +
	    		node->fY * node->fY +
	    		node->fZ * node->fZ);

            fprintf(fp, "%7.1f %7.1f %7.1f %10.3e %10.3e %10.3e %10.3e %7.1f %2d \n", 
		    node->x, node->y, node->z, 
		    fNor*velDir[X], fNor*velDir[Y],fNor*velDir[Z],
		    velocity, length, node->sgnv);
        }

        fclose(fp);

        return;
}
