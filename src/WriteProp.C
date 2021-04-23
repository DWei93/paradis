/****************************************************************************
 *
 *      Function:     WriteProp
 *      Description:  Depending on the input parameter, write a particular
 *                    time-dependent property out to that property's
 *                    diagnostic file, along with a timestamp.
 *
 *      NOTE: The contents of the file written below containing the
 *      dislocation density have changed over time.  The contents of
 *      the various versions are defined below.
 *
 *      Version 0 contents:
 *        - strain
 *        - dislocation density
 *        - deleted segment length
 *        - average velocity
 *        - std deviation of dislocation velocities
 *
 *      Version 1 contents: After the velocity std deviation, the
 *      following items were added:
 *        - version number
 *        - a set of values indicating the dislocation density for
 *          segments of specific groupings of burgers vectors.  These
 *          groupings differ for BCC and FCC mobility.
 *
 *          BCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 1] [-1-1-1]
 *            1         [-1 1 1] [ 1-1-1]
 *            2         [ 1-1 1] [-1 1-1]
 *            3         [ 1 1-1] [-1-1 1]
 *            4         [ 1 0 0] [-1 0 0]
 *                      [ 0 1 0] [ 0-1 0]
 *                      [ 0 0 1] [ 0 0-1]
 *
 *          FCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 0] [-1-1 0]
 *            1         [-1 1 0] [ 1-1 0]
 *            2         [ 1 0 1] [-1 0-1]
 *            3         [-1 0 1] [ 1 0-1]
 *            4         [ 0 1 1] [ 0-1-1]
 *            5         [ 0-1 1] [ 0 1-1]
 *            6         all others
 *
 *      Version 2 contents:
 *        - plastic strain added as first column of data
 *
 ****************************************************************************/

#include "Home.h"
#include "WriteProp.h"
#include "Util.h"
#include "Mobility.h"

#include <vector>
#include <string>

using namespace std;

#define DENSITY_FILE_VERSION 2

void WriteProp(Home_t *home, int property)
{
        int     i, numItems, tmpOffset, numBurgVectors;
        real8   al, am, an, amag, sigijk, pstnijk, dpstnijk, bl, bm, bn;
        real8   ptstnijk;
        real8   *localDensityVals  = (real8 *)NULL;
        real8   *globalDensityVals = (real8 *)NULL;
        real8   totDensityChange[14];
	    real8	TC;
        real8   dpstn[6], effPlaRate, sign;
        char    fileName[256], var[256];
        FILE    *fp;
        Param_t *param;
        int     loadDirID;

        param = home->param;
#ifdef VELOCITY_STATISTICS
        real8   nav_a = 0.0, nav_ia = 0.0;
        int     dissStartID = 3*param->nVelBins;
        int     dataLen = 3*param->nVelBins+3*param->nDissBins;
        real8   *locVals = (real8 *)calloc(dataLen, sizeof(real8));
        real8   *gloVals = (real8 *)calloc(dataLen, sizeof(real8));
#endif
		for (i = 0; i < 6; i++) {
            dpstn[i] = param->delpStrain[i];
        }
		
        effPlaRate = sqrt((dpstn[0]*dpstn[0]+ 
                           dpstn[1]*dpstn[1]+ 
                           dpstn[2]*dpstn[2]+ 
                       2.0*dpstn[3]*dpstn[3]+ 
                       2.0*dpstn[4]*dpstn[4]+ 
                       2.0*dpstn[5]*dpstn[5])*2.0/3.0)/param->lastDT;

	    TC = Sign(param->eRate);
        if(param->loadType!=6){
            al = param->edotdir[0];
            am = param->edotdir[1];
            an = param->edotdir[2];
            
            amag = sqrt(al*al+am*am+an*an);
            
            al /= amag;
            am /= amag;
            an /= amag;
            bl=al; bm=am; bn=an;
        }else{
            al = param->loadDir[0];
            am = param->loadDir[1];
            an = param->loadDir[2];
            bl = param->loadDir[3];
            bm = param->loadDir[4];
            bn = param->loadDir[5];

        }
        
        switch(property) {
        case DENSITY:
        
#ifdef PARALLEL
            numItems = param->numBurgGroups+1;
            tmpOffset = 1;
#ifdef _FEM
            numItems += 1;
            tmpOffset += 1;
#endif

            localDensityVals = (real8*)malloc((numItems)*sizeof(real8));
            globalDensityVals = (real8*)malloc((numItems)*sizeof(real8));

            for (i = 0; i < numItems; i++) {
                globalDensityVals[i] = 0.0;
            }

            localDensityVals[0] = param->delSegLength;
#ifdef _FEM
            localDensityVals[1] = param->fem_delSegLength;
#endif

            for (i = 0; i < param->numBurgGroups; i++) {
                localDensityVals[i+tmpOffset] = param->partialDisloDensity[i];
            }
        
            MPI_Allreduce(localDensityVals, globalDensityVals, numItems,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
            param->delSegLength = globalDensityVals[0];
#ifdef _FEM
            param->fem_delSegLength = globalDensityVals[1];
            param->fem_delSegLength *= param->burgVolFactor;
#endif

            for (i = 0; i < param->numBurgGroups; i++) {
                param->partialDisloDensity[i] = globalDensityVals[i+tmpOffset];
            }
        
            free(localDensityVals);
            free(globalDensityVals);
#endif

/*
 *          convert to density and write results to the file
 */
            param->delSegLength *= param->burgVolFactor;
        
            if (home->myDomain == 0) {
                snprintf(fileName, sizeof(fileName), "%s/density",
                         DIR_PROPERTIES);
                if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                    Fatal("can not open file %s", fileName);
                }

        
                pstnijk =  param->totpStn[0]*al*bl        +
                           param->totpStn[1]*am*bm        +
                           param->totpStn[2]*an*bn        +
                           param->totpStn[3]*(am*bn+an*bm)+
                           param->totpStn[4]*(an*bl+al*bn)+
                           param->totpStn[5]*(al*bm+bl*am);

                ptstnijk = param->totTwinpStn[0]*al*bl        +
                           param->totTwinpStn[1]*am*bm        +
                           param->totTwinpStn[2]*an*bn        +
                           param->totTwinpStn[3]*(am*bn+an*bm)+
                           param->totTwinpStn[4]*(an*bl+al*bn)+
                           param->totTwinpStn[5]*(al*bm+bl*am);
        
/*
 *              First print the standard stuff that's common
 *              regardless of which mobility (BCC, FCC, etc)
 *              is in use.
 */
                fprintf(fp, "%e %e %e %e %e %e %e",
                        pstnijk,param->eRate*param->timeNow, param->disloDensity,
                        param->delSegLength, param->vAverage,
                        param->vStDev,ptstnijk);

/*
 *              ****************************************************
 *              ***                IMPORTANT!                    ***
 *              ***   If you change the contents of the density  ***
 *              ***   file, increment the definition for         ***
 *              ***   DENSITY_FILE_VERSION above and document    ***
 *              ***   the new contents.  DO NOT REMOVE the       ***
 *              ***   description of the previous version!       ***
 *              ****************************************************
 *                  
 *              Now append to the line the file version number
 *              and the dislocation density for each of the
 *              groupings of burgers vectors appropriate to the
 *              mobility type.  Should be kept in sync with
 *              DeltaPlasticStrain().
 */
                fprintf(fp, " %d", DENSITY_FILE_VERSION);

                for (i = 0; i < param->numBurgGroups; i++) {
                    fprintf(fp, " %e", param->partialDisloDensity[i]);
                }

                fprintf(fp, "\n");
                fclose(fp);

#ifdef _FEM
                snprintf(fileName, sizeof(fileName), "%s/fem_density",
                         DIR_PROPERTIES);
                if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                    Fatal("can not open file %s", fileName);
                }
                fprintf(fp, "%d %e %e\n", home->cycle,
                        param->eRate*param->timeNow, param->fem_delSegLength);
                fclose(fp);
#endif
            }
        
/*
 *          Reinitialize accumulated length of deleted segments
 *          after writing
 */
            param->delSegLength = 0.0;

#ifdef _FEM
            param->fem_delSegLength = 0.0;
#endif
            break;
        
        case EPS:

#ifdef VELOCITY_STATISTICS
#ifdef PARALLEL
            
            for(i=0; i<dataLen; i++)gloVals[i] = 0.0;

            for(i=0; i<param->nVelBins; i++){
                locVals[3*i] = param->velTable[i][1];
                locVals[3*i+1] = param->velTable[i][2];
                locVals[3*i+2] = param->velTable[i][3];
            }
            for(i=0; i<param->nDissBins; i++){
                locVals[dissStartID+3*i] = param->dissTable[i][1];
                locVals[dissStartID+3*i+1] = param->dissTable[i][2];
                locVals[dissStartID+3*i+2] = param->dissTable[i][3];
            }

	        MPI_Reduce(locVals, gloVals, dataLen, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if(home->myDomain == 0){
                for(i=0; i<param->nVelBins; i++){
                    param->velTable[i][1] = gloVals[3*i];
                    param->velTable[i][2] = gloVals[3*i+1];
                    param->velTable[i][3] = gloVals[3*i+2];
                }
                for(i=0; i<param->nDissBins; i++){
                    param->dissTable[i][1] = gloVals[dissStartID+3*i];
                    param->dissTable[i][2] = gloVals[dissStartID+3*i+1];
                    param->dissTable[i][3] = gloVals[dissStartID+3*i+2];
                }
            }else{
                for(i=0; i<param->nVelBins; i++){
                    param->velTable[i][1] = 0.0;
                    param->velTable[i][2] = 0.0;
                    param->velTable[i][3] = 0.0;
                }
                for(i=0; i<param->nDissBins; i++){
                    param->dissTable[i][1] = 0.0;
                    param->dissTable[i][2] = 0.0;
                    param->dissTable[i][3] = 0.0;
                }
            }

            free(locVals);
            free(gloVals);
#endif
#endif

            if (home->myDomain == 0) {
                sigijk = param->appliedStress[0]*al*bl        +
                         param->appliedStress[1]*am*bm        +
                         param->appliedStress[2]*an*bn        +
                         param->appliedStress[3]*(am*bn+an*bm)+
                         param->appliedStress[4]*(an*bl+al*bn)+ 
                         param->appliedStress[5]*(al*bm+bl*am);
        
                pstnijk =  param->totpStn[0]*al*bl        +
                           param->totpStn[1]*am*bm        +
                           param->totpStn[2]*an*bn        +
                           param->totpStn[3]*(am*bn+an*bm)+
                           param->totpStn[4]*(an*bl+al*bn)+
                           param->totpStn[5]*(al*bm+bl*am);

                ptstnijk = param->totTwinpStn[0]*al*bl        +
                           param->totTwinpStn[1]*am*bm        +
                           param->totTwinpStn[2]*an*bn        +
                           param->totTwinpStn[3]*(am*bn+an*bm)+
                           param->totTwinpStn[4]*(an*bl+al*bn)+
                           param->totTwinpStn[5]*(al*bm+bl*am);
        
                snprintf(fileName, sizeof(fileName), "%s/time_Plastic_strain",
                         DIR_PROPERTIES);
                if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                    Fatal("can not open file %s", fileName);
                }

                fprintf(fp,"%e %e %e\n", param->timeNow, pstnijk, ptstnijk);
                fclose(fp);
        
                if ((param->loadType > 0)) {
        
                    snprintf(fileName, sizeof(fileName),
                             "%s/stress_Plastic_strain", DIR_PROPERTIES);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

                    fprintf(fp, "%e %e %e %e %e %e %e\n", pstnijk, param->appliedStress[0],
                            param->appliedStress[1], param->appliedStress[2], param->appliedStress[3],
                            param->appliedStress[4], param->appliedStress[5]);
                    fclose(fp);
#ifdef VELOCITY_STATISTICS
/*
 *                  dump the velocity problebility table
 */
                    nav_a = 0.0; nav_ia = 0.0;
                    for(i=0; i<param->nVelBins; i++){
                        param->velTable[i][1] /= param->cumulativeTime;
                        param->velTable[i][2] /= param->cumulativeTime;
                        param->velTable[i][3] /= param->cumulativeTime;
                        nav_a += param->velTable[i][1];
                        nav_a += param->velTable[i][2];
                        nav_ia += param->velTable[i][3];
                    }

                    snprintf(fileName, sizeof(fileName), "%s/unsigned-velocity-%d.plt", DIR_PROPERTIES, home->cycle);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

					fprintf(fp, "variables = velocity, width, p_unsigned_velocity_active, p_velocity_inactive\n"); 
		            fprintf(fp, "zone T = \"%d\", i = %d, j = 1, k = 1, SOLUTIONTIME = %d, F=POINT\n", 
                                home->cycle, param->nVelBins, home->cycle); 
                    fprintf(fp, "AUXDATA cumulativeTime = \"%f\"\n", param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_v_squal_active = \"%e\"\n", param->accumulation[0]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_v_squal_active = \"%e\"\n", param->accumulation[1]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA v_squal_active = \"%e\"\n", param->accumulation[2]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA vsqual_inactive = \"%e\"\n", param->accumulation[3]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_length_active = \"%e\"\n", param->aveLength[0]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_length_active = \"%e\"\n", param->aveLength[1]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA length_inactive = \"%e\"\n", param->aveLength[2]/param->cumulativeTime);
                    
                    i=0;
                    fprintf(fp, "%e %e %e %e\n", param->velTable[i][0], param->velMin,
                            (param->velTable[i][1]+param->velTable[i][2])/
                            (param->velMin*nav_a), param->velTable[i][3]/(param->velMin*nav_ia));

                    for(i=1; i<param->nVelBins; i++){
                        fprintf(fp, "%e %e %e %e\n", param->velTable[i][0], param->velMin*(param->baseVel-1)*pow(param->baseVel, ((real8)i)-1.0), 
                                (param->velTable[i][1]+param->velTable[i][2])/
                                (param->velMin*pow(param->baseVel,(real8)i)*(param->baseVel-1)*nav_a), 
                                param->velTable[i][3]/(param->velMin*pow(param->baseVel,(real8)i)*(param->baseVel-1)*nav_ia));
                    }

                    fclose(fp);                   

                    snprintf(fileName, sizeof(fileName), "%s/signed-velocity-%d.plt", DIR_PROPERTIES, home->cycle);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

					fprintf(fp, "variables = velocity, width, p_signed_velocity_active\n"); 
		            fprintf(fp, "zone T = \"%d\", i = %d, j = 1, k = 1, SOLUTIONTIME = %d, F=POINT\n", 
                                home->cycle, 2*param->nVelBins, home->cycle); 
                    fprintf(fp, "AUXDATA cumulativeTime = \"%f\"\n", param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_v_squal_active = \"%e\"\n", param->accumulation[0]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_v_squal_active = \"%e\"\n", param->accumulation[1]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA v_squal_active = \"%e\"\n", param->accumulation[2]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA vsqual_inactive = \"%e\"\n", param->accumulation[3]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_length_active = \"%e\"\n", param->aveLength[0]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_length_active = \"%e\"\n", param->aveLength[1]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA length_inactive = \"%e\"\n", param->aveLength[2]/param->cumulativeTime);
                    
                    for(i=param->nVelBins-1; i>0; i--){
                        fprintf(fp, "%e %e %e\n", -param->velTable[i][0], param->velMin*(param->baseVel-1)*pow(param->baseVel, ((real8)i)-1.0),
                                param->velTable[i][1]/
                                (param->velMin*pow(param->baseVel,(real8)i)*(param->baseVel-1)*nav_a)); 
                    }
                    i=0;
                    fprintf(fp, "%e %e %e\n", -param->velTable[i][0], param->velMin,
                            param->velTable[i][1]/
                            (param->velMin*nav_a));
                    fprintf(fp, "%e %e %e\n", param->velTable[i][0], param->velMin,
                            param->velTable[i][2]/
                            (param->velMin*nav_a));

                    for(i=1; i<param->nVelBins; i++){
                        fprintf(fp, "%e %e %e\n", param->velTable[i][0], param->velMin*(param->baseVel-1)*pow(param->baseVel, ((real8)i)-1.0),
                                param->velTable[i][2]/
                                (param->velMin*pow(param->baseVel,(real8)i)*(param->baseVel-1)*nav_a));
                    }

                    fclose(fp);                   

                    for(i=1; i<param->nVelBins; i++){
                        param->velTable[i][1] = 0.0; 
                        param->velTable[i][2] = 0.0;
                        param->velTable[i][3] = 0.0;
                    }
/*
 *                  dump the dissipation problebility table
 * */
                    nav_a = 0.0; nav_ia = 0.0;
                    for(i=0; i<param->nDissBins; i++){
                        param->dissTable[i][1] /= param->cumulativeTime;
                        param->dissTable[i][2] /= param->cumulativeTime;
                        param->dissTable[i][3] /= param->cumulativeTime;
                        nav_a += param->dissTable[i][1];
                        nav_a += param->dissTable[i][2];
                        nav_ia += param->dissTable[i][3];
                    }
                    
                    snprintf(fileName, sizeof(fileName), "%s/unsigned-dissipation-%d.plt", DIR_PROPERTIES, home->cycle);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

					fprintf(fp, "variables = dissipation, width, p_unsigned_disspation_active, p_disspation_inactive\n"); 
		            fprintf(fp, "zone T = \"%d\", i = %d, j = 1, k = 1, SOLUTIONTIME = %d, F=POINT\n", 
                                home->cycle, param->nDissBins, home->cycle); 
                    fprintf(fp, "AUXDATA cumulativeTime = \"%f\"\n", param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_diss_active = \"%e\"\n", param->accumulation[4]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_diss_active = \"%e\"\n", param->accumulation[5]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA diss_active = \"%e\"\n", param->accumulation[6]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA diss_inactive = \"%e\"\n", param->accumulation[7]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_length_active = \"%e\"\n", param->aveLength[0]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_length_active = \"%e\"\n", param->aveLength[1]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA length_inactive = \"%e\"\n", param->aveLength[2]/param->cumulativeTime);
                    
                    i=0;
                    fprintf(fp, "%e %e %e %e\n", param->dissTable[i][0], param->dissMin,
                            (param->dissTable[i][1]+param->dissTable[i][2])/
                            (param->dissMin*nav_a), param->dissTable[i][3]/(param->dissMin*nav_a));

                    for(i=1; i<param->nDissBins; i++){
                        fprintf(fp, "%e %e %e %e\n", param->dissTable[i][0], param->dissMin*(param->baseDiss-1)*pow(param->baseDiss, ((real8)i)-1.0),
                                (param->dissTable[i][1]+param->dissTable[i][2])/
                                (param->dissMin*pow(param->baseDiss,(real8)i)*(param->baseDiss-1)*nav_a), 
                                param->dissTable[i][3]/(param->dissMin*pow(param->baseDiss,(real8)i)*(param->baseDiss-1)*nav_ia));
                    }

                    fclose(fp);                   

                    snprintf(fileName, sizeof(fileName), "%s/signed-dissipation-%d.plt", DIR_PROPERTIES, home->cycle);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

					fprintf(fp, "variables = dissipation, width, p_signed_disspation_active\n"); 
		            fprintf(fp, "zone T = \"%d\", i = %d, j = 1, k = 1, SOLUTIONTIME = %d, F=POINT\n", 
                                home->cycle, 2*param->nDissBins, home->cycle); 
                    fprintf(fp, "AUXDATA cumulativeTime = \"%f\"\n", param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_diss_active = \"%e\"\n", param->accumulation[4]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_diss_active = \"%e\"\n", param->accumulation[5]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA diss_active = \"%e\"\n", param->accumulation[6]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA diss_inactive = \"%e\"\n", param->accumulation[7]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA neg_length_active = \"%e\"\n", param->aveLength[0]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA pos_length_active = \"%e\"\n", param->aveLength[1]/param->cumulativeTime);
                    fprintf(fp, "AUXDATA length_inactive = \"%e\"\n", param->aveLength[2]/param->cumulativeTime);
                    
                    for(i=param->nDissBins-1; i>0; i--){
                        fprintf(fp, "%e %e %e\n", -param->dissTable[i][0], param->dissMin*(param->baseDiss-1)*pow(param->baseDiss, ((real8)i)-1.0),
                                param->dissTable[i][1]/
                                (param->dissMin*pow(param->baseDiss,(real8)i)*(param->baseDiss-1)*nav_a)); 
                    }
                    i=0;
                    fprintf(fp, "%e %e %e\n", -param->dissTable[i][0], param->dissMin, 
                            param->dissTable[i][1]/
                            (param->dissMin*nav_a));
                    fprintf(fp, "%e %e %e\n", param->dissTable[i][0], param->dissMin, param->dissTable[i][2]/
                            (param->dissMin*nav_a));

                    for(i=1; i<param->nDissBins; i++){
                        fprintf(fp, "%e %e %e\n", param->dissTable[i][0], param->dissMin*(param->baseDiss-1)*pow(param->baseDiss, ((real8)i)-1.0),
                                param->dissTable[i][2]/
                                (param->dissMin*pow(param->baseDiss,(real8)i)*(param->baseDiss-1)*nav_a));
                    }

                    fclose(fp);                   
                    for(i=1; i<param->nDissBins; i++){
                        param->dissTable[i][1] = 0.0;
                        param->dissTable[i][2] = 0.0;
                        param->dissTable[i][3] = 0.0;
                    }
#endif
                    snprintf(fileName, sizeof(fileName),
                             "%s/stress_Total_strain.plt", DIR_PROPERTIES);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

                    sign  = 1;
                    switch (param->loadType){
                        case 1:
                            sign = Sign(param->eRate);
                            break;
                        case 6:
                            sign = Sign(param->eRate);
                            break;
                        case 7:
                            sign = Sign(param->strRate);
                            break;
                        default:
                            break;
                    }

					if(home->cycle == 1){
					    fprintf(fp, "variables = step, "); 

                        fprintf(fp, "strain, stress, plastic_strain, ");
                        if(param->polycrystal) fprintf(fp, "twinning_deformation, ");

#ifndef VELOCITY_STATISTICS
                        if(param->polycrystal) fprintf(fp, "disDen, disDen_TB\n");
                        else fprintf(fp, "disDen\n"); 
#else
                        fprintf(fp, "disDen, effPlaRate, neg_v_squa_active, pos_v_squa_active, v_squa_active," 
                                    "v_squa_inactive, neg_tauv_active, pos_tauv_active," 
                                    " tauv_active, tauv_inactive, cumulativeTime"
                                    ", ave_len_neg_active, ave_len_pos_active, ave_len_inactive\n");
#endif
					}
        
                    fprintf(fp,"%d ", home->cycle);
                    sigijk = param->appliedStress[0]*al*bl        +
                             param->appliedStress[1]*am*bm        +
                             param->appliedStress[2]*an*bn        +
                             param->appliedStress[3]*(am*bn+an*bm)+
                             param->appliedStress[4]*(an*bl+al*bn)+ 
                             param->appliedStress[5]*(al*bm+bl*am);
                    
                    pstnijk =  param->totpStn[0]*al*bl        +
                               param->totpStn[1]*am*bm        +
                               param->totpStn[2]*an*bn        +
                               param->totpStn[3]*(am*bn+an*bm)+
                               param->totpStn[4]*(an*bl+al*bn)+
                               param->totpStn[5]*(al*bm+bl*am);
                    
                    ptstnijk = param->totTwinpStn[0]*al*bl        +
                               param->totTwinpStn[1]*am*bm        +
                               param->totTwinpStn[2]*an*bn        +
                               param->totTwinpStn[3]*(am*bn+an*bm)+
                               param->totTwinpStn[4]*(an*bl+al*bn)+
                               param->totTwinpStn[5]*(al*bm+bl*am);
                    if (param->loadType==4){
                        fprintf(fp, "%e %e %e ",param->netCyclicStrain*1.0E2, sigijk/1.0E6, pstnijk*1.0E2);
                        if(param->polycrystal) fprintf(fp, "%e ", ptstnijk*1.0E2);
                        if(fabs(param->stainCycle) >= param->maxCycle)param->finish = 1;
                        
                    }else {
                        fprintf(fp, "%e %e %e ",sign*param->eRate*1.0E2*param->timeNow, sign*sigijk/1.0E6, sign*pstnijk*1.0E2);
                        if(param->polycrystal) fprintf(fp, "%e ", sign*ptstnijk*1.0E2);

                        if(fabs(param->eRate*param->timeNow) > param->maxStn){
                            if(home->myDomain==0)printf("Finish: max strain: %e\n", param->maxStn);
                            param->finish = 1;
                        }
                        if(fabs(pstnijk) > param->maxPst){
                            if(home->myDomain==0)printf("Finish: max strain: %e\n", param->maxPst);
                            param->finish = 1;
                        }
                    }


#ifndef VELOCITY_STATISTICS
                    if(param->polycrystal) fprintf(fp, "%e  %e\n",param->disloDensity, param->disloDensity_TB);
                    else fprintf(fp, "%e\n",param->disloDensity); 
#else
                    fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", param->disloDensity, effPlaRate, param->accumulation[0]/param->cumulativeTime, 
                                                       param->accumulation[1]/param->cumulativeTime, 
                                                       param->accumulation[2]/param->cumulativeTime, 
                                                       param->accumulation[3]/param->cumulativeTime, 
                                                       param->accumulation[4]/param->cumulativeTime, 
                                                       param->accumulation[5]/param->cumulativeTime, 
                                                       param->accumulation[6]/param->cumulativeTime, 
                                                       param->accumulation[7]/param->cumulativeTime, 
                                                       param->cumulativeTime,
                                                       param->aveLength[0]/param->cumulativeTime,
                                                       param->aveLength[1]/param->cumulativeTime,
                                                       param->aveLength[2]/param->cumulativeTime);
#endif
                    fclose(fp);

                    snprintf(fileName, sizeof(fileName), "Finish");
                    if((fp = fopen(fileName, "r")) != (FILE *)NULL){
                        fclose(fp);
                        param->finish = 1;
                    }
#if 0
                    switch (param->loadType) {
                    case 1|7:

					if(firstSS && home->cycle < 2 *  param->savepropfreq){
					   fprintf(fp, "variables = strain, s1, s2, s3, s4, s5, s6, pstrain, disDen\n"); 
					   firstSS = 0;
					}

        		    fprintf(fp, "%e %e %e %e %e %e %e %e %e\n", param->eRate*param->timeNow*TC*1.0e2,
        		            param->appliedStress[0],
        		            param->appliedStress[1], param->appliedStress[2], param->appliedStress[3],
        		            param->appliedStress[4], param->appliedStress[5], param->totpStn[3]*TC*1.0e2, param->disloDensity);
					if(pstnijk*TC*1.0e2 > 0.1){
					    static int first = 1;
					    if(first){
						first = 0;
					        printf("0.1%% pstn: %e %e\n",sigijk*TC/1.0e6,param->disloDensity);
					    }
					}
					if(param->eRate*param->timeNow*TC > param->maxStn)param->finish = 1;
			        if(pstnijk*TC > param->maxPst)param->finish = 1;
            
                        break;
                    case 4:
                        fprintf(fp, "%e %e %e %d\n",
                                param->netCyclicStrain, sigijk,
                                param->timeNow,param->numLoadCycle);
                        break;
                    }
#endif

                }
            }        
#ifdef VELOCITY_STATISTICS
            param->accumulation[0] = 0.0; 
            param->accumulation[1] = 0.0; 
            param->accumulation[2] = 0.0; 
            param->accumulation[3] = 0.0; 
            param->accumulation[4] = 0.0; 
            param->accumulation[5] = 0.0; 
            param->accumulation[6] = 0.0; 
            param->accumulation[7] = 0.0; 
            param->aveLength[0] = 0.0;
            param->aveLength[1] = 0.0;
            param->aveLength[2] = 0.0;
            param->cumulativeTime = 0.0;
#endif
            break;
               
        case ALL_EPS:
            if (home->myDomain == 0) {
                snprintf(fileName, sizeof(fileName), "%s/alleps",
                         DIR_PROPERTIES);
                if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                    Fatal("can not open file %s", fileName);
                }
                fprintf(fp,"%d %e %e %e %e %e %e %e %e\n",
                        home->cycle, param->timeNow,
                        param->totpStn[0],param->totpStn[1],param->totpStn[2],
                        param->totpStn[3],param->totpStn[4],param->totpStn[5],
                        param->disloDensity );
                fclose(fp);
            }
            break;
        
        case EPSDOT:
            if (home->myDomain == 0) {
                if (param->realdt > 1.e-20){
                    dpstnijk =  param->delpStrain[0]*al*bl     +
                                param->delpStrain[1]*am*bm     +
                                param->delpStrain[2]*an*bn     +
                                param->delpStrain[3]*(am*bn +bm*an)+
                                param->delpStrain[4]*(an*bl +bn*al)+
                                param->delpStrain[5]*(al*bm +bl*am);
        
                    snprintf(fileName, sizeof(fileName), "%s/epsdot",
                             DIR_PROPERTIES);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

                    fprintf(fp, "%e %e\n", param->timeNow,
                            fabs(dpstnijk/param->realdt));
                    fclose(fp);
                }
            }
            break;
        
        case DENSITY_DELTA:
/*
 *          Currently only supported for BCC mobilities
 */
            if (param->materialType == MAT_TYPE_BCC) {

/*
 *              We track 2 values (density gain and density loss) for
 *              each of 7 burgers vectors
 */
                numBurgVectors = 7;
                numItems = numBurgVectors * 2;
#ifdef PARALLEL

                for (i = 0; i < numItems; i++) {
                    totDensityChange[i] = 0.0;
                }

                MPI_Reduce(param->densityChange, totDensityChange, numItems,
                           MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
                for (i = 0; i < numItems; i++) {
                    totDensityChange[i] = param->densityChange[i];
                }
#endif
                if (home->myDomain == 0) {
                    real8 totGain=0.0, totLoss=0.0;

                    for (i = 0; i < numBurgVectors; i++) {
                        totGain += totDensityChange[i];
                        totLoss += totDensityChange[numBurgVectors+i];
                    }

                    snprintf(fileName, sizeof(fileName), "%s/density_delta",
                             DIR_PROPERTIES);
                    if((fp = fopen(fileName, "a")) == (FILE *)NULL){
                        Fatal("can not open file %s", fileName);
                    }

                    fprintf(fp, "%e %e %e ", param->eRate*param->timeNow,
                            totGain, totLoss);

                    for (i = 0; i < numItems; i++) {
                        fprintf(fp, "%e ", totDensityChange[i]);
                    }
                    fprintf(fp, "\n");

                    fclose(fp);
                }
/*
 *              Zero out the accumulated values since we've written the
 *              data to disk.
 */
                for (i = 0; i < numItems; i++) {
                    param->densityChange[i] = 0.0;
                }
            }
            break;

        default:
            Fatal("WriteProp: input parameter (property=%d) invalid", property);
            break;

        }  /* switch() */

        return;
}
