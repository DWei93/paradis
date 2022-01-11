/***************************************************************************
 *   
 *      Module:     DeltaPlasticStrain.
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will invoke the version of DeltaPlasticStrain*() 
 *                    appropriate to the type of material being simulated
 *
 ***************************************************************************/
#include "Home.h"
#include "Mobility.h"

void DeltaPlasticStrain(Home_t *home)
{
    Param_t     *param;
    int         totEleNum, i, j;
    real8       *locVals, *gloVals;

    param = home->param;
    totEleNum = param->meshNumX*param->meshNumY*param->meshNumZ;

    switch(home->param->materialType) {

        case MAT_TYPE_BCC:
            DeltaPlasticStrain_BCC(home);
            break;

        case MAT_TYPE_FCC:
            DeltaPlasticStrain_FCC(home);
            break;
    }


#ifdef PARALLEL
    locVals = (real8 *)calloc(totEleNum*NUM_COMPONENT_OF_PLASTIC_STRAIN, sizeof(real8));
    gloVals = (real8 *)calloc(totEleNum*NUM_COMPONENT_OF_PLASTIC_STRAIN, sizeof(real8));

    for (i=0; i<totEleNum; i++){
        for(j = 0; j<NUM_COMPONENT_OF_PLASTIC_STRAIN; j++){
            locVals[i*NUM_COMPONENT_OF_PLASTIC_STRAIN+j] = home->pstnEle[i][j];
        }
    }

	MPI_Reduce(locVals, gloVals, NUM_COMPONENT_OF_PLASTIC_STRAIN*totEleNum, MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);

    for (i=0; i<totEleNum; i++){
        for(j = 0; j<NUM_COMPONENT_OF_PLASTIC_STRAIN; j++){
            if(home->myDomain == 0){
                home->pstnEle[i][j] = gloVals[i*NUM_COMPONENT_OF_PLASTIC_STRAIN+j];
            }else{
                home->pstnEle[i][j] = 0.0;
            }
        }
    }

    free(locVals);
    free(gloVals);

#endif
#ifdef VELOCITY_STATISTICS

    real8   coeff = param->realdt/1.0E-9;

    if(param->disLength[0] > 1.0E-10){
        param->accumulation[0] += (param->delAccumulation[0]*param->realdt/param->disLength[0]/1.0E-9); 
        param->accumulation[4] += (param->delAccumulation[4]*param->realdt/param->disLength[0]/1.0E-9); 
    }

    if(param->disLength[1] > 1.0E-10){
        param->accumulation[1] += (param->delAccumulation[1]*param->realdt/param->disLength[1]/1.0E-9); 
        param->accumulation[5] += (param->delAccumulation[5]*param->realdt/param->disLength[1]/1.0E-9); 
    }

    if(param->disLength[0] > 1.0E-10 || param->disLength[1] > 1.0E-10){
        param->accumulation[2] += (param->delAccumulation[2]*param->realdt/(param->disLength[0]+param->disLength[1])/1.0E-9); 
        param->accumulation[6] += (param->delAccumulation[6]*param->realdt/(param->disLength[0]+param->disLength[1])/1.0E-9); 
    }

    if(param->disLength[2] > 1.0E-10){
        param->accumulation[3] += (param->delAccumulation[3]*param->realdt/param->disLength[2]/1.0E-9); 
        param->accumulation[7] += (param->delAccumulation[7]*param->realdt/param->disLength[2]/1.0E-9); 
    }

    param->cumulativeTime += (param->realdt/1.0E-9);
    param->aveLength[0] += (param->disLength[0]*coeff);
    param->aveLength[1] += (param->disLength[1]*coeff);
    param->aveLength[2] += (param->disLength[2]*coeff);

    int iBin;
    for(iBin=0; iBin<param->nDissBins; iBin++){
        param->dissTable[iBin][1] *= (param->realdt/1.0E-9);
        param->dissTable[iBin][2] *= (param->realdt/1.0E-9);
        param->dissTable[iBin][3] *= (param->realdt/1.0E-9);
    } 
    for(iBin=0; iBin<param->nVelBins; iBin++){
        param->velTable[iBin][1] *= (param->realdt/1.0E-9);
        param->velTable[iBin][2] *= (param->realdt/1.0E-9);
        param->velTable[iBin][3] *= (param->realdt/1.0E-9);
    } 
#endif
    return;
}
