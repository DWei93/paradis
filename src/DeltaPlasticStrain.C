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
    locVals = (real8 *)calloc(totEleNum*6, sizeof(real8));
    gloVals = (real8 *)calloc(totEleNum*6, sizeof(real8));

    for (i=0; i<totEleNum; i++){
        for(j = 0; j<6; j++){
            locVals[i*6+j] = home->pstnEle[i][j];
        }
    }

	MPI_Reduce(locVals, gloVals, 6*totEleNum, MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);

    for (i=0; i<totEleNum; i++){
        for(j = 0; j<6; j++){
            if(home->myDomain == 0){
                home->pstnEle[i][j] = gloVals[i*6+j];
            }else{
                home->pstnEle[i][j] = 0.0;
            }
        }
    }

    free(locVals);
    free(gloVals);

#endif

    param->accumulation[0] += (param->delAccumulation[0]*param->realdt/(param->disLength[0]+param->disLength[1])); 
    param->accumulation[1] += (param->delAccumulation[1]*param->realdt/param->disLength[0]); 
    param->accumulation[2] += (param->delAccumulation[2]*param->realdt/param->disLength[1]); 
    param->cumulativeTime += param->realdt;

    return;
}
