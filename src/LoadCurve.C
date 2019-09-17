/**************************************************************************
 *
 *      Author:  Moono Rhee
 *      Function: LoadCurve
 *
 *      Description: This subroutine defines the type of load curves.
 *                   Works only with the conventional x-y-z (global)
 *                   coordinate system.  If loading axis rotated, the
 *                   loading axis can be rotated, or One can rotate
 *                   the initial dislocation configuration to a
 *                   "laboratory" coordinate system.
 *
 *                   Types of load curves:
 *                      0  Creep
 *                      1  Constant strain test
 *                      2  Displacement-controlled
 *                      3  Junction unzipping jump test
 *                      4  Total strain controlled cyclic load
 *                      5  Plastic strain controlled cyclic load
 *                      6  Load-time curve 
 *
 *      Last Modified:  01/03/2001 - original version
 *                      03/13/2003 - M. Rhee Removed anisotropic elastic
 *                                   constants.  Modified to include
 *                                   isotropic Hooke's law for arbitray
 *                                   loading. 
 *                      11/11/2003 - MasatoH Implementation of loading axis
 *                                   rotation due to accumuration of
 *                                   material spin.  Instead of crystal
 *                                   system, lab frame is rotated in opposite
 *                                   way
 *                      06/23/2004 - M.Rhee Added strain decomposition and
 *                                   density flux decompostion.  Modified
 *                                   message passing calls for all decomposed
 *                                   strain/density info
 *                      07/12/2004 - Masato Strain contolled cyclic load
 *                                   is implemented.   
 *
 ***************************************************************************/
#include "Home.h"
#include "Util.h"
#include <stdio.h>
#include <math.h>
#include <WriteProp.h>


/*
 *      Ss(): Sine for small angle i.e. Ss~x
 *      Cs(): Cosine for small angle i.e. Cs~1-x^2/2
 */
#define Ss(a) ((a))
#define Cs(a) (1.0 - (0.5 * (a)*(a)))


/*
 *      Function:     SpinMatrix
 *      Description:  Small rotation matrix for accumulaed
 *                    rotations around axis 1, 2 and 3.
 *
 *                    Cs() = Cosine for small angle
 *                    Ss() = Sin for small angle
 */
static void SpinMatrix(real8 p1, real8 p2, real8 p3, real8 Rspin[3][3])
{   
        Rspin[0][0] =  Cs(p3)*Cs(p2);
        Rspin[1][1] =  Cs(p3)*Cs(p1) + Ss(p3)*Ss(p2)*Ss(p1);
        Rspin[2][2] =  Cs(p2)*Cs(p1);
        Rspin[0][1] = -Ss(p3)*Cs(p1) + Cs(p3)*Ss(p1)*Ss(p2);
        Rspin[1][2] = -Cs(p3)*Ss(p1) + Ss(p3)*Ss(p2)*Cs(p1);
        Rspin[2][0] = -Ss(p2);
        Rspin[0][2] =  Ss(p3)*Ss(p1) + Cs(p3)*Cs(p1)*Ss(p2);
        Rspin[1][0] =  Ss(p3)*Cs(p2);
        Rspin[2][1] =  Cs(p2)*Ss(p1);
        
        return;    
}


void LoadCurve(Home_t *home, real8 deltaStress[3][3])
{
        int     i, j, k, loadtype, indxerate;
        int     numLoadCycle, numLoadCycle2;       
        real8   youngs, erate, dtt;
        real8   shr, dstn[6];
        real8   modulus, dpl_stn, dStress, amag, al, am, an;
        real8   sigijk, stn_cut;
        real8   phi1, phi2, phi3;
        real8   Rspin[3][3];
        real8   tempedot[3], temppassedot[3];
        real8   pstnijk, eAmp, timeNow, cTimeOld;
        real8   dCyclicStrain;
        real8   totCyclicStrain, netCyclicStrain;
		real8   dpstn[6], dsig[6], destn[6];
        real8   tmpdSig[6];
		real8   (*C)[6], (*Sij)[6];
        real8   effPlaRate, stepStr;
        Param_t *param;
        int     maxDirID = -1;
        real8   maxDirVal = -1E50;
        char name[50];

        real8   displacement, vec[3], locMin, gloMin, py;
        Node_t  *node;
        int     locVal, gloVal;
        
        TimerStart(home, LOADCURVE);
        
        param     = home->param;
        loadtype  = param->loadType;
        shr       = param->shearModulus;
        youngs    = 2.0 * shr * (1.0+param->pois);
        erate     = param->eRate;
        dtt       = param->realdt;
        indxerate = param->indxErate;

        sigijk          = 0;
        totCyclicStrain = 0.0;
		
		for (i = 0; i < 6; i++) {
            dsig[i] = 0.0;
            destn[i] = 0.0;
            tmpdSig[i] = 0.0;
            dstn[i] = 0.0;
        }
        
/*
 *      Update young's modulus in Param_t so the correct value is written
 *      into the control file.  Value is written to control file for
 *      informational purposes only.
 */
        param->YoungsModulus = youngs;

/*
 *      for cyclic load
 */
        eAmp            = param->eAmp;
        timeNow         = param->timeNow;
        cTimeOld        = param->cTimeOld;
        numLoadCycle    = param->numLoadCycle;
        netCyclicStrain = param->netCyclicStrain;
        dCyclicStrain   = param->dCyclicStrain;
        
        deltaStress[0][0] = param->appliedStress[0];
        deltaStress[1][1] = param->appliedStress[1];
        deltaStress[2][2] = param->appliedStress[2];
        deltaStress[1][2] = param->appliedStress[3];
        deltaStress[2][0] = param->appliedStress[4];
        deltaStress[0][1] = param->appliedStress[5];
        
#if debug
        printf("cosinesmall in LC %e\n",Cs(0.001));
        SpinMatrix(0.1,0.1,0.1,Rspin);
        
        gnu_plot(home);
        
        if (home->cycle == 2) {
            Fatal("Stopping at DeltaPlasticStrain to debug");
        }
#endif
        
/*
 *      If we're including osmotic forces on dislocation segments, we
 *      need to use the delta plastic strain to adjust the vacancy
 *      concentration.  Basically, sum the diagonal of the delta plastic
 *      strain and add to the vacancy concentration.
 */
        if (param->vacancyConcEquilibrium > 0.0) {
            param->vacancyConc += (param->delpStrain[0] +
                                   param->delpStrain[1] +
                                   param->delpStrain[2]);
        }

/*
 *      Initialize the fluxtot values for later.  Unfortunately, different
 *      arrays and array sizes depending on material type.  Sigh.
 */
        switch (param->materialType) {
            case MAT_TYPE_BCC:
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 7; j++) {
                        param->fluxtot[i][j] = param->dfluxtot[i][j];
                    }
                }
                break;
            case MAT_TYPE_FCC:
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 7; j++) {
                        param->FCC_fluxtot[i][j] = param->FCC_dfluxtot[i][j];
                    }
                }
                break;
        }
        
        
        for (i = 0; i < 6; i++) {
            param->totpStn[i] += param->delpStrain[i];                       
            param->totTwinpStn[i] += param->delpTwinStrain[i];                       
            param->totpSpn[i] += param->delpSpin[i];                         
            param->totedgepStrain[i] += param->dedgepStrain[i];               
            param->totscrewpStrain[i] += param->dscrewpStrain[i];             
        }                                                                    
        
/*
 *      Part for Loading Axis Rotation due to Small Deformation Spin.
 *
 *      Some changes in rotation angles due to the deformation spins
 *      around x, y, and z axis
 */
        phi1 = - param->delpSpin[3];   
        phi2 =   param->delpSpin[4];   
        phi3 = - param->delpSpin[5];   
        
/*
 *      Matrix for (combined) rotation around x,y,and z in the sequence.
 *      This sequential rotation is correct only for small changes 
 *      in the angles since real material rotation occurs simultaneously.
 *      For counter-rotation, sign of phi is flipped.
 */
        SpinMatrix(phi1, phi2, phi3, Rspin);  
        
/*
 *      Compute nodal velocity : bug is fixed. Vector is address
 */
        tempedot[0] = param->edotdir[0];
        tempedot[1] = param->edotdir[1];
        tempedot[2] = param->edotdir[2];

        temppassedot[0] = 0.0;
        temppassedot[1] = 0.0;
        temppassedot[2] = 0.0;
        
        Matrix33Vector3Multiply(Rspin, tempedot, temppassedot);
        
        param->edotdir[0] = temppassedot[0];
        param->edotdir[1] = temppassedot[1];
        param->edotdir[2] = temppassedot[2];
        
#if 0
        printf("phi = %e %e %e \n", phi1, phi2, phi3);
        printf("Rspin = %e %e %e \n", Rspin[0][0], Rspin[0][1], Rspin[0][2]);
        printf("Rspin = %e %e %e \n", Rspin[1][0], Rspin[1][1], Rspin[1][2]);
        printf("Rspin = %e %e %e \n", Rspin[2][0], Rspin[2][1], Rspin[2][2]);
        printf("newl = %e %e %e \n", param->edotdir[0],
               param->edotdir[1],param->edotdir[2]);
        
        printf("param->timeNow = %e \n", param->timeNow);
        printf("youngs=%f\n", youngs);
        printf("erate=%f\n", erate);
        printf("  dtt=%e\n", dtt);
        printf("indxerate=%d\n", indxerate);
#endif
        
/*
 *      Arbitrary loading direction but keep in lab frame
 */
        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];

        amag = sqrt(al*al + am*am + an*an);

        al /= amag;
        am /= amag;
        an /= amag;
        
        switch(loadtype) {
/*
 *          creep - what we have been using so this should be default
 */
            case 0:
				
				Matrix33Init(deltaStress);
                break;
/*
 *          constant strain rate
 */
            case 1:
#if	0
/*
 *              Cover for specific loading direction also
 */
                dpl_stn=  param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs;
                } else {
                   modulus = 2.0 * shr;
                }

/*
 *              local in the [l m n] frame
 */
                dStress= modulus * (erate*dtt - dpl_stn);
        
/*#if 0
                printf("al=%e,am=%e,an=%e\n", al, am, an);
                printf("erate=%e,dtt=%e,dpl_stn=%e\n", erate, dtt, dpl_stn);
                printf("modulus=%e,dStress=%e\n", modulus, dStress);
*/
/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;
        
                param->totstraintensor[0] = erate * param->timeNow *  al*al;
                param->totstraintensor[1] = erate * param->timeNow *  am*am;
                param->totstraintensor[2] = erate * param->timeNow *  an*an;
                param->totstraintensor[3] = erate * param->timeNow *  an*am;
                param->totstraintensor[4] = erate * param->timeNow *  an*al;
                param->totstraintensor[5] = erate * param->timeNow *  al*am;
#endif
				
				for (i = 0; i < 6; i++) {
                    dpstn[i] = param->delpStrain[i];
                }
/*
 *              Set pointers to the precomputed 6x6 elastic constant
 *              matrix (C) and its inverse, the compliance matrix (Sij)
 */
                C   = param->elasticConstantMatrix;
                Sij = param->elasticConstantMatrixInv;
				
                effPlaRate = sqrt((dpstn[0]*dpstn[0]+ 
                                   dpstn[1]*dpstn[1]+ 
                                   dpstn[2]*dpstn[2]+ 
                               2.0*dpstn[3]*dpstn[3]+ 
                               2.0*dpstn[4]*dpstn[4]+ 
                               2.0*dpstn[5]*dpstn[5])*2.0/3.0)/param->lastDT;

                if(param->plaRateCtr <= 0.0 || 
                  (param->plaRateCtr > 0 && effPlaRate <= param->plaRateCtr)){
                    dstn[0] = param->loadDir[0]*dtt*erate;
                    dstn[1] = param->loadDir[1]*dtt*erate;
                    dstn[2] = param->loadDir[2]*dtt*erate;
                    dstn[3] = param->loadDir[3]*dtt*erate;
                    dstn[4] = param->loadDir[4]*dtt*erate;
                    dstn[5] = param->loadDir[5]*dtt*erate;
                }

                dsig[0] = (dstn[0] - dpstn[0]) /Sij[0][0];
                dsig[1] = (dstn[1] - dpstn[1]) /Sij[1][1];
                dsig[2] = (dstn[2] - dpstn[2]) /Sij[2][2];
                dsig[3] = 2.0 * (dstn[3] - dpstn[3]) /Sij[3][3];
                dsig[4] = 2.0 * (dstn[4] - dpstn[4]) /Sij[4][4];
                dsig[5] = 2.0 * (dstn[5] - dpstn[5]) /Sij[5][5];

				
				for(i=0; i<6; i++){
					for(j=0; j<6; j++){
						destn[i] += Sij[i][j] * dsig[j];
					}
					param->elasticStrainTensor[i] += destn[i];
				}
				
				for(i=0; i<6; i++){
					param->appliedStress[i] += dsig[i];
					
				}
				
                deltaStress[0][0] = dsig[0];
                deltaStress[1][1] = dsig[1];
                deltaStress[2][2] = dsig[2];

                deltaStress[1][2] = dsig[3];
                deltaStress[2][0] = dsig[4];
                deltaStress[0][1] = dsig[5];

                deltaStress[2][1] = deltaStress[1][2];
                deltaStress[0][2] = deltaStress[2][0];
                deltaStress[1][0] = deltaStress[0][1];

//                if(home->myDomain == 0)printf("Loading %e %e %e %e\n", erate*dtt, dpstn[3], deltaStress[2][1]/dtt, param->appliedStress[3]);
                break;
        
/*
 *          jump test
 */
            case 2:
                stn_cut = 5.e-12; 
                dStress = 1.e5;
                dpl_stn = param->delpStrain[indxerate-1];
        
                if ((dpl_stn > 0.0) && (dpl_stn < stn_cut)) {
                    param->appliedStress[indxerate-1] += dStress;
                    home->cycle = 1;
                } else if (dpl_stn < 0){
                } else {
                    if (home->cycle > 200000){
                         printf("Must be ok now printing the critical stress \n");
                         printf("Critical Stress =%e \n", sigijk);
                         printf("  minSeg   =%e \n", param->minSeg);
                         printf("  maxSeg   =%e \n", param->maxSeg);
                         Fatal("Doing clean terminate in LoadCurve");
                    }
               }
        
               home->cycle++;
        
               sigijk  =  param->appliedStress[indxerate-1];
               dpl_stn =  param->delpStrain[indxerate-1];
        
               if (home->cycle % 50 == 0) {
                   printf("sig=%e stn=%e dt=%e cyc=%d\n", sigijk,
                          dpl_stn, param->realdt, home->cycle);
               }
        
               break;
         
/*
 *          Junction unzipping jump test; not for general case yet */
            case 3:
#if 1
                stn_cut = 1.0e-10; /* lowbound strain cut for numerical */
                                   /* stability                         */
                dStress = 1.e5;
                dpl_stn = param->delpStrain[indxerate-1];
        
                if ((dpl_stn > 0.0) && (dpl_stn < stn_cut)) {
                    param->appliedStress[indxerate-1] += dStress;
                }

                sigijk =  param->appliedStress[indxerate-1];
                dpl_stn=  param->delpStrain[indxerate-1];
        
                if (home->cycle % 10 == 0) {
                    printf("sig=%e stn=%e dt=%e cyc=%d\n", sigijk,
                           dpl_stn, param->realdt, home->cycle);
                }
#endif
                break;
        
/*
 *          strain control cyclic load 
 *
 *              stainCycle    = current loading cycle
 *              eAmp          = strain amplitude for each side
 *              dCyclicStrain = change in the strain for each side
 *              acumStrain    = accumulated strain
 *              sgnLoad       = sign of load
 */
            case 4:
        
/*
 *              Cover for specific loading direction also
 */
                dpl_stn =  param->delpStrain[0]*al*al +
                           param->delpStrain[1]*am*am +
                           param->delpStrain[2]*an*an +
                           2.0*param->delpStrain[3]*am*an +
                           2.0*param->delpStrain[4]*an*al +
                           2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs; 
                } else {
                   modulus = 2.0 * shr; 
                }
        
                dCyclicStrain = erate*dtt;
                param->dCyclicStrain = dCyclicStrain;
        
/*
 *              local in the [l m n] frame
 */ 
                dStress= modulus * (dCyclicStrain - dpl_stn);
        
                totCyclicStrain = fabs(erate*timeNow); 
                numLoadCycle    = rint(0.5*totCyclicStrain/eAmp);
                numLoadCycle2   = rint(0.5*totCyclicStrain/eAmp-0.5);
        
                netCyclicStrain = fmod(totCyclicStrain, 2*eAmp); 
        
                if (fabs(netCyclicStrain) > eAmp) {
                    netCyclicStrain = 2*eAmp - netCyclicStrain;
                }
        
                netCyclicStrain = pow(-1, numLoadCycle2) *
                                  fabs(netCyclicStrain); 
        
                param->netCyclicStrain = netCyclicStrain;
        
#if 0
                printf("loading cycle 1&2 %d %d\n", numLoadCycle,
                       numLoadCycle2);
                printf("net strain %e\n", netCyclicStrain);
                printf("loading cycle %d\n", numLoadCycle);
                printf("Load Curve: dtt,totSt %e %e \n", dtt, totCyclicStrain);
                printf("Load Curve: dtt,dttOld dCyclic %e %e %e %e\n",
                       timeNow, cTimeOld, dCyclicStrain, eAmp);
#endif
        
                cTimeOld = timeNow;
                erate = fabs(erate)*pow(-1,numLoadCycle);
                dCyclicStrain = 0;
                param->cTimeOld = cTimeOld;
                param->numLoadCycle = numLoadCycle;
                param->eRate = erate;
        
/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am; 
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;
        
/*
                param->totstraintensor[0] = erate * param->timeNow *  al*al;
                param->totstraintensor[1] = erate * param->timeNow *  am*am;
                param->totstraintensor[2] = erate * param->timeNow *  an*an;
                param->totstraintensor[3] = erate * param->timeNow *  an*am;
                param->totstraintensor[4] = erate * param->timeNow *  an*al;
                param->totstraintensor[5] = erate * param->timeNow *  al*am; 
*/
                param->totstraintensor[0] = netCyclicStrain *  al*al;
        
                break;
        
        
/*
 *             Plastic strain control cyclic load 
 *                 stainCycle    = current loading cycle
 *                 eAmp          = strain amplitude for each side
 *                 dCyclicStrain = change in the strain for each side
 *                 acumStrain    = accumulated strain
 *                 sgnLoad       = sign of load
 */
            case 5:
              
/*
 *              Cover for specific loading direction als0
 */
                pstnijk = param->totpStn[0]*al*al +
                          param->totpStn[1]*am*am +
                          param->totpStn[2]*an*an +
                          2.0*param->totpStn[3]*am*an +
                          2.0*param->totpStn[4]*an*al +
                          2.0*param->totpStn[5]*al*am;
        
                dpl_stn = param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;
        
                if (indxerate <= 3) {
                   modulus = youngs; 
                } else  {
                   modulus = 2.0 * shr; 
                }
        
                dCyclicStrain = erate*dtt;
                param->dCyclicStrain = dCyclicStrain;
        
/*
 *              local in the [l m n] frame
 */ 
                dStress= modulus * (dCyclicStrain - dpl_stn);
        
                totCyclicStrain += fabs(dpl_stn); 
                numLoadCycle    = rint(0.5*pstnijk/eAmp);
                numLoadCycle2   = rint(0.5*pstnijk/eAmp-0.5);
        
                netCyclicStrain = fmod(pstnijk, 2*eAmp); 
        
                if (fabs(netCyclicStrain) > eAmp ) {
                    netCyclicStrain = 2*eAmp - netCyclicStrain;
                }
        
                netCyclicStrain = pow(-1, numLoadCycle2) *
                                  fabs(netCyclicStrain); 

                param->netCyclicStrain = netCyclicStrain;
        
#if 0
                printf("loading cycle 1&2 %d %d\n", numLoadCycle,
                       numLoadCycle2);
                printf("net strain %e\n", netCyclicStrain);
                printf("loading cycle %d\n", numLoadCycle);
                printf("Load Curve: dtt,totSt %e %e \n", dtt, totCyclicStrain);
                printf("Load Curve: dtt,dttOld dCyclic %e %e %e %e\n",
                       timeNow, cTimeOld, dCyclicStrain, eAmp);
#endif
                cTimeOld = timeNow;
                erate = fabs(erate) * pow(-1,numLoadCycle);
                dCyclicStrain = 0;
                param->cTimeOld = cTimeOld;
                param->numLoadCycle = numLoadCycle;
                param->eRate = erate;
        
/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;  
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;
        
/*   
 *               param->totstraintensor[0] = erate * param->timeNow *  al*al;
 *               param->totstraintensor[1] = erate * param->timeNow *  am*am;
 *               param->totstraintensor[2] = erate * param->timeNow *  an*an;
 *               param->totstraintensor[3] = erate * param->timeNow *  an*am;
 *               param->totstraintensor[4] = erate * param->timeNow *  an*al;
 *               param->totstraintensor[5] = erate * param->timeNow *  al*am; 
 */
                param->totstraintensor[0] = netCyclicStrain *  al*al;
                param->totstraintensor[1] = netCyclicStrain *  am*am;
                param->totstraintensor[2] = netCyclicStrain *  an*an;
                param->totstraintensor[3] = netCyclicStrain *  an*am;
                param->totstraintensor[4] = netCyclicStrain *  an*al;
                param->totstraintensor[5] = netCyclicStrain *  al*am; 
        
                break;

/*
 *          User defined load-time curve
 */
            case 6:
                break;
/*
 *          Stress loading with plastic strain control
 */
            case 7:
                C   = param->elasticConstantMatrix;
                Sij = param->elasticConstantMatrixInv;
    
                effPlaRate = sqrt((dpstn[0]*dpstn[0]+ 
                                   dpstn[1]*dpstn[1]+ 
                                   dpstn[2]*dpstn[2]+ 
                               2.0*dpstn[3]*dpstn[3]+ 
                               2.0*dpstn[4]*dpstn[4]+ 
                               2.0*dpstn[5]*dpstn[5])*2.0/3.0)/param->lastDT;

                if(param->plaRateCtr <= 0.0 || 
                  (param->plaRateCtr > 0 && effPlaRate <= param->plaRateCtr)){
                    for(i=0; i<6; i++){
                        dsig[i] = param->strRate*dtt*param->loadDir[i];

                        if(fabs(param->strRate*param->loadDir[i]) > maxDirVal){
                            maxDirID = i;
                            maxDirVal = fabs(param->strRate*param->loadDir[i]);
                        }
                    }
                
                }

                for(i=0; i<6; i++){
                    for(j=0; j<6; j++){
                        destn[i] += Sij[i][j] * dsig[j];
                    }
                    param->elasticStrainTensor[i] += destn[i];
                }
                
                for(i=0; i<6; i++){
                    param->appliedStress[i] += dsig[i];
                }
				
                if(param->maxStr > 0){
                    if(param->maxStr < fabs(param->appliedStress[maxDirID]) &&
                       (Sign(param->strRate*param->loadDir[maxDirID]) * 
                       Sign(param->appliedStress[maxDirID])) > 0){
                        param->strRate = -param->strRate;
#ifdef PARALLEL
                        MPI_Barrier(MPI_COMM_WORLD);
#endif
                        snprintf(name, 50, "-%05.0fMPa-%d", 
                                 param->appliedStress[maxDirID]/1.0E6, home->cycle);
                        GenerateOutput(home, STAGE_DEBUG, name);
                    }
                }

                deltaStress[0][0] = dsig[0];
                deltaStress[1][1] = dsig[1];
                deltaStress[2][2] = dsig[2];

                deltaStress[1][2] = dsig[3];
                deltaStress[2][0] = dsig[4];
                deltaStress[0][1] = dsig[5];

                deltaStress[2][1] = deltaStress[1][2];
                deltaStress[0][2] = deltaStress[2][0];
                deltaStress[1][0] = deltaStress[0][1];

                break;
        
/*
 *          Step Stress loading.
 */
            case 8:
				for (i = 0; i < 6; i++) {
                    dpstn[i] = param->delpStrain[i];
                }
                Sij = param->elasticConstantMatrixInv;
    
                effPlaRate = sqrt((dpstn[0]*dpstn[0]+ 
                                   dpstn[1]*dpstn[1]+ 
                                   dpstn[2]*dpstn[2]+ 
                               2.0*dpstn[3]*dpstn[3]+ 
                               2.0*dpstn[4]*dpstn[4]+ 
                               2.0*dpstn[5]*dpstn[5])*2.0/3.0)/param->lastDT;

                locVal = 1;
                gloVal = 0;

                locMin = param->position;
	            for(i=0;i<home->newNodeKeyPtr;i++){
                    if((node = home->nodeKeys[i]) == (Node_t *) NULL)continue;

                    VECTOR_ZERO(vec);
                    vec[1] = node->y - param->position;
                    ZImage(param, vec, vec+1, vec+2);
    
                    locVal &= int(fabs(vec[1]) >= 2000);

                    if(node->y > 0.25*param->Ly){
                        py = node->y - param->Ly;
                    }else{
                        py = node->y;
                    }
                    locMin = MIN(locMin, py);
                }

#ifdef PARALLEL
                MPI_Allreduce(&locVal, &gloVal, 1, MPI_INT,
                      MPI_MIN, MPI_COMM_WORLD);

                MPI_Allreduce(&locMin, &gloMin, 1, MPI_DOUBLE,
                      MPI_MIN, MPI_COMM_WORLD);
#else
                gloVal = locVal;
                gloMin = locMin;
#endif

                locVal = 0;
                if((fabs(effPlaRate) <= param->criticalEff && param->disRule == 0) 
                   || gloVal == 1){
                    locVal = 1;
                    
                    stepStr = param->targetStress/30.0;

                    dsig[0] = stepStr*param->loadDir[0];
                    dsig[1] = stepStr*param->loadDir[1];
                    dsig[2] = stepStr*param->loadDir[2];
                    dsig[3] = stepStr*param->loadDir[3];
                    dsig[4] = stepStr*param->loadDir[4];
                    dsig[5] = stepStr*param->loadDir[5];
                    
                    if(gloVal == 1){
                        vec[1] = param->position + 2000;
                        FoldBox(param, vec, vec+1, vec+2);
                        param->position = vec[1];
                        param->disRule = 1;
                    }else{
                        vec[1] = gloMin;
                        FoldBox(param, vec, vec+1, vec+2);
                        param->position = vec[1];
                    }

                }

#ifdef PARALLEL
                MPI_Allreduce(&locVal, &gloVal, 1, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);

#else
                gloVal = locVal;
#endif

                if(gloVal==1){
                    snprintf(name, 50, "-%07.2fMPa", param->appliedStress[3]/1.0E6);
                    GenerateOutput(home, STAGE_DEBUG, name);
                }

                if(param->appliedStress[3] > param->targetStress)param->finish = 1;

                for(i=0; i<6; i++){
                    for(j=0; j<6; j++){
                        destn[i] += Sij[i][j] * dsig[j];
                    }
                    param->elasticStrainTensor[i] += destn[i];
                }
                
                for(i=0; i<6; i++){
                    param->appliedStress[i] += dsig[i];
                }
				
                deltaStress[0][0] = dsig[0];
                deltaStress[1][1] = dsig[1];
                deltaStress[2][2] = dsig[2];

                deltaStress[1][2] = dsig[3];
                deltaStress[2][0] = dsig[4];
                deltaStress[0][1] = dsig[5];

                deltaStress[2][1] = deltaStress[1][2];
                deltaStress[0][2] = deltaStress[2][0];
                deltaStress[1][0] = deltaStress[0][1];

                break;
            default:
                Fatal("Load curves not defined. Stopping the program. \n");
            break;

        }  /* end: switch(loadtype) */
        
        TimerStop(home, LOADCURVE);
#if 0        
        deltaStress[0][0] = param->appliedStress[0] - deltaStress[0][0];
        deltaStress[1][1] = param->appliedStress[1] - deltaStress[1][1];
        deltaStress[2][2] = param->appliedStress[2] - deltaStress[2][2];
        deltaStress[1][2] = param->appliedStress[3] - deltaStress[1][2];
        deltaStress[2][0] = param->appliedStress[4] - deltaStress[2][0];
        deltaStress[0][1] = param->appliedStress[5] - deltaStress[0][1];
        deltaStress[2][1] = deltaStress[1][2];
        deltaStress[0][2] = deltaStress[2][0];
        deltaStress[1][0] = deltaStress[0][1];
#endif		
		
		
//		printf(" param->appliedStress[2] = %f\n",param->appliedStress[2]);
        return;
}
