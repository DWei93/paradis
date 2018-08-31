#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>

#include <mpi.h>

#include <stdlib.h>
#include <algorithm>
#include <numeric>
#define PI  3.1415926535897932385

using namespace std;

typedef struct {
    double  xBoundMin, xBoundMax;
    double  yBoundMin, yBoundMax;
    double  meshSizeX, meshSizeY;
    double  radius, burgMag;
    int     meshNumX, meshNumY;
} Param_t;


typedef struct {
   double   x,y;
} Point_t;


double randm(int *seed)
{  
        int ia, ib, i1, i2, i3;

        ia = *seed / 32768;
        ib = (*seed % 32768);
        i1 = ia * 30485;
        i2 = ib * 30485;
        i3 = ib * 48603;
        i1 = (i1 % 65536);
        i3 = (i3 % 65536);
        i1 = i1 + i3 + 13849 + i2 / 32768 + (ia % 2) * 32768;
        i2 = (i2 % 32768) + 12659;
        i1 = i1 + i2 / 32768;
        i2 = (i2 % 32768);
        i1 = (i1 % 65536);

        *seed = i1 * 32768 + i2;

        return(*seed * 4.65661287308E-10);
}


void Coor2EleNum(Param_t *param, double x, double y, int &nx, int &ny){

	double	eleLx,eleLy,eleLz;
    double   fnx, fny;

	eleLx = param->meshSizeX;
	eleLy = param->meshSizeY;

	fnx = floor((x - param->xBoundMin)/eleLx);
	fny = floor((y - param->xBoundMin)/eleLx);
    
    nx = (int)fnx;
    ny = (int)fny;

    return;
}
    
double Interpolation(double f00, double fx00, double fy00, double fxy00,
                     double f01, double fx01, double fy01, double fxy01,
                     double f10, double fx10, double fy10, double fxy10,
                     double f11, double fx11, double fy11, double fxy11,
                     double x, double y, double &bil)
{
    double a00, a01, a02, a03;
    double a10, a11, a12, a13;
    double a20, a21, a22, a23;
    double a30, a31, a32, a33;
    double xx[4], yy[4], temp[4];
    int     i;

    if(x<0.0 || x>1 || y<0.0 || y>1){
        printf("%f %f\n", x, y);
        exit (1);
    }

    xx[0] = 1.0;
    xx[1] = x;
    xx[2] = x*x;
    xx[3] = xx[2]*x;

    yy[0] = 1.0;
    yy[1] = y;
    yy[2] = y*y;
    yy[3] = yy[2]*y;

    a00 = f00;
    a01 = fy00;
    a02 = -3*f00 + 3*f01 - 2*fy00 - fy01;
    a03 = 2*f00 - 2*f01 + fy00 + fy01; 
    a10 = fx00; 
    a11 = fxy00;
    a12 = -3*fx00 + 3*fx01 - 2*fxy00 - fxy01;
    a13 = 2*fx00 - 2*fx01 + fxy00 + fxy01;
    a20 = -3*f00 + 3*f10 - 2*fx00 - fx10;
    a21 = -2*fxy00 - fxy10 - 3*fy00 + 3*fy10; 
    a22 = 9*f00 - 9*f01 - 9*f10 + 9*f11 + 6*fx00 - 6*fx01 + 3*fx10 - 
            3*fx11 + 4*fxy00 + 
            2*fxy01 + 2*fxy10 + fxy11 + 6*fy00 + 3*fy01 - 6*fy10 - 3*fy11;
    a23 = -6*f00 + 6*f01 + 6*f10 - 6*f11 - 4*fx00 + 4*fx01 - 2*fx10 + 
            2*fx11 - 2*fxy00 - 
            2*fxy01 - fxy10 - fxy11 - 3*fy00 - 3*fy01 + 3*fy10 + 3*fy11;
    a30 = 2*f00 - 2*f10 + fx00 + fx10; 
    a31 = fxy00 + fxy10 + 2*fy00 - 2*fy10; 
    a32 = -6*f00 + 6*f01 + 6*f10 - 6*f11 - 3*fx00 + 3*fx01 - 3*fx10 + 
            3*fx11 - 2*fxy00 - 
            fxy01 - 2*fxy10 - fxy11 - 4*fy00 - 2*fy01 + 4*fy10 + 2*fy11;
    a33 = 4*f00 - 4*f01 - 4*f10 + 4*f11 + 2*fx00 - 2*fx01 + 2*fx10 - 
            2*fx11 + fxy00 + 
            fxy01 + fxy10 + fxy11 + 2*fy00 + 2*fy01 - 2*fy10 - 2*fy11;

    temp[0] = (1.0-x)*f00 + x*f10;
    temp[1] = (1.0-x)*f01 + x*f11;

    bil = temp[0]*(1.0-y) + temp[1]*y;

    for(i=0;i<4;i++){
        temp[0] = xx[0]*a00 + xx[1]*a10 + xx[2]*a20 +xx[3]*a30;
        temp[1] = xx[0]*a01 + xx[1]*a11 + xx[2]*a21 +xx[3]*a31;
        temp[2] = xx[0]*a02 + xx[1]*a12 + xx[2]*a22 +xx[3]*a32;
        temp[3] = xx[0]*a03 + xx[1]*a13 + xx[2]*a23 +xx[3]*a33;
    }


    return(temp[0]*yy[0]+temp[1]*yy[1]+temp[2]*yy[2]+temp[3]*yy[3]);
}

int main(int argc,char** argv)
{
    double  xBoundMax, yBoundMax, xBoundMin, yBoundMin;    
    double  meshSizeX, meshSizeY, nXs, nYs;
    double  burgMag, (*potParams)[3];
    double  potExpection, potVariance, sigma;
    double  *localArray, *globalArray;
    double  radius, dx, dy, l1, l2;
    int     i, j, k, myDomain, nDomains;
    int     nPotParams, nMeshPoints;
    int     numX1, numY1;

    double  A, B, C, E, D, uni[2];
    double  potential, x0, y0, xj, yj, xi, yi;
    int     ii, jj, imin, imax, jmin, jmax;
    int     seed;
    double  e = 2.718281828459045;
    double  f, fx, fy, fxy;
    double  bli;
    
    string  str;
    char    name[20];

    Param_t param;
    
    seed = time(0);
    E = potExpection;
    D = potVariance;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myDomain);
    MPI_Comm_size(MPI_COMM_WORLD, &nDomains);


    param.xBoundMax = 10000.0;
    param.yBoundMax = 10000.0;
    xBoundMax = param.xBoundMax;
    yBoundMax = param.yBoundMax;

    param.xBoundMin = -10000.0;
    param.yBoundMin = -10000.0;
    xBoundMin = param.xBoundMin;
    yBoundMin = param.yBoundMin;

    param.meshSizeX = 5.0;
    param.meshSizeY = 5.0;
    param.radius = 100.0;
    meshSizeX = param.meshSizeX;
    meshSizeY = param.meshSizeY;
    radius = param.radius;

    nPotParams = 10000;
    potExpection = 0.04;
    potVariance = 0.0;
    sigma = 20.0;
    param.burgMag = 2.556e-10;
    burgMag = param.burgMag;

    nXs = floor((xBoundMax-xBoundMin)/meshSizeX);
    nYs = floor((yBoundMax-yBoundMin)/meshSizeY);

    param.meshNumX = (int)nXs;
    param.meshNumY = (int)nYs;

    nMeshPoints = ((int)nXs+1)*((int)nYs+1);
    numX1 = (int)nXs + 1;
    numY1 = (int)nYs + 1;

    meshSizeX = (xBoundMax-xBoundMin)/nXs;
    meshSizeY = (yBoundMax-yBoundMin)/nYs;
    param.meshSizeX = meshSizeX;
    param.meshSizeY = meshSizeY;

    localArray = (double *)malloc(((int)nXs+1)*((int)nYs+1)*4*sizeof(double ));
    globalArray = (double *)malloc(((int)nXs+1)*((int)nYs+1)*4*sizeof(double ));
    
    for(i=0; i<nMeshPoints*4; i++){
        localArray[i] = 0.0;
        globalArray[i] = 0.0;
    }

    seed = time(0) + myDomain;
    E = potExpection;
    D = potVariance;

    sprintf(name, "%d.plt", myDomain);
    str = name;
    ofstream out1(name);

    for(i=0; i<nPotParams/nDomains; i++){
        uni[0] = randm(&seed); 
        uni[1] = randm(&seed); 

        A = sqrt((-2)*log(uni[0]));
        B = 2 * PI*uni[1];
        C = A*cos(B);
        potential = (E + C*D);   
        x0 = randm(&seed)*(xBoundMax-xBoundMin) + xBoundMin;
        y0 = randm(&seed)*(yBoundMax-yBoundMin) + yBoundMin;

        Coor2EleNum(&param, x0-radius, y0-radius, imin, jmin);
        Coor2EleNum(&param, x0+radius, y0+radius, imax, jmax);

        out1 << x0 << " " << y0 << " " << potential << endl;
        
        for(j=imin; j<imax; j++){
            for(k=jmin; k<jmax; k++){
                xi = (((double)j)*param.meshSizeX + param.xBoundMin);
                yi = (((double)k)*param.meshSizeY + param.yBoundMin);
                dx = x0 - xi;
                dy = y0 - yi;
                if(dx*dx+dy*dy > radius*radius)continue;
            
                ii = j;
                jj = k;
                if(ii < 0) ii += param.meshNumX;
                if(ii >= param.meshNumX) ii -= param.meshNumX;        
                if(jj < 0) jj += param.meshNumY;
                if(jj >= param.meshNumY) jj -= param.meshNumY;        
            
                f = -potential*(yi-y0)/sigma/sigma*pow(e,-((xi-x0)*(xi-x0)+(yi-y0)*(yi-y0))/(2.0*sigma*sigma));
                fx = -(xi-x0)*f/(sigma*sigma);
                fy = -(yi-y0)*f/(sigma*sigma) + f/(yi-y0);
                fxy = (xi-x0)*(yi-y0)*f/(sigma*sigma*sigma*sigma) - (xi-x0)*f/(yi-y0)/(sigma*sigma);

                if(isnan(fx)) fx = 0.0;
                if(isnan(fy)) fy = 0.0;
                if(isnan(fxy)) fxy = 0.0;

                localArray[0+4*jj+4*numY1*ii] = f;
                localArray[1+4*jj+4*numY1*ii] = fx;
                localArray[2+4*jj+4*numY1*ii] = fy;
                localArray[3+4*jj+4*numY1*ii] = fxy;
                
            }
        }
    }
    out1.close(); 

    for(i=0; i<(int)nXs+1; i++){
        localArray[0+4*(int)nYs+4*numY1*i] = localArray[0+4*numY1*i];
        localArray[1+4*(int)nYs+4*numY1*i] = localArray[1+4*numY1*i];
        localArray[2+4*(int)nYs+4*numY1*i] = localArray[2+4*numY1*i];
        localArray[3+4*(int)nYs+4*numY1*i] = localArray[3+4*numY1*i];
    }
    for(i=0; i<(int)nYs+1; i++){
        localArray[0+4*i+4*numY1*((int)nXs)] = localArray[0+4*i];
        localArray[1+4*i+4*numY1*((int)nXs)] = localArray[1+4*i];
        localArray[2+4*i+4*numY1*((int)nXs)] = localArray[2+4*i];
        localArray[3+4*i+4*numY1*((int)nXs)] = localArray[3+4*i];
    }

#if 0

	MPI_Allreduce(localArray, globalArray, nMeshPoints*4, MPI_DOUBLE,
			   MPI_SUM, MPI_COMM_WORLD);
#else
	MPI_Reduce(localArray, globalArray, nMeshPoints*4, MPI_DOUBLE,
			   MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    free(localArray);

    if(myDomain == 0){
        ofstream outfile("Tec.plt");
        outfile << numY1*numX1*4 << " " << meshSizeX <<
                " " << meshSizeY << " " << numX1 << " " << numY1 <<endl;
        for(i=0; i<numX1; i++){
            for(j=0; j<numY1; j++){
                outfile << scientific << globalArray[0+4*j+4*numY1*i] << " ";
                outfile << scientific << globalArray[1+4*j+4*numY1*i] << " ";
                outfile << scientific << globalArray[2+4*j+4*numY1*i] << " ";
                outfile << scientific << globalArray[3+4*j+4*numY1*i] << endl;
            } 
        }
        outfile.close();
#if 0
        meshSizeX /= 40.0;
        meshSizeY /= 40.0;

        nXs = floor((xBoundMax - xBoundMin)/meshSizeX);
        nYs = floor((yBoundMax - yBoundMin)/meshSizeY);

        meshSizeX = (xBoundMax - xBoundMin)/nXs;
        meshSizeY = (yBoundMax - yBoundMin)/nYs;

        ofstream out2("refine.plt");
        for(i=0; (double)i<nXs+1; i++){
            for(j=0; (double)j<nYs+1; j++){
                x0 = ((double)i)*meshSizeX+xBoundMin;
                y0 = ((double)j)*meshSizeY+yBoundMin;
                Coor2EleNum(&param, x0, y0, imin, jmin);
                imax = imin + 1;
                jmax = jmin + 1;
                if(imax > param.meshNumX) imax -= param.meshNumX;
                if(jmax > param.meshNumY) jmax -= param.meshNumY;

                l1 = param.meshSizeX;
                l2 = param.meshSizeY;

                f = Interpolation(globalArray[0+jmin*4+imin*4*numY1], globalArray[1+jmin*4+imin*4*numY1]*l1,
                              globalArray[2+jmin*4+imin*4*numY1]*l2, globalArray[3+jmin*4+imin*4*numY1]*l1*l2,

                              globalArray[0+jmax*4+imin*4*numY1], globalArray[1+jmax*4+imin*4*numY1]*l1,
                              globalArray[2+jmax*4+imin*4*numY1]*l2, globalArray[3+jmax*4+imin*4*numY1]*l1*l2,

                              globalArray[0+jmin*4+imax*4*numY1], globalArray[1+jmin*4+imax*4*numY1]*l1,
                              globalArray[2+jmin*4+imax*4*numY1]*l2, globalArray[3+jmin*4+imax*4*numY1]*l1*l2,

                              globalArray[0+jmax*4+imax*4*numY1], globalArray[1+jmax*4+imax*4*numY1]*l1,
                              globalArray[2+jmax*4+imax*4*numY1]*l2, globalArray[3+jmax*4+imax*4*numY1]*l1*l2,

                              (x0-(param.xBoundMin+((double)imin*param.meshSizeX)))/param.meshSizeX,
                              (y0-(param.yBoundMin+((double)jmin*param.meshSizeY)))/param.meshSizeY, bli);

                out2 << x0 << " " << y0 << " " << f << " " << bli <<endl;
            }
        }
        out2.close();
     
#endif   
    }

    free(globalArray);
    MPI_Finalize();

    return (0);
}
