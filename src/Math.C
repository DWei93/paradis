/***************************************************************************
 *
 *      Module:       Math.C
 *      Description:  Contains various math functions
 *
 *      Includes:   
 *                  Arranege (int& n);
 *                  Arranege (int& n, int& m);
 *
 **************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "Math.h"

#include "Home.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

void FatalMath(const char *format, ...) 
{
        char    msg[512];
        va_list args;

        va_start(args, format);
        vsnprintf(msg, sizeof(msg)-1, format, args);
        msg[sizeof(msg)-1] = 0;
        va_end(args);
        printf("Fatal - Math: %s\n", msg);

#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        exit(1);
}

int Arrange(int n){
    if(n <= 0)return(1);

    int a = 1, i;
    for(i=n; i>0; i--){
        a *= i;
    }
    return(a);
}

int Arrange(int n, int m){
    if(n <= 0 || m <= 0)return(1);

    int nSubm = n - m;
    return(Arrange(n)/Arrange(nSubm));
}

int Combine(int n, int m){
    if(n <= 0 || m <= 0)return(1);
    
    return(Arrange(n,m)/Arrange(m));
}

#if 0
void Perm(int deep, int n, int m, bool *isUsed, int **results){
    int i;
    if (deep > n-m)return;

    for(i=0; i<n; i++){
        if(!isUsed[i]){
            isUsed[i] = 1;
            results[deep] = i;
            Perm(deep+1, n, m, isUsed, results);
            isUsed[i] = 0;
        }
    }
}
#endif

void GetArrange(int n, int m, int **arra){
    bool        *isUsed;
    double      *results;

    isUsed = (bool *)calloc(n,sizeof(bool));
    results = (double *)calloc(Arrange(n,m), sizeof(double));
}

void GetCombination(int n, int m, int **comb){

    bool *isUsed;
    isUsed = (bool *)calloc(n,sizeof(bool));
    
    
}

double Volume(double **p, int num){
    int         i, numComb;
    double      volume;
    
    if(num < 3)FatalMath("Volume: the shape has %d edges",num);
    numComb = Combine(num, 3);

//    for(i=0; i<; i++){
        
    
//    }

    return(0);   
}

int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2,
                           double *vec1, double *vec2)
{
    double  n1x = n1[0], n1y = n1[1], n1z = n1[2];
    double  n2x = n2[0], n2y = n2[1], n2z = n2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(sqrt(SQUARE(fabs(-(n1y*n2x) + n1x*n2y)) + 
       SQUARE(fabs(n1z*n2x - n1x*n2z)) + 
       SQUARE(fabs(-(n1z*n2y) + n1y*n2z))) < EPS1){
        if(fabs(n1x*(p1x - p2x) + n1y*(p1y - p2y) + n1z*(p1z - p2z)) < EPS2){
            VECTOR_COPY(vec1,n1);
            VECTOR_COPY(vec2,p1);
            return(PLANE_INTERSECTION);
        }
        return(NO_INTERSECTION);
    } 

    vec1[0] = -(n1z*n2y) + n1y*n2z;
    vec1[1] = n1z*n2x - n1x*n2z;
    vec1[2] = -(n1y*n2x) + n1x*n2y;
    NormalizeVec(vec1);

    if(sqrt(SQUARE(fabs(p1x - p2x)) + SQUARE(fabs(p1y - p2y)) + 
       SQUARE(fabs(p1z - p2z))) < EPS2){
        VECTOR_COPY(vec2,p2);
        return(LINE_INTERSECTION);
    }

    vec2[0] = ((-(SQUARE(n1y)*n2x) + n1x*n1y*n2y - n1z*(n1z*n2x - n1x*n2z))*
            (n2z*(n1x*p1x + n1y*p1y + n1z*p1z) - 
            n1z*(n2x*p2x + n2y*p2y + n2z*p2z)) + 
            (n1z*n2y - n1y*n2z)*((-(n1y*n2x) + n1x*n2y)*
            (n1x*p1x + n1y*p1y + n1z*p1z) + 
            (n1z*(n1z*(n2y*(p1x + p2x) - n2x*(p1y + p2y)) + 
            n1y*(-(n2z*(p1x + p2x)) + n2x*(p1z + p2z)) + 
            n1x*(n2z*(p1y + p2y) - n2y*(p1z + p2z))))/2.))/
            (n1z*(SQUARE(n1z)*(SQUARE(n2x) + SQUARE(n2y)) - 
            2*n1x*n1z*n2x*n2z - 2*n1y*n2y*(n1x*n2x + n1z*n2z) + 
            SQUARE(n1y)*(SQUARE(n2x) + SQUARE(n2z)) + 
            SQUARE(n1x)*(SQUARE(n2y) + SQUARE(n2z))));

    vec2[1] = (2*SQUARE(n1y)*(SQUARE(n2x) + SQUARE(n2z))*p1y + 
            SQUARE(n1x)*(2*n2x*n2y*(-p1x + p2x) + 2*SQUARE(n2y)*p2y + 
            SQUARE(n2z)*(p1y + p2y) + n2y*n2z*(-p1z + p2z)) + 
            SQUARE(n1z)*(n2x*n2y*(-p1x + p2x) + 
            SQUARE(n2x)*(p1y + p2y) + 
            2*n2y*(n2y*p2y + n2z*(-p1z + p2z))) + 
            n1y*n1z*(n2x*n2z*(p1x - p2x) + SQUARE(n2x)*(p1z - p2z) - 
            2*n2z*(n2y*(p1y + p2y) + n2z*(-p1z + p2z))) + 
            n1x*(-(n1z*(2*n2x*n2z*(p1y + p2y) + 
            n2y*(n2z*p1x + n2x*p1z - n2z*p2x - n2x*p2z))) + 
            n1y*(2*SQUARE(n2x)*(p1x - p2x) + 
            SQUARE(n2z)*(p1x - p2x) - 
            n2x*(2*n2y*(p1y + p2y) + n2z*(-p1z + p2z)))))/
            (2.*(SQUARE(n1z)*(SQUARE(n2x) + SQUARE(n2y)) - 
            2*n1x*n1z*n2x*n2z - 2*n1y*n2y*(n1x*n2x + n1z*n2z) + 
            SQUARE(n1y)*(SQUARE(n2x) + SQUARE(n2z)) + 
            SQUARE(n1x)*(SQUARE(n2y) + SQUARE(n2z))));

    vec2[2] = (2*SQUARE(n1z)*(SQUARE(n2x) + SQUARE(n2y))*p1z + 
            SQUARE(n1x)*(2*n2x*n2z*(-p1x + p2x) + n2y*n2z*(-p1y + p2y) + 
            2*SQUARE(n2z)*p2z + SQUARE(n2y)*(p1z + p2z)) + 
            SQUARE(n1y)*(n2x*n2z*(-p1x + p2x) + 
            SQUARE(n2x)*(p1z + p2z) + 
            2*n2z*(n2y*(-p1y + p2y) + n2z*p2z)) + 
            n1y*n1z*(n2x*n2y*(p1x - p2x) + SQUARE(n2x)*(p1y - p2y) + 
            2*n2y*(n2y*(p1y - p2y) - n2z*(p1z + p2z))) + 
            n1x*(-(n1y*(n2x*n2z*(p1y - p2y) + 
            n2y*(n2z*(p1x - p2x) + 2*n2x*(p1z + p2z)))) + 
            n1z*(2*SQUARE(n2x)*(p1x - p2x) + 
            SQUARE(n2y)*(p1x - p2x) + 
            n2x*(n2y*(p1y - p2y) - 2*n2z*(p1z + p2z)))))/
            (2.*(SQUARE(n1z)*(SQUARE(n2x) + SQUARE(n2y)) - 
            2*n1x*n1z*n2x*n2z - 2*n1y*n2y*(n1x*n2x + n1z*n2z) + 
            SQUARE(n1y)*(SQUARE(n2x) + SQUARE(n2z)) + 
            SQUARE(n1x)*(SQUARE(n2y) + SQUARE(n2z))));

    if(fabs(1-fabs(n1x)) < EPS1)vec2[0] = p1x;
    if(fabs(1-fabs(n2x)) < EPS1)vec2[0] = p2x;

    if(fabs(1-fabs(n1y)) < EPS1)vec2[1] = p1y;
    if(fabs(1-fabs(n2y)) < EPS1)vec2[1] = p2y;

    if(fabs(1-fabs(n1z)) < EPS1)vec2[2] = p1z;
    if(fabs(1-fabs(n2z)) < EPS1)vec2[2] = p2z;

    return(LINE_INTERSECTION);
} 

int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, 
                          double *vec1, double *vec2)
{
    double  dx = d[0], dy = d[1], dz = d[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(fabs(dx*nx + dy*ny + dz*nz) < EPS1*EPS1){
        if(fabs(nx*(p1x - p2x) + ny*(p1y - p2y) + nz*(p1z - p2z)) < EPS1*EPS2){
            VECTOR_COPY(vec1, d);
            VECTOR_COPY(vec2, p1);
            return(LINE_INTERSECTION);
        } 
        return(NO_INTERSECTION);
    }

    double  t = (-(nx*p1x) - ny*p1y - nz*p1z + nx*p2x + ny*p2y + nz*p2z)/
                (dx*nx + dy*ny + dz*nz);

    vec2[0] = p1x + t*dx;    
    vec2[1] = p1y + t*dy;    
    vec2[2] = p1z + t*dz;
    
    return(POINT_INTERSECTION);    
}

int PointLineIntersection(double *p1, double *d, double *p2,
                         double *vec, double *dis){

    double  p1x = p1[0], p1y=p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y=p2[1], p2z = p2[2];
    double  dx = d[0], dy=d[1], dz = d[2];

    VECTOR_ZERO(vec);
    if(SQUARE(dx) + SQUARE(dy) + SQUARE(dz) < EPS1*EPS1){
        FatalMath("PointLineIntersection: line direction is zero\n");
        return(NO_INTERSECTION);
    }
    
    vec[0] = p2x + (dx*(dx*(p1x - p2x) + dy*(p1y - p2y) + dz*(p1z - p2z)))/
             (SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
    vec[1] = p2y + (dy*(dx*(p1x - p2x) + dy*(p1y - p2y) + dz*(p1z - p2z)))/
             (SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
    vec[2] = (dz*(dx*(p1x - p2x) + dy*(p1y - p2y) + dz*(p1z - p2z)))/
             (SQUARE(dx) + SQUARE(dy) + SQUARE(dz)) + p2z;

    *dis = sqrt((SQUARE(dz)*(SQUARE(p1x) + SQUARE(p1y) - 2*p1x*p2x + 
       SQUARE(p2x) - 2*p1y*p2y + SQUARE(p2y)) - 
       2*dy*(p1y - p2y)*(dx*(p1x - p2x) + dz*(p1z - p2z)) - 
       2*dx*dz*(p1x - p2x)*(p1z - p2z) + 
       SQUARE(dy)*(SQUARE(p1x) + SQUARE(p1z) - 2*p1x*p2x + 
       SQUARE(p2x) - 2*p1z*p2z + SQUARE(p2z)) + 
       SQUARE(dx)*(SQUARE(p1y) + SQUARE(p1z) - 2*p1y*p2y + 
       SQUARE(p2y) - 2*p1z*p2z + SQUARE(p2z)))/
       (SQUARE(dx) + SQUARE(dy) + SQUARE(dz)));
    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}

int LineLineIntersection(double *d1, double *p1, double *d2, double *p2,
                         double *vec1, double *vec2, double *dis){
    double  d1x = d1[0], d1y = d1[1], d1z = d1[2];
    double  d2x = d2[0], d2y = d2[1], d2z = d2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);
    if(sqrt(SQUARE(fabs(-(d1y*d2x) + d1x*d2y)) + 
            SQUARE(fabs(d1z*d2x - d1x*d2z)) + 
            SQUARE(fabs(-(d1z*d2y) + d1y*d2z))) < EPS1){
        if(PointLineIntersection(p1, d2, p2, vec2, dis) == POINT_INTERSECTION){
            VECTOR_COPY(vec1, d1);
            *dis = fabs(*dis);
            return(LINE_INTERSECTION);
        }else{
            *dis = fabs(*dis);
            return(NO_INTERSECTION);
        }
    }

/*
 *  return the closed point on the first line
 */
    vec1[0] = p1x - (d1x*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/
              (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
              4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
              (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z)));        
    vec1[1] = p1y - (d1y*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/
              (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
              4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
              (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z)));
    vec1[2] = p1z - (d1z*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/
              (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
              4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
              (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))); 

/*
 *  return the closed point on the second line
 */
    vec2[0] = p2x + (d2x*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/
              (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
              2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
              SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
              SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

    vec2[1] = p2y + (d2y*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/
              (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
              2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
              SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
              SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

    vec2[2] = p2z + (d2z*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/
              (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
              2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
              SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
              SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

    *dis = sqrt(SQUARE(d1z*(d2y*(p1x - p2x) + d2x*(-p1y + p2y)) + 
           d1y*(-(d2z*p1x) + d2x*p1z + d2z*p2x - d2x*p2z) + 
           d1x*(d2z*p1y - d2y*p1z - d2z*p2y + d2y*p2z))/
           (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
           2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
           SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
           SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z))));

    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}


int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis)
{
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  temVec[3];

    if(nx*nx + ny*ny + nz*nz < EPS1*EPS1){
        printf("PointPlaneIntersection: normal is close to zero.\n");
        return(0);
    }

    *dis = sqrt(SQUARE(nx*p1x + ny*p1y + nz*p1z - nx*p2x - ny*p2y - nz*p2z)/
           (SQUARE(nx) + SQUARE(ny) + SQUARE(nz)));

    vec[0] = (SQUARE(ny)*p1x + SQUARE(nz)*p1x + SQUARE(nx)*p2x + 
             nx*ny*(-p1y + p2y) + nx*nz*(-p1z + p2z))/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    vec[1] = (SQUARE(nx)*p1y + SQUARE(nz)*p1y + nx*ny*(-p1x + p2x) + 
             SQUARE(ny)*p2y + ny*nz*(-p1z + p2z))/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    vec[2] = (SQUARE(nx)*p1z + SQUARE(ny)*p1z + nx*nz*(-p1x + p2x) + 
             ny*nz*(-p1y + p2y) + SQUARE(nz)*p2z)/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));

    temVec[0] = vec[0] - p1[0];
    temVec[1] = vec[1] - p1[1];
    temVec[2] = vec[2] - p1[2];
    NormalizeVec(temVec);
    if(DotProduct(temVec, n) < 0){
        *dis *= -1.0;
    }
    
    if(fabs(*dis) < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}

int PointPointIntersection(double *p1, double *p2, double *vec, double *dis){
    vec[0] = p2[0] - p1[0];
    vec[1] = p2[1] - p1[1];
    vec[2] = p2[2] - p1[2];

    *dis = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    
    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}
