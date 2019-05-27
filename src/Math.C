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

    return(0);   
}

int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2,
                           double *vec1, double *vec2)
{
    double  n1x = n1[0], n1y = n1[1], n1z = n1[2];
    double  n2x = n2[0], n2y = n2[1], n2z = n2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  dir[3], dir1[3], dir2[3], dot, dis, vec3[3], vec4[4], vec5[3], vec6[3];
    int     type;

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(sqrt(SQUARE(fabs(-(n1y*n2x) + n1x*n2y)) + 
       SQUARE(fabs(n1z*n2x - n1x*n2z)) + 
       SQUARE(fabs(-(n1z*n2y) + n1y*n2z))) < EPS1){
        if(fabs(n1x*(p1x - p2x) + n1y*(p1y - p2y) + n1z*(p1z - p2z)) < EPS2){
            VECTOR_COPY(vec1,n2);
            VECTOR_COPY(vec2,p2);
            return(PLANE_INTERSECTION);
        }

        VECTOR_COPY(vec1,n2);
        VECTOR_COPY(vec2,p2);
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

    cross(vec1, n1, vec3);
    if(Normal(vec3)< 1.0E-4){
        FatalMath("PlanePlaneIntersection in %s at %d\n", __FILE__, __LINE__);
    }
    NormalizeVec(vec3);

    type = LinePlaneIntersection(vec3, p1, n2, p2, vec4, vec2);
    if(type == NO_INTERSECTION){
        FatalMath("PlanePlaneIntersection in %s at %d\n", __FILE__, __LINE__);
    }

#if 0
    dir[0] = p2[0] - p1[0];
    dir[1] = p2[1] - p1[1];
    dir[2] = p2[2] - p1[2];
    if(Normal(dir) < 1.0E-2){
        VECTOR_COPY(vec2, p2);
        return(LINE_INTERSECTION);
    }
    NormalizeVec(dir);

    dot = DotProduct(dir, n1);
    dir1[0] = dir[0] - dot*n1[0];
    dir1[1] = dir[1] - dot*n1[1];
    dir1[2] = dir[2] - dot*n1[2];
    if(Normal(dir1) < 1.0E-4){
        VECTOR_COPY(vec2, p1);
        return(LINE_INTERSECTION);
    }

    NormalizeVec(dir1);

    dot = DotProduct(dir, n2);
    dir2[0] = dir[0] - dot*n2[0];
    dir2[1] = dir[1] - dot*n2[1];
    dir2[2] = dir[2] - dot*n2[2];
    if(Normal(dir2) < 1.0E-4){
        VECTOR_COPY(vec2, p2);
        return(LINE_INTERSECTION);
    }

    NormalizeVec(dir2);

    type = LineLineIntersection(dir1, p1, dir2, p2, vec5, vec2, &dis);

    if(fabs(dis) > 2){
        FormatVector(vec3, "p1");
        FormatVector(vec4, "p2");
        printf("Warning: PlanePlaneIntersection dis is %f\n", dis);
    }
#endif
    return(LINE_INTERSECTION);
} 

int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, 
                          double *vec1, double *vec2)
{
    double  dx = d[0], dy = d[1], dz = d[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  dis;

    VECTOR_ZERO(vec1);
    VECTOR_ZERO(vec2);

    if(fabs(dx*nx + dy*ny + dz*nz) < EPS1){
        VECTOR_COPY(vec1, d);
        PointPlaneIntersection(p1, n, p2, vec2, &dis);
        if(fabs(dis) < EPS2){
            VECTOR_COPY(vec2, p1);
            return(LINE_INTERSECTION);
        }else{ 
            return(NO_INTERSECTION);
        }
    }

    double  t = (-(nx*p1x) - ny*p1y - nz*p1z + nx*p2x + ny*p2y + nz*p2z)/
                (dx*nx + dy*ny + dz*nz);

    vec2[0] = p1x + t*dx;    
    vec2[1] = p1y + t*dy;    
    vec2[2] = p1z + t*dz;

    if(isnan(vec2[0]) || isinf(vec2[0]) ||
       isnan(vec2[1]) || isinf(vec2[1]) ||
       isnan(vec2[2]) || isinf(vec2[2])){
        printf(" LinePlaneIntersection: %f,%f,%f\n", 
                vec2[0], vec2[1], vec2[2]);
    }

    if(isnan(vec2[0]))vec2[0] = 0.0;
    if(isnan(vec2[1]))vec2[1] = 0.0;
    if(isnan(vec2[2]))vec2[2] = 0.0;
    
    if(isinf(vec2[0]) || isinf(vec2[1]) || isinf(vec2[2]) ){
        printf("d: {%f,%f,%f}\n", d[0], d[1], d[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("n: {%f,%f,%f}\n", n[0], n[1], n[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        exit(0);
    }
#if 0
    if(sqrt(SQUARE(vec2[2]) + SQUARE(vec2[2]) + SQUARE(vec2[2])) > 10000){
        printf("vec2: {%f,%f,%f}\n", vec2[0], vec2[1], vec2[2]);
        printf("d: {%f,%f,%f}\n", d[0], d[1], d[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("n: {%f,%f,%f}\n", n[0], n[1], n[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("PointPlaneIntersection\n");
    }
#endif

    return(POINT_INTERSECTION);    
}

int PointLineIntersection(double *p1, double *d, double *p2,
                         double *vec, double *dis){

    double  p1x = p1[0], p1y=p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y=p2[1], p2z = p2[2];
    double  dx = d[0], dy=d[1], dz = d[2];
    double  norm, norm2, invNorm, t;

    VECTOR_ZERO(vec);
    norm = sqrt(SQUARE(dx) + SQUARE(dy) + SQUARE(dz));
    invNorm = 1.0/norm;

    if(norm < EPS1){
        FatalMath("PointLineIntersection: line direction is zero\n");
        return(NO_INTERSECTION);
    }

    norm2 = sqrt(SQUARE(p1x-p2x) + SQUARE(p1y-p2y) + SQUARE(p1z-p2z));
    if(norm2 < EPS1){
        vec[0] = p2[0];
        vec[1] = p2[1];
        vec[2] = p2[2];
        *dis = 0.0;
        return(POINT_INTERSECTION);
    }

#if 0
    if(fabs(dx*(p1x-p2x) + dy*(p1y-p2y) + dz*(p1z-p2z)) > norm*norm2*(1-EPS1)){
        vec[0] = p1[0];
        vec[1] = p1[1];
        vec[2] = p1[2];
        *dis = 0.0;
        return(POINT_INTERSECTION);
    }
#endif
    t = (dx*(p1x - p2x) + dy*(p1y - p2y) + dz*(p1z - p2z))*(invNorm*invNorm);
    vec[0] = p2x + dx*t;
    vec[1] = p2y + dy*t;
    vec[2] = p2z + dz*t;

    *dis = sqrt(SQUARE(vec[0]-p1[0]) + SQUARE(vec[1]-p1[1]) + SQUARE(vec[2]-p1[2]) );

    if(isnan(*dis) || isinf(*dis)){
        printf(" PointLineIntersection: %f {%f,%f,%f}\n", 
               *dis, vec[0], vec[1], vec[2]);
        vec[0] = p1[0];
        vec[1] = p1[1];
        vec[2] = p1[2];
    }

    if(isnan(*dis))*dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;
#if 0
    if(*dis < EPS2 && sqrt(SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0])) > 10000){
        printf("vec: {%f,%f,%f}\n", vec[0], vec[1], vec[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("d: {%f,%f,%f}\n", d[0], d[1], d[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("PointLineIntersection");
    }
#endif
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

    if(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z) < EPS1*EPS1){
        FatalMath("LineLineIntersection: line direction is zero L1\n");
        return(NO_INTERSECTION);
    }
    if(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z) < EPS1*EPS1){
        FatalMath("LineLineIntersection: line direction is zero L2\n");
        return(NO_INTERSECTION);
    }

    if(fabs(1.0-fabs(d1x*d2x+d1y*d2y+d1z*d2z)) < EPS1){
        if(PointLineIntersection(p1, d2, p2, vec2, dis) == POINT_INTERSECTION){
            VECTOR_COPY(vec1, d1);
            *dis = fabs(*dis);
            return(LINE_INTERSECTION);
        }else{
            *dis = fabs(*dis);
            VECTOR_COPY(vec1, d2);
            VECTOR_COPY(vec2, p2);
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


    if(isnan(*dis) || isinf(*dis)){
        printf("LineLineIntersection: %f\n", *dis);
    }
    if(isnan(*dis)) *dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;

    if(fabs(*dis) < EPS2 && sqrt(SQUARE(vec2[0]-p1[0]) + SQUARE(vec2[0]-p1[0]) + SQUARE(vec2[0]-p1[0])) > 10000){
        printf("vec2: {%f,%f,%f}\n", vec2[0], vec2[1], vec2[2]);
        printf("d1: {%f,%f,%f}\n", d1[0], d1[1], d1[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("d2: {%f,%f,%f}\n", d2[0], d2[1], d2[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("LineLineIntersection\n");
    }

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

    if(isnan(*dis) || isinf(*dis)){
        printf(" PointPlaneIntersection: %f\n", *dis);
    }

    if(isnan(*dis)) *dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;

    vec[0] = (SQUARE(ny)*p1x + SQUARE(nz)*p1x + SQUARE(nx)*p2x + 
             nx*ny*(-p1y + p2y) + nx*nz*(-p1z + p2z))/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    vec[1] = (SQUARE(nx)*p1y + SQUARE(nz)*p1y + nx*ny*(-p1x + p2x) + 
             SQUARE(ny)*p2y + ny*nz*(-p1z + p2z))/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    vec[2] = (SQUARE(nx)*p1z + SQUARE(ny)*p1z + nx*nz*(-p1x + p2x) + 
             ny*nz*(-p1y + p2y) + SQUARE(nz)*p2z)/
             (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));

    temVec[0] = p1[0] - vec[0];
    temVec[1] = p1[1] - vec[1];
    temVec[2] = p1[2] - vec[2];
    NormalizeVec(temVec);
    if(DotProduct(temVec, n) < 0){
        *dis *= -1.0;
    }
    
    if(fabs(*dis) < EPS2 && sqrt(SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0]) + SQUARE(vec[0]-p1[0])) > 10000){
        printf("vec: {%f,%f,%f}\n", vec[0], vec[1], vec[2]);
        printf("p1: {%f,%f,%f}\n", p1[0], p1[1], p1[2]);
        printf("n: {%f,%f,%f}\n", n[0], n[1], n[2]);
        printf("p2: {%f,%f,%f}\n", p2[0], p2[1], p2[2]);
        printf("PointPlaneIntersection\n");
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
    if(isnan(*dis) || isinf(*dis)){
        printf(" PointPointIntersection: %f\n", *dis);
    }
    if(isnan(*dis)) *dis = 0.0;
    if(isinf(*dis)) *dis = 0.0;
    
    if(*dis < EPS2){
        return(POINT_INTERSECTION);
    }else{
        return(NO_INTERSECTION);
    }
}
