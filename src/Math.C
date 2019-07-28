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
#include <vector>
#include "Math.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef GSL
#include <gsl/gsl_sf_bessel.h>
#endif
/*
 *      Define some inline operations for handling 3-component vectors
 *
 * VECTOR_ADD:  add the components of the second vector to the first
 * VECTOR_COPY: copy the components of the second vector to the first
 * VECTOR_ZERO: Zero out the components of the vector
 */
#define VECTOR_ADD(a,b)  {(a)[0] += (b)[0]; (a)[1] += (b)[1]; (a)[2] += (b)[2];}
#define VECTOR_MINUS(a,b)  {(a)[0] -= (b)[0]; (a)[1] -= (b)[1]; (a)[2] -= (b)[2];}
#define VECTOR_COPY(a,b) {(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];}
#define VECTOR_ZERO(a)   {(a)[0] = 0; (a)[1]  = 0; (a)[2] = 0;}
#define VECTOR_MULTIPLY(a,x)	{(a)[0] *= (x);    (a)[1] *= (x);    (a)[2] *= (x);}
#define DotProduct(vec1,vec2)       \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))

/* normal of vector*/
extern double  Normal(double a[3]);
extern double  Normal(std::vector<double>& a);

/* Normalize the vector*/
extern void   Normalize(double *ax, double *ay, double *az);
extern void   NormalizeVec(double vec[3]);

/* cross production of vectors*/
extern void   cross(double a[3], double b[3], double c[3]);

/* Abort the program */
extern void Fatal(const char *format, ...);

#ifdef GSL
void GSLTest()
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  printf ("J0(%g) = %.18e\n", x, y);
  return ;
}
#endif


/** Check if the variable is numerical. */
inline bool Numerical(double v){return (!isnan(v) && !isinf(v));}
inline bool Numerical(double *vec){
    return (!isnan(vec[0]) && !isinf(vec[0]) &&
            !isnan(vec[1]) && !isinf(vec[1]) &&
            !isnan(vec[2]) && !isinf(vec[2]));
}

/** Determine the vector if is nomalized */
inline bool IsNormalized(double *vec, double error){
    double dot=DotProduct(vec, vec);
    return (dot>error && dot <2-error);
}
inline bool IsNormalized(double x, double y, double z, double error){
    double vec[3]={x,y,z};
    return IsNormalized(vec,error);
}

/** Determine if two directions are in the same direction. */
inline bool SameDirection(double *vec1, double *vec2, double error){
    double dot=DotProduct(vec1, vec2);
    return (dot>error && dot <2-error);
}
inline bool SameDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double vec1[3]={x1,y1,z1}, vec2[3]={x2,y2,z2};
    return SameDirection(vec1, vec2, error);
}

/** Determine if the two directions are aligned. */
inline bool AlignedDirection(double *vec1, double *vec2, double error){
    double dot=fabs(DotProduct(vec1, vec2));
    return (dot>error && dot <2-error);
}
inline bool AlignedDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double vec1[3]={x1,y1,z1}, vec2[3]={x2,y2,z2};
    return AlignedDirection(vec1, vec2, error);
}

/** Determine whether two directions are perpendicular to each other. */
inline bool PerpendicularDirection(double *vec1, double *vec2, double error){
    double dot=fabs(DotProduct(vec1, vec2));
    return dot<error;
}
inline bool PerpendicularDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double vec1[3]={x1,y1,z1}, vec2[3]={x2,y2,z2};
    return PerpendicularDirection(vec1, vec2, error);
}

/** Find intersection of two planes, plane2 is the default reference plane.
 *  \param[in] n1[3] normlized normal of plane1
 *  \param[in] p1[3] a point on plane1
 *  \param[in] n2[3] normlized normal of plane2
 *  \param[in] p2[3] a point on plane2
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \param[out] vec1[3] normalized vector that stores intersection info (if has)
 *  \param[out] vec2[3] second vector that stores intersection info (if has)
 *  \return intersection type: PLANE_INTERSECTION, NO_INTERSECTION and LINE_INTERSECTION */
int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2, double *vec1, double *vec2, double err1, double err2)
{
    double  n1x = n1[0], n1y = n1[1], n1z = n1[2];
    double  n2x = n2[0], n2y = n2[1], n2z = n2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  dir[3], dir1[3], dir2[3], dot, dis, vec3[3], vec4[4], vec5[3], vec6[3];
    int     type;

#ifdef DEBUG
    if(!IsNormalized(n1))Fatal("%s: n1 (%g,%g,%f) is not nomalized", __func__, n1x, n1y, n1z);
    if(!IsNormalized(n2))Fatal("%s: n2 (%g,%g,%f) is not nomalized", __func__, n2x, n2y, n2z);
#endif
    // if the two normals are aligned with each other, the intersection can not be line
    if(AlignedDirection(n1,n2,1-err1)){
        VECTOR_COPY(vec1,n2); VECTOR_COPY(vec2,p2);
        if(fabs(n1x*(p1x - p2x) + n1y*(p1y - p2y) + n1z*(p1z - p2z)) < err2) return(PLANE_INTERSECTION);
        else return(NO_INTERSECTION);
    } 

    // the direction of intersection line is the cross vector of the two normals
    vec1[0] = -(n1z*n2y) + n1y*n2z;
    vec1[1] = n1z*n2x - n1x*n2z;
    vec1[2] = -(n1y*n2x) + n1x*n2y;
    NormalizeVec(vec1);

    // if the two points are too closed, point2 is chosed to be the point on intersection line
    if(SQUARE(p1x - p2x) + SQUARE(p1y - p2y) + SQUARE(p1z - p2z) < err2*err2){
        VECTOR_COPY(vec2,p2); return(LINE_INTERSECTION);
    }

    // the glide direction from point1 to plane2 is defined by the cross production of normal1 and intersection line direction
    cross(vec1, n1, vec3);
#ifdef DEBUG 
    if(Normal(vec3)< 1.0E-4)Fatal("PlanePlaneIntersection in %s at %d\n", __FILE__, __LINE__);
#endif
    NormalizeVec(vec3);

    // the key point of inetrsection line can be defined by the intersection point of glide line and plane2
    type = LinePlaneIntersection(vec3, p1, n2, p2, vec4, vec2, err1, err2);
#ifdef DEBUG
    if(type == NO_INTERSECTION)Fatal("PlanePlaneIntersection in %s at %d\n", __FILE__, __LINE__);
#endif
    return(LINE_INTERSECTION);
} 

/** Find intersection of line and plane
 *  \param[in] d[3] direction of line
 *  \param[in] p1[3] a point on line
 *  \param[in] n[3] normlized normal of plane
 *  \param[in] p2[3] a point on plane
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \param[out] vec1[3] normalized vector that stores intersection info (if has)
 *  \param[out] vec2[3] second vector that stores intersection info (if has)
 *  \return intersection type: POINT_INTERSECTION, NO_INTERSECTION and LINE_INTERSECTION */
int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, double *vec1, double *vec2, double err1, double err2)
{
    double  dx = d[0], dy = d[1], dz = d[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  dis;

#ifdef DEBUG    
    if(!IsNormalized(d))Fatal("%s: line direction (%g,%g,%g) is not normalized.", __func__, dx, dy, dz);
    if(!IsNormalized(n))Fatal("%s: normal of plane (%g,%g,%g) is not normalized.", __func__, nx, ny, nz);
#endif
    // if the line direction is perpendicular to the normal of plane, the intersection can not be a point
    if(PerpendicularDirection(d, n, err1)){
        PointPlaneIntersection(p1, n, p2, vec2, &dis);
        VECTOR_COPY(vec1, d); VECTOR_COPY(vec2, p1);
        if(fabs(dis) < err2) return(LINE_INTERSECTION);
        else return(NO_INTERSECTION);
    }

    double t = (-(nx*p1x) - ny*p1y - nz*p1z + nx*p2x + ny*p2y + nz*p2z)/
                (dx*nx + dy*ny + dz*nz);

#ifdef DEBUG
    if(isinf(t) || isnan(t))Fatal("%s: %f",__func__,t);
#endif
    vec2[0] = p1x + t*dx;    
    vec2[1] = p1y + t*dy;    
    vec2[2] = p1z + t*dz;

#ifdef DEBUG
    if(isnan(vec2[0]) || isinf(vec2[0]) || isnan(vec2[1]) || isinf(vec2[1]) ||
       isnan(vec2[2]) || isinf(vec2[2])){
        Fatal("%s: %f,%f,%f\n",__func__,vec2[0], vec2[1], vec2[2]);
    }
#endif
    return(POINT_INTERSECTION);    
}

/** Find intersection of point and line, the line is reference
 *  \param[in] p1[3] coordinates of the point
 *  \param[in] d[3] direction of the line
 *  \param[in] p2[3] a point on the line
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \param[out] vec[3] coordinates that stores intersection info (if has)
 *  \param[out] *dis distance of point 1 and projection point
 *  \return intersection type: POINT_INTERSECTION and NO_INTERSECTION */
int PointLineIntersection(double *p1, double *d, double *p2, double *vec, double *dis, double err1, double err2){

    double  p1x = p1[0], p1y=p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y=p2[1], p2z = p2[2];
    double  dx = d[0], dy=d[1], dz = d[2];
    double  normal2, t;

#ifdef DEBUG
    if(!IsNormalized(d))Fatal("%s: direction (%g,%g,%g) is not normalized", __func__, dx,dy,dz);
#endif
    // if the two points are too close, point 2 is identified as the intersection point.
    normal2 = SQUARE(p1x-p2x) + SQUARE(p1y-p2y) + SQUARE(p1z-p2z);
    if(normal2 < err2*err2){
        VECTOR_COPY(vec, p2); *dis = sqrt(normal2);
        if(isnan(*dis))*dis=0.0;
        return(POINT_INTERSECTION);
    }

    normal2=dx*dx+dy*dy+dz*dz;
    t = (dx*(p1x - p2x) + dy*(p1y - p2y) + dz*(p1z - p2z))/normal2;
#ifdef DEBUG
    if(isinf(t) || isnan(t))Fatal("%s: %f",__func__,t);
#endif
    vec[0] = p2x + dx*t;
    vec[1] = p2y + dy*t;
    vec[2] = p2z + dz*t;

    *dis = SQUARE(vec[0]-p1[0]) + SQUARE(vec[1]-p1[1]) + SQUARE(vec[2]-p1[2]);

#ifdef DEBUG
    if(isnan(*dis) || isinf(*dis)) Fatal("%s: %f {%f,%f,%f}\n", __func__,*dis, vec[0], vec[1], vec[2]);
#endif
    
    *dis = sqrt(*dis); if(isnan(*dis))*dis=0.0;
    if(*dis < err2) return(POINT_INTERSECTION);
    else return(NO_INTERSECTION);
}

/** Find intersection of two lines, the second line is reference
 *  \param[in] d1[3] direction of line 1
 *  \param[in] p1[3] point on line 1
 *  \param[in] d2[3] direction of line 2
 *  \param[in] p2[3] point on line 2
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \param[out] vec1[3] intersection info. 1
 *  \param[out] vec2[3] intersection info. 2
 *  \param[out] *dis distance of the two lines
 *  \return intersection type: LINE_INTERSECTION, POINT_INTERSECTION, NO_INTERSECTION */
int LineLineIntersection(double *d1, double *p1, double *d2, double *p2, 
                         double *vec1, double *vec2, double *dis, double err1, double err2){
    double  d1x = d1[0], d1y = d1[1], d1z = d1[2];
    double  d2x = d2[0], d2y = d2[1], d2z = d2[2];
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];

#ifdef DEBUG
    if(!IsNormalized(d1))Fatal("%s: direction (%g,%g,%g) is not normalized", __func__, d1x,d1y,d1z);
    if(!IsNormalized(d2))Fatal("%s: direction (%g,%g,%g) is not normalized", __func__, d2x,d2y,d2z);
#endif

    // if the two lines are aligned with each other, the intersection can not be a point
    if(AlignedDirection(d1,d2,1-err1)){
        PointLineIntersection(p1, d2, p2, vec2, dis, err1, err2);
        if(*dis < err2){
            // If two parallel lines are within the err2, the intersection is line 2
            VECTOR_COPY(vec1, d2);
            return(LINE_INTERSECTION);
        }else{
            *dis = fabs(*dis);
            VECTOR_COPY(vec1, d2);
            return(NO_INTERSECTION);
        }
    }

    double denominator;
    // return closet point to line 2 on line 1 as intersection info. 1
    denominator = (4*SQUARE(d1x*d2x + d1y*d2y + d1z*d2z) - 
                  4*(SQUARE(d1x) + SQUARE(d1y) + SQUARE(d1z))*
                  (SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z)));        
#ifdef DEBUG
    if(isnan(denominator) || isinf(denominator))Fatal("%s: the denominator is %f",__func__,denominator);
#endif

    vec1[0] = p1x - (d1x*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/denominator;

    vec1[1] = p1y - (d1y*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/denominator;

    vec1[2] = p1z - (d1z*(-4*(SQUARE(d2x) + SQUARE(d2y) + SQUARE(d2z))*
              (d1x*(p1x - p2x) + d1y*(p1y - p2y) + d1z*(p1z - p2z)) + 
              4*(d1x*d2x + d1y*d2y + d1z*d2z)*
              (d2x*(p1x - p2x) + d2y*(p1y - p2y) + d2z*(p1z - p2z))))/denominator;

#ifdef DEBUG
    if(!Numerical(vec1))Fatal("%s: vec1 (%g,%g,%g) is not numerical", __func__, vec1[0], vec1[1], vec1[2]);
#endif
    // return closet point to line 1 on line 2 as intersection info. 2
    denominator = (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
                  2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
                  SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
                  SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));
#ifdef DEBUG
    if(isnan(denominator) || isinf(denominator))Fatal("%s: the denominator is %f",__func__,denominator);
#endif

    vec2[0] = p2x + (d2x*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/denominator;

    vec2[1] = p2y + (d2y*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/denominator;

    vec2[2] = p2z + (d2z*(SQUARE(d1z)*
              (d2x*p1x + d2y*p1y - d2x*p2x - d2y*p2y) + 
              SQUARE(d1y)*(d2x*(p1x - p2x) + d2z*(p1z - p2z)) + 
              d1x*d1z*(-(d2z*p1x) - d2x*p1z + d2z*p2x + d2x*p2z) + 
              SQUARE(d1x)*(d2y*p1y + d2z*p1z - d2y*p2y - d2z*p2z) + 
              d1y*(d1x*(-(d2y*p1x) - d2x*p1y + d2y*p2x + d2x*p2y) + 
              d1z*(-(d2z*p1y) - d2y*p1z + d2z*p2y + d2y*p2z))))/denominator;
#ifdef DEBUG
    if(!Numerical(vec2))Fatal("%s: vec1 (%g,%g,%g) is not numerical", __func__, vec2[0], vec2[1], vec2[2]);
#endif

    *dis = SQUARE(d1z*(d2y*(p1x - p2x) + d2x*(-p1y + p2y)) + 
           d1y*(-(d2z*p1x) + d2x*p1z + d2z*p2x - d2x*p2z) + 
           d1x*(d2z*p1y - d2y*p1z - d2z*p2y + d2y*p2z))/
           (SQUARE(d1z)*(SQUARE(d2x) + SQUARE(d2y)) - 
           2*d1x*d1z*d2x*d2z - 2*d1y*d2y*(d1x*d2x + d1z*d2z) + 
           SQUARE(d1y)*(SQUARE(d2x) + SQUARE(d2z)) + 
           SQUARE(d1x)*(SQUARE(d2y) + SQUARE(d2z)));

#ifdef DEBUG
    if(isnan(*dis) || isinf(*dis))Fatal("%s: *dis=%f",__func__, *dis);
#endif

    *dis = sqrt(*dis); if(isnan(*dis))*dis=0.0;
    if(*dis < err2) return(POINT_INTERSECTION);
    else return(NO_INTERSECTION);
}


/** Find intersection of point and plane, the plane is reference
 *  \param[in] p1[3] coordinates of the point
 *  \param[in] n[3] nromal of the line
 *  \param[in] p2[3] a point on the plane
 *  \param[out] vec[3] coordinates that stores intersection info (if has)
 *  \param[out] *dis distance of point 1 and projection point
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \return intersection type: POINT_INTERSECTION and NO_INTERSECTION */
int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis, double err1, double err2)
{
    double  p1x = p1[0], p1y = p1[1], p1z = p1[2];
    double  p2x = p2[0], p2y = p2[1], p2z = p2[2];
    double  nx = n[0], ny = n[1], nz = n[2];
    double  temVec[3], normal2;

#ifdef DEBUG
    if(!IsNormalized(n))Fatal("%s: the plane normal (%g,%g,%g) is not normalized", __func__, nx, ny, nz);
#endif

    normal2 = (SQUARE(nx) + SQUARE(ny) + SQUARE(nz));
    *dis = SQUARE(nx*p1x + ny*p1y + nz*p1z - nx*p2x - ny*p2y - nz*p2z)/normal2;

#ifdef DEBUG
    if(!Numerical(*dis))Fatal("%s: *dis=%f",__func__,*dis);
#endif 
    
    *dis = sqrt(*dis); if(isnan(*dis))*dis=0.0;

    vec[0] = (SQUARE(ny)*p1x + SQUARE(nz)*p1x + SQUARE(nx)*p2x + 
             nx*ny*(-p1y + p2y) + nx*nz*(-p1z + p2z))/normal2;
    vec[1] = (SQUARE(nx)*p1y + SQUARE(nz)*p1y + nx*ny*(-p1x + p2x) + 
             SQUARE(ny)*p2y + ny*nz*(-p1z + p2z))/normal2;
    vec[2] = (SQUARE(nx)*p1z + SQUARE(ny)*p1z + nx*nz*(-p1x + p2x) + 
             ny*nz*(-p1y + p2y) + SQUARE(nz)*p2z)/normal2;
#ifdef DEBUG
    if(!Numerical(vec))Fatal("the intersection point (%g,%g,%g)", vec[0],vec[1],vec[2]);
#endif

    temVec[0] = p1x - vec[0];
    temVec[1] = p1y - vec[1];
    temVec[2] = p1z - vec[2];
    NormalizeVec(temVec);
    if(DotProduct(temVec, n) < 0)*dis *= -1.0;
    
    if(fabs(*dis) < err2) return(POINT_INTERSECTION);
    else return(NO_INTERSECTION);
}

/** Find intersection of two points, the second point is reference
 *  \param[in] p1[3] coordinates of the point
 *  \param[in] p2[3] a point on the plane
 *  \param[out] vec[3] coordinates that stores intersection info (if has)
 *  \param[out] *dis distance of point 1 and projection point
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \return intersection type: POINT_INTERSECTION and NO_INTERSECTION */
int PointPointIntersection(double *p1, double *p2, double *vec, double *dis, double err1, double err2){

    VECTOR_COPY(vec, p2);
    *dis = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
#ifdef DEBUG
    if(!Numerical(*dis))Fatal("%s: distance %f is not numerical", *dis);
#endif
    *dis = sqrt(*dis); if(isnan(*dis))*dis=0.0;

    if(*dis < err2) return(POINT_INTERSECTION);
    else return(NO_INTERSECTION);
}

/** Find a intersection of segment and plane (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] nDir[3] normal of plane
 *  \param[in] point[3] coordinate of point on pkane
 *  \param[out] t normalized position of intersection pointi
 *  \return intersection type : NO_INTERSECTION and POINT_INTERSECTION */ 
int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double &t){

    double dot1 = nDir[0]*(p1[0]-point[0]) + nDir[1]*(p1[1]-point[1]) + nDir[2]*(p1[2]-point[2]);
    double dot2 = nDir[0]*(p2[0]-point[0]) + nDir[1]*(p2[1]-point[1]) + nDir[2]*(p2[2]-point[2]);
    
#ifdef DEBUG
    if(!IsNormalized(nDir))Fatal("%s: normal (%g,%g,%g) is not normalized");
    if(!Numerical(dot1))Fatal("%s: dot1=%f",__func__,dot1);
    if(!Numerical(dot2))Fatal("%s: dot2=%f",__func__,dot2);
#endif
    if(dot1*dot2 > 0)return(NO_INTERSECTION);
    double tot = fabs(dot1)+fabs(dot2);
    if(tot > 1.0E-80) t = fabs(dot1)/tot;
    else t = 0.0;
#ifdef DEBUG
    if(!Numerical(t))Fatal("%s: t=%f", __func__, t);
#endif
    return(POINT_INTERSECTION);
}

/** Find a intersection of segment and concave polygon (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] edges[6*nEdges] planes that containning concave polygon edges,
 *             each plane has normals and coodinates of point on the plane
 *  \param[in] nEdges # of edges of the concave polygon
 *  \param[out] interPos[3] the position of intersection point
 *  \param[out] interEdgeID index of edge which intersects with the segment
 *  \return intersection type : NO_INTERSECTION and POINT_INTERSECTION */ 
int SegmentConcavePolygonIntersection(double *p1, double *p2, double *edge, int nEdges, double *interPos, int &interEdgeID){
    int i, type; double t;
    for(i=0; i<nEdges; i++){
        if((type=SegmentPlaneIntersection(p1, p2, &(edge[i*6]), &(edge[i*6+3]), t))==POINT_INTERSECTION){
            interPos[0] = (1.0-t)*p1[0]+t*p2[0];
            interPos[1] = (1.0-t)*p1[1]+t*p2[1];
            interPos[2] = (1.0-t)*p1[2]+t*p2[2];
            break;
        }
    }
    interEdgeID = i;
    return type;
}
int SegmentConcavePolygonIntersection(double *p1, double *p2, double *edge, int nEdges, double *interPos){
    int interEdgeID;
    return SegmentConcavePolygonIntersection(p1, p2, edge, nEdges, interPos, interEdgeID);
}
int SegmentConcavePolygonIntersection(double *p1, double *p2, double *edge, int nEdges){
    int interEdgeID; double interPos[3];
    return SegmentConcavePolygonIntersection(p1, p2, edge, nEdges, interPos, interEdgeID);
}



