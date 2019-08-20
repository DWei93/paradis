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


#ifdef GSL
void GSLTest()
{
  double x = 1212.0;
  double y = gsl_sf_bessel_J0 (x);
  printf ("J0(%g) = %.18e\n", x, y);
  return ;
}
#endif


/** Find intersection of two planes, plane2 is the default reference plane.
 *  \param[in] n1[3] normlized normal of plane1
 *  \param[in] p1[3] a point on plane1
 *  \param[in] n2[3] normlized normal of plane2
 *  \param[in] p2[3] a point on plane2
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \param[out] vec1[3] normalized vector that stores intersection info (if has)
 *  \param[out] vec2[3] second vector that stores intersection info (if has)
 *  \return intersection type: PLANE_INTERSECTION, NO_INTERSECTION, LINE_INTERSECTION, OUTSIDE_SPACE and INSIDE_SPACE */
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
        double dot=n2x*(p1x - p2x) + n2y*(p1y - p2y) + n2z*(p1z - p2z);
        if(fabs(dot) < err2) return(PLANE_INTERSECTION);
        else if(dot > 0) return (NO_INTERSECTION|INSIDE_SPACE);
        else return (NO_INTERSECTION|OUTSIDE_SPACE);
    } 

    // the direction of intersection line is the cross product vector of the two normals
    vec1[0] = -(n1z*n2y) + n1y*n2z;
    vec1[1] = n1z*n2x - n1x*n2z;
    vec1[2] = -(n1y*n2x) + n1x*n2y;
    normalize_vec(vec1);

    // if the two points are too closed, point2 is chosed to be the point on intersection line
    if(SQUARE(p1x - p2x) + SQUARE(p1y - p2y) + SQUARE(p1z - p2z) < err2*err2){
        VECTOR_COPY(vec2,p2); return(LINE_INTERSECTION);
    }

    // the glide direction from point1 to plane2 is defined by the cross production of normal1 and intersection line direction
    cross_product(vec1, n1, vec3);
#ifdef DEBUG 
    if(vec3[0]*vec3[0]+vec3[1]*vec3[1]+vec3[2]*vec3[2]< 1.0E-8)Fatal("%s: cross porducet vector is zero\n", __func__);
#endif
    normalize_vec(vec3);

    // the key point of inetrsection line can be defined by the intersection point of glide line and plane2
    type = LinePlaneIntersection(vec3, p1, n2, p2, vec4, vec2, err1, err2);
#ifdef DEBUG
    if((type&NO_INTERSECTION)>0)Fatal("%s: can not find the intersection of line and plane", __func__);
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
 *  \return intersection type: POINT_INTERSECTION, NO_INTERSECTION, LINE_INTERSECTION, OUTSIDE_SPACE, and INSIDE_SPACE*/
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
        int type = PointPlaneIntersection(p1, n, p2, vec2, &dis);
        VECTOR_COPY(vec1, d); VECTOR_COPY(vec2, p1);
        if(fabs(dis) < err2) return(LINE_INTERSECTION);
        else return(type);
    }

    double t = (-(nx*p1x) - ny*p1y - nz*p1z + nx*p2x + ny*p2y + nz*p2z)/
                (dx*nx + dy*ny + dz*nz);

#ifdef DEBUG
    if(!Numerical(t))Fatal("%s: %f",__func__,t);
#endif
    vec2[0] = p1x + t*dx;    
    vec2[1] = p1y + t*dy;    
    vec2[2] = p1z + t*dz;

#ifdef DEBUG
    if(!Numerical(vec2)){
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
    if(!Numerical(denominator))Fatal("%s: the denominator is %f",__func__,denominator);
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
    if(!Numerical(denominator))Fatal("%s: the denominator is %f",__func__,denominator);
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
    if(!Numerical(*dis))Fatal("%s: *dis=%f",__func__, *dis);
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
 *  \return intersection type: POINT_INTERSECTION, NO_INTERSECTION, INSIDE_SPACE and OUTSIDE_SPACE */
int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis, double err1, double err2)
{
#if 1
    double distance;
    PointPlaneDistance(p1, n, p2, distance, vec);
    *dis = distance;
#else
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
    normalize_vec(temVec);
    if(DotProduct(temVec, n) < 0)*dis *= -1.0; 
#endif 
    if(fabs(*dis) < err2) return(POINT_INTERSECTION);
    else if(*dis>0)return(NO_INTERSECTION|INSIDE_SPACE);
    else return (NO_INTERSECTION|OUTSIDE_SPACE);
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

/** Get the distance between point and segment
 *  \param[in] p the coordinate of the point
 *  \param[in] p1,p2 the two end of the segment
 *  \return distance       */
double PointSegmentDistance(double *p, double *p1, double *p2){
    double lDir[3], vec2[3], dis, len;
    CoordMinus(p2,p1,lDir);
    len=normal(lDir);
#ifdef DEBUG
    if(len<LENGTH_ERROR)Fatal("%s: the segments is too short",__func__);
#endif
    normalize_vec(lDir);
    PointLineIntersection(p, lDir, p1, vec2, &dis);
    vec2[0]-=p1[0]; vec2[1]-=p1[1]; vec2[2]-=p1[2]; 
    
    if(Dot(vec2,lDir)<0)return PointPointDistance(p,p1);
    else if (normal(vec2)>len) return PointPointDistance(p,p2);
    else return dis;
}

/** Get the distance between point and polygon
 *  \param[in] p the coordinate of the point
 *  \param[in] v coordinates lis of vertices
 *  \param[in] polygon indexs of polygon vertices
 *  \param[in] startID index of start
 *  \return distance                    */
double PointConcavePolygonDistance(double *p, std::vector<double> &v, std::vector<int> &polygon, int startID){
    double p1[3], p2[3], dis, min=1E50;
    for(auto i=startID+1; i<startID+1+polygon[startID]; i++){
        Coordinate(v,polygon[i],p1);
        if(i<startID+polygon[startID])Coordinate(v,polygon[i+1],p2);
        else Coordinate(v,polygon[startID+1],p2);
        
        dis = PointSegmentDistance(p, p1, p2);
        min = (min>dis)?dis:min;
    }
    return min;
}

/** Find the intersection segment between plane and concave polygon`
 *  \param[in] n,p normal and point of the plane
 *  \param[in] vertices coodinates list
 *  \param[in] polygon index of polygon vertices
 *  \param[in] nPoints number polygon vertices
 *  \param[out]  intersections two coodinates of te intersection segment (if has)*/
bool PlaneConcavePolygonIntersection(double *n, double *p, double *vertices, int *polygon, int nPoints, double intersections[6])
{
    int next, j=0;
    double *v1, *v2, pos[3];
    for(auto i=0; i<nPoints; i++){
        next = (i<nPoints-1)?i+1:0;
        if(POINT_INTERSECTION==SegmentPlaneIntersection(Coordinate(vertices, polygon[i]), Coordinate(vertices,polygon[next]),n, p, pos)){
            intersections[3*j]=pos[0];intersections[3*j+1]=pos[1];intersections[3*j+2]=pos[2]; j++;
        }
        if(j==2)break;
    }
    if(j==2)return true;
    else return false;
}

/** Caculate centroid of the polygon
 *  \param[in] f_vert vertice indexs of polygons
 *  \param[in] v coodinates of all vertices
 *  \param[in] startID the index of polygon to be caculeted
 *  \param[out] cnetroid the centroid */
void CentroidOfPolygon(std::vector<int> &f_vert,std::vector<double> &v,int startID, double *centroid){

	int n=f_vert[startID++];

    double pos0[3], pos1[3], pos2[3], cpd[3], areas=0, area;
    VECTOR_ZERO(centroid); 
    for(int i=startID+1; i<startID+n-1; i++){
        pos1[0]=v[f_vert[i]*3]-v[f_vert[startID]*3];
        pos1[1]=v[f_vert[i]*3+1]-v[f_vert[startID]*3+1];
        pos1[2]=v[f_vert[i]*3+2]-v[f_vert[startID]*3+2];
        pos2[0]=v[f_vert[i+1]*3]-v[f_vert[startID]*3];
        pos2[1]=v[f_vert[i+1]*3+1]-v[f_vert[startID]*3+1];
        pos2[2]=v[f_vert[i+1]*3+2]-v[f_vert[startID]*3+2];

        pos0[0] = (pos1[0]+pos2[0])/3.0;
        pos0[1] = (pos1[1]+pos2[1])/3.0;
        pos0[2] = (pos1[2]+pos2[2])/3.0;
        cross_product(pos1, pos2, cpd);area = normal(cpd);
        areas += area;
        centroid[0] += area*pos0[0]; 
        centroid[1] += area*pos0[1]; 
        centroid[2] += area*pos0[2]; 
    }
    centroid[0] /= areas; centroid[1] /= areas; centroid[2] /= areas;
    centroid[0] += v[f_vert[startID]*3];
    centroid[1] += v[f_vert[startID]*3+1];
    centroid[2] += v[f_vert[startID]*3+2];
    return;
}

/** Caculate area of the polygon
 *  \param[in] f_vert vertice indexs of polygons
 *  \param[in] v coodinates of all vertices
 *  \param[in] startID the index of polygon to be caculeted
 *  \return area of polygon */
double AreaOfPolygon(std::vector<int> &f_vert,std::vector<double> &v,int startID){
	int n=f_vert[startID++];

    double pos0[3], pos1[3], pos2[3], cpd[3], areas=0, area;
    for(int i=startID+1; i<startID+n-1; i++){
        pos1[0]=v[f_vert[i]*3]-v[f_vert[startID]*3];
        pos1[1]=v[f_vert[i]*3+1]-v[f_vert[startID]*3+1];
        pos1[2]=v[f_vert[i]*3+2]-v[f_vert[startID]*3+2];
        pos2[0]=v[f_vert[i+1]*3]-v[f_vert[startID]*3];
        pos2[1]=v[f_vert[i+1]*3+1]-v[f_vert[startID]*3+1];
        pos2[2]=v[f_vert[i+1]*3+2]-v[f_vert[startID]*3+2];

        cross_product(pos1, pos2, cpd);area = normal(cpd);
        areas += area;
    }
    return 0.5*areas;
}

/** add and search vertices, used to assemble vertices list
 *  \param[in,out] vertices the coordinates of vertice list
 *  \param[in] x,y,z the coordinate to search or add
 *  \return the index of the coordinate */
int AddAndSearchVertice(std::vector<double> &vertices, double x, double y, double z)
{
    if(vertices.size() > 0){
        unsigned int i;
        for(i=0; i<vertices.size(); i+=3){
            if(fabs(vertices[i]-x)<1.0E-12 &&
               fabs(vertices[i+1]-y)<1.0E-12 &&
               fabs(vertices[i+2]-z)<1.0E-12)break;
        }
        if(i == vertices.size()){
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
        return int(i/3);
    }else{
        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);
        return(0);
    }
}

