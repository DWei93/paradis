#ifndef _Math_h
#define _Math_h

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <vector>

#define SQUARE(a) ((a)*(a))
#define Dot(vec1,vec2)       \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))

#define NO_INTERSECTION         0x00001
#define OUTSIDE_SPACE           0x00002
#define INSIDE_SPACE            0x00004
#define TANGENCY                0x00008
#define POINT_INTERSECTION      0x00010
#define LINE_INTERSECTION       0x00020
#define PLANE_INTERSECTION      0x00040
#define SEGMENT_INTERSECTION    0x00080

#define LENGTH_ERROR 1.0E-11
#define ONE_DEGREE 0.0174532925199432958
#define ANGLE_ERROR1 0.0003045864904521350  // Sin[ONE_DEGREE]^2
#define ANGLE_ERROR2 0.9996954135095479 // Cos[ONE_DEGREE]^2

#ifdef GSL
void GSLTest();
#endif

extern void Fatal(const char *format, ...);
extern double  Normal(double a[3]);
extern double  Normal(std::vector<double>& a);
extern void   Normalize(double *ax, double *ay, double *az);
extern void   NormalizeVec(double vec[3]);
extern void   cross(double a[3], double b[3], double c[3]);

/** Check if the variable is numerical. */
inline bool Numerical(double v){return (!isnan(v) && !isinf(v));}
inline bool Numerical(double *vec){
    return (!isnan(vec[0]) && !isinf(vec[0]) &&
            !isnan(vec[1]) && !isinf(vec[1]) &&
            !isnan(vec[2]) && !isinf(vec[2]));
}

/** Determine the vector if is nomalized */
inline bool IsNormalized(double *vec, double error=ANGLE_ERROR2);
inline bool IsNormalized(double x, double y, double z, double error=ANGLE_ERROR2);
inline bool IsNormalized(double *vec, double error){
    double dot=Dot(vec, vec);
    return (dot>error && dot <2-error);
}
inline bool IsNormalized(double x, double y, double z, double error){
    double vec[3]={x,y,z};
    return IsNormalized(vec,error);
}

/** Determine if two directions are in the same direction. */
inline bool SameDirection(double *vec1, double *vec2, double error=ANGLE_ERROR2);
inline bool SameDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR2);
inline bool SameDirection(double *vec1, double *vec2, double error){
    double dot=Dot(vec1, vec2);
    return (dot>error && dot <2-error);
}
inline bool SameDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double vec1[3]={x1,y1,z1}, vec2[3]={x2,y2,z2};
    return SameDirection(vec1, vec2, error);
}

/** Determine if the two directions are aligned. */
inline bool AlignedDirection(double *vec1, double *vec2, double error=ANGLE_ERROR2);
inline bool AlignedDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR2);
inline bool AlignedDirection(double *vec1, double *vec2, double error){
    double dot=fabs(Dot(vec1, vec2));
    return (dot>error && dot <2-error);
}
inline bool AlignedDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double vec1[3]={x1,y1,z1}, vec2[3]={x2,y2,z2};
    return AlignedDirection(vec1, vec2, error);
}

/** Determine whether two directions are perpendicular to each other. */
inline bool PerpendicularDirection(double *vec1, double *vec2, double error=ANGLE_ERROR1);
inline bool PerpendicularDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR1);
inline bool PerpendicularDirection(double *vec1, double *vec2, double error){
    double dot=fabs(Dot(vec1, vec2));
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
 *  \return intersection type: PLANE_INTERSECTION, NO_INTERSECTION, LINE_INTERSECTION, OUTSIDE_SPACE and INSIDE_SPACE */
int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2, double *vec1, 
                           double *vec2, double err1=ANGLE_ERROR1, double err2=LENGTH_ERROR);

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
int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, double *vec1, double *vec2, 
                          double err1=ANGLE_ERROR1, double err2=LENGTH_ERROR);

/** Find intersection of point and line, the line is reference
 *  \param[in] p1[3] coordinates of the point
 *  \param[in] d[3] direction of the line
 *  \param[in] p2[3] a point on the line
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \param[out] vec[3] coordinates that stores intersection info (if has)
 *  \param[out] *dis distance of point 1 and projection point
 *  \return intersection type: POINT_INTERSECTION and NO_INTERSECTION */
int PointLineIntersection(double *p1, double *d, double *p2, double *vec, double *dis,
                          double err1=ANGLE_ERROR1, double err2=LENGTH_ERROR);

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
int LineLineIntersection(double *d1, double *p1, double *d2, double *p2, double *vec1, double *vec2, 
                         double *dis, double err1=ANGLE_ERROR1, double err2=LENGTH_ERROR);

/** Find intersection of point and plane, the plane is reference
 *  \param[in] p1[3] coordinates of the point
 *  \param[in] n[3] nromal of the line
 *  \param[in] p2[3] a point on the plane
 *  \param[out] vec[3] coordinates that stores intersection info (if has)
 *  \param[out] *dis distance of point 1 and projection point
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \return intersection type: POINT_INTERSECTION, NO_INTERSECTION, INSIDE_SPACE and OUTSIDE_SPACE */
int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis,
                            double err1=ANGLE_ERROR1, double err2=LENGTH_ERROR);

/** Find intersection of two points, the second point is reference
 *  \param[in] p1[3] coordinates of the point
 *  \param[in] p2[3] a point on the plane
 *  \param[out] vec[3] coordinates that stores intersection info (if has)
 *  \param[out] *dis distance of point 1 and projection point
 *  \param[in] err1  error for angle
 *  \param[in] err2  error for length
 *  \return intersection type: POINT_INTERSECTION and NO_INTERSECTION */
int PointPointIntersection(double *p1, double *p2, double *vec, double *dis,
                            double err1=ANGLE_ERROR1, double err2=LENGTH_ERROR);
                           

/** Return the distance to plane 
 *  \param[in] p[3] coordinates of point
 *  \param[in] plane[6] normal and a point of the plane
 *  \param[out] dis distance form the plane to the point */
inline void PointPlaneDistance(double *p, double *plane, double &dis){
    double nx=plane[0], ny=plane[1], nz=plane[2];
#ifdef DEBUG
    if(!IsNormalized(nx,ny,nz))
        Fatal("in file %s at line %d: normal (%g,%g,%g) is not normalized", 
               __FILE__,__LINE__,nx,ny,nz);
#endif
    dis = (nx*(p[0]-plane[3])+ ny*(p[1]-plane[4]) + nz*(p[2]-plane[5]))/sqrt(nx*nx+ny*ny+nz*nz);
}
/** Return the distance to plane 
 *  \param[in] p[3] coordinates of point
 *  \param[in] normal[3] normal of the plane
 *  \param[in] point[3] point on the plane
 *  \param[out] dis distance form the plane to the point */
inline void PointPlaneDistance(double *p, double *normal, double *point, double &dis){
    double nx=normal[0], ny=normal[1], nz=normal[2];
#ifdef DEBUG
    if(!IsNormalized(nx,ny,nz))
        Fatal("in file %s at line %d: normal (%g,%g,%g) is not normalized", 
               __FILE__,__LINE__,nx,ny,nz);
#endif
    dis = (nx*(p[0]-point[0])+ ny*(p[1]-point[1]) + nz*(p[2]-point[2]))/sqrt(nx*nx+ny*ny+nz*nz);
}
/** Return the distance to plane and projection
 *  \param[in] p[3] coordinates of point
 *  \param[in] plane[6] normal and a point of the plane
 *  \param[out] dis distance form the plane to the point 
 *  \param[out] projection[3] projection point on the plane*/
inline void PointPlaneDistance(double *p, double *plane, double &dis, double *projection){
    double nx=plane[0], ny=plane[1], nz=plane[2];
#ifdef DEBUG
    if(!IsNormalized(nx,ny,nz))
        Fatal("in file %s at line %d: normal (%g,%g,%g) is not normalized", 
               __FILE__,__LINE__,nx,ny,nz);
#endif
    dis = (nx*(p[0]-plane[3])+ ny*(p[1]-plane[4]) + nz*(p[2]-plane[5]))/
          sqrt(nx*nx+ny*ny+nz*nz);
    projection[0] = p[0]-dis*nx;projection[1] = p[1]-dis*ny;projection[2] = p[2]-dis*nz;
}
/** Return the distance to plane and projection on plane
 *  \param[in] p[3] coordinates of point
 *  \param[in] normal[3] normal of the plane
 *  \param[in] point[3] point on the plane
 *  \param[out] projection[3] projection on the plane
 *  \param[out] dis distance form the plane to the point */
inline void PointPlaneDistance(double *p, double *normal, double *point,  double &dis, double *projection){
    double nx=normal[0], ny=normal[1], nz=normal[2];
#ifdef DEBUG
    if(!IsNormalized(nx,ny,nz))
        Fatal("in file %s at line %d: normal (%g,%g,%g) is not normalized", 
               __FILE__,__LINE__,nx,ny,nz);
#endif
    dis = (nx*(p[0]-point[0])+ ny*(p[1]-point[1]) + nz*(p[2]-point[2]))/sqrt(nx*nx+ny*ny+nz*nz);
    projection[0] = p[0]-dis*nx;projection[1] = p[1]-dis*ny;projection[2] = p[2]-dis*nz;
}

/** Find a intersection of segment and plane (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] nDir[3] normal of plane
 *  \param[in] point[3] coordinate of point on plane
 *  \param[in] plane[6] normal and a point of the plane
 *  \param[out] t normalized position of intersection point
 *  \param[out] intersection[3] position of intersection point
 *  \return intersection type: NO_INTERSECTION|INSIDE_SPACE, NO_INTERSECTION|OUTSIDE_SPACE, and POINT_INTERSECTION */ 
inline int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double &t){

    double dot1 = nDir[0]*(p1[0]-point[0]) + nDir[1]*(p1[1]-point[1]) + nDir[2]*(p1[2]-point[2]);
    double dot2 = nDir[0]*(p2[0]-point[0]) + nDir[1]*(p2[1]-point[1]) + nDir[2]*(p2[2]-point[2]);
#ifdef DEBUG
    if(!IsNormalized(nDir))Fatal("%s: normal (%g,%g,%g) is not normalized",__func__, nDir[0], nDir[1],nDir[2]);
    if(!Numerical(dot1))Fatal("%s: dot1=%f",__func__,dot1);
    if(!Numerical(dot2))Fatal("%s: dot2=%f",__func__,dot2);
#endif
    if(dot1 > 0 && dot2>0)return(NO_INTERSECTION|INSIDE_SPACE);
    if(dot1 < 0 && dot2<0)return(NO_INTERSECTION|OUTSIDE_SPACE);
    double tot = fabs(dot1)+fabs(dot2);
    if(tot > 1.0E-80) t = fabs(dot1)/tot;
    else t = 0.0;
#ifdef DEBUG
    if(!Numerical(t))Fatal("%s: t=%f", __func__, t);
#endif
    return(POINT_INTERSECTION);
}
inline int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double *intersection){
    int type; double t;
    if((type=SegmentPlaneIntersection(p1,p2,nDir,point,t))==POINT_INTERSECTION){
        intersection[0] = (1.0-t)*p1[0]+t*p2[0];
        intersection[1] = (1.0-t)*p1[1]+t*p2[1];
        intersection[2] = (1.0-t)*p1[2]+t*p2[2];
    }
    return type;
}
inline int SegmentPlaneIntersection(double *p1, double *p2, double *plane, double *intersection){
    return SegmentPlaneIntersection(p1,p2,plane,plane+3,intersection);
}

/** Find a intersection of segment and plane (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] nDir[3] normal of plane
 *  \param[in] point[3] coordinate of point on plane
 *  \param[in] plane[6] normal and a point of the plane
 *  \param[in] margin margin of the plane
 *  \param[out] t normalized position of intersection point
 *  \param[out] intersection[3] position of intersection point
 *  \return intersection type: NO_INTERSECTION|INSIDE_SPACE, NO_INTERSECTION|OUTSIDE_SPACE, and POINT_INTERSECTION */ 
inline int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double margin, double &t){
    double mPoint[3]={point[0]+margin*nDir[0],point[1]+margin*nDir[1],point[2]+margin*nDir[2]};
    return SegmentPlaneIntersection(p1, p2, nDir, mPoint, t);
}
inline int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double margin, double *intersection){
    double mPoint[3]={point[0]+margin*nDir[0],point[1]+margin*nDir[1],point[2]+margin*nDir[2]};
    return SegmentPlaneIntersection(p1, p2, nDir, mPoint, intersection);
}
inline int SegmentPlaneIntersection(double *p1, double *p2, double *plane, double margin, double *intersection){
    double mPoint[3]={plane[3]+margin*plane[0],plane[4]+margin*plane[1],plane[5]+margin*plane[2]};
    return SegmentPlaneIntersection(p1, p2, plane, mPoint, intersection);
}


/** Determines whether the node is inside a convex polyhedron
 *  \param[in] p[3] coordinates of the point
 *  \param[in] bounds boundary info. contains all plane normals and points on each plane
 *  \param[in] nBounds # of boundaries
 *  \param[in] margin margin of the boundaries
 *  \return true or false */
inline bool InConcavePolyhedron(double *p, double **bounds, int nBounds, double margin)
{
    double dis;
    for(auto i=0; i<nBounds; i++){
        PointPlaneDistance(p, bounds[i], dis);
        if(dis>margin)return false;
    }
    return true;
}
inline bool InConcavePolyhedron(double *p, double *bounds, int nBounds, double margin)
{
    double dis;
    for(auto i=0; i<nBounds; i++){
        PointPlaneDistance(p, &bounds[i*6], dis);
        if(dis>margin)return false;
    }
    return true;
}

/** Clip a segment in a concave polyhedron
 *  \param[in] p1[3] coordinates of point 1, this point is default in the polyhedron
 *  \param[in,out] p2[3] coordinates of point 2, and return the intersection point (if has)
 *  \param[in] bounds boundaries info. contains nromals and points of the polyhedron
 *  \param[in] nBounds # of boundaries
 *  \param[in] margin margin of the polyhedron
 *  \param[out] i ith boundary that has the intersection 
 *  \return NO_INTERSECTION|OUTSIDE_SPACE, POINT_INTERSECTION */
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, int &i)
{
    int type;
    for(i=0;i<nBounds;i++)if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,p2))==POINT_INTERSECTION)break;
    return type;
}
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, double margin, int &i)
{
    int type;
    for(i=0;i<nBounds;i++)if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,margin,p2))==POINT_INTERSECTION)break;
    return type;
}

/** Clip a segment in a concave polyhedron
 *  \param[in] p1[3] coordinates of point 1, this point is default in the polyhedron
 *  \param[in] p2[3] coordinates of point 2
 *  \param[in] bounds boundaries info. contains nromals and points of the polyhedron
 *  \param[in] nBounds # of boundaries
 *  \param[in] margin margin of the polyhedron
 *  \param[out] i ith boundary that has the intersection 
 *  \param[out] clipPos[3], clipping position
 *  \return NO_INTERSECTION|OUTSIDE_SPACE, POINT_INTERSECTION */
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, int &i, double *clipPos)
{
    int type;
    for(i=0;i<nBounds;i++)if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,clipPos))==POINT_INTERSECTION)break;
    return type;
}
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, double margin, int &i, double *clipPos)
{
    int type;
    for(i=0;i<nBounds;i++)if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,margin,clipPos))==POINT_INTERSECTION)break;
    return type;
}

#endif
