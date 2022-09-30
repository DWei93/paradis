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
#define ANGLE_ERROR1 0.00030458649045213  // Sin[ONE_DEGREE]^2
#define ANGLE_ERROR2 0.9996954135095479 // Cos[ONE_DEGREE]^2
#ifdef GSL
void GSLTest();
#endif

extern void Fatal(const char *format, ...);

/** normal of a vector */
inline double normal(double *a)
{return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);}
inline double normal(double x, double y, double z)
{return sqrt(x*x + y*y + z*z);}

/** Normalize vector 
 *  \param[in,out] vec vector to be normalized */
inline void normalize_vec(double *vec){
    double a2, a;
    a2 = (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (a2 > 0.0) {a = sqrt(a2);vec[0] /= a;vec[1] /= a;vec[2] /= a;}
    return;
}

/** Cross product of vector
 *  \param[in] a first vector 
 *  \param[in] b second vector
 *  \param[out] c coross product vector */
inline void cross_product(double *a, double *b, double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
inline void cross_product(double a0, double a1, double a2, double b0, double b1, double b2, double *c)
{
    c[0]=a1*b2-a2*b1;
    c[1]=a2*b0-a0*b2;
    c[2]=a0*b1-a1*b0;
}
#define CrossProduct(a,b,c) \
{\
    c[0]=a[1]*b[2]-a[2]*b[1]; \
    c[1]=a[2]*b[0]-a[0]*b[2]; \
    c[2]=a[0]*b[1]-a[1]*b[0];\
}

/* Dot production */
inline double dot_product(double a0, double a1, double a2, double b0, double b1, double b2){
    return a0*b0+a1*b1+a2*b2;
}

/* Negative vector*/
inline void NegativeVector(double *a, double *b){
    b[0] = -a[0];
    b[1] = -a[1];
    b[2] = -a[2];
}

/** Parallel component of vector a along the direction dir
 *  \param[in] a vector
 *  \param[in] dir direction (should be normalized)
 *  \param[out] vec parallel component  */
inline void ParallelComponent(double *a, double *dir, double *vec){
    double dot = Dot(a,dir); 
    vec[0] = dot*dir[0]; vec[1]=dot*dir[1]; vec[2]=dot*dir[2];
}

/** Normal component of vector a along the direction dir
 *  \param[in] a vector
 *  \param[in] dir direction (should be normalized)
 *  \param[out] vec normal component  */
inline void NormalComponent(double *a, double *dir, double *vec){
    double dot = Dot(a,dir); 
    vec[0] = a[0]-dot*dir[0]; vec[1]=a[1]-dot*dir[1]; vec[2]=a[2]-dot*dir[2];
}

/* Copy Vectors*/                
// Copying the values of the source vector into the destination vector.
template<class v>
inline void VectorCopy(v *dest, v *source){dest[0] = source[0];dest[1] = source[1];dest[2] = source[2];}
// Copying the values of s1, s2, and s3 into the vector dest.
template<class v>
inline void VectorCopy(v *dest, v s1, v s2, v s3){dest[0] = s1;dest[1] = s2;dest[2] = s3;}
// Copying the contents of source into dest.
template<class v>
inline void LongVectorCopy(v *dest, v *source, int num){memcpy(dest, source, num*sizeof(v));}
// Copying the first three elements of the source vector into the three destination vectors.
template<class v>
inline void VectorCopy(v &dest1, v &dest2, v &dest3, v *source){dest1 = source[0];dest2 = source[1];dest3 = source[2];}

/* Restore Vectors
 *  \param[in] dest bakup vector
 *  \param[out] s1 first component of the restored vector
 *  \param[out] s2 second component of the restored vector
 *  \param[out] s3 third component of the restored vector
 */
template<class v>
inline void VectorRestore(v *dest, v &s1, v &s2, v &s3){s1=dest[0];s2=dest[1];s3=dest[2];}


/* Two opearations of Coodinate*/
template<class v> 
inline void CoordAdd(v *a, v *b, v c[3]){c[0]=a[0]+b[0]; c[1]=a[1]+b[1]; c[2]=a[2]+b[2];}
template<class v, class q> 
inline void CoordAdd(q f1, v *a, q f2, v *b, v c[3]){
c[0]=f1*a[0]+f2*b[0]; c[1]=f1*a[1]+f2*b[1]; c[2]=f1*a[2]+f2*b[2];}
template<class v> 
inline void CoordMinus(v *a, v *b, v c[3]){c[0]=a[0]-b[0]; c[1]=a[1]-b[1]; c[2]=a[2]-b[2];}
template<class v>
inline void VectorAddSelf(v *a, v *b){a[0]+=b[0]; a[1]+=b[1]; a[2]+=b[2];}
template<class v>
inline void VectorMinusSelf(v *a, v *b){a[0]-=b[0]; a[1]-=b[1]; a[2]-=b[2];}
template<class v>
inline void VectorTimesSelf(v *a, v *b){a[0]*=b[0]; a[1]*=b[1]; a[2]*=b[2];}
template<class v>
inline void VectorDivideSelf(v *a, v *b){a[0]/=b[0]; a[1]/=b[1]; a[2]/=b[2];}

template<class v>
inline void VectorAddSelf(v *a, v b){a[0]+=b; a[1]+=b; a[2]+=b;}
template<class v>
inline void VectorMinusSelf(v *a, v b){a[0]-=b; a[1]-=b; a[2]-=b;}
template<class v>
inline void VectorTimesSelf(v *a, v b){a[0]*=b; a[1]*=b; a[2]*=b;}
template<class v>
inline void VectorDivideSelf(v *a, v b){a[0]/=b; a[1]/=b; a[2]/=b;}

/*return a coordinate from a vector*/
template<class v>
inline v *Coordinate(v *vector, int element){return (vector+3*element);}
template<class v>
inline void Coordinate(v *vector, int element, v a[3])
{a[0]=vector[3*element];a[1]=vector[3*element+1];a[2]=vector[3*element+2];}
template<class v>
inline void Coordinate(std::vector<v> &vector, int element, v p[3])
{p[0] = vector[3*element]; p[1] = vector[3*element+1]; p[2] = vector[3*element+2];}

template<class v>
inline void Oppsite(v a[3], v *b){b[0]=-a[0]; b[1]=-a[1]; b[2]=-a[2]; return;}

/** Check if the variable is numerical. */
inline bool Numerical(double v){return (!isnan(v) && !isinf(v));}
inline bool Numerical(double *vec){
    return (!isnan(vec[0]) && !isinf(vec[0]) &&
            !isnan(vec[1]) && !isinf(vec[1]) &&
            !isnan(vec[2]) && !isinf(vec[2]));
}
inline bool Numerical(double ax, double ay, double az){
    return (!isnan(ax) && !isinf(ax) &&
            !isnan(ay) && !isinf(ay) &&
            !isnan(az) && !isinf(az));

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
    double dot=x1*x2+y1*y2+z1*z2;
    return (dot>error && dot <2-error);
}

/** Determine if the two directions are aligned. */
inline bool AlignedDirection(double *vec1, double *vec2, double error=ANGLE_ERROR2);
inline bool AlignedDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR2);
inline bool AlignedDirection(double *vec1, double *vec2, double error){
    double dot=fabs(Dot(vec1, vec2));
    return (dot>error && dot <2-error);
}
inline bool AlignedDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double dot=fabs(x1*x2+y1*y2+z1*z2);
    return (dot>error && dot <2-error);
}

/** Determine whether two directions are perpendicular to each other. */
inline bool PerpendicularDirection(double *vec1, double *vec2, double error=ANGLE_ERROR1);
inline bool PerpendicularDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR1);
inline bool PerpendicularDirection(double *vec1, double *vec2, double error){
    double dot=fabs(Dot(vec1, vec2));
    return dot<error;
}
inline bool PerpendicularDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error){
    double dot=fabs(x1*x2+y1*y2+z1*z2);
    return dot<error;
}

/** Dertimne whether two vectors are aligned with each other and have same nagnitude */
inline bool AlignedVector(double *v1, double *v2, double error1=ANGLE_ERROR2, double error2=ANGLE_ERROR1){
    double dot=fabs(Dot(v1, v2)); double s1=Dot(v1,v1), s2=Dot(v2,v2);
    return (dot>error1 && dot <2-error1 && fabs(s1-s2) < error2);
}
inline bool AlignedVector(double x1, double y1, double z1, double x2, double y2, double z2, double error1=ANGLE_ERROR2, double error2=ANGLE_ERROR1){
    double dot=fabs(x1*x2+y1*y2+z1*z2); double s1=x1*x1+y1*y1+z1*z1; double s2=x2*x2+y2*y2+z2*z2;
    return (dot>error1 && dot <2-error1 && fabs(s1-s2) < error2);
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

/** Return the distance to plane 
 *  \param[in] p[3] coordinates of point
 *  \param[in] plane[6] normal and a point of the plane
 *  \param[in] point[3] point on the plane
 *  \param[out] projection[3] projection on the plane
 *  \param[out] dis distance form the plane to the point */
void PointPlaneDistance(double *p, double *plane, double &dis);
double PointPlaneDistance(double *p, double *normal, double *point);
void PointPlaneDistance(double *p, double *normal, double *point, double &dis);
void PointPlaneDistance(double *p, double *plane, double &dis, double *projection);
void PointPlaneDistance(double *p, double *normal, double *point,  double &dis, double *projection);

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
                           
/** Return the distance between two ponts
 *  \param[in] p1,p2 two node
 *  \return  distance */
inline double PointPointDistance(double *p1, double *p2)
{return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+ (p1[2]-p2[2])*(p1[2]-p2[2]));} 
inline double PointPointDistanceSquare(double *p1, double *p2)
{return (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+ (p1[2]-p2[2])*(p1[2]-p2[2]);}
inline double PointPointDistance(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z)
{return sqrt((p1x-p2x)*(p1x-p2x)+(p1y-p2y)*(p1y-p2y)+ (p1z-p2z)*(p1z-p2z)); }

/** Find a intersection of segment and plane (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] nDir[3] normal of plane
 *  \param[in] point[3] coordinate of point on plane
 *  \param[in] plane[6] normal and a point of the plane
 *  \param[out] t normalized position of intersection point
 *  \param[out] intersection[3] position of intersection point
 *  \param[out] dot1, dot2, signed distances between end and plane
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
    else {
#ifdef DEBUG
        printf("Warning in %s: dot1 %e, dot2 %e\n",__func__, dot1, dot2);
#endif
        t = 0.0;
    }
#ifdef DEBUG
    if(!Numerical(t))Fatal("%s: t=%f", __func__, t);
#endif
    return(POINT_INTERSECTION);
}
inline int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double *intersection){
    int type; double t;
    if((type=SegmentPlaneIntersection(p1,p2,nDir,point,t))==POINT_INTERSECTION){
        CoordAdd(1.0-t, p1, t, p2, intersection);
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
 *  \param[in] margin outer margin of the boundaries
 *  \return true or false */
inline bool InConcavePolyhedron(double *p, double **bounds, int nBounds, double margin)
{
    double dis;
    for(int i=0; i<nBounds; i++){
        PointPlaneDistance(p, bounds[i], dis);
        if(dis>margin)return false;
    }
    return true;
}
inline bool InConcavePolyhedron(double *p, double *bounds, int nBounds, double margin)
{
    double dis;
    for(int i=0; i<nBounds; i++){
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
 *  \param[in] margin outer margin of the polyhedron
 *  \param[out] edgeID the edge index that has the intersection 
 *  \return NO_INTERSECTION|OUTSIDE_SPACE, POINT_INTERSECTION */
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, int &edgeID)
{
    int i, type=-1, ftype=1, outspace=0;
    for(i=0;i<nBounds;i++){
        if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,p2))==POINT_INTERSECTION){edgeID=i; ftype=type;}
        else if((type&OUTSIDE_SPACE)>0)outspace=1;
    }
    if((ftype&NO_INTERSECTION)>0 && outspace==1)ftype|=(NO_INTERSECTION|OUTSIDE_SPACE);
    return ftype;
}
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, double margin, int &edgeID)
{
    int i, type=-1, ftype=-1, outspace=0;
    for(i=0;i<nBounds;i++){
        if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,margin,p2))==POINT_INTERSECTION){edgeID=i; ftype=type;}
        else if((type&OUTSIDE_SPACE)>0)outspace=1;
    }
    if((ftype&NO_INTERSECTION)>0 && outspace==1)ftype|=(NO_INTERSECTION|OUTSIDE_SPACE);
    return ftype;
}

/** Clip a segment in a concave polyhedron
 *  \param[in] p1[3] coordinates of point 1, this point is default in the polyhedron
 *  \param[in] p2[3] coordinates of point 2
 *  \param[in] bounds boundaries info. contains nromals and points of the polyhedron
 *  \param[in] nBounds # of boundaries
 *  \param[in] margin margin of the polyhedron
 *  \param[out] edgeID the edge index that has the intersection 
 *  \param[out] clipPos[3], clipping position
 *  \return NO_INTERSECTION|OUTSIDE_SPACE, POINT_INTERSECTION */
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, int &edgeID, double *clipPos)
{
    int i, type=-1, ftype=-1;
    for(i=0;i<nBounds;i++)if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,clipPos))==POINT_INTERSECTION){edgeID=i;ftype=type;}
    return ftype;
}
inline int ClipSegmentInConcavePolyhedron(double *p1, double *p2, double *bounds, int nBounds, double margin, int &edgeID, double *clipPos)
{
    int i, type=-1,ftype=-1;
    for(i=0;i<nBounds;i++)if((type=SegmentPlaneIntersection(p1,p2,bounds+6*i,margin,clipPos))==POINT_INTERSECTION){edgeID=i;ftype=type;}
    return ftype;
}

/** Get the distance between point and segment
 *  \param[in] p the coordinate of the point
 *  \param[in] p1,p2 the two end of the segment
 *  \return distance       */
double PointSegmentDistance(double *p, double *p1, double *p2);

/** Get the distance between point and polygon
 *  \param[in] p the coordinate of the point
 *  \param[in] v coordinates lis of vertices
 *  \param[in] polygon indexs of polygon vertices
 *  \param[in] startID index of start
 *  \return distance                    */
double PointConcavePolygonDistance(double *p, std::vector<double> &v, std::vector<int> &polygon, int startID);

/** find the intersection segment between plane and concave polygon
 *  \param[in] n,p normal and point of the plane
 *  \param[in] vertices coodinates list
 *  \param[in] polygon index of polygon vertices
 *  \param[in] nPoints number polygon vertices
 *  \param[out]  intersections two coodinates of te intersection segment (if has)*/
bool PlaneConcavePolygonIntersection(double *n, double *p, double *vertices, int *polygon, int nPoints, double intersections[6]);

/** Caculate centroid of the polygon
 *  \param[in] f_vert vertice indexs of polygons
 *  \param[in] v coodinates of all vertices
 *  \param[in] startID the index of polygon to be caculeted
 *  \param[out] cnetroid the centroid */
void CentroidOfPolygon(std::vector<int> &f_vert,std::vector<double> &v,int startID, double *centroid);

/** Caculate area of the polygon
 *  \param[in] f_vert vertice indexs of polygons
 *  \param[in] v coodinates of all vertices
 *  \param[in] startID the index of polygon to be caculeted
 *  \return area of polygon */
double AreaOfPolygon(std::vector<int> &f_vert,std::vector<double> &v,int startID);

/** add and search vertices, used to assemble vertices list
 *  \param[in,out] vertices the coordinates of vertice list
 *  \param[in] x,y,z the coordinate to search or add
 *  \return the index of the coordinate */
int AddAndSearchVertice(std::vector<double> &vertices, double x, double y, double z);

inline void BubbleSort2D(double *x, double *y, int n){
    double temp;
    for(int i=n-1; i>0; i--){
        for( int j=0; j<i; j++){
            if(x[j]>x[j+1]){
                temp=x[j]; x[j]=x[j+1]; x[j+1]=temp; 
                temp=y[j]; y[j]=y[j+1]; y[j+1]=temp; 
            }
        }
    }
    return;
}

bool FourPointCoplanar(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4,double y4, double z4, double err = LENGTH_ERROR);

/** Decompose segment in a box
 *  \param[in] x1,y1,z1 coordinate of first point on the segement
 *  \param[in] x2,y2,z2 coordinate of first point on the segement
 *  \param[in] box the range of box: xmin, xmax, ymin, ymax, zmin, zmax
 */
double DecomposeSegmentInBox(double *p1, double *p2, double *box);

bool inline InBox(double *p, double *box, double err=-LENGTH_ERROR){
    return( p[0]-box[0] > err && 
            box[1]-p[0] > err && 
            p[1]-box[2] > err && 
            box[3]-p[1] > err && 
            p[2]-box[4] > err && 
            box[5]-p[2] > err);
}

#endif
