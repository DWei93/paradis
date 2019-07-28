#ifndef _Math_h
#define _Math_h

#define SQUARE(a) ((a)*(a))

#define NO_INTERSECTION     0
#define POINT_INTERSECTION  1
#define LINE_INTERSECTION   2
#define PLANE_INTERSECTION  3

#define LENGTH_ERROR 1.0E-11
#define ONE_DEGREE 0.0174532925199432958
#define ANGLE_ERROR1 0.0003045864904521350
#define ANGLE_ERROR2 0.9996954135095479


#ifdef GSL
void GSLTest();
#endif


/** Check if the variable is numerical. */
inline bool Numerical(double v);
inline bool Numerical(double *vec);

/** Determine the vector if is nomalized */
inline bool IsNormalized(double *vec, double error=ANGLE_ERROR2);
inline bool IsNormalized(double x, double y, double z, double error=ANGLE_ERROR2);

/** Determine if two directions are in the same direction. */
inline bool SameDirection(double *vec1, double *vec2, double error=ANGLE_ERROR2);
inline bool SameDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR2);

/** Determine if the two directions are aligned. */
inline bool AlignedDirection(double *vec1, double *vec2, double error=ANGLE_ERROR2);
inline bool AlignedDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR2);

/** Determine whether two directions are perpendicular to each other. */
inline bool PerpendicularDirection(double *vec1, double *vec2, double error=ANGLE_ERROR1);
inline bool PerpendicularDirection(double x1, double y1, double z1, double x2, double y2, double z2, double error=ANGLE_ERROR1);

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
 *  \return intersection type: POINT_INTERSECTION, NO_INTERSECTION and LINE_INTERSECTION */
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
 *  \return intersection type: POINT_INTERSECTION and NO_INTERSECTION */
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
                           
/** Find a intersection of segment and plane (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] nDir[3] normal of plane
 *  \param[int] point[3] coordinate of point on pkane
 *  \param[out] t normalized position of intersection point */ 
int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double &t);

/** Find a intersection of segment and concave polygon (if has) 
 *  \param[in] p1[3] coordinates of the first end of segment
 *  \param[in] p2[3] coordinates of the second  end of segment
 *  \param[in] edges[6*nEdges] planes that containning concave polygon edges,
 *             each plane has normals and coodinates of point on the plane
 *  \param[in] nEdges # of edges of the concave polygon
 *  \param[out] interPos[3] the position of intersection point
 *  \param[out] interEdgeID index of edge which intersects with the segment
 *  \return intersection type : NO_INTERSECTION and POINT_INTERSECTION */ 
int SegmentConcavePolygonIntersection(double *p1, double *p2, double *edge, int nEdges, double *interPos, int &interEdgeID);
int SegmentConcavePolygonIntersection(double *p1, double *p2, double *edge, int nEdges, double *interPos);
int SegmentConcavePolygonIntersection(double *p1, double *p2, double *edge, int nEdges);

#endif
