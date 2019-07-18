
#ifndef _Math_h
#define _Math_h

#define SQUARE(a) ((a)*(a))

#define NO_INTERSECTION     0
#define POINT_INTERSECTION  1
#define LINE_INTERSECTION   2
#define PLANE_INTERSECTION  3
#define EPS1 1.0E-3
#define EPS2 1.0E-1

int Arrange(int n);
int Arrange(int n, int m);
double Volume(double **p, int num);

int PlanePlaneIntersection(double *n1, double *p1, double *n2, double*p2,
                           double *vec1, double *vec2);
int LinePlaneIntersection(double *d, double *p1, double *n, double *p2, 
                          double *vec1, double *vec2);
int LineLineIntersection(double *d1, double *p1, double *d2, double *p2,
                         double *vec1, double *vec2, double *dis);
int PointPlaneIntersection(double *p1, double *n, double *p2, double *vec, double *dis);
int PointLineIntersection(double *p1, double *d, double *p2, double *vec, double *dis);
int PointPointIntersection(double *p1, double *p2, double *vec, double *dis);
int SegmentPlaneIntersection(double *p1, double *p2, double *nDir, double *point, double &t);
int SegmentPolygonIntersection2D(double *p1, double *p2, double *edge, int nEdges, double *interPos);
int SegmentPolygonIntersection2D(double *p1, double *p2, double *edge, int nEdges, double *interPos, int &interEdgeID);;
#endif
