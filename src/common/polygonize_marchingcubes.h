#ifndef _POLYGONIZE_MARCHINGCUBES_H_
#define _POLYGONIZE_MARCHINGCUBES_H_
#include <vector>
using std::vector;

#ifndef EPSVERTEXDIST2
#define EPSVERTEXDIST2 1.0e-10
#endif

/* ************************************************************************** */
/** This class is a refactor of Paul Bourke's implementation.
   It closely follows the code 'example.c' contained in volexample.zip.
   For more details, see:
   http://paulbourke.net/geometry/polygonise/
   (last access: 2020-Jul-02)
*/
class PolygonizeMarchingCubes {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   typedef struct {
      double x,y,z;
   } XYZ;
   typedef struct {
      XYZ p[8];
      XYZ n[8];
      double val[8];
   } GRIDCELL;
   typedef struct {
      XYZ p[3];         /* Vertices */
      XYZ c;            /* Centroid */
      XYZ n[3];         /* Normal   */
   } TRIANGLE;
/* ************************************************************************** */
   PolygonizeMarchingCubes();
   ~PolygonizeMarchingCubes();
/** Given a grid cell and an isovalue, this function calculates the triangular
   facets requied to represent the isosurface through the cell.
   Returns the number of triangular facets, the array "triangles"
   (see private section, below) will be loaded up with the
   vertices at most 5 triangular facets.
   0 will be returned if the grid cell is either totally above
   of totally below the isovalue.
*/
   int PolygonizeVoxel(GRIDCELL &g,const double iso);
/** Returns the point between two points in the same ratio as
   isolevel is between valp1 and valp2.
*/
   XYZ VertexInterp(const double isolevel,\
         const XYZ &p1,const XYZ &p2,const double valp1,const double valp2);
   TRIANGLE GetTriangle(size_t idx) {return auxtriangles[idx];}
   static vector<vector<double> > GetTriangleVertices(const PolygonizeMarchingCubes::TRIANGLE &t);
   static double Distance(const PolygonizeMarchingCubes::XYZ &p1,const PolygonizeMarchingCubes::XYZ &p2);
   static double Distance2(const PolygonizeMarchingCubes::XYZ &p1,const PolygonizeMarchingCubes::XYZ &p2);
/* ************************************************************************** */
   const static int edgeTable[256];
   const static int triTable[256][16];
protected:
   TRIANGLE auxtriangles[12];
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _POLYGONIZE_MARCHINGCUBES_H_ */

