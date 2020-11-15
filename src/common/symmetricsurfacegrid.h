#ifndef _SYMMETRICSURFACEGRID_H_
#define _SYMMETRICSURFACEGRID_H_
#include <vector>
using std::vector;
#include "meshgrid.h"

/* ************************************************************************** */
class SymmetricSurfaceGrid : public MeshGrid {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   enum class Shape {
      NONE,//
      CUBOID,//Vertices 
      SPHERESQUARES, //A sphere whose vertices are determined through a projected cube
      SPHEREANGLES,//A sphere whose vertes are through angles theta and phi
      SPHEREICOSAHEDRON, //A sphere whose vertices are projected from an icosahedron
      CYLINDER
   };
   SymmetricSurfaceGrid();
   SymmetricSurfaceGrid(Shape sh);
   ~SymmetricSurfaceGrid();
   /** Projects the vertices of a cube into a unit sphere. If
    * it==0, then the mesh is a cube of edges equal to sqrt(2).
    * Each iteration divides all faces into 9 subfaces, i.e.,
    * each iteration divides a face into 3x3 smaller equal faces.  */
   void SetupSphereSquares(int it=0);
   /** Projects the vertices of the icosahedron onto a unit sphere. If
    * it==0, then the mesh is an icosahedron.  */
   void SetupSphereIcosahedron(const int it=0);
   void SubdivideAllSquares(int it=1);
   void SubdivideSquare(size_t idx);
   void SubdivideAllTriangles(int it=1);
   void SubdivideTriangle(size_t idx);
   void GenerateTrianglesFromSquares();
   void AddTrianglesFromSquare(size_t idx);
   void CheckTriangleOrientation();
   /** Returns the position of t within vertex. If t is not
    * in vertex, the function appends t to vertex and returns
    * vertex.size() before insertion [=vertex.size()-1 after insertion].  */
   void TrimFacesCentroidDotProdGreaterThanZero(const vector<double> &d);
   void RemoveFacesUsingVertices(const vector<size_t> &v2r);
   void ComputeCentroids();
   void RemoveUnusedVertices();
   void DisplayFaces();
/* ************************************************************************** */
   vector<vector<size_t> > sface;/*!< faces of squares.  */
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   Shape shape;
/* ************************************************************************** */
   void ClearArrays();
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _SYMMETRICSURFACEGRID_H_ */

