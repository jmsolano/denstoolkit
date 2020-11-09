#ifndef _SYMMETRICSURFACEGRID_H_
#define _SYMMETRICSURFACEGRID_H_
#include <vector>
using std::vector;

/* ************************************************************************** */
class SymmetricSurfaceGrid {
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
   /** Adds the vector t to all vertices.  */
   void Scale(const double f);
   void ScaleSingleVertex(const size_t idx,const double f);
   void Translate(const vector<double> &t);
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
   size_t VertexPosition(const vector<double> &t);
   void TrimFacesCentroidDotProdGreaterThanZero(const vector<double> &d);
   void RemoveFacesUsingVertices(const vector<size_t> &v2r);
   void TriangleCentroidDir(size_t tFIdx,vector<double> &c);
   void ComputeCentroids();
   //void EraseVertexAndAssociatedFaces(size_t idx);
   void RemoveUnusedVertices();
   void DisplayFaces();
   void UseNormals(bool un) { usenormals=un; }
   inline bool UseNormals() { return usenormals; }
/* ************************************************************************** */
   vector<vector<double> > vertex;
   vector<vector<double> > normal;
   vector<vector<double> > tcentroid;
   vector<double> center;
   vector<vector<size_t> > tface; /*!< faces of triangles.  */
   vector<vector<size_t> > sface;/*!< faces of squares.  */
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   Shape shape;
   bool usenormals;
/* ************************************************************************** */
   void ClearArrays();
   void ResizeMatrix(vector<vector<double> > &m,size_t nr,size_t nc=3);
   void ResizeMatrix(vector<vector<size_t> > &m,size_t nr,size_t nc=3);
   template <class MM> void ClearMatrix(vector<vector<MM> > &m) {
      for ( size_t i=0 ; i<m.size() ; ++i ) { m[i].clear(); }; m.clear(); }
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _SYMMETRICSURFACEGRID_H_ */

