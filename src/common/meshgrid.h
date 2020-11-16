#ifndef _MESHGRID_H_
#define _MESHGRID_H_
#include <vector>
using std::vector;

/* ************************************************************************** */
/** Instances of this class contains vertices, faces, centroids,
 * and normals of a triangular mesh. The class is intended to provide
 * data to generate povray's mesh2 objects and to serve as a base class
 * of triangular meshgrids. The class also has the most common
 * functions to manipulate the mesh (i.e. the vertices, faces, normals,
 * and centroids). This class assumes that each face is a triangle.
 * The MeshGrid object also has a vector to contain a scalar field
 * evalutated at each vertex (value).
 * The vectors/matrices (vectors of size 3/or arrays of size-3 vectors)
 * centroid, vertex, and normal are public, as well as the vector
 * value, in order to easy the data access. The user must be careful about
 * writing on this data objects.
 *   */
class MeshGrid {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   MeshGrid();
   ~MeshGrid();
   void Verbose(void) {verbose=true;}
   void Quiet(void) {verbose=false;}
   void UseNormals(bool un) { usenormals=un; }
   inline bool UseNormals() { return usenormals; }
   void ComputeCentroidsBase();
   virtual void ComputeCentroids() { return ComputeCentroidsBase(); }
   /** This translates triangles and centroids by t (i.e., it adds t to all coodinates).
    * This function also adds t to center.  */
   void TranslateBase(const vector<double> &t);
   virtual void Translate(const vector<double> &t) { return TranslateBase(t); }
   /** This function translate the vertices to the origin (i.e., it
    * moves the grid to the initial position), subsequently scales the vertices
    * by the factor f, and finally moves back the mesh to center.  */
   void ScaleVerticesBase(const double f);
   virtual void ScaleVertices(const double f) { return ScaleVerticesBase(f); }
   /** This function translates the vertes of index idx
    * to the origin, subsequently scales the vertex by the
    * factor f, and finally translates it back to center.  */
   void ScaleSingleVertex(const size_t idx,const double f);
   void NormaliseNormalsBase();
   virtual void NormaliseNormals() { return NormaliseNormalsBase(); }
   /** Returns the position of t within vertex. If t is not
    * in vertex, the function appends t to vertex and returns
    * vertex.size() before insertion [=vertex.size()-1 after insertion].  */
   size_t VertexPosition(const vector<double> &t);
   void TriangleCentroidDir(size_t tFIdx,vector<double> &c);
   void RemoveUnusedVerticesBase(vector<vector<size_t> > &f,\
         const size_t nvspf=3 /** number of vertices per face  */);
   virtual void RemoveUnusedVertices() { return RemoveUnusedVerticesBase(face,3); }
   virtual void DetermineEdges() {return DetermineEdgesBase(face,3); }
   void DetermineEdgesBase(const vector<vector<size_t> > &f,size_t nvspf);
   void FindAllVertexNeighbours();
/* ************************************************************************** */
   vector<vector<double> > vertex;/*!< Vertices container.  */
   vector<vector<double> > normal;/*!< Normals container.  */
   vector<vector<double> > centroid;/*!< Centroids container.  */
   vector<vector<size_t> > face; /*!< Faces container.  */
   vector<vector<size_t> > edge;/*!< Edges container.  */
   vector<vector<size_t> > vneigh2v;/*!< The neighbour vertices associated with each vertex.  */
   vector<double> value;/*!< The field values container (values at vertices).  */
   vector<double> cvalue;/*!< Another field values container (values at centroids).  */
   vector<double> center;/*!< The center of the mesh.  */
/* ************************************************************************** */
   template <class MM> void ClearMatrix(vector<vector<MM> > &m) {
      for ( size_t i=0 ; i<m.size() ; ++i ) { m[i].clear(); }; m.clear(); }
   /** Resizes the matrix m to the specified dimensions. This does not keep
    * the previous data. nr is the number of rows, and nc is the number of
    * columns.  */
   template <class MM> void ResizeMatrix(vector<vector<MM> > &m,size_t nr,size_t nc=3) {
       for ( size_t i=0 ; i<m.size() ; ++i ) { m[i].clear(); }
       m.resize(nr); for ( size_t i=0 ; i<nr ; ++i ) { m[i].resize(nc); }
   }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   bool verbose,usenormals;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _MESHGRID_H_ */

