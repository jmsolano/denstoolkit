/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/
#ifndef _MESHGRID_H_
#define _MESHGRID_H_
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;

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
   void ComputeNormal2FaceVectors();
   /** Compute the normal vector of a face. It assumes v_0, v_1 and v_2 (the vertices
    * v_i=face[idx][i]) produce outward normals. The function  < b> does not check < /b>
    * for idx<normal.size(); this should be checked elsewhere.  */
   void ComputeSingleNormal2Face(size_t idx);
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

/* ************************************************************************** */
class HelpersMeshGrid {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static void AddFaces2POVAsMesh(ofstream &ofil,\
         MeshGrid &grid,const double r, const double g, const double b,\
         const bool usenrmls,int usrntabs);
   static void AddFaces2POVAsMesh(ofstream &ofil,\
         MeshGrid &grid,const double gg,const bool usenrmls,int usrntabs) {
      return AddFaces2POVAsMesh(ofil,grid,gg,gg,gg,usenrmls,usrntabs); }
   static void AddVertices2POVAsSpheres(ofstream &ofil,\
         MeshGrid &grid,const double r,const double g,const double b,\
         const double sr,int usrntabs);
   static void AddVertices2POVAsSpheres(ofstream &ofil,\
         MeshGrid &grid,const double gg,const double sr,int usrntabs) {
      return AddVertices2POVAsSpheres(ofil,grid,gg,gg,gg,sr,usrntabs); }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */



#endif  /* _MESHGRID_H_ */

