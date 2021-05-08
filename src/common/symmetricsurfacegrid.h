/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
   void TrimFacesCentroidDotProdGreaterThanZero(const vector<double> &d) {return TrimFacesCentroidDotProdBetweenVals(d,0.0e0,1.0e0); }
   void TrimFacesCentroidDotProdBetweenVals(const vector<double> &d,const double val1,const double val2);
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

