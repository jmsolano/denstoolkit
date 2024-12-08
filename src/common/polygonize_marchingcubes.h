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
   int PolygonizeVoxelTetrahedrons(GRIDCELL &g,const double iso);
/** Returns the point between two points in the same ratio as
   isolevel is between valp1 and valp2.
*/
   int PolygoniseTri(GRIDCELL g,double iso,TRIANGLE *tri,int v0,int v1,int v2,int v3);
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

