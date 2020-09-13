#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_
#include <fstream>
using std::ofstream;
#include "gaussiancube.h"
#include "polygonize_marchingcubes.h"

/* ************************************************************************** */
class Isosurface {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Isosurface();
   ~Isosurface();
   bool Extract(const GaussianCube &g,const double isovalue) {return ExtractMarchingCubes(g,isovalue);}
   bool ExtractMarchingCubes(const GaussianCube &g,const double isoval);
   void Verbose(void) {verbose=true;}
   void Quiet(void) {verbose=false;}
   void SetRGB(double r,double g,double b) {rgb[0]=r; rgb[1]=g; rgb[2]=b;}
   size_t TrianglesSize() {return triangles.size();}
   PolygonizeMarchingCubes::TRIANGLE* GetTrianglePointer(size_t idx) {return &triangles[idx];}
/* ************************************************************************** */
   double rgb[3];
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void ResetTriangles();
   void ScaleTriangles(vector<double> dx);
   double isovalue;
/* ************************************************************************** */
   vector<PolygonizeMarchingCubes::TRIANGLE> triangles;
   bool verbose;
/* ************************************************************************** */
};
/* ************************************************************************** */
/* ************************************************************************** */
class HelpersIsosurface {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static PolygonizeMarchingCubes::GRIDCELL GetVoxel(const GaussianCube &g,\
         const size_t i,const size_t j,const size_t k);
   static bool AddIsosurfacePOV (ofstream &ofil,Isosurface &iso,int usrntabs=0);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _ISOSURFACE_H_ */

