#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_
#include <fstream>
using std::ofstream;
#include "gaussiancube.h"
#include "polygonize_marchingcubes.h"
#include "palette.h"

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
   void UseColorMap(bool uc) { usecolormap=uc; }
   bool UseColorMap(void) const { return usecolormap; }
   void ComputeCentroids();
   /** Copies the vector v to an internal vector. This function
    * sets: UseColorMap(true).  */
   bool SetProperty2Map(const vector<double> &p);
   inline void SetMinP2Map(const double mm) { pmin=mm; }
   inline void SetMaxP2Map(const double mm) { pmax=mm; }
   inline double MinP2Map() { return pmin; }
   inline double MaxP2Map() { return pmax; }
   void SearchMinAndMaxProp2Map();
   inline double Property2Map(size_t idx) { return prop2map[idx]; }
   void SetRGB(double r,double g,double b) {rgb[0]=r; rgb[1]=g; rgb[2]=b;}
   size_t TrianglesSize() {return triangles.size();}
   PolygonizeMarchingCubes::TRIANGLE* GetTrianglePointer(size_t idx) {return &triangles[idx];}
/* ************************************************************************** */
   double rgb[3];
   vector<vector<double> > centroid;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void ResetTriangles();
   void ScaleTriangles(vector<double> dx);
   double isovalue;
/* ************************************************************************** */
   vector<PolygonizeMarchingCubes::TRIANGLE> triangles;
   /** Stores values that are used to color the triangles. This
    * vector is passed through SetProperty2Map(...).  */
   vector<double> prop2map;
   double pmin; /*!< Minimum value of the prop2map to be used in color scale.  */
   double pmax; /*!< Maximum value of the prop2map to be used in color scale.  */
   bool verbose;
   bool usecolormap;
/* ************************************************************************** */
};
/* ************************************************************************** */
/* ************************************************************************** */
class HelpersIsosurface {
public:
/* ************************************************************************** */
   static PolygonizeMarchingCubes::GRIDCELL GetVoxel(const GaussianCube &g,\
         const size_t i,const size_t j,const size_t k);
   /** As the name suggests, it adds the triangles of the isosurface to the povfile 'ofil'.
    * This function uses the RGB color passed through Isosurface::SetRGB.  */
   static bool AddIsosurfacePOV(ofstream &ofil,Isosurface &iso,int usrntabs=0);
   /** As the name suggests, it adds the triangles of the isosurface to the povfile 'ofil'.
    * This function colors the isosurface according to the Isosurface::prop2map values.
    * However, it will override this behaviour in case the boolean Isosurface::UseColorMap()
    * is false, if this is the case, then the function calls AddIsosurfacePOV withouth the
    * palette. If the parameter pal==nullptr, then the function will setup the default
    * palette (pal->SetBlues()). Otherwise, the function will use the configuration 
    * of pal. */
   static bool AddIsosurfacePOV(ofstream &ofil,Isosurface &iso,shared_ptr<Palette> pal,int usrntabs=0);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};


#endif  /* _ISOSURFACE_H_ */

