#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_
#include <fstream>
using std::ofstream;
#include <memory>
using std::shared_ptr;
#include "meshgrid.h"
#include "gaussiancube.h"
#include "polygonize_marchingcubes.h"
#include "palette.h"

/* ************************************************************************** */
class Isosurface : public MeshGrid {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Isosurface();
   ~Isosurface();
   bool Extract(const GaussianCube &g,const double isovalue) {return ExtractMarchingCubes(g,isovalue);}
   bool ExtractMarchingCubes(const GaussianCube &g,const double isoval);
   void UseColorMap(bool uc) { usecolormap=uc; }
   bool UseColorMap(void) const { return usecolormap; }
   void UseTetrahedrons(bool ut) { usetetrahedrons=ut; }
   inline bool UseTetrahedrons() { return usetetrahedrons; }
   /** Copies the vector v to an internal vector. This function
    * sets: UseColorMap(true).  */
   bool SetProperty2Map(const vector<double> &p);
   inline void SetMinP2Map(const double mm) { pmin=mm; }
   inline void SetMaxP2Map(const double mm) { pmax=mm; }
   inline double MinP2Map() { return pmin; }
   inline double MaxP2Map() { return pmax; }
   void SearchMinAndMaxProp2Map();
   inline double Property2Map(size_t idx) { return value[idx]; }
   void SetRGB(double r,double g,double b) {rgb[0]=r; rgb[1]=g; rgb[2]=b;}
   PolygonizeMarchingCubes::TRIANGLE* GetTrianglePointer(size_t idx) {return &triangles[idx];}
   void SetNormal(const size_t triangIdx,const size_t normalIdx,const vector<double> &v);
   inline vector<double> Vertex(const size_t tIdx,const size_t pIdx) const
      { return vector<double>{triangles[tIdx].p[pIdx].x,triangles[tIdx].p[pIdx].y,triangles[tIdx].p[pIdx].z}; }
   inline vector<double> Normal(const size_t tIdx,const size_t nIdx) const
      { return vector<double>{triangles[tIdx].n[nIdx].x,triangles[tIdx].n[nIdx].y,triangles[tIdx].n[nIdx].z}; }
/* ************************************************************************** */
   double rgb[3];
   void GenerateMesh();
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
   vector<double> value;
   double pmin; /*!< Minimum value of the value to be used in color scale.  */
   double pmax; /*!< Maximum value of the value to be used in color scale.  */
   bool usecolormap;
   bool usetetrahedrons;
/* ************************************************************************** */
};
/* ************************************************************************** */
/* ************************************************************************** */
class HelpersIsosurface {
public:
/* ************************************************************************** */
   static PolygonizeMarchingCubes::GRIDCELL GetVoxel(const GaussianCube &g,\
         const size_t i,const size_t j,const size_t k);
   /** Adds the triangle mesh of the isosurface to the povfile 'ofil'.
    * This function uses the RGB color passed through Isosurface::SetRGB. The isosurface is
    * colored with a single RGB color.  */
   static bool AddIsosurfacePOVMeshNoNormals(ofstream &ofil,Isosurface &iso,int usrntabs=0);
   /** Adds the triangle mesh of the isosurface to the povfile 'ofil'.
    * This function uses RGB colors for each triangle using the Isosurface::value values.
    * Each triangle is coloured with a single color. The colour is chosen
    * through the palette pal and Isosurface::value. If pal==nullptr, then the function will create a default
    * palette.  */
   static bool AddIsosurfacePOVMeshNoNormals(ofstream &ofil,Isosurface &iso,shared_ptr<Palette> pal,int usrntabs=0);
   /** As the name suggests, it adds the triangle mesh of the isosurface to the povfile 'ofil'.
    * This function colors the isosurface according to the Isosurface::value values.
    * However, it will override this behaviour in case the boolean Isosurface::UseColorMap()
    * is false, if this is the case, then the function calls AddIsosurfacePOVNoNormals withouth the
    * palette. If the parameter pal==nullptr, then the function will setup the default
    * palette (pal->SetBlues()). Otherwise, the function will use the configuration 
    * of pal. */
   static bool AddIsosurfacePOVMeshWithNormals(ofstream &ofil,Isosurface &iso,shared_ptr<Palette> pal,int usrntabs=0);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};


#endif  /* _ISOSURFACE_H_ */

