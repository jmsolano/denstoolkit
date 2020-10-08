#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <memory>
using std::shared_ptr;
#include <algorithm>
#include "isosurface.h"
#include "povraytools.h"
#include "screenutils.h"
#include "palette.h"
#include "colorutils.h"

#ifndef BIGNUMBER4MINANDMAX
#define BIGNUMBER4MINANDMAX 1.0e0+100
#endif

Isosurface::Isosurface() {
   triangles.clear();
   prop2map.clear();
   isovalue=0.0e0;
   verbose=false;
   usecolormap=false;
   SetRGB(0.6,0.6,0.6);
   pmin=BIGNUMBER4MINANDMAX;
   pmax=-BIGNUMBER4MINANDMAX;
}
Isosurface::~Isosurface() {
   ResetTriangles();
}
void Isosurface::ResetTriangles() {
   triangles.clear();
}
bool Isosurface::ExtractMarchingCubes(const GaussianCube &g,const double isoval) {
   size_t nx=g.Nx()-2;
   size_t ny=g.Ny()-2;
   size_t nz=g.Nz()-2;
   if ( nx<1 || ny<1 || nz<1 ) {
      ScreenUtils::DisplayErrorMessage("Empty cube, cannot extract isosurface");
      return false;
   }
   isovalue=isoval;
   ResetTriangles();
   PolygonizeMarchingCubes::GRIDCELL voxel;
   PolygonizeMarchingCubes pol;
   int nt;
   for ( size_t i=0 ; i<nx ; ++i ) {
      if ( verbose && (i % (nx/10) ==0) ) {
         cout << "Slice " << i << "/" << nx << " [" <<
            size_t(int(100.0e0*double(i)/double(nx))) << "%]\r" << std::flush;
      }
      for ( size_t j=0 ; j<ny ; ++j ) {
         for ( size_t k=0 ; k<nz ; ++k ) {
            voxel=HelpersIsosurface::GetVoxel(g,i,j,k);
            nt=pol.PolygonizeVoxel(voxel,isovalue);
            for ( int l=0 ; l<nt ; ++l ) { triangles.push_back(pol.GetTriangle(l)); }
         }
      }
   }
   if ( verbose ) { cout << "Slice " << nx << "/" << nx << " [100%]" << endl; }
   ScaleTriangles(g.DX());
   return true;
}
void Isosurface::ScaleTriangles(vector<double> dx) {
   if ( dx.size()!=3 ) {
      ScreenUtils::DisplayErrorMessage("Not a 3D vector!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
   for ( size_t i=0 ; i<triangles.size() ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) {
         triangles[i].p[j].x*=dx[0];
         triangles[i].p[j].y*=dx[1];
         triangles[i].p[j].z*=dx[2];
      }
   }
}
void Isosurface::ComputeCentroids() {
   centroid.resize(triangles.size());
   for ( size_t i=0 ; i<centroid.size() ; ++i ) {
      centroid[i].resize(3);
      centroid[i][0]=triangles[i].p[0].x;
      centroid[i][1]=triangles[i].p[0].y;
      centroid[i][2]=triangles[i].p[0].z;
      centroid[i][0]+=triangles[i].p[1].x;
      centroid[i][1]+=triangles[i].p[1].y;
      centroid[i][2]+=triangles[i].p[1].z;
      centroid[i][0]+=triangles[i].p[2].x;
      centroid[i][1]+=triangles[i].p[2].y;
      centroid[i][2]+=triangles[i].p[2].z;
   }
   double oo3=1.0e0/3.0e0;
   for ( size_t i=0 ; i<centroid.size() ; ++i ) {
      centroid[i][0]*=oo3;
      centroid[i][1]*=oo3;
      centroid[i][2]*=oo3;
   }
   usecolormap=true;
}
bool Isosurface::SetProperty2Map(const vector<double> &p) {
   if ( p.size()!=triangles.size() ) {
      ScreenUtils::DisplayErrorMessage("Sizes of triangles and prop2map are different!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      usecolormap=false;
      return false;
   }
   prop2map=p;
   SearchMinAndMaxProp2Map();
   usecolormap=true;
   return true;
}
void Isosurface::SearchMinAndMaxProp2Map() {
   pmin=(*std::min_element(prop2map.begin(),prop2map.end()));
   pmax=(*std::max_element(prop2map.begin(),prop2map.end()));
}


/* ************************************************************************** */
/* ************************************************************************** */
PolygonizeMarchingCubes::GRIDCELL HelpersIsosurface::GetVoxel(\
      const GaussianCube &g,const size_t i,const size_t j,const size_t k) {
   PolygonizeMarchingCubes::GRIDCELL v;
   if ( /* i<0 || j<0 || k<0 || */ i>=(g.Nx()-1) || j>=(g.Ny()-1) || k>(g.Nz()-1) ) {
      ScreenUtils::DisplayErrorMessage("Out of bounds!");
      ScreenUtils::DisplayErrorMessage("Returning undefined voxel!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return v;
   }
   v.p[0].x = i;
   v.p[0].y = j;
   v.p[0].z = k;
   v.val[0] = g.Data(i,j,k);
   v.p[1].x = i+1;
   v.p[1].y = j;
   v.p[1].z = k; 
   v.val[1] = g.Data(i+1,j,k);
   v.p[2].x = i+1;
   v.p[2].y = j+1;
   v.p[2].z = k;
   v.val[2] = g.Data(i+1,j+1,k);
   v.p[3].x = i;
   v.p[3].y = j+1;
   v.p[3].z = k;
   v.val[3] = g.Data(i,j+1,k);
   v.p[4].x = i;
   v.p[4].y = j;
   v.p[4].z = k+1;
   v.val[4] = g.Data(i,j,k+1);
   v.p[5].x = i+1;
   v.p[5].y = j;
   v.p[5].z = k+1;
   v.val[5] = g.Data(i+1,j,k+1);
   v.p[6].x = i+1;
   v.p[6].y = j+1;
   v.p[6].z = k+1;
   v.val[6] = g.Data(i+1,j+1,k+1);
   v.p[7].x = i;
   v.p[7].y = j+1;
   v.p[7].z = k+1;
   v.val[7] = g.Data(i,j+1,k+1);
   return v;
}
bool HelpersIsosurface::AddIsosurfacePOV(ofstream &ofil,Isosurface &iso,int usrntabs) {
   int indlev=usrntabs;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "union {" << endl;
   PolygonizeMarchingCubes::TRIANGLE* t;
   for ( size_t i=0 ; i<iso.TrianglesSize() ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) {
         t=iso.GetTrianglePointer(i);
         HelpersPOVRay::WriteTriangle(ofil,indlev,\
               t->p[0].x,t->p[0].y,t->p[0].z,
               t->p[1].x,t->p[1].y,t->p[1].z,
               t->p[2].x,t->p[2].y,t->p[2].z,
               iso.rgb[0],iso.rgb[1],iso.rgb[2]);
      }
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;
   return true;
}
bool HelpersIsosurface::AddIsosurfacePOV(ofstream &ofil,Isosurface &iso,\
      shared_ptr<Palette> pal,int usrntabs) {
   if ( !iso.UseColorMap() ) { return AddIsosurfacePOV(ofil,iso,usrntabs); }
   int indlev=usrntabs;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "union {" << endl;
   PolygonizeMarchingCubes::TRIANGLE* t;
   if ( pal==nullptr ) { pal=std::make_shared<Palette>(); }
   pal->SetBlues();
   double pmin=iso.MinP2Map(),pmax=iso.MaxP2Map();
   double oodeltatimes255=255.0e0/(pmax-pmin);
   double tmp,r,g,b;
   size_t pos;
   for ( size_t i=0 ; i<iso.TrianglesSize() ; ++i ) {
      tmp=iso.Property2Map(i);
      if ( tmp<pmin ) { tmp=pmin; }
      if ( tmp>pmax ) { tmp=pmax; }
      tmp*=oodeltatimes255;
      if ( tmp>=0 ) {
         pos=size_t(tmp);
      } else {
         ScreenUtils::DisplayWarningMessage("tmp<0!");
         cout << __FILE__ << ", line: " << __LINE__ << '\n';
         pos=0;
      }
      pal->GetRGB(pos,r,g,b);
      //cout << pos << ": " << r << " " << g << " " << b << '\n';
      for ( size_t j=0 ; j<3 ; ++j ) {
         t=iso.GetTrianglePointer(i);
         HelpersPOVRay::WriteTriangle(ofil,indlev,\
               t->p[0].x,t->p[0].y,t->p[0].z,
               t->p[1].x,t->p[1].y,t->p[1].z,
               t->p[2].x,t->p[2].y,t->p[2].z,
               r,g,b);
      }
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;
   return true;
}

