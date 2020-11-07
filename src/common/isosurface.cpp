#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
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
   vertices.clear();
   faces.clear();
   isovalue=0.0e0;
   verbose=false;
   usecolormap=false;
   usetetrahedrons=false;
   SetRGB(0.6,0.6,0.6);
   pmin=BIGNUMBER4MINANDMAX;
   pmax=-BIGNUMBER4MINANDMAX;
}
Isosurface::~Isosurface() {
   ResetTriangles();
}
void Isosurface::ResetTriangles() {
   triangles.clear();
   prop2map.clear();
   faces.clear();
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
      if (usetetrahedrons) {
         for ( size_t j=0 ; j<ny ; ++j ) {
            for ( size_t k=0 ; k<nz ; ++k ) {
               voxel=HelpersIsosurface::GetVoxel(g,i,j,k);
               nt=pol.PolygonizeVoxelTetrahedrons(voxel,isovalue);
               for ( int l=0 ; l<nt ; ++l ) { triangles.push_back(pol.GetTriangle(l)); }
            }
         }
      } else {
         for ( size_t j=0 ; j<ny ; ++j ) {
            for ( size_t k=0 ; k<nz ; ++k ) {
               voxel=HelpersIsosurface::GetVoxel(g,i,j,k);
               nt=pol.PolygonizeVoxel(voxel,isovalue);
               for ( int l=0 ; l<nt ; ++l ) { triangles.push_back(pol.GetTriangle(l)); }
            }
         }
      }
   }
   if ( verbose ) { cout << "Slice " << nx << "/" << nx << " [100%]" << endl; }
   ScaleTriangles(g.DX());
   vector<double> t=g.X0();
   for ( size_t i=0 ; i<t.size() ; ++i ) { t[i]=-t[i]; }
   Translate(t);
   GenerateMesh();
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
   if ( p.size()!=triangles.size() && p.size()!=vertices.size() ) {
      ScreenUtils::DisplayErrorMessage("Sizes of triangles/vertices and prop2map are different!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      usecolormap=false;
      return false;
   }
   prop2map=p;
   SearchMinAndMaxProp2Map();
   usecolormap=true;
   return true;
}
void Isosurface::SetNormal(const size_t triangIdx,const size_t normalIdx,\
      const vector<double> &v) {
   triangles[triangIdx].n[normalIdx].x=v[0];
   triangles[triangIdx].n[normalIdx].y=v[1];
   triangles[triangIdx].n[normalIdx].z=v[2];
}
void Isosurface::SearchMinAndMaxProp2Map() {
   pmin=(*std::min_element(prop2map.begin(),prop2map.end()));
   pmax=(*std::max_element(prop2map.begin(),prop2map.end()));
}
void Isosurface::Translate(const vector<double> &t) {
   if ( t.size()!=3 ) {
      ScreenUtils::DisplayErrorMessage("The translation vector must have"
            "dimension 3!\n Nothing done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   size_t nn=triangles.size();
   for ( size_t i=0 ; i<nn ; ++i ) {
      triangles[i].p[0].x-=t[0];
      triangles[i].p[0].y-=t[1];
      triangles[i].p[0].z-=t[2];
      triangles[i].p[1].x-=t[0];
      triangles[i].p[1].y-=t[1];
      triangles[i].p[1].z-=t[2];
      triangles[i].p[2].x-=t[0];
      triangles[i].p[2].y-=t[1];
      triangles[i].p[2].z-=t[2];
   }
   nn=centroid.size();
   for ( size_t i=0 ; i<nn ; ++i ) {
      centroid[i][0]-=t[0];
      centroid[i][1]-=t[1];
      centroid[i][2]-=t[2];
   }
}
void Isosurface::NormalizeNormals() {
   double norm;
   PolygonizeMarchingCubes::XYZ* n;
   for ( size_t i=0 ; i<triangles.size() ; ++i ) {
      n=&(triangles[i].n[0]);
      norm=1.0e0/sqrt(((n->x)*(n->x))+((n->y)*(n->y))+((n->z)*(n->z)));
      (n->x)*=norm; (n->y)*=norm; (n->z)*=norm;
      n=&(triangles[i].n[1]);
      norm=1.0e0/sqrt(((n->x)*(n->x))+((n->y)*(n->y))+((n->z)*(n->z)));
      (n->x)*=norm; (n->y)*=norm; (n->z)*=norm;
      n=&(triangles[i].n[2]);
      norm=1.0e0/sqrt(((n->x)*(n->x))+((n->y)*(n->y))+((n->z)*(n->z)));
      (n->x)*=norm; (n->y)*=norm; (n->z)*=norm;
   }
}
void Isosurface::GenerateMesh() {
   size_t nt=triangles.size();
   if ( nt<1 ) { return; }
   //cout << "There were originally: " << nt*3 << "vertices." << '\n';
   faces.resize(nt);
   for ( size_t i=0 ; i<nt ; ++i ) { faces[i].resize(3); }
   vertices.clear();
   PolygonizeMarchingCubes::XYZ* v;
   vector<double> tmp(3);
   vector<size_t> face(3);
   bool matched;
   v=&(triangles[0].p[0]); face[0]=0;
   tmp[0]=v->x; tmp[1]=v->y; tmp[2]=v->z;
   vertices.push_back(tmp); 
   v=&(triangles[0].p[1]); face[1]=1;
   tmp[0]=v->x; tmp[1]=v->y; tmp[2]=v->z;
   vertices.push_back(tmp); 
   v=&(triangles[0].p[2]); face[2]=2;
   tmp[0]=v->x; tmp[1]=v->y; tmp[2]=v->z;
   vertices.push_back(tmp); 
   faces[0]=face;
   for ( size_t i=1 ; i<nt ; ++i ) {
      for ( size_t k=0 ; k<3 ; ++k ) {
         v=&(triangles[i].p[k]);
         tmp[0]=v->x; tmp[1]=v->y; tmp[2]=v->z;
         face[k]=vertices.size();   
         matched=false;
         for ( size_t j=0 ; j<face[k] ; ++j ) {
            if ( ( (tmp[0]-vertices[j][0])*(tmp[0]-vertices[j][0])\
                     +(tmp[1]-vertices[j][1])*(tmp[1]-vertices[j][1])\
                     +(tmp[2]-vertices[j][2])*(tmp[2]-vertices[j][2]))<EPSVERTEXDIST2 ) {
               face[k]=j;
               matched=true;
               break;
            }
         }
         if ( !matched ) { vertices.push_back(tmp); }
      }
      faces[i]=face;
   }
   //cout << "There are " << vertices.size() << " unique vertices." << '\n';
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
/*
bool HelpersIsosurface::AddIsosurfacePOVNoNormals(ofstream &ofil,Isosurface &iso,int usrntabs) {
   int indlev=usrntabs;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "union {" << endl;
   PolygonizeMarchingCubes::TRIANGLE* t;
   for ( size_t i=0 ; i<iso.TrianglesSize() ; ++i ) {
         t=iso.GetTrianglePointer(i);
         HelpersPOVRay::WriteTriangle(ofil,indlev,\
               t->p[0].x,t->p[0].y,t->p[0].z,
               t->p[1].x,t->p[1].y,t->p[1].z,
               t->p[2].x,t->p[2].y,t->p[2].z,
               iso.rgb[0],iso.rgb[1],iso.rgb[2]);
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;
   return true;
}
// */
bool HelpersIsosurface::AddIsosurfacePOVMeshNoNormals(ofstream &ofil,Isosurface &iso,int usrntabs) {
   int indlev=usrntabs;
   size_t nvm1=iso.vertices.size()-1;
   size_t nfm1=iso.faces.size()-1;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "mesh2 {" << endl;
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "vertex_vectors {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.vertices[i][0],iso.vertices[i][1],iso.vertices[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.vertices[nvm1][0],iso.vertices[nvm1][1],iso.vertices[nvm1][2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of vertex_vectors
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "face_indices {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nfm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nfm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.faces[i][0],iso.faces[i][1],iso.faces[i][2]);
      ofil << ",";
      if ( (i%5) == 4 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.faces[nfm1][0],iso.faces[nfm1][1],iso.faces[nfm1][2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of face_indices
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "pigment { rgb ";
   HelpersPOVRay::WriteVector(ofil,iso.rgb[0],iso.rgb[1],iso.rgb[2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of pigment
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl; //end of mesh2
   return true;
}
bool HelpersIsosurface::AddIsosurfacePOVMeshNoNormals(ofstream &ofil,Isosurface &iso,\
      shared_ptr<Palette> pal,int usrntabs) {
   if ( !iso.UseColorMap() ) {
      ScreenUtils::DisplayWarningMessage("UseColorMap()==false. Calling"\
            "\nHelpersIsosurface::AddIsosurfacePOVMeshNoNormals using single rgb.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return AddIsosurfacePOVMeshNoNormals(ofil,iso,usrntabs);
   }
   if ( pal==nullptr ) {
      pal=std::make_shared<Palette>();
      pal->SetBlues();
   }
   double pmin=iso.MinP2Map(),pmax=iso.MaxP2Map();
   double oodeltatimes255=255.0e0/(pmax-pmin);
   double tmp,r,g,b;
   size_t pos;
   int indlev=usrntabs;
   size_t nvm1=iso.vertices.size()-1;
   size_t nfm1=iso.faces.size()-1;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "mesh2 {" << endl;
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "vertex_vectors {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.vertices[i][0],iso.vertices[i][1],iso.vertices[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.vertices[nvm1][0],iso.vertices[nvm1][1],iso.vertices[nvm1][2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of vertex_vectors
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "texture_list {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   size_t nnvv=nvm1+1;
   ofil << thetabs << (nnvv) << ",\n";
   for ( size_t i=0 ; i<nnvv ; ++i ) {
      tmp=iso.Property2Map(i);
      if ( tmp<pmin ) { tmp=pmin; }
      if ( tmp>pmax ) { tmp=pmax; }
      tmp-=pmin;
      tmp*=oodeltatimes255;
      if ( tmp>=0 ) {
         pos=size_t(tmp);
      } else {
         ScreenUtils::DisplayWarningMessage("tmp<0!");
         cout << __FILE__ << ", line: " << __LINE__ << '\n';
         pos=0;
      }
      pal->GetRGB(pos,r,g,b);
      ofil << thetabs << "texture{pigment{rgb ";
      HelpersPOVRay::WriteVector(ofil,r,g,b);
      ofil << " }}" << (i<nvm1? ',' : ' ') << "\n";
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of texture_list
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "face_indices {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nfm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nfm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.faces[i][0],iso.faces[i][1],iso.faces[i][2]);
      ofil << ", " << i%(nvm1+1) << ", ";
      if ( (i%5) == 4 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.faces[nfm1][0],iso.faces[nfm1][1],iso.faces[nfm1][2]);
   ofil << ',' << nvm1 << '\n';
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of face_indices
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "pigment { rgb ";
   HelpersPOVRay::WriteVector(ofil,iso.rgb[0],iso.rgb[1],iso.rgb[2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of pigment
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl; //end of mesh2
   return true;
}
bool HelpersIsosurface::AddIsosurfacePOVMeshWithNormals(ofstream &ofil,Isosurface &iso,\
      shared_ptr<Palette> pal,int usrntabs) {
   if ( !iso.UseColorMap() ) {
      ScreenUtils::DisplayWarningMessage("UseColorMap()==false. Calling"\
            "\nHelpersIsosurface::AddIsosurfacePOVMeshNoNormals using single rgb.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return AddIsosurfacePOVMeshNoNormals(ofil,iso,usrntabs);
   }
   if ( pal==nullptr ) {
      pal=std::make_shared<Palette>();
      pal->SetBlues();
   }
   double pmin=iso.MinP2Map(),pmax=iso.MaxP2Map();
   double oodeltatimes255=255.0e0/(pmax-pmin);
   double tmp,r,g,b;
   size_t pos;
   int indlev=usrntabs;
   size_t nvm1=iso.vertices.size()-1;
   size_t nfm1=iso.faces.size()-1;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "mesh2 {" << endl;
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "vertex_vectors {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.vertices[i][0],iso.vertices[i][1],iso.vertices[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.vertices[nvm1][0],iso.vertices[nvm1][1],iso.vertices[nvm1][2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of vertex_vectors
   ofil << thetabs << "normal_vectors {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.normals[i][0],iso.normals[i][1],iso.normals[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.normals[nvm1][0],iso.normals[nvm1][1],iso.normals[nvm1][2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of normal_vectors
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "texture_list {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   size_t nnvv=nvm1+1;
   ofil << thetabs << (nnvv) << ",\n";
   for ( size_t i=0 ; i<nnvv ; ++i ) {
      tmp=iso.Property2Map(i);
      if ( tmp<pmin ) { tmp=pmin; }
      if ( tmp>pmax ) { tmp=pmax; }
      tmp-=pmin;
      tmp*=oodeltatimes255;
      if ( tmp>=0 ) {
         pos=size_t(tmp);
      } else {
         ScreenUtils::DisplayWarningMessage("tmp<0!");
         cout << __FILE__ << ", line: " << __LINE__ << '\n';
         pos=0;
      }
      pal->GetRGB(pos,r,g,b);
      ofil << thetabs << "texture{pigment{rgb ";
      HelpersPOVRay::WriteVector(ofil,r,g,b);
      ofil << " }}" << (i<nvm1? ',' : ' ') << "\n";
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of texture_list
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "face_indices {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nfm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nfm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,iso.faces[i][0],iso.faces[i][1],iso.faces[i][2]);
      ofil << ", " << iso.faces[i][0] << ',' << iso.faces[i][1] << ',' << iso.faces[i][2] << ", ";
      if ( (i%5) == 4 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,iso.faces[nfm1][0],iso.faces[nfm1][1],iso.faces[nfm1][2]);
   ofil << ", " <<  iso.faces[nfm1][0] << ',' << iso.faces[nfm1][1] << ',' << iso.faces[nfm1][2] << '\n';
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of face_indices
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "pigment { rgb ";
   HelpersPOVRay::WriteVector(ofil,iso.rgb[0],iso.rgb[1],iso.rgb[2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of pigment
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl; //end of mesh2
   return true;
}
/*
bool HelpersIsosurface::AddIsosurfacePOVWithNormals(ofstream &ofil,Isosurface &iso,int usrntabs) {
   int indlev=usrntabs;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "union {" << endl;
   PolygonizeMarchingCubes::TRIANGLE* t;
   vector<vector<double> > v(3),n(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      v[i].resize(3);
      n[i].resize(3);
   }
   for ( size_t i=0 ; i<iso.TrianglesSize() ; ++i ) {
         t=iso.GetTrianglePointer(i);
         v[0][0]=t->p[0].x; v[0][1]=t->p[0].x; v[0][2]=t->p[0].x;
         v[1][0]=t->p[1].y; v[1][1]=t->p[1].y; v[1][2]=t->p[1].y;
         v[2][0]=t->p[2].z; v[2][1]=t->p[2].z; v[2][2]=t->p[2].z;
         n[0][0]=t->n[0].x; n[0][1]=t->n[0].x; n[0][2]=t->n[0].x;
         n[1][0]=t->n[1].y; n[1][1]=t->n[1].y; n[1][2]=t->n[1].y;
         n[2][0]=t->n[2].z; n[2][1]=t->n[2].z; n[2][2]=t->n[2].z;
         HelpersPOVRay::WriteSmoothTriangle(ofil,indlev,\
               v,n,iso.rgb[0],iso.rgb[1],iso.rgb[2]);
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;
   return true;
}
// */
/*
bool HelpersIsosurface::AddIsosurfacePOVNoNormals(ofstream &ofil,Isosurface &iso,\
      shared_ptr<Palette> pal,int usrntabs) {
   if ( !iso.UseColorMap() ) { return AddIsosurfacePOVNoNormals(ofil,iso,usrntabs); }
   int indlev=usrntabs;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "union {" << endl;
   PolygonizeMarchingCubes::TRIANGLE* t;
   if ( pal==nullptr ) {
      pal=std::make_shared<Palette>();
      pal->SetBlues();
   }
   double pmin=iso.MinP2Map(),pmax=iso.MaxP2Map();
   double oodeltatimes255=255.0e0/(pmax-pmin);
   double tmp,r,g,b;
   size_t pos;
   for ( size_t i=0 ; i<iso.TrianglesSize() ; ++i ) {
      tmp=iso.Property2Map(i);
      if ( tmp<pmin ) { tmp=pmin; }
      if ( tmp>pmax ) { tmp=pmax; }
      tmp-=pmin;
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
      t=iso.GetTrianglePointer(i);
      HelpersPOVRay::WriteTriangle(ofil,indlev,\
            t->p[0].x,t->p[0].y,t->p[0].z,
            t->p[1].x,t->p[1].y,t->p[1].z,
            t->p[2].x,t->p[2].y,t->p[2].z,
            r,g,b);
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;
   return true;
}
// */
/*
bool HelpersIsosurface::AddIsosurfacePOVWithNormals(ofstream &ofil,Isosurface &iso,\
      shared_ptr<Palette> pal,int usrntabs) {
   if ( !iso.UseColorMap() ) { return AddIsosurfacePOVNoNormals(ofil,iso,usrntabs); }
   int indlev=usrntabs;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "merge {" << endl;
   PolygonizeMarchingCubes::TRIANGLE* t;
   if ( pal==nullptr ) {
      pal=std::make_shared<Palette>();
      pal->SetBlues();
   }
   vector<vector<double> > v(3),n(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      v[i].resize(3);
      n[i].resize(3);
   }
   double pmin=iso.MinP2Map(),pmax=iso.MaxP2Map();
   double oodeltatimes255=255.0e0/(pmax-pmin);
   double tmp,r,g,b;
   size_t pos;
   for ( size_t i=0 ; i<iso.TrianglesSize() ; ++i ) {
      tmp=iso.Property2Map(i);
      if ( tmp<pmin ) { tmp=pmin; }
      if ( tmp>pmax ) { tmp=pmax; }
      tmp-=pmin;
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
      t=iso.GetTrianglePointer(i);
      v[0][0]=t->p[0].x; v[0][1]=t->p[0].y; v[0][2]=t->p[0].z;
      v[1][0]=t->p[1].x; v[1][1]=t->p[1].y; v[1][2]=t->p[1].z;
      v[2][0]=t->p[2].x; v[2][1]=t->p[2].y; v[2][2]=t->p[2].z;
      n[0][0]=t->n[0].x; n[0][1]=t->n[0].y; n[0][2]=t->n[0].z;
      n[1][0]=t->n[1].x; n[1][1]=t->n[1].y; n[1][2]=t->n[1].z;
      n[2][0]=t->n[2].x; n[2][1]=t->n[2].y; n[2][2]=t->n[2].z;
      HelpersPOVRay::WriteSmoothTriangle(ofil,indlev,v,n,r,g,b);
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;
   return true;
}
// */
