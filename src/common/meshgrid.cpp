/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
#include <cstdlib>
#include <iostream>
using std::cout;
#include <string>
using std::string;
#include "meshgrid.h"
#include "matrixvectoroperations3d.h"
#include "screenutils.h"
#include "povraytools.h"
#include "mymath.h"

#ifndef MESHGRIDDEFCOORDSEPS2
#define  MESHGRIDDEFCOORDSEPS2 1.0e-12
#endif

MeshGrid::MeshGrid() :
   vertex(0),normal(0),centroid(0),face(0),edge(0),vneigh2v(0),\
      value(0),cvalue(0),center(3) {
   verbose=true;
   usenormals=false;
}
MeshGrid::~MeshGrid() {
   ClearMatrix(vertex);   
   ClearMatrix(normal);   
   ClearMatrix(centroid);   
   ClearMatrix(face);
   ClearMatrix(edge);
   ClearMatrix(vneigh2v);
   value.clear();
}
void MeshGrid::ComputeCentroidsBase() {
   size_t nf=face.size(),p0,p1,p2;
   const double oo3=1.0e0/3.0e0;
   ResizeMatrix(centroid,nf,3);
   for ( size_t i=0 ; i<nf ; ++i ) {
      p0=face[i][0]; p1=face[i][1]; p2=face[i][2];
      centroid[i][0] =(vertex[p0][0]); centroid[i][1] =(vertex[p0][1]); centroid[i][2] =(vertex[p0][2]);
      centroid[i][0]+=(vertex[p1][0]); centroid[i][1]+=(vertex[p1][1]); centroid[i][2]+=(vertex[p1][2]);
      centroid[i][0]+=(vertex[p2][0]); centroid[i][1]+=(vertex[p2][1]); centroid[i][2]+=(vertex[p2][2]);
      centroid[i][0]*=oo3;             centroid[i][1]*=oo3;             centroid[i][2]*=oo3;
   }
}
void MeshGrid::TranslateBase(const vector<double> &t) {
   size_t nv=vertex.size();
   for ( size_t i=0 ; i<nv ; ++i ) {
      vertex[i][0]+=t[0];
      vertex[i][1]+=t[1];
      vertex[i][2]+=t[2];
   }
   for ( size_t i=0 ; i<3 ; ++i ) { center[i]+=t[i]; }
}
void MeshGrid::ScaleVerticesBase(const double f) {
   size_t nv=vertex.size();
   double dist=MatrixVectorOperations3D::Norm(center);
   if ( dist>0.0e0 ) {
      for ( size_t i=0 ; i<nv ; ++i ) {
         vertex[i][0]-=center[0];
         vertex[i][1]-=center[1];
         vertex[i][2]-=center[2];
      }
   }
   for ( size_t i=0 ; i<nv ; ++i ) {
      vertex[i][0]*=f;
      vertex[i][1]*=f;
      vertex[i][2]*=f;
   }
   if ( dist>0.0e0 ) {
      for ( size_t i=0 ; i<nv ; ++i ) {
         vertex[i][0]+=center[0];
         vertex[i][1]+=center[1];
         vertex[i][2]+=center[2];
      }
   }
}
void MeshGrid::ScaleSingleVertex(const size_t idx,const double f) {
   double dist=MatrixVectorOperations3D::Norm(center);
   if ( dist>0.0e0 ) {
      vertex[idx][0]-=center[0];
      vertex[idx][1]-=center[1];
      vertex[idx][2]-=center[2];
   }
   vertex[idx][0]*=f;
   vertex[idx][1]*=f;
   vertex[idx][2]*=f;
   if ( dist>0.0e0 ) {
      vertex[idx][0]+=center[0];
      vertex[idx][1]+=center[1];
      vertex[idx][2]+=center[2];
   }
}
void MeshGrid::NormaliseNormalsBase() {
   size_t nn=normal.size();
   for ( size_t i=0 ; i<nn ; ++i ) {
      MatrixVectorOperations3D::Normalize(normal[i]);
   }
}
size_t MeshGrid::VertexPosition(const vector<double> &t) {
   for ( size_t i=0 ; i<vertex.size() ; ++i ) {
      if ( MatrixVectorOperations3D::Distance2(t,vertex[i]) <= MESHGRIDDEFCOORDSEPS2 ) {
         return i;
      }
   }
   vertex.push_back(t);
   return (vertex.size()-1);
}
void MeshGrid::TriangleCentroidDir(size_t tFIdx,vector<double> &c) {
   size_t p0,p1,p2;
   const double oo3=1.0e0/3.0e0;
   p0=face[tFIdx][0]; p1=face[tFIdx][1]; p2=face[tFIdx][2];
   c[0] =(vertex[p0][0]); c[1] =(vertex[p0][1]); c[2] =(vertex[p0][2]);
   c[0]+=(vertex[p1][0]); c[1]+=(vertex[p1][1]); c[2]+=(vertex[p1][2]);
   c[0]+=(vertex[p2][0]); c[1]+=(vertex[p2][1]); c[2]+=(vertex[p2][2]);
   c[0]*=oo3;             c[1]*=oo3;             c[2]*=oo3;
   c[0]-=(center[0]);     c[1]-=(center[1]);     c[2]-=(center[2]);
}
void MeshGrid::RemoveUnusedVerticesBase(vector<vector<size_t> > &f,const size_t nvspf) {
   size_t nf=f.size();//number of faces.
   size_t pos=0,tmp;
   bool iamused;
   size_t nv=vertex.size();
   while ( pos<nv ) {
      iamused=false;
      for ( size_t i=0 ; i<nf ; ++i ) {
         for ( size_t j=0 ; j<nvspf ; ++j ) {
            if ( f[i][j]==pos ) { iamused=true; }
         }
         if ( iamused ) { break; }
      }
      if ( !iamused ) {
         tmp=nv-1;
         for ( size_t i=0 ; i<3 ; ++i ) { vertex[pos][i]=vertex[tmp][i]; }
         for ( size_t i=0 ; i<nf ; ++i ) {
            for ( size_t j=0 ; j<nvspf ; ++j ) {
               if ( f[i][j]==tmp ) { f[i][j]=pos; }
            }
         }
         vertex.pop_back();
         --nv;
      } else {
         ++pos;
      }
   }
   ResizeMatrix(normal,vertex.size(),3);
   value.resize(vertex.size());
}
void MeshGrid::DetermineEdgesBase(const vector<vector<size_t> > &f,size_t nvspf) {
   vector<size_t> tmp(2);
   edge.reserve(2*vertex.size()+f.size());
   edge.resize(0);
   bool imnew;
   size_t nvspfm1=nvspf-1;
   size_t nf=f.size();
   for ( size_t i=0 ; i<nf ; ++i ) {
      for ( size_t k=0 ; k<nvspfm1 ; ++k ) {
         tmp[0]=f[i][k]; tmp[1]=f[i][k+1]; if ( tmp[0]>tmp[1] ) { std::swap(tmp[0],tmp[1]); }
         imnew=true;
         for ( size_t j=0 ; j<edge.size() ; ++j ) {
            if ( tmp[0]==edge[j][0] && tmp[1]==edge[j][1] ) { imnew=false; break; }
         }
         if ( imnew ) { edge.push_back(tmp); }
      }
      tmp[0]=f[i][nvspfm1]; tmp[1]=f[i][0]; if ( tmp[0]>tmp[1] ) { std::swap(tmp[0],tmp[1]); }
      imnew=true;
      for ( size_t j=0 ; j<edge.size() ; ++j ) {
         if ( tmp[0]==edge[j][0] && tmp[1]==edge[j][1] ) { imnew=false; break; }
      }
      if ( imnew ) { edge.push_back(tmp); }
   }
}
void MeshGrid::ComputeNormal2FaceVectors() {
   size_t nf=face.size();
   if ( normal.size() != nf ) { ResizeMatrix(normal,nf,3); }
   for ( size_t i=0 ; i<nf ; ++i ) { ComputeSingleNormal2Face(i); }
}
void MeshGrid::ComputeSingleNormal2Face(size_t idx) {
   double ra[3],rb[3],rc[3];
   size_t ia=face[idx][0],ib=face[idx][1],ic=face[idx][2];
   ra[0]=vertex[ia][0]; ra[1]=vertex[ia][1]; ra[2]=vertex[ia][2];
   rb[0]=vertex[ib][0]; rb[1]=vertex[ib][1]; rb[2]=vertex[ib][2];
   rc[0]=vertex[ic][0]; rc[1]=vertex[ic][1]; rc[2]=vertex[ic][2];
   double v1[3],v2[3];
   v1[0]=rb[0]-ra[0]; v1[1]=rb[1]-ra[1]; v1[2]=rb[2]-ra[2];
   v2[0]=rc[0]-ra[0]; v2[1]=rc[1]-ra[1]; v2[2]=rc[2]-ra[2];
   crossProductV3(v1,v2,rc);
   normal[idx][0]=rc[0]; normal[idx][1]=rc[1]; normal[idx][2]=rc[2];
}
void MeshGrid::FindAllVertexNeighbours() {
   if ( edge.size()==0 || edge.size()==string::npos ) { DetermineEdges(); }
   ClearMatrix(vneigh2v);
   vneigh2v.resize(vertex.size());
   for ( size_t i=0 ; i<vneigh2v.size() ; ++i ) { vneigh2v.reserve(12); }
   size_t ne=edge.size();
   for ( size_t i=0 ; i<ne ; ++i ) {
      vneigh2v[edge[i][0]].push_back(edge[i][1]);
      vneigh2v[edge[i][1]].push_back(edge[i][0]);
      //vneigh2v[edge[i][0]].push_back(i);
      //vneigh2v[edge[i][1]].push_back(i);
   }
   if ( verbose ) {
      cout << "After finding neighbours, there are " << vertex.size()
           << " vertices, " << edge.size() << " edges, and " << face.size()
           << " faces. Euler characteristic: " << (vertex.size()-edge.size()+face.size()) << '\n';
   }
   /*
   for ( size_t i=0 ; i<edge.size() ; ++i ) {
      cout << "Edge (" << i << "): ";
      for ( size_t j=0 ; j<edge[i].size() ; ++j ) {
         cout << ' ' << edge[i][j];
      }
      cout << '\n';
   }
   // */
}
/* ************************************************************************** */
/* ************************************************************************** */
void HelpersMeshGrid::AddFaces2POVAsMesh(ofstream &ofil,\
      MeshGrid &grid, const double r,const double g,const double b,\
      const bool usenrmls,int usrntabs) {
   int indlev=usrntabs;
   if ( grid.vertex.size() == 0 ) {
      ScreenUtils::DisplayErrorMessage("No vertices in the mesh!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   size_t nvm1=grid.vertex.size()-1;
   size_t nfm1=grid.face.size()-1;
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "mesh2 {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "vertex_vectors {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,grid.vertex[i][0],grid.vertex[i][1],grid.vertex[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,grid.vertex[nvm1][0],grid.vertex[nvm1][1],grid.vertex[nvm1][2]);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << '\n' << thetabs << "}\n";//end of vertex_vectors
   if ( usenrmls ) {
      if ( grid.normal.size()<1 ) {
         ScreenUtils::DisplayWarningMessage("There are no normals!");
         cout << __FILE__ << ", line: " << __LINE__ << '\n';
      } else {
         ofil << thetabs << "normal_vectors {\n";
         thetabs=HelpersPOVRay::IndTabsStr(indlev);
         ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
         for ( size_t i=0 ; i<nvm1 ; ++i ) {
            HelpersPOVRay::WriteVector(ofil,grid.normal[i][0],grid.normal[i][1],grid.normal[i][2]);
            ofil << ",";
            if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
         }
         HelpersPOVRay::WriteVector(ofil,grid.normal[nvm1][0],grid.normal[nvm1][1],grid.normal[nvm1][2]);
         thetabs=HelpersPOVRay::IndTabsStr(--indlev);
         ofil << thetabs << "}\n";//end of normal_vectors
      }
   }
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "texture_list {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   size_t nnvv=nvm1+1;
   ofil << thetabs << (nnvv) << ",\n";
   for ( size_t i=0 ; i<nnvv ; ++i ) {
      ofil << thetabs << "texture{pigment{rgb ";
      HelpersPOVRay::WriteVector(ofil,r,g,b);
      ofil << " } finish {ambient 0 emission  0.6}}" << (i<nvm1? ',' : ' ') << "\n";
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}\n";//end of texture_list
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "face_indices {\n";
   thetabs=HelpersPOVRay::IndTabsStr(indlev);
   ofil << thetabs << (nfm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nfm1 ; ++i ) {
      HelpersPOVRay::WriteVector(ofil,grid.face[i][0],grid.face[i][1],grid.face[i][2]);
      ofil << ", " << grid.face[i][0] << ',' << grid.face[i][1] << ',' << grid.face[i][2] << ", ";
      if ( (i%5) == 4 ) { ofil << '\n' << thetabs; }
   }
   HelpersPOVRay::WriteVector(ofil,grid.face[nfm1][0],grid.face[nfm1][1],grid.face[nfm1][2]);
   ofil << ", " <<  grid.face[nfm1][0] << ',' << grid.face[nfm1][1] << ',' << grid.face[nfm1][2] << '\n';
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}\n";//end of face_indices
   thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "pigment { rgb ";
   HelpersPOVRay::WriteVector(ofil,r,g,b);
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}\n";//end of pigment
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}\n"; //end of mesh2
}
void HelpersMeshGrid::AddVertices2POVAsSpheres(ofstream &ofil,\
         MeshGrid &grid,const double r,const double g,const double b,\
         const double sr,int usrntabs) {
   int indlev=usrntabs;
   if ( grid.vertex.size() == 0 ) {
      ScreenUtils::DisplayErrorMessage("No vertices in the mesh!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   size_t nv=grid.vertex.size();
   string thetabs=HelpersPOVRay::IndTabsStr(indlev++);
   ofil << thetabs << "union {\n";
   for ( size_t i=0 ; i<nv ; ++i ) {
      HelpersPOVRay::WriteSphere(ofil,indlev,\
            grid.vertex[i][0],grid.vertex[i][1],grid.vertex[i][2],sr,r,g,b);
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}\n";//end of union (vertices as spheres).
}

