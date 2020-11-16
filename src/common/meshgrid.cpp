#include <cstdlib>
#include <iostream>
using std::cout;
#include <string>
using std::string;
#include "meshgrid.h"
#include "matrixvectoroperations3d.h"

#ifndef MESHGRIDDEFCOORDSEPS2
#define  MESHGRIDDEFCOORDSEPS2 1.0e-12
#endif

MeshGrid::MeshGrid() :
   vertex(0),normal(0),centroid(0),face(0),value(0),cvalue(0),center(3),edge(0),vneigh2v(0) {
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
   cout << "After finding neighbours, there are " << vertex.size()
        << " vertices, " << edge.size() << " edges, and " << face.size()
        << " faces. Euler characteristic: " << (vertex.size()-edge.size()+face.size()) << '\n';
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

