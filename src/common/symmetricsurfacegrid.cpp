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
using std::endl;
using std::cerr;
#include <cmath>
#include <algorithm>
#include "screenutils.h"
#include "matrixvectoroperations3d.h"
#include "symmetricsurfacegrid.h"
#ifndef SYMMSURFGRIDDEFCOORDSEPS2
#define SYMMSURFGRIDDEFCOORDSEPS2 1.0e-12
#endif

SymmetricSurfaceGrid::SymmetricSurfaceGrid() : MeshGrid() {
   shape=Shape::NONE;
   centroid.clear();
   face.clear();
   sface.clear();
   center.resize(3);
   center[0]=center[1]=center[2]=0.0e0;
   usenormals=false;
}
SymmetricSurfaceGrid::SymmetricSurfaceGrid(Shape sh)
   : SymmetricSurfaceGrid() {
   shape=sh;
}
SymmetricSurfaceGrid::~SymmetricSurfaceGrid() {
   ClearArrays();
}
void SymmetricSurfaceGrid::ClearArrays() {
   ClearMatrix(vertex);
   ClearMatrix(normal);
   ClearMatrix(centroid);
   ClearMatrix(face);
   ClearMatrix(sface);
}
void SymmetricSurfaceGrid::SetupSphereSquares(int it) {
   if ( vertex.size()!=0 ) {
      ClearArrays();
      shape=Shape::NONE;
   }
   vertex.resize(8);
   static const double ooSq3=0.577350269189625764509149; //1/sqrt(3)
   vertex[0]={ -ooSq3, -ooSq3, -ooSq3 };
   vertex[1]={  ooSq3, -ooSq3, -ooSq3 };
   vertex[2]={  ooSq3,  ooSq3, -ooSq3 };
   vertex[3]={ -ooSq3,  ooSq3, -ooSq3 };
   vertex[4]={ -ooSq3, -ooSq3,  ooSq3 };
   vertex[5]={  ooSq3, -ooSq3,  ooSq3 };
   vertex[6]={  ooSq3,  ooSq3,  ooSq3 };
   vertex[7]={ -ooSq3,  ooSq3,  ooSq3 };
   normal.resize(8);
   for ( size_t i=0 ; i<8 ; ++i ) { normal[i].resize(3); }
   sface.resize(6);
   for ( size_t i=0 ; i<6 ; ++i ) { sface[i].resize(3); }
   sface[0]={1,2,6,5};
   sface[1]={2,3,7,6};
   sface[2]={0,4,7,3};
   sface[3]={0,1,5,4};
   sface[4]={4,5,6,7};
   sface[5]={0,3,2,1};
   for ( int i=1 ; i<it ; ++i ) { SubdivideAllSquares(); }
   GenerateTrianglesFromSquares();
   CheckTriangleOrientation();
   shape=Shape::SPHERESQUARES;
}
void SymmetricSurfaceGrid::SetupSphereIcosahedron(const int it) {
   if ( vertex.size()!=0 ) {
      ClearArrays();
      shape=Shape::NONE;
   }
   static const double V0=0.0e0;
   static const double V5=2.0e0/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
   static const double V8=(1.0e0+sqrt(5.0e0))/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
   vertex.resize(12);
   vertex[ 0]={  V0, -V5,  V8 };
   vertex[ 1]={  V8,  V0,  V5 };
   vertex[ 2]={  V8,  V0, -V5 };
   vertex[ 3]={ -V8,  V0, -V5 };
   vertex[ 4]={ -V8,  V0,  V5 };
   vertex[ 5]={ -V5,  V8,  V0 };
   vertex[ 6]={  V5,  V8,  V0 };
   vertex[ 7]={  V5, -V8,  V0 };
   vertex[ 8]={ -V5, -V8,  V0 };
   vertex[ 9]={  V0, -V5, -V8 };
   vertex[10]={  V0,  V5, -V8 };
   vertex[11]={  V0,  V5,  V8 };
   normal.resize(12);
   for ( size_t i=0 ; i<12 ; ++i ) { normal[i].resize(3); }
   face.resize(20);
   for ( size_t i=0 ; i<20 ; ++i ) { face[i].resize(3); }
   face[ 0]={6,1,2};
   face[ 1]={2,1,7};
   face[ 2]={5,3,4};
   face[ 3]={8,4,3};
   face[ 4]={11,6,5};
   face[ 5]={10,5,6};
   face[ 6]={2,9,10};
   face[ 7]={3,10,9};
   face[ 8]={9,7,8};
   face[ 9]={0,8,7};
   face[10]={1,11,0};
   face[11]={4,0,11};
   face[12]={10,6,2};
   face[13]={11,1,6};
   face[14]={10,3,5};
   face[15]={11,5,4};
   face[16]={9,2,7};
   face[17]={0,7,1};
   face[18]={8,3,9};
   face[19]={0,4,8};
   for ( int i=1 ; i<it ; ++i ) { SubdivideAllTriangles(); }
   CheckTriangleOrientation();
   shape=Shape::SPHEREICOSAHEDRON;
}
void SymmetricSurfaceGrid::SubdivideAllSquares(int it) {
   for ( int j=0 ; j<it ; ++j ) {
      size_t ns=sface.size();
      for ( size_t i=0 ; i<ns ; ++i ) { SubdivideSquare(i); }
   }
   for ( size_t i=0 ; i<normal.size() ; ++i ) { normal[i].clear(); }
   normal.resize(vertex.size());
   for ( size_t i=0 ; i<normal.size() ; ++i ) { normal[i].resize(3); }
   GenerateTrianglesFromSquares();
   CheckTriangleOrientation();
}
void SymmetricSurfaceGrid::SubdivideAllTriangles(int it) {
   for ( int j=0 ; j<it ; ++j ) {
      size_t nt=face.size();
      for ( size_t i=0 ; i<nt ; ++i ) { SubdivideTriangle(i); }
   }
   for ( size_t i=0 ; i<normal.size() ; ++i ) { normal[i].clear(); }
   normal.resize(face.size());
   for ( size_t i=0 ; i<normal.size() ; ++i ) { normal[i].resize(3); }
}
void SymmetricSurfaceGrid::SubdivideTriangle(size_t idx) {
   vector<double> a(3),b(3),c(3),d(3),e(3),f(3);
   size_t iA,iB,iC,iD,iE,iF;
   iA=face[idx][0];
   iB=face[idx][1];
   iC=face[idx][2];
   a=vertex[iA];
   b=vertex[iB];
   c=vertex[iC];
   MatrixVectorOperations3D::Add(a,b,d);
   MatrixVectorOperations3D::Add(b,c,e);
   MatrixVectorOperations3D::Add(a,c,f);
   MatrixVectorOperations3D::Normalize(e);
   MatrixVectorOperations3D::Normalize(d);
   MatrixVectorOperations3D::Normalize(f);
   iD=VertexPosition(d);
   iE=VertexPosition(e);
   iF=VertexPosition(f);
   face[idx][0]=iA;
   face[idx][1]=iD;
   face[idx][2]=iF;
   vector<size_t> h(3);
   h[0]=iD; h[1]=iB; h[2]=iE;
   face.push_back(h);
   h[0]=iF; h[1]=iD; h[2]=iE;
   face.push_back(h);
   h[0]=iF; h[1]=iE; h[2]=iC;
   face.push_back(h);
}
void SymmetricSurfaceGrid::SubdivideSquare(size_t idx) {
   vector<double> A(3),B(3),C(3),D(3);
   vector<double> E(3),F(3),G(3),H(3);
   vector<double> I(3),J(3),K(3),L(3);
   vector<double> M(3),N(3),O(3),P(3);
   vector<double> tmp(3);
   double oo3=1.0e0/3.0e0;
   size_t iA,iB,iC,iD;
   size_t iE,iF,iG,iH;
   size_t iI,iJ,iK,iL;
   size_t iM,iN,iO,iP;
   iA=sface[idx][0]; iB=sface[idx][1]; iC=sface[idx][2]; iD=sface[idx][3];
   A=vertex[iA]; B=vertex[iB]; C=vertex[iC]; D=vertex[iD];

   MatrixVectorOperations3D::AminusB(B,A,tmp); MatrixVectorOperations3D::Scale(oo3,tmp);
   MatrixVectorOperations3D::Add(A,tmp,E);
   MatrixVectorOperations3D::Add(E,tmp,F);
   MatrixVectorOperations3D::AminusB(C,B,tmp); MatrixVectorOperations3D::Scale(oo3,tmp);
   MatrixVectorOperations3D::Add(B,tmp,G);
   MatrixVectorOperations3D::Add(G,tmp,H);
   MatrixVectorOperations3D::AminusB(D,C,tmp); MatrixVectorOperations3D::Scale(oo3,tmp);
   MatrixVectorOperations3D::Add(C,tmp,I);
   MatrixVectorOperations3D::Add(I,tmp,J);
   MatrixVectorOperations3D::AminusB(A,D,tmp); MatrixVectorOperations3D::Scale(oo3,tmp);
   MatrixVectorOperations3D::Add(D,tmp,K);
   MatrixVectorOperations3D::Add(K,tmp,L);
   MatrixVectorOperations3D::AminusB(J,E,tmp); MatrixVectorOperations3D::Scale(oo3,tmp);
   MatrixVectorOperations3D::Add(E,tmp,M);
   MatrixVectorOperations3D::Add(M,tmp,P);
   MatrixVectorOperations3D::AminusB(I,F,tmp); MatrixVectorOperations3D::Scale(oo3,tmp);
   MatrixVectorOperations3D::Add(F,tmp,N);
   MatrixVectorOperations3D::Add(N,tmp,O);

   MatrixVectorOperations3D::Normalize(E);
   MatrixVectorOperations3D::Normalize(F);
   MatrixVectorOperations3D::Normalize(G);
   MatrixVectorOperations3D::Normalize(H);
   MatrixVectorOperations3D::Normalize(I);
   MatrixVectorOperations3D::Normalize(J);
   MatrixVectorOperations3D::Normalize(K);
   MatrixVectorOperations3D::Normalize(L);
   MatrixVectorOperations3D::Normalize(M);
   MatrixVectorOperations3D::Normalize(N);
   MatrixVectorOperations3D::Normalize(O);
   MatrixVectorOperations3D::Normalize(P);
   iE=VertexPosition(E);
   iF=VertexPosition(F);
   iG=VertexPosition(G);
   iH=VertexPosition(H);
   iI=VertexPosition(I);
   iJ=VertexPosition(J);
   iK=VertexPosition(K);
   iL=VertexPosition(L);
   iM=VertexPosition(M);
   iN=VertexPosition(N);
   iO=VertexPosition(O);
   iP=VertexPosition(P);

   sface[idx][0]=iA;
   sface[idx][1]=iE;
   sface[idx][2]=iM;
   sface[idx][3]=iL;
   vector<size_t> h(4);
   h[0]=iL; h[1]=iM; h[2]=iP; h[3]=iK;
   sface.push_back(h);
   h[0]=iK; h[1]=iP; h[2]=iJ; h[3]=iD;
   sface.push_back(h);

   h[0]=iE; h[1]=iF; h[2]=iN; h[3]=iM;
   sface.push_back(h);
   h[0]=iM; h[1]=iN; h[2]=iO; h[3]=iP;
   sface.push_back(h);
   h[0]=iP; h[1]=iO; h[2]=iI; h[3]=iJ;
   sface.push_back(h);

   h[0]=iF; h[1]=iB; h[2]=iG; h[3]=iN;
   sface.push_back(h);
   h[0]=iN; h[1]=iG; h[2]=iH; h[3]=iO;
   sface.push_back(h);
   h[0]=iO; h[1]=iH; h[2]=iC; h[3]=iI;
   sface.push_back(h);
}
void SymmetricSurfaceGrid::CheckTriangleOrientation() {
   size_t nf=face.size();
   if ( nf==0 || nf==std::string::npos ) {
      ScreenUtils::DisplayWarningMessage("No triangles to check!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   vector<double> c(3),x1mx0(3),x2mx0(3),n(3);
   size_t k[3];
   double dot;
   for ( size_t i=0 ; i<nf ; ++i ) {
      c[0]=c[1]=c[2]=0.0e0;
      for ( size_t j=0 ; j<3 ; ++j ) {
         k[j]=face[i][j];
         c[0]+=vertex[k[j]][0];
         c[1]+=vertex[k[j]][1];
         c[2]+=vertex[k[j]][2];
      }
      x1mx0[0]=vertex[k[1]][0]-vertex[k[0]][0];
      x1mx0[1]=vertex[k[1]][1]-vertex[k[0]][1];
      x1mx0[2]=vertex[k[1]][2]-vertex[k[0]][2];
      x2mx0[0]=vertex[k[2]][0]-vertex[k[0]][0];
      x2mx0[1]=vertex[k[2]][1]-vertex[k[0]][1];
      x2mx0[2]=vertex[k[2]][2]-vertex[k[0]][2];
      n=MatrixVectorOperations3D::CrossProduct(x1mx0,x2mx0);
      dot=MatrixVectorOperations3D::InnerProduct(c,n);
      if ( dot<=0.0e0 ) {
         ScreenUtils::DisplayWarningMessage("Inwards normal!");
         cout << "tIdx: " << i << '\n';
         cout << __FILE__ << ", line: " << __LINE__ << '\n';
      }
   }
}
void SymmetricSurfaceGrid::GenerateTrianglesFromSquares() {
   for ( size_t i=0 ; i<face.size() ; ++i ) { face[i].clear(); }
   face.clear();
   for ( size_t i=0 ; i<sface.size() ; ++i ) {
      AddTrianglesFromSquare(i);
   }
}
void SymmetricSurfaceGrid::AddTrianglesFromSquare(size_t idx) {
   vector<double> a(3),b(3),c(3),d(3);
   size_t iA,iB,iC,iD;
   iA=sface[idx][0];
   iB=sface[idx][1];
   iC=sface[idx][2];
   iD=sface[idx][3];
   a=vertex[iA];
   b=vertex[iB];
   c=vertex[iC];
   d=vertex[iD];
   vector<size_t> h(3);
   h[0]=iA; h[1]=iB; h[2]=iC;
   face.push_back(h);
   h[0]=iA; h[1]=iC; h[2]=iD;
   face.push_back(h);
}
void SymmetricSurfaceGrid::TrimFacesCentroidDotProdBetweenVals(const vector<double> &d,const double val1,const double val2) {
   //cout << "Entering TrimFacesCentroidDotProdGreaterThanZero" << '\n';
   size_t count=0,nf=face.size();
   double dot,oodnorm,cnorm;
   oodnorm=1.0e0/MatrixVectorOperations3D::Norm(d);
   vector<double> c(3);
   vector<vector<size_t> > kf;
   ResizeMatrix(kf,nf,3);
   for ( size_t i=0 ; i<nf ; ++i ) {
      TriangleCentroidDir(i,c);
      cnorm=MatrixVectorOperations3D::Norm(c);
      dot=MatrixVectorOperations3D::InnerProduct(c,d)*oodnorm/cnorm;
      if ( dot>val1 && dot<=val2 ) {
         for ( size_t j=0 ; j<3 ; ++j ) { kf[count][j]=face[i][j]; }
         ++count;
      }
   }
   ResizeMatrix(face,count);
   for ( size_t i=0 ; i<count ; ++i ) {
      face[i][0]=kf[i][0];
      face[i][1]=kf[i][1];
      face[i][2]=kf[i][2];
   }
   ClearMatrix(kf);
   RemoveUnusedVertices();
}
void SymmetricSurfaceGrid::RemoveFacesUsingVertices(const vector<size_t> &v2r) {
   //cout << "Entering RemoveFacesUsingVertices" << '\n';
   size_t count=0,nf=face.size();
   bool rmme;
   vector<vector<size_t> > kf;
   ResizeMatrix(kf,nf,3);
   //cout << "Faces to be removed:" << '\n';
   for ( size_t i=0 ; i<nf ; ++i ) {
      rmme=false;
      for ( size_t j=0 ; j<v2r.size() ; ++j ) {
         for ( size_t k=0 ; k<3 ; ++k ) {
            if ( face[i][k]==v2r[j] ) { rmme=true; }
         }
         if ( rmme ) { break; }
      }
      if ( !rmme ) {
         for ( size_t k=0 ; k<3 ; ++k ) { kf[count][k]=face[i][k]; }
         ++count;
      } else {
         if ( verbose ) {
            cout << "Removing face " << i << ':' << '\n';
            for ( size_t k=0 ; k<3 ; ++k ) {
               cout << face[i][k] << ' ';
            }
            cout << '\n';
         }
      }
   }
   ResizeMatrix(face,count,3);
   for ( size_t i=0 ; i<count ; ++i ) {
      face[i][0]=kf[i][0];
      face[i][1]=kf[i][1];
      face[i][2]=kf[i][2];
   }
   ClearMatrix(kf);
   RemoveUnusedVertices();
   //cout << "Leaving RemoveFacesUsingVertices" << '\n';
}
void SymmetricSurfaceGrid::RemoveUnusedVertices() {
   //ScreenUtils::DisplayWarningMessage("SymmetricSurfaceGrid::RemoveUnusedVertices is not tested.");
   //cout << __FILE__ << ", line: " << __LINE__ << '\n';
   //size_t nvspf=0; //number of vertices per face.
   //size_t nf=0;//number of faces.
   //vector<vector<size_t> > *f=nullptr;
   //vector<vector<double> > *ctd=nullptr;
   if ( shape==Shape::SPHEREICOSAHEDRON ) {
      //f=&face;
      //ctd=&centroid;
      //nf=face.size();
      //nvspf=3;
      RemoveUnusedVerticesBase(face,size_t(3));
   } else if ( shape==Shape::SPHERESQUARES ) {
      //f=&sface;
      //nf=sface.size();
      //nvspf=4;
      cout << "This function may need adjustments." << '\n';
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      RemoveUnusedVerticesBase(sface,size_t(4));
   } else {
      ScreenUtils::DisplayWarningMessage("Non valid shape!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   /*
   cout << "Testing face selection" << '\n';
   cout << (*f)[0][0] << ' ' << (*f)[0][1] << ' ' << (*f)[0][2] << '\n';
   cout << face[0][0] << ' ' << face[0][1] << ' ' << face[0][2] << '\n';
   // */
   /*
   size_t pos=0,tmp;
   bool iamused;
   size_t nv=vertex.size();
   while ( pos<nv ) {
      iamused=false;
      for ( size_t i=0 ; i<nf ; ++i ) {
         for ( size_t j=0 ; j<nvspf ; ++j ) {
            if ( (*f)[i][j]==pos ) { iamused=true; }
         }
         if ( iamused ) { break; }
      }
      if ( !iamused ) {
         tmp=nv-1;
         for ( size_t i=0 ; i<3 ; ++i ) { vertex[pos][i]=vertex[tmp][i]; }
         for ( size_t i=0 ; i<nf ; ++i ) {
            for ( size_t j=0 ; j<nvspf ; ++j ) {
               if ( (*f)[i][j]==tmp ) { (*f)[i][j]=pos; }
            }
         }
         vertex.pop_back();
         --nv;
      } else {
         ++pos;
      }
   }
   //if ( ctd!=nullptr ) { ResizeMatrix((*ctd),(*f).size(),3); }
   ResizeMatrix(normal,vertex.size(),3);
   value.resize(vertex.size());
   // */
}
void SymmetricSurfaceGrid::DisplayFaces() {
   if ( shape==Shape::SPHEREICOSAHEDRON ) {
      for ( size_t i=0 ; i<face.size() ; ++i ) {
         for ( size_t j=0 ; j<3 ; ++j ) { cout << face[i][j] << ' '; }
         cout << '\n';
      }
   } else if ( shape==Shape::SPHERESQUARES ) {
      for ( size_t i=0 ; i<sface.size() ; ++i ) {
         for ( size_t j=0 ; j<4 ; ++j ) { cout << sface[i][j] << ' '; }
         cout << '\n';
      }
   }
}
void SymmetricSurfaceGrid::ComputeCentroids() {
   if ( shape!=Shape::SPHEREICOSAHEDRON ) {
      ScreenUtils::DisplayWarningMessage("ComputeCentroids is not implemented for this mesh.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   return ComputeCentroidsBase();
   /*
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
   // */
}

