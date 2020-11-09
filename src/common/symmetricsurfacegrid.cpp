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

SymmetricSurfaceGrid::SymmetricSurfaceGrid() {
   shape=Shape::NONE;
   vertex.clear();
   normal.clear();
   tcentroid.clear();
   tface.clear();
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
   ClearMatrix(tcentroid);
   ClearMatrix(tface);
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
   tface.resize(20);
   for ( size_t i=0 ; i<20 ; ++i ) { tface[i].resize(3); }
   tface[ 0]={6,1,2};
   tface[ 1]={2,1,7};
   tface[ 2]={5,3,4};
   tface[ 3]={8,4,3};
   tface[ 4]={11,6,5};
   tface[ 5]={10,5,6};
   tface[ 6]={2,9,10};
   tface[ 7]={3,10,9};
   tface[ 8]={9,7,8};
   tface[ 9]={0,8,7};
   tface[10]={1,11,0};
   tface[11]={4,0,11};
   tface[12]={10,6,2};
   tface[13]={11,1,6};
   tface[14]={10,3,5};
   tface[15]={11,5,4};
   tface[16]={9,2,7};
   tface[17]={0,7,1};
   tface[18]={8,3,9};
   tface[19]={0,4,8};
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
      size_t nt=tface.size();
      for ( size_t i=0 ; i<nt ; ++i ) { SubdivideTriangle(i); }
   }
   for ( size_t i=0 ; i<normal.size() ; ++i ) { normal[i].clear(); }
   normal.resize(tface.size());
   for ( size_t i=0 ; i<normal.size() ; ++i ) { normal[i].resize(3); }
}
size_t SymmetricSurfaceGrid::VertexPosition(const vector<double> &t) {
   for ( size_t i=0 ; i<vertex.size() ; ++i ) {
      if ( MatrixVectorOperations3D::Distance2(t,vertex[i]) <= SYMMSURFGRIDDEFCOORDSEPS2 ) {
         return i;
      }
   }
   vertex.push_back(t);
   return (vertex.size()-1);
}
void SymmetricSurfaceGrid::SubdivideTriangle(size_t idx) {
   vector<double> a(3),b(3),c(3),d(3),e(3),f(3);
   size_t iA,iB,iC,iD,iE,iF;
   iA=tface[idx][0];
   iB=tface[idx][1];
   iC=tface[idx][2];
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
   tface[idx][0]=iA;
   tface[idx][1]=iD;
   tface[idx][2]=iF;
   vector<size_t> h(3);
   h[0]=iD; h[1]=iB; h[2]=iE;
   tface.push_back(h);
   h[0]=iF; h[1]=iD; h[2]=iE;
   tface.push_back(h);
   h[0]=iF; h[1]=iE; h[2]=iC;
   tface.push_back(h);
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
   size_t nf=tface.size();
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
         k[j]=tface[i][j];
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

void SymmetricSurfaceGrid::Translate(const vector<double> &t) {
   size_t nv=vertex.size();
   for ( size_t i=0 ; i<nv ; ++i ) {
      vertex[i][0]+=t[0];
      vertex[i][1]+=t[1];
      vertex[i][2]+=t[2];
   }
   for ( size_t i=0 ; i<3 ; ++i ) { center[i]+=t[i]; }
}
void SymmetricSurfaceGrid::Scale(const double f) {
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
void SymmetricSurfaceGrid::ScaleSingleVertex(const size_t idx,const double f) {
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
void SymmetricSurfaceGrid::GenerateTrianglesFromSquares() {
   for ( size_t i=0 ; i<tface.size() ; ++i ) { tface[i].clear(); }
   tface.clear();
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
   tface.push_back(h);
   h[0]=iA; h[1]=iC; h[2]=iD;
   tface.push_back(h);
}
void SymmetricSurfaceGrid::TrimFacesCentroidDotProdGreaterThanZero(const vector<double> &d) {
   //cout << "Entering TrimFacesCentroidDotProdGreaterThanZero" << '\n';
   size_t count=0,nf=tface.size();
   double dot;
   vector<double> c(3);
   vector<vector<size_t> > kf;
   ResizeMatrix(kf,nf,3);
   for ( size_t i=0 ; i<nf ; ++i ) {
      TriangleCentroidDir(i,c);
      dot=MatrixVectorOperations3D::InnerProduct(c,d);
      if ( dot>0.0e0 ) {
         for ( size_t j=0 ; j<3 ; ++j ) { kf[count][j]=tface[i][j]; }
         ++count;
      }
   }
   ResizeMatrix(tface,count);
   for ( size_t i=0 ; i<count ; ++i ) {
      tface[i][0]=kf[i][0];
      tface[i][1]=kf[i][1];
      tface[i][2]=kf[i][2];
   }
   ClearMatrix(kf);
   RemoveUnusedVertices();
}
void SymmetricSurfaceGrid::RemoveFacesUsingVertices(const vector<size_t> &v2r) {
   //cout << "Entering RemoveFacesUsingVertices" << '\n';
   size_t count=0,nf=tface.size();
   bool rmme;
   vector<vector<size_t> > kf;
   ResizeMatrix(kf,nf,3);
   //cout << "Faces to be removed:" << '\n';
   for ( size_t i=0 ; i<nf ; ++i ) {
      rmme=false;
      for ( size_t j=0 ; j<v2r.size() ; ++j ) {
         for ( size_t k=0 ; k<3 ; ++k ) {
            if ( tface[i][k]==v2r[j] ) { rmme=true; }
         }
         if ( rmme ) { break; }
      }
      if ( !rmme ) {
         for ( size_t k=0 ; k<3 ; ++k ) { kf[count][k]=tface[i][k]; }
         ++count;
      } else {
         cout << "Removing face " << i << ':' << '\n';
         for ( size_t k=0 ; k<3 ; ++k ) {
            cout << tface[i][k] << ' ';
         }
         cout << '\n';
      }
   }
   ResizeMatrix(tface,count,3);
   for ( size_t i=0 ; i<count ; ++i ) {
      tface[i][0]=kf[i][0];
      tface[i][1]=kf[i][1];
      tface[i][2]=kf[i][2];
   }
   ClearMatrix(kf);
   RemoveUnusedVertices();
   //cout << "Leaving RemoveFacesUsingVertices" << '\n';
}
void SymmetricSurfaceGrid::ResizeMatrix(vector<vector<double> > &m,size_t nr,size_t nc) {
   size_t neor=m.size();
   for ( size_t i=0 ; i<neor ; ++i ) { m[i].clear(); }
   if ( neor!=nr ) { m.resize(nr); }
   for ( size_t i=0 ; i<nr ; ++i ) { m[i].resize(nc); }
}
void SymmetricSurfaceGrid::ResizeMatrix(vector<vector<size_t> > &m,size_t nr,size_t nc) {
   size_t neor=m.size();
   for ( size_t i=0 ; i<neor ; ++i ) { m[i].clear(); }
   if ( neor!=nr ) { m.resize(nr); }
   for ( size_t i=0 ; i<nr ; ++i ) { m[i].resize(nc); }
}
void SymmetricSurfaceGrid::TriangleCentroidDir(size_t tFIdx,vector<double> &c) {
   size_t p0,p1,p2;
   const double oo3=1.0e0/3.0e0;
   p0=tface[tFIdx][0]; p1=tface[tFIdx][1]; p2=tface[tFIdx][2];
   c[0] =(vertex[p0][0]); c[1] =(vertex[p0][1]); c[2] =(vertex[p0][2]);
   c[0]+=(vertex[p1][0]); c[1]+=(vertex[p1][1]); c[2]+=(vertex[p1][2]);
   c[0]+=(vertex[p2][0]); c[1]+=(vertex[p2][1]); c[2]+=(vertex[p2][2]);
   c[0]*=oo3;             c[1]*=oo3;             c[2]*=oo3;
   c[0]-=(center[0]);     c[1]-=(center[1]);     c[2]-=(center[2]);
}
/*
void SymmetricSurfaceGrid::EraseVertexAndAssociatedFaces(size_t idx) {
   size_t pos,pm1;
   size_t npos=0;
   --npos;
   bool iusethisvertex=false;
   if ( shape==Shape::SPHEREICOSAHEDRON ) {
      pos=tface.size()-1;
      pm1=tface.size()-1;
      while ( pos!=npos ) {
         iusethisvertex=false;
         for ( size_t j=0 ; j<3 ; ++j ) {
            if ( tface[pos][j]==idx ) { iusethisvertex=true; }
         }
         if ( !iusethisvertex ) {
            for ( size_t j=0 ; j<3 ; ++j ) { tface[pos][j]=tface[pm1][j]; }
            tface.pop_back();
            --pm1;
         }
         --pos;
      }
      pm1=vertex.size()-1;
      for ( size_t i=0 ; i<3 ; ++i ) { vertex[idx][i]=vertex[pm1][i]; }
      vertex.pop_back();
      for ( size_t i=0 ; i<tface.size() ; ++i ) {
         for ( size_t j=0 ; j<3 ; ++j ) {
            if ( tface[i][j]==pm1 ) { tface[i][j]=idx; }
         }
      }
   } else {
      ScreenUtils::DisplayWarningMessage("EraseVertexAndAssociatedFaces does not work for this shape!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
}
// */
void SymmetricSurfaceGrid::RemoveUnusedVertices() {
   //ScreenUtils::DisplayWarningMessage("SymmetricSurfaceGrid::RemoveUnusedVertices is not tested.");
   //cout << __FILE__ << ", line: " << __LINE__ << '\n';
   size_t nvspf=0; //number of vertices per face.
   size_t nf=0;//number of faces.
   vector<vector<size_t> > *f=nullptr;
   vector<vector<double> > *ctd=nullptr;
   if ( shape==Shape::SPHEREICOSAHEDRON ) {
      f=&tface;
      ctd=&tcentroid;
      nf=tface.size();
      nvspf=3;
   } else if ( shape==Shape::SPHERESQUARES ) {
      f=&sface;
      nf=sface.size();
      nvspf=4;
   } else {
      ScreenUtils::DisplayWarningMessage("Non valid shape!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   /*
   cout << "Testing face selection" << '\n';
   cout << (*f)[0][0] << ' ' << (*f)[0][1] << ' ' << (*f)[0][2] << '\n';
   cout << tface[0][0] << ' ' << tface[0][1] << ' ' << tface[0][2] << '\n';
   // */
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
}
void SymmetricSurfaceGrid::DisplayFaces() {
   if ( shape==Shape::SPHEREICOSAHEDRON ) {
      for ( size_t i=0 ; i<tface.size() ; ++i ) {
         for ( size_t j=0 ; j<3 ; ++j ) { cout << tface[i][j] << ' '; }
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
   size_t nf=tface.size(),p0,p1,p2;
   const double oo3=1.0e0/3.0e0;
   ResizeMatrix(tcentroid,nf,3);
   for ( size_t i=0 ; i<nf ; ++i ) {
      p0=tface[i][0]; p1=tface[i][1]; p2=tface[i][2];
      tcentroid[i][0] =(vertex[p0][0]); tcentroid[i][1] =(vertex[p0][1]); tcentroid[i][2] =(vertex[p0][2]);
      tcentroid[i][0]+=(vertex[p1][0]); tcentroid[i][1]+=(vertex[p1][1]); tcentroid[i][2]+=(vertex[p1][2]);
      tcentroid[i][0]+=(vertex[p2][0]); tcentroid[i][1]+=(vertex[p2][1]); tcentroid[i][2]+=(vertex[p2][2]);
      tcentroid[i][0]*=oo3;             tcentroid[i][1]*=oo3;             tcentroid[i][2]*=oo3;
   }
}

