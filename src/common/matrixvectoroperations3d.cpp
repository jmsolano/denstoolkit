#include <cstdlib>
#include <cassert>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cmath>
#include "matrixvectoroperations3d.h"


vector<vector<double> > MatrixVectorOperations3D::GetMatrixToAlignXToV(vector<double> &ix,vector<double> &iv) {
#if DEBUG
   if ( ix.size()!=3 || iv.size()!=3 ) {
      cerr << "Vectors have incorrect sizes!" << endl;
      cerr << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
#endif
   vector<double> z=iv;
   vector<double> vv=CrossProduct(ix,z);
   double v0=vv[0],v1=vv[1],v2=vv[2];
   double s=Norm(vv);
   double c=InnerProduct(z,ix);
   vector<vector<double> > vm=Zeros();
   vm[0][0]=0.0e0; vm[0][1]=-v2;   vm[0][2]=v1;
   vm[1][0]=v2;    vm[1][1]=0.0e0; vm[1][2]=-v0;
   vm[2][0]=-v1;   vm[2][1]=v0;    vm[2][2]=0.0e0;
   vector<vector<double> > vm2=MatrixProduct(vm,vm);
   vector<vector<double> > mm=UnitMatrix();
   double fact=(1.0e0-c)/(s*s);
   for ( int i=0 ; i<3 ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) {
         mm[i][j]+=(vm[i][j]+fact*vm2[i][j]);
      }
   }
   return mm;
}
vector<double> MatrixVectorOperations3D::CrossProduct(const vector<double> &a,
      const vector<double> &b) {
   vector<double> c(3);
   c[0]=a[1]*b[2]-a[2]*b[1];
   c[1]=a[2]*b[0]-a[0]*b[2];
   c[2]=a[0]*b[1]-a[1]*b[0];
   return c;
}
double MatrixVectorOperations3D::InnerProduct(const vector<double> &a,
      const vector<double> &b) {
#if DEBUG
   if ( a.size()!=3 || b.size()!=3 ) {
      cerr << "Vectors have incorrect sizes!" << endl;
      cerr << __FILE__ << ", line: " << __LINE__ << endl;
      return 0.0e0;
   }
#endif
   double res=0.0e0;
   for ( int i=0 ; i<3 ; ++i ) { res+=(a[i]*b[i]); }
   return res;
}
double MatrixVectorOperations3D::Distance2(const vector<double> &a,const vector<double> &b) {
   double sum=0.0e0;
   for ( size_t i=0 ; i<3 ; ++i ) { sum+=((a[i]-b[i])*(a[i]-b[i])); }
   return sum;
}
void MatrixVectorOperations3D::Normalize(vector<double> &a) {
   double oomag=1.0e0/Norm(a);
   for ( int i=0 ; i<3 ; ++i ) { a[i]*=oomag; }
}
double MatrixVectorOperations3D::Determinant(const vector<vector<double> > &M) {
   double res=M[0][0]*(M[1][1]*M[2][2]-M[2][1]*M[1][2]);
   res-=M[0][1]*(M[1][0]*M[2][2]-M[2][0]*M[1][2]);
   res+=M[0][2]*(M[1][0]*M[2][1]-M[2][0]*M[1][1]);
   return res;
}
vector<vector<double> > MatrixVectorOperations3D::UnitMatrix() {
   vector<vector<double> > m=Zeros();
   for ( int i=0 ; i<3 ; ++i ) { m[i][i]=1.0e0; }
   return m;
}
vector<vector<double> > MatrixVectorOperations3D::Zeros() {
   vector<vector<double> > m(3);
   for ( int i=0 ; i<3 ; ++i ) {
      m[i].resize(3);
      for ( int j=0 ; j<3 ; ++j ) { m[i][j]=0.0e0; }
   }
   return m;
}
vector<double> MatrixVectorOperations3D::UnitVector(size_t idx) {
   vector<double> res={0.0e0,0.0e0,0.0e0};
   res[idx]=1.0e0;
   return res;
}
vector<vector<double> > MatrixVectorOperations3D::MatrixProduct(const vector<vector<double> > &a,
         const vector<vector<double> > &b) {
   vector<vector<double> > c(3);
   for ( int i=0 ; i<3 ; ++i ) {
      c[i].resize(3);
      for ( int j=0 ; j<3 ; ++j ) { c[i][j]=0.0e0; }
   }
#if DEBUG
   if ( a.size()!=3 || b.size()!=3 ) {
      cerr << "Vectors have incorrect sizes!" << endl;
      cerr << __FILE__ << ", line: " << __LINE__ << endl;
      return c;
   }
#endif
   for ( int i=0 ; i<3 ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) { c[i][j]=0.0e0; }
   }
   for ( int i=0 ; i<3 ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) {
         for ( int k=0 ; k<3 ; ++k ) { c[i][j]+=(a[i][k]*b[k][j]); }
      }
   }
   return c;
}
vector<vector<double> > MatrixVectorOperations3D::GetMatrixToAlignVToZ(vector<double> &v) {
   vector<double> z={0.0e0,0.0e0,1.0e0};
   return GetMatrixToAlignXToV(v,z);
}
vector<double> MatrixVectorOperations3D::MatrixVectorProduct(
      const vector<vector<double> > &m,vector<double> &v) {
   vector<double> res=v;
   size_t dim=v.size();
   if ( dim!=m[0].size() ) {
      cerr << "Warning: very likely, the matrix.vector cannot be performed." << endl
           << "   The sizes of the vector and matrix do not match!" << endl;
      cerr << __FILE__ << ", line: " << __LINE__ << endl;
   }
   for ( size_t i=0 ; i<dim ; ++i ) {
      res[i]=0.0e0;
      for ( size_t j=0 ; j<dim ; ++j ) {
         res[i]+=(m[i][j]*v[j]);
      }
   }
   return res;
}
/* Check https://en.wikipedia.org/wiki/Rotation_matrix  */
void MatrixVectorOperations3D::RotateAroundXAxis(vector<double> &v,double angle) {
   double ct=cos(angle),st=sin(angle);
   vector<vector<double> > m=Zeros();
   m[0][0]=1.0e0;
   m[1][1]=ct; m[1][2]=-st;
   m[2][1]=st; m[2][2]=ct;
   vector<double> tmp=v;
   v=MatrixVectorProduct(m,tmp);
}
/* Check https://en.wikipedia.org/wiki/Rotation_matrix  */
void MatrixVectorOperations3D::RotateAroundYAxis(vector<double> &v,double angle) {
   double ct=cos(angle),st=sin(angle);
   vector<vector<double> > m=Zeros();
   m[0][0]=ct; m[0][2]=st;
   m[1][1]=1.0e0;
   m[2][0]=-st; m[2][2]=ct;
   vector<double> tmp=v;
   v=MatrixVectorProduct(m,tmp);
}
/* Check https://en.wikipedia.org/wiki/Rotation_matrix  */
void MatrixVectorOperations3D::RotateAroundZAxis(vector<double> &v,double angle) {
   double ct=cos(angle),st=sin(angle);
   vector<vector<double> > m=Zeros();
   m[0][0]=ct; m[0][1]=-st;
   m[1][0]=st; m[1][1]=ct;
   m[2][2]=1.0e0;
   vector<double> tmp=v;
   v=MatrixVectorProduct(m,tmp);
}
/* In this library, we use the convention Z1X2Z3, which matches the convention
 * used in Goldstein's Classical Mechanics book (2nd Ed.). See also
 * https://en.wikipedia.org/wiki/Euler_angles
 * section "Conversion to other orientation representations"/"Rotation matrix".
 * Look for the table and see column "Proper Euler angles"  */
vector<vector<double> > MatrixVectorOperations3D::GetEulerRotationMatrix(const double alpha,\
      const double beta,const double gamma) {
   double c1=cos(alpha), c2=cos(beta), c3=cos(gamma);
   double s1=sin(alpha), s2=sin(beta), s3=sin(gamma);
   vector<vector<double> > R=Zeros();
   // --------------------------------------------------------------------
   R[0][0]= c1*c3-c2*s1*s3; R[0][1]=-c1*s3-c2*c3*s1; R[0][2]= s1*s2;
   R[1][0]= c3*s1+c1*c2*s3; R[1][1]= c1*c2*c3-s1*s3; R[1][2]=-c1*s2;
   R[2][0]= s2*s3;          R[2][1]= c3*s2;          R[2][2]=    c2;
   // --------------------------------------------------------------------
   return R;
}
vector<vector<double> > MatrixVectorOperations3D::GetEulerRotationMatrix() {
   double alpha=2.0e0*M_PI*(double(rand())/double(RAND_MAX));
   double beta=2.0e0*M_PI*(double(rand())/double(RAND_MAX));
   double gamma=2.0e0*M_PI*(double(rand())/double(RAND_MAX));
   //cout << "alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << endl;
   return GetEulerRotationMatrix(alpha,beta,gamma);
}
vector<double> MatrixVectorOperations3D::RandomUnitVector() {
   vector<double> res=UnitVector(0);
   for ( size_t i=0 ; i<3 ; ++i ) { res[i]=-1.0e0+2.0e0*(double(rand())/double(RAND_MAX)); }
   Normalize(res);
   return res;
}
vector<double> MatrixVectorOperations3D::RandomVector(const double mag) {
   vector<double> res=RandomUnitVector();
   for ( size_t i=0 ; i<3 ; ++i ) { res[i]*=mag; }
   return res;
}
vector<double> MatrixVectorOperations3D::RandomVector() {
   double mag=1.0e0+2.0e0*(double(rand())/double(RAND_MAX));
   return RandomVector(mag);
}
