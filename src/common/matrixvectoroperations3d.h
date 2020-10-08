#ifndef _MATRIXVECTOROPERATIONS3D_H_
#define _MATRIXVECTOROPERATIONS3D_H_
#include <vector>
using std::vector;
#include <cmath>
#include <random>
#include <chrono>

#ifndef COORDSEPSILON
#define COORDSEPSILON 5.0e-01
#endif

/* ************************************************************************** */
class MatrixVectorOperations3D {
/* ************************************************************************** */
public:
   MatrixVectorOperations3D() {srand(long(std::chrono::high_resolution_clock::now().time_since_epoch().count()));}
   static vector<vector<double> > GetMatrixToAlignXToV(vector<double> &ix,vector<double> &iv);
   static vector<vector<double> > GetMatrixToAlignVToZ(vector<double> &v);
   static vector<double> CrossProduct(const vector<double> &a,const vector<double> &b);
   static double InnerProduct(const vector<double> &a,const vector<double> &b);
   inline static double Norm(const vector<double> &a) {return sqrt(InnerProduct(a,a));}
   static double Distance2(const vector<double> &a,const vector<double> &b);
   inline static double Distance(const vector<double> &a,const vector<double> &b) {return sqrt(Distance2(a,b));}
   static void Normalize(vector<double> &a);
   static double Determinant(const vector<vector<double> > &M);
   static vector<vector<double> > UnitMatrix();
   static vector<vector<double> > Zeros();
   /** Returns a unit vector, v, with v[idx]=1.  */
   static vector<double> UnitVector(size_t idx);
   static vector<double> X() {return UnitVector(0);}
   static vector<double> Y() {return UnitVector(1);}
   static vector<double> Z() {return UnitVector(2);}
   /** Returns a unit vector with random components.  */
   vector<double> RandomUnitVector();
   vector<double> RandomVector(const double mag);
   /** Returns a vector with random direction, and random magnitude (between 1 and 3).  */
   vector<double> RandomVector();
   static vector<vector<double> > MatrixProduct(const vector<vector<double> > &a,
         const vector<vector<double> > &b);
   static vector<double> MatrixVectorProduct(const vector<vector<double> > &m,vector<double> &v);
   vector<vector<double> > GetEulerRotationMatrix(const double alpha, const double beta,const double gamma);
   /** This will generate a rotation matrix with random Euler angles.  */
   vector<vector<double> > GetEulerRotationMatrix(void);
   static void RotateAroundXAxis(vector<double> &v,double angle=M_PI);
   static void RotateAroundYAxis(vector<double> &v,double angle=M_PI);
   static void RotateAroundZAxis(vector<double> &v,double angle=M_PI);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _MATRIXVECTOROPERATIONS3D_H_ */

