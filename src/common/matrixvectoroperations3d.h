/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
   static vector<vector<double> > GetRotationMatrixAroundX(const double angle);
   static vector<vector<double> > GetRotationMatrixAroundY(const double angle);
   static vector<vector<double> > GetRotationMatrixAroundZ(const double angle);
   static vector<double> CrossProduct(const vector<double> &a,const vector<double> &b);
   static double InnerProduct(const vector<double> &a,const vector<double> &b);
   inline static double Norm(const vector<double> &a) {return sqrt(InnerProduct(a,a));}
   static double Distance2(const vector<double> &a,const vector<double> &b);
   inline static double Distance(const vector<double> &a,const vector<double> &b) {return sqrt(Distance2(a,b));}
   static void Normalize(vector<double> &a);
   static double Determinant(const vector<vector<double> > &M);
   static vector<vector<double> > UnitMatrix();
   static vector<vector<double> > Zeros();
   static vector<double> ZeroVector() {return vector<double>(3);}
   static void Add(const vector<double> &a,const vector<double> &b,vector<double> &c);
   static void AminusB(const vector<double> &a,const vector<double> &b,vector<double> &c);
   inline static void Scale(const double f,vector<double> &v) { v[0]*=f; v[1]*=f; v[2]*=f; };
   inline static void Scale(vector<double> &v,const double f) { Scale(f,v); };
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
   static void Transpose(vector<vector<double> > &m);
   static vector<vector<double> > Transpose(const vector<vector<double> > &m);
   /** The result matrix should be used to rotate the particles of the system. */
   static vector<vector<double> > GetRotationMatrix2AlignActive(const vector<double> &A,
         const vector<double> &B,const vector<double> &C);
   /** The result matrix should be used to rotate the frame reference (it also works
    * for rotating position and view/right/up vectors of povray cameras.  */
   static vector<vector<double> > GetRotationMatrix2AlignPassive(const vector<double> &A,
         const vector<double> &B,const vector<double> &C);
   static void TransformByMatrixMultiplication(const vector<vector<double> > &M,vector<double> &v);
   static vector<vector<double> > GetRotationMatrixAroundAxis(const vector<double> &omega,const double angle);
   /** Creates the unit vectors to form a Cartesian system, with the following definitions:
    * Y=(A-B); Z=(C-B)xY, X=YxZ. Here X, Y, and Z are NOT aligned with the standard x, y, and z, but
    * are aligned according to A, B, and C. The function DOES NOT check vector sizes of A, B, or C.  */
   static void GetCartesianSystemFrom3Vectors(const vector<double> &A,const vector<double> &B,const vector<double> &C,\
         vector<double> &xx,vector<double> &yy,vector<double> &zz);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _MATRIXVECTOROPERATIONS3D_H_ */

