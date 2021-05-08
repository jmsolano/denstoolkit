/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
/*
 Eigen-decomposition code for SYMMETRIC 2x2, 3x3, 4x4, 5x5, and 6x6 matrices.
 
 The code below is a modification from the code found in 
 
 http://barnesc.blogspot.ca/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
 
 In turn, the above modification has the following text:
 
 Eigen decomposition code for symmetric 3x3 matrices, copied from the public
 domain Java Matrix library JAMA.
 
 
 */

#ifndef _EIGENDECOMPOSITIONJAMA_H_
#define _EIGENDECOMPOSITIONJAMA_H_
#include <cmath>
#include <vector>
using std::vector;
#define N2 2
#define N3 3
#define N4 4

/* ************************************************************************** */
class EigenDecompositionJAMA {
/* ************************************************************************** */
/* Symmetric matrix A => eigenvectors in columns of V, corresponding
 eigenvalues in d. */
public:
/* ************************************************************************** */
   /** Evaluates eigenvalues (d) and eigenvectors (V) of the matrix A.
    * Eigenvectors are row-ordered, i.e. V[i] is the i-th eigenvector.
    * Eigenvalues are ordered from the least to the greatest. */
   static void EigenDecomposition2(vector<vector<double> > &A,vector<vector<double> > &V,vector<double> &d);
   static void EigenDecomposition3(vector<vector<double> > &A,vector<vector<double> > &V,vector<double> &d);
   static void EigenDecomposition4(vector<vector<double> > &A,vector<vector<double> > &V,vector<double> &d);
/* ************************************************************************** */
   /** Evaluates eigenvalues (d) and eigenvectors (V) of the matrix A.
    * Eigenvectors are column-ordered, i.e. 
    * the i-th eigenvector is {V[0][i], V[1][i], V[2][i]}.
    * Eigenvalues are ordered from the least to the greatest. */
   static void EigenDecomposition2(double (&A)[N2][N2], double (&V)[N2][N2], double (&d)[N2]);
   static void EigenDecomposition3(double (&A)[N3][N3], double (&V)[N3][N3], double (&d)[N3]);
   static void EigenDecomposition4(double (&A)[N4][N4], double (&V)[N4][N4], double (&d)[N4]);
   static inline double hypot2(double x, double y) { return sqrt(x*x+y*y); }
/* ************************************************************************** */
protected:
   static void tred22(double (&V)[N2][N2], double (&d)[N2], double (&e)[N2]);
   static void tql22(double (&V)[N2][N2], double (&d)[N2], double (&e)[N2]);
   static void tred23(double (&V)[N3][N3], double (&d)[N3], double (&e)[N3]);
   static void tql23(double (&V)[N3][N3], double (&d)[N3], double (&e)[N3]);
   static void tred24(double (&V)[N4][N4], double (&d)[N4], double (&e)[N4]);
   static void tql24(double (&V)[N4][N4], double (&d)[N4], double (&e)[N4]);
/* ************************************************************************** */
};

#endif//_EIGENDECOMPOSITIONJAMA_H_

