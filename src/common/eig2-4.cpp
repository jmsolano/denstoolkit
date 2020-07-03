/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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
 
 Eigen-decomposition code for symmetric 3x3 and 4x4 matrices.
 
 The code below is a modification from the code found in 
 
 http://barnesc.blogspot.ca/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
 
 In turn, the above modification has the following text:
 
 Eigen decomposition code for symmetric 3x3 matrices, copied from the public
 domain Java Matrix library JAMA.
 
 The i-th eigenvector is {V[0][i], V[1][i], V[2][i]}
 
 */


#ifndef _EIG2_4_CPP_
#define _EIG2_4_CPP_

#include "eig2-4.h"
#include <math.h>

#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

#define N2 2
#define N3 3
#define N4 4


//*************************************************************************************************
static solreal hypot2(solreal x, solreal y) {
   return sqrt(x*x+y*y);
}
//*************************************************************************************************
// Symmetric Householder reduction to tridiagonal form.
//*************************************************************************************************
static void tred23(solreal V[N3][N3], solreal d[N3], solreal e[N3]) {
   
   //  This is derived from the Algol procedures tred2 by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
   for (int j = 0; j < N3; j++) {
      d[j] = V[N3-1][j];
   }
   
   // Householder reduction to tridiagonal form.
   
   for (int i = N3-1; i > 0; i--) {
      
      // Scale to avoid under/overflow.
      
      solreal scale = 0.0;
      solreal h = 0.0;
      for (int k = 0; k < i; k++) {
         scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
         e[i] = d[i-1];
         for (int j = 0; j < i; j++) {
            d[j] = V[i-1][j];
            V[i][j] = 0.0;
            V[j][i] = 0.0;
         }
      } else {
         
         // Generate Householder vector.
         
         for (int k = 0; k < i; k++) {
            d[k] /= scale;
            h += d[k] * d[k];
         }
         solreal f = d[i-1];
         solreal g = sqrt(h);
         if (f > 0) {
            g = -g;
         }
         e[i] = scale * g;
         h = h - f * g;
         d[i-1] = f - g;
         for (int j = 0; j < i; j++) {
            e[j] = 0.0;
         }
         
         // Apply similarity transformation to remaining columns.
         
         for (int j = 0; j < i; j++) {
            f = d[j];
            V[j][i] = f;
            g = e[j] + V[j][j] * f;
            for (int k = j+1; k <= i-1; k++) {
               g += V[k][j] * d[k];
               e[k] += V[k][j] * f;
            }
            e[j] = g;
         }
         f = 0.0;
         for (int j = 0; j < i; j++) {
            e[j] /= h;
            f += e[j] * d[j];
         }
         solreal hh = f / (h + h);
         for (int j = 0; j < i; j++) {
            e[j] -= hh * d[j];
         }
         for (int j = 0; j < i; j++) {
            f = d[j];
            g = e[j];
            for (int k = j; k <= i-1; k++) {
               V[k][j] -= (f * e[k] + g * d[k]);
            }
            d[j] = V[i-1][j];
            V[i][j] = 0.0;
         }
      }
      d[i] = h;
   }
   
   // Accumulate transformations.
   
   for (int i = 0; i < N3-1; i++) {
      V[N3-1][i] = V[i][i];
      V[i][i] = 1.0;
      solreal h = d[i+1];
      if (h != 0.0) {
         for (int k = 0; k <= i; k++) {
            d[k] = V[k][i+1] / h;
         }
         for (int j = 0; j <= i; j++) {
            solreal g = 0.0;
            for (int k = 0; k <= i; k++) {
               g += V[k][i+1] * V[k][j];
            }
            for (int k = 0; k <= i; k++) {
               V[k][j] -= g * d[k];
            }
         }
      }
      for (int k = 0; k <= i; k++) {
         V[k][i+1] = 0.0;
      }
   }
   for (int j = 0; j < N3; j++) {
      d[j] = V[N3-1][j];
      V[N3-1][j] = 0.0;
   }
   V[N3-1][N3-1] = 1.0;
   e[0] = 0.0;
} 

//*************************************************************************************************
// Symmetric tridiagonal QL algorithm.
//*************************************************************************************************

static void tql23(solreal V[N3][N3], solreal d[N3], solreal e[N3]) {
   
   //  This is derived from the Algol procedures tql2, by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
   for (int i = 1; i < N3; i++) {
      e[i-1] = e[i];
   }
   e[N3-1] = 0.0;
   
   solreal f = 0.0;
   solreal tst1 = 0.0;
   solreal eps = pow(2.0,-52.0);
   for (int l = 0; l < N3; l++) {
      
      // Find small subdiagonal element
      
      tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < N3) {
         if (fabs(e[m]) <= eps*tst1) {
            break;
         }
         m++;
      }
      
      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.
      
      if (m > l) {
         int iter = 0;
         do {
            iter = iter + 1;  // (Could check iteration count here.)
            
            // Compute implicit shift
            
            solreal g = d[l];
            solreal p = (d[l+1] - g) / (2.0 * e[l]);
            solreal r = hypot2(p,1.0);
            if (p < 0) {
               r = -r;
            }
            d[l] = e[l] / (p + r);
            d[l+1] = e[l] * (p + r);
            solreal dl1 = d[l+1];
            solreal h = g - d[l];
            for (int i = l+2; i < N3; i++) {
               d[i] -= h;
            }
            f = f + h;
            
            // Implicit QL transformation.
            
            p = d[m];
            solreal c = 1.0;
            solreal c2 = c;
            solreal c3 = c;
            solreal el1 = e[l+1];
            solreal s = 0.0;
            solreal s2 = 0.0;
            for (int i = m-1; i >= l; i--) {
               c3 = c2;
               c2 = c;
               s2 = s;
               g = c * e[i];
               h = c * p;
               r = hypot2(p,e[i]);
               e[i+1] = s * r;
               s = e[i] / r;
               c = p / r;
               p = c * d[i] - s * g;
               d[i+1] = h + s * (c * g + s * d[i]);
               
               // Accumulate transformation.
               
               for (int k = 0; k < N3; k++) {
                  h = V[k][i+1];
                  V[k][i+1] = s * V[k][i] + c * h;
                  V[k][i] = c * V[k][i] - s * h;
               }
            }
            p = -s * s2 * c3 * el1 * e[l] / dl1;
            e[l] = s * p;
            d[l] = c * p;
            
            // Check for convergence.
            
         } while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
   }
   
   // Sort eigenvalues and corresponding vectors.
   
   for (int i = 0; i < N3-1; i++) {
      int k = i;
      solreal p = d[i];
      for (int j = i+1; j < N3; j++) {
         if (d[j] < p) {
            k = j;
            p = d[j];
         }
      }
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (int j = 0; j < N3; j++) {
            p = V[j][i];
            V[j][i] = V[j][k];
            V[j][k] = p;
         }
      }
   }
}
//*************************************************************************************************
void eigen_decomposition3(solreal A[N3][N3], solreal V[N3][N3], solreal d[N3]) {
   solreal e[N3];
   for (int i = 0; i < N3; i++) {
      for (int j = 0; j < N3; j++) {
         V[i][j] = A[i][j];
      }
   }
   tred23(V, d, e);
   tql23(V, d, e);
}

//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
// Symmetric Householder reduction to tridiagonal form.
//*************************************************************************************************
static void tred24(solreal V[N4][N4], solreal d[N4], solreal e[N4]) {
   
   //  This is derived from the Algol procedures tred2 by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
   for (int j = 0; j < N4; j++) {
      d[j] = V[N4-1][j];
   }
   
   // Householder reduction to tridiagonal form.
   
   for (int i = N4-1; i > 0; i--) {
      
      // Scale to avoid under/overflow.
      
      solreal scale = 0.0;
      solreal h = 0.0;
      for (int k = 0; k < i; k++) {
         scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
         e[i] = d[i-1];
         for (int j = 0; j < i; j++) {
            d[j] = V[i-1][j];
            V[i][j] = 0.0;
            V[j][i] = 0.0;
         }
      } else {
         
         // Generate Householder vector.
         
         for (int k = 0; k < i; k++) {
            d[k] /= scale;
            h += d[k] * d[k];
         }
         solreal f = d[i-1];
         solreal g = sqrt(h);
         if (f > 0) {
            g = -g;
         }
         e[i] = scale * g;
         h = h - f * g;
         d[i-1] = f - g;
         for (int j = 0; j < i; j++) {
            e[j] = 0.0;
         }
         
         // Apply similarity transformation to remaining columns.
         
         for (int j = 0; j < i; j++) {
            f = d[j];
            V[j][i] = f;
            g = e[j] + V[j][j] * f;
            for (int k = j+1; k <= i-1; k++) {
               g += V[k][j] * d[k];
               e[k] += V[k][j] * f;
            }
            e[j] = g;
         }
         f = 0.0;
         for (int j = 0; j < i; j++) {
            e[j] /= h;
            f += e[j] * d[j];
         }
         solreal hh = f / (h + h);
         for (int j = 0; j < i; j++) {
            e[j] -= hh * d[j];
         }
         for (int j = 0; j < i; j++) {
            f = d[j];
            g = e[j];
            for (int k = j; k <= i-1; k++) {
               V[k][j] -= (f * e[k] + g * d[k]);
            }
            d[j] = V[i-1][j];
            V[i][j] = 0.0;
         }
      }
      d[i] = h;
   }
   
   // Accumulate transformations.
   
   for (int i = 0; i < N4-1; i++) {
      V[N4-1][i] = V[i][i];
      V[i][i] = 1.0;
      solreal h = d[i+1];
      if (h != 0.0) {
         for (int k = 0; k <= i; k++) {
            d[k] = V[k][i+1] / h;
         }
         for (int j = 0; j <= i; j++) {
            solreal g = 0.0;
            for (int k = 0; k <= i; k++) {
               g += V[k][i+1] * V[k][j];
            }
            for (int k = 0; k <= i; k++) {
               V[k][j] -= g * d[k];
            }
         }
      }
      for (int k = 0; k <= i; k++) {
         V[k][i+1] = 0.0;
      }
   }
   for (int j = 0; j < N4; j++) {
      d[j] = V[N4-1][j];
      V[N4-1][j] = 0.0;
   }
   V[N4-1][N4-1] = 1.0;
   e[0] = 0.0;
} 
//*************************************************************************************************
// Symmetric tridiagonal QL algorithm.
//*************************************************************************************************
static void tql24(solreal V[N4][N4], solreal d[N4], solreal e[N4]) {
   
   //  This is derived from the Algol procedures tql2, by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
   for (int i = 1; i < N4; i++) {
      e[i-1] = e[i];
   }
   e[N4-1] = 0.0;
   
   solreal f = 0.0;
   solreal tst1 = 0.0;
   solreal eps = pow(2.0,-52.0);
   for (int l = 0; l < N4; l++) {
      
      // Find small subdiagonal element
      
      tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < N4) {
         if (fabs(e[m]) <= eps*tst1) {
            break;
         }
         m++;
      }
      
      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.
      
      if (m > l) {
         int iter = 0;
         do {
            iter = iter + 1;  // (Could check iteration count here.)
            
            // Compute implicit shift
            
            solreal g = d[l];
            solreal p = (d[l+1] - g) / (2.0 * e[l]);
            solreal r = hypot2(p,1.0);
            if (p < 0) {
               r = -r;
            }
            d[l] = e[l] / (p + r);
            d[l+1] = e[l] * (p + r);
            solreal dl1 = d[l+1];
            solreal h = g - d[l];
            for (int i = l+2; i < N4; i++) {
               d[i] -= h;
            }
            f = f + h;
            
            // Implicit QL transformation.
            
            p = d[m];
            solreal c = 1.0;
            solreal c2 = c;
            solreal c3 = c;
            solreal el1 = e[l+1];
            solreal s = 0.0;
            solreal s2 = 0.0;
            for (int i = m-1; i >= l; i--) {
               c3 = c2;
               c2 = c;
               s2 = s;
               g = c * e[i];
               h = c * p;
               r = hypot2(p,e[i]);
               e[i+1] = s * r;
               s = e[i] / r;
               c = p / r;
               p = c * d[i] - s * g;
               d[i+1] = h + s * (c * g + s * d[i]);
               
               // Accumulate transformation.
               
               for (int k = 0; k < N4; k++) {
                  h = V[k][i+1];
                  V[k][i+1] = s * V[k][i] + c * h;
                  V[k][i] = c * V[k][i] - s * h;
               }
            }
            p = -s * s2 * c3 * el1 * e[l] / dl1;
            e[l] = s * p;
            d[l] = c * p;
            
            // Check for convergence.
            
         } while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
   }
   
   // Sort eigenvalues and corresponding vectors.
   
   for (int i = 0; i < N4-1; i++) {
      int k = i;
      solreal p = d[i];
      for (int j = i+1; j < N4; j++) {
         if (d[j] < p) {
            k = j;
            p = d[j];
         }
      }
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (int j = 0; j < N4; j++) {
            p = V[j][i];
            V[j][i] = V[j][k];
            V[j][k] = p;
         }
      }
   }
}
//*************************************************************************************************
void eigen_decomposition4(solreal A[N4][N4], solreal V[N4][N4], solreal d[N4]) {
   solreal e[N4];
   for (int i = 0; i < N4; i++) {
      for (int j = 0; j < N4; j++) {
         V[i][j] = A[i][j];
      }
   }
   tred24(V, d, e);
   tql24(V, d, e);
}
//*************************************************************************************************
//*************************************************************************************************
// Symmetric Householder reduction to tridiagonal form.
//*************************************************************************************************
static void tred22(solreal V[N2][N2], solreal d[N2], solreal e[N2]) {
   
   //  This is derived from the Algol procedures tred2 by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
   for (int j = 0; j < N2; j++) {
      d[j] = V[N2-1][j];
   }
   
   // Householder reduction to tridiagonal form.
   
   for (int i = N2-1; i > 0; i--) {
      
      // Scale to avoid under/overflow.
      
      solreal scale = 0.0;
      solreal h = 0.0;
      for (int k = 0; k < i; k++) {
         scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
         e[i] = d[i-1];
         for (int j = 0; j < i; j++) {
            d[j] = V[i-1][j];
            V[i][j] = 0.0;
            V[j][i] = 0.0;
         }
      } else {
         
         // Generate Householder vector.
         
         for (int k = 0; k < i; k++) {
            d[k] /= scale;
            h += d[k] * d[k];
         }
         solreal f = d[i-1];
         solreal g = sqrt(h);
         if (f > 0) {
            g = -g;
         }
         e[i] = scale * g;
         h = h - f * g;
         d[i-1] = f - g;
         for (int j = 0; j < i; j++) {
            e[j] = 0.0;
         }
         
         // Apply similarity transformation to remaining columns.
         
         for (int j = 0; j < i; j++) {
            f = d[j];
            V[j][i] = f;
            g = e[j] + V[j][j] * f;
            for (int k = j+1; k <= i-1; k++) {
               g += V[k][j] * d[k];
               e[k] += V[k][j] * f;
            }
            e[j] = g;
         }
         f = 0.0;
         for (int j = 0; j < i; j++) {
            e[j] /= h;
            f += e[j] * d[j];
         }
         solreal hh = f / (h + h);
         for (int j = 0; j < i; j++) {
            e[j] -= hh * d[j];
         }
         for (int j = 0; j < i; j++) {
            f = d[j];
            g = e[j];
            for (int k = j; k <= i-1; k++) {
               V[k][j] -= (f * e[k] + g * d[k]);
            }
            d[j] = V[i-1][j];
            V[i][j] = 0.0;
         }
      }
      d[i] = h;
   }
   
   // Accumulate transformations.
   
   for (int i = 0; i < N2-1; i++) {
      V[N2-1][i] = V[i][i];
      V[i][i] = 1.0;
      solreal h = d[i+1];
      if (h != 0.0) {
         for (int k = 0; k <= i; k++) {
            d[k] = V[k][i+1] / h;
         }
         for (int j = 0; j <= i; j++) {
            solreal g = 0.0;
            for (int k = 0; k <= i; k++) {
               g += V[k][i+1] * V[k][j];
            }
            for (int k = 0; k <= i; k++) {
               V[k][j] -= g * d[k];
            }
         }
      }
      for (int k = 0; k <= i; k++) {
         V[k][i+1] = 0.0;
      }
   }
   for (int j = 0; j < N2; j++) {
      d[j] = V[N2-1][j];
      V[N2-1][j] = 0.0;
   }
   V[N2-1][N2-1] = 1.0;
   e[0] = 0.0;
} 

//*************************************************************************************************
// Symmetric tridiagonal QL algorithm.
//*************************************************************************************************

static void tql22(solreal V[N2][N2], solreal d[N2], solreal e[N2]) {
   
   //  This is derived from the Algol procedures tql2, by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.
   
   for (int i = 1; i < N2; i++) {
      e[i-1] = e[i];
   }
   e[N2-1] = 0.0;
   
   solreal f = 0.0;
   solreal tst1 = 0.0;
   solreal eps = pow(2.0,-52.0);
   for (int l = 0; l < N2; l++) {
      
      // Find small subdiagonal element
      
      tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < N2) {
         if (fabs(e[m]) <= eps*tst1) {
            break;
         }
         m++;
      }
      
      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.
      
      if (m > l) {
         int iter = 0;
         do {
            iter = iter + 1;  // (Could check iteration count here.)
            
            // Compute implicit shift
            
            solreal g = d[l];
            solreal p = (d[l+1] - g) / (2.0 * e[l]);
            solreal r = hypot2(p,1.0);
            if (p < 0) {
               r = -r;
            }
            d[l] = e[l] / (p + r);
            d[l+1] = e[l] * (p + r);
            solreal dl1 = d[l+1];
            solreal h = g - d[l];
            for (int i = l+2; i < N2; i++) {
               d[i] -= h;
            }
            f = f + h;
            
            // Implicit QL transformation.
            
            p = d[m];
            solreal c = 1.0;
            solreal c2 = c;
            solreal c3 = c;
            solreal el1 = e[l+1];
            solreal s = 0.0;
            solreal s2 = 0.0;
            for (int i = m-1; i >= l; i--) {
               c3 = c2;
               c2 = c;
               s2 = s;
               g = c * e[i];
               h = c * p;
               r = hypot2(p,e[i]);
               e[i+1] = s * r;
               s = e[i] / r;
               c = p / r;
               p = c * d[i] - s * g;
               d[i+1] = h + s * (c * g + s * d[i]);
               
               // Accumulate transformation.
               
               for (int k = 0; k < N2; k++) {
                  h = V[k][i+1];
                  V[k][i+1] = s * V[k][i] + c * h;
                  V[k][i] = c * V[k][i] - s * h;
               }
            }
            p = -s * s2 * c3 * el1 * e[l] / dl1;
            e[l] = s * p;
            d[l] = c * p;
            
            // Check for convergence.
            
         } while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
   }
   
   // Sort eigenvalues and corresponding vectors.
   
   for (int i = 0; i < N2-1; i++) {
      int k = i;
      solreal p = d[i];
      for (int j = i+1; j < N2; j++) {
         if (d[j] < p) {
            k = j;
            p = d[j];
         }
      }
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (int j = 0; j < N2; j++) {
            p = V[j][i];
            V[j][i] = V[j][k];
            V[j][k] = p;
         }
      }
   }
}
//*************************************************************************************************
void eigen_decomposition2(solreal A[N2][N2], solreal V[N2][N2], solreal d[N2]) {
   solreal e[N2];
   for (int i = 0; i < N2; i++) {
      for (int j = 0; j < N2; j++) {
         V[i][j] = A[i][j];
      }
   }
   tred22(V, d, e);
   tql22(V, d, e);
}
//*************************************************************************************************
//*************************************************************************************************
#endif//_EIG2_4_CPP_


