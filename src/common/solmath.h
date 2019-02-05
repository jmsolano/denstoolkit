/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.3.0
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */

/*
 solmath.h
 
 
 Created by Juan Manuel Solano Altamirano on 11/03/13.
 e-mail: jmsolanoalt@gmail.com
 This program was developed in The University of Guelph,
 50 Stone Road West, Guelph
 ON, N1G 2W1, Canada
 
 */

#ifndef _SOLMATH_H_
#define _SOLMATH_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef EPSTODIVIDE
#define EPSTODIVIDE (1.0e-60)
#endif

#ifndef _SOL_USE_SAFE_CHECKS_
#define _SOL_USE_SAFE_CHECKS_ 1
#endif

/* ************************************************************************** */
/* ************************************************************************** */
int factorial(const int n); //This function only works up to n=16.
/* ************************************************************************** */
int doubfact(const int n);  //This function only works up to n=19.
/* ************************************************************************** */
solreal magV3(solreal (&v)[3]);
/* ************************************************************************** */
void normalizeV3(solreal (&v)[3]);
/* ************************************************************************** */
void crossProductV3(solreal (&a)[3],solreal (&b)[3],solreal (&c)[3]);
/*
   This function calculates the cross product: c=axb
*/
/* ************************************************************************** */
solreal dotProductV3(solreal (&a)[3],solreal (&b)[3]);
/* ************************************************************************** */
solreal detM3x3(solreal (&a)[3][3]);
/* ************************************************************************** */
void invertM3x3(solreal (&oM)[3][3],solreal (&rM)[3][3]);
/* ************************************************************************** */
/* ************************************************************************** */
inline void sort2intmin2max(int &a,int&b)
{
   if (a<=b) {
      return;
   } else {
      int c=a;
      a=b;
      b=c;
   }
   return;
}
/* ************************************************************************** */
solreal rfactorial(const int n);
/* ************************************************************************** */
solreal BoysFunction(const int m,solreal x);
/* ************************************************************************** */
solreal RecursiveBoysFunction(const int m,const solreal x);
/* ************************************************************************** */
solreal TabBoysFunction(const int m,const solreal x);
/* ************************************************************************** */
void BoysFunction(const solreal x,const int nmax,solreal (&fn)[7]);
/* ************************************************************************** */
void BoysFunctionRec(const solreal x,const int nmax,solreal (&fn)[7]);
/* ************************************************************************** */
void BoysFunctionTab(const solreal x,const int nmax,solreal (&fn)[7]);
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
#endif//_SOLMATH_H_

