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
 solmath.h
 
 
 Created by Juan Manuel Solano Altamirano on 11/03/13.
 e-mail: jmsolanoalt@gmail.com
 This program was developed in The University of Guelph,
 50 Stone Road West, Guelph
 ON, N1G 2W1, Canada
 
 */

#ifndef _MYMATH_H_
#define _MYMATH_H_
#ifndef EPSTODIVIDE
#define EPSTODIVIDE (1.0e-60)
#endif
/* ************************************************************************** */
/* ************************************************************************** */
int factorial(const int n); //This function only works up to n=16.
/* ************************************************************************** */
int doubfact(const int n);  //This function only works up to n=19.
/* ************************************************************************** */
double magV3(double (&v)[3]);
/* ************************************************************************** */
void normalizeV3(double (&v)[3]);
/* ************************************************************************** */
void crossProductV3(double (&a)[3],double (&b)[3],double (&c)[3]);
/*
   This function calculates the cross product: c=axb
*/
/* ************************************************************************** */
double dotProductV3(double (&a)[3],double (&b)[3]);
/* ************************************************************************** */
double detM3x3(double (&a)[3][3]);
/* ************************************************************************** */
void invertM3x3(double (&oM)[3][3],double (&rM)[3][3]);
/* ************************************************************************** */
/* ************************************************************************** */
inline void sort2intmin2max(int &a,int&b) {
   if (a<=b) {return; }
   else {
      int c=a;
      a=b;
      b=c;
   }
   return;
}
/* ************************************************************************** */
double rfactorial(const int n);
/* ************************************************************************** */
double BoysFunction(const int m,double x);
/* ************************************************************************** */
double RecursiveBoysFunction(const int m,const double x);
/* ************************************************************************** */
double TabBoysFunction(const int m,const double x);
/* ************************************************************************** */
void BoysFunction(const double x,const int nmax,double (&fn)[7]);
/* ************************************************************************** */
void BoysFunctionRec(const double x,const int nmax,double (&fn)[7]);
/* ************************************************************************** */
void BoysFunctionTab(const double x,const int nmax,double (&fn)[7]);
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
#endif//_MYMATH_H_

