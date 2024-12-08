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
#ifndef _INTERPOLATORS_H_
#define _INTERPOLATORS_H_
#include <vector>
using std::vector;
#include "spline.h"

/* ************************************************************************** */
class Interpolator {
/* ************************************************************************** */
public:
   Interpolator();
/* ************************************************************************** */
   static double Linear(const double x1,const double y1,const double x2,const double y2,const double xi);
   static double Linear(const vector<double> &x,const vector<double> &y,const double xi);
   double Spline(const vector<double> &x,const vector<double> &y,const double xi);
   inline double Spline(const double xi) {return sp(xi);} //Unsafe, but faster (in principle).
   /** Returns true if v[0]<=vi<=v[n-1]. It saves the left (l) and right (r) positions that encloses the value vi in v. If
    * there is a v[k] element that exactly matches vi (i.e. vi==v[k]), then it sets
    * l=k, and r=k+1.  */
   static bool BinarySearch(const vector<double> &v,size_t &ul,size_t &ur,const double vi);
   void SetPoints(const vector<double> &x,const vector<double> &y);
   /** This will prepare the interpolator object to be ready for a new interpolation
    * set. Usefor for reusing the same Interpolator instance with different
    * interpolation vectors.  */
   void Reset() { spIsSetup=false; }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   tk::spline sp;
   bool spIsSetup;
   vector<double> xx,yy;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTERPOLATORS_H_ */

