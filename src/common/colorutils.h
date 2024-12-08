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
#ifndef _COLORUTILS_H_
#define _COLORUTILS_H_
#include <vector>
using std::vector;
#include <cmath>

/* ************************************************************************** */
class ColorUtils {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   /* Calculate hsl from rgb
      Hue is in degrees
      Lightness is between 0 and 1
      Saturation is between 0 and 1
      This code is a refactor of the code taken from:
      http://paulbourke.net/miscellaneous/colourspace/ */
   static void rgb2hsl(const double r,const double g,const double b,\
         double &h,double &s,double &l);
   static inline double rgbDouble(int v) {return (double(v)*oo255);}
   static inline int rgbInt(const double v) { return v<=1.0e0 ? round(255.0e0*v) : 1.0e0;} 
   /** Calculate rgb from hsl, reverse of rgb2hsl(...)
      Hue is in degrees
      Lightness is between 0 and 1
      Saturation is between 0 and 1
      This code is a refactor of the code taken from:
      http://paulbourke.net/miscellaneous/colourspace/ */
   static void hsl2rgb(const double h,const double s,const double l,\
         double &r,double &g,double &b);
/* ************************************************************************** */
   static constexpr double oo255=1.0e0/255.0e0;
   static constexpr double oo60=1.0e0/60.0e0;
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _COLORUTILS_H_ */

