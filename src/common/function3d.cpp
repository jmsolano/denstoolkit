/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cmath>
#include "function3d.h"

Function3D::Function3D() {
   x1.resize(3);
   x2.resize(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      x1[i]=-1.75e0;
      x2[i]=1.24e0;
   }
}
double Function3D::testF1(const vector<double> &x) {
   double res1=0.0e0,res2=0.0e0;;
   for ( size_t i=0 ; i<x.size() ; ++i ) {
      res1+=((x[i]-x1[i])*(x[i]-x1[i]));
      res2+=((x[i]-x2[i])*(x[i]-x2[i]));
   }
   return (exp(-res1)+exp(-res2));
}
double Function3D::analytF1Res() {
   return 1.74685581658e0; //int _{3/2}^{3/2}testF1
}


