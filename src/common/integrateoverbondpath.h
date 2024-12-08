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

#ifndef _INTEGRATEOVERBONDPATH_H_
#define _INTEGRATEOVERBONDPATH_H_
#include "fldtypesdef.h"
#include <fstream>
using std::ofstream;

/* ************************************************************************** */
class IntegrateOverBondPath {
/* ************************************************************************** */
public:
   IntegrateOverBondPath(class GaussWaveFunction &ugwf,class CritPtNetWork &ucpn,\
         ScalarFieldType utp=DENS);
   ~IntegrateOverBondPath();
   void ComputeAllIntegralsOverBondPaths(void);
   double GetBondPathIntegral();
   void WriteIntegralValuesToFile(ofstream &ofil,double globalEnergy=0.0e0);
   inline double GetBondPathIntegral(int idx) {return integralValue[idx];}
   inline string GetFieldTypeLabelShort() {return GetFieldTypeKeyShort(myCharFieldType);}
   inline string GetFieldTypeLabelLong() {return GetFieldTypeKeyLong(myCharFieldType);}
/* ************************************************************************** */
protected:
   class GaussWaveFunction *wf;
   class CritPtNetWork *cp;
   /* Functions  */
   void Init(void);
   void GetIntermediateCoordinatesAndDistanceBetweenPoints(int startIdx,double** (&arr),double &h,\
         double (&x0)[3],double (&x1)[3],double (&x2)[3],double (&x3)[3]);
   void ComputeScalarFunctionValuesAtIntermediatePoints(double (&x0)[3],double (&x1)[3],\
         double (&x2)[3],double (&x3)[3],double (&f)[4]);
   double ComputeBondPathIntegral(int bpIdx);
   void ComputeBondPathIntegrals(void);
   bool AllocateAuxiliaryArrays(void);
   double TheFunction(double (&x)[3]);
/* ************************************************************************** */
   static constexpr double oo3=1.0e0/3.0e0,threeo8=3.0e0/8.0e0;
   static constexpr double threeo4=3.0e0/4.0e0;
   static constexpr double threeopi=3.0e0/3.14159265358979323846264e0;
   int nbgp;
   ScalarFieldType myFieldType;
   char myCharFieldType;
   double *integralValue;
/* ************************************************************************** */
   IntegrateOverBondPath(); //Default constructor is not allowed!
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATEOVERBONDPATH_H_ */

