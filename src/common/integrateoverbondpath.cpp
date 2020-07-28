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

#ifndef _INTEGRATEOVERBONDPATH_CPP_
#define _INTEGRATEOVERBONDPATH_CPP_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "gausswavefunction.h"
#include "critptnetwork.h"
#include "integrateoverbondpath.h"
#include "screenutils.h"
#include "mymemory.h"
#include "mymath.h"
#include "fileutils.h"

void IntegrateOverBondPath::Init(void) {
   wf=NULL;
   cp=NULL;
   integralValue=NULL;
   nbgp=0;
   myFieldType=NONE;
}
IntegrateOverBondPath::IntegrateOverBondPath(GaussWaveFunction &ugwf,critPtNetWork &ucpn,\
      ScalarFieldType utp) {
   Init();
   wf=&ugwf;
   if ( !ucpn.iKnowBCPs() ) {
      ScreenUtils::DisplayErrorMessage("First seek the critical points!");
      ScreenUtils::DisplayWarningMessage("critPtNetWork pointer is set to null!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   if ( !ucpn.iKnowBGPs() ) {
      ScreenUtils::DisplayErrorMessage("First compute the bond paths!");
      ScreenUtils::DisplayWarningMessage("critPtNetWork pointer is set to null!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   cp=&ucpn;
   nbgp=cp->nBGP;
   if ( !AllocateAuxiliaryArrays() ) {
      ScreenUtils::DisplayErrorMessage("AuxiliaryArrys could not be allocated!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
   myFieldType=utp;
   myCharFieldType=ConvertScalarFieldType2Char(myFieldType);
}
bool IntegrateOverBondPath::AllocateAuxiliaryArrays(void) {
   return MyMemory::Alloc1DRealArray("integralValue",nbgp,integralValue);
}
IntegrateOverBondPath::~IntegrateOverBondPath() {
   MyMemory::Dealloc1DRealArray(integralValue);
   wf=NULL;
   cp=NULL;
}
double IntegrateOverBondPath::ComputeBondPathIntegral(int bpIdx) {
   double xi[3],xip1[3],xip2[3],xip3[3],h,f[4];
   int nReg=cp->conBCP[bpIdx][2]-1;
   double res=0.0e0;
   for ( int i=0 ; i<nReg ; ++i ) {
      GetIntermediateCoordinatesAndDistanceBetweenPoints(i,cp->RBGP[bpIdx],h,\
            xi,xip1,xip2,xip3);
      ComputeScalarFunctionValuesAtIntermediatePoints(xi,xip1,xip2,xip3,f);
      res+=(threeo8*h*(f[0]+(3.0e0*f[1])+(3.0e0*f[2])+f[3]));
   }
   return res;
}
void IntegrateOverBondPath::GetIntermediateCoordinatesAndDistanceBetweenPoints(int startIdx,\
      double** (&arr),double &h,double (&x0)[3],double (&x1)[3],\
      double (&x2)[3],double (&x3)[3]) {
   double x3mx0[3];
   for ( int i=0 ; i<3 ; ++i ) {
      x0[i]=arr[startIdx][i];
      x3[i]=arr[startIdx+1][i];
      x3mx0[i]=x3[i]-x0[i];
   }
   h=magV3(x3mx0)*oo3;
   for ( int i=0 ; i<3 ; ++i ) {
      x1[i]=x0[i]+oo3*x3mx0[i];
      x2[i]=x1[i]+oo3*x3mx0[i];
   }
}
void IntegrateOverBondPath::ComputeScalarFunctionValuesAtIntermediatePoints(\
      double (&x0)[3],double (&x1)[3],double (&x2)[3],double (&x3)[3],\
      double (&f)[4]) {
   f[0]=TheFunction(x0);
   f[1]=TheFunction(x1);
   f[2]=TheFunction(x2);
   f[3]=TheFunction(x3);
}
void IntegrateOverBondPath::ComputeBondPathIntegrals(void) {
   for ( int i=0 ; i<nbgp ; ++i ) {
      integralValue[i]=ComputeBondPathIntegral(i);
   }
}
double IntegrateOverBondPath::GetBondPathIntegral(void) {
   double res=0.0e0;
   for ( int i=0 ; i<nbgp ; ++i ) { res+=integralValue[i]; }
   return res;
}
void IntegrateOverBondPath::ComputeAllIntegralsOverBondPaths(void) {
   ComputeBondPathIntegrals();
}
double IntegrateOverBondPath::TheFunction(double (&x)[3]) {
   double trho;
   switch ( myFieldType ) {
      case DENS :
         return wf->EvalDensity(x[0],x[1],x[2]);
         break;
      case KEDK :
         return wf->EvalKineticEnergyK(x[0],x[1],x[2]);
         break;
      case KEDG :
         return wf->EvalKineticEnergyG(x[0],x[1],x[2]);
         break;
      case REDG :
         return wf->EvalReducedDensityGradient(x[0],x[1],x[2]);
         break;
      case MEPD :
         return wf->EvalMolElecPot(x[0],x[1],x[2]);
         break;
      case EDFTA :
         trho=wf->EvalDensity(x[0],x[1],x[2]);
         return (-threeo4*pow((threeopi*trho),1.0e0/3.0e0));
         break;
      default :
         break;
   }
   return 0.0e0;
}
void IntegrateOverBondPath::WriteIntegralValuesToFile(ofstream &ofil,double globalEnergy) {
   bool haveGE=false;
   double factor=1.0e0;
   if ( globalEnergy!=0.0e0 ) {
      haveGE=true;
      factor=globalEnergy/GetBondPathIntegral();
   }
   FileUtils::WriteScrStarLine(ofil);
   FileUtils::WriteCenteredString(ofil,GetFieldTypeLabelLong());
   FileUtils::WriteScrStarLine(ofil);
   ofil << "Total integral: " << GetBondPathIntegral();
   if ( haveGE ) {
      ofil << " " << globalEnergy << endl;
   }
   ofil << endl;
   for ( int i=0 ; i<nbgp ; ++i ) {
      ofil << cp->lblBCP[i] << " bond path integral: " << GetBondPathIntegral(i);
      if ( haveGE ) {
         ofil << " " << (GetBondPathIntegral(i)*factor);
      }
      ofil << endl;
   }
}
#endif  /* _INTEGRATEOVERBONDPATH_CPP_ */

