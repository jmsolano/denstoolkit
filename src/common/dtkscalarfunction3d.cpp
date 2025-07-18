/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
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
#include "dtkscalarfunction3d.h"
#include "screenutils.h"
#include "fldtypesdef.h"

bool DTKScalarFunction::prntunknfld1st=true;
DTKScalarFunction::DTKScalarFunction() : Function3D() {
   wf=nullptr;
   sft=ScalarFieldType::NONE;
   selfnc=nullptr;
}
DTKScalarFunction::DTKScalarFunction(GaussWaveFunction &ugwf) : DTKScalarFunction() {
   wf=&ugwf;
   sft=ScalarFieldType::DENS;
}
void DTKScalarFunction::SetScalarFunction(char t) {
   sft=Char2ScalarFieldType(t);
   SelectScalarFunctionPtr(t);
}
void DTKScalarFunction::SelectScalarFunctionPtr(const char prop) {
   return SelectScalarFunctionPtr(Char2ScalarFieldType(prop));
}
void DTKScalarFunction::SelectScalarFunctionPtr(const ScalarFieldType ft) {
   switch ( ft ) {
      case DENS :
         selfnc=&GaussWaveFunction::EvalDensity;
         break;
      case DENSM :
         selfnc=&GaussWaveFunction::EvalFTDensity;
         break;
      case MGRD :
         selfnc=&GaussWaveFunction::EvalMagGradRho;
         break;
      case LAPD :
         selfnc=&GaussWaveFunction::EvalLapRho;
         break;
      case ELFD :
         selfnc=&GaussWaveFunction::EvalELF;
         break;
      case SENT :
         selfnc=&GaussWaveFunction::EvalShannonEntropy;
         break;
      case SENTM :
         selfnc=&GaussWaveFunction::EvalMomentumShannonEntropy;
         break;
      case KEDK :
         selfnc=&GaussWaveFunction::EvalKineticEnergyK;
         break;
      case KEDKM :
         selfnc=&GaussWaveFunction::EvalFTKineticEnergy;
         break;
      case KEDG :
         selfnc=&GaussWaveFunction::EvalKineticEnergyG;
         break;
      case MGLD :
         selfnc=&GaussWaveFunction::EvalMagGradLOL;
         break;
      case MEPD :
         selfnc=&GaussWaveFunction::EvalMolElecPot;
         break;
      case MLED :
         selfnc=&GaussWaveFunction::EvalMagLED;
         break;
      case REDG :
         selfnc=&GaussWaveFunction::EvalReducedDensityGradient;
         break;
      case ROSE :
         selfnc=&GaussWaveFunction::EvalRoSE;
         break;
      case VPED :
         selfnc=&GaussWaveFunction::EvalVirialPotentialEnergyDensity;
         break;
      case NCIS :
         selfnc=&GaussWaveFunction::EvalNCIs;
         break;
      case NCIL :
         selfnc=&GaussWaveFunction::EvalNCILambda;
         break;
      case ELLPY :
         selfnc=&GaussWaveFunction::EvalEllipticity;
         break;
      case DORI :
         selfnc=&GaussWaveFunction::EvalDORI;
         break;
      case SPND :
         selfnc=&GaussWaveFunction::EvalSpinDensity;
         break;
      case RHO2 :
         selfnc=&GaussWaveFunction::EvalOneElecDisequilibrium;
         break;
      case RHO2M:
         selfnc=&GaussWaveFunction::EvalFTOneElecDisequilibrium;
         break;
      case SCFD :
         selfnc=&GaussWaveFunction::EvalCustomScalarField;
         break;
      default :
         ScreenUtils::DisplayErrorMessage("Field not setup to be choosable!\n"
               "Setting up electron deensity as the function.");
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         selfnc=&GaussWaveFunction::EvalDensity;
         break;
   }
}
void DTKScalarFunction::NumericalGradient(double x,double y,double z,double (&g)[3]) {
   double hh=0.02;
   double oo2h=0.5e0/hh;
   double ffpp=f(x+hh,y,z)-f(x-hh,y,z);
   g[0]=oo2h*ffpp;
   ffpp=f(x,y+hh,z)-f(x,y-hh,z);
   g[1]=oo2h*ffpp;
   ffpp=f(x,y,z+hh)-f(x,y,z-hh);
   g[2]=oo2h*ffpp;
}
