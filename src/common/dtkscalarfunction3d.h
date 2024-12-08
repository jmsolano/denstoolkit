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
#ifndef _DTKSCALARFUNCTION3D_H_
#define _DTKSCALARFUNCTION3D_H_
#include <vector>
using std::vector;
#include "../common/fldtypesdef.h"
#include "function3d.h"
#include "../common/gausswavefunction.h"

/* ************************************************************************** */
class DTKScalarFunction : public Function3D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   DTKScalarFunction();
   DTKScalarFunction(GaussWaveFunction &ugwf);
   void SetScalarFunction(char t);
   inline double f(const vector<double> &x) { return f(x[0],x[1],x[2]); }
   inline double f(const double x,const double y,const double z) {
      double trho;
      switch ( sft ) {
         case ScalarFieldType::DENS :
            return wf->EvalDensity(x,y,z);
            break;
         case ScalarFieldType::DENSM :
            return wf->EvalFTDensity(x,y,z);
            break;
         case ScalarFieldType::MGRD :
            return wf->EvalMagGradRho(x,y,z);
            break;
         case ScalarFieldType::LAPD :
            return wf->EvalLapRho(x,y,z);
            break;
         case ScalarFieldType::LOLD :
            return wf->EvalLOL(x,y,z);
            break;
         case ScalarFieldType::ELFD :
            return wf->EvalELF(x,y,z);
            break;
         case ScalarFieldType::SENT :
            return wf->EvalShannonEntropy(x,y,z);
            break;
         case ScalarFieldType::SENTM :
            return wf->EvalMomentumShannonEntropy(x,y,z);
            break;
         case ScalarFieldType::KEDK :
            return wf->EvalKineticEnergyK(x,y,z);
            break;
         case ScalarFieldType::KEDKM :
            return wf->EvalFTKineticEnergy(x,y,z);
            break;
         case ScalarFieldType::KEDG :
            return wf->EvalKineticEnergyG(x,y,z);
            break;
         case ScalarFieldType::MGLD :
            return wf->EvalMagGradLOL(x,y,z);
            break;
         case ScalarFieldType::MEPD :
            return wf->EvalMolElecPot(x,y,z);
            break;
         case ScalarFieldType::MLED :
            return wf->EvalMagLED(x,y,z);
            break;
         case ScalarFieldType::REDG :
            return wf->EvalReducedDensityGradient(x,y,z);
            break;
         case ScalarFieldType::ROSE :
            return wf->EvalRoSE(x,y,z);
            break;
         case ScalarFieldType::VPED :
            return wf->EvalVirialPotentialEnergyDensity(x,y,z);
            break;
         case ScalarFieldType::NCIS :
            return wf->EvalNCIs(x,y,z);
            break;
         case ScalarFieldType::NCIL :
            return wf->EvalNCILambda(x,y,z);
            break;
         case ScalarFieldType::EDFTA :
            trho=wf->EvalDensity(x,y,z);
            return (-0.75e0*pow((0.9549296585513720146e0*trho),1.0e0/3.0e0));
            break;
         case ScalarFieldType::ELLPY :
            return wf->EvalEllipticity(x,y,z);
            break;
         case ScalarFieldType::DORI :
            return wf->EvalDORI(x,y,z);
            break;
         case ScalarFieldType::RHO2 :
            return wf->EvalOneElecDisequilibrium(x,y,z);
            break;
         case ScalarFieldType::RHO2M:
            return wf->EvalFTOneElecDisequilibrium(x,y,z);
            break;
         case ScalarFieldType::SPND :
            return wf->EvalSpinDensity(x,y,z);
            break;
         case ScalarFieldType::SCFD :
            return wf->EvalCustomScalarField(x,y,z);
            break;
         default :
            if ( prntunknfld1st ) {
               printf("\033[31mUnknown field!\n");
               printf("%s, line: %d\033[0m\n",__FILE__,__LINE__);
               prntunknfld1st=false;
            }
            return 0.0e0;
            break;
      }
      return 0.0e0;
   }
   inline void gradf(double x,double y,double z,double (&gs)[3]) {
      double trho;
      switch ( sft ) {
         case ScalarFieldType::DENS :
            wf->EvalRhoGradRho(x,y,z,trho,gs);
            break;
         case ScalarFieldType::LOLD :
            double xx[3],h[3][3];
            xx[0]=x; xx[1]=y; xx[2]=z;
            wf->EvalHessLOL(xx,trho,gs,h);
            break;
         case ScalarFieldType::REDG :
         case ScalarFieldType::NCIS :
            wf->EvalGradReducedDensityGradient(x,y,z,gs);
            break;
         case ScalarFieldType::GLOL :
         case ScalarFieldType::LEDV :
         case ScalarFieldType::VCFD :
            if ( prntunknfld1st ) {
               printf("\033[31mThis gradient field cannot be computed!\n");
               printf("In the current version, only the gradients of Rho, ReducedGradientDensity,\n"
                     " and LOL are implemented...");
               printf("%s, line: %d\033[0m\n",__FILE__,__LINE__);
               prntunknfld1st=false;
            }
            gs[0]=gs[1]=gs[2]=0.0e0;
            break;
         default :
            NumericalGradient(x,y,z,gs);
            break;
      }
   }
   double fptr(double x,double y,double z) {return (wf->*selfnc)(x,y,z);}
   void NumericalGradient(double x,double y,double z,double (&g)[3]);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void SelectScalarFunctionPtr(const char prop);
   void SelectScalarFunctionPtr(const ScalarFieldType ft);
   ScalarFieldType sft;
   GaussWaveFunction* wf;
   static bool prntunknfld1st; /*!< Whether to print error messages.  */
   double (GaussWaveFunction::*selfnc)(double,double,double);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DTKSCALARFUNCTION3D_H_ */

