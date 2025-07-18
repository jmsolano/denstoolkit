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
#include "gausswavefunction.h"
/* ************************************************************************************ */
/* In this file, double is a mask for a double, i.e., a double is a double.
 * This name is used in order to keep compatibility with future implementations
 * requiring floats rather than doubles (old GPUs).  */
double GaussWaveFunction::EvalCustomScalarField(double x,double y,double z) {
   /* 
    * double rho=EvalDensity(x,y,z) //Electron density [ED]
    * double maggrho=EvalMagGradRho(x,y,z); //Magnitude of the Gradient of ED
    * double lap=EvalLapRho(x,y,z); //Laplacian of ED
    * double lol=EvalLOL(x,y,z); //Localized Orbital Locator
    * double elf=EvalELF(x,y,z); //Returns the Electron Localized Function
    * double shent=EvalShannonEntropy(x,y,z); //Shannon entropy density
    * double mssent=EvalMomentumShannonEntropy(px,py,pz); //momentum space shannon
    *                                     //entropy at the momentum-point (px,py,pz)
    * double keG=EvalKineticEnergyG(x,y,z); //self descriptive
    * double keK=EvalKineticEnergyK(x,y,z); //self descriptive
    * double ftrho=EvalFTDensity(px,py,pz); //the momentum-space electron density
    *                                        // at the momentum-point (px,py,pz)
    * double mglol=EvalMagGradLOL(x,y,z); // Magnitude of grad(LOL)
    * double mep=EvalMolElecPot(x,y,z); //Molecular Electrostatic Potential
    * double magled=EvalMagLED(x,y,z); //Magnitude of LED
    * double rose=EvalRoSE(x,y,z); //Region of Slow Electrons
    * double s=EvalReducedDensityGradient(x,y,z); //self descriptive
    *
    * */

   /* What follows is an example of how to implement the field rho^2
    *
    * for other fields, you can choose one or more of the above enlisted fields.
    * */
   double eve[3][3],eva[3],h[3][3];
   EvalHessian(x,y,z,h);
   EigenDecompositionJAMA::EigenDecomposition3(h,eve,eva);
   double ellip=(eva[0]/eva[1])-1.0e0;
   return ellip;
}
/* ************************************************************************************ */
void GaussWaveFunction::EvalCustomVectorField(double x,double y,double z,\
      double (&v)[3]) {
   /* 
    * double rho;
    * evalRhoGradRho(x,y,z,rho,v);//Stores the Electron density [ED] in rho, and the 
    *                             //gradient in v
    *
    *
    * double led[3],xx[3];
    * xx[0]=x; xx[1]=y; xx[2]=z;
    * evalLED(xx,led);  //stores the vector LED in led
    *
    *
    *
    * */

   /* The following example implements grad(rho)/rho, it may not have a meaning, 
    * and it only has the purpose of showing its implementation  */
   static const double USRFLD_EPS_DEF=1.0e-10; //avoid division by zero
   double rho,gr[3];
   EvalRhoGradRho(x,y,z,rho,gr);
   if ( rho<USRFLD_EPS_DEF ) {rho=USRFLD_EPS_DEF;}
   for ( int i=0 ; i<3 ; i++ ) {v[i]=gr[i]/rho;}
   return;
}
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */


