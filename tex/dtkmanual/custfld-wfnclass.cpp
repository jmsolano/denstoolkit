/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
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



#ifndef _CUSTFLD_WFNCLASS_CPP_
#define _CUSTFLD_WFNCLASS_CPP_

/* ************************************************************************************ */
/* In this file, double is a mask for a double, i.e., a double is a double.
 * This name is used in order to keep compatibility with future implementations
 * requiring floats rather than doubles (old GPUs).  */
double GaussWaveFunction::evalCustomScalarField(double x,double y,double z)
{
   /* 
    * double rho=evalDensity(x,y,z) //Electron density [ED]
    * double maggrho=evalMagGradRho(x,y,z); //Magnitude of the Gradient of ED
    * double lap=evalLapRho(x,y,z); //Laplacian of ED
    * double lol=evalLOL(x,y,z); //Localized Orbital Locator
    * double elf=evalELF(x,y,z); //Returns the Electron Localized Function
    * double shent=evalShannonEntropy(x,y,z); //Shannon entropy density
    * double mssent=evalMomentumShannonEntropy(px,py,pz); //momentum space shannon
    *                                     //entropy at the momentum-point (px,py,pz)
    * double keG=evalKineticEnergyG(x,y,z); //self descriptive
    * double keK=evalKineticEnergyK(x,y,z); //self descriptive
    * double ftrho=evalFTDensity(px,py,pz); //the momentum-space electron density
    *                                        // at the momentum-point (px,py,pz)
    * double mglol=evalMagGradLOL(x,y,z); // Magnitude of grad(LOL)
    * double mep=evalMolElecPot(x,y,z); //Molecular Electrostatic Potential
    * double magled=evalMagLED(x,y,z); //Magnitude of LED
    * double rose=evalRoSE(x,y,z); //Region of Slow Electrons
    * double s=evalReducedDensityGradient(x,y,z); //self descriptive
    *
    * */

   /* What follows is an example of how to implement the field rho^2
    *
    * for other fields, you can choose one or more of the above enlisted fields.
    * */
   double rho=evalDensity(x,y,z);
   return (rho*rho);
}
/* ************************************************************************************ */
void GaussWaveFunction::evalCustomVectorField(double x,double y,double z,\
      double (&v)[3])
{
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
   evalRhoGradRho(x,y,z,rho,gr);
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


#endif  /* _CUSTFLD_WFNCLASS_CPP_ */

