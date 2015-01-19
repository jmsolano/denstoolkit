#ifndef _CUSTFLD_WFNCLASS_CPP_
#define _CUSTFLD_WFNCLASS_CPP_

/* ************************************************************************************ */
/* In this file, solreal is a mask for a double, i.e., a solreal is a double.
 * This name is used in order to keep compatibility with future implementations
 * requiring floats rather than doubles (old GPUs).  */
solreal gaussWaveFunc::evalCustomScalarField(solreal x,solreal y,solreal z)
{
   /* 
    * solreal rho=evalDensity(x,y,z) //Electron density [ED]
    * solreal maggrho=evalMagGradRho(x,y,z); //Magnitude of the Gradient of ED
    * solreal lap=evalLapRho(x,y,z); //Laplacian of ED
    * solreal lol=evalLOL(x,y,z); //Localized Orbital Locator
    * solreal elf=evalELF(x,y,z); //Returns the Electron Localized Function
    * solreal shent=evalShannonEntropy(x,y,z); //Shannon entropy density
    * solreal mssent=evalMomentumShannonEntropy(px,py,pz); //momentum space shannon
    *                                     //entropy at the momentum-point (px,py,pz)
    * solreal keG=evalKineticEnergyG(x,y,z); //self descriptive
    * solreal keK=evalKineticEnergyK(x,y,z); //self descriptive
    * solreal ftrho=evalFTDensity(px,py,pz); //the momentum-space electron density
    *                                        // at the momentum-point (px,py,pz)
    * solreal mglol=evalMagGradLOL(x,y,z); // Magnitude of grad(LOL)
    * solreal mep=evalMolElecPot(x,y,z); //Molecular Electrostatic Potential
    * solreal magled=evalMagLED(x,y,z); //Magnitude of LED
    * solreal rose=evalRoSE(x,y,z); //Region of Slow Electrons
    * solreal s=evalReducedDensityGradient(x,y,z); //self descriptive
    *
    * */

   /* What follows is an example of how to implement the field rho^2
    *
    * for other fields, you can choose one or more of the above enlisted fields.
    * */
   solreal rho=evalDensity(x,y,z);
   return (rho*rho);
}
/* ************************************************************************************ */
void gaussWaveFunc::evalCustomVectorField(solreal x,solreal y,solreal z,\
      solreal (&v)[3])
{
   /* 
    * solreal rho;
    * evalRhoGradRho(x,y,z,rho,v);//Stores the Electron density [ED] in rho, and the 
    *                             //gradient in v
    *
    *
    * solreal led[3],xx[3];
    * xx[0]=x; xx[1]=y; xx[2]=z;
    * evalLED(xx,led);  //stores the vector LED in led
    *
    *
    *
    * */

   /* The following example implements grad(rho)/rho, it may not have a meaning, 
    * and it only has the purpose of showing its implementation  */
   static const USRFLD_EPS_DEF=1.0e-10; //avoid division by zero
   solreal rho,gr[3];
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

