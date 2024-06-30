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
