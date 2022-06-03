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
         case ScalarFieldType::REDG :
            return wf->EvalReducedDensityGradient(x,y,z);
            break;
         case ScalarFieldType::ROSE :
            return wf->EvalRoSE(x,y,z);
            break;
         case ScalarFieldType::VPED :
            return wf->EvalVirialPotentialEnergyDensity(x,y,z);
            break;
         case ScalarFieldType::EDFTA :
            trho=wf->EvalDensity(x,y,z);
            return (-0.75e0*pow((0.9549296585513720146e0*trho),1.0e0/3.0e0));
            break;
         case ScalarFieldType::ELLPY :
            return wf->EvalEllipticity(x,y,z);
            break;
         case ScalarFieldType::SCFD :
            return wf->EvalCustomScalarField(x,y,z);
            break;
         default :
            if ( prntunknfld1st ) {
               printf("Unknown field!\n");
               printf("%s, line: %d\n",__FILE__,__LINE__);
               prntunknfld1st=false;
            }
            return 0.0e0;
            break;
      }
      return 0.0e0;
   }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   ScalarFieldType sft;
   GaussWaveFunction* wf;
   static bool prntunknfld1st; /*!< Whether to print error messages.  */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DTKSCALARFUNCTION3D_H_ */

