#ifndef _DTKSCALARFUNCTION6D_H_
#define _DTKSCALARFUNCTION6D_H_
#include <vector>
using std::vector;
#include "../common/fldtypesdef.h"
#include "../common/gausswavefunction.h"

/* ************************************************************************** */
class DTKScalarFunction6D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   DTKScalarFunction6D();
   DTKScalarFunction6D(GaussWaveFunction &ugwf);
   void SetScalarFunction(char t);
   inline double f(const vector<double> &x) { return f(x[0],x[1],x[2],x[3],x[4],x[5]); }
   inline double f(const double x,const double y,const double z,\
                   const double xp,const double yp,const double zp) {
      switch ( sft ) {
         case ScalarField6DType::DM1 :
            return wf->EvalDensityMatrix1(x,y,z,xp,yp,zp);
            break;
         default :
            if ( prntunknfld1st ) {
               printf("\033[31mUnknown or unimplemented 6D field!\n");
               printf("%s, line: %d\033[0m\n",__FILE__,__LINE__);
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
   ScalarField6DType sft;
   GaussWaveFunction* wf;
   static bool prntunknfld1st; /*!< Whether to print error messages.  */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DTKSCALARFUNCTION6D_H_ */

