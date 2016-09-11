#ifndef _INTEGRATEOVERBONDPATH_H_
#define _INTEGRATEOVERBONDPATH_H_
#include "fldtypesdef.h"
#include <fstream>
using std::ofstream;

/* ************************************************************************** */
class IntegrateOverBondPath {
/* ************************************************************************** */
public:
   IntegrateOverBondPath(class GaussWaveFunction &ugwf,class critPtNetWork &ucpn,\
         ScalarFieldType utp=DENS);
   ~IntegrateOverBondPath();
   void ComputeAllIntegralsOverBondPaths(void);
   solreal GetBondPathIntegral();
   void WriteIntegralValuesToFile(ofstream &ofil);
   inline solreal GetBondPathIntegral(int idx) {return integralValue[idx];}
   inline string GetFieldTypeLabelShort() {return getFieldTypeKeyShort(myCharFieldType);}
   inline string GetFieldTypeLabelLong() {return getFieldTypeKeyLong(myCharFieldType);}
/* ************************************************************************** */
protected:
   class GaussWaveFunction *wf;
   class critPtNetWork *cp;
   /* Functions  */
   void init(void);
   void GetIntermediateCoordinatesAndDistanceBetweenPoints(int startIdx,solreal** (&arr),solreal &h,\
         solreal (&x0)[3],solreal (&x1)[3],solreal (&x2)[3],solreal (&x3)[3]);
   void ComputeScalarFunctionValuesAtIntermediatePoints(solreal (&x0)[3],solreal (&x1)[3],\
         solreal (&x2)[3],solreal (&x3)[3],solreal (&f)[4]);
   solreal ComputeBondPathIntegral(int bpIdx);
   void ComputeBondPathIntegrals(void);
   bool AllocateAuxiliaryArrays(void);
   solreal theFunction(solreal (&x)[3]);
/* ************************************************************************** */
   static constexpr solreal oo3=1.0e0/3.0e0,threeo8=3.0e0/8.0e0;
   static constexpr solreal threeo4=3.0e0/4.0e0;
   static constexpr solreal threeopi=3.0e0/3.14159265358979323846264e0;
   int nbgp;
   ScalarFieldType myFieldType;
   char myCharFieldType;
   solreal *integralValue;
/* ************************************************************************** */
   IntegrateOverBondPath(); //Default constructor is not allowed!
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATEOVERBONDPATH_H_ */

