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
#include "solscrutils.h"
#include "solmemhand.h"
#include "solmath.h"
#include "solfileutils.h"

void IntegrateOverBondPath::init(void) {
   wf=NULL;
   cp=NULL;
   integralValue=NULL;
   nbgp=0;
   myFieldType=NONE;
}
IntegrateOverBondPath::IntegrateOverBondPath(GaussWaveFunction &ugwf,critPtNetWork &ucpn,\
      ScalarFieldType utp) {
   init();
   wf=&ugwf;
   if ( !ucpn.iKnowBCPs() ) {
      displayErrorMessage("First seek the critical points!");
      displayWarningMessage("critPtNetWork pointer is set to null!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   if ( !ucpn.iKnowBGPs() ) {
      displayErrorMessage("First compute the bond paths!");
      displayWarningMessage("critPtNetWork pointer is set to null!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   cp=&ucpn;
   nbgp=cp->nBGP;
   if ( !AllocateAuxiliaryArrays() ) {
      displayErrorMessage("AuxiliaryArrys could not be allocated!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
   myFieldType=utp;
   myCharFieldType=convertScalarFieldType2Char(myFieldType);
}
bool IntegrateOverBondPath::AllocateAuxiliaryArrays(void) {
   return alloc1DRealArray("integralValue",nbgp,integralValue);
}
IntegrateOverBondPath::~IntegrateOverBondPath() {
   dealloc1DRealArray(integralValue);
   wf=NULL;
   cp=NULL;
}
solreal IntegrateOverBondPath::ComputeBondPathIntegral(int bpIdx) {
   solreal xi[3],xip1[3],xip2[3],xip3[3],h,f[4];
   int nReg=cp->conBCP[bpIdx][2]-1;
   solreal res=0.0e0;
   for ( int i=0 ; i<nReg ; ++i ) {
      GetIntermediateCoordinatesAndDistanceBetweenPoints(i,cp->RBGP[bpIdx],h,\
            xi,xip1,xip2,xip3);
      ComputeScalarFunctionValuesAtIntermediatePoints(xi,xip1,xip2,xip3,f);
      res+=(threeo8*h*(f[0]+(3.0e0*f[1])+(3.0e0*f[2])+f[3]));
   }
   return res;
}
void IntegrateOverBondPath::GetIntermediateCoordinatesAndDistanceBetweenPoints(int startIdx,\
      solreal** (&arr),solreal &h,solreal (&x0)[3],solreal (&x1)[3],\
      solreal (&x2)[3],solreal (&x3)[3]) {
   solreal x3mx0[3];
   for ( int i=0 ; i<3 ; ++i ) {
      x0[i]=arr[startIdx][i];
      x3[i]=arr[startIdx+1][i];
      x3mx0[i]=x3[i]-x0[i];
   }
   h=magV3(x3mx0)*oo3;
   normalizeV3(x3mx0);
   for ( int i=0 ; i<3 ; ++i ) {
      x1[i]=x0[i]+(h*x3mx0[i]);
      x2[i]=x1[i]+(h*x3mx0[i]);
   }
}
void IntegrateOverBondPath::ComputeScalarFunctionValuesAtIntermediatePoints(\
      solreal (&x0)[3],solreal (&x1)[3],solreal (&x2)[3],solreal (&x3)[3],\
      solreal (&f)[4]) {
   f[0]=theFunction(x0);
   f[1]=theFunction(x1);
   f[2]=theFunction(x2);
   f[3]=theFunction(x3);
}
void IntegrateOverBondPath::ComputeBondPathIntegrals(void) {
   for ( int i=0 ; i<nbgp ; ++i ) {
      integralValue[i]=ComputeBondPathIntegral(i);
   }
}
solreal IntegrateOverBondPath::GetBondPathIntegral(void) {
   solreal res=0.0e0;
   for ( int i=0 ; i<nbgp ; ++i ) { res+=integralValue[i]; }
   return res;
}
void IntegrateOverBondPath::ComputeAllIntegralsOverBondPaths(void) {
   ComputeBondPathIntegrals();
}
solreal IntegrateOverBondPath::theFunction(solreal (&x)[3]) {
   solreal trho;
   switch ( myFieldType ) {
      case DENS :
         return wf->evalDensity(x[0],x[1],x[2]);
         break;
      case KEDK :
         return wf->evalKineticEnergyK(x[0],x[1],x[2]);
         break;
      case KEDG :
         return wf->evalKineticEnergyG(x[0],x[1],x[2]);
         break;
      case REDG :
         return wf->evalReducedDensityGradient(x[0],x[1],x[2]);
         break;
      case MEPD :
         return wf->evalMolElecPot(x[0],x[1],x[2]);
         break;
      case EDFTA :
         trho=wf->evalDensity(x[0],x[1],x[2]);
         return (-threeo4*pow((threeopi*trho),1.0e0/3.0e0));
         break;
      default :
         break;
   }
   return 0.0e0;
}
void IntegrateOverBondPath::WriteIntegralValuesToFile(ofstream &ofil) {
   writeCommentedScrStarLine(ofil);
   centerCommentedString(GetFieldTypeLabelLong(),ofil);
   writeCommentedScrStarLine(ofil);
   ofil << "Total integral: " << GetBondPathIntegral() << endl;
   for ( int i=0 ; i<nbgp ; ++i ) {
      ofil << cp->lblBCP[i] << " bond path integral: " << GetBondPathIntegral(i) << endl;
   }
}

#endif  /* _INTEGRATEOVERBONDPATH_CPP_ */

