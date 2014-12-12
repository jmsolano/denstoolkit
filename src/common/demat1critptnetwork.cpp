#ifndef _DEMAT1CRITPTNETWORK_CPP_
#define _DEMAT1CRITPTNETWORK_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <cmath>
#include "demat1critptnetwork.h"
#include "wavefunctionclass.h"
#include "solmemhand.h"
#include "solscrutils.h"
#include "solstringtools.h"


#ifndef MAXSIZECPARRAYS
#define MAXSIZECPARRAYS 20
#endif

#ifndef EPSFABSDIFFCOORD
#define EPSFABSDIFFCOORD (1.0e-06)
#endif

/* ************************************************************************************ */
DeMat1CriticalPointNetwork::DeMat1CriticalPointNetwork()
{
   init();
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::init()
{
   wf=NULL;
   ata=atb=0;
   nACP=nSCP=nRCP=0;
   for ( int i=0 ; i<3 ; i++ ) {x1[i]=x2[i]=x2mx1[i]=x2pmx1p[i]=0.0e0;}
   lenline=0.0e0;
   asACP=asSCP=asRCP=MAXSIZECPARRAYS;
   alloc2DRealArray("RACP",asACP,6,RACP);
   alloc2DRealArray("RSCP",asSCP,6,RSCP);
   alloc2DRealArray("RRCP",asRCP,6,RRCP);
   computePolygonVertices();
   return;
}
/* ************************************************************************************ */
DeMat1CriticalPointNetwork::DeMat1CriticalPointNetwork(gaussWaveFunc *usrwf,\
      int at1,int at2)
{
   init();
   wf=usrwf;
   if ( at1==at2 || at1<0 || at1>=(wf->nNuc) || at2<0 || at2>=(wf->nNuc) ) {
      displayErrorMessage(string("Non valid atom indices! at1: "+getStringFromInt(at1)\
               +", at2: "+getStringFromInt(at2)));
      wf=NULL;
      return;
   }
   ata=at1;
   atb=at2;
   for ( int i=0 ; i<3 ; i++ ) {
      x1[i]=wf->R[3*ata+i];
      x2[i]=wf->R[3*atb+i];
      x2mx1[i]=x2[i]-x1[i];
      x2pmx1p[i]=x2mx1[i];
   }
   lenline=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {lenline+=(x2mx1[i]*x2mx1[i]);}
   lenline=sqrt(lenline);
   oolenline=1.0e0/lenline;
}
/* ************************************************************************************ */
DeMat1CriticalPointNetwork::~DeMat1CriticalPointNetwork()
{
   wf=NULL;
   dealloc2DRealArray(RACP,asACP);
   dealloc2DRealArray(RSCP,asRCP);
   dealloc2DRealArray(RRCP,asSCP);
}
/* ************************************************************************************ */
solreal DeMat1CriticalPointNetwork::PolyV[nPolyV][2];
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::computePolygonVertices(void)
{
   solreal pio6=4.0e0*atan(1.0e0)/6.0e0;
   solreal alpha;
   for ( int i=0 ; i<2 ; i++ ) {PolyV[0][i]=0.0e0;}
   for ( int i=1 ; i<nPolyV ; i++ ) {
      alpha=solreal(i-1)*pio6;
      PolyV[i][0]=cos(alpha);
      PolyV[i][1]=sin(alpha);
   }
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::getXCoordinatesFromUV(solreal uu,solreal vv,\
      solreal (&xx)[3],solreal (&xp)[3])
{
   if ( uu>lenline || uu<0.0e0 ) {
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
   }
   if ( vv>lenline || vv<0.0e0 ) {
      displayErrorMessage("v out of range! v should be: 0<=v<=1");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
   }
   for ( int i=0 ; i<3 ; i++ ) {
      xx[i]=x1[i]+x2mx1[i]*vv;
      xp[i]=x1[i]+x2mx1[i]*uu;
   }
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::evalUVGrad(solreal uu,solreal vv,
      solreal &gamm,solreal (&uvg)[2])
{
   solreal xx[3],xp[3],gg[3],gp[3];
   getXCoordinatesFromUV(uu,vv,xx,xp);
   wf->evalGradDensityMatrix1(xx[0],xx[1],xx[2],xp[0],xp[1],xp[2],gamm,gg,gp);
   solreal sum=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {sum+=x2mx1[i]*gg[i];}
   uvg[0]=sum;
   sum=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {sum+=x2mx1[i]*gp[i];}
   uvg[1]=sum;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::evalUVHessian(solreal uu,solreal vv,solreal &gamm,\
         solreal (&uvg)[2],solreal (&uvh)[2][2])
{
   solreal xx[3],xp[3],gg[3],gp[3],hhhh[3][3],hphh[3][3],hphp[3][3];
   getXCoordinatesFromUV(uu,vv,xx,xp);
   wf->evalHessDensityMatrix1(xx,xp,gamm,gg,gp,hhhh,hphh,hphp);
   solreal sum;
   /* ---------------------------------------------------  */
   sum=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {
      for ( int j=0 ; j<3 ; j++ ) {
         sum+=(x2mx1[i]*x2pmx1p[j]*hhhh[i][j]);
      }
   }
   uvh[0][0]=sum;
   /* ---------------------------------------------------  */
   sum=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {
      for ( int j=0 ; j<3 ; j++ ) {
         sum+=(x2mx1[i]*x2pmx1p[j]*hphp[i][j]);
      }
   }
   uvh[1][1]=sum;
   /* ---------------------------------------------------  */
   sum=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {
      for ( int j=0 ; j<3 ; j++ ) {
         sum+=(x2mx1[i]*x2pmx1p[j]*hphh[i][j]);
      }
   }
   uvh[0][1]=uvh[1][0]=sum;
   /* ---------------------------------------------------  */
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
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */



#endif  /* _DEMAT1CRITPTNETWORK_CPP_ */

