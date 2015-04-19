/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.1
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
#include "eig2-4.h"


#ifndef MAXSIZECPARRAYS
#define MAXSIZECPARRAYS 20
#endif

#ifndef DEMAT1EPSEIGENVALUECPSEARCH
#define DEMAT1EPSEIGENVALUECPSEARCH (1.0e-08)
#endif

#ifndef DEMAT1MAXSTEPSIZEACPRHOSEARCH
#define DEMAT1MAXSTEPSIZEACPRHOSEARCH (0.08)
#endif

#ifndef DEMAT1MAXSTEPSIZESCPRHOSEARCH
#define DEMAT1MAXSTEPSIZESCPRHOSEARCH (0.04)
#endif

#ifndef DEMAT1MAXSTEPSIZERCPRHOSEARCH
#define DEMAT1MAXSTEPSIZERCPRHOSEARCH (0.02)
#endif

#ifndef SIGNF(a)                                                                                          
#define SIGNF(a) ((a)>=0?(1):(-1))
#endif

#ifndef DEMAT1EPSGRADMAG
#define DEMAT1EPSGRADMAG 1.0e-05
#endif

#ifndef DEMAT1EPSFABSDIFFCOORD
#define DEMAT1EPSFABSDIFFCOORD (1.0e-04)
#endif

#ifndef DEMAT1MINGAMMSIGNIFICATIVEVAL
#define DEMAT1MINGAMMSIGNIFICATIVEVAL (1.0e-04)
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
   alloc2DRealArray("RACP",asACP,6,RACP,1.0e+50);
   alloc2DRealArray("RSCP",asSCP,6,RSCP,1.0e+50);
   alloc2DRealArray("RRCP",asRCP,6,RRCP,1.0e+50);
   alloc1DStringArray("lblACP",asACP,lblACP);
   alloc1DStringArray("lblSCP",asSCP,lblSCP);
   alloc1DStringArray("lblRCP",asRCP,lblRCP);
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
   dealloc2DRealArray(RSCP,asSCP);
   dealloc2DRealArray(RRCP,asRCP);
   dealloc1DStringArray(lblACP);
   dealloc1DStringArray(lblSCP);
   dealloc1DStringArray(lblRCP);
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
   if ( uu>1.0e0 || uu<0.0e0 ) {
      //displayErrorMessage("u out of range! u should be: 0<=u<=1");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
   }
   if ( vv>1.0e0 || vv<0.0e0 ) {
      //displayErrorMessage("v out of range! v should be: 0<=v<=1");
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
   for ( int i=0 ; i<3 ; i++ ) {sum+=x2mx1[i]*gg[i];}
   uvg[0]=sum;
   /* ---------------------------------------------------  */
   sum=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {sum+=x2mx1[i]*gp[i];}
   uvg[1]=sum;
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
void DeMat1CriticalPointNetwork::getACPStep(solreal (&g)[2],solreal (&hess)[2][2],\
      solreal (&hh)[2],int &sig)

{
   static solreal eive[2][2],b[2],F[2];
   eigen_decomposition2(hess, eive, b);
   for (int i=0; i<2; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<2; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   static solreal h3[3][3],m3[3][3],v3[3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         h3[i][j]=0.0e0;
         m3[i][j]=0.0e0;
      }
      v3[i]=0.0e0;
   }
   h3[0][0]=b[0]; h3[1][1]=b[1];
   h3[0][2]=h3[2][0]=F[0];
   h3[1][2]=h3[2][1]=F[1];
   eigen_decomposition3(h3, m3, v3);
   solreal lp=v3[2];
#if DEBUG
   if (lp<=0.0e0) {
      displayWarningMessage(string("lp<=0!: "+getStringFromReal(lp)));
      for ( int i=0 ; i<3 ; i++ ) {cout << v3[i] << " ";}
      cout << endl;
      printM2x2Comp("hess:\n",hess);
   }
#endif
   if ( fabs(lp)<DEMAT1EPSEIGENVALUECPSEARCH ) {lp=DEMAT1EPSEIGENVALUECPSEARCH;}
   hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<2; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp);
      hh[1]-=eive[1][j]*F[j]/(b[j]-lp);
   }
   for (int i=0; i<2; i++) {
      if (fabs(hh[i])>DEMAT1MAXSTEPSIZEACPRHOSEARCH) {
         hh[i]=SIGNF(hh[i])*DEMAT1MAXSTEPSIZEACPRHOSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<2; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=3;
         break;
      }
   }
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::getSCPStep(solreal (&g)[2],solreal (&hess)[2][2],\
      solreal (&hh)[2],int &sig)

{
   solreal eive[2][2],b[2],F[2];
   eigen_decomposition2(hess, eive, b);
   for (int i=0; i<2; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<2; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   //cout << "b: " << b[0] << " " << b[1] << endl;
   //cout << "F: " << F[0] << " " << F[1] << endl;
   solreal lp=0.5e0*(b[0]+sqrt(b[0]*b[0]+4.0e0*F[0]*F[0]));
   solreal ln=0.5e0*(b[1]-sqrt(b[1]*b[1]+4.0e0*F[1]*F[1]));
#if DEBUG
   if (lp<=0.0e0) {
      displayWarningMessage(string("lp<=0!: "+getStringFromReal(lp)));
      cout << endl;
   }
   if (ln>=0.0e0) {
      displayWarningMessage(string("ln>=0!: "+getStringFromReal(ln)));
      cout << endl;
   }
   if (b[0]==lp) {
      displayWarningMessage(string("b[0]==lp!"));
      cout << endl;
   }
#endif
   if ( fabs(lp)<DEMAT1EPSEIGENVALUECPSEARCH ) {lp=DEMAT1EPSEIGENVALUECPSEARCH;}
   if ( fabs(ln)<DEMAT1EPSEIGENVALUECPSEARCH ) {ln=DEMAT1EPSEIGENVALUECPSEARCH;}
   hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<2; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln);
   }
   for (int i=0; i<2; i++) {
      if (fabs(hh[i])>DEMAT1MAXSTEPSIZEACPRHOSEARCH) {
         hh[i]=SIGNF(hh[i])*DEMAT1MAXSTEPSIZEACPRHOSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<2; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=-3;
         break;
      }
   }
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::getRCPStep(solreal (&g)[2],solreal (&hess)[2][2],\
      solreal (&hh)[2],int &sig)

{
   static solreal eive[2][2],b[2],F[2];
   eigen_decomposition2(hess, eive, b);
   for (int i=0; i<2; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<2; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   static solreal h3[3][3],m3[3][3],v3[3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         h3[i][j]=0.0e0;
         m3[i][j]=0.0e0;
      }
      v3[i]=0.0e0;
   }
   h3[0][0]=b[0]; h3[1][1]=b[1];
   h3[0][2]=h3[2][0]=F[0];
   h3[1][2]=h3[2][1]=F[1];
   eigen_decomposition3(h3, m3, v3);
   solreal ln=v3[0];
#if DEBUG
   if (ln<=0.0e0) {
      displayWarningMessage(string("ln<=0!: "+getStringFromReal(ln)));
      for ( int i=0 ; i<3 ; i++ ) {cout << v3[i] << " ";}
      cout << endl;
      printM2x2Comp("hess:\n",hess);
   }
#endif
   if ( fabs(ln)<DEMAT1EPSEIGENVALUECPSEARCH ) {ln=DEMAT1EPSEIGENVALUECPSEARCH;}
   hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<2; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-ln);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln);
   }
   for (int i=0; i<2; i++) {
      if (fabs(hh[i])>DEMAT1MAXSTEPSIZERCPRHOSEARCH) {
         hh[i]=SIGNF(hh[i])*DEMAT1MAXSTEPSIZERCPRHOSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<2; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=3;
         break;
      }
   }
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekGammaACP(solreal (&x)[2],solreal &gamm2ret,solreal (&g)[2],\
         int &sig,int maxit)
{
   static solreal gamm,gr[2],hr[2][2],dx[2];
   //static int sig;
   evalUVHessian(x[0],x[1],gamm,gr,hr);
   //cout << "gfromACP: " << g[0] << " " << gr[1] << endl;
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]);
   solreal magh=magd;
   if (magd<=DEMAT1EPSGRADMAG) {
      magd=0.001e0;
      x[0]+=0.001e0;
   }
   int count=0;
   while (((magd>DEMAT1EPSGRADMAG)&&(magh>DEMAT1EPSGRADMAG))&&(count<maxit)) {
      getACPStep(gr,hr,dx,sig);
      //cout << " gfromACP: " << gr[0] << " " << gr[1] << endl;
      //cout << "dxfromACP: " << dx[0] << " " << dx[1] << endl;
      for (int i=0; i<2; i++) {x[i]+=dx[i];}
      evalUVHessian(x[0],x[1],gamm,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
      count++;
   }
   //cout << "iter: " << count << endl;
   //cout << "ACPSig: " << sig << endl;
   for ( int k=0 ; k<2 ; k++ ) {g[k]=gr[k];}
   gamm2ret=gamm;
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekGammaSCP(solreal (&x)[2],solreal &gamm2ret,solreal (&g)[2],\
         int &sig,int maxit)
{
   static solreal gamm,gr[2],hr[2][2],dx[2];
   //static int sig;
   evalUVHessian(x[0],x[1],gamm,gr,hr);
   //cout << "gfromSCP: " << g[0] << " " << gr[1] << endl;
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]);
   solreal magh=magd;
   /*
   if (magd<=DEMAT1EPSGRADMAG) {
      magd=0.001e0;
      x[0]+=0.001e0;
   }
   // */
   int count=0;
   while (((magd>DEMAT1EPSGRADMAG)&&(magh>DEMAT1EPSGRADMAG))&&(count<maxit)) {
      getSCPStep(gr,hr,dx,sig);
      //cout << " gfromSCP: " << gr[0] << " " << gr[1] << endl;
      //cout << "dxfromSCP: " << dx[0] << " " << dx[1] << endl;
      for (int i=0; i<2; i++) {x[i]+=dx[i];}
      evalUVHessian(x[0],x[1],gamm,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
      //cout << "   x: " << x[0] << " " << x[1] << endl;
      //cout << "magd: " << magd << ", magh: " << magh << endl;
      //cout << " eps: " << DEMAT1EPSGRADMAG << endl;
      //cout << (((magd>DEMAT1EPSGRADMAG)&&(magh>DEMAT1EPSGRADMAG))&&(count<maxit) ? "yes" : "no") << endl;
      count++;
   }
   //cout << "iter: " << count << endl;
   //cout << "SCPSig: " << sig << endl;
   for ( int k=0 ; k<2 ; k++ ) {g[k]=gr[k];}
   gamm2ret=gamm;
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekGammaRCP(solreal (&x)[2],solreal &gamm2ret,solreal (&g)[2],\
         int &sig,int maxit)
{
   static solreal gamm,gr[2],hr[2][2],dx[2];
   //static int sig;
   evalUVHessian(x[0],x[1],gamm,gr,hr);
   //cout << "gfromACP: " << g[0] << " " << gr[1] << endl;
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]);
   solreal magh=magd;
   if (magd<=DEMAT1EPSGRADMAG) {
      getRCPStep(gr,hr,dx,sig);
      //magd=0.001e0;
      //x[0]+=0.001e0;
      if (sig==2) {
         cout << "RCPSig: " << sig << ", |g|: " << magd
            << ", it: -1" << endl;
      }
      return;
   }
   int count=0;
   while (((magd>DEMAT1EPSGRADMAG)&&(magh>DEMAT1EPSGRADMAG))&&(count<maxit)) {
      getRCPStep(gr,hr,dx,sig);
      //cout << " gfromRCP: " << gr[0] << " " << gr[1] << endl;
      //cout << "dxfromRCP: " << dx[0] << " " << dx[1] << endl;
      for (int i=0; i<2; i++) {x[i]+=dx[i];}
      evalUVHessian(x[0],x[1],gamm,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
      count++;
   }
   //cout << "iter: " << count << endl;
   /*
   if (sig==2) {
      cout << "RCPSig: " << sig << ", |g|: " << magd
           << ", it: " << count << endl;
   }
   // */
   for ( int k=0 ; k<2 ; k++ ) {g[k]=gr[k];}
   gamm2ret=gamm;
   return;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::setGammaACP(int ndivs)
{
   solreal xx[2],ddx=1.0e0/solreal(ndivs-1); 
   int ntot=ndivs*ndivs,count=1;
   cout << "Processing seed: " << std::flush;
   for ( int i=0 ; i<ndivs ; i++ ) {
         xx[0]=solreal(i)*ddx;
      for ( int j=0 ; j<ndivs ; j++ ) {
         xx[1]=solreal(j)*ddx;
         cout << "\rProcessing seed: " << (ndivs*i+j+1) << "/" << ntot << std::flush;
         seekGammaACPsAroundAPoint(xx,(ddx*0.45e0));
         ++count;
      }
   }
   cout << endl;
  return; 
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::setGammaSCP(int ndivs)
{
   solreal xx[2],ddx=1.0e0/solreal(ndivs-1); 
   int ntot=ndivs*ndivs,count=1;
   cout << "Processing seed: " << std::flush;
   for ( int i=0 ; i<ndivs ; i++ ) {
         xx[0]=solreal(i)*ddx;
      for ( int j=0 ; j<ndivs ; j++ ) {
         xx[1]=solreal(j)*ddx;
         cout << "\rProcessing seed: " << (ndivs*i+j+1) << "/" << ntot << std::flush;
         seekGammaSCPsAroundAPoint(xx,(ddx*0.3e0));
         ++count;
      }
   }
   cout << endl;
  return; 
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::setGammaRCP(int ndivs)
{
   solreal xx[2],ddx=1.0e0/solreal(ndivs-1); 
   int ntot=ndivs*ndivs,count=1;
   cout << "Processing seed: " << std::flush;
   for ( int i=0 ; i<ndivs ; i++ ) {
         xx[0]=solreal(i)*ddx;
      for ( int j=0 ; j<ndivs ; j++ ) {
         xx[1]=solreal(j)*ddx;
         cout << "\rProcessing seed: " << std::setw(4) << (ndivs*i+j+1) << "/" << ntot << std::flush;
         seekGammaRCPsAroundAPoint(xx,(ddx*0.1e0));
         ++count;
      }
   }
   cout << endl;
  return; 
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::addGammaACP(solreal (&x)[2],string lbl)
{
   //if ((x[0]<0.0e0)||(x[0]>1.0e0)) {cout << "Out of box (u)...\n"; return;}
   //if ((x[1]<0.0e0)||(x[1]>1.0e0)) {cout << "Out of box (v)...\n"; return;}
   if ( nACP==0 ) {
      for ( int i=0 ; i<2 ; i++ ) {RACP[0][i]=x[i];}
      lblACP[0]=lbl;
      ++nACP;
      return;
   }
   size_t pos;
   if ( imNew(x,asACP,RACP,pos) ) {
      for ( int i=0 ; i<2 ; i++ ) {RACP[pos][i]=x[i];}
      lblACP[pos]=lbl;
      ++nACP;
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::addGammaSCP(solreal (&x)[2],string lbl)
{
   //if ((x[0]<0.0e0)||(x[0]>1.0e0)) {cout << "Out of box (u)...\n"; return;}
   //if ((x[1]<0.0e0)||(x[1]>1.0e0)) {cout << "Out of box (v)...\n"; return;}
   if ( nSCP==0 ) {
      for ( int i=0 ; i<2 ; i++ ) {RSCP[0][i]=x[i];}
      lblSCP[0]=lbl;
      ++nSCP;
      return;
   }
   size_t pos;
   if ( imNew(x,asSCP,RSCP,pos) ) {
      for ( int i=0 ; i<2 ; i++ ) {RSCP[pos][i]=x[i];}
      lblSCP[pos]=lbl;
      ++nSCP;
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::addGammaRCP(solreal (&x)[2],string lbl)
{
   //if ((x[0]<0.0e0)||(x[0]>1.0e0)) {cout << "Out of box (u)...\n"; return;}
   //if ((x[1]<0.0e0)||(x[1]>1.0e0)) {cout << "Out of box (v)...\n"; return;}
   if ( nRCP==0 ) {
      for ( int i=0 ; i<2 ; i++ ) {RRCP[0][i]=x[i];}
      lblRCP[0]=lbl;
      ++nRCP;
      return;
   }
   size_t pos;
   if ( imNew(x,asRCP,RRCP,pos) ) {
      for ( int i=0 ; i<2 ; i++ ) {RRCP[pos][i]=x[i];}
      lblRCP[pos]=lbl;
      ++nRCP;
   }
}
/* ************************************************************************************ */
bool DeMat1CriticalPointNetwork::imNew(solreal (&x)[2],int dim,solreal ** (&arr),\
      size_t &pos)
{
   solreal ee;
   int firstzeropos=0,k=1;
   ee=arr[0][0]*arr[0][0]+arr[0][1]*arr[0][1];
   while ((ee<1.0e49)&&(k<dim)) {
      ++firstzeropos;
      ee=arr[k][0]*arr[k][0]+arr[k][1]*arr[k][1];
      ++k;
   }
   if (k==dim) {
      cout << "Warning: end of the array reached, perhaps you need a larger array...\n";
      cout << "dim: " << dim << ". Returning false...\n";
      return false;
   }
   k=0;
   while (k<firstzeropos) {
      if ((fabs(x[0]-arr[k][0])<DEMAT1EPSFABSDIFFCOORD)&&
          (fabs(x[1]-arr[k][1])<DEMAT1EPSFABSDIFFCOORD)) {
         pos=k;
         return false;
      }
      ++k;
   }
   pos=firstzeropos;
   return true;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekGammaACPsAroundAPoint(solreal (&oo)[2],solreal ddxx)
{
   solreal xxx[2],gamm,ggg[2];
   string lbl;
   for ( int i=0 ; i<nPolyV ; i++ ) {
      for ( int j=0 ; j<2 ; j++ ) {xxx[j]=oo[j]+PolyV[i][j]*ddxx;}
      lbl="ACP"+getStringFromInt(nACP+1);
      seekSingleGammaACP(xxx,gamm,ggg,lbl);
   }
   
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekGammaSCPsAroundAPoint(solreal (&oo)[2],solreal ddxx)
{
   solreal xxx[2],gamm,ggg[2];
   string lbl;
   for ( int i=0 ; i<nPolyV ; i++ ) {
      for ( int j=0 ; j<2 ; j++ ) {xxx[j]=oo[j]+PolyV[i][j]*ddxx;}
      lbl="SCP"+getStringFromInt(nSCP+1);
      seekSingleGammaSCP(xxx,gamm,ggg,lbl);
   }
   
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekGammaRCPsAroundAPoint(solreal (&oo)[2],solreal ddxx)
{
   solreal xxx[2],gamm,ggg[2];
   string lbl;
   for ( int i=0 ; i<nPolyV ; i++ ) {
      for ( int j=0 ; j<2 ; j++ ) {xxx[j]=oo[j]+PolyV[i][j]*ddxx;}
      lbl="RCP"+getStringFromInt(nRCP+1);
      seekSingleGammaRCP(xxx,gamm,ggg,lbl);
   }
   
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekSingleGammaACP(solreal (&xs)[2],solreal &gamm,\
      solreal (&gg)[2],string &lbl)
{
   /*
   if ( xs[0]<0.0e0 || xs[1]<0.0e0 || xs[0]>1.0e0 || xs[1]>1.0e0 ) {
      gg[0]=gg[1]=1.0e0;
      return;
   }
   // */
   int sig;
   seekGammaACP(xs,gamm,gg,sig);
   solreal magg=getV2Norm(gg);
   if ( magg<=DEMAT1EPSGRADMAG && gamm>DEMAT1MINGAMMSIGNIFICATIVEVAL && (sig==-2) ) {
      //cout << "Possible ACP found: " << xs[0] << " " << xs[1] << endl;
      addGammaACP(xs,lbl);
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekSingleGammaSCP(solreal (&xs)[2],solreal &gamm,\
      solreal (&gg)[2],string &lbl)
{
   /*
   if ( xs[0]<0.0e0 || xs[1]<0.0e0 || xs[0]>1.0e0 || xs[1]>1.0e0 ) {
      gg[0]=gg[1]=1.0e0;
      return;
   }
   // */
   int sig;
   seekGammaSCP(xs,gamm,gg,sig);
   solreal magg=getV2Norm(gg);
   if ( magg<=DEMAT1EPSGRADMAG && gamm>DEMAT1MINGAMMSIGNIFICATIVEVAL && (sig==0) ) {
      //cout << "Possible SCP found: " << xs[0] << " " << xs[1] << "; g: " \
           << gg[0] << " " << gg[1] << ", gm: " << gamm << endl;
      addGammaSCP(xs,lbl);
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::seekSingleGammaRCP(solreal (&xs)[2],solreal &gamm,\
      solreal (&gg)[2],string &lbl)
{
   /*
   if ( xs[0]<0.0e0 || xs[1]<0.0e0 || xs[0]>1.0e0 || xs[1]>1.0e0 ) {
      gg[0]=gg[1]=1.0e0;
      return;
   }
   // */
   int sig;
   seekGammaRCP(xs,gamm,gg,sig);
   solreal magg=getV2Norm(gg);
   if ( magg<=DEMAT1EPSGRADMAG && gamm>DEMAT1MINGAMMSIGNIFICATIVEVAL && (sig==2) ) {
      //cout << "Possible RCP found: " << xs[0] << " " << xs[1] << "; g: " \
           << gg[0] << " " << gg[1] << ", gm: " << gamm << endl;
      addGammaRCP(xs,lbl);
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::displayACPsInfo(void)
{
   if ( nACP==0 ) {
      cout << "No ACPs found!" << endl;
      return;
   }
   solreal gg[2],gamm;
   for ( int i=0 ; i<nACP ; i++ ) {
      evalUVGrad(RACP[i][0],RACP[i][1],gamm,gg);
      cout << "ACP(" << std::setw(3) << std::setfill('0') << i
           << ").\tx:" << RACP[i][0] << " " << RACP[i][1]
           << "; g: " << gg[0] << " " << gg[1] << endl;
   }
   
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::displaySCPsInfo(void)
{
   if ( nSCP==0 ) {
      cout << "No SCPs found!" << endl;
      return;
   }
   solreal gg[2],gamm;
   for ( int i=0 ; i<nSCP ; i++ ) {
      evalUVGrad(RSCP[i][0],RSCP[i][1],gamm,gg);
      cout << "SCP(" << std::setw(3) << std::setfill('0') << i
           << ").\tx:" << RSCP[i][0] << " " << RSCP[i][1]
           << "; g: " << gg[0] << " " << gg[1] << endl;
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::displayRCPsInfo(void)
{
   if ( nRCP==0 ) {
      cout << "No RCPs found!" << endl;
      return;
   }
   solreal gg[2],gamm;
   for ( int i=0 ; i<nRCP ; i++ ) {
      evalUVGrad(RRCP[i][0],RRCP[i][1],gamm,gg);
      cout << "RCP(" << std::setw(3) << std::setfill('0') << i
           << ").\tx:" << RRCP[i][0] << " " << RRCP[i][1]
           << "; g: " << gg[0] << " " << gg[1] << endl;
   }
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::displayCPsInfo(void)
{
   displayACPsInfo();
   displaySCPsInfo();
   displayRCPsInfo();
   cout << "nACP-nSCP+nRCP: " << (nACP-nSCP+nRCP) << endl;
}
/* ************************************************************************************ */
void DeMat1CriticalPointNetwork::setGammaCriticalPoints(void)
{
   cout << "Scanning for ACPs..." << endl;
   setGammaACP();
   cout << "Scanning for SCPs..." << endl;
   setGammaSCP();
   cout << "Scanning for RCPs..." << endl;
   setGammaRCP();
}
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */



#endif  /* _DEMAT1CRITPTNETWORK_CPP_ */

