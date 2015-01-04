/*
 *  critptnetwork.cpp
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 05/06/13.
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef _CRITPTNETWORK_CPP_
#define _CRITPTNETWORK_CPP_

#include "solscrutils.h"
#include "solfileutils.h"
#include "critptnetwork.h"
#include "bondnetwork.h"
#include "eig2-4.h"
#include "iofuncts-cpx.h"
#include "atomradiicust.h"
#include "solstringtools.h"
// The first 94 atomic radii are given,
//  the rest are set to be 0.80e0
//
//#include "atomcolschcust.h" //Choose this for the palette defined by JMHP
#include "atomcolschjmol.h" //Choose this for the palette used in JMol

/* ********************************************************************************* */
/* ********************************************************************************* */
/* ********************************************************************************* */

/* ************************************************************************************* */
critPtNetWork::critPtNetWork()
{
   //publics:
   dACP=dBCP=dRCP=dCCP=0;
   nACP=nBCP=nRCP=nCCP=0;
   nBGP=0;
   atBCP=NULL;
   RACP=RBCP=RRCP=RCCP=NULL;
   RBGP=NULL;
   lblACP=lblBCP=lblRCP=lblCCP=NULL;
   for (int i=0; i<3; i++) {centMolecVec[i]=0.0e0;}
   //privates:
   iknowacps=iknowbcps=iknowrcps=iknowccps=false;
   iknowallcps=false;
   iknowbgps=false;
   mycptype=NONE;
   drawNuc=false;
   drawBnd=true;
   drawBGPs=false;
   tubeBGPStyle=false;
   //maxiterracp=MAXITERATIONACPSEARCH;
   //maxiterrbcp=MAXITERATIONBCPSEARCH;
   //maxiterrrcp=MAXITERATIONRCPSEARCH;
   //maxiterrccp=MAXITERATIONCCPSEARCH;
   return;
}
critPtNetWork::~critPtNetWork()
{
   dealloc3DRealArray(RBGP,dBCP,ARRAYSIZEGRADPATH);
   dealloc2DRealArray(RACP,dACP);
   dealloc1DStringArray(lblACP);
   dealloc2DRealArray(RBCP,dBCP);
   dealloc1DStringArray(lblBCP);
   dealloc2DIntArray(atBCP,dBCP);
   dealloc2DRealArray(RRCP,dRCP);
   dealloc1DStringArray(lblRCP);
   dealloc2DRealArray(RCCP,dCCP);
   dealloc1DStringArray(lblCCP);
}
/* ************************************************************************************* */
solreal critPtNetWork::V0=0.0e0;
solreal critPtNetWork::V5=2.0e0/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
solreal critPtNetWork::V8=(1.0e0+sqrt(5.0e0))/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
solreal critPtNetWork::IHV[nIHV][3]={
   {  V0,  V0,  V0 },
   {  V0,  V0,  2.0e0*EPSFABSDIFFCOORD },
   {  V0, -V5,  V8 },
   {  V8,  V0,  V5 },
   {  V8,  V0, -V5 },
   { -V8,  V0, -V5 },
   { -V8,  V0,  V5 },
   { -V5,  V8,  V0 },
   {  V5,  V8,  V0 },
   {  V5, -V8,  V0 },
   { -V5, -V8,  V0 },
   {  V0, -V5, -V8 },
   {  V0,  V5, -V8 },
   {  V0,  V5,  V8 }
};
//The first point is the origin plus a small epsilon, not a vertex.
//Added for looking at the atomic center and for easy the
//addition of the labels and coordinates in
//functions add*?CP
/* ************************************************************************************* */
void critPtNetWork::setCriticalPoints(bondNetWork &bn,gaussWaveFunc &wf,ScalarFieldType ft)
                            //ft\in{NONE,DENS,MGRD,LAPD,LOLD,ELFD,SENT,KEDK,KEDG}
{
   mycptype=ft;
   if (!bn.imstp()) {
      displayErrorMessage("Trying to use a non set up bond network object!");
      exit(1);
   }
   switch (ft) {
      case DENS:
         cout << "Scanning for Density Critical Points." << endl;
         dACP=bn.nNuc*MAXRHOACPSPERATOM;
         break;
      case LOLD:
         cout << "Scanning for LOL Critical Points." << endl;
         dACP=bn.nNuc*MAXLOLACPSPERATOM;
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   //int eACP=nACP+bn.countAtomsOfAtomicNumber(1);
   alloc2DRealArray(string("RACP"),dACP,3,RACP,1.0e+50);
   alloc1DStringArray("lblACP",dACP,lblACP);
   cout << "Looking for Attractor Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowacps=setRhoACPs(bn,wf);
         //dBCP=(eACP*(eACP-1))/2;
         dBCP=(nACP*(nACP-1))/2;
         break;
      case LOLD:
         iknowacps=setLOLACPs(bn,wf);
         dBCP=(nACP*(nACP-1))/2;
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if (iknowacps) {
      alloc2DRealArray(string("RBCP"),dBCP,3,RBCP,1.0e+50);
      alloc2DIntArray(string("atBCP"),dBCP,3,atBCP,-1);
      alloc1DStringArray("lblBCP",dBCP,lblBCP);
   } else {
      displayErrorMessage("First look for ACPs...\nQuiting...");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
   }
#if DEBUG
   cout << "nACP: " << nACP << ", dACP: " << dACP << endl;
#endif
   cout << "Looking for Bond Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowbcps=setRhoBCPs(bn,wf);
         break;
      case LOLD:
         iknowbcps=setLOLBCPs(bn,wf);
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if (iknowbcps) {
      dRCP=(nBCP*(nBCP-1))/2;
      alloc2DRealArray(string("RRCP"),dRCP,3,RRCP,1.0e+50);
      alloc1DStringArray("lblRCP",dRCP,lblRCP);
   }
#if DEBUG
   cout << "nBCP: " << nBCP << ", dBCP: " << dBCP << endl;
#endif
   cout << "Looking for Ring Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowrcps=setRhoRCPs(bn,wf);
         break;
      case LOLD:
         iknowrcps=setLOLRCPs(bn,wf);
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if (iknowrcps) {
      dCCP=(nRCP*(nRCP-1))/2;
      alloc2DRealArray(string("RCCP"),dCCP,3,RCCP,1.0e+50);
      alloc1DStringArray("lblCCP",dCCP,lblCCP);
   }
#if DEBUG
   cout << "nRCP: " << nRCP << ", dRCP: " << dRCP << endl;
#endif
   cout << "Looking for Cage Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowccps=setRhoCCPs(bn,wf);
         break;
      case LOLD:
         iknowccps=setLOLCCPs(bn,wf);
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
#if DEBUG
   cout << "nCCP: " << nCCP << ", dCCP: " << dCCP << endl;
#endif
   iknowallcps=(iknowacps&&iknowbcps&&iknowrcps&&iknowccps);
   if (iknowallcps) {
      printScrCharLine('*');
      cout << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
      printScrCharLine('*');
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekRhoACP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],\
      gaussWaveFunc &wf,int maxit)
{
   static solreal rho,gr[3],hr[3][3],dx[3];
   static int sig;
   wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
   solreal magh=magd;
   if (magd<=EPSGRADMAG) {
      magd=0.1e0;
      x[0]+=0.01e0;
   }
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&(count<maxit)) {
      getACPStep(gr,hr,dx,sig);
      //if (sig!=-3) {
      //   x[0]+=0.1e0;
      //   magd=0.1e0;
      //   count++;
      //   continue;
      //}
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      count++;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekRhoBCP(solreal (&x)[3],gaussWaveFunc &wf)
{
   static solreal rho,gr[3],hr[3][3],dx[3];
   static int sig;
   wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
   if (magd<=EPSGRADMAG) {
      magd=0.1e0;
      x[0]+=0.05e0;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&
          (count<MAXITERATIONBCPSEARCH)&&(rho>MINRHOSIGNIFICATIVEVAL)) {
      getBCPStep(gr,hr,dx,sig);
      //if (sig!=-1) {
      //   x[0]+=0.1e0;
      //   magd=0.1e0;
      //   count++;
      //   continue;
      //}
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      count++;
   }
   //cout << "IterBCP: " << count << endl;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekRhoRCP(solreal (&x)[3],gaussWaveFunc &wf)
{
   static solreal rho,gr[3],hr[3][3],dx[3];
   static int sig;
   wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
   if (magd<=EPSGRADMAG) {
      magd=0.1e0;
      x[0]+=0.05e0;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&
          (count<MAXITERATIONRCPSEARCH)&&(rho>MINRHOSIGNIFICATIVEVAL)) {
      getRCPStep(gr,hr,dx,sig);
      //cout << "rho: " << rho << ", dx: " << dx[0] << " " << dx[1] << " " << dx[2]  << endl;
      //if (sig!=1) {
      //   cout << "sigRCP: " << sig << endl;
      //   x[0]+=0.1e0;
      //   magd=0.1e0;
      //   count++;
      //   continue;
      //}
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
      //cout << "rho: " << rho << ", R: " << x[0] << " " << x[1] << " " << x[2]  << endl;
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      count++;
   }
   //cout << "IterRCP: " << count << endl;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekRhoCCP(solreal (&x)[3],gaussWaveFunc &wf)
{
   static solreal rho,gr[3],hr[3][3],dx[3];
   static int sig;
   wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
   solreal magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
   if (magd<=EPSGRADMAG) {
      magd=0.1e0;
      x[0]+=0.05e0;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&
          (count<MAXITERATIONRCPSEARCH)&&(rho>MINRHOSIGNIFICATIVEVAL)) {
      getCCPStep(gr,hr,dx,sig);
      //if (sig!=3) {
      //   x[0]+=0.1e0;
      //   magd=0.1e0;
      //   count++;
      //   continue;
      //}
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf.evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      count++;
   }
   //cout << "IterRCP: " << count << endl;
   //displayWarningMessage("seekRhoCCP(...) under construction.");
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekLOLACP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf)
{
   static solreal lol,gl[3],hl[3][3],dx[3];
   static int sig;
   wf.evalHessLOL(x,lol,gl,hl);
   solreal magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
   solreal magh=magd;
   if (magd<=EPSGRADMAG) {
      magd=0.1e0;
      x[0]+=0.05e0;
   }
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&(count<MAXITERATIONACPSEARCH)) {
      getACPStep(gl,hl,dx,sig);
      //if (sig!=-3) {
      //   x[0]+=0.1e0;
      //   magd=0.1e0;
      //   count++;
      //   continue;
      //}
      for (int i=0; i<3; i++) {
         if (fabs(dx[i])>MAXSTEPSIZEACPLOLSEARCH) {
            dx[i]=SIGNF(dx[i])*MAXSTEPSIZEACPLOLSEARCH;
         }
      }
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf.evalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      count++;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   //cout << "LOLACP: " << x[0] << " " << x[1] << " " << x[2] << endl;
   //cout << "   lol: " << lol << endl;
   //cout << "    gl: " << gl[0] << " " << gl[1] << " " << gl[2] << endl;
   //cout << "  Iter: " << count << endl;
   //return;
   //displayWarningMessage("seekLOLACP(...) under construction.");
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekLOLBCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf)
{
   static solreal lol,gl[3],hl[3][3],dx[3];
   static int sig;
   wf.evalHessLOL(x,lol,gl,hl);
   solreal magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
   solreal magh=magd;
   if (magd<=EPSGRADMAG) {
      magd=0.1e0;
      x[0]+=0.05e0;
   }
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&(count<MAXITERATIONACPSEARCH)) {
      getBCPStep(gl,hl,dx,sig);
      //if (sig!=-1) {
      //   x[0]+=0.1e0;
      //   magd=0.1e0;
      //   count++;
      //   continue;
      //}
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEBCPLOLSEARCH*/);}
      wf.evalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      count++;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   //cout << "LOLBCP: " << x[0] << " " << x[1] << " " << x[2] << endl;
   //cout << "   lol: " << lol << endl;
   //cout << "    gl: " << gl[0] << " " << gl[1] << " " << gl[2] << endl;
   //cout << "  Iter: " << count << endl;
   //displayWarningMessage("seekLOLBCP(...) under construction.");
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekLOLRCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf)
{
   displayWarningMessage("seekLOLRCP(...) under construction.");
   return;
}
/* ************************************************************************************* */
void critPtNetWork::seekLOLCCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf)
{
   displayWarningMessage("seekLOLCCP(...) under construction.");
   return;
}
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
void critPtNetWork::getACPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   static solreal eive[3][3],b[3],F[3];
   eigen_decomposition3(hess, eive, b);
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   static solreal h4[4][4],m4[4][4],v4[4];
   for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
         h4[i][j]=0.0e0;
         m4[i][j]=0.0e0;
      }
      v4[i]=0.0e0;
   }
   h4[0][0]=b[0]; h4[1][1]=b[1]; h4[2][2]=b[2];
   h4[0][3]=h4[3][0]=F[0];
   h4[1][3]=h4[3][1]=F[1];
   h4[2][3]=h4[3][2]=F[2];
   eigen_decomposition4(h4, m4, v4);
   solreal lp=v4[3];
   if ( fabs(lp)<EPSEIGENVALUECPSEARCH ) {lp=EPSEIGENVALUECPSEARCH;}
#if DEBUG
   if (lp<=0.0e0) {
      displayWarningMessage(string("lp<=0!: "+getStringFromReal(lp)));
      for ( int i=0 ; i<4 ; i++ ) {cout << v4[i] << " ";}
      cout << endl;
      printM3x3Comp("hess:\n",hess);
   }
#endif
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp+EPSGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-lp+EPSGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-lp+EPSGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZEACPRHOSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZEACPRHOSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<3; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=0;
         break;
      }
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getBCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   static solreal eive[3][3],b[3];
   eigen_decomposition3(hess, eive, b);
   static solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   static solreal h3[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {h3[i][j]=0.00000e0;}
   }
   solreal ln=0.5e0*(b[2]-sqrt(b[2]*b[2]+4.0e0*F[2]*F[2]));
   h3[0][0]=b[0];
   h3[1][1]=b[1];
   h3[2][0]=h3[0][2]=F[0];
   h3[2][1]=h3[1][2]=F[1];
   static solreal m3[3][3],vv[3];
   eigen_decomposition3(h3, m3, vv);
   solreal lp=vv[2];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp+EPSGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-lp+EPSGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln+EPSGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZEBCPSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZEBCPSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<3; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=0;
         break;
      }
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getRCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   static solreal eive[3][3],b[3];
   eigen_decomposition3(hess, eive, b);
   static solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   static solreal h3[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {h3[i][j]=0.00000e0;}
   }
   solreal lp=0.5e0*(b[0]+sqrt(b[0]*b[0]+4.0e0*F[0]*F[0]));
   h3[0][0]=b[1];
   h3[1][1]=b[2];
   h3[2][0]=h3[0][2]=F[1];
   h3[2][1]=h3[1][2]=F[2];
   static solreal m3[3][3],vv[3];
   eigen_decomposition3(h3, m3, vv);
   solreal ln=vv[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp+EPSGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln+EPSGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln+EPSGRADMAG);
   }
   //cout << "hh: " << hh[0] << " " << hh[1] << " " << hh[2] <<endl;
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZERCPSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZERCPSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<3; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=0;
         break;
      }
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getCCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   static solreal eive[3][3],b[3];
   eigen_decomposition3(hess, eive, b);
   static solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   static solreal h4[4][4],m4[4][4],v4[4];
   for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
         h4[i][j]=0.0e0;
         m4[i][j]=0.0e0;
      }
      v4[i]=0.0e0;
   }
   h4[0][0]=b[0]; h4[1][1]=b[1]; h4[2][2]=b[2];
   h4[0][3]=h4[3][0]=F[0];
   h4[1][3]=h4[3][1]=F[1];
   h4[2][3]=h4[3][2]=F[2];
   eigen_decomposition4(h4, m4, v4);
   solreal ln=v4[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-ln+EPSGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln+EPSGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln+EPSGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZECCPSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZECCPSEARCH;
      }
   }
   sig=0;
   for (int i=0; i<3; i++) {
      sig+=SIGNF(b[i]);
      if (b[i]==0.0e0) {
         sig=0;
         break;
      }
   }
   return;
}
/* ************************************************************************************* */
bool critPtNetWork::setRhoACPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   solreal rad,x[3],rho,g[3],magg;
   string lbl;
   nACP=0;
   //nACP=bn.nNuc;
   for (int i=0; i<bn.nNuc; i++) {
      //for (int j=0; j<3; j++) {RACP[i][j]=bn.R[i][j];}
      for ( int j=0 ; j<3 ; j++ ) {x[j]=bn.R[i][j];}
      seekRhoACP(x,rho,g,wf);
      magg=0.0e0;
      for ( int k=0 ; k<3 ; k++ ) {magg+=(g[k]*g[k]);}
      magg=sqrt(magg);
      //cout << "magg: " << magg << " (" << bn.atLbl[i] << ")" << endl;
      if ( (rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSRHOACPGRADMAG) ) {
         lbl=bn.atLbl[i];
         addRhoACP(x,lbl,bn);
      }
      //lblACP[i]=bn.atLbl[i];
   }
   solreal regnACP=nACP;
   cout << "Found " << regnACP << " regular ACPs" << endl;
   int nhyd,noth;
   nhyd=bn.countAtomsOfAtomicNumber(1);
   noth=bn.nNuc-nhyd;
   solreal perchyd,percoth;
   perchyd=100.0e0*solreal(nhyd)/solreal(noth);
   percoth=100.0e0-perchyd;
   cout << "Looking around atoms..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   int itmp;
   for (int i=0; i<bn.nNuc; i++) {
      //if (bn.atNum[i]!=0) {continue;}
      rad=0.2e0;//atomicRadius[bn.atNum[i]]*0.1e0;
      lbl=(bn.atLbl[i]+string("nn"));
      itmp=0;
      for (int j=0; j<6; j++) {
         for (int k=0; k<3; k++) {x[k]=bn.R[i][k]+IHV[j][k]*rad;}
         seekRhoACP(x,rho,g,wf);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if ((rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSRHOACPGRADMAG)) {
            lbl+=getStringFromInt(itmp);
            addRhoACP(x,lbl,bn);
            ++itmp;
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((bn.nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   cout << "Looking for possible extra ACPs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal extbd,maxbd,dd;
   maxbd=bn.maxBondDist;
   extbd=maxbd*1.1e0;
   //cout << "extd: " << extbd << endl;
   for (int i=0; i<bn.nNuc; i++) {
      for (int j=i+1; j<bn.nNuc; j++) {
         dd=0.0e0;
         for (int k=0; k<3; k++) {
            dd+=((bn.R[i][k]-bn.R[j][k])*(bn.R[i][k]-bn.R[j][k]));
         }
         dd=sqrt(dd);
         if (dd<=extbd) {
            //cout << "Possible extended ACP Seed found!" << endl;
            for (int k=0; k<3; k++) {x[k]=0.5e0*(bn.R[i][k]+bn.R[j][k]);}
            wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            magg=0.0e0;
            for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
            magg=sqrt(magg);
            if ((rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSGRADMAG)) {
               addRhoACP(x,lbl,bn);
            }
            seekRhoACP(x,rho,g,wf,(MAXITERATIONACPSEARCH));
            wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            magg=0.0e0;
            for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
            magg=sqrt(magg);
            lbl="nnn";
            if ((rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSGRADMAG)) {
               addRhoACP(x,lbl,bn);
            }
}
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setRhoBCPs(bondNetWork &bn,gaussWaveFunc &wf)
{
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   int ata,atb,k,mypos;
   solreal x[3],g[3],rho;
   string lbl;
   for (int i=0; i<bn.nNuc; i++) {
      ata=i;
      atb=bn.bNet[ata][0];
      k=1;
      while ((k<MAXBONDINGATOMS)&&(atb>0)) {
         for (int j=0; j<3; j++) {x[j]=0.5e0*(bn.R[ata][j]+bn.R[atb][j]);}
         seekRhoBCP(x,wf);
         //cout << "BCP: " << x[0] << " " << x[1] << " " << x[2] << endl;
         //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         if (ata<atb) {
            lbl=bn.atLbl[ata]+string("-")+bn.atLbl[atb];
         } else {
            lbl=bn.atLbl[atb]+string("-")+bn.atLbl[ata];
         }
         mypos=addRhoBCP(x,lbl,bn);
         if (mypos>=0&&mypos<dBCP) {
            atBCP[mypos][0]=ata;
            atBCP[mypos][1]=atb;
         }
         atb=bn.bNet[ata][k];
         k++;
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
   int normalbcp=nBCP;
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   //------------------------------------------------------
   cout << nBCP << " BCPs found.\n";
   solreal extbd,dd,magg;
   string ll;
   cout << "Including non-nuclear attractors..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for ( int i=0 ; i<nACP ; i++ ) {
      for ( int j=bn.nNuc ; j<nACP ; j++ ) {
         dd=0.0e0;
         for (int k=0; k<3; k++) {
            dd+=((RACP[i][k]-RACP[j][k])*(RACP[i][k]-RACP[j][k]));
         }
         dd=sqrt(dd);
         if (dd<=bn.maxBondDist) {
            //cout << "Possible extended BCP Seed found!" << endl;
            //cout << bn.atLbl[i] << "-" << bn.atLbl[j] << endl;
            //cx=0.5e0*(bn.R[i][0]+bn.R[j][0]);
            //cy=0.5e0*(bn.R[i][1]+bn.R[j][1]);
            //cz=0.5e0*(bn.R[i][2]+bn.R[j][2]);
            for (int k=0; k<3; k++) {x[k]=0.5e0*(RACP[i][k]+RACP[j][k]);}
            if ((wf.evalDensity(x[0],x[1],x[2])>MINRHOSIGNIFICATIVEVAL)) {
               seekRhoBCP(x,wf);
               if (i<j) {
                  lbl=string("*")+lblACP[i]+string("-")+lblACP[j];
               } else {
                  lbl=string("*")+lblACP[j]+string("-")+lblACP[i];
               }
               wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
               magg=0.0e0;
               for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
               magg=sqrt(magg);
               if ((rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSGRADMAG)) {
                  mypos=addRhoBCP(x,lbl,bn);
                  if ((mypos>=0)&&(mypos<dBCP)) {
                     atBCP[mypos][0]=ata;
                     atBCP[mypos][1]=atb;
                  }
                  //cout << x[0] << " " << x[1] << " " << x[2] << endl;
                  //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
                  //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
               }
            }
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
#endif
   cout << endl << "Looking for possible extra BCPs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   extbd=bn.maxBondDist*EXTENDEDBONDDISTFACTOR;
   //cout << "extd: " << extbd << endl;
   for (int i=0; i<bn.nNuc; i++) {
      for (int j=i+1; j<bn.nNuc; j++) {
         dd=0.0e0;
         for (int k=0; k<3; k++) {
            dd+=((bn.R[i][k]-bn.R[j][k])*(bn.R[i][k]-bn.R[j][k]));
         }
         dd=sqrt(dd);
         if ((dd>=bn.maxBondDist)&&(dd<=extbd)) {
            //cout << "Possible extended BCP Seed found!" << endl;
            //cout << bn.atLbl[i] << "-" << bn.atLbl[j] << endl;
            //cx=0.5e0*(bn.R[i][0]+bn.R[j][0]);
            //cy=0.5e0*(bn.R[i][1]+bn.R[j][1]);
            //cz=0.5e0*(bn.R[i][2]+bn.R[j][2]);
            for (int k=0; k<3; k++) {x[k]=0.5e0*(bn.R[i][k]+bn.R[j][k]);}
            if ((wf.evalDensity(x[0],x[1],x[2])>MINRHOSIGNIFICATIVEVAL)) {
               seekRhoBCP(x,wf);
               if (i<j) {
                  lbl=string("*")+bn.atLbl[i]+string("-")+bn.atLbl[j];
               } else {
                  lbl=string("*")+bn.atLbl[j]+string("-")+bn.atLbl[i];
               }
               wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
               magg=0.0e0;
               for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
               magg=sqrt(magg);
               if ((rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSGRADMAG)) {
                  mypos=addRhoBCP(x,lbl,bn);
                  if ((mypos>=0)&&(mypos<dBCP)) {
                     atBCP[mypos][0]=ata;
                     atBCP[mypos][1]=atb;
                  }
                  //cout << x[0] << " " << x[1] << " " << x[2] << endl;
                  //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
                  //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
               }
            }
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   if (normalbcp==nBCP) {
      cout << "No more BCPs found." << endl;
   } else {
      solreal tbcp[3];
      for ( int i=normalbcp ; i<nBCP ; i++ ) {
         for ( int k=0 ; k<3 ; k++ ) {tbcp[k]=RBCP[i][k];}
         findTwoClosestAtoms(tbcp,wf,ata,atb);
         atBCP[i][0]=ata;
         atBCP[i][1]=atb;
         lblBCP[i]="*"+wf.atLbl[ata]+"-"+wf.atLbl[atb];
      }
      cout << (nBCP-normalbcp) << " new BCP";
      if ((nBCP-normalbcp)>1) {cout << "s";}
      cout << " found! Total number of BCPs: " << nBCP << endl;
   }
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setRhoRCPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   solreal x[3],magg,rho,g[3];
   string lbl;
   for (int i=0; i<nBCP; i++) {
      for (int j=(i+1); j<nBCP; j++) {
         magg=0.0e0;
         for (int k=0; k<3; k++) {
            x[k]=0.5e0*(RBCP[i][k]+RBCP[j][k]);
            magg+=(RBCP[i][k]-RBCP[j][k])*(RBCP[i][k]-RBCP[j][k]);
         }
         magg=sqrt(magg);
         if (magg>(bn.maxBondDist*2.0)) {continue;}
         seekRhoRCP(x,wf);
         lbl=(lblBCP[i]+string("-")+lblBCP[j]);
         wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if (/*(rho>MINRHOSIGNIFICATIVEVAL)&&*/(magg<EPSGRADMAG)) {
            addRhoRCP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
            //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nBCP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   //cout << "Checkpoint..." << endl;
   for (int i=0; i<nRCP; i++) {
      removeRedundInLabel(lblRCP[i]); //cout << "i: " << i << endl;
   }
   cout << nRCP << " RCPs found.\n";
   //displayWarningMessage("setRhoRCPs(...) under construction.");
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setRhoCCPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   solreal x[3],magg,rho,g[3];
   string lbl;
   for (int i=0; i<nRCP; i++) {
      for (int j=(i+1); j<nRCP; j++) {
         magg=0.0e0;
         for (int k=0; k<3; k++) {
            x[k]=0.5e0*(RRCP[i][k]+RRCP[j][k]);
            magg+=(RRCP[i][k]-RRCP[j][k])*(RRCP[i][k]-RRCP[j][k]);
         }
         magg=sqrt(magg);
         if (magg>(bn.maxBondDist*BONDISTEXTCCPSEARCHFACTOR)) {continue;}
         seekRhoCCP(x,wf);
         lbl=(lblRCP[i]+string("-")+lblRCP[j]);
         wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if (/*(rho>MINRHOSIGNIFICATIVEVAL)&&*/(magg<EPSGRADMAG)) {
            addRhoCCP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
            //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nRCP))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   for (int i=0; i<nCCP; i++) {
      removeRedundInLabel(lblCCP[i]);
      cout << lblCCP[i] << endl;
   }
   cout << nCCP << " CCPs found.\n";
   //displayWarningMessage("setRhoCCPs(...) under construction.");
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setLOLACPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   //string cl="a";
   solreal x[3],magg,lol,g[3],rad;
   string lbl;
   int nhyd,noth;
   nhyd=bn.countAtomsOfAtomicNumber(1);
   noth=bn.nNuc-nhyd;
   solreal perchyd,percoth;
   perchyd=100.0e0*solreal(nhyd)/solreal(noth);
   percoth=100.0e0-perchyd;
   //cout << "There are " << nhyd << " hydrogens" << endl;
   cout << "Looking around Hydrogen atoms..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for (int i=0; i<bn.nNuc; i++) {
      if (bn.atNum[i]!=0) {continue;}
      rad=0.1e0;//atomicRadius[bn.atNum[i]]*0.1e0;
      lbl=(bn.atLbl[i]+string("a"));
      for (int j=0; j<nIHV; j++) {
         for (int k=0; k<3; k++) {x[k]=bn.R[i][k]+IHV[j][k]*rad;}
         seekLOLACP(x,lol,g,wf);
         //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         //cout << "glolnuc: " << g[0] << " " << g[1] << " " << g[2] << endl;
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if (/*(wf.evalDensity(x[0],x[1],x[2])>MINRHOSIGNIFICATIVEVAL)&&*/(magg<EPSGRADMAG)) {
            addLOLACP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((bn.nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   cout << "Looking around other atoms..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for (int i=0; i<bn.nNuc; i++) {
      if (bn.atNum[i]==0) {continue;}
      //rad=atomicRadius[bn.atNum[i]]*2.0e0/3.0e0;
      rad=0.01e0;
      lbl=(bn.atLbl[i]+string("a"));
      for (int j=0; j<nIHV; j++) {
         for (int k=0; k<3; k++) {x[k]=bn.R[i][k]+IHV[j][k]*rad;}
         seekLOLACP(x,lol,g,wf);
         //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if (/*(wf.evalDensity(x[0],x[1],x[2])>MINRHOSIGNIFICATIVEVAL)&&*/(magg<EPSGRADMAG)) {
            addLOLACP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
            //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         }
         //*
         for (int k=0; k<3; k++) {x[k]=bn.R[i][k]+IHV[j][k]*rad*40.0e0;}
         seekLOLACP(x,lol,g,wf);
         //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if (/*(wf.evalDensity(x[0],x[1],x[2])>MINRHOSIGNIFICATIVEVAL)&&*/(magg<EPSGRADMAG)) {
            addLOLACP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
            //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         }
         // */
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((bn.nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   cout << "Looking between atom pairs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   
   //-----------------------------------------------
   int ata,atb,k;
   for (int i=0; i<bn.nNuc; i++) {
      ata=i;
      atb=bn.bNet[ata][0];
      k=1;
      while ((k<MAXBONDINGATOMS)&&(atb>=0)) {
         for (int j=0; j<3; j++) {x[j]=0.5e0*(bn.R[ata][j]+bn.R[atb][j]);}
         seekLOLACP(x,lol,g,wf);
         //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         //cout << "BCP: " << x[0] << " " << x[1] << " " << x[2] << endl;
         //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         if (ata<atb) {
            lbl=bn.atLbl[ata]+string("_")+bn.atLbl[atb];
         } else if (ata==atb) {
            lbl=bn.atLbl[ata]+string("a");
         } else {
            lbl=bn.atLbl[atb]+string("_")+bn.atLbl[ata];
         }
         addLOLACP(x,lbl,bn);
         atb=bn.bNet[ata][k];
         k++;
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((bn.nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   cout << nACP << " ACPs found.\n";
   //return true;
   //displayWarningMessage("setLOLACPs(...) under construction.");
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setLOLBCPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   solreal x[3],magg,lol,g[3];
   string lbl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for (int i=0; i<nACP; i++) {
      for (int j=(i+1); j<nACP; j++) {
         magg=0.0e0;
         for (int k=0; k<3; k++) {
            x[k]=0.5e0*(RACP[i][k]+RACP[j][k]);
            magg+=(RACP[i][k]-RACP[j][k])*(RACP[i][k]-RACP[j][k]);
         }
         magg=sqrt(magg);
         //if (magg>bn.maxBondDist) {continue;}
         seekLOLBCP(x,lol,g,wf);
         lbl=(lblACP[i]+string("-")+lblACP[j]);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if (/*(rho>MINRHOSIGNIFICATIVEVAL)&&*/(magg<EPSGRADMAG)) {
            addLOLBCP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
            //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
   //for (int i=0; i<nBCP; i++) {removeRedundInLabel(lblRCP[i]);}
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   /*
   cout << "Looking around ACPs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal rad;
   for (int i=0; i<bn.nNuc; i++) {
      //rad=atomicRadius[bn.atNum[i]]*2.0e0/3.0e0;
      rad=0.1e0;
      lbl=(bn.atLbl[i]+string("a"));
      for (int j=0; j<nIHV; j++) {
         for (int k=0; k<3; k++) {x[k]=bn.R[i][k]+IHV[j][k]*rad;}
         seekLOLBCP(x,lol,g,wf);
         //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
         magg=0.0e0;
         for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
         magg=sqrt(magg);
         if ((magg<EPSGRADMAG)) {
            addLOLBCP(x,lbl,bn);
            //cout << x[0] << " " << x[1] << " " << x[2] << endl;
            //wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
            //cout << "     " << g[0] << " " << g[1] << " " << g[2] << endl;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(solreal(i)/solreal((bn.nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   // */
   cout << nBCP << " BCPs found.\n";
   //displayWarningMessage("setLOLBCPs(...) under construction.");
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setLOLRCPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   displayWarningMessage("No LOL RCP will be look for...");
   displayWarningMessage("(critPtNetWork::setLOLRCPs(...) under construction.)");
   return false;
}
/* ************************************************************************************* */
bool critPtNetWork::setLOLCCPs(bondNetWork &bn,gaussWaveFunc &wf)
{
   displayWarningMessage("No LOL CCP will be look for...");
   displayWarningMessage("(critPtNetWork::setLOLCCPs(...) under construction.)");
   return false;
}
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
void critPtNetWork::displayACPCoords(void)
{
   if (!iknowacps) {
      displayErrorMessage("First seek the ACPs.\nNo ACP to display its coordinates.");
      return;
   }
   cout << scientific << setprecision(12);
   printScrCharLine('+');
   centerString("Coordinates of Attractor Critical Points");
   printScrCharLine('-');
   for (int i=0; i<nACP; i++) {
      cout << "ACP[" << (i+1) << "](" << lblACP[i] << "): ";
      for (int j=0; j<3; j++) {cout << RACP[i][j] << " ";}
      cout << endl;
   }
   printScrCharLine('+');
   cout.unsetf(ios::scientific);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::displayBCPCoords(void)
{
   if (!iknowbcps) {
      displayErrorMessage("First seek the BCPs.\nNo BCP to display its coordinates.");
      return;
   }
   cout << scientific << setprecision(12);
   printScrCharLine('+');
   centerString("Coordinates of Bond Critical Points");
   printScrCharLine('-');
   for (int i=0; i<nBCP; i++) {
      cout << "BCP[" << (i+1) << "](" << lblBCP[i] << "): ";
      for (int j=0; j<3; j++) {cout << RBCP[i][j] << " ";}
      cout << endl;
   }
   printScrCharLine('+');
   cout.unsetf(ios::scientific);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::displayRCPCoords(void)
{
   if (!iknowrcps) {
      displayErrorMessage("First seek the RCPs.\nNo RCP to display its coordinates.");
      return;
   }
   cout << scientific << setprecision(11);
   printScrCharLine('+');
   centerString("Coordinates of Ring Critical Points");
   printScrCharLine('-');
   for (int i=0; i<nRCP; i++) {
      cout << "RCP[" << (i+1) << "](" << lblRCP[i] << "): ";
      for (int j=0; j<3; j++) {cout << RRCP[i][j] << " ";}
      cout << endl;
   }
   printScrCharLine('+');
   cout.unsetf(ios::scientific);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::displayCCPCoords(void)
{
   if (!iknowccps) {
      displayErrorMessage("First seek the CCPs.\nNo CCP to display its coordinates.");
      return;
   }
   cout << scientific << setprecision(12);
   printScrCharLine('+');
   centerString("Coordinates of Cage Critical Points");
   printScrCharLine('-');
   for (int i=0; i<nCCP; i++) {
      cout << "CCP[" << (i+1) << "](" << lblCCP[i] << "): ";
      for (int j=0; j<3; j++) {cout << RCCP[i][j] << " ";}
      cout << endl;
   }
   printScrCharLine('+');
   cout.unsetf(ios::scientific);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::removeRedundInLabel(string &lbl)
{
   /*
   string tl,sl,fl;
   tl=lbl+"-";
   sl=fl="";
   cout << tl << endl;
   size_t pos=tl.find_first_of("-")+1,len=tl.length(),clen;
#if DEBUG
   int count=0;
#endif
   while ((pos!=string::npos)&&(pos<len)&&(count<50)) {
      if ((pos!=string::npos)&&(pos<len)) {sl=tl.substr(0,pos);} else {sl="";}
      clen=sl.length(); pos=0;
      cout << "clen=" << clen << "; sl='" << sl << "'" << endl;
      while (pos!=string::npos&&pos<len) {
         if (tl.length()==0||clen==0) {break;}
         tl.erase(pos,clen);
         pos=tl.find(sl);
      }
      fl+=sl;
      //cout << tl << endl;
      pos=tl.find_first_of("-")+1;
      if (pos==tl.length()) {fl+=tl.substr(0,tl.length()); break;}
      count++;
   }
   clen=fl.length()-1;
   if (fl[clen]=='-') {fl.erase(clen,1);}
   cout << fl << endl;
   lbl=fl;
   // */
   string tl=lbl,chunk,tchk,tl2,fl="";
   chunk=getFirstChunkOfLabel(tl);
   while (chunk.length()>0) {
      fl+=(chunk+"-");
      tl.erase(0,(chunk.length()+1));
      tchk=getFirstChunkOfLabel(tl);
      tl2="";
      while (tchk.length()>0) {
         if (tchk!=chunk) {
            tl2+=(tchk+"-");
         }
         tl.erase(0,(tchk.length()+1));
         tchk=getFirstChunkOfLabel(tl);
      }
      tl=tl2;
      chunk=getFirstChunkOfLabel(tl);
   }
   size_t clen=fl.length()-1;
   if (fl[clen]=='-') {fl.erase(clen,1);}
   //cout << fl << endl;
   lbl=fl;
   return;
}
/* ************************************************************************************* */
string critPtNetWork::getFirstChunkOfLabel(string &lbl)
{
   if (lbl.length()==0) {
      return string("");
   }
   size_t pos=lbl.find_first_of("-");
   return lbl.substr(0,pos);
}
/* ************************************************************************************* */
void critPtNetWork::addRhoACP(solreal (&x)[3],string &lbl,bondNetWork &bn)
{
   if (nACP==0) {
      for (int i=0; i<3; i++) {RACP[0][i]=x[i];}
      lblACP[0]=lbl;
      (lbl[lbl.length()-1])++; //Check what is this for...
      nACP++;
      return;
   }
   if ((x[0]<bn.bbmin[0])||(x[0]>bn.bbmax[0])) {cout << "Out of box (x)...\n"; return;}
   if ((x[1]<bn.bbmin[1])||(x[1]>bn.bbmax[1])) {cout << "Out of box (y)...\n"; return;}
   if ((x[2]<bn.bbmin[2])||(x[2]>bn.bbmax[2])) {cout << "Out of box (z)...\n"; return;}
   size_t pos;
   if (imNew(x,dACP,RACP,pos)) {
      for (int i=0; i<3; i++) {RACP[pos][i]=x[i];}
      lblACP[pos]=lbl;
      (lbl[lbl.length()-1])++;
      nACP++;
   }
   return;
}
/* ************************************************************************************* */
int critPtNetWork::addRhoBCP(solreal (&x)[3],string &lbl,bondNetWork &bn)
{
   if (nBCP==0) {
      for (int i=0; i<3; i++) {RBCP[0][i]=x[i];}
      lblBCP[0]=lbl;
      nBCP++;
      return 0;
   }
   if ((x[0]<bn.bbmin[0])||(x[0]>bn.bbmax[0])) {/*cout << "Out of box (x)...\n";*/ return -1;}
   if ((x[1]<bn.bbmin[1])||(x[1]>bn.bbmax[1])) {/*cout << "Out of box (y)...\n";*/ return -1;}
   if ((x[2]<bn.bbmin[2])||(x[2]>bn.bbmax[2])) {/*cout << "Out of box (z)...\n";*/ return -1;}
   size_t pos;
   if (imNew(x,dBCP,RBCP,pos)) {
      for (int i=0; i<3; i++) {RBCP[pos][i]=x[i];}
      lblBCP[pos]=lbl;
      nBCP++;
   } else {
      return -1;
   }
   return int(pos);
}
/* ************************************************************************************* */
void critPtNetWork::addRhoRCP(solreal (&x)[3],string &lbl,bondNetWork &bn)
{
   if (nRCP==0) {
      for (int i=0; i<3; i++) {RRCP[0][i]=x[i];}
      lblRCP[0]=lbl;
      nRCP++;
      return;
   }
   if ((x[0]<bn.bbmin[0])||(x[0]>bn.bbmax[0])) {/*cout << "\nOut of box (x)... " << lbl;*/ return;}
   if ((x[1]<bn.bbmin[1])||(x[1]>bn.bbmax[1])) {/*cout << "\nOut of box (y)... " << lbl;*/ return;}
   if ((x[2]<bn.bbmin[2])||(x[2]>bn.bbmax[2])) {/*cout << "\nOut of box (z)... " << lbl;*/ return;}
   size_t pos;
   if (imNew(x,dRCP,RRCP,pos)) {
      for (int i=0; i<3; i++) {RRCP[pos][i]=x[i];}
      lblRCP[pos]=lbl;
      nRCP++;
   } else {
      //cout << pos << " " << endl;
      //cout << "lbl from addRhoRCP: " << lbl << endl;
      if ((int(pos)<nRCP)&&(pos!=string::npos)) {
         lblRCP[pos]+=(string("-")+lbl);
         //cout << "added: " << lbl << endl;
      }
      //lblRCP[pos]+=(string("-")+lbl);
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::addRhoCCP(solreal (&x)[3],string &lbl,bondNetWork &bn)
{
   if (nCCP==0) {
      for (int i=0; i<3; i++) {RCCP[0][i]=x[i];}
      lblCCP[0]=lbl;
      nCCP++;
      return;
   }
   if ((x[0]<bn.bbmin[0])||(x[0]>bn.bbmax[0])) {/*cout << "\nOut of box (x)... " << lbl;*/ return;}
   if ((x[1]<bn.bbmin[1])||(x[1]>bn.bbmax[1])) {/*cout << "\nOut of box (y)... " << lbl;*/ return;}
   if ((x[2]<bn.bbmin[2])||(x[2]>bn.bbmax[2])) {/*cout << "\nOut of box (z)... " << lbl;*/ return;}
   size_t pos;
   if (imNew(x,dCCP,RCCP,pos)) {
      for (int i=0; i<3; i++) {RCCP[pos][i]=x[i];}
      lblCCP[pos]=lbl;
      nCCP++;
   } else {
      //cout << pos << " " << endl;
      //cout << "lbl from addRhoRCP: " << lbl << endl;
      if ((int(pos)<nCCP)&&(pos!=string::npos)) {
         lblCCP[pos]+=(string("-")+lbl);
         //cout << "added: " << lbl << endl;
      }
      //lblRCP[pos]+=(string("-")+lbl);
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::addLOLACP(solreal (&x)[3],string &lbl,bondNetWork &bn)
{
   if (nACP==0) {
      for (int i=0; i<3; i++) {RACP[0][i]=x[i];}
      lblACP[0]=lbl;
      (lbl[lbl.length()-1])++;
      nACP++;
      return;
   }
   if ((x[0]<bn.bbmin[0])||(x[0]>bn.bbmax[0])) {/*cout << "Out of box (x)...\n";*/ return;}
   if ((x[1]<bn.bbmin[1])||(x[1]>bn.bbmax[1])) {/*cout << "Out of box (y)...\n";*/ return;}
   if ((x[2]<bn.bbmin[2])||(x[2]>bn.bbmax[2])) {/*cout << "Out of box (z)...\n";*/ return;}
   size_t pos;
   if (imNew(x,dACP,RACP,pos)) {
      for (int i=0; i<3; i++) {RACP[pos][i]=x[i];}
      lblACP[pos]=lbl;
      (lbl[lbl.length()-1])++;
      nACP++;
   }
   return;
}
/* ************************************************************************************* */
void critPtNetWork::addLOLBCP(solreal (&x)[3],string &lbl,bondNetWork &bn)
{
   if (nBCP==0) {
      for (int i=0; i<3; i++) {RBCP[0][i]=x[i];}
      lblBCP[0]=lbl;
      nBCP++;
      return;
   }
   if ((x[0]<bn.bbmin[0])||(x[0]>bn.bbmax[0])) {/*cout << "Out of box (x)...\n";*/ return;}
   if ((x[1]<bn.bbmin[1])||(x[1]>bn.bbmax[1])) {/*cout << "Out of box (y)...\n";*/ return;}
   if ((x[2]<bn.bbmin[2])||(x[2]>bn.bbmax[2])) {/*cout << "Out of box (z)...\n";*/ return;}
   size_t pos;
   if (imNew(x,dBCP,RBCP,pos)) {
      for (int i=0; i<3; i++) {RBCP[pos][i]=x[i];}
      lblBCP[pos]=lbl;
      nBCP++;
   }
   return;
}
/* ************************************************************************************* */
bool critPtNetWork::imNew(solreal (&x)[3],int dim,solreal ** (&arr),size_t &pos)
{
   solreal ee;
   int firstzeropos=0,k=1;
   ee=arr[0][0]*arr[0][0]+arr[0][1]*arr[0][1]+arr[0][2]*arr[0][2];
   while ((ee<1.0e49)&&(k<dim)) {
      firstzeropos++;
      ee=arr[k][0]*arr[k][0]+arr[k][1]*arr[k][1]+arr[k][2]*arr[k][2];
      k++;
   }
   if (k==dim) {
      cout << "Warning: end of the array reached, perhaps you need a larger array...\n";
      cout << "Returning false...\n";
      return false;
   }
   k=0;
   while (k<firstzeropos) {
      if ((fabs(x[0]-arr[k][0])<EPSFABSDIFFCOORD)&&
          (fabs(x[1]-arr[k][1])<EPSFABSDIFFCOORD)&&
          (fabs(x[2]-arr[k][2])<EPSFABSDIFFCOORD)) {
         pos=k;
#if DEBUG
         //cout << "Already exists!" << endl;
         //DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
         return false;
      }
      k++;
   }
   pos=firstzeropos;
   return true;
}
/* ************************************************************************************* */
void critPtNetWork::displayIHVCoords(void)
{
   cout << scientific << setprecision(8);
   printScrCharLine('+');
   centerString("Coordinates of Icosahedron Vertices");
   printScrCharLine('-');
   for (int i=0; i<nIHV; i++) {
      cout << "V[" << (i+1) << "]: ";
      for (int j=0; j<3; j++) {cout << IHV[i][j] << " ";}
      cout << endl;
   }
   printScrCharLine('+');
   cout.unsetf(ios::scientific);
   return;
}
/* ********************************************************************************* */
void critPtNetWork::printAllFieldProperties(solreal &x,solreal &y,solreal &z,gaussWaveFunc &wf)
{
   static solreal rho,lol,xx[3],g[3],hess[3][3];
   static solreal eivec[3][3],eival[3];
   xx[0]=x;
   xx[1]=y;
   xx[2]=z;
   wf.evalRhoGradRho(x,y,z,rho,g);
   wf.evalHessian(xx[0],xx[1],xx[2],hess);
   eigen_decomposition3(hess, eivec, eival);
   cout << "  R:           " << x << " " << y << " " << z << endl;
   cout << "  Rho:         " << rho << endl;
   cout << "  GradRho:     " << g[0] << " " << g[1] << " " << g[2] << endl;
   cout << "  HessRho:     " << hess[0][0] << " " << hess[0][1] << " " << hess[0][2] << endl;
   cout << "               " << hess[1][0] << " " << hess[1][1] << " " << hess[1][2] << endl;
   cout << "               " << hess[2][0] << " " << hess[2][1] << " " << hess[2][2] << endl;
   cout << "  EigVal Hess: " << eival[0] << " " << eival[1] << " " << eival[2] << endl;
   cout << "  EigVec Hess: " << eivec[0][0] << " " << eivec[1][0] << " " << eivec[2][0] << endl;
   cout << "               " << eivec[0][1] << " " << eivec[1][1] << " " << eivec[2][1] << endl;
   cout << "               " << eivec[0][2] << " " << eivec[1][2] << " " << eivec[2][2] << endl;
   cout << "  LapRho:      " << (hess[0][0]+hess[1][1]+hess[2][2]) << endl;
   wf.evalHessLOL(xx,lol,g,hess);//(x,lol,gl,hl)
   eigen_decomposition3(hess, eivec, eival);
   cout << "  LOL:         " << lol << endl;
   cout << "  gradLOL:     " << g[0] << " " << g[1] << " " << g[2] << endl;
   cout << "  HessLOL:     " << hess[0][0] << " " << hess[0][1] << " " << hess[0][2] << endl;
   cout << "               " << hess[1][0] << " " << hess[1][1] << " " << hess[1][2] << endl;
   cout << "               " << hess[2][0] << " " << hess[2][1] << " " << hess[2][2] << endl;
   cout << "  EigVal HLOL: " << eival[0] << " " << eival[1] << " " << eival[2] << endl;
   cout << "  ELF:         " << wf.evalELF(x,y,z) << endl;
   cout << "  K.E. G.:     " << wf.evalKineticEnergyG(x,y,z) << endl;
   cout << "  K.E. K.:     " << wf.evalKineticEnergyK(x,y,z) << endl;
   cout << "  Shann. Ent.: " << wf.evalShannonEntropy(x,y,z) << endl;
   return;
}
/* ********************************************************************************* */
void critPtNetWork::writeAllFieldProperties(ofstream &ofil,solreal &x,solreal &y,solreal &z,gaussWaveFunc &wf)
{
   static solreal rho,lol,xx[3],g[3],hess[3][3];
   static solreal eivec[3][3],eival[3];
   xx[0]=x;
   xx[1]=y;
   xx[2]=z;
   wf.evalRhoGradRho(x,y,z,rho,g);
   wf.evalHessian(xx[0],xx[1],xx[2],hess);
   eigen_decomposition3(hess, eivec, eival);
   ofil << scientific << setprecision(12);
   ofil << "  R:           ";
   for (int i=0; i<3; i++) {
      ofil.width(20);
      ofil << xx[i];
   }
   ofil << endl;
   ofil << "  Rho:         "; ofil.width(20); ofil  << rho << endl;
   ofil << "  GradRho:     ";
   for (int i=0; i<3; i++) {
      ofil.width(20);
      ofil << g[i];
   }
   ofil << endl;
   for (int i=0; i<3; i++) {
      if (i==0) {
         ofil << "  HessRho:     ";
      } else {
         ofil << "               ";
      }
      for (int j=0; j<3; j++) {
         ofil.width(20);
         ofil << hess[i][j];
      }
      ofil << endl;
   }
   ofil << "  EigVal Hess: ";
   for (int i=0; i<3; i++) {
      ofil.width(20);
      ofil << eival[i];
   }
   ofil << endl;
   for (int i=0; i<3; i++) {
      if (i==0) {
         ofil << "  EigVec Hess: ";
      } else {
         ofil << "               ";
      }
      for (int j=0; j<3; j++) {ofil.width(20); ofil << eivec[j][i];}
      ofil << endl;
   }
   ofil << "  LapRho:      "; ofil.width(20);
   ofil << (hess[0][0]+hess[1][1]+hess[2][2]) << endl;
   wf.evalHessLOL(xx,lol,g,hess);//(x,lol,gl,hl)
   eigen_decomposition3(hess, eivec, eival);
   ofil << "  LOL:         "; ofil.width(20); ofil << lol << endl;
   ofil << "  gradLOL:     ";
   for (int i=0; i<3; i++) {
      ofil.width(20);
      ofil << g[i];
   }
   ofil << endl;
   for (int i=0; i<3; i++) {
      if (i==0) {
         ofil << "  HessLOL:     ";
      } else {
         ofil << "               ";
      }
      for (int j=0; j<3; j++) {
         ofil.width(20);
         ofil << hess[i][j];
      }
      ofil << endl;
   }
   ofil << "  EigVal HLOL: ";
   for (int i=0; i<3; i++) {
      ofil.width(20);
      ofil << eival[i];
   }
   ofil << endl;
   ofil << "  ELF:         "; ofil.width(20); ofil << wf.evalELF(x,y,z) << endl;
   ofil << "  K.E. G.:     "; ofil.width(20); ofil << wf.evalKineticEnergyG(x,y,z) << endl;
   ofil << "  K.E. K.:     "; ofil.width(20); ofil << wf.evalKineticEnergyK(x,y,z) << endl;
   ofil << "  Shann. Ent.: "; ofil.width(20); ofil << wf.evalShannonEntropy(x,y,z) << endl;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::printCPProps(gaussWaveFunc &wf)
{
   solreal x,y,z;//,gx,gy,gz,rho,hxx,hyy,hzz,hxy,hxz,hyz;
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   string cptp;
   switch (mycptype) {
      case DENS:
         cptp="Density ";
         break;
      case LOLD:
         cptp="LOL ";
      default:
         break;
   }
   cout << cptp << "Attractor Critical Points Information" << endl;
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   if (iknowacps) {
      for (int i=0; i<80; i++) {cout << '-';}
      cout << endl;
      for (int i=0; i<nACP; i++) {
         cout << "ACP(" << (i+1) << ") [" << lblACP[i] << "]:" << endl;
         x=RACP[i][0];
         y=RACP[i][1];
         z=RACP[i][2];
         printAllFieldProperties(x,y,z,wf);
         for (int i=0; i<80; i++) {cout << '-';}
         cout << endl;
      }
   } else {
      cout << "No ACP information available!" << endl;
   }
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   cout << cptp << "Bond Critical Points Information" << endl;
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   if (iknowbcps) {
      for (int i=0; i<80; i++) {cout << '-';}
      cout << endl;
      for (int i=0; i<nBCP; i++) {
         cout << "BCP(" << (i+1) << ") [" << lblBCP[i] << "]:" << endl;
         x=RBCP[i][0];
         y=RBCP[i][1];
         z=RBCP[i][2];
         printAllFieldProperties(x,y,z,wf);
         for (int i=0; i<80; i++) {cout << '-';}
         cout << endl;
      }
   } else {
      cout << "No BCP information available!" << endl;
   }
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   cout << cptp << "Ring Critical Points Information" << endl;
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   if (iknowrcps) {
      for (int i=0; i<80; i++) {cout << '-';}
      cout << endl;
      if (nRCP==0) {
         cout << "No Ring Critical Points have been found." << endl;
         for (int i=0; i<80; i++) {cout << '-';}
         cout << endl;
      }
      for (int i=0; i<nRCP; i++) {
         cout << "RCP(" << (i+1) << ") [" << lblRCP[i] << "]:" << endl;
         x=RRCP[i][0];
         y=RRCP[i][1];
         z=RRCP[i][2];
         printAllFieldProperties(x,y,z,wf);
         for (int i=0; i<80; i++) {cout << '-';}
         cout << endl;
      }
   } else {
      cout << "No RCP information available!" << endl;
   }
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   cout << cptp << "Cage Critical Points Information" << endl;
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   if (iknowccps) {
      for (int i=0; i<80; i++) {cout << '-';}
      cout << endl;
      if (nCCP==0) {
         cout << "No Cage Critical Points have been found." << endl;
         for (int i=0; i<80; i++) {cout << '-';}
         cout << endl;
      }
      for (int i=0; i<nCCP; i++) {
         cout << "CCP(" << (i+1) << ") [" << lblCCP[i] << "]:" << endl;
         x=RCCP[i][0];
         y=RCCP[i][1];
         z=RCCP[i][2];
         printAllFieldProperties(x,y,z,wf);
         for (int i=0; i<80; i++) {cout << '-';}
         cout << endl;
      }
   } else {
      cout << "No CCP information available!" << endl;
   }
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   cout << cptp << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
   for (int i=0; i<80; i++) {cout << '*';}
   cout << endl;
   return;
}
/* ************************************************************************************* */
bool critPtNetWork::makePOVFile(string pnam,bondNetWork &bn,povRayConfProp &pvp,
                                int campos)
{
   ofstream pof;
   pof.open(pnam.c_str(),ios::out);
   if (!(pof.good())) {
      cout << "Error: File " << pnam << "could not be opened...\n";
#if DEBUG
      cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
      return false;
   }
   string tmpbool;
#if CHOOSEPOVVERSION36
   pof << "#version 3.6; //Unless you know what you are doing, do not modify this line..." << endl;
#endif
   pof << "#include \"colors.inc\"" << endl;
   writeScrCharLine(pof,'/');
   pof << "//" << endl;
#if DEBUG
   pof << "//Code generated by the class critPtNetWork." << endl;
#endif
   pof << "//Below you can find some options to be parsed to povray" << endl;
   pof << "//set your custom values." << endl;
   pof << "//You can reconstruct the image using the script dtkpov2png" << endl;
   pof << "//" << endl;
   writeScrCharLine(pof,'/');
   pof << "#declare GNUPlotAngle1=0;" << endl;
   pof << "#declare GNUPlotAngle2=0;" << endl;
   pof << "#declare YAngle=0;" << endl;
   if (drawNuc) {tmpbool="true";} else {tmpbool="false";}
   pof << "#declare DrawAtomTranspSpheres=" << tmpbool << ";" << endl;
   if (drawBnd) {tmpbool="true";} else {tmpbool="false";}
   pof << "#declare DrawStandardBonds=" << tmpbool << ";" << endl;
   pof << "#declare DrawAttractorCriticalPoints=true;" << endl;
   pof << "#declare DrawBondCriticalPoints=true;" << endl;
   pof << "#declare DrawRingCriticalPoints=true;" << endl;
   pof << "#declare DrawCageCriticalPoints=true;" << endl;
   if (drawBGPs&&(!tubeBGPStyle)) {tmpbool="true";} else {tmpbool="false";}
   pof << "#declare DrawGradientPathSpheres=" << tmpbool << ";" << endl;
   if (drawBGPs&&tubeBGPStyle) {tmpbool="true";} else {tmpbool="false";}
   pof << "#declare DrawGradientPathTubes=" << tmpbool << ";" << endl;
   pof << "//  Activation of \"DrawGradientPathSpheres\" requires deactivation of "
   << "\n//  \"DrawGradientPathTubes\", and vice versa." << endl;
   solreal allcprad=bn.drawAtSize*ATOMCRITICALPOINTSIZEFACTOR;
   pof << "#declare RadiusAllCriticalPoints=" << allcprad << ";" << endl;
   pof << "#declare ColorACP=rgb <0.0,0.0,0.0>;" << endl;
   pof << "#declare RadiusACP=RadiusAllCriticalPoints;" << endl;
   //pof << "#declare ColorBCP=rgb <0.3,0.3,0.3>;" << endl;
   pof << "#declare ColorBCP=rgb <0.0,0.6,1.0>;" << endl;
   pof << "#declare RadiusBCP=RadiusAllCriticalPoints;" << endl;
   //pof << "#declare ColorRCP=rgb <0.6,0.6,0.6>;" << endl;
   pof << "#declare ColorRCP=rgb <1.0,1.0,0.0>;" << endl;
   pof << "#declare RadiusRCP=RadiusAllCriticalPoints;" << endl;
   //pof << "#declare ColorCCP=rgb <0.9,0.9,0.9>;" << endl;
   pof << "#declare ColorCCP=rgb <1.0,0.0,0.0>;" << endl;
   pof << "#declare RadiusCCP=RadiusAllCriticalPoints;" << endl;
   pof << "#declare ColorABGradPath=rgb <0.0,1.0,0.0>;" << endl;
   pof << "#default { finish { specular 0.3 roughness 0.03 phong .1 } }" << endl;
   writeScrCharLine(pof,'/');
   pof << "//For the colors, instead of rgb <...>, you may want to try Red, Yellow, ..." << endl;
   pof << "//  or any of the colors defined in \"colors.inc\"" << endl;
   pof << "//" << endl;
   writeScrCharLine(pof,'/');
   pof << "// END OF CUSTOM OPTIONS" << endl;
   writeScrCharLine(pof ,'/');
   if (!(bn.imstp())) {bn.setUpBNW();}
   if (!(bn.ballAndStickMode)) {bn.drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
#if DEBUG
   displayWarningMessage("In this version, calling critPtNetWork::makePovFILE(...)\n\
will overwrite the original coordinates of the critical points\n\
and the coordinates on the bondnetwork object as well.");
#endif
   centerMolecule(bn);
   bn.calcViewRadius();
   //cout << "rView: " << bn.rView << endl;
   solreal camdist=2.5e0;
   for (int i=0; i<3; i++) {pvp.locCam[i]=0.0e0;}
   switch (campos) {
      case 1:
         pvp.locCam[2]=camdist;
         break;
      case 10:
         pvp.locCam[1]=camdist;
         break;
      case 100:
         pvp.locCam[0]=camdist;
         break;
      case 11:
         pvp.locCam[2]=pvp.locCam[1]=camdist;
         break;
      case 110:
         pvp.locCam[0]=pvp.locCam[1]=camdist;
         break;
      case 101:
         pvp.locCam[0]=pvp.locCam[2]=camdist;
         break;
      case 111:
         for (int i=0; i<3; i++) {
            pvp.locCam[i]=camdist;
         }
         break;
      default:
         displayWarningMessage("Not a valid direction. Using 001.");
         pvp.locCam[2]=camdist;
         break;
   }
   for (int i=0; i<3; i++) {
      pvp.locCam[i]*=bn.rView;
      for (int j=0; j<2; j++) {
         pvp.lightSource[j][i]*=(bn.rView*2.0e0);
      }
   }
   pvp.inccolors=false;
   //pvp.writeHeader(pof,false);
   pof << "global_settings { ambient_light White }" << endl;
   pof << "\nbackground { color < 0, 0.5, 0.7 > }\n" << endl;
   solreal zsep=0.5e0;
   pvp.lightSource[1][0]=zsep;
   pvp.lightSource[1][1]=zsep;
   pvp.lightSource[1][2]=1.0e0;
   pvp.addLightSource(zsep,-zsep,1.0e0);
   pvp.addLightSource(-zsep,zsep,1.0e0);
   pvp.addLightSource(-zsep,-zsep,1.0e0);
   for (int i=1; i<pvp.nLightSources; i++) {
      for (int j=0; j<3; j++) {pvp.lightSource[i][j]*=(bn.rView*4.0e0);}
   }
   for (int i=0; i<pvp.nLightSources; i++) {
      pvp.writeLightSource(pof,i,0.5,"  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   }
   pof << "camera {" << endl;
   pof << "  up < 0, 1, 0 >" << endl;
   pof << "  right < -4/3, 0, 0 >" << endl;
   pof << "  location ";
   writePoVVector(pof,pvp.locCam[0],pvp.locCam[1],pvp.locCam[2]);
   pof << endl << "  look_at ";
   writePoVVector(pof,pvp.lookAtCam[0],pvp.lookAtCam[1],pvp.lookAtCam[2]);
   pof << endl << "  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >";
   pof << endl << "}" << endl;
   if (bn.spaceFillingMode) {
      pof << indTabsStr(pvp.currIndLev) << "merge{" << endl;
      pvp.currIndLev++;
   }
   writeScrCharLine(pof,'/');
   pof << "#if(DrawAtomTranspSpheres)" << endl;
   putNuclei(pof,bn);
   pof << "#end\n//end if DrawAtomTranspSpheres" << endl;
   writeScrCharLine(pof,'/');
   pof << "#if(DrawStandardBonds)" << endl;
   putBonds(pof,bn);
   pof << "#end\n//end if DrawStandardBonds" << endl;
   writeScrCharLine(pof,'/');
   if (iknowacps) {
      //int indacp;
      pof << "#if(DrawAttractorCriticalPoints)" << endl;
      //solreal acprad;
      for (int i=0; i<nACP; i++) {
         //acprad=bn.drawAtSize*ATOMCRITICALPOINTSIZEFACTOR;
         //writePOVSphere(pof,0,RACP[i][0],RACP[i][1],RACP[i][2],acprad,
         //               0.0e0,0.0e0,0.0e0);
         writePOVSphere(pof,0,RACP[i][0],RACP[i][1],RACP[i][2],"RadiusACP",
                        "ColorACP");
      }
      pof << "#end\n//end if DrawAttractorCriticalPoints" << endl;
      writeScrCharLine(pof,'/');
   }
   if (iknowbcps) {
      //int indbcp;
      pof << "#if(DrawBondCriticalPoints)" << endl;
      //solreal bcprad;
      for (int i=0; i<nBCP; i++) {
         //bcprad=bn.drawAtSize*ATOMCRITICALPOINTSIZEFACTOR;
         writePOVSphere(pof,0,RBCP[i][0],RBCP[i][1],RBCP[i][2],"RadiusBCP",
                        "ColorBCP");
      }
      pof << "#end\n//end if DrawBondCriticalPoints" << endl;
      writeScrCharLine(pof,'/');
   }
   if (iknowrcps) {
      //int indrcp;
      pof << "#if(DrawRingCriticalPoints)" << endl;
      solreal rcprad;
      for (int i=0; i<nRCP; i++) {
         rcprad=bn.drawAtSize*ATOMCRITICALPOINTSIZEFACTOR;
         writePOVSphere(pof,0,RRCP[i][0],RRCP[i][1],RRCP[i][2],"RadiusRCP",
                        "ColorRCP");
      }
      pof << "#end\n//end if DrawRingCriticalPoints" << endl;
      writeScrCharLine(pof,'/');
   }
   if (iknowccps) {
      //int indccp;
      pof << "#if(DrawCageCriticalPoints)" << endl;
      solreal ccprad;
      for (int i=0; i<nCCP; i++) {
         ccprad=bn.drawAtSize*ATOMCRITICALPOINTSIZEFACTOR;
         writePOVSphere(pof,0,RCCP[i][0],RCCP[i][1],RCCP[i][2],"RadiusCCP",
                        "ColorCCP");
      }
      pof << "#end\n//end if DrawCageCriticalPoints" << endl;
      writeScrCharLine(pof,'/');
   }
   if (iknowbgps) {
      solreal gprad=0.05;
      int npts;
      pof << "#if(DrawGradientPathSpheres)" << endl;
      pof << "union {" << endl;
      for (int i=0; i<nBCP; i++) {
         npts=atBCP[i][2];
         writePOVSphere(pof,1,RBGP[i][0][0],RBGP[i][0][1],RBGP[i][0][2], \
                        gprad,"ColorABGradPath");
         for (int j=1; j<npts; j++) {
            writePOVSphere(pof,1,RBGP[i][j][0],RBGP[i][j][1],RBGP[i][j][2], \
                           gprad,"ColorABGradPath");
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathSpheres" << endl;
   }
   if (iknowbgps) {
      solreal gprad=0.05;
      int npts;
      pof << "#if(DrawGradientPathTubes)" << endl;
      pof << "union {" << endl;
      for (int i=0; i<nBCP; i++) {
         npts=atBCP[i][2];
         writePOVSphere(pof,1,RBGP[i][0][0],RBGP[i][0][1],RBGP[i][0][2], \
                        gprad,"ColorABGradPath");
         for (int j=1; j<npts; j++) {
            writePOVSphere(pof,1,RBGP[i][j][0],RBGP[i][j][1],RBGP[i][j][2], \
                           gprad,"ColorABGradPath");
            //if (tubeBGPStyle) {
               writePOVCylinder(pof,1,
                                RBGP[i][j][0],RBGP[i][j][1],RBGP[i][j][2], \
                                RBGP[i][j-1][0],RBGP[i][j-1][1],RBGP[i][j-1][2], \
                                gprad,"ColorABGradPath");
            //}
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathTubes" << endl;
   }
   pof.close();
   return true;
}
/* ************************************************************************************* */
void critPtNetWork::putBonds(ofstream &pof,bondNetWork &bn)
{
   pof << "union{" << endl;
   int k=0,atni,atnk;
   solreal startpt[3],frak1;
   for (int i=0; i<bn.nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bn.bNet[i][j];
         atni=bn.atNum[i];
         atnk=bn.atNum[k];
         //frak1=atomicRadius[atni]/(atomicRadius[atni]+atomicRadius[atnk]);
         frak1=getAtomicVDWRadius(atni)/(getAtomicVDWRadius(atni)+getAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=bn.R[i][l]*(1.0e0-frak1)+bn.R[k][l]*frak1;
         }
         if (k>0) {
            writePOVCylinder(pof,1,
                             bn.R[i][0],bn.R[i][1],bn.R[i][2],
                             startpt[0],startpt[1],startpt[2],
                             bn.drawStickSize*ATOMCRITICALPOINTSIZEFACTOR,
                             getAtomicRColorReal(atni),getAtomicGColorReal(atni),
                             getAtomicBColorReal(atni));
            writePOVCylinder(pof,1,
                             startpt[0],startpt[1],startpt[2],
                             bn.R[k][0],bn.R[k][1],bn.R[k][2],
                             bn.drawStickSize*ATOMCRITICALPOINTSIZEFACTOR,
                             getAtomicRColorReal(atnk),getAtomicGColorReal(atnk),
                             getAtomicBColorReal(atnk));
         }
      }
      writePOVSphere(pof,0,bn.R[i][0],bn.R[i][1],bn.R[i][2],
                     bn.drawStickSize*ATOMCRITICALPOINTSIZEFACTOR,
                     getAtomicRColorReal(atni),getAtomicGColorReal(atni),
                     getAtomicBColorReal(atni));
   }
   pof << "}" << endl;
   return;
}
/* ********************************************************************************* */
void critPtNetWork::putNuclei(ofstream & pof,bondNetWork &bn)
{
   int atomn;
   solreal atrad;
   for (int i=0; i<bn.nNuc; i++) {
      atomn=bn.atNum[i];
      atrad=bn.drawAtSize;
      writePOVTransparentSphere(pof,0,bn.R[i][0],bn.R[i][1],bn.R[i][2],atrad,
                                getAtomicRColorReal(atomn),getAtomicGColorReal(atomn),
                                getAtomicBColorReal(atomn),0.7);
   }
   return;
}
/* ********************************************************************************* */
void critPtNetWork::centerMolecule(bondNetWork &bn)
{
   solreal trn[3];
   for (int i=0; i<3; i++) {
      trn[i]=0.5e0*(bn.rmax[i]+bn.rmin[i]);
   }
   if (iknowacps) {
      for (int i=0; i<nACP; i++) {
         for (int j=0; j<3; j++) {
            RACP[i][j]-=trn[j];
         }
      }
   }
   if (iknowbcps) {
      for (int i=0; i<nBCP; i++) {
         for (int j=0; j<3; j++) {
            RBCP[i][j]-=trn[j];
         }
      }
   }
   if (iknowrcps) {
      for (int i=0; i<nRCP; i++) {
         for (int j=0; j<3; j++) {
            RRCP[i][j]-=trn[j];
         }
      }
   }
   if (iknowccps) {
      for (int i=0; i<nCCP; i++) {
         for (int j=0; j<3; j++) {
            RCCP[i][j]-=trn[j];
         }
      }
   }
   if (iknowbgps) {
      for (int i=0; i<nBCP; i++) {
         for (int j=0; j<atBCP[i][2]; j++) {
            for (int k=0; k<3; k++) {
               RBGP[i][j][k]-=trn[k];
            }
         }
      }
   }
   bn.centerMolecule();
   for (int i=0; i<3; i++) {centMolecVec[i]=trn[i];}
   return;
}
/* ************************************************************************************* */
void critPtNetWork::writeCPProps(string &ofnam,string &wfnnam,gaussWaveFunc &wf)
{
   ofstream ofil;
   ofil.open(ofnam.c_str());
   solreal x,y,z,gx,gy,gz,rho,hxx,hyy,hzz,hxy,hxz,hyz;
   string cptp;
   switch (mycptype) {
      case DENS:
         cptp="Density ";
         break;
      case LOLD:
         cptp="LOL ";
      default:
         break;
   }
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   ofil << "Wave function from file: " << wfnnam <<endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   ofil << cptp << "Attractor Critical Points Information" << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   if (iknowacps) {
      for (int i=0; i<80; i++) {ofil << '-';}
      ofil << endl;
      for (int i=0; i<nACP; i++) {
         ofil << "ACP(" << (i+1) << ") [" << lblACP[i] << "]:" << endl;
         x=RACP[i][0];
         y=RACP[i][1];
         z=RACP[i][2];
         writeAllFieldProperties(ofil,x,y,z,wf);
         for (int i=0; i<80; i++) {ofil << '-';}
         ofil << endl;
      }
   } else {
      ofil << "No ACP information available!" << endl;
   }
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   ofil << cptp << "Bond Critical Points Information" << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   if (iknowbcps) {
      for (int i=0; i<80; i++) {ofil << '-';}
      ofil << endl;
      for (int i=0; i<nBCP; i++) {
         ofil << "BCP(" << (i+1) << ") [" << lblBCP[i] << "]:" << endl;
         x=RBCP[i][0];
         y=RBCP[i][1];
         z=RBCP[i][2];
         writeAllFieldProperties(ofil,x,y,z,wf);
         for (int i=0; i<80; i++) {ofil << '-';}
         ofil << endl;
      }
   } else {
      ofil << "No BCP information available!" << endl;
   }
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   ofil << cptp << "Ring Critical Points Information" << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   if (iknowrcps) {
      for (int i=0; i<80; i++) {ofil << '-';}
      ofil << endl;
      if (nRCP==0) {
         ofil << "No Ring Critical Points have been found." << endl;
         for (int i=0; i<80; i++) {ofil << '-';}
         ofil << endl;
      }
      for (int i=0; i<nRCP; i++) {
         ofil << "RCP(" << (i+1) << ") [" << lblRCP[i] << "]:" << endl;
         x=RRCP[i][0];
         y=RRCP[i][1];
         z=RRCP[i][2];
         writeAllFieldProperties(ofil,x,y,z,wf);
         for (int i=0; i<80; i++) {ofil << '-';}
         ofil << endl;
      }
   } else {
      ofil << "No RCP information available!" << endl;
   }
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   ofil << cptp << "Cage Critical Points Information" << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   if (iknowccps) {
      for (int i=0; i<80; i++) {ofil << '-';}
      ofil << endl;
      if (nCCP==0) {
         ofil << "No Cage Critical Points have been found." << endl;
         for (int i=0; i<80; i++) {ofil << '-';}
         ofil << endl;
      }
      for (int i=0; i<nCCP; i++) {
         ofil << "CCP(" << (i+1) << ") [" << lblCCP[i] << "]:" << endl;
         x=RCCP[i][0];
         y=RCCP[i][1];
         z=RCCP[i][2];
         writeAllFieldProperties(ofil,x,y,z,wf);
         for (int i=0; i<80; i++) {ofil << '-';}
         ofil << endl;
      }
   } else {
      ofil << "No CCP information available!" << endl;
   }
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   ofil << cptp << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   /*
   if (iknowbgps) {
      for (int i=0; i<80; i++) {ofil << '*';}
      ofil << endl;
      ofil << cptp << " Gradient Path information." << endl;
      for (int i=0; i<80; i++) {ofil << '*';}
      ofil << endl;
      for (int i=0; i<nBCP; i++) {
         ofil << "#" << cptp << " Bond Gradient Path between atoms " << ((atBCP[i][0])+1)
              << " and " << ((atBCP[i][1])+1) << " (" << lblBCP[i] << "). Number of points:" << endl;
         ofil << "% " << atBCP[i][2] << endl;
         for (int j=0; j<atBCP[i][2]; j++) {
            for (int k=0; k<3; k++) {ofil << RBGP[i][j][k] << " ";}
            ofil << endl;
         }
      }
   }
   // */
   ofil.close();
   return;
}
/* ************************************************************************************* */
void critPtNetWork::drawNuclei(bool dn)
{
   drawNuc=dn;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::drawBonds(bool db)
{
   drawBnd=db;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::drawBondGradPaths(bool dbg)
{
   drawBGPs=dbg;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::tubeStyleBGP(bool stl)
{
   tubeBGPStyle=stl;
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getNextPointInGradientPathSimple(solreal (&xn)[3],solreal &stepsize,solreal &mgg,gaussWaveFunc &wf)
{
   static solreal rho,g[3],magg;
   wf.evalRhoGradRho(xn[0],xn[1],xn[2],rho,g);
   magg=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   for (int i=0; i<3; i++) {
      xn[i]+=(stepsize*(g[i]/magg));
   }
   mgg=magg;
   //cout << "mgg: " << mgg << endl;
   return;
}
/* ************************************************************************************* */
int critPtNetWork::findGradientPathSimple(int iacp1,int iacp2,int ibcp,gaussWaveFunc &wf)
{
   static solreal ro[3],rn[3],rho,g[3],h[3][3],eive[3][3],eival[3],hstep,dist,maggrad;
   for (int i=0; i<3; i++) {
      ro[i]=RBCP[ibcp][i];
   }
   wf.evalHessian(ro[0],ro[1],ro[2],rho,g,h);
   eigen_decomposition3(h,eive,eival);
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {
      maggrad+=(eive[i][2]*eive[i][2]);
   }
   maggrad=sqrt(maggrad);
   hstep=DEFAULTHGRADIENTPATHS;
   for (int i=0; i<3; i++) {rn[i]=ro[i]+hstep*eive[i][2]/maggrad;}
   wf.evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   //cout << "maggrad: " << maggrad << endl;
   for (int i=0; i<3; i++) {
      RGP[0][i]=rn[i];
   }
   bool iminacp;
   iminacp=false;
   int count=1;
   int maxit=ARRAYSIZEGRADPATH/2;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)&&(maggrad<MAXGRADMAG)) {
      getNextPointInGradientPathSimple(rn,hstep,maggrad,wf);
      for (int i=0; i<3; i++) {RGP[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp1][i])*(rn[i]-RACP[iacp1][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp2][i])*(rn[i]-RACP[iacp2][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      //cout << "rn(" << count << "): " << rn[0] << " " << rn[1] << " " << rn[2] << endl;
      count++;
   }
   invertOrderBGPPoints(count);
   for (int i=0; i<3; i++) {rn[i]=ro[i]-hstep*eive[i][2];}
   wf.evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   //cout << "maggrad: " << maggrad << endl;
   for (int i=0; i<3; i++) {
      RGP[count][i]=rn[i];
   }
   //int cbase=count;
   maxit+=count;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)&&(maggrad<MAXGRADMAG)) {
      getNextPointInGradientPathSimple(rn,hstep,maggrad,wf);
      for (int i=0; i<3; i++) {RGP[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp1][i])*(rn[i]-RACP[iacp1][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp2][i])*(rn[i]-RACP[iacp2][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      //cout << "rn(" << count << "): " << rn[0] << " " << rn[1] << " " << rn[2] << endl;
      count++;
   }
   return count-1;
   //wf.displayAllFieldProperties(rn[0],rn[1],rn[2]);
}
/* ************************************************************************************* */
int critPtNetWork::findGradientPathRK5(int iacp1,int iacp2,int ibcp,gaussWaveFunc &wf)
{
   static solreal ro[3],rn[3],rho,g[3],h[3][3],eive[3][3],eival[3],hstep,dist,maggrad;
   for (int i=0; i<3; i++) {
      ro[i]=RBCP[ibcp][i];
   }
   wf.evalHessian(ro[0],ro[1],ro[2],rho,g,h);
   eigen_decomposition3(h,eive,eival);
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {
      maggrad+=(eive[i][2]*eive[i][2]);
   }
   maggrad=sqrt(maggrad);
   hstep=DEFAULTHGRADIENTPATHS;
   for (int i=0; i<3; i++) {rn[i]=ro[i]+hstep*eive[i][2]/maggrad;}
   wf.evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   //cout << "maggrad: " << maggrad << endl;
   for (int i=0; i<3; i++) {
      RGP[0][i]=rn[i];
   }
   bool iminacp;
   iminacp=false;
   int count=1;
   int maxit=ARRAYSIZEGRADPATH/2;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)&&(maggrad<MAXGRADMAG)) {
      getNextPointInGradientPathRK5(rn,hstep,maggrad,wf);
      for (int i=0; i<3; i++) {RGP[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp1][i])*(rn[i]-RACP[iacp1][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp2][i])*(rn[i]-RACP[iacp2][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      //cout << "rn(" << count << "): " << rn[0] << " " << rn[1] << " " << rn[2] << endl;
      count++;
   }
   invertOrderBGPPoints(count);
   for (int i=0; i<3; i++) {rn[i]=ro[i]-hstep*eive[i][2];}
   wf.evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   //cout << "maggrad: " << maggrad << endl;
   for (int i=0; i<3; i++) {
      RGP[count][i]=rn[i];
   }
   //int cbase=count;
   maxit+=count;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)&&(maggrad<MAXGRADMAG)) {
      getNextPointInGradientPathRK5(rn,hstep,maggrad,wf);
      for (int i=0; i<3; i++) {RGP[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp1][i])*(rn[i]-RACP[iacp1][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-RACP[iacp2][i])*(rn[i]-RACP[iacp2][i]));}
      dist=sqrt(dist);
      if (dist<=0.01) {iminacp=true;}
      //cout << "rn(" << count << "): " << rn[0] << " " << rn[1] << " " << rn[2] << endl;
      count++;
   }
   return count-1;
   //wf.displayAllFieldProperties(rn[0],rn[1],rn[2]);
}
/* ************************************************************************************* */
void critPtNetWork::addGradientPathToPOVFile(int npts,ofstream &ofil)
{
   //int indgp;
   solreal gprad;
   for (int i=0; i<npts; i++) {
      gprad=0.05;
      writePOVSphere(ofil,0,RGP[i][0],RGP[i][1],RGP[i][2],gprad,0.0e0,1.0e0,0.0e0);
   }
}
/* ************************************************************************************* */
void critPtNetWork::setBondPaths(gaussWaveFunc &wf)
{
   if (!iknowbcps) {
      displayErrorMessage("Please look first for the BCPs...\nNothing to be done!");
      return;
   }
   int npts;
   alloc3DRealArray(string("RBGP"),dBCP,ARRAYSIZEGRADPATH,3,RBGP);
   cout << "Calculating Bond Gradient Paths..." << endl;
   //cout << "nBCP: " << nBCP << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal hstep,rseed[3];
   hstep=DEFAULTHGRADIENTPATHS;
   int arrsize=ARRAYSIZEGRADPATH;
   int at1,at2;
   for (int i=0; i<nBCP; i++) {
      for (int k=0; k<3; k++) {rseed[k]=RBCP[i][k];}
      at1=atBCP[i][0];
      at2=atBCP[i][1];
      //cout << "BCP[" << (i+1) << "]: " << "at1: " << at1 << "(" << wf.atLbl[at1] << "), at2: "\
           << at2 << "(" << wf.atLbl[at2] << ")" << endl;
      //npts=findGradientPathRK5(atBCP[i][0],atBCP[i][1],i,wf);
      npts=findSingleRhoGradientPathRK5(at1,at2,hstep,arrsize,RBGP[i],rseed,wf);
      atBCP[i][2]=npts;
      if (npts>0) {nBGP++;}
      //for (int j=0; j<npts; j++) {
      //   for (int k=0; k<3; k++) {RBGP[i][j][k]=RGP[j][k];}
      //}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nBCP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   if (nBGP!=nBCP) {displayWarningMessage("For some unknown reason nBGP!=nBCP...");}
   iknowbgps=true;
}
/* ************************************************************************************* */
void critPtNetWork::getNextPointInGradientPathRK5(solreal (&xn)[3],
                                                  solreal &stepsize,solreal &mgg,gaussWaveFunc &wf)
{
   static const solreal b[15]={0.2e0, 3.0e0/40.0e0, 9.0e0/40.0e0, 0.3e0, -0.9e0, 1.2e0, \
      -11.0e0/54.0e0, 2.5e0, -70.0e0/27.0e0, 35.0e0/27.0e0, 1631.0e0/55296.0e0, \
      175.0e0/512.0e0, 575.0e0/13824.0e0, 44275.0e0/110592.0e0, 253.0e0/4096.0e0};
	static const solreal c1=37.0e0/378.0e0;
	static const solreal c3=250.0e0/621.0e0;
	static const solreal c4=125.0e0/594.0e0;
	static const solreal c6=512.0e0/1771.0e0;
	solreal k[6][3],rho,g[3],xt[3],maggrad;
	int offset;
	for(int i=0; i<6; i++) {
		offset=((i*(i-1))>>1);
		for(int j=0; j<3; j++) {
			xt[j] = xn[j];
			for(int l=0; l<i; l++) {
				xt[j] += b[l+offset]*k[l][j];
         }
		}
		wf.evalRhoGradRho(xt[0],xt[1],xt[2],rho,g);
      maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
		for(int l=0; l<3; l++) {
         k[i][l] = stepsize*g[l]/maggrad;
      }
	}
	for(int i=0; i<3; i++) {
		xn[i]+=(c1*k[0][i]+c3*k[2][i]+c4*k[3][i]+c6*k[5][i]);
   }
   wf.evalRhoGradRho(xn[0],xn[1],xn[2],rho,g);
   mgg=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::invertOrderBGPPoints(int dim)
{
   int numop=((dim+1)>>1);
   
   solreal tmp;
   for (int i=0; i<numop; i++) {
      for (int j=0; j<3; j++) {
         tmp=RGP[i][j];
         RGP[i][j]=RGP[dim-i-1][j];
         RGP[dim-i-1][j]=tmp;
      }
   }
}
/* ************************************************************************************* */
void critPtNetWork::invertOrderBGPPoints(int dim,solreal** (&arr))
{
#if DEBUG
   if (arr==NULL) {
      displayErrorMessage("The pointer for this array is null!!!");
      DISPLAYDEBUGINFOFILELINE;
      return;
   }
#endif
   int numop=((dim+1)>>1);
   solreal tmp;
   for (int i=0; i<numop; i++) {
      for (int j=0; j<3; j++) {
         tmp=arr[i][j];
         arr[i][j]=arr[dim-i-1][j];
         arr[dim-i-1][j]=tmp;
      }
   }
}
/* ************************************************************************************* */
bool critPtNetWork::seekSingleRhoBCP(int ata,int atb,gaussWaveFunc &wf,solreal (&x)[3])
{
   
   int idxa,idxb;
   idxa=3*ata;
   idxb=3*atb;
   //cout << "ata: " << ata << ", atb: " << atb << endl;
   for (int j=0; j<3; j++) {x[j]=0.5e0*(wf.R[idxa+j]+wf.R[idxb+j]);}
   seekRhoBCP(x,wf);
   solreal rho,g[3];
   wf.evalRhoGradRho(x[0],x[1],x[2],rho,g);
   solreal magg=0.0e0;
   for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
   magg=sqrt(magg);
   if (!((rho>MINRHOSIGNIFICATIVEVAL)&&(magg<EPSGRADMAG))) {
      displayErrorMessage(string("The chosen atoms ("+wf.atLbl[ata]\
               +","+wf.atLbl[atb]+") do not lead to a BCP!"));
      //wf.displayAllFieldProperties(x[0],x[1],x[2]);
      return false;
   }
   return true;
}
/* ************************************************************************************* */
int critPtNetWork::findSingleRhoGradientPathRK5(int at1,int at2,solreal hstep,
                                                int dima,solreal** (&arbgp),solreal (&ro)[3],
                                                gaussWaveFunc &wf)
{
   solreal rn[3],rho,g[3],h[3][3],eive[3][3],eival[3],dist,maggrad;
   seekSingleRhoBCP(at1,at2,wf,ro);
   //cout << "ro: " << ro[0] << " " << ro[1] << " " << ro[2] << endl;
   wf.evalHessian(ro[0],ro[1],ro[2],rho,g,h);
   eigen_decomposition3(h,eive,eival);
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {maggrad+=(eive[i][2]*eive[i][2]);}
   maggrad=sqrt(maggrad);
   //hstep=DEFAULTHGRADIENTPATHS;
   for (int i=0; i<3; i++) {rn[i]=ro[i]+hstep*eive[i][2]/maggrad;}
   wf.evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   //cout << "maggrad: " << maggrad << endl;
   for (int i=0; i<3; i++) {
      arbgp[0][i]=ro[i];
      arbgp[1][i]=rn[i];
   }
   bool iminacp;
   iminacp=false;
   int count=2;
   int maxit=dima/2;
   int iacp1,iacp2;
   iacp1=3*at1;
   iacp2=3*at2;
   solreal loopdist;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)) {
      getNextPointInGradientPathRK5(rn,hstep,maggrad,wf);
      for (int i=0; i<3; i++) {arbgp[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf.R[iacp1+i])*(rn[i]-wf.R[iacp1+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {
         iminacp=true;
         count++;
         for (int i=0; i<3; i++) {arbgp[count][i]=wf.R[iacp1+i];}
      }
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf.R[iacp2+i])*(rn[i]-wf.R[iacp2+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {
         iminacp=true;
         count++;
         for (int i=0; i<3; i++) {arbgp[count][i]=wf.R[iacp2+i];}
      }
      if (count>5) {
         loopdist=0.0e0;
         for (int i=0; i<3; i++) {loopdist+=((rn[i]-arbgp[count-2][i])*(rn[i]-arbgp[count-2][i]));}
         loopdist=sqrt(loopdist);
         if ((5.0e0*loopdist)<hstep) {
            iminacp=true;
            count--;
         }
      }
      //cout << "rn(" << count << "): " << rn[0] << " " << rn[1] << " " << rn[2] << endl;
      count++;
      if (count==dima) {
         displayWarningMessage("You need a bigger array for storing the BGP coordinates!");
         return count-1;
      }
   }
   invertOrderBGPPoints(count,arbgp);
   //cout << "Checkpoint... count: " << count  << endl;
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {maggrad+=(eive[i][2]*eive[i][2]);}
   maggrad=sqrt(maggrad);
   for (int i=0; i<3; i++) {rn[i]=ro[i]-hstep*eive[i][2]/maggrad;}
   wf.evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   //cout << "maggrad: " << maggrad << endl;
   for (int i=0; i<3; i++) {
      arbgp[count][i]=rn[i];
   }
   count++;
   //int cbase=count;
   maxit+=count;
   iminacp=false;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)) {
      getNextPointInGradientPathRK5(rn,hstep,maggrad,wf);
      for (int i=0; i<3; i++) {arbgp[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf.R[iacp1+i])*(rn[i]-wf.R[iacp1+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {iminacp=true;}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf.R[iacp2+i])*(rn[i]-wf.R[iacp2+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {iminacp=true;}
      if (count>5) {
         loopdist=0.0e0;
         for (int i=0; i<3; i++) {loopdist+=((rn[i]-arbgp[count-2][i])*(rn[i]-arbgp[count-2][i]));}
         loopdist=sqrt(loopdist);
         if ((5.0e0*loopdist)<hstep) {
            iminacp=true;
            count--;
         }
      }
      //cout << "rn(" << count << "): " << rn[0] << " " << rn[1] << " " << rn[2] << endl;
      count++;
      if (count==dima) {
         displayWarningMessage("You need a bigger array for storing the BGP coordinates!");
         return count-1;
      }
   }
   //count++;
   if (count==dima) {
      displayWarningMessage("You need a bigger array for storing the BGP coordinates!");
      return count-1;
   }
   dist=0.0e0;
   for (int i=0; i<3; i++) {dist+=((arbgp[count-1][i]-wf.R[iacp1+i])*(arbgp[count-1][i]-wf.R[iacp1+i]));}
   dist=sqrt(dist);
   if (dist<=hstep) {
      for (int i=0; i<3; i++) {arbgp[count][i]=wf.R[iacp1+i];}
      count++;
   }
   dist=0.0e0;
   for (int i=0; i<3; i++) {dist+=((arbgp[count-1][i]-wf.R[iacp2+i])*(arbgp[count-1][i]-wf.R[iacp2+i]));}
   dist=sqrt(dist);
   if (dist<=hstep) {
      for (int i=0; i<3; i++) {arbgp[count][i]=wf.R[iacp2+i];}
      count++;
   }
   //cout << "Checkpoint2... count: " << count  << endl;
   return count;
   //wf.displayAllFieldProperties(rn[0],rn[1],rn[2]);
}
/* ************************************************************************************* */
void critPtNetWork::findTwoClosestAtoms(solreal (&xo)[3],gaussWaveFunc &wf,int &idx1st,\
      int &idx2nd)
{
   if ( wf.nNuc<2 ) {idx1st=0; idx2nd=0; return;}
   solreal xmagt=0.0e0,xmag1=0.0e0,xmag2=0.0e0;
   int ii1,ii2,iit;
   for ( int k=0 ; k<3 ; k++ ) {xmag1+=((xo[k]-wf.getR(0,k))*(xo[k]-wf.getR(0,k)));}
   ii1=0;
   for ( int k=0 ; k<3 ; k++ ) {xmag2+=((xo[k]-wf.getR(1,k))*(xo[k]-wf.getR(1,k)));}
   ii2=1;
   if ( xmag1>xmag2 ) {
      xmagt=xmag1; xmag1=xmag2; xmag2=xmagt;
      ii1=1;
      ii2=0;
   }
   for ( int i=2 ; i<wf.nNuc ; i++ ) {
      xmagt=0.0e0;
      for ( int k=0 ; k<3 ; k++ ) {xmagt+=((xo[k]-wf.getR(i,k))*(xo[k]-wf.getR(i,k)));}
      if ( xmagt<xmag2 ) {xmag2=xmagt; ii2=i;}
      if ( xmag1>xmag2 ) {
         xmagt=xmag1; xmag1=xmag2; xmag2=xmagt;
         iit=ii1;     ii1=ii2;     ii2=iit;
      }
   }
   idx1st=ii1;
   idx2nd=ii2;
#if DEBUG
      if ( idx1st==idx2nd ) {
         displayWarningMessage("Identical atoms!");
         DISPLAYDEBUGINFOFILELINE;
      }
#endif /* ( DEBUG ) */
}
/* ************************************************************************************* */
bool critPtNetWork::iKnowACPs(void)
{
   return iknowacps;
}
/* ************************************************************************************* */
bool critPtNetWork::iKnowBCPs(void)
{
   return iknowbcps;
}
/* ************************************************************************************* */
bool critPtNetWork::iKnowRCPs(void)
{
   return iknowrcps;
}
/* ************************************************************************************* */
bool critPtNetWork::iKnowCCPs(void)
{
   return iknowccps;
}
/* ************************************************************************************* */
ScalarFieldType critPtNetWork::myCPType()
{
   return mycptype;
}
/* ************************************************************************************* */
bool critPtNetWork::iKnowBGPs(void)
{
   return iknowbgps;
}
/* ************************************************************************************* */
bool critPtNetWork::readFromFile(string inname)
{
   string tmps;
   tmps=inname.substr(inname.length()-3,inname.length());
   if (tmps!="cpx") {
      displayErrorMessage("Not a valid file!");
      return false;
   }
   ifstream cfil;
   cfil.open(inname.c_str(),ios::in);
   if (!(cfil.good())) {
      displayErrorMessage("This file could not be opened...");
      return false;
   }
   cout << "Loading Critical Points information from file:\n   '" << inname << "'..." << endl;
   mycptype=cpxGetCriticalPointFieldType(cfil);
   nACP=cpxGetNOfACPs(cfil);
   switch (mycptype) {
      case DENS:
         dACP=nACP*MAXRHOACPSPERATOM;
         break;
      case LOLD:
         dACP=nACP*MAXLOLACPSPERATOM;
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if (dACP<0) {
      displayWarningMessage(string("ACPs are not given in file "+inname));
      iknowacps=false;
   } else {
      alloc2DRealArray(string("RACP"),dACP,3,RACP,1.0e+50);
      alloc1DStringArray("lblACP",dACP,lblACP);
      if (nACP>=0) {
         cpxGetACPCartCoordFromFile(cfil,nACP,RACP);
         cpxGetACPLabelsFromFile(cfil,nACP,lblACP);
         iknowacps=true;
      } else {
         iknowacps=false;
      }
   }
   if (iknowacps) {
      nBCP=cpxGetNOfBCPs(cfil);
      dBCP=(nACP*(nACP-1))/2;
      alloc2DRealArray(string("RBCP"),dBCP,3,RBCP,1.0e+50);
      alloc2DIntArray(string("atBCP"),dBCP,3,atBCP,-1);
      alloc1DStringArray("lblBCP",dBCP,lblBCP);
      if (nBCP>=0) {
         cpxGetBCPCartCoordFromFile(cfil,nBCP,RBCP);
         cpxGetBCPConnectivityFromFile(cfil,nBCP,atBCP);
         cpxGetBCPLabelsFromFile(cfil,nBCP,lblBCP);
         iknowbcps=true;
      } else {
         iknowbcps=false;
      }
   } else {
      displayErrorMessage("First look for ACPs...\nQuiting...");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
   }
   if (iknowbcps) {
      nRCP=cpxGetNOfRCPs(cfil);
      dRCP=(nBCP*(nBCP-1))/2;
      alloc2DRealArray(string("RRCP"),dRCP,3,RRCP,1.0e+50);
      alloc1DStringArray("lblRCP",dRCP,lblRCP);
      if (nRCP>=0) {
         cpxGetRCPCartCoordFromFile(cfil,nRCP,RRCP);
         cpxGetRCPLabelsFromFile(cfil,nRCP,lblRCP);
         iknowrcps=true;
      } else {
         iknowrcps=false;
      }
   }
   if (iknowrcps) {
      nCCP=cpxGetNOfCCPs(cfil);
      dCCP=(nRCP*(nRCP-1))/2;
      alloc2DRealArray(string("RCCP"),dCCP,3,RCCP,1.0e+50);
      alloc1DStringArray("lblCCP",dCCP,lblCCP);
      if (nCCP>=0) {
         cpxGetCCPCartCoordFromFile(cfil,nCCP,RCCP);
         cpxGetCCPLabelsFromFile(cfil,nCCP,lblCCP);
         iknowccps=true;
      } else {
         iknowccps=false;
      }
   }
   iknowallcps=(iknowacps&&iknowbcps&&iknowrcps&&iknowccps);
   nBGP=cpxGetNOfBondPaths(cfil);
   if (dBCP>0) {alloc3DRealArray(string("RBGP"),dBCP,ARRAYSIZEGRADPATH,3,RBGP);}
   if (nBGP>=0&&atBCP!=NULL) {
      cpxGetNOfPtsPerBondPath(cfil,nBGP,atBCP);
      cpxGetBondPathData(cfil,nBGP,atBCP,RBGP);
      iknowbgps=true;
   } else {iknowbgps=false;}
   cout << "Critical Point State Loaded!" << endl;
   displayStatus(true);
   return true;
}
/* ************************************************************************************* */
void critPtNetWork::displayStatus(bool lngdesc)
{
   printScrCharLine('+');
   if (lngdesc) {
      if (atBCP==NULL) {displayWarningMessage("atBCP is not allocated.");}
      if (RACP==NULL) {displayWarningMessage("RACP is not allocated.");}
      if (RBCP==NULL) {displayWarningMessage("RBCP is not allocated.");}
      if (RRCP==NULL) {displayWarningMessage("RRCP is not allocated.");}
      if (RCCP==NULL) {displayWarningMessage("RCCP is not allocated.");}
      if (RBGP==NULL) {displayWarningMessage("RBGP is not allocated.");}
      if (lblACP==NULL) {displayWarningMessage("lblACP is not allocated.");}
      if (lblBCP==NULL) {displayWarningMessage("lblBCP is not allocated.");}
      if (lblRCP==NULL) {displayWarningMessage("lblRCP is not allocated.");}
      if (lblCCP==NULL) {displayWarningMessage("lblCCP is not allocated.");}
   }
   printScrCharLine('-');
   if (iknowacps) {cout << nACP << " ACPs found." << endl;}
   else {displayWarningMessage("ACPs search needed...");}
   if (iknowbcps) {cout << nBCP << " BCPs found." << endl;}
   else {displayWarningMessage("BCPs search needed...");}
   if (iknowrcps) {cout << nRCP << " RCPs found." << endl;}
   else {displayWarningMessage("RCPs search needed...");}
   if (iknowccps) {cout << nCCP << " CCPs found." << endl;}
   else {displayWarningMessage("CCPs search needed...");}
   if (iknowallcps) {
      printScrCharLine('-');
      cout << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
      printScrCharLine('-');
   }
   if (iknowbgps) {
      if (nBGP!=nBCP) {displayWarningMessage("For some unknown reason nBGP!=nBCP");}
      else {cout << nBGP << " Bond Gradient Paths found." << endl;}
   } else {displayWarningMessage("BGPs search needed...");}
   printScrCharLine('+');
   return;
}
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
/* ************************************************************************************* */
#endif//_CRITPTNETWORK_CPP_

