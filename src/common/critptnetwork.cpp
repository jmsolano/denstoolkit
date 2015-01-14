
/*
 *  critptnetwork.cpp
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 05/06/13.
 *  Copyright 2013-2015. All rights reserved.
 *  This is a new version, several functions and some
 *  redesign of the class is done in 03-04/Jan/2015.
 *
 */

#ifndef _CRITPTNETWORK_CPP_
#define _CRITPTNETWORK_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

#include "solscrutils.h"
#include "solfileutils.h"
#include "critptnetwork.h"
#include "wavefunctionclass.h"
#include "bondnetwork.h"
#include "eig2-4.h"
#include "iofuncts-cpx.h"
#include "atomradiicust.h"
#include "solstringtools.h"
#include "solpovtools.h"
// The first 94 atomic radii are given,
//  the rest are set to be 0.80e0
//
//#include "atomcolschcust.h" //Choose this for the palette defined by JMHP
#include "atomcolschjmol.h" //Choose this for the palette used in JMol

/* ************************************************************************************ */
//                    Macros
/* ************************************************************************************ */
#define CPNW_EPSFABSDIFFCOORD (1.0e-06)
#define CPNW_MINRHOSIGNIFICATIVEVAL (5.0e-05)
#define CPNW_MINLOLSIGNIFICATIVEVAL (5.0e-04)
#define CPNW_EXTENDEDBONDDISTFACTOR (3.5e0)
#define CPNW_BONDISTEXTCCPSEARCHFACTOR (2.0e0)
#define CPNW_ATOMCRITICALPOINTSIZEFACTOR (0.2e0)

#ifndef CPNW_EPSEIGENVALUECPSEARCH
#define CPNW_EPSEIGENVALUECPSEARCH (1.0e-08)
#endif

#ifndef CPNW_MAXSTEPSIZEACPRHOSEARCH
#define CPNW_MAXSTEPSIZEACPRHOSEARCH (0.1e0)
#endif

#ifndef CPNW_MAXSTEPSIZEBCPSEARCH
#define CPNW_MAXSTEPSIZEBCPSEARCH 0.4
#endif

#ifndef CPNW_MAXSTEPSIZERCPSEARCH
#define CPNW_MAXSTEPSIZERCPSEARCH 0.35
#endif

#ifndef CPNW_MAXSTEPSIZECCPSEARCH
#define CPNW_MAXSTEPSIZECCPSEARCH 0.3
#endif

#ifndef CPNW_MAXSTEPSIZEACPLOLSEARCH
#define CPNW_MAXSTEPSIZEACPLOLSEARCH 0.01
#endif

#ifndef CPNW_DEFAULTGRADIENTPATHS
#define CPNW_DEFAULTGRADIENTPATHS 0.1
#endif

#ifndef CPNW_MAXITERATIONACPSEARCH
#define CPNW_MAXITERATIONACPSEARCH 20
#endif

#ifndef CPNW_MAXITERATIONBCPSEARCH
#define CPNW_MAXITERATIONBCPSEARCH 40
#endif

#ifndef CPNW_MAXITERATIONRCPSEARCH
#define CPNW_MAXITERATIONRCPSEARCH 40
#endif

#ifndef CPNW_MAXITERATIONCCPSEARCH
#define CPNW_MAXITERATIONCCPSEARCH 120
#endif

//#ifndef CPNW_CPNW_EPSRHOACPGRADMAG
//#define CPNW_CPNW_EPSRHOACPGRADMAG (1.0e-12)
//#endif

#ifndef CPNW_EPSRHOACPGRADMAG
#define CPNW_EPSRHOACPGRADMAG (5.0e-10)
#endif

#ifndef CPNW_EPSLOLACPGRADMAG
#define CPNW_EPSLOLACPGRADMAG (5.0e-04)
#endif

#ifndef CPNW_MAXRHOACPSPERATOM
#define CPNW_MAXRHOACPSPERATOM (4)
#endif

#ifndef CPNW_MAXLOLACPSPERATOM
#define CPNW_MAXLOLACPSPERATOM (40)
#endif

#ifndef CPNW_MINARRAYSIZE
#define CPNW_MINARRAYSIZE (8)
#endif

/* ************************************************************************************ */
void critPtNetWork::init()
{
   //publics
   nACP=nBCP=nRCP=nCCP=0;
   normalbcp=0;
   nBGP=0;
   atBCP=NULL;
   RACP=RBCP=RRCP=RCCP=NULL;
   RBGP=NULL;
   lblACP=lblBCP=lblRCP=lblCCP=NULL;
   for (int i=0; i<3; i++) {centMolecVec[i]=0.0e0;}
   //privates
   dACP=dBCP=dRCP=dCCP=0;
   maxItACP=CPNW_MAXITERATIONACPSEARCH;
   maxItBCP=CPNW_MAXITERATIONBCPSEARCH;
   maxItRCP=CPNW_MAXITERATIONRCPSEARCH;
   maxItCCP=CPNW_MAXITERATIONCCPSEARCH;
   stepSizeACP=CPNW_MAXSTEPSIZEACPRHOSEARCH;
   stepSizeBCP=CPNW_MAXSTEPSIZEBCPSEARCH;
   stepSizeRCP=CPNW_MAXSTEPSIZERCPSEARCH;
   stepSizeCCP=CPNW_MAXSTEPSIZECCPSEARCH;
   iknowacps=iknowbcps=iknowrcps=iknowccps=false;
   iknowallcps=false;
   iknowbgps=false;
   mycptype=NONE;
   drawNuc=false;
   drawBnd=true;
   drawBGPs=false;
   tubeBGPStyle=false;
   wf=NULL;
   bn=NULL;
}
/* ************************************************************************************ */
critPtNetWork::critPtNetWork(gaussWaveFunc &uwf,bondNetWork &ubn)
{
   init();
   wf=&uwf;
   bn=&ubn;
}
/* ************************************************************************************ */
critPtNetWork::~critPtNetWork()
{
   dealloc3DRealArray(RBGP,dBCP,CPNW_ARRAYSIZEGRADPATH);
   dealloc2DRealArray(RACP,dACP);
   dealloc1DStringArray(lblACP);
   dealloc2DRealArray(RBCP,dBCP);
   dealloc1DStringArray(lblBCP);
   dealloc2DIntArray(atBCP,dBCP);
   dealloc2DRealArray(RRCP,dRCP);
   dealloc1DStringArray(lblRCP);
   dealloc2DRealArray(RCCP,dCCP);
   dealloc1DStringArray(lblCCP);
   wf=NULL;
   bn=NULL;
}
/* ************************************************************************************ */
solreal critPtNetWork::V0=0.0e0;
solreal critPtNetWork::V5=2.0e0/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
solreal critPtNetWork::V8=(1.0e0+sqrt(5.0e0))/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
solreal critPtNetWork::IHV[nIHV][3]={
   {  V0,  V0,  V0 },
   {  2.0e0*CPNW_EPSFABSDIFFCOORD,  V0,  V0 },
   {  V0,  2.0e0*CPNW_EPSFABSDIFFCOORD,  V0 },
   {  V0,  V0,  2.0e0*CPNW_EPSFABSDIFFCOORD },
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
/* ************************************************************************************ */
void critPtNetWork::setCriticalPoints(ScalarFieldType ft)
{
   mycptype=ft;
   if (!bn->imstp()) {
      displayErrorMessage("Trying to use a non set up bond network object!");
      return;
   }
   switch (ft) {
      case DENS:
         cout << "Scanning for Density Critical Points." << endl;
         dACP=(bn->nNuc)*CPNW_MAXRHOACPSPERATOM;
         break;
      case LOLD:
         cout << "Scanning for LOL Critical Points." << endl;
         dACP=(bn->nNuc)*CPNW_MAXLOLACPSPERATOM;
         stepSizeACP=CPNW_MAXSTEPSIZEACPLOLSEARCH;
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   alloc2DRealArray(string("RACP"),dACP,3,RACP,1.0e+50);
   alloc1DStringArray("lblACP",dACP,lblACP);
   cout << "Looking for Attractor Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowacps=setRhoACPs();
         break;
      case LOLD:
         iknowacps=setLOLACPs();
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if (iknowacps) {
      dBCP=(nACP*(nACP-1))/2;
      if ( dBCP<CPNW_MINARRAYSIZE ) {dBCP=CPNW_MINARRAYSIZE;}
      alloc2DRealArray(string("RBCP"),dBCP,3,RBCP,1.0e+50);
      alloc2DIntArray(string("atBCP"),dBCP,3,atBCP,-1);
      alloc1DStringArray("lblBCP",dBCP,lblBCP);
   } else {
      displayErrorMessage("First look for ACPs...\n");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
      return;
   }
   cout << "Looking for Bond Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowbcps=setRhoBCPs();
         break;
      case LOLD:
         iknowbcps=setLOLBCPs();
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         //exit(1);
         break;
   }
#if DEBUG
   cout << "nBCP: " << nBCP << ", dBCP: " << dBCP << endl;
#endif
   if (iknowbcps) {
      dRCP=(nBCP*(nBCP-1))/2;
      if ( dRCP<CPNW_MINARRAYSIZE ) {dRCP=CPNW_MINARRAYSIZE;}
      alloc2DRealArray(string("RRCP"),dRCP,3,RRCP,1.0e+50);
      alloc1DStringArray("lblRCP",dRCP,lblRCP);
   }
   cout << "Looking for Ring Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowrcps=setRhoRCPs();
         break;
      case LOLD:
         iknowrcps=setLOLRCPs();
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
#if DEBUG
   cout << "nRCP: " << nRCP << ", dRCP: " << dRCP << endl;
#endif
   if (iknowrcps) {
      dCCP=(nRCP*(nRCP-1))/2;
      if ( dCCP<CPNW_MINARRAYSIZE ) {dCCP=CPNW_MINARRAYSIZE;}
      alloc2DRealArray(string("RCCP"),dCCP,3,RCCP,1.0e+50);
      alloc1DStringArray("lblCCP",dCCP,lblCCP);
   }
   cout << "Looking for Cage Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowccps=setRhoCCPs();
         break;
      case LOLD:
         iknowccps=setLOLCCPs();
         break;
      default:
         displayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
#if DEBUG
   cout << "nCCP: " << nCCP << ", dCCP: " << dCCP << endl;
#endif
   if ( mkextsearch ) {extendedSearchCPs();}
   iknowallcps=(iknowacps&&iknowbcps&&iknowrcps&&iknowccps);
   if (iknowallcps) {
      printScrCharLine('*');
      cout << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
      printScrCharLine('*');
   }
}
/* ************************************************************************************ */
void critPtNetWork::setBondPaths()
{
   if (!iknowbcps) {
      displayErrorMessage("Please look first for the BCPs...\nNothing to be done!");
      return;
   }
   int npts;
   alloc3DRealArray(string("RBGP"),dBCP,CPNW_ARRAYSIZEGRADPATH,3,RBGP);
   cout << "Calculating Bond Gradient Paths..." << endl;
   //cout << "nBCP: " << nBCP << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal hstep,rseed[3];
   hstep=CPNW_DEFAULTGRADIENTPATHS;
   int arrsize=CPNW_ARRAYSIZEGRADPATH;
   int at1,at2;
   for (int i=0; i<nBCP; i++) {
      for (int k=0; k<3; k++) {rseed[k]=RBCP[i][k];}
      at1=atBCP[i][0];
      at2=atBCP[i][1];
      npts=findSingleRhoGradientPathRK5(at1,at2,hstep,arrsize,RBGP[i],rseed);
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
/* ************************************************************************************ */
bool critPtNetWork::setRhoACPs()
{
   solreal x[3],rho,g[3],magg;
   string lbl;
   int sig;
   bool tmpbool=false;
   for (int i=0; i<(bn->nNuc); i++) {
      for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
      seekRhoACP(x,rho,g,sig);
      magg=computeMagnitudeV3(g);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         lbl=bn->atLbl[i];
         addRhoACP(x,sig,lbl);
      } else {
         for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
         lbl=bn->atLbl[i]+"+";
         addRhoACP(x,sig,lbl);
         tmpbool=true;
      }
   }
   if ( tmpbool ) {
      displayWarningMessage("Some ACPs presented |nabla rho|!= 0... Look for '+' in labels.");
   }
   if ( nACP!=(bn->nNuc) ) {
      displayWarningMessage("Number of regular ACPs is different from number of Nuclei!");
   }
   int regnACP=nACP;
   cout << "Found " << regnACP << " regular ACPs (" << bn->nNuc << " nuc)" << endl;
   cout << "Looking between atom pairs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal xs[3],rad;
   for ( int i=0 ; i<(bn->nNuc) ; i++ ) {
      for ( int j=(i+1) ; j<(bn->nNuc) ; j++ ) {
         for ( int k=0 ; k<3 ; k++ ) {xs[k]=(bn->R[i][k]-bn->R[j][k]);}
         rad=computeMagnitudeV3(xs);
         if ( rad<=(bn->maxBondDist) ) {
            for ( int k=0 ; k<3 ; k++ ) {xs[k]=0.5e0*(bn->R[i][k]+bn->R[j][k]);}
            lbl="NNACP";
            rad*=0.3e0;
            seekRhoACPsAroundAPoint(xs,rad,lbl,8);
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((bn->nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   if (regnACP==nACP) {
      cout << "No more ACPs found." << endl;
   } else {
      cout << (nACP-regnACP) << " new ACP";
      if ((nACP-regnACP)>1) {cout << "s";}
      cout << " found! Total number of ACPs: " << nACP << endl;
   }
   return true;
}
/* ************************************************************************************ */
bool critPtNetWork::setRhoBCPs(void)
{
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   int ata,atb,k,sig;
   int mypos;
   solreal x[3],g[3],rho;
   string lbl;
   for (int i=0; i<(bn->nNuc); i++) {
      ata=i;
      atb=bn->bNet[ata][0];
      k=1;
      while ((k<MAXBONDINGATOMS)&&(atb>0)) {
         for (int j=0; j<3; j++) {x[j]=0.5e0*(bn->R[ata][j]+bn->R[atb][j]);}
         seekRhoBCP(x,rho,g,sig);
         if (ata<atb) {
            lbl=bn->atLbl[ata]+string("-")+bn->atLbl[atb];
         } else {
            lbl=bn->atLbl[atb]+string("-")+bn->atLbl[ata];
         }
         if (rho>CPNW_MINRHOSIGNIFICATIVEVAL&&addRhoBCP(x,sig,lbl,mypos)&&mypos>=0) {
            atBCP[mypos][0]=ata;
            atBCP[mypos][1]=atb;
         }
         atb=bn->bNet[ata][k];
         ++k;
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
   normalbcp=nBCP;
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   /* ------------------------------------------------------ */
   cout << nBCP << " BCPs found.\n";
   string ll;
   cout << "Including non-nuclear attractors..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for ( int i=(bn->nNuc) ; i<nACP ; i++ ) {
      seekRhoBCPWithExtraACP(i,(bn->maxBondDist));
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   cout << "Looking for possible extra BCPs (Hydrogen bonds)..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal extbd=(bn->maxBondDist*CPNW_EXTENDEDBONDDISTFACTOR),dd=0.0e0;
   for (int i=0; i<nACP; i++) {
      for (int j=(i+1); j<nACP; j++) {
         for (int k=0; k<3; k++) {x[k]=((RACP[i][k]-RACP[j][k]));}
         dd=computeMagnitudeV3(x);
         if ((dd>=(bn->maxBondDist*0.9e0))&&(dd<=extbd)) {
            for (int k=0; k<3; k++) {x[k]=0.5e0*(RACP[i][k]+RACP[j][k]);}
            if ((wf->evalDensity(x[0],x[1],x[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
               seekRhoBCP(x,rho,g,sig);
               if (i<j) {ata=i; atb=j;} else {ata=j; atb=i;}
               lbl=string("*")+lblACP[ata]+string("-")+lblACP[atb];
               wf->evalRhoGradRho(x[0],x[1],x[2],rho,g);
               if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(computeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
                  addRhoBCP(x,sig,lbl,mypos);
                  if ((mypos>=0)&&(mypos<dBCP)) {
                     atBCP[mypos][0]=ata;
                     atBCP[mypos][1]=atb;
                  }
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
         findTwoClosestACPs(tbcp,ata,atb);
         atBCP[i][0]=ata;
         atBCP[i][1]=atb;
         lblBCP[i]="*"+lblACP[ata]+"-"+lblACP[atb];
      }
      cout << (nBCP-normalbcp) << " new BCP";
      if ((nBCP-normalbcp)>1) {cout << "s";}
      cout << " found! Total number of BCPs: " << nBCP << endl;
   }
   return true;
}
/* ************************************************************************************ */
bool critPtNetWork::setRhoRCPs(void)
{
   solreal x[3],rho,g[3];
   string lbl;
   int sig,mypos;
   for (int i=0; i<nBCP; i++) {
      for (int j=(i+1); j<nBCP; j++) {
         for (int k=0; k<3; k++) {x[k]=(RBCP[i][k]-RBCP[j][k]);}
         if (computeMagnitudeV3(x)>(bn->maxBondDist*2.0e0)) {continue;}
         for ( int k=0 ; k<3 ; k++ ) {x[k]=0.5e0*(RBCP[i][k]+RBCP[j][k]);}
         seekRhoRCP(x,rho,g,sig);
         lbl=(lblBCP[i]+string("-")+lblBCP[j]);
         //wf->evalRhoGradRho(x[0],x[1],x[2],rho,g);
         if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(computeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
            addRhoRCP(x,sig,lbl,mypos);
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
/* ************************************************************************************ */
bool critPtNetWork::setRhoCCPs(void)
{
   solreal x[3],rho,g[3];
   string lbl;
   int sig,mypos;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for (int i=0; i<nRCP; i++) {
      for (int j=(i+1); j<nRCP; j++) {
         for (int k=0; k<3; k++) {x[k]=(RRCP[i][k]-RRCP[j][k]);}
         if (computeMagnitudeV3(x)>(bn->maxBondDist*CPNW_BONDISTEXTCCPSEARCHFACTOR)) {continue;}
         for (int k=0; k<3; k++) {x[k]=0.5e0*(RRCP[i][k]+RRCP[j][k]);}
         seekRhoCCP(x,rho,g,sig);
         lbl=(lblRCP[i]+string("-")+lblRCP[j]);
         wf->evalRhoGradRho(x[0],x[1],x[2],rho,g);
         if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(computeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
            addRhoCCP(x,sig,lbl,mypos);
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
      //cout << lblCCP[i] << endl;
   }
   cout << nCCP << " CCPs found.\n";
   //displayWarningMessage("setRhoCCPs(...) under construction.");
   return true;
}
/* ************************************************************************************ */
bool critPtNetWork::setLOLACPs()
{
   solreal x[3],rho,g[3],magg;
   string lbl;
   int sig;
   bool tmpbool=false;
   for (int i=0; i<(bn->nNuc); i++) {
      for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
      seekLOLACP(x,rho,g,sig);
      magg=computeMagnitudeV3(g);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSLOLACPGRADMAG) ) {
         lbl=bn->atLbl[i];
         addRhoACP(x,sig,lbl);
      } else {
         for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
         lbl=bn->atLbl[i]+"+";
         addRhoACP(x,sig,lbl);
         tmpbool=true;
      }

   }
   if ( tmpbool ) {
      displayWarningMessage("Some ACPs presented |nabla rho|!= 0... Look for '+' in labels.");
   }
   if ( nACP!=(bn->nNuc) ) {
      displayWarningMessage("Number of regular ACPs is different from number of Nuclei!");
   }
   int regnACP=nACP;
   cout << "Found " << regnACP << " regular ACPs (" << bn->nNuc << " nuc)" << endl;
   cout << "Looking around nuclei..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal xs[3],rad;
   for ( int i=0 ; i<(bn->nNuc) ; i++ ) {
         for ( int k=0 ; k<3 ; k++ ) {xs[k]=(bn->R[i][k]);}
         rad=computeMagnitudeV3(xs);
         lbl="NNLOLACP"+getStringFromInt(i);
         rad*=0.01e0;
         seekLOLACPsAroundAPoint(xs,rad,lbl,-1);
#if USEPROGRESSBAR
         printProgressBar(int(100.0e0*solreal(i)/solreal((bn->nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   /* -------------------------------------------------------------------------------  */
   //*
   cout << "Looking between atom pairs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   //solreal xs[3],rad;
   for ( int i=0 ; i<(bn->nNuc) ; i++ ) {
      for ( int j=(i+1) ; j<(bn->nNuc) ; j++ ) {
         for ( int k=0 ; k<3 ; k++ ) {xs[k]=(bn->R[i][k]-bn->R[j][k]);}
         rad=computeMagnitudeV3(xs);
         if ( rad<=(2.0e0*bn->maxBondDist) ) {
            for ( int k=0 ; k<3 ; k++ ) {xs[k]=0.5e0*(bn->R[i][k]+bn->R[j][k]);}
            lbl="NNLOLACP"+getStringFromInt((i+1));
            rad*=0.1e0;
            seekLOLACPsAroundAPoint(xs,rad,lbl,-1);
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((bn->nNuc-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   // */
   if (regnACP==nACP) {
      cout << "No more ACPs found." << endl;
   } else {
      cout << (nACP-regnACP) << " new ACP";
      if ((nACP-regnACP)>1) {cout << "s";}
      cout << " found! Total number of ACPs: " << nACP << endl;
   }
   return true;
}
/* ************************************************************************************ */
bool critPtNetWork::setLOLBCPs(void)
{
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   int ata,atb,k,sig;
   int mypos;
   solreal x[3],g[3],rho;
   string lbl;
   for (int i=0; i<(bn->nNuc); i++) {
      ata=i;
      atb=bn->bNet[ata][0];
      k=1;
      while ((k<MAXBONDINGATOMS)&&(atb>0)) {
         for (int j=0; j<3; j++) {x[j]=0.5e0*(bn->R[ata][j]+bn->R[atb][j]);}
         seekLOLBCP(x,rho,g,sig);
         if (ata<atb) {
            lbl=bn->atLbl[ata]+string("-")+bn->atLbl[atb];
         } else {
            lbl=bn->atLbl[atb]+string("-")+bn->atLbl[ata];
         }
         if (addRhoBCP(x,sig,lbl,mypos)&&mypos>=0) {
            atBCP[mypos][0]=ata;
            atBCP[mypos][1]=atb;
         }
         atb=bn->bNet[ata][k];
         ++k;
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
   normalbcp=nBCP;
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   /* ------------------------------------------------------ */
   cout << nBCP << " BCPs found.\n";
   string ll;
   cout << "Including non-nuclear attractors..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for ( int i=(bn->nNuc) ; i<nACP ; i++ ) {
      seekLOLBCPWithExtraACP(i,(bn->maxBondDist));
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nACP-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   /*
   cout << "Looking for possible extra BCPs (Hydrogen bonds)..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   solreal extbd=(bn->maxBondDist*CPNW_EXTENDEDBONDDISTFACTOR),dd=0.0e0;
   for (int i=0; i<nACP; i++) {
      for (int j=(i+1); j<nACP; j++) {
         for (int k=0; k<3; k++) {x[k]=((RACP[i][k]-RACP[j][k]));}
         dd=computeMagnitudeV3(x);
         if ((dd>=(bn->maxBondDist*0.9e0))&&(dd<=extbd)) {
            for (int k=0; k<3; k++) {x[k]=0.5e0*(RACP[i][k]+RACP[j][k]);}
            if ((wf->evalDensity(x[0],x[1],x[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
               seekRhoBCP(x,rho,g,sig);
               if (i<j) {ata=i; atb=j;} else {ata=j; atb=i;}
               lbl=string("*")+lblACP[ata]+string("-")+lblACP[atb];
               wf->evalRhoGradRho(x[0],x[1],x[2],rho,g);
               if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(computeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
                  addRhoBCP(x,sig,lbl,mypos);
                  if ((mypos>=0)&&(mypos<dBCP)) {
                     atBCP[mypos][0]=ata;
                     atBCP[mypos][1]=atb;
                  }
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
         findTwoClosestACPs(tbcp,ata,atb);
         atBCP[i][0]=ata;
         atBCP[i][1]=atb;
         lblBCP[i]="*"+lblACP[ata]+"-"+lblACP[atb];
      }
      cout << (nBCP-normalbcp) << " new BCP";
      if ((nBCP-normalbcp)>1) {cout << "s";}
      cout << " found! Total number of BCPs: " << nBCP << endl;
   }
   // */
   return true;
}
/* ************************************************************************************* */
bool critPtNetWork::setLOLRCPs(void)
{
   displayWarningMessage("No LOL RCP will be look for...");
   displayWarningMessage("(critPtNetWork::setLOLRCPs(...) under construction.)");
   return false;
}
/* ************************************************************************************* */
bool critPtNetWork::setLOLCCPs(void)
{
   displayWarningMessage("No LOL CCP will be look for...");
   displayWarningMessage("(critPtNetWork::setLOLCCPs(...) under construction.)");
   return false;
}
/* ************************************************************************************ */
bool critPtNetWork::addRhoACP(solreal (&x)[3],int sig,string &lbl)
{
   if ( sig!=-3 ) {return false;}
   if (nACP==0) {
      for (int i=0; i<3; i++) {RACP[0][i]=x[i];}
      lblACP[0]=lbl;
      //(lbl[lbl.length()-1])++; //Check what is this for...
      ++nACP;
      return true;
   }
   //if ((x[0]<bn->bbmin[0])||(x[0]>bn->bbmax[0])) {cout << "Out of box (x)...\n"; return false;}
   //if ((x[1]<bn->bbmin[1])||(x[1]>bn->bbmax[1])) {cout << "Out of box (y)...\n"; return false;}
   //if ((x[2]<bn->bbmin[2])||(x[2]>bn->bbmax[2])) {cout << "Out of box (z)...\n"; return false;}
   size_t pos;
   if (imNew(x,dACP,RACP,pos)) {
      for (int i=0; i<3; i++) {RACP[pos][i]=x[i];}
      lblACP[pos]=lbl;
      //(lbl[lbl.length()-1])++;
      ++nACP;
      return true;
   }
   return false;
}
/* ************************************************************************************ */
bool critPtNetWork::addRhoBCP(solreal (&x)[3],int sig,string &lbl,int &pos)
{
   if ( sig!=-1 ) {
      pos=-1;
      return false;
   }
   if (nBCP==0) {
      for (int i=0; i<3; i++) {RBCP[0][i]=x[i];}
      lblBCP[0]=lbl;
      ++nBCP;
      pos=0;
      return true;
   }
   pos=-1;
   //if ((x[0]<bn->bbmin[0])||(x[0]>bn->bbmax[0])) {cout << "Out of box (x)...\n"; return false;}
   //if ((x[1]<bn->bbmin[1])||(x[1]>bn->bbmax[1])) {cout << "Out of box (y)...\n"; return false;}
   //if ((x[2]<bn->bbmin[2])||(x[2]>bn->bbmax[2])) {cout << "Out of box (z)...\n"; return false;}
   size_t ttpos;
   if (imNew(x,dBCP,RBCP,ttpos)) {
      for (int i=0; i<3; i++) {RBCP[ttpos][i]=x[i];}
      lblBCP[ttpos]=lbl;
      ++nBCP;
      pos=int(ttpos);
      return true;
   }
   return false;
}
/* ************************************************************************************ */
bool critPtNetWork::addRhoRCP(solreal (&x)[3],int sig,string &lbl,int &pos)
{
   if ( sig!=1 ) {
      pos=-1;
      return false;
   }
   if (nRCP==0) {
      for (int i=0; i<3; i++) {RRCP[0][i]=x[i];}
      lblRCP[0]=lbl;
      ++nRCP;
      pos=0;
      return true;
   }
   pos=-1;
   //if ((x[0]<bn->bbmin[0])||(x[0]>bn->bbmax[0])) {cout << "Out of box (x)...\n"; return false;}
   //if ((x[1]<bn->bbmin[1])||(x[1]>bn->bbmax[1])) {cout << "Out of box (y)...\n"; return false;}
   //if ((x[2]<bn->bbmin[2])||(x[2]>bn->bbmax[2])) {cout << "Out of box (z)...\n"; return false;}
   size_t ttpos;
   if (imNew(x,dRCP,RRCP,ttpos)) {
      for (int i=0; i<3; i++) {RRCP[ttpos][i]=x[i];}
      lblRCP[ttpos]=lbl;
      ++nRCP;
      pos=int(ttpos);
      return true;
   } else {
      if ((int(ttpos)<nRCP)&&(ttpos!=string::npos)) {
         lblRCP[ttpos]+=(string("-")+lbl);
      }
   }
   return false;
}
/* ************************************************************************************ */
bool critPtNetWork::addRhoCCP(solreal (&x)[3],int sig,string &lbl,int &pos)
{
   if ( sig!=3 ) {
      pos=-1;
      return false;
   }
   if (nCCP==0) {
      for (int i=0; i<3; i++) {RCCP[0][i]=x[i];}
      lblCCP[0]=lbl;
      ++nCCP;
      pos=0;
      return true;
   }
   pos=-1;
   //if ((x[0]<bn->bbmin[0])||(x[0]>bn->bbmax[0])) {cout << "Out of box (x)...\n"; return false;}
   //if ((x[1]<bn->bbmin[1])||(x[1]>bn->bbmax[1])) {cout << "Out of box (y)...\n"; return false;}
   //if ((x[2]<bn->bbmin[2])||(x[2]>bn->bbmax[2])) {cout << "Out of box (z)...\n"; return false;}
   size_t ttpos;
   if (imNew(x,dCCP,RCCP,ttpos)) {
      for (int i=0; i<3; i++) {RCCP[ttpos][i]=x[i];}
      lblCCP[ttpos]=lbl;
      ++nCCP;
      pos=int(ttpos);
      return true;
   } else {
      if ((int(ttpos)<nCCP)&&(ttpos!=string::npos)) {
         lblCCP[ttpos]+=(string("-")+lbl);
      }
   }
   return false;
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoACPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
      string const &blbl,int uunvrt)
{
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   solreal xx[3],rho,gg[3],magg;
   int sig;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      seekRhoACP(xx,rho,gg,sig);
      magg=computeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (addRhoACP(xx,sig,lbl)) {++lbl[(lbl.length()-1)];}
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoBCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
      string const &blbl,int uunvrt)
{
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   solreal xx[3],rho,gg[3],magg;
   int sig,pos;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      seekRhoBCP(xx,rho,gg,sig);
      magg=computeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (addRhoBCP(xx,sig,lbl,pos)) {++lbl[(lbl.length()-1)];}
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoRCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
      string const &blbl,int uunvrt)
{
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   solreal xx[3],rho,gg[3],magg;
   int sig,pos;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      seekRhoRCP(xx,rho,gg,sig);
      magg=computeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (addRhoRCP(xx,sig,lbl,pos)) {++lbl[(lbl.length()-1)];}
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoCCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
      string const &blbl,int uunvrt)
{
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   solreal xx[3],rho,gg[3],magg;
   int sig,pos;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      seekRhoCCP(xx,rho,gg,sig);
      magg=computeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (addRhoCCP(xx,sig,lbl,pos)) {++lbl[(lbl.length()-1)];}
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoBCPWithExtraACP(int acppos,solreal maxrad)
{
   solreal xx[3],rho,g[3];
   int sig;
   string lbl="";
   int mypos=-1;
   for ( int i=0 ; i<nACP ; i++ ) {
      if ( i==acppos ) {continue;}
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=(RACP[acppos][k]-RACP[i][k]);}
      if ( computeMagnitudeV3(xx)<=maxrad ) {
         for (int k=0; k<3; k++) {xx[k]=0.5e0*(RACP[i][k]+RACP[acppos][k]);}
         if ((wf->evalDensity(xx[0],xx[1],xx[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
            seekRhoBCP(xx,rho,g,sig);
            if (i<acppos) {
               lbl=string("*")+lblACP[i]+string("-")+lblACP[acppos];
            } else {
               lbl=string("*")+lblACP[acppos]+string("-")+lblACP[i];
            }
            wf->evalRhoGradRho(xx[0],xx[1],xx[2],rho,g);
            if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&
                  (computeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
               addRhoBCP(xx,sig,lbl,mypos);
               if ((mypos>=0)&&(mypos<dBCP)) {
                  atBCP[mypos][0]=i;
                  atBCP[mypos][1]=acppos;
               }
            }
         }
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLBCPWithExtraACP(int acppos,solreal maxrad)
{
   solreal xx[3],lol,g[3],hl[3][3];
   int sig;
   string lbl="";
   int mypos=-1;
   for ( int i=0 ; i<nACP ; i++ ) {
      if ( i==acppos ) {continue;}
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=(RACP[acppos][k]-RACP[i][k]);}
      if ( computeMagnitudeV3(xx)<=maxrad ) {
         for (int k=0; k<3; k++) {xx[k]=0.5e0*(RACP[i][k]+RACP[acppos][k]);}
         if ((wf->evalDensity(xx[0],xx[1],xx[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
            seekLOLBCP(xx,lol,g,sig);
            if (i<acppos) {
               lbl=string("*")+lblACP[i]+string("-")+lblACP[acppos];
            } else {
               lbl=string("*")+lblACP[acppos]+string("-")+lblACP[i];
            }
            wf->evalHessLOL(xx,lol,g,hl);
            if (computeMagnitudeV3(g)<CPNW_EPSLOLACPGRADMAG) {
               addRhoBCP(xx,sig,lbl,mypos);
               if ((mypos>=0)&&(mypos<dBCP)) {
                  atBCP[mypos][0]=i;
                  atBCP[mypos][1]=acppos;
               }
            }
         }
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLACPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
      string const &blbl,int uunvrt)
{
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   solreal xx[3],rho,gg[3],magg;
   int sig;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      seekLOLACP(xx,rho,gg,sig);
      magg=computeMagnitudeV3(gg);
      if ( /*(rho>CPNW_MINLOLSIGNIFICATIVEVAL)&&*/(magg<CPNW_EPSLOLACPGRADMAG) ) {
         if (addRhoACP(xx,sig,lbl)) {++lbl[(lbl.length()-1)];}
      }
   }
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLACPsOnASphere(int atIdx,int nDivR,int nDivT,int nDivP,\
      solreal radmin,solreal radmax)
{
   static const solreal ppii=4.0e0*atan(1.0e0);
   if ( atIdx<0 || atIdx>=bn->nNuc ) {displayErrorMessage("Requesting a non existent atom!");}
   if ( nDivR<=0 || nDivT<=0 || nDivP<=0 ) {
      displayErrorMessage("number of points must be greater than zero!");
   }
   cout << "Looking for LOL ACPs around atom " << (atIdx+1) << "(" 
        << wf->atLbl[atIdx] << ")" << endl;
   solreal xx[3],lol,gl[3];
   int sig,nseeds,count;
   string blbl="LOLACPSp"+getStringFromInt(atIdx+1)+"*",lbl;
   solreal cc[3],dr,dt,dp,currr,currt,currp;
   dr=(radmax-radmin)/solreal(nDivR-1);
   dt=ppii/solreal(nDivT-1);
   dp=2.0e0*ppii/solreal(nDivP);
   nseeds=nDivR*nDivT*nDivP;
   count=0;
   int origacp=nACP,idxLbl=1;
   lbl=blbl+getStringFromInt(idxLbl);
   for ( int i=0 ; i<3 ; i++ ) {cc[i]=bn->R[atIdx][i];}
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for ( int ir=0 ; ir<nDivR ; ir++ ) {
      currr=radmin+solreal(ir)*dr;
      for ( int it=0 ; it<nDivT ; it++ ) {
         currt=(dt*solreal(it));
         for ( int ip=0 ; ip<nDivP ; ip++ ) {
            currp=dp*solreal(ip);
            xx[0]=cc[0]+currr*sin(currt)*cos(currp);
            xx[1]=cc[1]+currr*sin(currt)*sin(currp);
            xx[2]=cc[2]+currr*cos(currt);
            seekLOLACP(xx,lol,gl,sig);
            if ( computeMagnitudeV3(gl)<CPNW_EPSLOLACPGRADMAG ) {
               if (addRhoACP(xx,sig,lbl)){lbl=blbl+getStringFromInt(++idxLbl);}
            }
            ++count;
         }
      }
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(count)/solreal(nseeds-1)));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   if ( origacp<nACP ) {cout << "Found " << (nACP-origacp) << " new ACPs!..." << endl;}
}
/* ************************************************************************************ */
void critPtNetWork::extendedSearchCPs(void)
{
   solreal xs[3],dx;
   dx=bn->maxBondDist*0.1e0;
   string lbl="";
   int count=0;
   int initacp=nACP,initbcp=nBCP,initrcp=nRCP,initccp=nCCP;
   cout << "Looking around extraBCPs..." << endl;
   for ( int i=normalbcp ; i<nBCP ; i++ ) {
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
      for ( int k=0 ; k<3 ; k++ ) {xs[k]=RBCP[i][k];}
      lbl="extACPb";
      seekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extRCPb";
      seekRhoRCPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extCCPb";
      seekRhoCCPsAroundAPoint(xs,dx,lbl,nIHV);
      ++count;
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(count)/solreal(nBCP-normalbcp-1)));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   count=0;
   cout << "Looking around RCPs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for ( int i=0 ; i<nRCP ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xs[k]=RRCP[i][k];}
      lbl="extACPr";
      seekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extBCPr";
      seekRhoBCPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extCCPr";
      seekRhoCCPsAroundAPoint(xs,dx,lbl,nIHV);
      ++count;
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(count)/solreal(nBCP-normalbcp-1)));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   count=0;
   cout << "Looking around CCPs..." << endl;
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   for ( int i=0 ; i<nCCP ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xs[k]=RCCP[i][k];}
      lbl="extACPc";
      seekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extBCPc";
      seekRhoBCPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extRCPc";
      seekRhoRCPsAroundAPoint(xs,dx,lbl,nIHV);
      ++count;
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(count)/solreal(nBCP-normalbcp-1)));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   for (int i=0; i<nRCP; i++) {removeRedundInLabel(lblRCP[i]);}
   for (int i=0; i<nCCP; i++) {removeRedundInLabel(lblCCP[i]);}
   bool foundnewcps=false;
   if ( initacp<nACP ) {
      cout << "Found " << (nACP-initacp) << " new ACPs." << endl;
      foundnewcps=true;
   }
   if ( initbcp<nBCP ) {
      cout << "Found " << (nBCP-initbcp) << " new BCPs." << endl;
      foundnewcps=true;
   }
   if ( initrcp<nRCP ) {
      cout << "Found " << (nRCP-initrcp) << " new RCPs." << endl;
      foundnewcps=true;
   }
   if ( initccp<nCCP ) {
      cout << "Found " << (nCCP-initccp) << " new CCPs." << endl;
      foundnewcps=true;
   }
   if ( foundnewcps ) {
      printBetweenStarLines(string("nACP-nBCP+nRCP-nCCP = "+\
               getStringFromInt(nACP-nBCP+nRCP-nCCP)));
   }
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
void critPtNetWork::removeRedundInLabel(string &lbl)
{
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
/* ************************************************************************************ */
string critPtNetWork::getFirstChunkOfLabel(string &lbl)
{
   if (lbl.length()==0) {
      return string("");
   }
   size_t pos=lbl.find_first_of("-");
   return lbl.substr(0,pos);
}
/* ************************************************************************************ */
void critPtNetWork::printAllFieldProperties(solreal x,solreal y,solreal z)
{
   wf->displayAllFieldProperties(x,y,z);
}
/* ************************************************************************************ */
void critPtNetWork::writeAllFieldProperties(solreal x,solreal y,solreal z,ofstream &ofil)
{
   wf->writeAllFieldProperties(x,y,z,ofil);
}
/* ************************************************************************************ */
void critPtNetWork::getACPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   solreal eive[3][3],b[3],F[3];
   eigen_decomposition3(hess, eive, b);
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   solreal h4[4][4],m4[4][4],v4[4];
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
   if ( fabs(lp)<CPNW_EPSEIGENVALUECPSEARCH ) {lp=CPNW_EPSEIGENVALUECPSEARCH;}
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
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp+CPNW_EPSRHOACPGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-lp+CPNW_EPSRHOACPGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-lp+CPNW_EPSRHOACPGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>stepSizeACP) {
         hh[i]=SIGNF(hh[i])*stepSizeACP;
      }
   }
   sig=computeSignature(b);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getBCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   solreal eive[3][3],b[3];
   eigen_decomposition3(hess, eive, b);
   solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   solreal h3[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {h3[i][j]=0.00000e0;}
   }
   solreal ln=0.5e0*(b[2]-sqrt(b[2]*b[2]+4.0e0*F[2]*F[2]));
   h3[0][0]=b[0];
   h3[1][1]=b[1];
   h3[2][0]=h3[0][2]=F[0];
   h3[2][1]=h3[1][2]=F[1];
   solreal m3[3][3],vv[3];
   eigen_decomposition3(h3, m3, vv);
   solreal lp=vv[2];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp+CPNW_EPSRHOACPGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-lp+CPNW_EPSRHOACPGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln+CPNW_EPSRHOACPGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>stepSizeBCP) {
         hh[i]=SIGNF(hh[i])*stepSizeBCP;
      }
   }
   sig=computeSignature(b);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getRCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   solreal eive[3][3],b[3];
   eigen_decomposition3(hess, eive, b);
   solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   solreal h3[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {h3[i][j]=0.00000e0;}
   }
   solreal lp=0.5e0*(b[0]+sqrt(b[0]*b[0]+4.0e0*F[0]*F[0]));
   h3[0][0]=b[1];
   h3[1][1]=b[2];
   h3[2][0]=h3[0][2]=F[1];
   h3[2][1]=h3[1][2]=F[2];
   solreal m3[3][3],vv[3];
   eigen_decomposition3(h3, m3, vv);
   solreal ln=vv[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp+CPNW_EPSRHOACPGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln+CPNW_EPSRHOACPGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln+CPNW_EPSRHOACPGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>stepSizeRCP) {
         hh[i]=SIGNF(hh[i])*stepSizeRCP;
      }
   }
   sig=computeSignature(b);
   return;
}
/* ************************************************************************************* */
void critPtNetWork::getCCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig)
{
   solreal eive[3][3],b[3];
   eigen_decomposition3(hess, eive, b);
   solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   solreal h4[4][4],m4[4][4],v4[4];
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
      hh[0]-=eive[0][j]*F[j]/(b[j]-ln+CPNW_EPSRHOACPGRADMAG);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln+CPNW_EPSRHOACPGRADMAG);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln+CPNW_EPSRHOACPGRADMAG);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>stepSizeCCP) {
         hh[i]=SIGNF(hh[i])*stepSizeCCP;
      }
   }
   sig=computeSignature(b);
   return;
}
/* ************************************************************************************ */
int critPtNetWork::computeSignature(solreal (&ev)[3])
{
   int res=0;
   for ( int i=0 ; i<3 ; i++ ) {
      res+=SIGNF(ev[i]);
      if ( ev[i]==0.0e0 ) {
         res=0;
#if DEBUG
         displayWarningMessage("Zero signature!");
         printV3Comp("EigenValues: ",ev);
#endif
         break;
      }
   }
   return res;
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoACP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig)
{
   solreal rho,gr[3],hr[3][3],dx[3];
   wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=computeSignature(hr);
   solreal magd=computeMagnitudeV3(gr);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-3 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItACP)) {
      getACPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoBCP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig)
{
   solreal rho,gr[3],hr[3][3],dx[3];
   wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=computeSignature(hr);
   solreal magd=computeMagnitudeV3(gr);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-1 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItBCP)) {
      getBCPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoRCP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig)
{
   solreal rho,gr[3],hr[3][3],dx[3];
   wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=computeSignature(hr);
   solreal magd=computeMagnitudeV3(gr);
   if ( magd<=CPNW_EPSRHOACPGRADMAG && sig==1 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItRCP)) {
      getRCPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekRhoCCP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig)
{
   solreal rho,gr[3],hr[3][3],dx[3];
   wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=computeSignature(hr);
   solreal magd=computeMagnitudeV3(gr);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==3 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItCCP)) {
      getCCPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->evalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLACP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig)
{
   solreal lol,gl[3],hl[3][3],dx[3];
   wf->evalHessLOL(x,lol,gl,hl);
   sig=computeSignature(hl);
   solreal magd=computeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-3 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItACP)) {
      getACPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->evalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLBCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig)
{
   solreal lol,gl[3],hl[3][3],dx[3];
   wf->evalHessLOL(x,lol,gl,hl);
   sig=computeSignature(hl);
   solreal magd=computeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-1 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItBCP)) {
      getBCPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->evalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLRCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig)
{
   solreal lol,gl[3],hl[3][3],dx[3];
   wf->evalHessLOL(x,lol,gl,hl);
   sig=computeSignature(hl);
   solreal magd=computeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==1 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItRCP)) {
      getRCPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->evalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::seekLOLCCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig)
{
   solreal lol,gl[3],hl[3][3],dx[3];
   wf->evalHessLOL(x,lol,gl,hl);
   sig=computeSignature(hl);
   solreal magd=computeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==3 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   solreal magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItRCP)) {
      getCCPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->evalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
/* ************************************************************************************ */
int critPtNetWork::computeSignature(solreal (&hh)[3][3])
{
   solreal eive[3][3],b[3];
   eigen_decomposition3(hh, eive, b);
   return computeSignature(b);
}
/* ************************************************************************************ */
bool critPtNetWork::imNew(solreal (&x)[3],int dim,solreal ** (&arr),size_t &pos)
{
   solreal ee;
   int firstzeropos=0,k=1;
   ee=arr[0][0]*arr[0][0]+arr[0][1]*arr[0][1]+arr[0][2]*arr[0][2];
   while ((ee<1.0e49)&&(k<dim)) {
      firstzeropos++;
      ee=arr[k][0]*arr[k][0]+arr[k][1]*arr[k][1]+arr[k][2]*arr[k][2];
      ++k;
   }
   if (k==dim) {
      cout << "Warning: end of the array reached, perhaps you need a larger array...\n";
      cout << "Returning false...\n";
      return false;
   }
   k=0;
   while (k<firstzeropos) {
      if ((fabs(x[0]-arr[k][0])<CPNW_EPSFABSDIFFCOORD)&&
            (fabs(x[1]-arr[k][1])<CPNW_EPSFABSDIFFCOORD)&&
            (fabs(x[2]-arr[k][2])<CPNW_EPSFABSDIFFCOORD)) {
         pos=k;
         return false;
      }
      ++k;
   }
   pos=firstzeropos;
   return true;
}
/* ************************************************************************************ */
void critPtNetWork::displayXCPCoords(char cpt)
{
   bool cpknow=false;
   string cplbl="";
   string longcplbl="";
   string *lblArray=NULL;
   int theNCP=0;
   solreal **RArray=NULL;
   switch ( cpt ) {
      case 'a' :
      case 'A' :
         cpknow=iknowacps;
         cplbl="ACP";
         longcplbl="Attractor";
         lblArray=lblACP;
         theNCP=nACP;
         RArray=RACP;
         break;
      case 'b' :
      case 'B' :
         cpknow=iknowbcps;
         cplbl="BCP";
         longcplbl="Bond";
         lblArray=lblBCP;
         theNCP=nBCP;
         RArray=RBCP;
         break;
      case 'r' :
      case 'R' :
         cpknow=iknowrcps;
         cplbl="RCP";
         longcplbl="Ring";
         lblArray=lblRCP;
         theNCP=nRCP;
         RArray=RRCP;
         break;
      case 'c' :
      case 'C' :
         cpknow=iknowccps;
         cplbl="CCP";
         longcplbl="Cage";
         lblArray=lblCCP;
         theNCP=nCCP;
         RArray=RCCP;
         break;
      default :
         break;
   }
   if (!cpknow) {
      displayErrorMessage(string("First seek the "+cplbl+"s.\nNo "\
               +cplbl+" to display its coordinates."));
      return;
   }
   cout << scientific << setprecision(12);
   printScrCharLine('+');
   centerString(string("Coordinates of "+longcplbl+" Critical Points"));
   printScrCharLine('-');
   for (int i=0; i<theNCP; i++) {
      cout << cplbl << "[" << (i+1) << "](" << lblArray[i] << "): ";
      for (int j=0; j<3; j++) {cout << RArray[i][j] << " ";}
      cout << endl;
   }
   printScrCharLine('+');
   cout.unsetf(ios::scientific);
   return;
}
/* ************************************************************************************ */
void critPtNetWork::displayAllCPCoords(void)
{
   displayXCPCoords('a');
   displayXCPCoords('b');
   displayXCPCoords('r');
   displayXCPCoords('c');
}
/* ************************************************************************************ */
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
void critPtNetWork::findTwoClosestAtoms(solreal (&xo)[3],int &idx1st,int &idx2nd)
{
   if ( wf->nNuc<2 ) {idx1st=0; idx2nd=0; return;}
   solreal xmagt=0.0e0,xmag1=0.0e0,xmag2=0.0e0;
   int ii1,ii2,iit;
   for ( int k=0 ; k<3 ; k++ ) {xmag1+=((xo[k]-wf->getR(0,k))*(xo[k]-wf->getR(0,k)));}
   ii1=0;
   for ( int k=0 ; k<3 ; k++ ) {xmag2+=((xo[k]-wf->getR(1,k))*(xo[k]-wf->getR(1,k)));}
   ii2=1;
   if ( xmag1>xmag2 ) {
      xmagt=xmag1; xmag1=xmag2; xmag2=xmagt;
      ii1=1;
      ii2=0;
   }
   for ( int i=2 ; i<wf->nNuc ; i++ ) {
      xmagt=0.0e0;
      for ( int k=0 ; k<3 ; k++ ) {xmagt+=((xo[k]-wf->getR(i,k))*(xo[k]-wf->getR(i,k)));}
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
   return;
}
/* ************************************************************************************ */
void critPtNetWork::findTwoClosestACPs(solreal (&xo)[3],int &idx1st,int &idx2nd)
{
   if ( nACP<2 ) {idx1st=0; idx2nd=0; return;}
   solreal xmagt=0.0e0,xmag1=0.0e0,xmag2=0.0e0;
   int ii1,ii2,iit;
   for ( int k=0 ; k<3 ; k++ ) {xmag1+=((xo[k]-RACP[0][k])*(xo[k]-RACP[0][k]));}
   ii1=0;
   for ( int k=0 ; k<3 ; k++ ) {xmag2+=((xo[k]-RACP[1][k])*(xo[k]-RACP[1][k]));}
   ii2=1;
   if ( xmag1>xmag2 ) {
      xmagt=xmag1; xmag1=xmag2; xmag2=xmagt;
      ii1=1;
      ii2=0;
   }
   for ( int i=2 ; i<nACP ; i++ ) {
      xmagt=0.0e0;
      for ( int k=0 ; k<3 ; k++ ) {xmagt+=((xo[k]-RACP[i][k])*(xo[k]-RACP[i][k]));}
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
      displayWarningMessage("Identical ACPs!");
      DISPLAYDEBUGINFOFILELINE;
   }
#endif /* ( DEBUG ) */
   return;
}
/* ************************************************************************************ */
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
/* ************************************************************************************ */
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
/* ************************************************************************************ */
void critPtNetWork::writeCPProps(string &ofnam,string &wfnnam)
{
   ofstream ofil;
   ofil.open(ofnam.c_str());
   solreal x,y,z;
   //solreal gx,gy,gz,rho,hxx,hyy,hzz,hxy,hxz,hyz;
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
   ofil << "nACP: " << nACP << "\nnBCP: " << nBCP
        << "\nnRCP: " << nRCP << "\nnCCP: " << nCCP << endl;
   ofil << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
   writeScrStarLine(ofil);
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
         writeAllFieldProperties(x,y,z,ofil);
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
         writeAllFieldProperties(x,y,z,ofil);
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
         writeAllFieldProperties(x,y,z,ofil);
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
         writeAllFieldProperties(x,y,z,ofil);
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
   ofil.close();
   return;
}
/* ************************************************************************************ */
bool critPtNetWork::makePOVFile(string pnam,povRayConfProp &pvp,int campos)
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
   solreal allcprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
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
   if (!(bn->imstp())) {bn->setUpBNW();}
   if (!(bn->ballAndStickMode)) {bn->drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
#if DEBUG
   displayWarningMessage("In this version, calling critPtNetWork::makePovFILE(...)\n\
         will overwrite the original coordinates of the critical points\n\
         and the coordinates on the bondnetwork object as well.");
#endif
   centerMolecule();
   bn->calcViewRadius();
   //cout << "rView: " << bn->rView << endl;
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
      pvp.locCam[i]*=bn->rView;
      for (int j=0; j<2; j++) {
         pvp.lightSource[j][i]*=(bn->rView*2.0e0);
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
      for (int j=0; j<3; j++) {pvp.lightSource[i][j]*=(bn->rView*4.0e0);}
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
   if (bn->spaceFillingMode) {
      pof << indTabsStr(pvp.currIndLev) << "merge{" << endl;
      pvp.currIndLev++;
   }
   writeScrCharLine(pof,'/');
   pof << "#if(DrawAtomTranspSpheres)" << endl;
   putNuclei(pof);
   pof << "#end\n//end if DrawAtomTranspSpheres" << endl;
   writeScrCharLine(pof,'/');
   pof << "#if(DrawStandardBonds)" << endl;
   putBonds(pof);
   pof << "#end\n//end if DrawStandardBonds" << endl;
   writeScrCharLine(pof,'/');
   if (iknowacps) {
      //int indacp;
      pof << "#if(DrawAttractorCriticalPoints)" << endl;
      //solreal acprad;
      for (int i=0; i<nACP; i++) {
         //acprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
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
         //bcprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
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
         rcprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
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
         ccprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
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
            writePOVCylinder(pof,1,\
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
void critPtNetWork::putBonds(ofstream &pof)
{
   pof << "union{" << endl;
   int k=0,atni,atnk;
   solreal startpt[3],frak1;
   for (int i=0; i<bn->nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bn->bNet[i][j];
         atni=bn->atNum[i];
         atnk=bn->atNum[k];
         //frak1=atomicRadius[atni]/(atomicRadius[atni]+atomicRadius[atnk]);
         frak1=getAtomicVDWRadius(atni)/(getAtomicVDWRadius(atni)+getAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=bn->R[i][l]*(1.0e0-frak1)+bn->R[k][l]*frak1;
         }
         if (k>0) {
            writePOVCylinder(pof,1,
                  bn->R[i][0],bn->R[i][1],bn->R[i][2],
                  startpt[0],startpt[1],startpt[2],
                  bn->drawStickSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR,
                  getAtomicRColorReal(atni),getAtomicGColorReal(atni),
                  getAtomicBColorReal(atni));
            writePOVCylinder(pof,1,
                  startpt[0],startpt[1],startpt[2],
                  bn->R[k][0],bn->R[k][1],bn->R[k][2],
                  bn->drawStickSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR,
                  getAtomicRColorReal(atnk),getAtomicGColorReal(atnk),
                  getAtomicBColorReal(atnk));
         }
      }
      writePOVSphere(pof,0,bn->R[i][0],bn->R[i][1],bn->R[i][2],
            bn->drawStickSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR,
            getAtomicRColorReal(atni),getAtomicGColorReal(atni),
            getAtomicBColorReal(atni));
   }
   pof << "}" << endl;
   return;
}
/* ************************************************************************************ */
void critPtNetWork::putNuclei(ofstream & pof)
{
   int atomn;
   solreal atrad;
   for (int i=0; i<bn->nNuc; i++) {
      atomn=bn->atNum[i];
      atrad=bn->drawAtSize;
      writePOVTransparentSphere(pof,0,bn->R[i][0],bn->R[i][1],bn->R[i][2],atrad,
            getAtomicRColorReal(atomn),getAtomicGColorReal(atomn),
            getAtomicBColorReal(atomn),0.7);
   }
   return;
}
/* ************************************************************************************ */
void critPtNetWork::centerMolecule(void)
{
   solreal trn[3];
   for (int i=0; i<3; i++) {
      trn[i]=0.5e0*(bn->rmax[i]+bn->rmin[i]);
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
   bn->centerMolecule();
   for (int i=0; i<3; i++) {centMolecVec[i]=trn[i];}
   return;
}
/* ************************************************************************************ */
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
         dACP=nACP*CPNW_MAXRHOACPSPERATOM;
         break;
      case LOLD:
         dACP=nACP*CPNW_MAXLOLACPSPERATOM;
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
   if (dBCP>0) {alloc3DRealArray(string("RBGP"),dBCP,CPNW_ARRAYSIZEGRADPATH,3,RBGP);}
   if (nBGP>=0&&atBCP!=NULL) {
      cpxGetNOfPtsPerBondPath(cfil,nBGP,atBCP);
      cpxGetBondPathData(cfil,nBGP,atBCP,RBGP);
      iknowbgps=true;
   } else {iknowbgps=false;}
   cout << "Critical Point State Loaded!" << endl;
   displayStatus(true);
   return true;
}
/* ************************************************************************************ */
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
/* ************************************************************************************ */
int critPtNetWork::findSingleRhoGradientPathRK5(int at1,int at2,solreal hstep,\
      int dima,solreal** (&arbgp),solreal (&ro)[3])
{
   solreal rn[3],rho,g[3],h[3][3],eive[3][3],eival[3],dist,maggrad;
   seekSingleRhoBCP(at1,at2,ro);
   //cout << "ro: " << ro[0] << " " << ro[1] << " " << ro[2] << endl;
   wf->evalHessian(ro[0],ro[1],ro[2],rho,g,h);
   eigen_decomposition3(h,eive,eival);
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {maggrad+=(eive[i][2]*eive[i][2]);}
   maggrad=sqrt(maggrad);
   //hstep=DEFAULTHGRADIENTPATHS;
   for (int i=0; i<3; i++) {rn[i]=ro[i]+hstep*eive[i][2]/maggrad;}
   wf->evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
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
      getNextPointInGradientPathRK5(rn,hstep,maggrad);
      for (int i=0; i<3; i++) {arbgp[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf->R[iacp1+i])*(rn[i]-wf->R[iacp1+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {
         iminacp=true;
         count++;
         for (int i=0; i<3; i++) {arbgp[count][i]=wf->R[iacp1+i];}
      }
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf->R[iacp2+i])*(rn[i]-wf->R[iacp2+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {
         iminacp=true;
         count++;
         for (int i=0; i<3; i++) {arbgp[count][i]=wf->R[iacp2+i];}
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
   wf->evalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
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
      getNextPointInGradientPathRK5(rn,hstep,maggrad);
      for (int i=0; i<3; i++) {arbgp[count][i]=rn[i];}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf->R[iacp1+i])*(rn[i]-wf->R[iacp1+i]));}
      dist=sqrt(dist);
      if (dist<=hstep) {iminacp=true;}
      dist=0.0e0;
      for (int i=0; i<3; i++) {dist+=((rn[i]-wf->R[iacp2+i])*(rn[i]-wf->R[iacp2+i]));}
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
   for (int i=0; i<3; i++) {dist+=((arbgp[count-1][i]-wf->R[iacp1+i])*(arbgp[count-1][i]-wf->R[iacp1+i]));}
   dist=sqrt(dist);
   if (dist<=hstep) {
      for (int i=0; i<3; i++) {arbgp[count][i]=wf->R[iacp1+i];}
      count++;
   }
   dist=0.0e0;
   for (int i=0; i<3; i++) {dist+=((arbgp[count-1][i]-wf->R[iacp2+i])*(arbgp[count-1][i]-wf->R[iacp2+i]));}
   dist=sqrt(dist);
   if (dist<=hstep) {
      for (int i=0; i<3; i++) {arbgp[count][i]=wf->R[iacp2+i];}
      count++;
   }
   //cout << "Checkpoint2... count: " << count  << endl;
   return count;
   //wf->displayAllFieldProperties(rn[0],rn[1],rn[2]);
}
/* ************************************************************************************ */
bool critPtNetWork::seekSingleRhoBCP(int ata,int atb,solreal (&x)[3])
{
   int idxa,idxb;
   idxa=3*ata;
   idxb=3*atb;
   //cout << "ata: " << ata << ", atb: " << atb << endl;
   for (int j=0; j<3; j++) {x[j]=0.5e0*(wf->R[idxa+j]+wf->R[idxb+j]);}
   solreal rho,g[3];
   int sig;
   seekRhoBCP(x,rho,g,sig);
   //wf->evalRhoGradRho(x[0],x[1],x[2],rho,g);
   //solreal magg=0.0e0;
   //for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
   //magg=sqrt(magg);
   if (!((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(computeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG))) {
      displayErrorMessage(string("The chosen atoms ("+wf->atLbl[ata]\
               +","+wf->atLbl[atb]+") do not lead to a BCP!"));
      //wf->displayAllFieldProperties(x[0],x[1],x[2]);
      return false;
   }
   return true;
}
/* ************************************************************************************ */
void critPtNetWork::getNextPointInGradientPathRK5(solreal (&xn)[3],solreal &stepsize,\
      solreal &mgg)
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
      wf->evalRhoGradRho(xt[0],xt[1],xt[2],rho,g);
      maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      for(int l=0; l<3; l++) {
         k[i][l] = stepsize*g[l]/maggrad;
      }
   }
   for(int i=0; i<3; i++) {
      xn[i]+=(c1*k[0][i]+c3*k[2][i]+c4*k[3][i]+c6*k[5][i]);
   }
   wf->evalRhoGradRho(xn[0],xn[1],xn[2],rho,g);
   mgg=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   return;
}
/* ************************************************************************************ */




#endif  /* _CRITPTNETWORK_CPP_ */

