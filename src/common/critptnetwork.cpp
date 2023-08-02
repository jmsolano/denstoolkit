/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/
#ifndef _CRITPTNETWORK_CPP_
#define _CRITPTNETWORK_CPP_
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setprecision;
using std::scientific;
#include "screenutils.h"
#include "fileutils.h"
#include "critptnetwork.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "eigendecompositionjama.h"
#include "iofuncts-cpx.h"
#include "atomradiicust.h"
#include "stringtools.h"
#include "povraytools.h"
#include "mymath.h"
// The first 94 atomic radii are given,
//  the rest are set to be 0.80e0
//
//#include "atomcolschcust.h" //Choose this for the palette defined by JMHP
#include "atomcolschjmol.h" //Choose this for the palette used in JMol

#define CPNW_EPSFABSDIFFCOORD (1.0e-06)
#define CPNW_MINRHOSIGNIFICATIVEVAL (5.0e-05)
#define CPNW_MINLOLSIGNIFICATIVEVAL (5.0e-04)
#define CPNW_EXTENDEDBONDDISTFACTOR (3.5e0)
#define CPNW_BONDISTEXTCCPSEARCHFACTOR (2.0e0)

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

#ifndef CPNW_MAXITERATIONRINGPATHBISECT
#define CPNW_MAXITERATIONRINGPATHBISECT 60
#endif

#ifndef CPNW_MINDISTFORPINTINPATH
#define CPNW_MINDISTFORPINTINPATH (0.6e0)
#endif

//#ifndef CPNW_EPSRHOACPGRADMAG
//#define CPNW_EPSRHOACPGRADMAG (1.0e-12)
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

void CritPtNetWork::Init() {
   //publics
   nACP=nBCP=nRCP=nCCP=0;
   normalbcp=0;
   nBGP=0;
   conBCP=NULL;
   conRCP=conCCP=NULL;
   RACP=RBCP=RRCP=RCCP=NULL;
   RGP=NULL;
   RBGP=NULL;
   RRGP=RCGP=NULL;
   lblACP=lblBCP=lblRCP=lblCCP=NULL;
   for (int i=0; i<3; i++) {centMolecVec[i]=0.0e0;}
   //privates
   dACP=dBCP=dRCP=dCCP=0;
   maxItACP=CPNW_MAXITERATIONACPSEARCH;
   maxItBCP=CPNW_MAXITERATIONBCPSEARCH;
   maxItRCP=CPNW_MAXITERATIONRCPSEARCH;
   maxItCCP=CPNW_MAXITERATIONCCPSEARCH;
   maxGradPathNPts=CPNW_ARRAYSIZEGRADPATH;
   stepSizeACP=CPNW_MAXSTEPSIZEACPRHOSEARCH;
   stepSizeBCP=CPNW_MAXSTEPSIZEBCPSEARCH;
   stepSizeRCP=CPNW_MAXSTEPSIZERCPSEARCH;
   stepSizeCCP=CPNW_MAXSTEPSIZECCPSEARCH;
   stepSizeBGP=CPNW_DEFAULTGRADIENTPATHS;
   iknowacps=iknowbcps=iknowrcps=iknowccps=false;
   iknowallcps=false;
   iknowbgps=iknowrgps=iknowcgps=false;
   iknowallgps=false;
   mycptype=NONE;
   maxBondDist=maxBCPACPDist=-1.0e+50;
   drawNuc=false;
   drawBnd=true;
   drawBGPs=drawRGPs=drawCGPs=false;
   tubeBGPStyle=false;
   wf=NULL;
   bn=NULL;
}
CritPtNetWork::CritPtNetWork(GaussWaveFunction &uwf,BondNetWork &ubn) {
   Init();
   wf=&uwf;
   bn=&ubn;
}
CritPtNetWork::~CritPtNetWork() {
   MyMemory::Dealloc4DRealArray(RCGP,dCCP,CPNW_MAXRCPSCONNECTEDTOCCP,\
         maxGradPathNPts);
   MyMemory::Dealloc4DRealArray(RRGP,dRCP,CPNW_MAXBCPSCONNECTEDTORCP,\
         maxGradPathNPts);
   MyMemory::Dealloc2DRealArray(RGP,maxGradPathNPts);
   MyMemory::Dealloc3DRealArray(RBGP,dBCP,maxGradPathNPts);
   MyMemory::Dealloc2DRealArray(RACP,dACP);
   MyMemory::Dealloc1DStringArray(lblACP);
   MyMemory::Dealloc2DRealArray(RBCP,dBCP);
   MyMemory::Dealloc1DStringArray(lblBCP);
   MyMemory::Dealloc2DIntArray(conBCP,dBCP);
   MyMemory::Dealloc2DRealArray(RRCP,dRCP);
   MyMemory::Dealloc1DStringArray(lblRCP);
   MyMemory::Dealloc3DIntArray(conRCP,dRCP,2);
   MyMemory::Dealloc3DIntArray(conCCP,dCCP,2);
   MyMemory::Dealloc2DRealArray(RCCP,dCCP);
   MyMemory::Dealloc1DStringArray(lblCCP);
   wf=NULL;
   bn=NULL;
}
double CritPtNetWork::V0=0.0e0;
double CritPtNetWork::V5=2.0e0/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
double CritPtNetWork::V8=(1.0e0+sqrt(5.0e0))/(sqrt(4.0e0+(1.0e0+sqrt(5.0e0))*(1.0e0+sqrt(5.0e0))));
double CritPtNetWork::IHV[nIHV][3]={
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
void CritPtNetWork::SetCriticalPoints(ScalarFieldType ft) {
   if (!bn->ImStp()) {
      ScreenUtils::DisplayErrorMessage("Trying to use a non set up bond network object!");
      return;
   }
   SetupACPs(ft);
   SetACPs(ft);
   SetupBCPs(ft);
   SetBCPs(ft);
   if (iknowbcps) {
      dRCP=(nBCP*(nBCP-1))/2;
      if ( dRCP<CPNW_MINARRAYSIZE ) {dRCP=CPNW_MINARRAYSIZE;}
      MyMemory::Alloc2DRealArray(string("RRCP"),dRCP,3,RRCP,1.0e+50);
      MyMemory::Alloc1DStringArray("lblRCP",dRCP,lblRCP);
      MyMemory::Alloc3DIntArray("conRCP",dRCP,2,CPNW_MAXBCPSCONNECTEDTORCP,conRCP,-1);
   }
   cout << "Looking for Ring Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowrcps=SetRhoRCPs();
         break;
      case LOLD:
         iknowrcps=SetLOLRCPs();
         break;
      default:
         ScreenUtils::DisplayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
#if DEBUG
   cout << "nRCP: " << nRCP << ", dRCP: " << dRCP << endl;
#endif
   if (iknowrcps) {
      dCCP=(nRCP*(nRCP-1))/2;
      if ( dCCP<CPNW_MINARRAYSIZE ) {dCCP=CPNW_MINARRAYSIZE;}
      MyMemory::Alloc2DRealArray(string("RCCP"),dCCP,3,RCCP,1.0e+50);
      MyMemory::Alloc1DStringArray("lblCCP",dCCP,lblCCP);
      MyMemory::Alloc3DIntArray("conCCP",dCCP,2,CPNW_MAXRCPSCONNECTEDTOCCP,conCCP,-1);
   }
   cout << "Looking for Cage Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowccps=SetRhoCCPs();
         break;
      case LOLD:
         iknowccps=SetLOLCCPs();
         break;
      default:
         ScreenUtils::DisplayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
#if DEBUG
   cout << "nCCP: " << nCCP << ", dCCP: " << dCCP << endl;
#endif
   if ( mkextsearch ) {ExtendedSearchCPs();}
   iknowallcps=(iknowacps&&iknowbcps&&iknowrcps&&iknowccps);
   if (iknowallcps) {
      ScreenUtils::PrintScrCharLine('*');
      cout << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
      ScreenUtils::PrintScrCharLine('*');
   }
}
void CritPtNetWork::SetupACPs(ScalarFieldType ft) {
   mycptype=ft;
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
         ScreenUtils::DisplayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if ( RGP!=NULL ) {
      ScreenUtils::DisplayWarningMessage("RGP already allocated! nothing to do.");
   }
   if ( RACP!=NULL ) {
      ScreenUtils::DisplayWarningMessage("RACP already allocated! Nothing to do.");
      return;
   }
   MyMemory::Alloc2DRealArray(string("RGP"),maxGradPathNPts,3,RGP,1.0e+50);
   MyMemory::Alloc2DRealArray(string("RACP"),dACP,3,RACP,1.0e+50);
   MyMemory::Alloc1DStringArray("lblACP",dACP,lblACP);
}
void CritPtNetWork::SetACPs(ScalarFieldType ft) {
   cout << "Looking for Attractor Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowacps=SetRhoACPs();
         break;
      case LOLD:
         iknowacps=SetLOLACPs();
         break;
      default:
         ScreenUtils::DisplayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
}
void CritPtNetWork::SetupBCPs(ScalarFieldType ft) {
   if ( mycptype!=ft ) {
      ScreenUtils::DisplayWarningMessage("Change of field type is not allowed, using previous type!");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
   }
   if (iknowacps) {
      if ( RBCP!=NULL ) {
         ScreenUtils::DisplayWarningMessage("RBCP already allocated! Nothing to do.");
         return;
      }
      dBCP=(nACP*(nACP-1))/2;
      if ( dBCP<CPNW_MINARRAYSIZE ) {dBCP=CPNW_MINARRAYSIZE;}
      MyMemory::Alloc2DRealArray(string("RBCP"),dBCP,3,RBCP,1.0e+50);
      MyMemory::Alloc2DIntArray(string("conBCP"),dBCP,3,conBCP,-1);
      MyMemory::Alloc1DStringArray("lblBCP",dBCP,lblBCP);
   } else {
      ScreenUtils::DisplayErrorMessage("First look for ACPs...\n");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif
   }
   return;
}
void CritPtNetWork::SetBCPs(ScalarFieldType ft) {
   cout << "Looking for Bond Critical Points..." << endl;
   switch (ft) {
      case DENS:
         iknowbcps=SetRhoBCPs();
         break;
      case LOLD:
         iknowbcps=SetLOLBCPs();
         break;
      default:
         ScreenUtils::DisplayWarningMessage("This field has not been implemented!");
         break;
   }
#if DEBUG
   cout << "nBCP: " << nBCP << ", dBCP: " << dBCP << endl;
#endif
}
void CritPtNetWork::SetupBondPaths(void) {
   MyMemory::Alloc3DRealArray(string("RBGP"),dBCP,maxGradPathNPts,3,RBGP);
}
void CritPtNetWork::SetBondPaths() {
   if (!iknowbcps) {
      ScreenUtils::DisplayErrorMessage("Please look first for the BCPs...\nNothing to be done!");
      return;
   }
   SetupBondPaths();
   int npts;
   cout << "Calculating Bond Gradient Paths..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double hstep,rseed[3];
   hstep=stepSizeBGP;
   int at1,at2;
   for (int i=0; i<nBCP; i++) {
      for (int k=0; k<3; k++) {rseed[k]=RBCP[i][k];}
      at1=conBCP[i][0];
      at2=conBCP[i][1];
      npts=FindSingleRhoBondGradientPathRK5(at1,at2,hstep,maxGradPathNPts,RBGP[i],rseed);
      conBCP[i][2]=npts;
      if (npts>0) {nBGP++;}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nBCP>1) ? (nBCP-1) : 1 )));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   if (nBGP!=nBCP) {ScreenUtils::DisplayWarningMessage("For some unknown reason nBGP!=nBCP...");}
   iknowbgps=true;
   iknowallgps=(iknowbgps&&iknowrgps&&iknowcgps);
   FindMaxBondDist();
}
bool CritPtNetWork::SetRhoACPs() {
   double x[3],rho,g[3],magg;
   string lbl;
   int sig;
   bool tmpbool=false;
   for (int i=0; i<(bn->nNuc); i++) {
      for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
      SeekRhoACP(x,rho,g,sig);
      magg=ComputeMagnitudeV3(g);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         lbl=bn->atLbl[i];
         AddRhoACP(x,sig,lbl);
      } else {
         for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
         lbl=bn->atLbl[i]+"+";
         AddRhoACP(x,sig,lbl);
         tmpbool=true;
      }
   }
   if ( tmpbool ) {
      ScreenUtils::DisplayWarningMessage("Some ACPs presented |nabla rho|!= 0... Look for '+' in labels.");
   }
   if ( nACP!=(bn->nNuc) ) {
      ScreenUtils::DisplayWarningMessage("Number of regular ACPs is different from number of Nuclei!");
   }
   int regnACP=nACP;
   cout << "Found " << regnACP << " regular ACPs (" << bn->nNuc << " nuc)" << endl;
   cout << "Looking between atom pairs..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double xs[3],rad;
   for ( int i=0 ; i<(bn->nNuc) ; i++ ) {
      for ( int j=(i+1) ; j<(bn->nNuc) ; j++ ) {
         for ( int k=0 ; k<3 ; k++ ) {xs[k]=(bn->R[i][k]-bn->R[j][k]);}
         rad=ComputeMagnitudeV3(xs);
         if ( rad<=(bn->maxBondDist) ) {
            for ( int k=0 ; k<3 ; k++ ) {xs[k]=0.5e0*(bn->R[i][k]+bn->R[j][k]);}
            lbl="NNACP";
            rad*=0.3e0;
            SeekRhoACPsAroundAPoint(xs,rad,lbl,nIHV);
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((bn->nNuc > 1)? (bn->nNuc-1) : 1)));
#endif
   }
   if (regnACP==nACP) {
      cout << "\nNo more ACPs found so far..." << endl;
   } else {
      cout << (nACP-regnACP) << " new ACP";
      if ((nACP-regnACP)>1) {cout << "s";}
      cout << " found! Total number of ACPs: " << nACP << endl;
   }
   regnACP=nACP;
   cout << "Looking between atom triads..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   int ppp=bn->nNuc;
   for ( int i=0 ; i<(ppp) ; i++ ) {
      for ( int j=(i+1) ; j<(ppp) ; j++ ) {
         for ( int l=0 ; l<3 ; ++l ) { xs[l]=(bn->R[i][l])-(bn->R[j][l]); }
         rad=ComputeMagnitudeV3(xs);
         if ( fabs(rad-(bn->maxBondDist))>1.0e-4 ) { continue; }
         for ( int k=(j+1) ; k<(ppp) ; ++k ) {
            for ( int l=0 ; l<3 ; ++l ) { xs[l]=(bn->R[j][l])-(bn->R[k][l]); }
            rad=ComputeMagnitudeV3(xs);
            if ( fabs(rad-(bn->maxBondDist))>1.0e-4 ) { continue; }
            for ( int l=0 ; l<3 ; ++l ) { xs[l]=(bn->R[i][l]); }
            for ( int l=0 ; l<3 ; ++l ) { xs[l]+=(bn->R[j][l]); }
            for ( int l=0 ; l<3 ; ++l ) { xs[l]+=(bn->R[k][l]); xs[l]/=3.0e0; }
            lbl="NNACP";
            rad=(bn->maxBondDist)*0.01e0;
            stepSizeACP=0.01;
            SeekRhoACPsAroundAPoint(xs,rad,lbl,nIHV);
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((bn->nNuc > 1)? (bn->nNuc-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
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
bool CritPtNetWork::SetRhoBCPs(void) {
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   int ata,atb,k,sig;
   int mypos;
   double x[3],g[3],rho;
   string lbl;
   for (int i=0; i<(bn->nNuc); i++) {
      ata=i;
      atb=bn->bNet[ata][0];
      k=1;
      while ((k<MAXBONDINGATOMS)&&(atb>0)) {
         for (int j=0; j<3; j++) {x[j]=0.5e0*(bn->R[ata][j]+bn->R[atb][j]);}
         SeekRhoBCP(x,rho,g,sig);
         if (ata<atb) {
            lbl=bn->atLbl[ata]+string("-")+bn->atLbl[atb];
         } else {
            lbl=bn->atLbl[atb]+string("-")+bn->atLbl[ata];
         }
         if (rho>CPNW_MINRHOSIGNIFICATIVEVAL&&AddRhoBCP(x,sig,lbl,mypos)&&mypos>=0) {
            conBCP[mypos][0]=ata;
            conBCP[mypos][1]=atb;
         }
         atb=bn->bNet[ata][k];
         ++k;
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nACP>1)? (nACP-1) : 1)));
#endif
   }
   normalbcp=nBCP;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   /* ------------------------------------------------------ */
   cout << nBCP << " BCPs found.\n";
   string ll;
   cout << "Including non-nuclear attractors..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   for ( int i=(bn->nNuc) ; i<nACP ; i++ ) {
      SeekRhoBCPWithExtraACP(i,(bn->maxBondDist));
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nACP>1)? (nACP-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   cout << "Looking for possible extra BCPs (Hydrogen bonds)..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double extbd=(bn->maxBondDist*CPNW_EXTENDEDBONDDISTFACTOR),dd=0.0e0;
   for (int i=0; i<nACP; i++) {
      for (int j=(i+1); j<nACP; j++) {
         for (int k=0; k<3; k++) {x[k]=((RACP[i][k]-RACP[j][k]));}
         dd=ComputeMagnitudeV3(x);
         if ((dd>=(bn->maxBondDist*0.9e0))&&(dd<=extbd)) {
            for (int k=0; k<3; k++) {x[k]=0.5e0*(RACP[i][k]+RACP[j][k]);}
            if ((wf->EvalDensity(x[0],x[1],x[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
               SeekRhoBCP(x,rho,g,sig);
               if (i<j) {ata=i; atb=j;} else {ata=j; atb=i;}
               lbl=string("*")+lblACP[ata]+string("-")+lblACP[atb];
               wf->EvalRhoGradRho(x[0],x[1],x[2],rho,g);
               if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(ComputeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
                  AddRhoBCP(x,sig,lbl,mypos);
                  if ((mypos>=0)&&(mypos<dBCP)) {
                     conBCP[mypos][0]=ata;
                     conBCP[mypos][1]=atb;
                  }
               }
            }
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nACP>1)? (nACP-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   if (normalbcp==nBCP) {
      cout << "No more BCPs found." << endl;
   } else {
      double tbcp[3];
      for ( int i=normalbcp ; i<nBCP ; i++ ) {
         for ( int k=0 ; k<3 ; k++ ) {tbcp[k]=RBCP[i][k];}
         FindTwoClosestACPs(tbcp,ata,atb);
         conBCP[i][0]=ata;
         conBCP[i][1]=atb;
         lblBCP[i]="*"+lblACP[ata]+"-"+lblACP[atb];
      }
      cout << (nBCP-normalbcp) << " new BCP";
      if ((nBCP-normalbcp)>1) {cout << "s";}
      cout << " found! Total number of BCPs: " << nBCP << endl;
   }
   return true;
}
bool CritPtNetWork::SetRhoRCPs(void) {
   double x[3],rho,g[3];
   string lbl;
   int sig,mypos;
   for (int i=0; i<nBCP; i++) {
      for (int j=(i+1); j<nBCP; j++) {
         for (int k=0; k<3; k++) {x[k]=(RBCP[i][k]-RBCP[j][k]);}
         if (ComputeMagnitudeV3(x)>(bn->maxBondDist*2.0e0)) {continue;}
         for ( int k=0 ; k<3 ; k++ ) {x[k]=0.5e0*(RBCP[i][k]+RBCP[j][k]);}
         SeekRhoRCP(x,rho,g,sig);
         lbl=(lblBCP[i]+string("-")+lblBCP[j]);
         //wf->EvalRhoGradRho(x[0],x[1],x[2],rho,g);
         if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(ComputeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
            AddRhoRCP(x,sig,lbl,mypos);
            if (mypos>=0) {
               AddBCP2ConRCP(mypos,i);
               AddBCP2ConRCP(mypos,j);
            }  
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nBCP>1)? (nBCP-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   //cout << "Checkpoint..." << endl;
   for (int i=0; i<nRCP; i++) {
      RemoveRedundInLabel(lblRCP[i]); //cout << "i: " << i << endl;
   }
   cout << nRCP << " RCPs found.\n";
   //ScreenUtils::DisplayWarningMessage("SetRhoRCPs(...) under construction.");
   return true;
}
bool CritPtNetWork::SetRhoCCPs(void) {
   double x[3],rho,g[3];
   string lbl;
   int sig,mypos;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   for (int i=0; i<nRCP; i++) {
      for (int j=(i+1); j<nRCP; j++) {
         for (int k=0; k<3; k++) {x[k]=(RRCP[i][k]-RRCP[j][k]);}
         if (ComputeMagnitudeV3(x)>(bn->maxBondDist*CPNW_BONDISTEXTCCPSEARCHFACTOR)) {continue;}
         for (int k=0; k<3; k++) {x[k]=0.5e0*(RRCP[i][k]+RRCP[j][k]);}
         SeekRhoCCP(x,rho,g,sig);
         lbl=(lblRCP[i]+string("-")+lblRCP[j]);
         wf->EvalRhoGradRho(x[0],x[1],x[2],rho,g);
         if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(ComputeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
            AddRhoCCP(x,sig,lbl,mypos);
            if (mypos>=0) {
               AddRCP2ConCCP(mypos,i);
               AddRCP2ConCCP(mypos,j);
            } 
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nRCP>0)? (nRCP) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   for (int i=0; i<nCCP; i++) {
      RemoveRedundInLabel(lblCCP[i]);
      //cout << lblCCP[i] << endl;
   }
   cout << nCCP << " CCPs found.\n";
   //ScreenUtils::DisplayWarningMessage("SetRhoCCPs(...) under construction.");
   return true;
}
bool CritPtNetWork::SetLOLACPs() {
   double x[3],rho,g[3],magg;
   string lbl;
   int sig;
   bool tmpbool=false;
   for (int i=0; i<(bn->nNuc); i++) {
      for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
      SeekLOLACP(x,rho,g,sig);
      magg=ComputeMagnitudeV3(g);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSLOLACPGRADMAG) ) {
         lbl=bn->atLbl[i];
         AddRhoACP(x,sig,lbl);
      } else {
         for ( int j=0 ; j<3 ; j++ ) {x[j]=bn->R[i][j];}
         lbl=bn->atLbl[i]+"+";
         AddRhoACP(x,sig,lbl);
         tmpbool=true;
      }
   }
   if ( tmpbool ) {
      ScreenUtils::DisplayWarningMessage("Some ACPs presented |nabla rho|!= 0... Look for '+' in labels.");
   }
   if ( nACP!=(bn->nNuc) ) {
      ScreenUtils::DisplayWarningMessage("Number of regular ACPs is different from number of Nuclei!");
   }
   int regnACP=nACP;
   cout << "Found " << regnACP << " regular ACPs (" << bn->nNuc << " nuc)" << endl;
   cout << "Looking around nuclei..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double xs[3],rad;
   for ( int i=0 ; i<(bn->nNuc) ; i++ ) {
         for ( int k=0 ; k<3 ; k++ ) {xs[k]=(bn->R[i][k]);}
         rad=ComputeMagnitudeV3(xs);
         lbl="NNLOLACP"+StringTools::GetStringFromInt(i);
         rad*=0.01e0;
         SeekLOLACPsAroundAPoint(xs,rad,lbl,-1);
#if USEPROGRESSBAR
         ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
                  double((bn->nNuc > 1)? (bn->nNuc-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   /* -------------------------------------------------------------------------------  */
   //*
   cout << "Looking between atom pairs..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   //double xs[3],rad;
   for ( int i=0 ; i<(bn->nNuc) ; i++ ) {
      for ( int j=(i+1) ; j<(bn->nNuc) ; j++ ) {
         for ( int k=0 ; k<3 ; k++ ) {xs[k]=(bn->R[i][k]-bn->R[j][k]);}
         rad=ComputeMagnitudeV3(xs);
         if ( rad<=(2.0e0*bn->maxBondDist) ) {
            for ( int k=0 ; k<3 ; k++ ) {xs[k]=0.5e0*(bn->R[i][k]+bn->R[j][k]);}
            lbl="NNLOLACP"+StringTools::GetStringFromInt((i+1));
            rad*=0.1e0;
            SeekLOLACPsAroundAPoint(xs,rad,lbl,-1);
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((bn->nNuc > 1)? (bn->nNuc-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
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
bool CritPtNetWork::SetLOLBCPs(void) {
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   int ata,atb,k,sig;
   int mypos;
   double x[3],g[3],rho;
   string lbl;
   for (int i=0; i<(bn->nNuc); i++) {
      ata=i;
      atb=bn->bNet[ata][0];
      k=1;
      while ((k<MAXBONDINGATOMS)&&(atb>0)) {
         for (int j=0; j<3; j++) {x[j]=0.5e0*(bn->R[ata][j]+bn->R[atb][j]);}
         SeekLOLBCP(x,rho,g,sig);
         if (ata<atb) {
            lbl=bn->atLbl[ata]+string("-")+bn->atLbl[atb];
         } else {
            lbl=bn->atLbl[atb]+string("-")+bn->atLbl[ata];
         }
         if (AddRhoBCP(x,sig,lbl,mypos)&&mypos>=0) {
            conBCP[mypos][0]=ata;
            conBCP[mypos][1]=atb;
         }
         atb=bn->bNet[ata][k];
         ++k;
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nACP>1)? (nACP-1) : 1)));
#endif
   }
   normalbcp=nBCP;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   /* ------------------------------------------------------ */
   cout << nBCP << " BCPs found.\n";
   string ll;
   cout << "Including non-nuclear attractors..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   for ( int i=(bn->nNuc) ; i<nACP ; i++ ) {
      SeekLOLBCPWithExtraACP(i,(bn->maxBondDist));
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/\
               double((nACP>1)? (nACP-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   /*
   cout << "Looking for possible extra BCPs (Hydrogen bonds)..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double extbd=(bn->maxBondDist*CPNW_EXTENDEDBONDDISTFACTOR),dd=0.0e0;
   for (int i=0; i<nACP; i++) {
      for (int j=(i+1); j<nACP; j++) {
         for (int k=0; k<3; k++) {x[k]=((RACP[i][k]-RACP[j][k]));}
         dd=ComputeMagnitudeV3(x);
         if ((dd>=(bn->maxBondDist*0.9e0))&&(dd<=extbd)) {
            for (int k=0; k<3; k++) {x[k]=0.5e0*(RACP[i][k]+RACP[j][k]);}
            if ((wf->EvalDensity(x[0],x[1],x[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
               SeekRhoBCP(x,rho,g,sig);
               if (i<j) {ata=i; atb=j;} else {ata=j; atb=i;}
               lbl=string("*")+lblACP[ata]+string("-")+lblACP[atb];
               wf->EvalRhoGradRho(x[0],x[1],x[2],rho,g);
               if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(ComputeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
                  AddRhoBCP(x,sig,lbl,mypos);
                  if ((mypos>=0)&&(mypos<dBCP)) {
                     conBCP[mypos][0]=ata;
                     conBCP[mypos][1]=atb;
                  }
               }
            }
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((nACP-1))));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   if (normalbcp==nBCP) {
      cout << "No more BCPs found." << endl;
   } else {
      double tbcp[3];
      for ( int i=normalbcp ; i<nBCP ; i++ ) {
         for ( int k=0 ; k<3 ; k++ ) {tbcp[k]=RBCP[i][k];}
         FindTwoClosestACPs(tbcp,ata,atb);
         conBCP[i][0]=ata;
         conBCP[i][1]=atb;
         lblBCP[i]="*"+lblACP[ata]+"-"+lblACP[atb];
      }
      cout << (nBCP-normalbcp) << " new BCP";
      if ((nBCP-normalbcp)>1) {cout << "s";}
      cout << " found! Total number of BCPs: " << nBCP << endl;
   }
   // */
   return true;
}
bool CritPtNetWork::SetLOLRCPs(void) {
   ScreenUtils::DisplayWarningMessage("No LOL RCP will be look for...");
   ScreenUtils::DisplayWarningMessage("(CritPtNetWork::SetLOLRCPs(...) under construction.)");
   return false;
}
bool CritPtNetWork::SetLOLCCPs(void) {
   ScreenUtils::DisplayWarningMessage("No LOL CCP will be look for...");
   ScreenUtils::DisplayWarningMessage("(CritPtNetWork::SetLOLCCPs(...) under construction.)");
   return false;
}
bool CritPtNetWork::AddRhoACP(double (&x)[3],int sig,string &lbl) {
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
   if (ImNew(x,dACP,RACP,pos)) {
      for (int i=0; i<3; i++) {RACP[pos][i]=x[i];}
      lblACP[pos]=lbl;
      //(lbl[lbl.length()-1])++;
      ++nACP;
      return true;
   }
   return false;
}
bool CritPtNetWork::AddRhoBCP(double (&x)[3],int sig,string &lbl,int &pos) {
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
   if (ImNew(x,dBCP,RBCP,ttpos)) {
      for (int i=0; i<3; i++) {RBCP[ttpos][i]=x[i];}
      lblBCP[ttpos]=lbl;
      ++nBCP;
      pos=int(ttpos);
      return true;
   }
   return false;
}
bool CritPtNetWork::AddRhoRCP(double (&x)[3],int sig,string &lbl,int &pos) {
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
   if (ImNew(x,dRCP,RRCP,ttpos)) {
      for (int i=0; i<3; i++) {RRCP[ttpos][i]=x[i];}
      lblRCP[ttpos]=lbl;
      ++nRCP;
      pos=int(ttpos);
      return true;
   } else {
      if ((int(ttpos)<nRCP)&&(ttpos!=string::npos)) {
         lblRCP[ttpos]+=(string("-")+lbl);
         pos=int(ttpos);
      }
   }
   return false;
}
bool CritPtNetWork::AddRhoCCP(double (&x)[3],int sig,string &lbl,int &pos) {
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
   if (ImNew(x,dCCP,RCCP,ttpos)) {
      for (int i=0; i<3; i++) {RCCP[ttpos][i]=x[i];}
      lblCCP[ttpos]=lbl;
      ++nCCP;
      pos=int(ttpos);
      return true;
   } else {
      if ((int(ttpos)<nCCP)&&(ttpos!=string::npos)) {
         lblCCP[ttpos]+=(string("-")+lbl);
         pos=int(ttpos);
      }
   }
   return false;
}
void CritPtNetWork::SeekRhoACPsAroundAPoint(double const (&oo)[3],double const ddxx,\
      string const &blbl,int uunvrt) {
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   double xx[3],rho,gg[3],magg;
   int sig;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      SeekRhoACP(xx,rho,gg,sig);
      magg=ComputeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (AddRhoACP(xx,sig,lbl)) {++lbl[(lbl.length()-1)];}
      }
   }
}
void CritPtNetWork::SeekRhoBCPsAroundAPoint(double const (&oo)[3],double const ddxx,\
      string const &blbl,int uunvrt) {
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   double xx[3],rho,gg[3],magg;
   int sig,pos;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      SeekRhoBCP(xx,rho,gg,sig);
      magg=ComputeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (AddRhoBCP(xx,sig,lbl,pos)) {++lbl[(lbl.length()-1)];}
      }
   }
}
void CritPtNetWork::SeekRhoRCPsAroundAPoint(double const (&oo)[3],double const ddxx,\
      string const &blbl,int uunvrt) {
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   double xx[3],rho,gg[3],magg;
   int sig,pos;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      SeekRhoRCP(xx,rho,gg,sig);
      magg=ComputeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (AddRhoRCP(xx,sig,lbl,pos)) {++lbl[(lbl.length()-1)];}
      }
   }
}
void CritPtNetWork::SeekRhoCCPsAroundAPoint(double const (&oo)[3],double const ddxx,\
      string const &blbl,int uunvrt) {
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   double xx[3],rho,gg[3],magg;
   int sig,pos;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      SeekRhoCCP(xx,rho,gg,sig);
      magg=ComputeMagnitudeV3(gg);
      if ( (rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(magg<CPNW_EPSRHOACPGRADMAG) ) {
         if (AddRhoCCP(xx,sig,lbl,pos)) {++lbl[(lbl.length()-1)];}
      }
   }
}
void CritPtNetWork::SeekRhoBCPWithExtraACP(int acppos,double maxrad) {
   double xx[3],rho,g[3];
   int sig;
   string lbl="";
   int mypos=-1;
   for ( int i=0 ; i<nACP ; i++ ) {
      if ( i==acppos ) {continue;}
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=(RACP[acppos][k]-RACP[i][k]);}
      if ( ComputeMagnitudeV3(xx)<=maxrad ) {
         for (int k=0; k<3; k++) {xx[k]=0.5e0*(RACP[i][k]+RACP[acppos][k]);}
         if ((wf->EvalDensity(xx[0],xx[1],xx[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
            SeekRhoBCP(xx,rho,g,sig);
            if (i<acppos) {
               lbl=string("*")+lblACP[i]+string("-")+lblACP[acppos];
            } else {
               lbl=string("*")+lblACP[acppos]+string("-")+lblACP[i];
            }
            wf->EvalRhoGradRho(xx[0],xx[1],xx[2],rho,g);
            if ((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&
                  (ComputeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG)) {
               AddRhoBCP(xx,sig,lbl,mypos);
               if ((mypos>=0)&&(mypos<dBCP)) {
                  conBCP[mypos][0]=i;
                  conBCP[mypos][1]=acppos;
               }
            }
         }
      }
   }
}
void CritPtNetWork::SeekLOLBCPWithExtraACP(int acppos,double maxrad) {
   double xx[3],lol,g[3],hl[3][3];
   int sig;
   string lbl="";
   int mypos=-1;
   for ( int i=0 ; i<nACP ; i++ ) {
      if ( i==acppos ) {continue;}
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=(RACP[acppos][k]-RACP[i][k]);}
      if ( ComputeMagnitudeV3(xx)<=maxrad ) {
         for (int k=0; k<3; k++) {xx[k]=0.5e0*(RACP[i][k]+RACP[acppos][k]);}
         if ((wf->EvalDensity(xx[0],xx[1],xx[2])>CPNW_MINRHOSIGNIFICATIVEVAL)) {
            SeekLOLBCP(xx,lol,g,sig);
            if (i<acppos) {
               lbl=string("*")+lblACP[i]+string("-")+lblACP[acppos];
            } else {
               lbl=string("*")+lblACP[acppos]+string("-")+lblACP[i];
            }
            wf->EvalHessLOL(xx,lol,g,hl);
            if (ComputeMagnitudeV3(g)<CPNW_EPSLOLACPGRADMAG) {
               AddRhoBCP(xx,sig,lbl,mypos);
               if ((mypos>=0)&&(mypos<dBCP)) {
                  conBCP[mypos][0]=i;
                  conBCP[mypos][1]=acppos;
               }
            }
         }
      }
   }
}
void CritPtNetWork::SeekLOLACPsAroundAPoint(double const (&oo)[3],double const ddxx,\
      string const &blbl,int uunvrt) {
   int nvrt=uunvrt;
   if ( uunvrt>nIHV || uunvrt<1 ) {nvrt=nIHV;}
   double xx[3],rho,gg[3],magg;
   int sig;
   string lbl=blbl;
   lbl.append("1");
   for ( int i=0 ; i<nvrt ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xx[k]=oo[k]+IHV[i][k]*ddxx;}
      SeekLOLACP(xx,rho,gg,sig);
      magg=ComputeMagnitudeV3(gg);
      if ( /*(rho>CPNW_MINLOLSIGNIFICATIVEVAL)&&*/(magg<CPNW_EPSLOLACPGRADMAG) ) {
         if (AddRhoACP(xx,sig,lbl)) {++lbl[(lbl.length()-1)];}
      }
   }
}
void CritPtNetWork::SeekLOLACPsOnASphere(int atIdx,int nDivR,int nDivT,int nDivP,\
      double radmin,double radmax) {
   static const double ppii=4.0e0*atan(1.0e0);
   if ( atIdx<0 || atIdx>=bn->nNuc ) {ScreenUtils::DisplayErrorMessage("Requesting a non existent atom!");}
   if ( nDivR<=0 || nDivT<=0 || nDivP<=0 ) {
      ScreenUtils::DisplayErrorMessage("number of points must be greater than zero!");
   }
   cout << "Looking for LOL ACPs around atom " << (atIdx+1) << "(" 
        << wf->atLbl[atIdx] << ")" << endl;
   double xx[3],lol,gl[3];
   int sig,nseeds,count;
   string blbl="LOLACPSp"+StringTools::GetStringFromInt(atIdx+1)+"*",lbl;
   double cc[3],dr,dt,dp,currr,currt,currp;
   dr=(radmax-radmin)/double(nDivR-1);
   dt=ppii/double(nDivT-1);
   dp=2.0e0*ppii/double(nDivP);
   nseeds=nDivR*nDivT*nDivP;
   count=0;
   int origacp=nACP,idxLbl=1;
   lbl=blbl+StringTools::GetStringFromInt(idxLbl);
   for ( int i=0 ; i<3 ; i++ ) {cc[i]=bn->R[atIdx][i];}
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   for ( int ir=0 ; ir<nDivR ; ir++ ) {
      currr=radmin+double(ir)*dr;
      for ( int it=0 ; it<nDivT ; it++ ) {
         currt=(dt*double(it));
         for ( int ip=0 ; ip<nDivP ; ip++ ) {
            currp=dp*double(ip);
            xx[0]=cc[0]+currr*sin(currt)*cos(currp);
            xx[1]=cc[1]+currr*sin(currt)*sin(currp);
            xx[2]=cc[2]+currr*cos(currt);
            SeekLOLACP(xx,lol,gl,sig);
            if ( ComputeMagnitudeV3(gl)<CPNW_EPSLOLACPGRADMAG ) {
               if (AddRhoACP(xx,sig,lbl)){lbl=blbl+StringTools::GetStringFromInt(++idxLbl);}
            }
            ++count;
         }
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(count)/\
               double((nseeds>1) ? (nseeds-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   if ( origacp<nACP ) {cout << "Found " << (nACP-origacp) << " new ACPs!..." << endl;}
}
void CritPtNetWork::ExtendedSearchCPs(void) {
   double xs[3],dx;
   dx=bn->maxBondDist*0.1e0;
   string lbl="";
   int count=0;
   int initacp=nACP,initbcp=nBCP,initrcp=nRCP,initccp=nCCP;
   cout << "Looking around extraBCPs..." << endl;
   for ( int i=normalbcp ; i<nBCP ; i++ ) {
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
      for ( int k=0 ; k<3 ; k++ ) {xs[k]=RBCP[i][k];}
      lbl="extACPb";
      SeekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extRCPb";
      SeekRhoRCPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extCCPb";
      SeekRhoCCPsAroundAPoint(xs,dx,lbl,nIHV);
      ++count;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(count)/\
               double(((nBCP-normalbcp)>1)? (nBCP-normalbcp-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   count=0;
   cout << "Looking around RCPs..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   for ( int i=0 ; i<nRCP ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xs[k]=RRCP[i][k];}
      lbl="extACPr";
      SeekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extBCPr";
      SeekRhoBCPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extCCPr";
      SeekRhoCCPsAroundAPoint(xs,dx,lbl,nIHV);
      ++count;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(count)/\
               double((nRCP>1)? (nRCP-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   count=0;
   cout << "Looking around CCPs..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   for ( int i=0 ; i<nCCP ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {xs[k]=RCCP[i][k];}
      lbl="extACPc";
      SeekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extBCPc";
      SeekRhoBCPsAroundAPoint(xs,dx,lbl,nIHV);
      lbl="extRCPc";
      SeekRhoRCPsAroundAPoint(xs,dx,lbl,nIHV);
      ++count;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(count)/\
               double(((nBCP-normalbcp)>1)? (nBCP-normalbcp-1) : 1)));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   for (int i=0; i<nRCP; i++) {RemoveRedundInLabel(lblRCP[i]);}
   for (int i=0; i<nCCP; i++) {RemoveRedundInLabel(lblCCP[i]);}
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
      ScreenUtils::PrintBetweenStarLines(string("nACP-nBCP+nRCP-nCCP = "+\
               StringTools::GetStringFromInt(nACP-nBCP+nRCP-nCCP)));
   }
}
void CritPtNetWork::CustomSearchTwoACPs(int acpIdx1,int acpIdx2) {
   cout << "Looking around a custom seed (two acps)..." << endl;
   string str1=string("This acp does not exist! There are ");
   str1+=StringTools::GetStringFromInt(nACP);
   str1+=" ACPs in this molecule.";
   string str2="Requested acp: ";
   string str3=". Nothing to do...";
   if ( acpIdx1>=nACP ) {
      ScreenUtils::DisplayWarningMessage(str1);
      ScreenUtils::DisplayWarningMessage(str2+StringTools::GetStringFromInt(acpIdx1+1)+str3);
      return;
   }
   if ( acpIdx2>=nACP ) {
      ScreenUtils::DisplayWarningMessage(str1);
      ScreenUtils::DisplayWarningMessage(str2+StringTools::GetStringFromInt(acpIdx2+1)+str3);
      return;
   }
   double xs[3];
   for ( int i=0 ; i<3 ; ++i ) {
      xs[i]=0.5e0*(RACP[acpIdx1][i]+RACP[acpIdx2][i]);
   }
   CustomSearchCPs(xs);
}
void CritPtNetWork::CustomSearchCPs(double (&xs)[3]) {
   double dx;
   dx=bn->maxBondDist*0.05e0;
   string lbl="";
   int initacp=nACP,initbcp=nBCP,initrcp=nRCP,initccp=nCCP;
   cout << "Looking for ACPs..." << endl;
   lbl="custACP";
   SeekRhoACPsAroundAPoint(xs,dx,lbl,nIHV);
   cout << "Looking for BCPs..." << endl;
   lbl="custBCP";
   SeekRhoBCPsAroundAPoint(xs,dx,lbl,nIHV);
   cout << "Looking for RCPs..." << endl;
   lbl="custRCP";
   SeekRhoRCPsAroundAPoint(xs,dx,lbl,nIHV);
   cout << "Looking for CCPs..." << endl;
   lbl="custCCP";
   SeekRhoCCPsAroundAPoint(xs,dx,lbl,nIHV);
   for (int i=0; i<nRCP; i++) {RemoveRedundInLabel(lblRCP[i]);}
   for (int i=0; i<nCCP; i++) {RemoveRedundInLabel(lblCCP[i]);}
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
      ScreenUtils::PrintBetweenStarLines(string("nACP-nBCP+nRCP-nCCP = "+\
               StringTools::GetStringFromInt(nACP-nBCP+nRCP-nCCP)));
   }
}
void CritPtNetWork::RemoveRedundInLabel(string &lbl) {
   string tl=lbl,chunk,tchk,tl2,fl="";
   chunk=GetFirstChunkOfLabel(tl);
   while (chunk.length()>0) {
      fl+=(chunk+"-");
      tl.erase(0,(chunk.length()+1));
      tchk=GetFirstChunkOfLabel(tl);
      tl2="";
      while (tchk.length()>0) {
         if (tchk!=chunk) {
            tl2+=(tchk+"-");
         }
         tl.erase(0,(tchk.length()+1));
         tchk=GetFirstChunkOfLabel(tl);
      }
      tl=tl2;
      chunk=GetFirstChunkOfLabel(tl);
   }
   size_t clen=fl.length()-1;
   if (fl[clen]=='-') {fl.erase(clen,1);}
   //cout << fl << endl;
   lbl=fl;
   return;
}
string CritPtNetWork::GetFirstChunkOfLabel(string &lbl) {
   if (lbl.length()==0) {
      return string("");
   }
   size_t pos=lbl.find_first_of("-");
   return lbl.substr(0,pos);
}
void CritPtNetWork::PrintAllFieldProperties(double x,double y,double z) {
   wf->DisplayAllFieldProperties(x,y,z);
}
void CritPtNetWork::WriteAllFieldProperties(double x,double y,double z,ofstream &ofil) {
   wf->WriteAllFieldProperties(x,y,z,ofil);
}
void CritPtNetWork::GetACPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig) {
   double eive[3][3],b[3],F[3];
   EigenDecompositionJAMA::EigenDecomposition3(hess, eive, b);
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   double h4[4][4],m4[4][4],v4[4];
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
   EigenDecompositionJAMA::EigenDecomposition4(h4, m4, v4);
   double lp=v4[3];
   if ( fabs(lp)<CPNW_EPSEIGENVALUECPSEARCH ) {lp=CPNW_EPSEIGENVALUECPSEARCH;}
#if DEBUG
   if (lp<=0.0e0) {
      ScreenUtils::DisplayWarningMessage(string("lp<=0!: "+StringTools::GetStringFromReal(lp)));
      for ( int i=0 ; i<4 ; i++ ) {cout << v4[i] << " ";}
      cout << endl;
      ScreenUtils::PrintM3x3Comp("hess:\n",hess);
   }
#endif
   double dd;
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      dd=b[j]-lp; if ( fabs(dd) < 1.0e-15 ) { dd+=1.0e-15; }
      hh[0]-=eive[0][j]*F[j]/dd;
      hh[1]-=eive[1][j]*F[j]/dd;
      hh[2]-=eive[2][j]*F[j]/dd;
   }
   dd=ComputeMagnitudeV3(hh);
   if (dd>stepSizeACP) {
      for (int i=0; i<3; i++) {
         hh[i]*=(stepSizeACP/dd);
      }
   }
   sig=ComputeSignature(b);
   return;
}
void CritPtNetWork::GetBCPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig) {
   double eive[3][3],b[3];
   EigenDecompositionJAMA::EigenDecomposition3(hess, eive, b);
   double F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   double h3[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {h3[i][j]=0.00000e0;}
   }
   double ln=0.5e0*(b[2]-sqrt(b[2]*b[2]+4.0e0*F[2]*F[2]));
   h3[0][0]=b[0];
   h3[1][1]=b[1];
   h3[2][0]=h3[0][2]=F[0];
   h3[2][1]=h3[1][2]=F[1];
   double m3[3][3],vv[3];
   EigenDecompositionJAMA::EigenDecomposition3(h3, m3, vv);
   double dd,lp=vv[2];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      dd=b[j]-lp; if ( fabs(dd) < 1.0e-15 ) { dd+=1.0e-15; }
      hh[0]-=eive[0][j]*F[j]/dd;
      hh[1]-=eive[1][j]*F[j]/dd;
      dd=b[j]-ln; if ( fabs(dd) < 1.0e-15 ) { dd+=1.0e-15; }
      hh[2]-=eive[2][j]*F[j]/dd;
   }
   dd=ComputeMagnitudeV3(hh);
   if (dd>stepSizeBCP) {
      for (int i=0; i<3; i++) {
         hh[i]*=(stepSizeBCP/dd);
      }
   }
   sig=ComputeSignature(b);
   return;
}
void CritPtNetWork::GetRCPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig) {
   double eive[3][3],b[3];
   EigenDecompositionJAMA::EigenDecomposition3(hess, eive, b);
   double F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   double h3[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {h3[i][j]=0.00000e0;}
   }
   double lp=0.5e0*(b[0]+sqrt(b[0]*b[0]+4.0e0*F[0]*F[0]));
   h3[0][0]=b[1];
   h3[1][1]=b[2];
   h3[2][0]=h3[0][2]=F[1];
   h3[2][1]=h3[1][2]=F[2];
   double m3[3][3],vv[3];
   EigenDecompositionJAMA::EigenDecomposition3(h3, m3, vv);
   double dd,ln=vv[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      dd=b[j]-lp; if ( fabs(dd) < 1.0e-15 ) { dd+=1.0e-15; }
      hh[0]-=eive[0][j]*F[j]/dd;
      dd=b[j]-ln; if ( fabs(dd) < 1.0e-15 ) { dd+=1.0e-15; }
      hh[1]-=eive[1][j]*F[j]/dd;
      hh[2]-=eive[2][j]*F[j]/dd;
   }
   dd=ComputeMagnitudeV3(hh);
   if (dd>stepSizeRCP) {
      for (int i=0; i<3; i++) {
         hh[i]*=(stepSizeRCP/dd);
      }
   }
   sig=ComputeSignature(b);
   return;
}
void CritPtNetWork::GetCCPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig) {
   double eive[3][3],b[3];
   EigenDecompositionJAMA::EigenDecomposition3(hess, eive, b);
   double F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   double h4[4][4],m4[4][4],v4[4];
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
   EigenDecompositionJAMA::EigenDecomposition4(h4, m4, v4);
   double dd,ln=v4[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      dd=b[j]-ln; if ( fabs(dd) < 1.0e-15 ) { dd+=1.0e-15; }
      hh[0]-=eive[0][j]*F[j]/dd;
      hh[1]-=eive[1][j]*F[j]/dd;
      hh[2]-=eive[2][j]*F[j]/dd;
   }
   dd=ComputeMagnitudeV3(hh);
   if (dd>stepSizeCCP) {
      for (int i=0; i<3; i++) {
         hh[i]*=(stepSizeCCP/dd);
      }
   }
   sig=ComputeSignature(b);
   return;
}
int CritPtNetWork::ComputeSignature(double (&ev)[3]) {
   int res=0;
   for ( int i=0 ; i<3 ; i++ ) {
      res+=SIGNF(ev[i]);
      if ( ev[i]==0.0e0 ) {
         res=0;
#if DEBUG
         ScreenUtils::DisplayWarningMessage("Zero signature!");
         ScreenUtils::PrintV3Comp("EigenValues: ",ev);
#endif
         break;
      }
   }
   return res;
}
void CritPtNetWork::SeekRhoACP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig) {
   double rho,gr[3],hr[3][3],dx[3];
   wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=ComputeSignature(hr);
   double magd=ComputeMagnitudeV3(gr);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-3 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItACP)) {
      GetACPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
void CritPtNetWork::SeekRhoBCP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig) {
   double rho,gr[3],hr[3][3],dx[3];
   wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=ComputeSignature(hr);
   double magd=ComputeMagnitudeV3(gr);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-1 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItBCP)) {
      GetBCPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
void CritPtNetWork::SeekRhoRCP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig) {
   double rho,gr[3],hr[3][3],dx[3];
   wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=ComputeSignature(hr);
   double magd=ComputeMagnitudeV3(gr);
   if ( magd<=CPNW_EPSRHOACPGRADMAG && sig==1 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItRCP)) {
      GetRCPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
void CritPtNetWork::SeekRhoCCP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig) {
   double rho,gr[3],hr[3][3],dx[3];
   wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
   sig=ComputeSignature(hr);
   double magd=ComputeMagnitudeV3(gr);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==3 ) {
      rho2ret=rho;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gr[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItCCP)) {
      GetCCPStep(gr,hr,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=dx[i];}
      wf->EvalHessian(x[0],x[1],x[2],rho,gr,hr);
      magd=sqrt(gr[0]*gr[0]+gr[1]*gr[1]+gr[2]*gr[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for ( int k=0 ; k<3 ; k++ ) {g[k]=gr[k];}
   rho2ret=rho;
   return;
}
void CritPtNetWork::SeekLOLACP(double (&x)[3],double &ll,double (&g)[3],int &sig) {
   double lol,gl[3],hl[3][3],dx[3];
   wf->EvalHessLOL(x,lol,gl,hl);
   sig=ComputeSignature(hl);
   double magd=ComputeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-3 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItACP)) {
      GetACPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->EvalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
void CritPtNetWork::SeekLOLBCP(double (&x)[3],double &ll,double (&g)[3],int &sig) {
   double lol,gl[3],hl[3][3],dx[3];
   wf->EvalHessLOL(x,lol,gl,hl);
   sig=ComputeSignature(hl);
   double magd=ComputeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==-1 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItBCP)) {
      GetBCPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->EvalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
void CritPtNetWork::SeekLOLRCP(double (&x)[3],double &ll,double (&g)[3],int &sig) {
   double lol,gl[3],hl[3][3],dx[3];
   wf->EvalHessLOL(x,lol,gl,hl);
   sig=ComputeSignature(hl);
   double magd=ComputeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==1 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItRCP)) {
      GetRCPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->EvalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
void CritPtNetWork::SeekLOLCCP(double (&x)[3],double &ll,double (&g)[3],int &sig) {
   double lol,gl[3],hl[3][3],dx[3];
   wf->EvalHessLOL(x,lol,gl,hl);
   sig=ComputeSignature(hl);
   double magd=ComputeMagnitudeV3(gl);
   if ( magd<CPNW_EPSRHOACPGRADMAG && sig==3 ) {
      ll=lol;
      for ( int i=0 ; i<3 ; i++ ) {g[i]=gl[i];}
      return;
   }
   double magh=magd;
   int count=0;
   while (((magd>CPNW_EPSRHOACPGRADMAG)&&(magh>CPNW_EPSFABSDIFFCOORD))&&(count<maxItRCP)) {
      GetCCPStep(gl,hl,dx,sig);
      for (int i=0; i<3; i++) {x[i]+=(dx[i]/*MAXSTEPSIZEACPLOLSEARCH*/);}
      wf->EvalHessLOL(x,lol,gl,hl);
      magd=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      magh=sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
      ++count;
   }
   for (int i=0; i<3; i++) {g[i]=gl[i];}
   ll=lol;
   return;
}
int CritPtNetWork::ComputeSignature(double (&hh)[3][3]) {
   double eive[3][3],b[3];
   EigenDecompositionJAMA::EigenDecomposition3(hh, eive, b);
   return ComputeSignature(b);
}
bool CritPtNetWork::ImNew(double (&x)[3],int dim,double ** (&arr),size_t &pos) {
   double ee;
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
void CritPtNetWork::DisplayXCPCoords(char cpt) {
   bool cpknow=false;
   string cplbl="";
   string longcplbl="";
   string *lblArray=NULL;
   int theNCP=0;
   double **RArray=NULL;
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
      ScreenUtils::DisplayErrorMessage(string("First seek the "+cplbl+"s.\nNo "\
               +cplbl+" to display its coordinates."));
      return;
   }
   cout << scientific << setprecision(12);
   ScreenUtils::PrintScrCharLine('+');
   ScreenUtils::CenterString(string("Coordinates of "+longcplbl+" Critical Points"));
   ScreenUtils::PrintScrCharLine('-');
   for (int i=0; i<theNCP; i++) {
      cout << cplbl << "[" << (i+1) << "](" << lblArray[i] << "): ";
      for (int j=0; j<3; j++) {cout << RArray[i][j] << " ";}
      cout << endl;
   }
   ScreenUtils::PrintScrCharLine('+');
   cout.unsetf(std::ios::scientific);
   return;
}
void CritPtNetWork::DisplayAllCPCoords(void) {
   DisplayXCPCoords('a');
   DisplayXCPCoords('b');
   DisplayXCPCoords('r');
   DisplayXCPCoords('c');
}
void CritPtNetWork::DisplayIHVCoords(void) {
   cout << scientific << setprecision(8);
   ScreenUtils::PrintScrCharLine('+');
   ScreenUtils::CenterString("Coordinates of Icosahedron Vertices");
   ScreenUtils::PrintScrCharLine('-');
   for (int i=0; i<nIHV; i++) {
      cout << "V[" << (i+1) << "]: ";
      for (int j=0; j<3; j++) {cout << IHV[i][j] << " ";}
      cout << endl;
   }
   ScreenUtils::PrintScrCharLine('+');
   cout.unsetf(std::ios::scientific);
   return;
}
void CritPtNetWork::FindTwoClosestAtoms(double (&xo)[3],int &idx1st,int &idx2nd) {
   if ( wf->nNuc<2 ) {idx1st=0; idx2nd=0; return;}
   double xmagt=0.0e0,xmag1=0.0e0,xmag2=0.0e0;
   int ii1,ii2,iit;
   for ( int k=0 ; k<3 ; k++ ) {xmag1+=((xo[k]-wf->GetR(0,k))*(xo[k]-wf->GetR(0,k)));}
   ii1=0;
   for ( int k=0 ; k<3 ; k++ ) {xmag2+=((xo[k]-wf->GetR(1,k))*(xo[k]-wf->GetR(1,k)));}
   ii2=1;
   if ( xmag1>xmag2 ) {
      xmagt=xmag1; xmag1=xmag2; xmag2=xmagt;
      ii1=1;
      ii2=0;
   }
   for ( int i=2 ; i<wf->nNuc ; i++ ) {
      xmagt=0.0e0;
      for ( int k=0 ; k<3 ; k++ ) {xmagt+=((xo[k]-wf->GetR(i,k))*(xo[k]-wf->GetR(i,k)));}
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
      ScreenUtils::DisplayWarningMessage("Identical atoms!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
#endif /* ( DEBUG ) */
   return;
}
void CritPtNetWork::FindTwoClosestAtomsToBCP(const int bcpIdx,int &at1Idx,int&at2Idx) {
   if ( wf->nNuc<2 ) {at1Idx=0; at2Idx=0; return;}
#if DEBUG
   if ( !iknowbcps ) {
      ScreenUtils::DisplayErrorMessage("First look for BCPs!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      at1Idx=at2Idx=0;
      return;
   }
#endif /* ( DEBUG ) */
   if ( bcpIdx>=nBCP || bcpIdx<0 ) {
      ScreenUtils::DisplayErrorMessage("Requesting a non-existing BCP!");
      at1Idx=at2Idx=0;
      return;
   }
   double xo[3];
   for ( int i=0 ; i<3 ; ++i ) { xo[i]=RBCP[bcpIdx][i]; }
   return FindTwoClosestAtoms(xo,at1Idx,at2Idx);
}
void CritPtNetWork::FindTwoClosestACPs(double (&xo)[3],int &idx1st,int &idx2nd) {
   if ( nACP<2 ) {idx1st=0; idx2nd=0; return;}
   double xmagt=0.0e0,xmag1=0.0e0,xmag2=0.0e0;
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
      ScreenUtils::DisplayWarningMessage("Identical ACPs!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
#endif /* ( DEBUG ) */
   return;
}
void CritPtNetWork::InvertOrderBGPPoints(int dim) {
   int numop=((dim+1)>>1);

   double tmp;
   for (int i=0; i<numop; i++) {
      for (int j=0; j<3; j++) {
         tmp=RGP[i][j];
         RGP[i][j]=RGP[dim-i-1][j];
         RGP[dim-i-1][j]=tmp;
      }
   }
}
void CritPtNetWork::InvertOrderBGPPoints(int dim,double** (&arr)) {
#if DEBUG
   if (arr==NULL) {
      ScreenUtils::DisplayErrorMessage("The pointer for this array is null!!!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return;
   }
#endif
   int numop=((dim+1)>>1);
   double tmp;
   for (int i=0; i<numop; i++) {
      for (int j=0; j<3; j++) {
         tmp=arr[i][j];
         arr[i][j]=arr[dim-i-1][j];
         arr[dim-i-1][j]=tmp;
      }
   }
}
void CritPtNetWork::WriteCPProps(string &ofnam,string &wfnnam) {
   ofstream ofil;
   ofil.open(ofnam.c_str());
   double x,y,z;
   //double gx,gy,gz,rho,hxx,hyy,hzz,hxy,hxz,hyz;
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
   FileUtils::WriteScrStarLine(ofil,false);
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
         WriteAllFieldProperties(x,y,z,ofil);
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
         WriteAllFieldProperties(x,y,z,ofil);
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
         WriteAllFieldProperties(x,y,z,ofil);
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
         WriteAllFieldProperties(x,y,z,ofil);
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
bool CritPtNetWork::MakePOVFile(string pnam,POVRayConfiguration &pvp,int campos) {
   ofstream pof;
   pof.open(pnam.c_str(),std::ios::out);
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
   FileUtils::WriteScrCharLine(pof,'/',false);
   pof << "//" << endl;
#if DEBUG
   pof << "//Code generated by the class CritPtNetWork." << endl;
#endif
   pof << "//Below you can find some options to be parsed to povray" << endl;
   pof << "//set your custom values." << endl;
   pof << "//You can reconstruct the image using the script dtkpov2png" << endl;
   pof << "//" << endl;
   FileUtils::WriteScrCharLine(pof,'/',false);
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
   double allcprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
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
   pof << "#declare ColorABGradPath=rgb <0.0,0.2,1.0>;" << endl;
   pof << "#declare ColorARGradPath=rgb <0.0,0.8,0.0>;" << endl;
   pof << "#declare ColorACGradPath=rgb <1.0,0.5,0.0>;" << endl;
   pof << "#declare TransmitAtomSphere=0.7;" << endl;
   pof << "#declare TransmitStdBondCylinder=0.7;" << endl;
   pof << "#default { finish { specular 0.3 roughness 0.03 phong .1 } }" << endl;
   FileUtils::WriteScrCharLine(pof,'/',false);
   pof << "//For the colors, instead of rgb <...>, you may want to try Red, Yellow, ..." << endl;
   pof << "//  or any of the colors defined in \"colors.inc\"" << endl;
   pof << "//" << endl;
   FileUtils::WriteScrCharLine(pof,'/',false);
   pof << "// END OF CUSTOM OPTIONS" << endl;
   FileUtils::WriteScrCharLine(pof ,'/',false);
   if (!(bn->ImStp())) {bn->SetUpBNW();}
   if (!(bn->ballAndStickMode)) {bn->drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
#if DEBUG
   ScreenUtils::DisplayWarningMessage("In this version, calling CritPtNetWork::makePovFILE(...)\n\
         will overwrite the original coordinates of the critical points\n\
         and the coordinates on the bondnetwork object as well.");
#endif
   CenterMolecule();
   bn->CalcViewRadius();
   //cout << "rView: " << bn->rView << endl;
   double camdist=2.5e0;
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
         ScreenUtils::DisplayWarningMessage("Not a valid direction. Using 001.");
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
   double zsep=0.5e0;
   pvp.lightSource[1][0]=zsep;
   pvp.lightSource[1][1]=zsep;
   pvp.lightSource[1][2]=1.0e0;
   pvp.AddLightSource(zsep,-zsep,1.0e0);
   pvp.AddLightSource(-zsep,zsep,1.0e0);
   pvp.AddLightSource(-zsep,-zsep,1.0e0);
   for (int i=1; i<pvp.nLightSources; i++) {
      for (int j=0; j<3; j++) {pvp.lightSource[i][j]*=(bn->rView*4.0e0);}
   }
   for (int i=0; i<pvp.nLightSources; i++) {
      pvp.WriteLightSource(pof,i,0.5,"  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   }
   pof << "camera {" << endl;
   pof << "  up < 0, 1, 0 >" << endl;
   pof << "  right < -4/3, 0, 0 >" << endl;
   pof << "  location ";
   HelpersPOVRay::WriteVector(pof,pvp.locCam[0],pvp.locCam[1],pvp.locCam[2]);
   pof << endl << "  look_at ";
   HelpersPOVRay::WriteVector(pof,pvp.lookAtCam[0],pvp.lookAtCam[1],pvp.lookAtCam[2]);
   pof << endl << "  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >";
   pof << endl << "}" << endl;
   if (bn->spaceFillingMode) {
      pof << HelpersPOVRay::IndTabsStr(pvp.currIndLev) << "merge{" << endl;
      pvp.currIndLev++;
   }
   FileUtils::WriteScrCharLine(pof,'/',false);
   pof << "#if(DrawAtomTranspSpheres)" << endl;
   PutNuclei(pof);
   pof << "#end\n//end if DrawAtomTranspSpheres" << endl;
   FileUtils::WriteScrCharLine(pof,'/',false);
   pof << "#if(DrawStandardBonds)" << endl;
   PutBonds(pof);
   pof << "#end\n//end if DrawStandardBonds" << endl;
   FileUtils::WriteScrCharLine(pof,'/',false);
   if (iknowacps) {
      //int indacp;
      pof << "#if(DrawAttractorCriticalPoints)" << endl;
      //double acprad;
      for (int i=0; i<nACP; i++) {
         //acprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
         //HelpersPOVRay::WriteSphere(pof,0,RACP[i][0],RACP[i][1],RACP[i][2],acprad,
               //               0.0e0,0.0e0,0.0e0);
         HelpersPOVRay::WriteSphere(pof,0,RACP[i][0],RACP[i][1],RACP[i][2],"RadiusACP",
               "ColorACP");
      }
      pof << "#end\n//end if DrawAttractorCriticalPoints" << endl;
      FileUtils::WriteScrCharLine(pof,'/',false);
   }
   if (iknowbcps) {
      //int indbcp;
      pof << "#if(DrawBondCriticalPoints)" << endl;
      //double bcprad;
      for (int i=0; i<nBCP; i++) {
         //bcprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
         HelpersPOVRay::WriteSphere(pof,0,RBCP[i][0],RBCP[i][1],RBCP[i][2],"RadiusBCP",
               "ColorBCP");
      }
      pof << "#end\n//end if DrawBondCriticalPoints" << endl;
      FileUtils::WriteScrCharLine(pof,'/',false);
   }
   if (iknowrcps) {
      //int indrcp;
      pof << "#if(DrawRingCriticalPoints)" << endl;
      //double rcprad;
      for (int i=0; i<nRCP; i++) {
         //rcprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
         HelpersPOVRay::WriteSphere(pof,0,RRCP[i][0],RRCP[i][1],RRCP[i][2],"RadiusRCP",
               "ColorRCP");
      }
      pof << "#end\n//end if DrawRingCriticalPoints" << endl;
      FileUtils::WriteScrCharLine(pof,'/',false);
   }
   if (iknowccps) {
      //int indccp;
      pof << "#if(DrawCageCriticalPoints)" << endl;
      //double ccprad;
      for (int i=0; i<nCCP; i++) {
         //ccprad=bn->drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
         HelpersPOVRay::WriteSphere(pof,0,RCCP[i][0],RCCP[i][1],RCCP[i][2],"RadiusCCP",
               "ColorCCP");
      }
      pof << "#end\n//end if DrawCageCriticalPoints" << endl;
      FileUtils::WriteScrCharLine(pof,'/',false);
   }
   if (iknowbgps) {
      double gprad=0.06;
      int npts;
      pof << "#if(DrawGradientPathSpheres)" << endl;
      pof << "union {" << endl;
      for (int i=0; i<nBCP; i++) {
         npts=conBCP[i][2];
         HelpersPOVRay::WriteSphere(pof,1,RBGP[i][0][0],RBGP[i][0][1],RBGP[i][0][2], \
               gprad,"ColorABGradPath");
         for (int j=1; j<npts; j++) {
            HelpersPOVRay::WriteSphere(pof,1,RBGP[i][j][0],RBGP[i][j][1],RBGP[i][j][2], \
                  gprad,"ColorABGradPath");
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathSpheres" << endl;
   }
   if (iknowbgps) {
      double gprad=0.06;
      int npts;
      pof << "#if(DrawGradientPathTubes)" << endl;
      pof << "union {" << endl;
      for (int i=0; i<nBCP; i++) {
         npts=conBCP[i][2];
         HelpersPOVRay::WriteSphere(pof,1,RBGP[i][0][0],RBGP[i][0][1],RBGP[i][0][2], \
               gprad,"ColorABGradPath");
         for (int j=1; j<npts; j++) {
            HelpersPOVRay::WriteSphere(pof,1,RBGP[i][j][0],RBGP[i][j][1],RBGP[i][j][2], \
                  gprad,"ColorABGradPath");
            //if (tubeBGPStyle) {
            HelpersPOVRay::WriteCylinder(pof,1,\
                  RBGP[i][j][0],RBGP[i][j][1],RBGP[i][j][2], \
                  RBGP[i][j-1][0],RBGP[i][j-1][1],RBGP[i][j-1][2], \
                  gprad,"ColorABGradPath");
            //}
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathTubes" << endl;
   }
   if ( iknowrgps ) {
      double gprad=0.045;
      int npts,currBcpPos;
      pof << "#if(DrawGradientPathTubes)" << endl;
      pof << "union {" << endl;
      for (int rcpIdx=0; rcpIdx<nRCP; ++rcpIdx) {
         currBcpPos=0;
         while ( conRCP[rcpIdx][1][currBcpPos]>0 ) {
            npts=conRCP[rcpIdx][1][currBcpPos];
            HelpersPOVRay::WriteSphere(pof,1,RRGP[rcpIdx][currBcpPos][0][0],\
                  RRGP[rcpIdx][currBcpPos][0][1],\
                  RRGP[rcpIdx][currBcpPos][0][2],gprad,"ColorARGradPath");
            for ( int j=1 ; j<npts ; ++j ) {
               HelpersPOVRay::WriteSphere(pof,1,RRGP[rcpIdx][currBcpPos][j][0],\
                     RRGP[rcpIdx][currBcpPos][j][1],\
                     RRGP[rcpIdx][currBcpPos][j][2],gprad,"ColorARGradPath");
               if ( j%2 ==0 ) {
               HelpersPOVRay::WriteCylinder(pof,1,\
                     RRGP[rcpIdx][currBcpPos][j][0],\
                     RRGP[rcpIdx][currBcpPos][j][1],\
                     RRGP[rcpIdx][currBcpPos][j][2],\
                     RRGP[rcpIdx][currBcpPos][j-1][0],\
                     RRGP[rcpIdx][currBcpPos][j-1][1],\
                     RRGP[rcpIdx][currBcpPos][j-1][2],\
                     gprad,"ColorARGradPath");
               }
            }
            ++currBcpPos;
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathTubes" << endl;
      //-------------------
      pof << "#if(DrawGradientPathSpheres)" << endl;
      pof << "union {" << endl;
      for (int rcpIdx=0; rcpIdx<nRCP; ++rcpIdx) {
         currBcpPos=0;
         while ( conRCP[rcpIdx][1][currBcpPos]>0 ) {
            npts=conRCP[rcpIdx][1][currBcpPos];
            for ( int j=0 ; j<npts ; ++j ) {
               HelpersPOVRay::WriteSphere(pof,1,RRGP[rcpIdx][currBcpPos][j][0],\
                     RRGP[rcpIdx][currBcpPos][j][1],\
                     RRGP[rcpIdx][currBcpPos][j][2],gprad,"ColorARGradPath");
            }
            ++currBcpPos;
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathSpheres" << endl;
 
   }
   if ( iknowcgps ) {
      double gprad=0.045;
      int npts,currRcpPos;
      pof << "#if(DrawGradientPathTubes)" << endl;
      pof << "union {" << endl;
      for (int ccpIdx=0; ccpIdx<nCCP; ++ccpIdx) {
         currRcpPos=0;
         while ( conCCP[ccpIdx][0][currRcpPos]>=0 ) {
            npts=conCCP[ccpIdx][1][currRcpPos];
            if ( npts<2 ) { ++currRcpPos; continue; }
            HelpersPOVRay::WriteSphere(pof,1,RCGP[ccpIdx][currRcpPos][0][0],\
                  RCGP[ccpIdx][currRcpPos][0][1],\
                  RCGP[ccpIdx][currRcpPos][0][2],gprad,"ColorACGradPath");
            for ( int j=1 ; j<npts ; ++j ) {
               if (j%3 != 1) {
               HelpersPOVRay::WriteSphere(pof,1,RCGP[ccpIdx][currRcpPos][j][0],\
                     RCGP[ccpIdx][currRcpPos][j][1],\
                     RCGP[ccpIdx][currRcpPos][j][2],gprad,"ColorACGradPath");
               }
               if (j%3 ==0) {
               HelpersPOVRay::WriteCylinder(pof,1,\
                     RCGP[ccpIdx][currRcpPos][j][0],\
                     RCGP[ccpIdx][currRcpPos][j][1],\
                     RCGP[ccpIdx][currRcpPos][j][2],\
                     RCGP[ccpIdx][currRcpPos][j-1][0],\
                     RCGP[ccpIdx][currRcpPos][j-1][1],\
                     RCGP[ccpIdx][currRcpPos][j-1][2],\
                     gprad,"ColorACGradPath");
               }
            }
            ++currRcpPos;
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathTubes" << endl;
      //-------------------
      pof << "#if(DrawGradientPathSpheres)" << endl;
      pof << "union {" << endl;
      for (int ccpIdx=0; ccpIdx<nCCP; ++ccpIdx) {
         currRcpPos=0;
         while ( conCCP[ccpIdx][0][currRcpPos]>=0 ) {
            npts=conCCP[ccpIdx][1][currRcpPos];
            for ( int j=0 ; j<npts ; ++j ) {
               HelpersPOVRay::WriteSphere(pof,1,RCGP[ccpIdx][currRcpPos][j][0],\
                     RCGP[ccpIdx][currRcpPos][j][1],\
                     RCGP[ccpIdx][currRcpPos][j][2],gprad,"ColorACGradPath");
            }
            ++currRcpPos;
         }
      }
      pof << "}" << endl;
      pof << "#end\n//end if DrawGradientPathSpheres" << endl;
 
   }
   pof.close();
   return true;
}
void CritPtNetWork::PutBonds(ofstream &pof) {
   string pigmstr="transmit TransmitStdBondCylinder";
   pof << "union{" << endl;
   int k=0,atni,atnk;
   double startpt[3],frak1;
   for (int i=0; i<bn->nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bn->bNet[i][j];
         atni=bn->atNum[i];
         atnk=bn->atNum[k];
         //frak1=atomicRadius[atni]/(atomicRadius[atni]+atomicRadius[atnk]);
         frak1=GetAtomicVDWRadius(atni)/(GetAtomicVDWRadius(atni)+GetAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=bn->R[i][l]*(1.0e0-frak1)+bn->R[k][l]*frak1;
         }
         if (k>0) {
            HelpersPOVRay::WriteCylinder(pof,1,
                  bn->R[i][0],bn->R[i][1],bn->R[i][2],
                  startpt[0],startpt[1],startpt[2],
                  bn->drawStickSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR,
                  GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
                  GetAtomicBColorReal(atni),pigmstr);
            HelpersPOVRay::WriteCylinder(pof,1,
                  startpt[0],startpt[1],startpt[2],
                  bn->R[k][0],bn->R[k][1],bn->R[k][2],
                  bn->drawStickSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR,
                  GetAtomicRColorReal(atnk),GetAtomicGColorReal(atnk),
                  GetAtomicBColorReal(atnk),pigmstr);
         }
      }
      HelpersPOVRay::WriteSphere(pof,0,bn->R[i][0],bn->R[i][1],bn->R[i][2],
            bn->drawStickSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR,
            GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
            GetAtomicBColorReal(atni));
   }
   pof << "}" << endl;
   return;
}
void CritPtNetWork::PutNuclei(ofstream & pof) {
   int atomn;
   double atrad;
   string transmStr="TransmitAtomSphere";
   for (int i=0; i<bn->nNuc; i++) {
      atomn=bn->atNum[i];
      atrad=bn->drawAtSize;
      HelpersPOVRay::WriteTransparentSphere(pof,0,bn->R[i][0],bn->R[i][1],bn->R[i][2],atrad,
            GetAtomicRColorReal(atomn),GetAtomicGColorReal(atomn),
            GetAtomicBColorReal(atomn),transmStr);
   }
   return;
}
void CritPtNetWork::CenterMolecule(void) {
   double trn[3];
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
         for (int j=0; j<conBCP[i][2]; j++) {
            for (int k=0; k<3; k++) {
               RBGP[i][j][k]-=trn[k];
            }
         }
      }
   }
   if ( iknowrgps ) {
      int kn,kp;
      for ( int rcpIdx=0 ; rcpIdx<nRCP ; ++rcpIdx ) {
         kn=GetNofRingPathsOfRCP(rcpIdx);
         for ( int j=0 ; j<kn ; ++j ) {
            kp=conRCP[rcpIdx][1][j];
            for ( int k=0 ; k<kp ; ++k ) {
               for ( int l=0 ; l<3 ; ++l ) {RRGP[rcpIdx][j][k][l]-=trn[l];}
            }
         }
      }
   }
   if ( iknowcgps ) {
      int kn,kp;
      for ( int ccpIdx=0 ; ccpIdx<nCCP ; ++ccpIdx ) {
         kn=GetNofCagePathsOfCCP(ccpIdx);
         for ( int j=0 ; j<kn ; ++j ) {
            kp=conCCP[ccpIdx][1][j];
            for ( int k=0 ; k<kp ; ++k ) {
               for ( int l=0 ; l<3 ; ++l ) {RCGP[ccpIdx][j][k][l]-=trn[l];}
            }
         }
      }
   }
   bn->CenterMolecule();
   for (int i=0; i<3; i++) {centMolecVec[i]=trn[i];}
   return;
}
bool CritPtNetWork::ReadFromFile(string inname) {
   string tmps;
   tmps=inname.substr(inname.length()-3,inname.length());
   if (tmps!="cpx") {
      ScreenUtils::DisplayErrorMessage("Not a valid file!");
      return false;
   }
   ifstream cfil;
   cfil.open(inname.c_str(),std::ios::in);
   if (!(cfil.good())) {
      ScreenUtils::DisplayErrorMessage("This file could not be opened...");
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
         ScreenUtils::DisplayWarningMessage("This field has not been implemented!");
         exit(1);
         break;
   }
   if (dACP<0) {
      ScreenUtils::DisplayWarningMessage(string("ACPs are not given in file "+inname));
      iknowacps=false;
   } else {
      MyMemory::Alloc2DRealArray(string("RACP"),dACP,3,RACP,1.0e+50);
      MyMemory::Alloc1DStringArray("lblACP",dACP,lblACP);
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
      MyMemory::Alloc2DRealArray(string("RBCP"),dBCP,3,RBCP,1.0e+50);
      MyMemory::Alloc2DIntArray(string("conBCP"),dBCP,3,conBCP,-1);
      MyMemory::Alloc1DStringArray("lblBCP",dBCP,lblBCP);
      if (nBCP>=0) {
         cpxGetBCPCartCoordFromFile(cfil,nBCP,RBCP);
         cpxGetBCPConnectivityFromFile(cfil,nBCP,conBCP);
         cpxGetBCPLabelsFromFile(cfil,nBCP,lblBCP);
         iknowbcps=true;
      } else {
         iknowbcps=false;
      }
   } else {
      ScreenUtils::DisplayErrorMessage("First look for ACPs...\nQuiting...");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif
   }
   if (iknowbcps) {
      nRCP=cpxGetNOfRCPs(cfil);
      dRCP=(nBCP*(nBCP-1))/2;
      if ( dRCP<CPNW_MINARRAYSIZE ) {dRCP=CPNW_MINARRAYSIZE;}
      MyMemory::Alloc2DRealArray(string("RRCP"),dRCP,3,RRCP,1.0e+50);
      MyMemory::Alloc1DStringArray("lblRCP",dRCP,lblRCP);
      MyMemory::Alloc3DIntArray("conRCP",dRCP,2,CPNW_MAXBCPSCONNECTEDTORCP,conRCP,-1);
      if (nRCP>=0) {
         cpxGetRCPCartCoordFromFile(cfil,nRCP,RRCP);
         cpxGetRCPLabelsFromFile(cfil,nRCP,lblRCP);
         cpxGetRCPConnectivityFromFile(cfil,nRCP,conRCP);
         iknowrcps=true;
      } else {
         iknowrcps=false;
      }
   }
   if (iknowrcps) {
      nCCP=cpxGetNOfCCPs(cfil);
      dCCP=(nRCP*(nRCP-1))/2;
      if ( dCCP<CPNW_MINARRAYSIZE ) {dCCP=CPNW_MINARRAYSIZE;}
      MyMemory::Alloc2DRealArray(string("RCCP"),dCCP,3,RCCP,1.0e+50);
      MyMemory::Alloc1DStringArray("lblCCP",dCCP,lblCCP);
      MyMemory::Alloc3DIntArray("conCCP",dCCP,2,CPNW_MAXRCPSCONNECTEDTOCCP,conCCP,-1);
      if (nCCP>=0) {
         cpxGetCCPCartCoordFromFile(cfil,nCCP,RCCP);
         cpxGetCCPLabelsFromFile(cfil,nCCP,lblCCP);
         cpxGetCCPConnectivityFromFile(cfil,nCCP,conCCP);
         iknowccps=true;
      } else {
         iknowccps=false;
      }
   }
   iknowallcps=(iknowacps&&iknowbcps&&iknowrcps&&iknowccps);
   nBGP=cpxGetNOfBondPaths(cfil);
   if (dBCP>0) {MyMemory::Alloc3DRealArray(string("RBGP"),dBCP,maxGradPathNPts,3,RBGP);}
   if (nBGP>=0&&conBCP!=NULL) {
      cpxGetNOfPtsPerBondPath(cfil,nBGP,conBCP);
      cpxGetBondPathData(cfil,nBGP,conBCP,RBGP);
      iknowbgps=true;
   } else {iknowbgps=false;}
   nRGP=cpxGetNOfRingPaths(cfil);
   if (dRCP>0) {
      MyMemory::Alloc4DRealArray(string("RRGP"),dRCP,CPNW_MAXBCPSCONNECTEDTORCP,\
            maxGradPathNPts,3,RRGP);
   }
   if ( nRGP>=0 && conRCP!=NULL ) {
      cpxGetNOfPtsPerRingPath(cfil,nRCP,conRCP);
      cpxGetRingPathData(cfil,nRCP,conRCP,RRGP);
      iknowrgps=true;
   } else { iknowrgps=false; }
   nCGP=cpxGetNOfCagePaths(cfil);
   if (dCCP>0) {
      MyMemory::Alloc4DRealArray(string("RCGP"),dCCP,CPNW_MAXRCPSCONNECTEDTOCCP,\
            maxGradPathNPts,3,RCGP);
   }
   if ( nCGP>=0 && conCCP!=NULL ) {
      cpxGetNOfPtsPerCagePath(cfil,nCCP,conCCP);
      cpxGetCagePathData(cfil,nCCP,conCCP,RCGP);
      iknowcgps=true;
   } else { iknowcgps=false; }
   iknowallgps=iknowbgps&&iknowrgps&&iknowcgps;
   cout << "Critical Point State Loaded!" << endl;
   DisplayStatus(true);
   return true;
}
void CritPtNetWork::DisplayStatus(bool lngdesc) {
   ScreenUtils::PrintScrCharLine('+');
   if (lngdesc) {
      if (conBCP==NULL) {ScreenUtils::DisplayWarningMessage("conBCP is not allocated.");}
      if (RACP==NULL) {ScreenUtils::DisplayWarningMessage("RACP is not allocated.");}
      if (RBCP==NULL) {ScreenUtils::DisplayWarningMessage("RBCP is not allocated.");}
      if (RRCP==NULL) {ScreenUtils::DisplayWarningMessage("RRCP is not allocated.");}
      if (RCCP==NULL) {ScreenUtils::DisplayWarningMessage("RCCP is not allocated.");}
      if (RBGP==NULL) {ScreenUtils::DisplayWarningMessage("RBGP is not allocated.");}
      if (lblACP==NULL) {ScreenUtils::DisplayWarningMessage("lblACP is not allocated.");}
      if (lblBCP==NULL) {ScreenUtils::DisplayWarningMessage("lblBCP is not allocated.");}
      if (lblRCP==NULL) {ScreenUtils::DisplayWarningMessage("lblRCP is not allocated.");}
      if (lblCCP==NULL) {ScreenUtils::DisplayWarningMessage("lblCCP is not allocated.");}
   }
   ScreenUtils::PrintScrCharLine('-');
   if (iknowacps) {cout << nACP << " ACPs found." << endl;}
   else {ScreenUtils::DisplayWarningMessage("ACPs search needed...");}
   if (iknowbcps) {cout << nBCP << " BCPs found." << endl;}
   else {ScreenUtils::DisplayWarningMessage("BCPs search needed...");}
   if (iknowrcps) {cout << nRCP << " RCPs found." << endl;}
   else {ScreenUtils::DisplayWarningMessage("RCPs search needed...");}
   if (iknowccps) {cout << nCCP << " CCPs found." << endl;}
   else {ScreenUtils::DisplayWarningMessage("CCPs search needed...");}
   if (iknowallcps) {
      ScreenUtils::PrintScrCharLine('-');
      cout << "nACP-nBCP+nRCP-nCCP: " << (nACP-nBCP+nRCP-nCCP) << endl;
      ScreenUtils::PrintScrCharLine('-');
   }
   if (iknowbgps) {
      if (nBGP!=nBCP) {ScreenUtils::DisplayWarningMessage("For some unknown reason nBGP!=nBCP");}
      else {cout << nBGP << " Bond Gradient Paths found." << endl;}
   } else {ScreenUtils::DisplayWarningMessage("BGPs search needed...");}
   ScreenUtils::PrintScrCharLine('+');
   return;
}
int CritPtNetWork::FindSingleRhoBondGradientPathRK5(int at1,int at2,double hstep,\
      int dima,double** (&arbgp),double (&ro)[3]) {
   double rn[3],rho,g[3],h[3][3],eive[3][3],eival[3],dist,maggrad;
   SeekSingleRhoBCP(at1,at2,ro);
   wf->EvalHessian(ro[0],ro[1],ro[2],rho,g,h);
   EigenDecompositionJAMA::EigenDecomposition3(h,eive,eival);
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {maggrad+=(eive[i][2]*eive[i][2]);}
   maggrad=sqrt(maggrad);
   for (int i=0; i<3; i++) {rn[i]=ro[i]+hstep*eive[i][2]/maggrad;}
   wf->EvalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
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
   double loopdist;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)) {
      GetNextPointInGradientPathRK5UpHill(rn,hstep,maggrad);
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
      count++;
      if (count==dima) {
         ScreenUtils::DisplayWarningMessage("You need a bigger array for storing the BGP coordinates!");
         return count-1;
      }
   }
   InvertOrderBGPPoints(count,arbgp);
   maggrad=0.0e0;
   for (int i=0; i<3; i++) {maggrad+=(eive[i][2]*eive[i][2]);}
   maggrad=sqrt(maggrad);
   for (int i=0; i<3; i++) {rn[i]=ro[i]-hstep*eive[i][2]/maggrad;}
   wf->EvalRhoGradRho(rn[0],rn[1],rn[2],rho,g);
   maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   for (int i=0; i<3; i++) {
      arbgp[count][i]=rn[i];
   }
   count++;
   maxit+=count;
   iminacp=false;
   while ((!iminacp)&&(count<maxit)&&(maggrad>EPSGRADMAG)) {
      GetNextPointInGradientPathRK5UpHill(rn,hstep,maggrad);
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
         ScreenUtils::DisplayWarningMessage("You need a bigger array for storing the BGP coordinates!");
         return count-1;
      }
   }
   if (count==dima) {
      ScreenUtils::DisplayWarningMessage("You need a bigger array for storing the BGP coordinates!");
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
   for (int i=0; i<3; i++) {dist+=((arbgp[count-1][i]-wf->R[iacp2+i])\
         *(arbgp[count-1][i]-wf->R[iacp2+i]));}
   dist=sqrt(dist);
   if (dist<=hstep) {
      for (int i=0; i<3; i++) {arbgp[count][i]=wf->R[iacp2+i];}
      count++;
   }
   dist=0.0e0;
   for (int i=0; i<3; i++) {dist+=((arbgp[0][i]-wf->R[iacp1+i])\
         *(arbgp[0][i]-wf->R[iacp1+i]));}
   dist=sqrt(dist);
   if ( dist>(2.0e0*stepSizeBGP) ) {
      InvertOrderBGPPoints(count,arbgp);
   }
   return count;
}
int CritPtNetWork::FindSingleRhoRingGradientPathRK5(int rcpIdx,\
      int bcpIdxInRRGP,double hstep,int dima,\
      double** (&arrgp)) {
   int bcpGlobIdx=conRCP[rcpIdx][0][bcpIdxInRRGP];
#if DEBUG
   if ( bcpIdxInRRGP>CPNW_MAXBCPSCONNECTEDTORCP || bcpIdxInRRGP < 0 ) {
      ScreenUtils::DisplayErrorMessage("Out of conRCP bounds!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   if ( rcpIdx>nRCP ) {
      ScreenUtils::DisplayErrorMessage("rcpIdx>nRCP");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   if ( bcpGlobIdx>nBCP || bcpGlobIdx<0 ) {
      ScreenUtils::DisplayErrorMessage("Non existent bcp!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
#endif
   int count;
   double xm[3],xb[3],xr[3],xn[3],xrmxb[3]; //,xmmxr[3];
   double magd=0.0e0,maxalllen=bn->maxBondDist*1.5e0;
   for ( int i=0 ; i<3 ; ++i ) {
      xb[i]=RBCP[bcpGlobIdx][i];
      xr[i]=RRCP[rcpIdx][i];
      xrmxb[i]=xr[i]-xb[i];
      xn[i]=xrmxb[i];
      magd+=(xn[i]*xn[i]);
   }
   magd=sqrt(magd);
   for ( int i=0 ; i<3 ; ++i ) {
      xn[i]/=magd;
      xn[i]*=hstep;
      xn[i]+=xb[i];
   }
   bool imatrcp=WalkGradientPathRK5ToEndPoint(xb,xn,xr,xm,magd,hstep,\
         dima,arrgp,count,maxalllen,false /* uphilldir=false  */);
   if ( imatrcp ) {return count;}
   //cout << "rcpIdx: " << rcpIdx << ", bcpGlobIdx: " << bcpGlobIdx \
   //<< ", count: " << count << ", magd: " << magd << '\n';
   //if ( magd>maxBCPACPDist ) { return -1; }
   double ux[3],uy[3]; //,xm0mxr[3],uz[3];
   //for ( int i=0 ; i<3 ; ++i ) { xm0mxr[i]=xm[i]-xr[i]; }
   double hess[3][3],eivec[3][3],tmpv[3];
   wf->EvalHessian(xb[0],xb[1],xb[2],magd,ux,hess);
   if ( magV3(ux)>CPNW_EPSRHOACPGRADMAG ) {
      ScreenUtils::DisplayWarningMessage(string("Not at a BCP? (")+StringTools::GetStringFromReal(magV3(ux))\
            +string(")"));
   }
   EigenDecompositionJAMA::EigenDecomposition3(hess,eivec,tmpv);
   for ( int i=0 ; i<3 ; ++i ) {
      ux[i]=eivec[i][0];
      uy[i]=eivec[i][1];
      //uz[i]=eivec[i][2];
      //xmmxr[i]=xm[i]-xr[i];
   }
   double dir1[3],dir2[3];
   for ( int i=0 ; i<3 ; ++i ) {
      dir1[i]=xrmxb[i];
      dir2[i]=xm[i]-xb[i];
   }
   double alpha,beta,delta,gamma,currdmin;
   alpha=atan2(dotProductV3(dir1,uy),dotProductV3(dir1,ux));
   beta=atan2(dotProductV3(dir2,uy),dotProductV3(dir2,ux));
   delta=alpha-beta;
   gamma=alpha;
   int sign0=(delta >= 0.0e0 ? 1 : -1),currSign;
   int mymaxiter=fabs(int(6.28318530717959e0/delta))+10;
   bool samedir=true;
   int iter=0;
   magd=1.0e0;
   do {
      gamma+=delta;
      for ( int i=0 ; i<3 ; ++i ) {
         dir2[i]=(cos(gamma)*ux[i]+sin(gamma)*uy[i]);
         xn[i]=dir2[i];
      }
      normalizeV3(xn);
      for ( int i=0 ; i<3 ; ++i ) {
         xn[i]*=0.5e0*hstep;
         xn[i]+=xb[i];
      }
      imatrcp=WalkGradientPathRK5ToEndPoint(xb,xn,xr,xm,currdmin,hstep,\
            dima,arrgp,count,maxalllen,false);
      if ( imatrcp ) {
         //displayGreenMessage("imatrcp!");
         for ( int i=0 ; i<3 ; ++i ) { tmpv[i]=xm[i]-xb[i]; }
         magd=atan2(dotProductV3(tmpv,uy),dotProductV3(tmpv,ux));
         //cout << "iter: " << iter << ", delta: " << magd << endl;
         return count;
      }
      for ( int i=0 ; i<3 ; ++i ) { tmpv[i]=xm[i]-xb[i]; }
      magd=atan2(dotProductV3(tmpv,uy),dotProductV3(tmpv,ux));
      currSign=((alpha-magd)>=0? 1 : -1);
      if ( currSign!=sign0 ) {
         samedir=false;
         break;
      } else {
         for ( int i=0 ; i<3 ; ++i ) { dir1[i]=dir2[i]; }
      }
      ++iter;
   } while (samedir && iter<mymaxiter);
   iter=0;
   do {
      for ( int i=0 ; i<3 ; ++i ) { xn[i]=0.5e0*(dir1[i]+dir2[i]); }
      normalizeV3(xn);
      for ( int i=0 ; i<3 ; ++i ) {
         xn[i]*=hstep;
         xn[i]+=xb[i];
      }
      imatrcp=WalkGradientPathRK5ToEndPoint(xb,xn,xr,xm,currdmin,hstep,\
            dima,arrgp,count,maxalllen,false);
      if ( imatrcp ) { return count; }
      for ( int i=0 ; i<3 ; ++i ) { tmpv[i]=xm[i]-xb[i]; }
      magd=atan2(dotProductV3(tmpv,uy),dotProductV3(tmpv,ux));
      currSign=((alpha-magd)>=0? 1 : -1);
      if ( currSign == sign0 ) {
         for ( int i=0; i<3; ++i ) { dir1[i]+=dir2[i]; dir1[i]*=0.5e0; }
      } else {
         for ( int i=0; i<3; ++i ) { dir2[i]+=dir1[i]; dir2[i]*=0.5e0; }
      }
      ++iter;
   } while (iter<100 && fabs(magd)>CPNW_EPSRHOACPGRADMAG);
   //cout << "iter: " << iter << ", delta: " << magd << endl;
   for ( int i=0 ; i<3 ; ++i ) {tmpv[i]=xn[i]-xr[i];}
   if ( magV3(tmpv)<(5.0e0*hstep) ) {
      return count;
   } else {
      //cout << "iter: " << iter << ", |tmpv|: " << magV3(tmpv) << endl;
      return -1;
   }
}
int CritPtNetWork::FindSingleRhoCageGradientPathRK5(int ccpIdx,\
      int rcpIdxInRCGP,double hstep,int dima,\
      double** (&arrgp)) {
   int rcpGlobIdx=conCCP[ccpIdx][0][rcpIdxInRCGP];
#if DEBUG
   if ( rcpIdxInRCGP>CPNW_MAXRCPSCONNECTEDTOCCP || rcpIdxInRCGP < 0 ) {
      ScreenUtils::DisplayErrorMessage("Out of conCCP bounds!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   if ( ccpIdx>nCCP ) {
      ScreenUtils::DisplayErrorMessage("ccpIdx>nCCP");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   if ( rcpGlobIdx>nRCP || rcpGlobIdx<0 ) {
      ScreenUtils::DisplayErrorMessage("Non existent rcp!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
#endif
   double xn[3],xr[3],xc[3],xm[3];
   double magd;
   for ( int i=0 ; i<3 ; ++i ) {
      xr[i]=RRCP[rcpGlobIdx][i];
      xc[i]=RCCP[ccpIdx][i];
   }
   double hess[3][3],eivec[3][3],tmpv[3],dir2min[3];
   wf->EvalHessian(xr[0],xr[1],xr[2],magd,tmpv,hess); //here tmpv is gradRho
   if ( magV3(tmpv)>CPNW_EPSRHOACPGRADMAG ) {
      ScreenUtils::DisplayWarningMessage(string("Not at an RCP? (")+\
            StringTools::GetStringFromReal(magV3(tmpv))+string(")"));
   }
   EigenDecompositionJAMA::EigenDecomposition3(hess,eivec,tmpv);
   for ( int i=0 ; i<3 ; ++i ) {
      dir2min[i]=eivec[i][0];
      xn[i]=xr[i]+hstep*dir2min[i];
   }
   int count;
   double maxalllen=bn->maxBondDist;
   bool imatccp=WalkGradientPathRK5ToEndPoint(xr,xn,xc,xm,magd,hstep,\
         dima,arrgp,count,maxalllen,false); //uphill=false
   if ( imatccp ) {return count;}
   for ( int i=0 ; i<3 ; ++i ) {
      dir2min[i]=eivec[i][0];
      xn[i]=xr[i]-hstep*dir2min[i];
   }
   imatccp=WalkGradientPathRK5ToEndPoint(xr,xn,xc,xm,magd,hstep,\
         dima,arrgp,count,maxalllen,false); //uphill=false
   if ( imatccp ) {return count;}
   if ( magd<0.5e0 ) {
      ScreenUtils::DisplayErrorMessage("Unknown error!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      cout << "count: " << count << ", magd: " << magd << ", maxalllen: " << maxalllen << '\n';
      wf->DisplayAllFieldProperties(xn[0],xn[1],xn[2]);
   }
#if DEBUG
   cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   wf->DisplayAllFieldProperties(xn[0],xn[1],xn[2]);
#endif /* ( DEBUG ) */
   return -1;
}
bool CritPtNetWork::WalkGradientPathRK5ToEndPoint(\
      double (&xi)[3],double (&x1)[3],\
      double (&xe)[3],double (&xm)[3],double &dm,double hstep,int dima,\
      double** (&arrgp), int &npia,double maxlen,bool uphilldir) {
   double xn[3],xmin[3],dmin=0.0e0,magd=0.0e0;
   for ( int i=0 ; i<3 ; ++i ) {
      RGP[0][i]=xi[i];
      RGP[1][i]=x1[i];
      xn[i]=x1[i];
      xmin[i]=xn[i]-xe[i];
      dmin+=(xmin[i]*xmin[i]);
      magd+=((xn[i]-x1[i])*(xn[i]-x1[i]));
   }
   double mxlen2=maxlen*maxlen,pathlength=magd;
   double epsd2=hstep*hstep;
   //double epsd2=CPNW_EPSFABSDIFFCOORD*CPNW_EPSFABSDIFFCOORD;
   double maggrad=1.0e0;
   bool imatend=false;
   int count=2;
   while ((!imatend)&&(count<dima)&&(pathlength<mxlen2)) {
      if ( uphilldir ) {
         GetNextPointInGradientPathRK5UpHill(xn,hstep,maggrad);
      } else {
         GetNextPointInGradientPathRK5DownHill(xn,hstep,maggrad);
      }
      magd=maggrad=0.0e0;
      for (int i=0; i<3; ++i) {
         magd+=((xn[i]-xe[i])*(xn[i]-xe[i]));
         RGP[count][i]=xn[i];
         maggrad+=((xn[i]-RGP[count-1][i])*(xn[i]-RGP[count-1][i]));
      }
      pathlength+=maggrad;
      if ( magd<=epsd2 || maggrad<1.0e-6 ) {
         imatend=true;
      }
      if ( magd<=dmin ) {
         dmin=magd;
         for ( int i=0 ; i<3 ; ++i ) { xmin[i]=xn[i]; }
      } else {
         ++count;
         break;
      }
      ++count;
   }
   if ( count==maxGradPathNPts ) {
      ScreenUtils::DisplayWarningMessage("End or array reached, need larger array?");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
   }
   for ( int i=0 ; i<3 ; ++i ) {xm[i]=xmin[i];}
   dm=sqrt(dmin);
   //if ( imatend ) {
      CopyRGP2Array(arrgp,count);
      npia=count;
   //} else {
      //CopyRGP2Array(arrgp,count);
      //npia=count;
   //}
   return imatend;
}
bool CritPtNetWork::SeekSingleRhoBCP(int ata,int atb,double (&x)[3]) {
   int idxa,idxb;
   idxa=3*ata;
   idxb=3*atb;
   //cout << "ata: " << ata << ", atb: " << atb << endl;
   for (int j=0; j<3; j++) {x[j]=0.5e0*(wf->R[idxa+j]+wf->R[idxb+j]);}
   double rho,g[3];
   int sig;
   SeekRhoBCP(x,rho,g,sig);
   //wf->EvalRhoGradRho(x[0],x[1],x[2],rho,g);
   //double magg=0.0e0;
   //for (int n=0; n<3; n++) {magg+=(g[n]*g[n]);}
   //magg=sqrt(magg);
   if (!((rho>CPNW_MINRHOSIGNIFICATIVEVAL)&&(ComputeMagnitudeV3(g)<CPNW_EPSRHOACPGRADMAG))) {
      ScreenUtils::DisplayErrorMessage(string("The chosen atoms ("+wf->atLbl[ata]\
               +","+wf->atLbl[atb]+") probably do not lead to a BCP!"));
#if DEBUG
      wf->DisplayAllFieldProperties(x[0],x[1],x[2]);
#endif
      return false;
   }
   return true;
}
void CritPtNetWork::GetNextPointInGradientPathRK5UpHill(double (&xn)[3],double &stepsize,\
      double &mgg) {
   /*
   static const double b[15]={0.2e0, 3.0e0/40.0e0, 9.0e0/40.0e0, 0.3e0, -0.9e0, 1.2e0, \
      -11.0e0/54.0e0, 2.5e0, -70.0e0/27.0e0, 35.0e0/27.0e0, 1631.0e0/55296.0e0, \
         175.0e0/512.0e0, 575.0e0/13824.0e0, 44275.0e0/110592.0e0, 253.0e0/4096.0e0};
   static const double c1=37.0e0/378.0e0;
   static const double c3=250.0e0/621.0e0;
   static const double c4=125.0e0/594.0e0;
   static const double c6=512.0e0/1771.0e0;
   double k[6][3],rho,g[3],xt[3],maggrad;
   int offset;
   for(int i=0; i<6; i++) {
      offset=((i*(i-1))>>1);
      for(int j=0; j<3; j++) {
         xt[j] = xn[j];
         for(int l=0; l<i; l++) {
            xt[j] += b[l+offset]*k[l][j];
         }
      }
      wf->EvalRhoGradRho(xt[0],xt[1],xt[2],rho,g);
      maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      for(int l=0; l<3; l++) {
         k[i][l] = stepsize*g[l]/maggrad;
      }
   }
   // */
   static const double b[21]={1.0e0/5.0e0, \
      3.0e0/40.0e0, 9.0e0/40.0e0, \
      44.0e0/45.0e0, -56.0e0/15.0e0, 32.0e0/9.0e0, \
      19372.0e0/6561.0e0, -25360.0e0/2187.0e0, 64448.0e0/6561.0e0, -212.0e0/729.0e0,\
      9017.0e0/3168.0e0, -355.0e0/33.0e0, 46732.0e0/5247.0e0, 49.0e0/176.0e0, -5103.0e0/18656.0e0,\
      35.0e0/384.0e0, 0.0e0, 500.0e0/1113.0e0, 125.0e0/192.0e0, -2187.0e0/6784.0e0,11.0e0/84.0e0};
   static const double c1=35.0e0/384.0e0;
   static const double c3=500.0e0/1113.0e0;
   static const double c4=125.0e0/192.0e0;
   static const double c5=-2187.0e0/6784.0e0;
   static const double c6=11.0e0/84.0e0;
   double k[7][3],rho,g[3],xt[3],maggrad;
   int offset;
   for(int i=0; i<7; i++) {
      offset=((i*(i-1))>>1);
      for(int j=0; j<3; j++) {
         xt[j] = xn[j];
         for(int l=0; l<i; l++) {
            xt[j] += b[l+offset]*k[l][j];
         }
      }
      wf->EvalRhoGradRho(xt[0],xt[1],xt[2],rho,g);
      maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      for(int l=0; l<3; l++) {
         k[i][l] = stepsize*g[l]/maggrad;
      }
   }
   for(int i=0; i<3; i++) {
      xn[i]+=(c1*k[0][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
   }
   wf->EvalRhoGradRho(xn[0],xn[1],xn[2],rho,g);
   mgg=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   return;
}
void CritPtNetWork::GetNextPointInGradientPathRK5DownHill(double (&xn)[3],\
      double &stepsize,double &mgg) {
   static const double b[21]={1.0e0/5.0e0, \
      3.0e0/40.0e0, 9.0e0/40.0e0, \
         44.0e0/45.0e0, -56.0e0/15.0e0, 32.0e0/9.0e0, \
         19372.0e0/6561.0e0, -25360.0e0/2187.0e0, \
         64448.0e0/6561.0e0, -212.0e0/729.0e0,\
         9017.0e0/3168.0e0, -355.0e0/33.0e0, 46732.0e0/5247.0e0, \
         49.0e0/176.0e0, -5103.0e0/18656.0e0,\
         35.0e0/384.0e0, 0.0e0, 500.0e0/1113.0e0, 125.0e0/192.0e0, \
         -2187.0e0/6784.0e0,11.0e0/84.0e0};
   static const double c1=35.0e0/384.0e0;
   static const double c3=500.0e0/1113.0e0;
   static const double c4=125.0e0/192.0e0;
   static const double c5=-2187.0e0/6784.0e0;
   static const double c6=11.0e0/84.0e0;
   double k[7][3],rho,g[3],xt[3],maggrad;
   int offset;
   for(int i=0; i<7; i++) {
      offset=((i*(i-1))>>1);
      for(int j=0; j<3; j++) {
         xt[j] = xn[j];
         for(int l=0; l<i; l++) {
            xt[j] += b[l+offset]*k[l][j];
         }
      }
      wf->EvalRhoGradRho(xt[0],xt[1],xt[2],rho,g);
      maggrad=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      for(int l=0; l<3; l++) {
         k[i][l] = -stepsize*g[l]/maggrad;
      }
   }
   for(int i=0; i<3; i++) {
      xn[i]+=(c1*k[0][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
   }
   wf->EvalRhoGradRho(xn[0],xn[1],xn[2],rho,g);
   mgg=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   return;
}
void CritPtNetWork::FindMaxBondDist() {
   if ( !iknowbgps ) {
      ScreenUtils::DisplayErrorMessage("First you need to find the Bond Gradient Paths!");
      return;
   }
   if ( wf->nNuc==1 ) { maxBondDist=0.0e0; return; }
   double magd;
   int bcp1,bcp2;
   for ( int i=0 ; i<nBCP ; ++i ) {
      bcp1=conBCP[i][0];
      bcp2=conBCP[i][1];
      magd=0.0e0;
      for ( int j=0 ; j<3 ; ++j ) {
         magd+=((RACP[bcp1][j]-RACP[bcp2][j])*(RACP[bcp1][j]-RACP[bcp2][j]));
      }
      if ( magd>maxBondDist ) { maxBondDist=magd; }
      magd=0.0e0;
      for ( int j=0 ; j<3 ; ++j ) {
         magd+=((RBCP[i][j]-RACP[bcp1][j])*(RBCP[i][j]-RACP[bcp1][j]));
      }
      if ( magd>maxBCPACPDist ) { maxBCPACPDist=magd; }
      magd=0.0e0;
      for ( int j=0 ; j<3 ; ++j ) {
         magd+=((RBCP[i][j]-RACP[bcp2][j])*(RBCP[i][j]-RACP[bcp2][j]));
      }
      if ( magd>maxBCPACPDist ) { maxBCPACPDist=magd; }
   }
   maxBondDist=sqrt(maxBondDist);
   maxBCPACPDist=sqrt(maxBCPACPDist);
}
void CritPtNetWork::AddBCP2ConRCP(const int rcpIdx,const int bcpIdx) {
#if DEBUG
   if ( conRCP==NULL ) {
      ScreenUtils::DisplayErrorMessage("conRCP is not allocated!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
#endif /* ( DEBUG ) */
   int k=0;
   while ( k<CPNW_MAXBCPSCONNECTEDTORCP && conRCP[rcpIdx][0][k]!=bcpIdx ) {
      if ( conRCP[rcpIdx][0][k]<0 ) {conRCP[rcpIdx][0][k]=bcpIdx; return;}
      ++k;
   }
   if ( k==CPNW_MAXBCPSCONNECTEDTORCP ) {
      ScreenUtils::DisplayWarningMessage("Perhaps you need a larger array for conRCP!");
#if DEBUG
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
   }
}
void CritPtNetWork::AddRCP2ConCCP(const int ccpIdx,const int rcpIdx) {
#if DEBUG
   if ( conCCP==NULL ) {
      ScreenUtils::DisplayErrorMessage("conCCP is not allocated!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
#endif /* ( DEBUG ) */
   int k=0;
   while ( k<CPNW_MAXRCPSCONNECTEDTOCCP && conCCP[ccpIdx][0][k]!=rcpIdx ) {
      if ( conCCP[ccpIdx][0][k]<0 ) {conCCP[ccpIdx][0][k]=rcpIdx; return;}
      ++k;
   }
   if ( k==CPNW_MAXRCPSCONNECTEDTOCCP ) {
      ScreenUtils::DisplayWarningMessage("Perhaps you need a larger array for conCCP!");
#if DEBUG
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
   }
}
void CritPtNetWork::ForceBCPConnectivity(int bcpIdx,int acpIdx1,int acpIdx2) {
#if DEBUG
   if ( conBCP==NULL ) {
      ScreenUtils::DisplayErrorMessage("conBCP is not allocated!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
#endif
   if ( bcpIdx>=nBCP ) {
      ScreenUtils::DisplayWarningMessage(string("The BCP of index ")+StringTools::GetStringFromInt(bcpIdx)+string(" does not exist!"));
      ScreenUtils::DisplayWarningMessage("Nothingo to do...");
      return;
   }
   if ( acpIdx1>=nACP ) {
      ScreenUtils::DisplayWarningMessage(string("The ACP of index ")+StringTools::GetStringFromInt(acpIdx1)+string(" does not exist!"));
      ScreenUtils::DisplayWarningMessage("Nothingo to do...");
      return;
   }
   if ( acpIdx2>=nACP ) {
      ScreenUtils::DisplayWarningMessage(string("The ACP of index ")+StringTools::GetStringFromInt(acpIdx2)+string(" does not exist!"));
      ScreenUtils::DisplayWarningMessage("Nothingo to do...");
      return;
   }
   conBCP[bcpIdx][0]=acpIdx1;
   conBCP[bcpIdx][1]=acpIdx2;
   lblBCP[bcpIdx]=lblACP[acpIdx1]+"-"+lblACP[acpIdx2];
}
void CritPtNetWork::CorrectRCPConnectivity(void) {
   FindMaxBondDist();
   //double maxallwdd=1.414213562373095e0*maxBCPACPDist;
   double maxallwdd=0.9*maxBondDist; //This must be CritPtNetWork::maxBondDist.
   double dd;
   int j,bcpIdx;
   for ( int i=0 ; i<nRCP ; ++i ) {
      j=0;
      while ( j<CPNW_MAXBCPSCONNECTEDTORCP && conRCP[i][0][j]>=0 ) {
         bcpIdx=conRCP[i][0][j];
         dd=0.0e0;
         for ( int k=0; k<3 ; ++k ) {
            dd+=((RRCP[i][k]-RBCP[bcpIdx][k])*(RRCP[i][k]-RBCP[bcpIdx][k]));
         }
         dd=sqrt(dd);
         //cout << "d: " << dd << "; " << lblBCP[bcpIdx] << endl;
         if ( dd<=maxallwdd ) {
            ++j;
         } else {
            RemoveFromConRCP(i,j);
         }
      }
   }
   double tmpv[3];
   for ( int i=0 ; i<nRCP ; ++i ) {
      for ( int j=0 ; j<nBCP ; ++j ) {
         for ( int k=0 ; k<3 ; ++k ) {tmpv[k]=RRCP[i][k]-RBCP[j][k];}
         if ( magV3(tmpv)<maxBCPACPDist ) {AddToConRCP(i,j);}
      }
   }
}
void CritPtNetWork::RemoveFromConRCP(const int rcpIdx,const int pos2rem) {
   int total=0;
   while ( conRCP[rcpIdx][0][total]>=0 ) {++total;}
   --total;
   if ( pos2rem==total ) {
      conRCP[rcpIdx][0][total]=conRCP[rcpIdx][1][total]=-1;
      return;
   }
   int tmp=conRCP[rcpIdx][0][total];
   conRCP[rcpIdx][0][pos2rem]=tmp;
   tmp=conRCP[rcpIdx][1][total];
   conRCP[rcpIdx][1][pos2rem]=tmp;
   conRCP[rcpIdx][0][total]=conRCP[rcpIdx][1][total]=-1;
   /*
   int total=0;
   while ( conRCP[rcpIdx][0][total]>=0 ) {++total;}
   if ( conRCP[rcpIdx][0][total-1]==bcpIdx ) {
      conRCP[rcpIdx][0][total-1]=conRCP[rcpIdx][1][total-1]=-1;
      return;
   }
   int mypos=0;
   while ( conRCP[rcpIdx][0][mypos]!=bcpIdx && mypos<total ) {++mypos;}
   --total;
   int tmp=conRCP[rcpIdx][0][total];
   conRCP[rcpIdx][0][mypos]=tmp;
   tmp=conRCP[rcpIdx][1][total];
   conRCP[rcpIdx][1][mypos];
   conRCP[rcpIdx][0][total]=conRCP[rcpIdx][1][total]=-1;
   // */
}
void CritPtNetWork::AddToConRCP(const int rcpIdx,const int bcpIdx) {
   int mypos=0;
   while (conRCP[rcpIdx][0][mypos]>=0) {
      if ( conRCP[rcpIdx][0][mypos]==bcpIdx ) {return;}
      ++mypos;
      if ( mypos==CPNW_MAXBCPSCONNECTEDTORCP ) {
         ScreenUtils::DisplayWarningMessage("End of array reached, need bigger array!");
#if DEBUG
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
         return;
      }
   }
   conRCP[rcpIdx][0][mypos]=bcpIdx;
}
void CritPtNetWork::AddToConCCP(const int ccpIdx,const int rcpIdx) {
   int mypos=0;
   while (conCCP[ccpIdx][0][mypos]>=0) {
      if ( conCCP[ccpIdx][0][mypos]==rcpIdx ) {return;}
      ++mypos;
      if ( mypos==CPNW_MAXRCPSCONNECTEDTOCCP ) {
         ScreenUtils::DisplayWarningMessage("End of array reached, need bigger array!");
#if DEBUG
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
         return;
      }
   }
   conCCP[ccpIdx][0][mypos]=rcpIdx;
}
void CritPtNetWork::SetRingPaths() {
   if (!iknowrcps) {
      ScreenUtils::DisplayErrorMessage("Please look first for the RCPs...\nNothing to be done!");
      return;
   }
   MyMemory::Alloc4DRealArray(string("RRGP"),dRCP,CPNW_MAXBCPSCONNECTEDTORCP,\
         maxGradPathNPts,3,RRGP);
   CorrectRCPConnectivity();
   cout << "Calculating Ring Gradient Paths..." << endl;
#if DEBUG
   cout << "nRCP: " << nRCP << endl;
#endif /* ( DEBUG ) */
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double hstep; //,rseed[3];
   hstep=stepSizeBGP;
   int currBcpPos,npts;
#if DEBUG
   int bcpIdx;
#endif
   for (int rcpIdx=0; rcpIdx<nRCP; rcpIdx++) {
      currBcpPos=0;
      while ( conRCP[rcpIdx][0][currBcpPos]>=0 ) {
#if DEBUG
         bcpIdx=conRCP[rcpIdx][0][currBcpPos];
         if ( bcpIdx>=nBCP || bcpIdx<0 ) {
            ScreenUtils::DisplayErrorMessage("BCP out of bounds!");
            cout << "bcpIdx: " << bcpIdx << " (nBCP=" << nBCP << ')' << '\n';
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         }
#endif
         //for (int k=0; k<3; k++) {rseed[k]=RBCP[bcpIdx][k];}
         npts=FindSingleRhoRingGradientPathRK5(rcpIdx,\
               currBcpPos,hstep,maxGradPathNPts,RRGP[rcpIdx][currBcpPos]);
#if DEBUG
         if ( npts<0 ) {
            ScreenUtils::DisplayWarningMessage("Catched -1");
            cout << "bcpIdx: " << bcpIdx << " (nBCP=" << nBCP << ')' << '\n';
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         }
#endif /* ( DEBUG ) */
         conRCP[rcpIdx][1][currBcpPos]=npts;
         ++currBcpPos;
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(rcpIdx)/\
               double( (nRCP>1) ? (nRCP-1) : 1 )));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   iknowrgps=true;
   iknowallgps=(iknowbgps&&iknowrgps&&iknowcgps);
}
void CritPtNetWork::SetCagePaths(void) {
  if (!iknowccps) {
      ScreenUtils::DisplayErrorMessage("Please look first for the CCPs...\nNothing to be done!");
      return;
   }
   MyMemory::Alloc4DRealArray(string("RCGP"),dCCP,CPNW_MAXRCPSCONNECTEDTOCCP,\
         maxGradPathNPts,3,RCGP);
   cout << "Calculating Cage Gradient Paths..." << endl;
#if DEBUG
   cout << "nCCP: " << nCCP << endl;
#endif /* ( DEBUG ) */
   /* Remember this is just a trial function,  */
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   double hstep; //,rseed[3];
   hstep=stepSizeBGP;
   //int rcpIdx;
   int currRcpPos,npts;
   for (int ccpIdx=0; ccpIdx<nCCP; ccpIdx++) {
      currRcpPos=0;
      while ( conCCP[ccpIdx][0][currRcpPos]>=0 ) {
         //rcpIdx=conCCP[ccpIdx][0][currRcpPos];
#if DEBUG
         if ( ccpIdx>=nCCP || ccpIdx<0 ) {
            ScreenUtils::DisplayErrorMessage("CCP out of bounds!");
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         }
#endif
         //for (int k=0; k<3; k++) {rseed[k]=RCCP[ccpIdx][k];}
         npts=FindSingleRhoCageGradientPathRK5(ccpIdx,\
               currRcpPos,hstep,maxGradPathNPts,RCGP[ccpIdx][currRcpPos]);
#if DEBUG
         if ( npts<0 ) {
            ScreenUtils::DisplayWarningMessage("Catched -1");
            cout << "ccpIdx: " << ccpIdx << " (nCCP=" << nCCP << ')' << '\n';
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         }
#endif /* ( DEBUG ) */
         conCCP[ccpIdx][1][currRcpPos]=npts;
         ++currRcpPos;
      }
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(ccpIdx)/\
               double( (nCCP>1) ? (nCCP-1) : 1 )));
#endif
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   iknowcgps=true;
   iknowallgps=(iknowbgps&&iknowrgps&&iknowcgps);
}
int CritPtNetWork::GetNofRingPathsOfRCP(int rcpIdx) {
   if ( !iknowrgps ) {
      ScreenUtils::DisplayWarningMessage("First seek for Ring Gradient Paths! Returning 0.");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
   }
   int res=0;
   while ( conRCP[rcpIdx][0][res]>=0 ) {++res;}
   return res;
}
int CritPtNetWork::GetNofCagePathsOfCCP(int ccpIdx) {
   if ( !iknowcgps ) {
      ScreenUtils::DisplayWarningMessage("First seek for Cage Gradient Paths! Returning 0.");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
   }
   int res=0;
   while ( conCCP[ccpIdx][0][res]>=0 ) {++res;}
   return res;
}
int CritPtNetWork::GetTotalNofRingPaths(void) {
#if DEBUG
   if ( !iknowrgps ) {
      ScreenUtils::DisplayWarningMessage("First seek for Ring Gradient Paths! Returning -1.");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return -1;
   }
#endif /* ( DEBUG ) */
   int res=0;
   for ( int i=0 ; i<nRCP ; ++i ) {res+=GetNofRingPathsOfRCP(i);}
   return res;
}
int CritPtNetWork::GetTotalNofCagePaths(void) {
#if DEBUG
   if ( !iknowcgps ) {
      ScreenUtils::DisplayWarningMessage("First seek for Cage Gradient Paths! Returning -1.");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return -1;
   }
#endif /* ( DEBUG ) */
   int res=0;
   for ( int i=0 ; i<nCCP ; ++i ) {res+=GetNofCagePathsOfCCP(i);}
   return res;
}
void CritPtNetWork::CopyRGP2Array(double** (&thearr),int nn) {
   for ( int i=0 ; i<nn ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) {thearr[i][j]=RGP[i][j];}
   }
}

#endif  /* _CRITPTNETWORK_CPP_ */

