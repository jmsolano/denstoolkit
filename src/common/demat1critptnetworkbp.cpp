#ifndef _DEMAT1CRITPTNETWORKBP_CPP_
#define _DEMAT1CRITPTNETWORKBP_CPP_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "demat1critptnetworkbp.h"
#include "gausswavefunction.h"
#include "critptnetwork.h"
#include "solscrutils.h"
#include "solmemhand.h"

#ifndef DM1CPNWBP_DEFAULTSTEPSIZEBGP
#define DM1CPNWBP_DEFAULTSTEPSIZEBGP 0.025e0
#endif

#ifndef DM1CPNWBP_DEFAULTARRAYSIZEBGP
#define DM1CPNWBP_DEFAULTARRAYSIZEBGP 300
#endif

DeMat1CriticalPointNetworkBP::DeMat1CriticalPointNetworkBP(\
      GaussWaveFunction &usrwf,bondNetWork &usrbn) {
   init();
   if ( !usrwf.imldd ) { 
      cout << "Error: First load a Gaussian Wave Function!" << endl;
      return;
   }
   wf=&usrwf;
   bn=&usrbn;
   imsetup=InitSafetyChecks();
   if ( !imsetup ) {
      displayErrorMessage("DeMat1CriticalPointNetworkBP object could not be setut!");
      return;
   }
   cpn=new critPtNetWork(usrwf,usrbn);
   imsetup=imsetup&&SetupCPN();
}
DeMat1CriticalPointNetworkBP::~DeMat1CriticalPointNetworkBP() {
   destroy();
}
void DeMat1CriticalPointNetworkBP::init() {
   wf=NULL;
   cpn=NULL;
   gCICP=NULL;
   hCICP=NULL;
   imsetup=false;
   nCICP=0;
}
void DeMat1CriticalPointNetworkBP::destroy(void) {
   wf=NULL;
   bn=NULL;
   if ( cpn!=NULL ) { delete cpn; cpn=NULL; }
   dealloc1DRealArray(gCICP);
   dealloc2DRealArray(hCICP,nCICP);
   imsetup=false;
}
bool DeMat1CriticalPointNetworkBP::InitSafetyChecks(void) {
   if ( wf->nNuc<2 ) {
      displayErrorMessage("The wavefunction must have at least two atoms!");
      return false;
   }
   return true;
}
bool DeMat1CriticalPointNetworkBP::SetupCPN(void) {
   cpn->setMaxGradPathNPts(DM1CPNWBP_DEFAULTARRAYSIZEBGP);
   cpn->setStepSizeBGP(DM1CPNWBP_DEFAULTSTEPSIZEBGP);
   cpn->setCriticalPoints(DENS);
   cpn->setBondPaths();
   nCICP=cpn->nBGP;
   return true;
}
void DeMat1CriticalPointNetworkBP::AllocAuxArrays(void) {
   if ( !cpn->iKnowBGPs() ) {
      displayErrorMessage("First setup the critPtNetWork object!");
   }
   nCICP=cpn->nBGP;
   alloc1DRealArray(string("gCICP"),nCICP,gCICP);
   alloc2DRealArray(string("hCICP"),nCICP,3,hCICP);
}
void DeMat1CriticalPointNetworkBP::ComputeCoreInteractionCPs(void) {
   if ( !imsetup ) {
      displayErrorMessage("DeMat1CriticalPointNetworkBP object is not setup!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
      return;
   }
   if ( !CPSafetyChecks() ) {
      displayErrorMessage("DeMat1CriticalPointNetworkBP's internal critPtNetWork is not ready!");
      return;
   }
   //for ( int i=0 ; i<nCICP ; ++i ) { ComputeSingleCICP(i); }
   ComputeSingleCICP(0);
}
bool DeMat1CriticalPointNetworkBP::CPSafetyChecks(void) {
   return cpn->iKnowBGPs();
}
void DeMat1CriticalPointNetworkBP::ComputeSingleCICP(int idx) {
#if DEBUG
   if ( idx>=nCICP ) {
      displayErrorMessage("Trying to compute non existent CICP!");
      cout << "nCICP: " << nCICP << endl;
      DISPLAYDEBUGINFOFILELINE;
      return;
   }
#endif
   solreal tmp;
   int np=cpn->conBCP[idx][2];
   int j=cpn->conBCP[idx][0];
   for ( int i=0 ; i<3 ; ++i ) { cout << cpn->RACP[j][i] << " "; }
   cout << endl;
   for ( int i=0 ; i<np ; ++i ) {
      tmp=((cpn->RBGP[idx][i][0])*(cpn->RBGP[idx][i][0]))+\
          ((cpn->RBGP[idx][i][1])*(cpn->RBGP[idx][i][1]))+\
          ((cpn->RBGP[idx][i][2])*(cpn->RBGP[idx][i][2]));
      cout << cpn->RBGP[idx][i][0] << " " << cpn->RBGP[idx][i][1] << " "\
           << cpn->RBGP[idx][i][2] << " " << sqrt(tmp) << endl;
   }
   j=cpn->conBCP[idx][1];
   for ( int i=0 ; i<3 ; ++i ) { cout << cpn->RACP[j][i] << " "; }
   cout << endl;
}


#endif  /* _DEMAT1CRITPTNETWORKBP_CPP_ */

