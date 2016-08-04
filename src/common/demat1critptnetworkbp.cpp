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

#ifndef DM1CPNWBP_DEFAULTSTEPSIZEBGP
#define DM1CPNWBP_DEFAULTSTEPSIZEBGP 0.025e0
#endif

/* ************************************************************************** */
DeMat1CriticalPointNetworkBP::DeMat1CriticalPointNetworkBP(\
      GaussWaveFunction &usrwf,bondNetWork &usrbn) {
   init();
   if ( !usrwf.imldd ) { 
      cout << "Error: First load a Gaussian Wave Function!" << endl;
      return;
   }
   wf=&usrwf;
   bn=&usrbn;
   cpn=new critPtNetWork(usrwf,usrbn);
   imsetup=SafetyChecks();
}
/* ************************************************************************** */
DeMat1CriticalPointNetworkBP::~DeMat1CriticalPointNetworkBP() {
   destroy();
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::init() {
   wf=NULL;
   cpn=NULL;
   imsetup=false;
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::destroy(void) {
   wf=NULL;
   bn=NULL;
   if ( cpn!=NULL ) { delete cpn; cpn=NULL; }
   imsetup=false;
}
/* ************************************************************************** */
bool DeMat1CriticalPointNetworkBP::SafetyChecks(void) {
   if ( wf->nNuc<2 ) {
      displayErrorMessage("The wavefunction must have at least two atoms!");
      return false;
   }
   cpn->setCriticalPoints(DENS);
   cpn->setStepSizeBGP(DM1CPNWBP_DEFAULTSTEPSIZEBGP);
   cpn->setBondPaths();
   return true;
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::SeekBondPath(int ata,int atb) {
   at1=ata; at2=atb;
   solreal xbcp[3];
   cpn->seekSingleRhoBCP(at1,at2,xbcp);
   double hstep=DM1CPNWBP_DEFAULTSTEPSIZEBGP;
   int arrsize=CPNW_ARRAYSIZEGRADPATH;
   int npts=cpn->findSingleRhoBondGradientPathRK5(at1,at2,hstep,arrsize,cpn->RBGP[0],xbcp);
   cout << "npts: " << npts << endl;
}
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_CPP_ */

