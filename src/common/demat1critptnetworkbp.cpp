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

/* ************************************************************************** */
DeMat1CriticalPointNetworkBP::DeMat1CriticalPointNetworkBP(\
      GaussWaveFunction &usrwf,critPtNetWork &usrcp) {
   init();
   if ( !usrwf.imldd ) { 
      cout << "Error: First load a Gaussian Wave Function!" << endl;
      return;
   }
   wf=&usrwf;
   cp=&usrcp;
   SafetyChecks();
   imsetup=true;
}
/* ************************************************************************** */
DeMat1CriticalPointNetworkBP::~DeMat1CriticalPointNetworkBP() {
   destroy();
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::init() {
   wf=NULL;
   cp=NULL;
   imsetup=false;
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::destroy(void) {
   wf=NULL;
   cp=NULL;
   imsetup=false;
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::SafetyChecks(void) {
   imsetup=false;
   if ( wf->nNuc<2 ) {
      displayErrorMessage("The wavefunction must have at least two atoms!");
   }
   if ( !cp->iKnowACPs() ) {
      displayErrorMessage("Please find the ACPs before setting up "
            "DeMat1CriticalPointNetworkBP object!");
   }
   imsetup=true;
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::SetupBondPath(int ata,int atb) {
   at1=ata; at2=atb;
}
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_CPP_ */

