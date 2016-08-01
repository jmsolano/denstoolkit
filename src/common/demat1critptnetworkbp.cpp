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

/* ************************************************************************** */
DeMat1CriticalPointNetworkBP::DeMat1CriticalPointNetworkBP(\
      GaussWaveFunction &usrwf,critPtNetWork &usrcp) {
   init();
   if ( !usrwf.imldd ) { 
      cout << "Error: First load a Gaussian Wave Function!" << endl;
      return;
   }
   if ( !usrcp.iKnowBGPs() ) {
      cout << "Error: First seek for Bond Gradient Paths!" << endl;
      return;
   }
   wf=&usrwf;
   cp=&usrcp;
}
/* ************************************************************************** */
DeMat1CriticalPointNetworkBP::~DeMat1CriticalPointNetworkBP() {
   destroy();
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::init() {
   wf=NULL;
   cp=NULL;
}
/* ************************************************************************** */
void DeMat1CriticalPointNetworkBP::destroy(void) {
   wf=NULL;
   cp=NULL;
}
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_CPP_ */

