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
#include "solmath.h"
#include "eig2-4.h"
#include "eig6.h"

#ifndef DM1CPNWBP_DEFAULTSTEPSIZEBGP
#define DM1CPNWBP_DEFAULTSTEPSIZEBGP 0.025e0
#endif

#ifndef DM1CPNWBP_DEFAULTARRAYSIZEBGP
#define DM1CPNWBP_DEFAULTARRAYSIZEBGP 300
#endif

#ifndef DM1CPNWBP_DEFAULTEPSDISTANCE
#define DM1CPNWBP_DEFAULTEPSDISTANCE 1.0e-04
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
void DeMat1CriticalPointNetworkBP::ComputeCoreInteractionCPs2D(void) {
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
   for ( int i=0 ; i<nCICP ; ++i ) { ComputeSingleCICP2D(i); }
   //ComputeSingleCICP2D(0);
}
void DeMat1CriticalPointNetworkBP::ComputeCoreInteractionCPs6D(void) {
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
   for ( int i=0 ; i<nCICP ; ++i ) { ComputeSingleCICP6D(i); }
   //ComputeSingleCICP6D(0);
}
bool DeMat1CriticalPointNetworkBP::CPSafetyChecks(void) {
   return cpn->iKnowBGPs();
}
void DeMat1CriticalPointNetworkBP::ComputeSingleCICP2D(int idx) {
#if DEBUG
   if ( idx>=nCICP ) {
      displayErrorMessage("Trying to compute non existent CICP!");
      cout << "nCICP: " << nCICP << endl;
      DISPLAYDEBUGINFOFILELINE;
      return;
   }
#endif
   int acp1Idx=cpn->conBCP[idx][0],acp2Idx=cpn->conBCP[idx][1];
   solreal xx[3],xxp[3],gamm,gg[3],gp[3],hh[3][3],hph[3][3],hp[3][3];
   for ( int i=0 ; i<3 ; ++i ) {
      xx[i]=cpn->RACP[acp1Idx][i];
      xxp[i]=cpn->RACP[acp2Idx][i];
   }
   wf->evalHessDensityMatrix1(xx,xxp,gamm,gg,gp,hh,hph,hp);
   solreal e1[3],e2[3];
   GetTangentialVectors(idx,e1,e2);
   solreal huv[2][2],tmp;
   huv[0][0]=huv[1][1]=huv[0][1]=0.0e0;
   for ( int i=0 ; i<3 ; ++i ) {
      tmp=0.0e0;
      for ( int j=0 ; j<3 ; ++j ) { tmp+=(hh[i][j]*e1[j]); }
      huv[0][0]+=tmp*e1[i];
      tmp=0.0e0;
      for ( int j=0 ; j<3 ; ++j ) { tmp+=hph[i][j]*e1[j]; }
      huv[0][1]+=tmp*e2[i];
      tmp=0.0e0;
      for ( int j=0 ; j<3 ; ++j ) { tmp+=hp[i][j]*e2[j]; }
      huv[1][1]+=tmp*e2[i];
   }
   huv[1][0]=huv[0][1];
   solreal eivec[2][2],eival[2];
   eigen_decomposition2(huv,eivec,eival);
   cout << cpn->lblBCP[idx] << " " << eival[0] << " " << eival[1] << " (signature: "
        << GetSignature(eival) << ")" << endl;
}
void DeMat1CriticalPointNetworkBP::ComputeSingleCICP6D(int idx) {
#if DEBUG
   if ( idx>=nCICP ) {
      displayErrorMessage("Trying to compute non existent CICP!");
      cout << "nCICP: " << nCICP << endl;
      DISPLAYDEBUGINFOFILELINE;
      return;
   }
#endif
   int acp1Idx=cpn->conBCP[idx][0],acp2Idx=cpn->conBCP[idx][1];
   solreal xx[3],xxp[3],gamm,gg[3],gp[3],hh[3][3],hph[3][3],hp[3][3];
   for ( int i=0 ; i<3 ; ++i ) {
      xx[i]=cpn->RACP[acp1Idx][i];
      xxp[i]=cpn->RACP[acp2Idx][i];
   }
   wf->evalHessDensityMatrix1(xx,xxp,gamm,gg,gp,hh,hph,hp);
   solreal hess[6][6],eivec[6][6],eival[6];
   /* ************************************************************************** */
   hess[0][0]=hh[0][0]; hess[0][1]=hh[0][1]; hess[0][2]=hh[0][2]; hess[0][3]=hph[0][0]; hess[0][4]=hph[0][1]; hess[0][5]=hph[0][2];
   hess[1][0]=hh[1][0]; hess[1][1]=hh[1][1]; hess[1][2]=hh[1][2]; hess[1][3]=hph[1][0]; hess[1][4]=hph[1][1]; hess[1][5]=hph[1][2];
   hess[2][0]=hh[2][0]; hess[2][1]=hh[2][1]; hess[2][2]=hh[2][2]; hess[2][3]=hph[2][0]; hess[2][4]=hph[2][1]; hess[2][5]=hph[2][2];
   /* ************************************************************************** */
   hess[3][0]=hph[0][0]; hess[3][1]=hph[1][0]; hess[3][2]=hph[2][0]; hess[3][3]=hp[0][0]; hess[3][4]=hp[0][1]; hess[3][5]=hp[0][2];
   hess[4][0]=hph[0][1]; hess[4][1]=hph[1][1]; hess[4][2]=hph[2][1]; hess[4][3]=hp[1][0]; hess[4][4]=hp[1][1]; hess[4][5]=hp[1][2];
   hess[5][0]=hph[0][2]; hess[5][1]=hph[1][2]; hess[5][2]=hph[2][2]; hess[5][3]=hp[2][0]; hess[5][4]=hp[2][1]; hess[5][5]=hp[2][2];
   /* ************************************************************************** */
   eigen_decomposition6(hess,eivec,eival);
   cout << cpn->lblBCP[idx];
   for ( int i=0 ; i<6 ; ++i ) { cout << " " << eival[i]; }
   cout << " (signature: "
        << GetSignature(eival) << ")" << endl;
}
void DeMat1CriticalPointNetworkBP::GetTangentialVectors(const int bcpIdx,solreal (&e1)[3],\
      solreal (&e2)[3]) {
   int npbgp=cpn->conBCP[bcpIdx][2];
   solreal xi[3],xip1[3];
   for ( int i=0 ; i<3 ; ++i ) {
      xi[i]=cpn->RBGP[bcpIdx][0][i];
      xip1[i]=cpn->RBGP[bcpIdx][1][i];
   }
   for ( int i=0 ; i<3 ; ++i ) { e1[i]=xip1[i]-xi[i]; }
   for ( int i=0 ; i<3 ; ++i ) {
      xi[i]=cpn->RBGP[bcpIdx][npbgp-2][i];
      xip1[i]=cpn->RBGP[bcpIdx][npbgp-1][i];
   }
   for ( int i=0 ; i<3 ; ++i ) { e2[i]=xip1[i]-xi[i]; }
   normalizeV3(e1);
   normalizeV3(e2);
}
int DeMat1CriticalPointNetworkBP::GetSignature(solreal (&v)[2]) {
   int s=0;
   for ( int i=0 ; i<2 ; ++i ) { v[i] >=0.0e0 ? ++s : --s; }
   return s;
}
int DeMat1CriticalPointNetworkBP::GetSignature(solreal (&v)[6]) {
   int s=0;
   for ( int i=0 ; i<6 ; ++i ) { v[i] >=0.0e0 ? ++s : --s; }
   return s;
}

#endif  /* _DEMAT1CRITPTNETWORKBP_CPP_ */


