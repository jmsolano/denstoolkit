/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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
#include "bondnetwork.h"
#include "screenutils.h"
#include "mymemory.h"
#include "mymath.h"
#include "eigendecompositionjama.h"
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
      ScreenUtils::DisplayErrorMessage("DeMat1CriticalPointNetworkBP object could not be setut!");
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
   eivalNN6D=eivalCICP2D=eivalCICP6D=NULL;
   sigNN6D=sigCICP2D=sigCICP6D=NULL;
   imsetup=false;
   nCICP=0;
}
void DeMat1CriticalPointNetworkBP::destroy(void) {
   wf=NULL;
   bn=NULL;
   if ( cpn!=NULL ) { delete cpn; cpn=NULL; }
   MyMemory::Dealloc2DRealArray(eivalCICP2D,nCICP);
   MyMemory::Dealloc2DRealArray(eivalCICP6D,nCICP);
   MyMemory::Dealloc2DRealArray(eivalNN6D,nCICP);
   MyMemory::Dealloc1DIntArray(sigCICP2D);
   MyMemory::Dealloc1DIntArray(sigCICP6D);
   MyMemory::Dealloc1DIntArray(sigNN6D);
   imsetup=false;
}
bool DeMat1CriticalPointNetworkBP::InitSafetyChecks(void) {
   if ( wf->nNuc<2 ) {
      ScreenUtils::DisplayErrorMessage("The wavefunction must have at least two atoms!");
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
   AllocAuxArrays();
   return true;
}
bool DeMat1CriticalPointNetworkBP::AllocAuxArrays(void) {
   bool res;
   if ( !(res=cpn->iKnowBGPs()) ) {
      ScreenUtils::DisplayErrorMessage("First setup the critPtNetWork object!");
      return res;
   }
   nCICP=cpn->nBGP;
   res=res&&MyMemory::Alloc2DRealArray(string("eivalCICP2D"),nCICP,2,eivalCICP2D);
   res=res&&MyMemory::Alloc2DRealArray(string("eivalCICP6D"),nCICP,6,eivalCICP6D);
   res=res&&MyMemory::Alloc2DRealArray(string("eivalNN6D"),nCICP,6,eivalNN6D);
   res=res&&MyMemory::Alloc1DIntArray(string("sigCICP2D"),nCICP,sigCICP2D);
   res=res&&MyMemory::Alloc1DIntArray(string("sigCICP6D"),nCICP,sigCICP6D);
   res=res&&MyMemory::Alloc1DIntArray(string("sigNN6D"),nCICP,sigNN6D);
   return res;
}
void DeMat1CriticalPointNetworkBP::ComputeCoreInteractionCPs2D(void) {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("DeMat1CriticalPointNetworkBP object is not setup!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
      return;
   }
   if ( !CPSafetyChecks() ) {
      ScreenUtils::DisplayErrorMessage("DeMat1CriticalPointNetworkBP's internal critPtNetWork is not ready!");
      return;
   }
   for ( int i=0 ; i<nCICP ; ++i ) { ComputeSingleCICP2D(i); }
   //ComputeSingleCICP2D(0);
}
void DeMat1CriticalPointNetworkBP::ComputeCoreInteractionCPs6D(void) {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("DeMat1CriticalPointNetworkBP object is not setup!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
      return;
   }
   if ( !CPSafetyChecks() ) {
      ScreenUtils::DisplayErrorMessage("DeMat1CriticalPointNetworkBP's internal critPtNetWork is not ready!");
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
      ScreenUtils::DisplayErrorMessage("Trying to compute non existent CICP!");
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
   EigenDecompositionJAMA::EigenDecomposition2(huv,eivec,eival);
   sigCICP2D[idx]=GetSignature(eival);
   for ( int i=0 ; i<2 ; ++i ) { eivalCICP2D[idx][i]=eival[i]; }
#if DEBUG
   for ( int i=0 ; i<2 ; ++i ) { cout << " " << eival[i]; }
   cout << " (signature cicp: "
        << sigCICP2D[idx] << ")" << endl;
#endif /* ( DEBUG ) */
}
void DeMat1CriticalPointNetworkBP::ComputeSingleCICP6D(int idx) {
#if DEBUG
   if ( idx>=nCICP ) {
      ScreenUtils::DisplayErrorMessage("Trying to compute non existent CICP!");
      cout << "nCICP: " << nCICP << endl;
      DISPLAYDEBUGINFOFILELINE;
      return;
   }
#endif
   /* acp-acp (cicp) determination  */
   int acp1Idx=cpn->conBCP[idx][0],acp2Idx=cpn->conBCP[idx][1];
   solreal xx[3],xxp[3],gamm,gg[3],gp[3],hh[3][3],hph[3][3],hp[3][3];
   for ( int i=0 ; i<3 ; ++i ) {
      xx[i]=cpn->RACP[acp1Idx][i];
      xxp[i]=cpn->RACP[acp2Idx][i];
   }
   wf->evalHessDensityMatrix1(xx,xxp,gamm,gg,gp,hh,hph,hp);
   solreal hess[6][6],eivec[6][6],eival[6];
   assignHessian6D(hh,hph,hp,hess);
   eigen_decomposition6(hess,eivec,eival);
   sigCICP6D[idx]=GetSignature(eival);
   for ( int i=0 ; i<6 ; ++i ) { eivalCICP6D[idx][i]=eival[i]; }
   /* nuc-nuc correlation determination  */
   cpn->findTwoClosestAtomsToBCP(idx,acp1Idx,acp2Idx);
   for ( int i=0 ; i<3 ; ++i ) {
      xx[i]=bn->R[acp1Idx][i];
      xxp[i]=bn->R[acp2Idx][i];
   }
   wf->evalHessDensityMatrix1(xx,xxp,gamm,gg,gp,hh,hph,hp);
   assignHessian6D(hh,hph,hp,hess);
   eigen_decomposition6(hess,eivec,eival);
   sigNN6D[idx]=GetSignature(eival);
   for ( int i=0 ; i<6 ; ++i ) { eivalNN6D[idx][i]=eival[i]; }
#if DEBUG
   cout << cpn->lblBCP[idx] << " (eival cicp): ";
   for ( int i=0 ; i<6 ; ++i ) { cout << " " << eivalCICP6D[idx][i]; }
   cout << " (signature cicp: "
        << sigCICP6D[idx] << ")" << endl;
   cout << "(eival nuc-nuc [" << bn->atLbl[acp1Idx] << "," 
        << bn->atLbl[acp2Idx] << "] ): ";
   for ( int i=0 ; i<6 ; ++i ) { cout << " " << eivalNN6D[idx][i]; }
   cout << " (signature nuc-nuc: "
        << sigNN6D[idx] << ")" << endl;
#endif /* ( DEBUG ) */
}
void DeMat1CriticalPointNetworkBP::assignHessian6D(solreal (&hh)[3][3],\
      solreal (&hph)[3][3],solreal (&hp)[3][3],solreal (&hess)[6][6]) {
   /* ************************************************************************** */
   hess[0][0]=hh[0][0]; hess[0][1]=hh[0][1]; hess[0][2]=hh[0][2]; hess[0][3]=hph[0][0]; hess[0][4]=hph[0][1]; hess[0][5]=hph[0][2];
   hess[1][0]=hh[1][0]; hess[1][1]=hh[1][1]; hess[1][2]=hh[1][2]; hess[1][3]=hph[1][0]; hess[1][4]=hph[1][1]; hess[1][5]=hph[1][2];
   hess[2][0]=hh[2][0]; hess[2][1]=hh[2][1]; hess[2][2]=hh[2][2]; hess[2][3]=hph[2][0]; hess[2][4]=hph[2][1]; hess[2][5]=hph[2][2];
   /* ************************************************************************** */
   hess[3][0]=hph[0][0]; hess[3][1]=hph[1][0]; hess[3][2]=hph[2][0]; hess[3][3]=hp[0][0]; hess[3][4]=hp[0][1]; hess[3][5]=hp[0][2];
   hess[4][0]=hph[0][1]; hess[4][1]=hph[1][1]; hess[4][2]=hph[2][1]; hess[4][3]=hp[1][0]; hess[4][4]=hp[1][1]; hess[4][5]=hp[1][2];
   hess[5][0]=hph[0][2]; hess[5][1]=hph[1][2]; hess[5][2]=hph[2][2]; hess[5][3]=hp[2][0]; hess[5][4]=hp[2][1]; hess[5][5]=hp[2][2];
   /* ************************************************************************** */
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
bool DeMat1CriticalPointNetworkBP::differentSignaturesCICPvsNN(void) {
   bool res=false;
   for ( int i=0 ; i<nCICP ; ++i ) {
      if ( sigNN6D[i]!=sigCICP6D[i] ) {
         res=true;
         break;
      }
   }
   return res;
}
#endif  /* _DEMAT1CRITPTNETWORKBP_CPP_ */


