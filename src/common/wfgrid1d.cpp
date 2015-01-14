//
//  wfgrid1d.cpp
//  
//
//  Created by Juan Manuel Solano on 2013-06-03.
//
//

#ifndef _WFGRID1D_CPP_
#define _WFGRID1D_CPP_

#include "wfgrid1d.h"
#include "solmemhand.h"
#include "solmath.h"

/* ********************************************************************************** */
/* ********************************************************************************** */
waveFunctionGrid1D::waveFunctionGrid1D()
{
   for (int i=0; i<3; i++) {
      Ca[i]=0.0e0;
      Cb[i]=0.0e0;
   }
   npts=DEFAULTPOINTSPERDIRECTION;
   dx=1.0e0/solreal(DEFAULTPOINTSPERDIRECTION);
   comments=string("#");
   prop1d=NULL;
   prop2plot=NONE;
   imsetup=false;
}
/* ********************************************************************************** */
waveFunctionGrid1D::~waveFunctionGrid1D()
{
   dealloc1DRealArray(prop1d);
}
/* ********************************************************************************** */
void waveFunctionGrid1D::setNPts(int nn)
{
   npts=nn;
   return;
}
/* ********************************************************************************** */
int waveFunctionGrid1D::getNPts(void)
{
   return npts;
}
/* ********************************************************************************** */
void waveFunctionGrid1D::setUpSimpleLine(bondNetWork &bn,int na,int nb)
{
   if (!(bn.imstp())) {
      cout << "Error: Trying to use a non set-up bondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   if ((na>=bn.nNuc)||(nb>=bn.nNuc)) {
      setScrRedBoldFont();
      cout << "Error: The atom you requested does not exist!" << endl;
      cout << "The grid could not be set up." << endl;
      setScrNormalFont();
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      exit(1);
      return;
   }
   solreal va[3],vb[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i];
      vb[i]=bn.R[nb][i];
   }
   setUpSimpleLine(bn,va,vb);
   dx=2.0e0/solreal(npts-1);
   alloc1DRealArray("prop1d",npts,prop1d);
   imsetup=true;
   return;
}
/* ********************************************************************************** */
void waveFunctionGrid1D::setUpSimpleLine(bondNetWork &bn,solreal (&ta)[3],solreal (&tb)[3])
{
   solreal eta[3],orig[3];
   for (int i=0; i<3; i++) {
      eta[i]=tb[i]-ta[i];
      orig[i]=0.5e0*(ta[i]+tb[i]);
   }
   normalizeV3(eta);
   solreal tmp;
   maxdim=tmp=0.0e0;
   for (int i=0; i<3; i++) {
      tmp+=(bn.bbmin[i]*bn.bbmin[i]);
      maxdim+=(bn.bbmax[i]*bn.bbmax[i]);
   }
   if (tmp>maxdim) {
      maxdim=sqrt(tmp);
   } else {
      maxdim=sqrt(maxdim);
   }
   maxdim+=(EXTRASPACELINEFACTOR*bn.maxBondDist);
   //cout << "maxdim: " << maxdim <<endl;
   for (int i=0; i<3; i++) {
      Ca[i]=(maxdim*(-eta[i]))+orig[i];
      Cb[i]=(maxdim*(eta[i]))+orig[i];
   }
   //cout << "Ca: " << Ca[0] << " " << Ca[1] << " " << Ca[2] << endl;
   //cout << "Cb: " << Cb[0] << " " << Cb[1] << " " << Cb[2] << endl;
   //cout << "Or: " << orig[0] << " " << orig[1] << " " << orig[2] << endl;
   return;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatRho(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalDensity(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
void waveFunctionGrid1D::setUpSimpleLine(bondNetWork &bn,int na)
{
   if (!(bn.imstp())) {
      cout << "Error: Trying to use a non set-up bondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   if (na>=bn.nNuc) {
      cout << "Error: The atom you requested does not exist!" << endl;
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   solreal va[3],vb[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i]-1.0e0;
      vb[i]=va[i]+2.0e0;
   }
   setUpSimpleLine(bn,va,vb);
   dx=2.0e0/solreal(npts-1);
   alloc1DRealArray("prop1d",npts,prop1d);
   imsetup=true;
   return;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatLapRho(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalLapRho(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatELF(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalELF(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatLOL(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalLOL(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatMagGradLOL(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3],gl[3],hl[3][3],lol;
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      wf.evalHessLOL(xx,lol,gl,hl);
      prop1d[i]=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      //prop1d[i]=wf.evalLOL(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatShannonEntropy(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalShannonEntropy(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatMagGradRho(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalMagGradRho(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatKinetEnerDensG(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalKineticEnergyG(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatKinetEnerDensK(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalKineticEnergyK(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatMolElecPot(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalMolElecPot(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatMagLED(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalMagLED(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatRedDensGrad(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalReducedDensityGradient(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatRoSE(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalRoSE(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
bool waveFunctionGrid1D::writeLineDatScalarCustFld(ofstream &ofil,gaussWaveFunc &wf)
{
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << solreal(i) << " " << 2.0e0*solreal(i) << " " << 3.2e0*solreal(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   solreal delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=solreal(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.evalCustomScalarField(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   solreal e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
/* ********************************************************************************** */
void waveFunctionGrid1D::makeDat(string &onam,gaussWaveFunc &wf,ScalarFieldType ft)
{
   if (!wf.imldd) {
      cout << "Error: trying to use a non loaded wave function object!\nNothing done!\n";
      return;
   }
   char cft=convertScalarFieldType2Char(ft);
   comments+=string("Property: ");
   comments+=getFieldTypeKeyLong(cft);
   ofstream ofil;
   ofil.open(onam.c_str());
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   switch (ft) {
      case DENS:
         writeLineDatRho(ofil,wf);
         break;
      case LAPD:
         writeLineDatLapRho(ofil,wf);
         break;
      case ELFD:
         writeLineDatELF(ofil,wf);
         break;
      case SENT:
         writeLineDatShannonEntropy(ofil,wf);
         break;
      case MGRD:
         writeLineDatMagGradRho(ofil,wf);
         break;
      case LOLD:
         writeLineDatLOL(ofil,wf);
         break;
      case MGLD:
         writeLineDatMagGradLOL(ofil,wf);
         break;
      case KEDG:
         writeLineDatKinetEnerDensG(ofil,wf);
         break;
      case KEDK:
         writeLineDatKinetEnerDensK(ofil,wf);
         break;
      case MEPD:
         writeLineDatMolElecPot(ofil,wf);
         break;
      case MLED :
         writeLineDatMagLED(ofil,wf);
         break;
      case ROSE :
         writeLineDatRoSE(ofil,wf);
         break;
      case REDG :
         writeLineDatRedDensGrad(ofil,wf);
         break;
      case SCFD :
         writeLineDatScalarCustFld(ofil,wf);
         break;
      default:
         cout << "Error: Field type not known!\n dat file incomplete!" << endl;
         break;
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   ofil.close();
   return;
}
/* ********************************************************************************** */

#endif//_WFGRID1D_CPP_
