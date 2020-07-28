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

/* wfgrid1d.cpp
   Created by Juan Manuel Solano on 2013-06-03.
*/
#ifndef _WFGRID1D_CPP_
#define _WFGRID1D_CPP_
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include "wfgrid1d.h"
#include "mymemory.h"
#include "mymath.h"

WaveFunctionGrid1D::WaveFunctionGrid1D() {
   for (int i=0; i<3; i++) {
      Ca[i]=0.0e0;
      Cb[i]=0.0e0;
   }
   npts=DEFAULTPOINTSPERDIRECTION;
   dx=1.0e0/double(DEFAULTPOINTSPERDIRECTION);
   comments=string("#");
   prop1d=NULL;
   prop2plot=NONE;
   imsetup=false;
}
WaveFunctionGrid1D::~WaveFunctionGrid1D() {
   MyMemory::Dealloc1DRealArray(prop1d);
}
void WaveFunctionGrid1D::SetNPts(int nn) {
   npts=nn;
   return;
}
int WaveFunctionGrid1D::GetNPts(void) {
   return npts;
}
void WaveFunctionGrid1D::SetUpSimpleLine(bondNetWork &bn,int na,int nb) {
   if (!(bn.imstp())) {
      cout << "Error: Trying to use a non set-up bondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   if ((na>=bn.nNuc)||(nb>=bn.nNuc)) {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: The atom you requested does not exist!" << endl;
      cout << "The grid could not be set up." << endl;
      ScreenUtils::SetScrNormalFont();
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      exit(1);
      return;
   }
   double va[3],vb[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i];
      vb[i]=bn.R[nb][i];
   }
   SetUpSimpleLine(bn,va,vb);
   dx=2.0e0/double(npts-1);
   MyMemory::Alloc1DRealArray("prop1d",npts,prop1d);
   imsetup=true;
   return;
}
void WaveFunctionGrid1D::SetUpSimpleLine(bondNetWork &bn,double (&ta)[3],double (&tb)[3]) {
   double eta[3],orig[3];
   for (int i=0; i<3; i++) {
      eta[i]=tb[i]-ta[i];
      orig[i]=0.5e0*(ta[i]+tb[i]);
   }
   normalizeV3(eta);
   double tmp;
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
bool WaveFunctionGrid1D::WriteLineDatRho(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalDensity(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
void WaveFunctionGrid1D::SetUpSimpleLine(bondNetWork &bn,int na) {
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
   double va[3],vb[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i]-1.0e0;
      vb[i]=va[i]+2.0e0;
   }
   SetUpSimpleLine(bn,va,vb);
   dx=2.0e0/double(npts-1);
   MyMemory::Alloc1DRealArray("prop1d",npts,prop1d);
   imsetup=true;
   return;
}
bool WaveFunctionGrid1D::WriteLineDatLapRho(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalLapRho(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatELF(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalELF(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatLOL(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalLOL(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatMagGradLOL(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3],gl[3],hl[3][3],lol;
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      wf.EvalHessLOL(xx,lol,gl,hl);
      prop1d[i]=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
      //prop1d[i]=wf.EvalLOL(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatShannonEntropy(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalShannonEntropy(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatMagGradRho(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalMagGradRho(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalKineticEnergyG(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalKineticEnergyK(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatMolElecPot(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalMolElecPot(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatMagLED(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalMagLED(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatRedDensGrad(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalReducedDensityGradient(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatRoSE(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalRoSE(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatScalarCustFld(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalCustomScalarField(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalVirialPotentialEnergyDensity(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
bool WaveFunctionGrid1D::WriteLineDatEllipticity(ofstream &ofil,GaussWaveFunction &wf) {
   if (!imsetup) {
      cout << "Error: the grid has not been set up!" << endl;
      cout << "No output will be written." << endl;
      for (int j=0; j<4; j++) {
         for (int i=0; i<4; i++) {
            ofil << double(i) << " " << 2.0e0*double(i) << " " << 3.2e0*double(i) << endl;
         }
         ofil << endl;
      }
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return false;
   }
   double delta[3],xx[3];
   for (int i=0; i<3; i++) {
      delta[i]=(Cb[i]-Ca[i]);
      xx[i]=Ca[i];
   }
   for (int i=0; i<3; i++) {delta[i]/=double(npts-1);}
   for (int i=0; i<npts; i++) {
      prop1d[i]=wf.EvalEllipticity(xx[0],xx[1],xx[2]);
      for (int j=0; j<3; j++) {xx[j]+=delta[j];}
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
   }
   //cout << "xx: " << xx[0] << " " << xx[1] << " " << xx[2] << endl;
   double e2=-1.0e0*maxdim;
   for (int i=0; i<npts; i++) {
      ofil << e2 << " " << prop1d[i] << endl;
      e2+=dx*maxdim;
   }
   return true;
}
void WaveFunctionGrid1D::MakeDat(string &onam,GaussWaveFunction &wf,ScalarFieldType ft) {
   if (!wf.imldd) {
      cout << "Error: trying to use a non loaded wave function object!\nNothing done!\n";
      return;
   }
   char cft=ConvertScalarFieldType2Char(ft);
   comments+=string("Property: ");
   comments+=GetFieldTypeKeyLong(cft);
   ofstream ofil;
   ofil.open(onam.c_str());
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch (ft) {
      case DENS:
         WriteLineDatRho(ofil,wf);
         break;
      case LAPD:
         WriteLineDatLapRho(ofil,wf);
         break;
      case ELFD:
         WriteLineDatELF(ofil,wf);
         break;
      case SENT:
         WriteLineDatShannonEntropy(ofil,wf);
         break;
      case MGRD:
         WriteLineDatMagGradRho(ofil,wf);
         break;
      case LOLD:
         WriteLineDatLOL(ofil,wf);
         break;
      case MGLD:
         WriteLineDatMagGradLOL(ofil,wf);
         break;
      case KEDG:
         WriteLineDatKinetEnerDensG(ofil,wf);
         break;
      case KEDK:
         WriteLineDatKinetEnerDensK(ofil,wf);
         break;
      case MEPD:
         WriteLineDatMolElecPot(ofil,wf);
         break;
      case MLED :
         WriteLineDatMagLED(ofil,wf);
         break;
      case ROSE :
         WriteLineDatRoSE(ofil,wf);
         break;
      case REDG :
         WriteLineDatRedDensGrad(ofil,wf);
         break;
      case SCFD :
         WriteLineDatScalarCustFld(ofil,wf);
         break;
      case VPED :
         WriteLineDatVirialPotentialEnergyDensity(ofil,wf);
         break;
      case ELLPY :
         WriteLineDatEllipticity(ofil,wf);
         break;
      default:
         cout << "Error: Field type not known!\n dat file incomplete!" << endl;
         break;
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   ofil.close();
   return;
}

#endif//_WFGRID1D_CPP_
