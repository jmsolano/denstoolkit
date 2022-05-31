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
/* wfgrid2d.cpp
    Created by Juan Manuel Solano Altamirano on 10/05/13.
    Copyright 2013. All rights reserved.
*/
/* Check for errors: Na=Nb=Nc */
#ifndef _WFGRID2D_CPP_
#define _WFGRID2D_CPP_
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::scientific;
using std::setprecision;
#include "wfgrid2d.h"
#include "mymemory.h"
#include "mymath.h"

WaveFunctionGrid2D::WaveFunctionGrid2D() {
   for (int i=0; i<3; i++) {
      Ca[i]=0.0e0;
      Cb[i]=0.0e0;
      Cc[i]=0.0e0;
      Cd[i]=0.0e0;
      dircos1[i]=0.0e0;
      dircos2[i]=0.0e0;
   }
   for (int i=0; i<2; i++) {
      npts[i]=DEFAULTPOINTSPERDIRECTION;
      dx[i]=1.0e0/double(DEFAULTPOINTSPERDIRECTION);
   }
   comments=string("#");
   prop1d=NULL;
   prop2d=NULL;
   prop2plot=NONE;
   imsetup=false;
}
WaveFunctionGrid2D::~WaveFunctionGrid2D() {
   MyMemory::Dealloc1DRealArray(prop1d);
   MyMemory::Dealloc2DRealArray(prop2d,npts[1]);
}
void WaveFunctionGrid2D::SetNPts(int nx,int ny) {
   npts[0]=nx;
   npts[1]=ny;
   return;
}
void WaveFunctionGrid2D::SetNPts(int nn) {
   for (int i=0; i<2; i++) {npts[i]=nn;}
   return;
}
int WaveFunctionGrid2D::GetNPts(int ii) {
   return npts[ii];   
}
void WaveFunctionGrid2D::SetUpSimplePlane(BondNetWork &bn,int na,int nb,int nc) {
   if (!(bn.ImStp())) {
      cout << "Error: Trying to use a non set-up BondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   if ((na>=bn.nNuc)||(nb>=bn.nNuc)||(nc>=bn.nNuc)) {
      cout << "Error: The atom you requested does not exist!" << endl;
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   double va[3],vb[3],vc[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i];
      vb[i]=bn.R[nb][i];
      vc[i]=bn.R[nc][i];
   }
   SetUpSimplePlane(bn,va,vb,vc);
   for (int i=0; i<2; i++) {dx[i]=2.0e0/double(npts[i]-1);}
   MyMemory::Alloc1DRealArray("prop1d",npts[1],prop1d);
   MyMemory::Alloc2DRealArray(string("prop2d"),npts[1],2,prop2d);
   imsetup=true;
   return;
}
void WaveFunctionGrid2D::SetUpSimplePlane(BondNetWork &bn,int na,int nb) {
   if (!(bn.ImStp())) {
      cout << "Error: Trying to use a non set-up BondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   if ((na>=bn.nNuc)||(nb>=bn.nNuc)) {
      cout << "Error: The atom you requested does not exist!" << endl;
      cout << "The grid could not be set up." << endl;
#if DEBUG
      cout << "From: " << __FILE__ << " at line " << __LINE__ << endl;
#endif
      return;
   }
   double va[3],vb[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i];
      vb[i]=bn.R[nb][i];
   }
   SetUpSimplePlane(bn,va,vb);
   for (int i=0; i<2; i++) {dx[i]=2.0e0/double(npts[i]-1);}
   MyMemory::Alloc1DRealArray("prop1d",npts[1],prop1d);
   MyMemory::Alloc2DRealArray(string("prop2d"),npts[1],2,prop2d);
   imsetup=true;
   return;
}
void WaveFunctionGrid2D::SetUpSimplePlane(BondNetWork &bn,int na) {
   if (!(bn.ImStp())) {
      cout << "Error: Trying to use a non set-up BondNetWork object!\n";
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
   double va[3];
   for (int i=0; i<3; i++) {
      va[i]=bn.R[na][i];
   }
   SetUpSimplePlane(bn,va);
   for (int i=0; i<2; i++) {dx[i]=2.0e0/double(npts[i]-1);}
   MyMemory::Alloc1DRealArray("prop1d",npts[1],prop1d);
   MyMemory::Alloc2DRealArray(string("prop2d"),npts[1],2,prop2d);
   imsetup=true;
   return;
}
void WaveFunctionGrid2D::SetUpSimplePlane(BondNetWork &bn,
                                          double (&ta)[3],double (&tb)[3],double (&tc)[3]) {
   /*
    Three special cases to consider:
    --> One single atom
    --> Two atoms
    --> Thre colineal atoms: (\vec{A}\cdot\vec{B})/(|\vec{A}||\vec{B}|)=1.0e0
    */
   double n[3],eta[3],xi[3];//,orig[3];
   for (int i=0; i<3; i++) {
      eta[i]=ta[i]-tb[i];
      xi[i]=ta[i]-tc[i];
      orig[i]=(ta[i]+tb[i]+tc[i])/3.0e0;
   }
   double meta,mxi;
   meta=magV3(eta);
   mxi=magV3(xi);
   normalizeV3(eta);
   normalizeV3(xi);
   double tmp=0.0e0;
   for (int i=0; i<3; i++) {tmp+=eta[i]*xi[i];}
   if (tmp<(1.0e0-COLLINEAREPS)) {
      crossProductV3(eta,xi,n);
      crossProductV3(n,eta,xi);
      normalizeV3(n);
      normalizeV3(eta);
      normalizeV3(xi);
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
      maxdim+=(EXTRASPACEPLANEFACTOR*bn.maxBondDist);
      for (int i=0; i<3; i++) {
         Ca[i]=(maxdim*(-eta[i]-xi[i]))+orig[i];
         Cb[i]=(maxdim*(eta[i]-xi[i]))+orig[i];
         Cc[i]=(maxdim*(eta[i]+xi[i]))+orig[i];
         Cd[i]=(maxdim*(-eta[i]+xi[i]))+orig[i];
         dircos1[i]=eta[i];
         dircos2[i]=xi[i];
      }
   } else {
      if (meta>mxi) {
         SetUpSimplePlane(bn,ta,tb);
      } else {
         SetUpSimplePlane(bn,ta,tc);
      }
   }
   return;
}
void WaveFunctionGrid2D::SetUpSimplePlane(BondNetWork &bn,double (&ta)[3],double (&tb)[3]) {
   double n[3],eta[3],xi[3];//,orig[3];
   for (int i=0; i<3; i++) {
      eta[i]=tb[i]-ta[i];
      orig[i]=0.5e0*(ta[i]+tb[i]);
   }
   int posmin;
   double amin=1.0e+50;
   for (int i=0; i<3; i++) {
      if (fabs(eta[i])<amin) {
         amin=fabs(eta[i]);
         posmin=i;
      }
      xi[i]=0.0e0;
   }
   xi[posmin]=1.0e0;
   /*
   int posmax;
   double amax=-1.0e+50;
   for (int i=0; i<3; i++) {
      if (fabs(eta[i])>amax) {
         amax=fabs(eta[i]);
         posmax=i;
      }
      xi[i]=1.0e0;
   }
   xi[posmax]=0.0e0;
   // */
   crossProductV3(eta,xi,n);
   crossProductV3(n,eta,xi);
   normalizeV3(n);
   normalizeV3(eta);
   normalizeV3(xi);
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
   maxdim+=(EXTRASPACEPLANEFACTOR*bn.maxBondDist);
   for (int i=0; i<3; i++) {
      Ca[i]=(maxdim*(-eta[i]-xi[i]))+orig[i];
      Cb[i]=(maxdim*(eta[i]-xi[i]))+orig[i];
      Cc[i]=(maxdim*(eta[i]+xi[i]))+orig[i];
      Cd[i]=(maxdim*(-eta[i]+xi[i]))+orig[i];
      dircos1[i]=eta[i];
      dircos2[i]=xi[i];
   }
   /*
   cout << Ca[0] << " " << Ca[1] << " " << Ca[2] << endl;
   cout << Cb[0] << " " << Cb[1] << " " << Cb[2] << endl;
   cout << Cc[0] << " " << Cc[1] << " " << Cc[2] << endl;
   cout << Cd[0] << " " << Cd[1] << " " << Cd[2] << endl;
   // */
   return;
}
void WaveFunctionGrid2D::SetUpSimplePlane(BondNetWork &bn,double (&ta)[3]) {
   double eta[3],xi[3];
   for (int i=0; i<3; i++) {
      eta[i]=0.0e0;
      xi[i]=0.0e0;
   }
   eta[0]=1.0e0;
   xi[1]=1.0e0;
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
   maxdim+=(EXTRASPACEPLANEFACTOR*bn.maxBondDist*1.5e0);
   for (int i=0; i<3; i++) {
      Ca[i]=(maxdim*(-eta[i]-xi[i]))+ta[i];
      Cb[i]=(maxdim*(eta[i]-xi[i]))+ta[i];
      Cc[i]=(maxdim*(eta[i]+xi[i]))+ta[i];
      Cd[i]=(maxdim*(-eta[i]+xi[i]))+ta[i];
      dircos1[i]=eta[i];
      dircos2[i]=xi[i];
   }
   return;
}
bool WaveFunctionGrid2D::WritePlaneTsvRho(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalDensity(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return true;
}
void WaveFunctionGrid2D::MakeTsv(string &onam,GaussWaveFunction &wf,ScalarFieldType ft) {
   if (!wf.imldd) {
      cout << "Error: trying to use a non loaded wave function object!\nNothing done!\n";
      return;
   }
   char cft=ConvertScalarFieldType2Char(ft);
   comments+=string("Property: ");
   comments+=GetFieldTypeKeyLong(cft);
   ofstream ofil;
   ofil.open(onam.c_str());
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
      ofil.close();
      return;
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch (ft) {
      case DENS:
         WritePlaneTsvRho(ofil,wf);
         break;
      case LAPD:
         WritePlaneTsvLapRho(ofil,wf);
         break;
      case ELFD:
         WritePlaneTsvELF(ofil,wf);
         break;
      case SENT:
         WritePlaneTsvShannonEntropy(ofil,wf);
         break;
      case MGRD:
         WritePlaneTsvMagGradRho(ofil,wf);
         break;
      case LOLD:
         WritePlaneTsvLOL(ofil,wf);
         break;
      case MGLD:
         WritePlaneTsvMagGradLOL(ofil,wf);
         break;
      case GLOL:
         WritePlaneTsvGradLOL(ofil,wf);
         break;
      case KEDG:
         WritePlaneTsvKinetEnerDensG(ofil,wf);
         break;
      case KEDK:
         WritePlaneTsvKinetEnerDensK(ofil,wf);
         break;
      case MEPD:
         WritePlaneTsvMolElecPot(ofil,wf);
         break;
      case LEDV :
         WritePlaneTsvLED(ofil,wf);
         break;
      case MLED :
         WritePlaneTsvMagLED(ofil,wf);
         break;
      case ROSE :
         WritePlaneTsvRoSE(ofil,wf);
         break;
      case REDG :
         WritePlaneTsvRedDensMag(ofil,wf);
         break;
      case SCFD :
         WritePlaneTsvScalarCustFld(ofil,wf);
         break;
      case VCFD :
         WritePlaneTsvVectorCustFld(ofil,wf);
         break;
      case VPED :
         WritePlaneTsvVirialPotentialEnergyDensity(ofil,wf);
         break;
      case ELLPY:
         WritePlaneTsvEllipticity(ofil,wf);
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
bool WaveFunctionGrid2D::WritePlaneTsvLapRho(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalLapRho(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvELF(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalELF(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvShannonEntropy(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalShannonEntropy(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvMagGradRho(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalMagGradRho(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvLOL(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalLOL(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvMagGradLOL(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3],gl[3],hl[3][3],lol;
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         wf.EvalHessLOL(xx,lol,gl,hl);
         prop1d[j]=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalKineticEnergyG(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalKineticEnergyK(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvGradLOL(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3],gl[3],hl[3][3],lol;
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         wf.EvalHessLOL(xx,lol,gl,hl);
         //prop1d[j]=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
         prop2d[j][0]=prop2d[j][1]=0.0e0;
         for (int p=0; p<3; p++) {
            prop2d[j][0]+=(dircos1[p]*gl[p]);
            prop2d[j][1]+=(dircos2[p]*gl[p]);
         }
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop2d[j][0] << "\t" << prop2d[j][1] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvMolElecPot(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalMolElecPot(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvLED(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3],led[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         wf.EvalLED(xx,led);
         //prop1d[j]=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
         prop2d[j][0]=prop2d[j][1]=0.0e0;
         for (int p=0; p<3; p++) {
            prop2d[j][0]+=(dircos1[p]*led[p]);
            prop2d[j][1]+=(dircos2[p]*led[p]);
         }
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop2d[j][0] << "\t" << prop2d[j][1] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvMagLED(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalMagLED(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvRedDensMag(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalReducedDensityGradient(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvRoSE(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalRoSE(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvScalarCustFld(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalCustomScalarField(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvVectorCustFld(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3],vcf[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         wf.EvalCustomVectorField(xx[0],xx[1],xx[2],vcf);
         //prop1d[j]=sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
         prop2d[j][0]=prop2d[j][1]=0.0e0;
         for (int p=0; p<3; p++) {
            prop2d[j][0]+=(dircos1[p]*vcf[p]);
            prop2d[j][1]+=(dircos2[p]*vcf[p]);
         }
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop2d[j][0] << "\t" << prop2d[j][1] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalVirialPotentialEnergyDensity(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}
bool WaveFunctionGrid2D::WritePlaneTsvEllipticity(ofstream &ofil,GaussWaveFunction &wf) {
   ofil << scientific << setprecision (10);
   double e1,e2,xx[3];
   e1=-1.0e0;
   for (int i=0; i<npts[0]; i++) {
      e2=-1.0e0;
      for (int j=0; j<npts[1]; j++) {
         for (int k=0; k<3; k++) {
            xx[k]=Ca[k]*(1.0e0-e1)*(1.0e0-e2);
            xx[k]+=Cb[k]*(1.0e0+e1)*(1.0e0-e2);
            xx[k]+=Cc[k]*(1.0e0+e1)*(1.0e0+e2);
            xx[k]+=Cd[k]*(1.0e0-e1)*(1.0e0+e2);
            xx[k]*=0.25e0;
         }
         prop1d[j]=wf.EvalEllipticity(xx[0],xx[1],xx[2]);
         e2+=dx[1];
      }
      e2=-1.0e0*maxdim;
      for (int j=0; j<npts[1]; j++) {
         ofil << e1*maxdim << "\t" << e2 << "\t" << prop1d[j] << endl;
         e2+=dx[1]*maxdim;
      }
      e1+=dx[0];
      ofil << endl;
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return false;
}

#endif//_WFGRID2D_CPP_
