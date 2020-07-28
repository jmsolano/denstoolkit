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
/* wfgrid3d.cpp
    Created by Juan Manuel Solano Altamirano on 06/05/13.
    Copyright 2013. All rights reserved.
*/
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::scientific;
#include "wfgrid3d.h"
#include "solcubetools.h"
#include "bondnetwork.h"

#ifndef PARALLELIZEDTK
#define PARALLELIZEDTK 0
#endif

WaveFunctionGrid3D::WaveFunctionGrid3D() {
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         dx[i][j]=0.0e0;
      }
      npts[i]=DEFAULTPOINTSPERDIRECTION;
      xin[i]=0.0e0;
   }
   comments=string("#");
   prop1d=NULL;
   prop2plot=NONE;
   imsetup=false;
}
WaveFunctionGrid3D::~WaveFunctionGrid3D() {
   MyMemory::Dealloc1DRealArray(prop1d);
}
void WaveFunctionGrid3D::SetUpSimpleGrid(GaussWaveFunction &wf,BondNetWork &bn) {
   if (!(bn.ImStp())) {
      cout << "Error: Trying to use a non set-up BondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
      return;
   }
   for (int i=0; i<3; i++) {
      xin[i]=bn.bbmin[i]-(EXTRASPACECUBEFACTOR*bn.maxBondDist);
      dx[i][i]=(bn.bbmax[i]-bn.bbmin[i]+(2.0e0*EXTRASPACECUBEFACTOR)*bn.maxBondDist)/double(npts[i]-1);
   }
   MyMemory::Alloc1DRealArray("prop1d",npts[2],prop1d);
   imsetup=true;
   return;
}
void WaveFunctionGrid3D::SetUpSmartCuboidGrid(GaussWaveFunction &wf,BondNetWork &bn,int nmx) {
   if (!(bn.ImStp())) {
      cout << "Error: Trying to use a non set-up BondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
      return;
   }
   double xsp=(2.0e0*EXTRASPACECUBEFACTOR)*bn.maxBondDist;
   double maxdim=-1.0e+50,tmp;
   for (int i=0; i<3; i++) {
      tmp=bn.bbmax[i]-bn.bbmin[i]+xsp;
      if (tmp>maxdim) {
         maxdim=tmp;
      }
   }
   for (int i=0; i<3; i++) {
      tmp=bn.bbmax[i]-bn.bbmin[i]+xsp;
      npts[i]=floor((tmp/maxdim)*double(nmx));
      xin[i]=bn.bbmin[i]-(EXTRASPACECUBEFACTOR*bn.maxBondDist);
      dx[i][i]=(bn.bbmax[i]-bn.bbmin[i]+xsp)/double(npts[i]-1);
   }
   MyMemory::Alloc1DRealArray("prop1d",npts[2],prop1d);
   imsetup=true;
   return;
}
void WaveFunctionGrid3D::SetUpCenteredGrid(GaussWaveFunction &wf,BondNetWork &bn,\
      const int at1,const int at2,const double len,const int nmx) {
   if (!(bn.ImStp())) {
      cout << "Error: Trying to use a non set-up BondNetWork object!\n";
      cout << "The grid could not be set up." << endl;
      return;
   }
   double x0[3];
   for ( int i=0 ; i<3 ; ++i ) { x0[i]=0.5e0*(bn.R[at1-1][i]+bn.R[at2-1][i]); }
   for (int i=0; i<3; i++) {
      xin[i]=x0[i]-0.5*len;
      dx[i][i]=len/double(npts[i]-1);
   }
   MyMemory::Alloc1DRealArray("prop1d",npts[2],prop1d);
   imsetup=true;
   return;
}
void WaveFunctionGrid3D::SetNPts(int nx,int ny,int nz) {
   npts[0]=nx;
   npts[1]=ny;
   npts[2]=nz;
   return;
}
void WaveFunctionGrid3D::SetNPts(int nn) {
   for (int i=0; i<3; i++) {npts[i]=nn;}
   return;
}
int WaveFunctionGrid3D::GetNPts(int ii) {
   return npts[ii];
}
void WaveFunctionGrid3D::WriteCubeRho(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalDensity(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeLapRho(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalLapRho(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeELF(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalELF(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeShannonEntropy(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalShannonEntropy(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeMagGradRho(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalMagGradRho(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeLOL(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
//#if PARALLELIZEDTK
//            prop1d[k]=wf.EvalLOLNew(xx,yy,zz);
//#else
            prop1d[k]=wf.EvalLOL(xx,yy,zz);
//#endif
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalKineticEnergyG(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalKineticEnergyK(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeMagGradLOL(ofstream &ofil,GaussWaveFunction &wf) {
   //double xx,yy,zz;
   //xx=xin[0];
   //yy=xin[1];
   //zz=xin[2];
   double xl[3],gl[3],hl[3][3],lol,glm;
   xl[0]=xin[0];
   for (int i=0; i<npts[0]; i++) {
      xl[1]=xin[1];
      for (int j=0; j<npts[1]; j++) {
         xl[2]=xin[2];
         for (int k=0; k<npts[2]; k++) {
            wf.EvalHessLOL(xl,lol,gl,hl);
            glm=0.0e0;
            for (int m=0; m<3; m++) {glm+=gl[m]*gl[m];}
            prop1d[k]=sqrt(glm);
            //prop1d[k]=wf.EvalMagGradRho(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            xl[2]+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         xl[1]+=dx[1][1];
      }
      xl[0]+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeMolElecPot(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalMolElecPot(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeMagLED(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalMagLED(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeRedDensGrad(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalReducedDensityGradient(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeRoSE(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalRoSE(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeScalarCustFld(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalCustomScalarField(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::MakeCube(string &onam,GaussWaveFunction &wf,ScalarFieldType ft) {
   if (!wf.imldd) {
      cout << "Error: trying to use a non loaded wave function object!\nNothing done!\n";
      return;
   }
   char cft=ConvertScalarFieldType2Char(ft);
   comments+=string("Property: ");
   comments+=GetFieldTypeKeyLong(cft);
   ofstream ofil;
   ofil.open(onam.c_str());
   WriteCubeHeader(ofil,wf.title[0],comments,npts,xin,dx,wf.nNuc,wf.atCharge,wf.R);
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch (ft) {
      case DENS:
         WriteCubeRho(ofil,wf);
         break;
      case LAPD:
         WriteCubeLapRho(ofil,wf);
         break;
      case ELFD:
         WriteCubeELF(ofil,wf);
         break;
      case SENT:
         WriteCubeShannonEntropy(ofil,wf);
         break;
      case MGRD:
         WriteCubeMagGradRho(ofil,wf);
         break;
      case LOLD:
         WriteCubeLOL(ofil,wf);
         break;
      case MGLD:
         WriteCubeMagGradLOL(ofil,wf);
         break;
      case KEDG:
         WriteCubeKinetEnerDensG(ofil,wf);
         break;
      case KEDK:
         WriteCubeKinetEnerDensK(ofil,wf);
         break;
      case MEPD:
         WriteCubeMolElecPot(ofil,wf);
         break;
      case MLED :
         WriteCubeMagLED(ofil,wf);
         break;
      case REDG :
         WriteCubeRedDensGrad(ofil,wf);
         break;
      case ROSE :
         WriteCubeRoSE(ofil,wf);
         break;
      case SCFD :
         WriteCubeScalarCustFld(ofil,wf);
         break;
      case VPED :
         WriteCubeVirialPotentialEnergyDensity(ofil,wf);
         break;
      case ELLPY :
         WriteCubeEllipticity(ofil,wf);
         break;
      case NCIS :
         WriteCubeNCIRedDensGrad(ofil,wf);
         break;
      case NCIL :
         WriteCubeNCIRho(ofil,wf);
         break;
      default:
         cout << "Error: Field type not known!\n Cube incomplete!" << endl;
         break;
   }
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   ofil.close();
   return;
}
void WaveFunctionGrid3D::WriteCubeEllipticity(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalEllipticity(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalVirialPotentialEnergyDensity(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeNCIRedDensGrad(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalNCIs(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}
void WaveFunctionGrid3D::WriteCubeNCIRho(ofstream &ofil,GaussWaveFunction &wf) {
   double xx,yy,zz;
   xx=xin[0];
   yy=xin[1];
   zz=xin[2];
   for (int i=0; i<npts[0]; i++) {
      yy=xin[1];
      for (int j=0; j<npts[1]; j++) {
         zz=xin[2];
         for (int k=0; k<npts[2]; k++) {
            prop1d[k]=wf.EvalNCILambda(xx,yy,zz);
            //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
            zz+=dx[2][2];
         }
         WriteCubeProp(ofil,npts[2],prop1d);
         yy+=dx[1][1];
      }
      xx+=dx[0][0];
#if USEPROGRESSBAR
      ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts[0]-1))));
#endif
   }
   return;
}

