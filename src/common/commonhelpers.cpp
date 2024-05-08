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
#include <cstdlib>
#include <iostream>
using std::cout;
#include <iomanip>
using std::scientific;
using std::setprecision;
#include "commonhelpers.h"
#include "atomcolschjmol.h"
#include "atomradiicust.h"
#include "mymath.h"
#include "matrixvectoroperations3d.h"

void CommonHelpers::PutNuclei(ofstream &ofil,BondNetWork &bn,int ntbs,\
      const string trnsmat,bool cpkview) {
   int atomn;
   double atrad;
   int nt=ntbs;
   ofil << "//          Nuclei of molecule." << '\n';
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union {" << '\n';
   for (int i=0; i<bn.nNuc; i++) {
      atomn=bn.atNum[i];
      if ( cpkview ) {
         atrad=GetAtomicVDWRadius(atomn)*AUTOMATICSPACEFILLINGRATIO;
      } else {
         atrad=bn.drawAtSize;
      }
      HelpersPOVRay::WriteTransparentSphere(ofil,nt,bn.R[i][0],bn.R[i][1],bn.R[i][2],atrad,
            GetAtomicRColorReal(atomn),GetAtomicGColorReal(atomn),
            GetAtomicBColorReal(atomn),trnsmat);
   }
   ofil << HelpersPOVRay::IndTabsStr(--nt) << "}\n";
}
void CommonHelpers::PutSpecialSpheres(ofstream &ofil,int ntbs,\
      const vector<vector<double> > &sp,const string trnsmat) {
   if ( sp[0].size()<7 ) {
      ScreenUtils::DisplayErrorMessage("CommonHelpers::PutSpecialSpheres needs 7 componets per sp[i] vector!");
      return;
   }
   double atrad;
   int nt=ntbs;
   ofil << "//          Special additional points." << '\n';
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union {" << '\n';
   for (size_t i=0; i<sp.size(); i++) {
      atrad=sp[i][3];
      HelpersPOVRay::WriteTransparentSphere(ofil,nt,sp[i][0],sp[i][1],sp[i][2],atrad,
            sp[i][4],sp[i][5],sp[i][6],trnsmat);
   }
   ofil << HelpersPOVRay::IndTabsStr(--nt) << "}\n";
}
void CommonHelpers::PutBonds(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmbnd) {
   string pigmstr="transmit ";
   pigmstr+=trnsmbnd;
   int nt=ntbs;
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union{" << '\n';
   int k=0,atni,atnk;
   double startpt[3],frak1;
   for (int i=0; i<bn.nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bn.bNet[i][j];
         atni=bn.atNum[i];
         atnk=bn.atNum[k];
         frak1=GetAtomicVDWRadius(atni)/(GetAtomicVDWRadius(atni)+GetAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=bn.R[i][l]*(1.0e0-frak1)+bn.R[k][l]*frak1;
         }
         if (k>0) {
            HelpersPOVRay::WriteCylinder(ofil,nt,
                  bn.R[i][0],bn.R[i][1],bn.R[i][2],
                  startpt[0],startpt[1],startpt[2],
                  bn.drawStickSize,
                  GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
                  GetAtomicBColorReal(atni),pigmstr);
            HelpersPOVRay::WriteCylinder(ofil,nt,
                  startpt[0],startpt[1],startpt[2],
                  bn.R[k][0],bn.R[k][1],bn.R[k][2],
                  bn.drawStickSize,
                  GetAtomicRColorReal(atnk),GetAtomicGColorReal(atnk),
                  GetAtomicBColorReal(atnk),pigmstr);
         }
      }
      HelpersPOVRay::WriteSphere(ofil,nt,bn.R[i][0],bn.R[i][1],bn.R[i][2],
            bn.drawStickSize,
            GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
            GetAtomicBColorReal(atni));
   }
   ofil << HelpersPOVRay::IndTabsStr(--nt) << "}" << '\n';
}
void CommonHelpers::WriteAngleDeclarations(ofstream &ofil,POVRayConfiguration &pvc) {
   ofil << "#declare GNUPlotAngle1=" << pvc.vecAngView[0]
        << "; // Equivalent to XAngle." << '\n';
   ofil << "#declare GNUPlotAngle2=" << pvc.vecAngView[2]
        << "; // Equivalent to ZAngle." << '\n';
   ofil << "#declare YAngle=" << pvc.vecAngView[1] << ";" << '\n';
}
void CommonHelpers::RenderPovfile(const string &povname,bool verbose) {
   string cmd="dtkpov2png "+povname;
   if ( !verbose ) { cmd+=" 2>/dev/null"; }
   system(cmd.c_str());
   if (verbose) {
      cout << "Rendering done." << '\n';
   }
}
void CommonHelpers::RotateCameraAroundLocCam(POVRayConfiguration &pvc,double angle) {
   vector<double> u(3);
   u[0]=pvc.locCam[0]; u[1]=pvc.locCam[1]; u[2]=pvc.locCam[2];
   double angrad=angle*M_PI/180.0e0;
   vector<vector<double> > R=MatrixVectorOperations3D::GetRotationMatrixAroundAxis(u,angrad);
   pvc.ApplyRotationMatrixToCameraAndLightSources(R);
}
void CommonHelpers::RotateCameraAroundUp(POVRayConfiguration &pvc,double angle) {
   vector<double> u(3);
   u[0]=pvc.vecUp[0]; u[1]=pvc.vecUp[1]; u[2]=pvc.vecUp[2];
   double angrad=angle*M_PI/180.0e0;
   vector<vector<double> > R=MatrixVectorOperations3D::GetRotationMatrixAroundAxis(u,angrad);
   pvc.ApplyRotationMatrixToCameraAndLightSources(R);
}
void CommonHelpers::RotateCameraAroundRight(POVRayConfiguration &pvc,double angle) {
   vector<double> u(3);
   u[0]=pvc.vecRight[0]; u[1]=pvc.vecRight[1]; u[2]=pvc.vecRight[2];
   double angrad=angle*M_PI/180.0e0;
   vector<vector<double> > R=MatrixVectorOperations3D::GetRotationMatrixAroundAxis(u,angrad);
   pvc.ApplyRotationMatrixToCameraAndLightSources(R);
}
double CommonHelpers::DetermineMonoatomicIntegralLimit(GaussWaveFunction &wf,\
      const char ft) {
   double rt=0.5e0;
   double rmx=-1.0e+50;
   bool momsp=((ft == 'm') || (ft == 'T') || (ft == 'k'));
   if ( momsp ) {
      while ( wf.EvalFTDensity(rt,0.0e0,0.0e0) >= 1.0e-14 ) {
         rt+=0.5e0;
      }
      if ( rmx<fabs(rt) ) { rmx=fabs(rt); }
      //cout << "rmxx: " << rt << '\n';
      rt=0.5e0;
      while ( wf.EvalFTDensity(0.0e0,rt,0.0e0) >= 1.0e-14 ) {
         rt+=0.5e0;
      }
      if ( rmx<fabs(rt) ) { rmx=fabs(rt); }
      //cout << "rmxy: " << rt << '\n';
      rt=0.5e0;
      while ( wf.EvalFTDensity(0.0e0,0.0e0,rt) >= 1.0e-14 ) {
         rt+=0.5e0;
      }
      if ( rmx<fabs(rt) ) { rmx=fabs(rt); }
      //cout << "rmxz: " << rt << '\n';
   } else {
      while ( wf.EvalDensity(0.0e0,0.0e0,rt) >= 1.0e-14 ) {
         rt+=0.5e0;
      }
      rmx=fabs(rt);
   }
   IsSphericallySymmetric(wf,(!momsp));
   return rmx;
}
void CommonHelpers::DetermineDiatomicIntegralLimits(GaussWaveFunction &wf,\
         const char ft, const vector<double> &r0, const vector<double> &r1,\
         const double a0,const double a1,vector<double> &r0mx,\
         vector<double> &r1mx,vector<double> &rmid) {
   double dir[3];
   for ( size_t i=0 ; i<3 ; ++i ) {
      dir[i]=r1[i]-r0[i];
      r0mx[i]=r0[i];
      r1mx[i]=r1[i];
   }
   double magdir=magV3(dir);
   normalizeV3(dir);
   for ( size_t i=0 ; i<3 ; ++i ) {
      rmid[i]=r0mx[i]+(dir[i])*(a0-0.5e0*(a0+a1-magdir));
   }
   if ( (ft == 'm') || (ft == 'T') || (ft == 'k') ) {
      ScreenUtils::DisplayWarningMessage("You should be using spherical integration\n"
            "domain!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      while ( wf.EvalFTDensity(r1mx[0],r1mx[1],r1mx[2]) >= 1.0e-14 ) {
         for (int i=0; i<3; ++i) { r1mx[i] += (0.10e0*dir[i]); }
      }
      while ( wf.EvalFTDensity(r0mx[0],r0mx[1],r0mx[2]) >= 1.0e-14 ) {
         for (int i=0; i<3; ++i) { r0mx[i] -= (0.10e0*dir[i]); }
      }
   } else {
      while ( wf.EvalDensity(r1mx[0],r1mx[1],r1mx[2]) >= 1.0e-14 ) {
         for (int i=0; i<3; ++i) { r1mx[i] += (0.10e0*dir[i]); }
      }
      while ( wf.EvalDensity(r0mx[0],r0mx[1],r0mx[2]) >= 1.0e-14 ) {
         for (int i=0; i<3; ++i) { r0mx[i] -= (0.10e0*dir[i]); }
      }
   }
}
bool CommonHelpers::AtomsAreZAligned(GaussWaveFunction &wf) {
   if ( wf.nNuc!=2 ) {
      ScreenUtils::DisplayErrorMessage("This function can be used for diatomic molecules only!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   double xx[3];
   for ( int i=0 ; i<3 ; ++i ) { xx[i]=wf.GetR(1,i)-wf.GetR(0,i); }
   if ( fabs(magV3(xx)-fabs(xx[2]))>0.000001 ) {
      return false;
   }
   return true;
}
bool CommonHelpers::IsSphericallySymmetric(GaussWaveFunction &wf,bool possp) {
   double rho[3];
   double rt=GetAtomicVDWRadiusAtomicUnits(wf.atCharge[0]-1);;
   string msg=" (pos): ";
   if ( possp ) {
      rho[0]=wf.EvalDensity(rt,0.0e0,0.0e0);
      rho[1]=wf.EvalDensity(0.0e0,rt,0.0e0);
      rho[2]=wf.EvalDensity(0.0e0,0.0e0,rt);
   } else {
      msg=" (mom): ";
      rt=0.50;
      rho[0]=wf.EvalFTDensity(rt,0.0e0,0.0e0);
      rho[1]=wf.EvalFTDensity(0.0e0,rt,0.0e0);
      rho[2]=wf.EvalFTDensity(0.0e0,0.0e0,rt);
   }
   bool symm=((fabs(rho[0]-rho[1])<1.0e-12) && (fabs(rho[0]-rho[2])<1.0e-12) \
           && (fabs(rho[1]-rho[2])<1.0e-12) );
   if ( !symm ) {
      ScreenUtils::DisplayWarningMessage("This wavefunction is not spherically symetric!");
      cout << scientific;
      cout << "Dxy" << msg << fabs(rho[0]-rho[1]) << '\n';
      cout << "Dxz" << msg << fabs(rho[0]-rho[2]) << '\n';
      cout << "Dyz" << msg << fabs(rho[1]-rho[2]) << '\n';
   }
   return symm;
}

