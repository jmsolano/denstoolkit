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
using std::endl;
using std::cerr;
#include <cmath>
#include <fstream>
using std::ofstream;
#include <iomanip>
#include "circledots3d.h"
#include "screenutils.h"
#include "stringtools.h"
#include "mymemory.h"
#include "mymath.h"
#include "fileutils.h"


double const CircleDots3D::twoPi=6.28318530717958647692529e0;
CircleDots3D::CircleDots3D() {
   Init();
}
CircleDots3D::~CircleDots3D() {
   MyMemory::Dealloc2DRealArray(xx_,npts_);
}
void CircleDots3D::Init() {
   npts_=0;
   for ( int i=0 ; i<3 ; ++i ) {oo_[i]=ue1_[i]=ue2_[i]=ue3_[i]=0.0e0;}
   radius_=1.0e0;
   xx_=NULL;
   imsetup=havee1=havee2=havee3=false;
}
double CircleDots3D::GetCartCoord(const int i,const int j) {
   if ( (xx_!=NULL) && (i<npts_) && (j<3) ) {
      return xx_[i][j];
   } else {
      if ( xx_==NULL ) {
         ScreenUtils::DisplayErrorMessage("First setup this object");
#if DEBUG
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
         return 0.0e0;
      }
      if ( i>=npts_ ) {
         ScreenUtils::DisplayErrorMessage(string("Only ")+StringTools::GetStringFromInt(i)\
               +string(" points in the circle!"));
#if DEBUG
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
            return 0.0e0;
      }
      if ( j>2 ) {
         ScreenUtils::DisplayErrorMessage("Non valid cartesian index!");
#if DEBUG
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
           return 0.0e0;
      }
   }
   return 0.0e0;
}
double CircleDots3D::GetPhi(const int i) {
   if ( (xx_!=NULL) && (i<npts_) ) {
      return xx_[i][3];
   } else {
      if ( xx_==NULL ) {
         ScreenUtils::DisplayErrorMessage("First setup this object");
#if DEBUG
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
         return 0.0e0;
      }
      if ( i>=npts_ ) {
         ScreenUtils::DisplayErrorMessage(string("Only ")+StringTools::GetStringFromInt(i)\
               +string(" points in the circle!"));
#if DEBUG
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
            return 0.0e0;
      }
   }
   return 0.0e0;
}
void CircleDots3D::SetE1(const double x,const double y,const double z) {
   ue1_[0]=x; ue1_[1]=y; ue1_[2]=z;
   if ( sqrt(ue1_[0]*ue1_[0]+ue1_[1]*ue1_[1]+ue1_[2]*ue1_[2])<=0.0e0 ) {
      ScreenUtils::DisplayErrorMessage("Please provide a non-zero vector!");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
      for ( int i=0 ; i<3 ; ++i ) { ue1_[i]=0.0e0; }
      return; 
   } else { cout << "magE1: " << magV3(ue1_) << endl;}
   havee1=true;
}
void CircleDots3D::SetE2(const double x,const double y,const double z) {
   ue2_[0]=x; ue2_[1]=y; ue2_[2]=z;
   if ( sqrt(ue2_[0]*ue2_[0]+ue2_[1]*ue2_[1]+ue2_[2]*ue2_[2])<=0.0e0 ) {
      ScreenUtils::DisplayErrorMessage("Please provide a non-zero vector!");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
      for ( int i=0 ; i<3 ; ++i ) { ue2_[i]=0.0e0; }
      return; 
   }
   havee2=true;
}
void CircleDots3D::SetE1AndE2(const double (&ee1)[3],const double (&ee2)[3]) {
   for ( int i=0 ; i<3 ; ++i ) {
      ue1_[i]=ee1[i];
      ue2_[i]=ee2[i];
   }
   havee1=havee2=true;
}
void CircleDots3D::ComputeUE3(void) {
   if ( !(havee1&&havee2) ) {
      ScreenUtils::DisplayErrorMessage("First set e1 and e2!");
#if DEBUG
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
         return;
   }
   crossProductV3(ue1_,ue2_,ue3_);
   havee3=true;
}
void CircleDots3D::SetupCircle(void) {
   ComputeUE3();
   if ( !havee3 ) {
      ScreenUtils::DisplayErrorMessage("First setup the e3 vector!");
      return;
   }
   if ( npts_<=0 ) {
      ScreenUtils::DisplayErrorMessage("First set the number of points for the discretized"
            "circle!");
      return;
   }
   for ( int i=0 ; i<3 ; ++i ) {
      e1_[i]=ue1_[i]; normalizeV3(e1_);
      e3_[i]=ue3_[i]; normalizeV3(e3_);
      crossProductV3(e3_,e1_,e2_); normalizeV3(e2_);
   }
   dphi_=twoPi/double(npts_-1);
   MyMemory::Alloc2DRealArray("xx_",npts_,4,xx_);
   double phi,cp,sp;
   for ( int i=0 ; i<npts_ ; ++i ) {
      phi=double(i)*dphi_;
      cp=cos(phi)*radius_; sp=sin(phi)*radius_;
      for ( int j=0 ; j<3 ; ++j ) {
         xx_[i][j] =(cp*e1_[j]);
         xx_[i][j]+=(sp*e2_[j]);
      }
      xx_[i][3]=phi;
   }
   imsetup=true;
}
void CircleDots3D::DisplayCoordinates(void) {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("The circle is not setup!");
#if DEBUG
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif /* ( DEBUG ) */
      return;
   }
   cout << std::scientific << std::setprecision(10);
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::CenterString("Coordinates of circle.");
   ScreenUtils::PrintScrCharLine('-');
   ScreenUtils::PrintV3Comp("e1: ",e1_);
   ScreenUtils::PrintV3Comp("e2: ",e2_);
   ScreenUtils::PrintV3Comp("e3: ",e3_);
   ScreenUtils::PrintScrStarLine();
   int nnww=floor(log10(double(npts_-1)))+1;
   for ( int i=0 ; i<npts_ ; ++i ) {
      cout << "P[" << std::setw(nnww) << std::setfill('0') << i << "]: ";
      for ( int j=0 ; j<4 ; ++j ) { cout << std::setw(18) << xx_[i][j]; }
      cout << endl;
   }
   cout << "width: " << nnww << endl;
   ScreenUtils::PrintScrStarLine();
}
void CircleDots3D::WriteCoordinates(const string &oname,bool wrtoo) {
   ofstream ofil(oname.c_str());
   FileUtils::WriteScrStarLine(ofil);
   FileUtils::WriteCenteredString(ofil,"Coordinates of circle, centered at");
   ofil << std::scientific << std::setprecision(12);
   ofil << "#\n# " << oo_[0] << " " << oo_[1] << " " << oo_[2] << "\n#" << endl;
   ofil << "#The circle will be draw in the plane spanned by: " << endl;
   FileUtils::WriteV3Components(ofil,"#v1: ",ue1_);
   FileUtils::WriteV3Components(ofil,"#v2: ",ue2_);
   ofil << "#Final unit vectors:" << endl;
   FileUtils::WriteV3Components(ofil,"#e1: ",e1_);
   FileUtils::WriteV3Components(ofil,"#e2: ",e2_);
   FileUtils::WriteV3Components(ofil,"#e3: ",e3_);
   if ( wrtoo ) {
      ofil << "#The first point is the center of the circle." << endl;
      FileUtils::WriteV3Components(ofil,oo_);
   }
   for ( int i=0 ; i<npts_ ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) {ofil << xx_[i][j] << " ";}
      ofil << xx_[i][3] << endl;
   }
   ofil.close();
}
void CircleDots3D::SetOrigin(const double x,const double y,const double z) {
   oo_[0]=x; oo_[1]=y; oo_[2]=z;
}

