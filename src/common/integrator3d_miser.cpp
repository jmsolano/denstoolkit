/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.1
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
#include <algorithm>
using std::max;
using std::min;
#include <cmath>
#include <string>
using std::string;
#include "integrator3d_miser.h"

#define MNPT 15
#define MNBS (60)
#define MMCI_TINY 1.0e-30
#define MMCI_BIG 1.0e+30
#define PFAC 0.1e0


size_t Integrator3DMiser::neval=0;
Integrator3DMiser::Integrator3DMiser() : Integrator3D() {
   variance=-1.0e0;
   dith=0.0e0;
   domaintype=Integrator3D::DomainType::CUBOIDBOX;
   xa.resize(3);
   xb.resize(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      xa[i]=-1.5e0;
      xb[i]=1.5e0;
   }
   totalnpts=DEFAULTMISERINTEGNPTS;
}
Integrator3DMiser::Integrator3DMiser(shared_ptr<Function3D> i) : Integrator3DMiser() {
   SetIntegrand(i);
}
void Integrator3DMiser::SetXMin(const vector<double> & xm) {
   for ( size_t i=0 ; i<3 ; ++i ) { xa[i]=xm[i]; }
}
void Integrator3DMiser::SetXMax(const vector<double> & xm) {
   for ( size_t i=0 ; i<3 ; ++i ) { xb[i]=xm[i]; }
}
void Integrator3DMiser::ComputeIntegral() {
   neval=0;
   double vol=1.0e0;
   vector<double> len(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      len[i]=xb[i]-xa[i];
      vol*=len[i];
   }
   double mean,var;
   Miser(xa,xb,totalnpts,mean,var);
   variance=vol*vol*var;
   result=(vol*mean);
}
void Integrator3DMiser::Miser(const vector<double> &xa,const vector<double> &xb,\
      const size_t npts,double &mean,double &var) {
   double sum=0.0e0,summ=0.0e0,summ2=0.0e0,fval=0.0e0;
   size_t npre;
   vector<double> x(3),len(3);
   for ( size_t i=0 ; i<3 ; ++i ) { len[i]=xb[i]-xa[i]; }
   if ( npts<MNBS ) {
      summ=summ2=0.0e0;
      for ( size_t i=0 ; i<npts ; ++i ) {
         for ( size_t j=0 ; j<3 ; ++j ) { x[j]=(xa[j]+len[j]*rg.rannr()); }
         fval=integrand->f(x);
         summ+=fval;
         summ2+=(fval*fval);
         ++neval;
      }
      mean=summ/double(npts);
      var=max(MMCI_TINY,(summ2-summ*summ/double(npts))/double(npts*npts));
   } else {
      npre=max(size_t(double(npts)*PFAC),size_t(MNPT)); 
      double s;
      vector<double> rmid(3),fminl(3),fminr(3),fmaxl(3),fmaxr(3);
      for ( size_t j=0 ; j<3 ; ++j ) {
         s=copysign(dith,(rg.rannr()-0.5e0));
         rmid[j]=(0.5e0+s)*xa[j]+(0.5e0-s)*xb[j];
         fminl[j]=fminr[j]=MMCI_BIG;
         fmaxl[j]=fmaxr[j]=-MMCI_BIG;
      }
      for ( size_t i=0 ; i<npre ; ++i ) {
         for ( size_t j=0 ; j<3 ; ++j ) { x[j]=(xa[j]+len[j]*rg.rannr()); }
         fval=integrand->f(x);
         for ( size_t j=0 ; j<3 ; ++j ) {
            if ( x[j]<=rmid[j] ) {
               fminl[j]=min(fminl[j],fval);
               fmaxl[j]=max(fmaxl[j],fval);
            } else {
               fminr[j]=min(fminr[j],fval);
               fmaxr[j]=max(fmaxr[j],fval);
            }
         }
         ++neval;
      }
      double sumb=MMCI_BIG;
      size_t jb=string::npos;
      double siglb=1.0e0,sigrb=1.0e0,sigl,sigr;
      for ( size_t j=0 ; j<3 ; ++j ) {
         if ( (fmaxl[j]>fminl[j])&&(fmaxr[j]>fminr[j]) ) {
            sigl=max(MMCI_TINY,pow((fmaxl[j]-fminl[j]),(2.0e0/3.0e0)));
            sigr=max(MMCI_TINY,pow((fmaxr[j]-fminr[j]),(2.0e0/3.0e0)));
            sum=sigl+sigr;
            if ( sum<=sumb ) {
               sumb=sum;
               jb=j;
               siglb=sigl;
               sigrb=sigr;
            }
         }
      }
      if ( jb==string::npos ) { jb=size_t(3.0e0*rg.rannr()); }
      double rgl=xa[jb],rgm=rmid[jb],rgr=xb[jb],fracl=fabs((rgm-rgl)/(rgr-rgl));
      size_t nptl=size_t(MNPT+(npts-npre-2*MNPT)*fracl*siglb/(fracl*siglb+(1.0e0-fracl)*sigrb));
      size_t nptr=npts-npre-nptl;
      double meanl,varl,meanr,varr;
      vector<double> xal=xa,xbl=xb,xar=xa,xbr=xb;
      xbl[jb]=rmid[jb];
      Miser(xal,xbl,nptl,meanl,varl);
      xar[jb]=rmid[jb];
      Miser(xar,xbr,nptr,meanr,varr);
      mean=fracl*meanl+(1.0e0-fracl)*meanr;
      var=fracl*fracl*varl+(1.0e0-fracl)*(1.0e0-fracl)*varr;
   }
}
void Integrator3DMiser::DisplayProperties() {
   if ( verbosity>0 ) {
      BaseDisplayProperties();
   }
   cout << "Lower,back,left corner: " << xa[0] << ' ' << xa[1] << ' ' << xa[2] << '\n';
   cout << "Upwer,front,right corner: " << xb[0] << ' ' << xb[1] << ' ' << xb[2] << '\n';
   cout << "Property at (lower,back,left) corner: " << integrand->f(xa) << '\n';
   cout << "Property at (upper,front,right) corner: " << integrand->f(xb) << '\n';
   cout << "Number of evaluations: " << NumberOfEvaluations() << '\n';
}
void Integrator3DMiser::DisplayResults() {
   cout << "Integral: " << result << '\n';
   cout << "Variance: " << variance << '\n';
}
void Integrator3DMiser::WriteProperties(ofstream &ofil) {
   ofil << "#Integrator properties:\n";
   FileUtils::WriteScrStarLine(ofil);
   ofil << "Integrator type: Monte Carlo Miser\n";
   ofil << "Lower,back,left corner: " << xa[0] << ' ' << xa[1] << ' ' << xa[2] << '\n';
   ofil << "Upwer,front,right corner: " << xb[0] << ' ' << xb[1] << ' ' << xb[2] << '\n';
   ofil << "Property at (lower,back,left) corner: " << integrand->f(xa) << '\n';
   ofil << "Property at (upper,front,right) corner: " << integrand->f(xb) << '\n';
   ofil << "Number of evaluations: " << NumberOfEvaluations() << '\n';
   FileUtils::WriteScrStarLine(ofil);
}
void Integrator3DMiser::WriteResults(ofstream &ofil) {
   ofil << "#Results:\n";
   FileUtils::WriteScrStarLine(ofil);
   ofil << "Integral: " << result << '\n';
   ofil << "Variance: " << variance << '\n';
   FileUtils::WriteScrStarLine(ofil);
}

