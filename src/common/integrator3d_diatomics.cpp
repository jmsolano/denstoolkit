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
#include <cmath>
#include "mymath.h"
#include "screenutils.h"
#include "integrator3d_diatomics.h"
#include "basegausslegendre.h"
#include "basesphtdesign.h"

Integrator3DDiatomics::Integrator3DDiatomics() : Integrator3D() {
   for ( size_t i=0 ; i<xt.size() ; ++i ) { xt[i].clear(); }
   xt.clear();
   wt.clear();
   xrc.clear();
   xrv.clear();
   xpc.clear();
   xpv.clear();
   wrc.clear();
   wrv.clear();
   wpc.clear();
   wpv.clear();
   firstisupper=true;
   SetDomainType(DomainType::DIATOMIC);
   SetIntegratorType(IntegratorType::CUBATURE);
   imsetup=false;
}
Integrator3DDiatomics::Integrator3DDiatomics(shared_ptr<Function3D> i)
   : Integrator3DDiatomics() {
   SetIntegrand(i);
}
Integrator3DDiatomics::~Integrator3DDiatomics() {
   for ( size_t i=0 ; i<xt.size() ; ++i ) { xt[i].clear(); }
   xt.clear();
   wt.clear();
   xrc.clear();
   xrv.clear();
   xpc.clear();
   xpv.clear();
   wrc.clear();
   wrv.clear();
   wpc.clear();
   wpv.clear();
}
void Integrator3DDiatomics::ComputeIntegral() {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("Data is not setup!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   BuildUpperCubature();
   double sum=0.0e0;
   size_t nn=xt.size();
   for ( size_t i=0 ; i<nn ; ++i ) {
      sum+=(wt[i]*(integrand->f(xt[i])));
   }
   BuildLowerCubature();
   for ( size_t i=0 ; i<nn ; ++i ) {
      sum+=(wt[i]*(integrand->f(xt[i])));
   }
   result=2.0e0*M_PI*sum;
}
bool Integrator3DDiatomics::SetupCubature(const vector<double> &xx0,const vector<double> &xx1,\
         const vector<double> &xxc,const double rr0,const double rr1,const double Rout,\
         const int radord,const int phiord) {
   bool havevecs=true;
   int glro=radord;
   if ( radord%2 == 1 ) {
      ++glro;
      ScreenUtils::DisplayWarningMessage("Due to technical reasons, it is preferred to\n"
            "use an even number of points for the radial integration.");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   int glpo=phiord;
   if ( phiord%2 == 1 ) {
      ++glpo;
      ScreenUtils::DisplayWarningMessage("Due to technical reasons, it is preferred to\n"
            "use an even number of points for the angular integration.");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   BaseGaussLegendre::GetWeightsAndAbscissas(wrc,xrc,glro);
   BaseGaussLegendre::GetWeightsAndAbscissas(wrv,xrv,glro/2);
   BaseGaussLegendre::GetWeightsAndAbscissas(wpc,xpc,glpo);
   BaseGaussLegendre::GetWeightsAndAbscissas(wpv,xpv,glpo/2);
   wt.resize((xrc.size())*(xpc.size())+2*((xrv.size())*(xpv.size())),0.0e0);
   xt.resize((xrc.size())*(xpc.size())+2*((xrv.size())*(xpv.size())));
   for ( size_t i=0 ; i<xt.size() ; ++i ) { xt[i].resize(3,0e0); }
   for ( int i=0 ; i<3 ; ++i ) { xdir[i]=(xx1[i]-xx0[i]); }
   for ( int i=0 ; i<3 ; ++i ) { xc[i]=xxc[i]; }
   normalizeV3(xdir);
   for ( int i=0 ; i<3 ; ++i ) { x0[i]=xx0[i]; }
   for ( int i=0 ; i<3 ; ++i ) { x1[i]=xx1[i]; }
   r0=rr0;
   r1=rr1;
   double zz[3]={0.0e0,0.0e0,1.0e0};
   if ( dotProductV3(xdir,zz)<0.0e0 ) {
      if ( verbosity>0 ) {
         cout << "The first atom (wfn) is at the top cuadrant." << '\n';
      }
      firstisupper=true;
   } else {
      if ( verbosity>0 ) {
         cout << "The first atom (wfn) is at the bottom cuadrant." << '\n';
      }
      firstisupper=false;
   }
   R=Rout;
   imsetup=true;
   return havevecs;
}
void Integrator3DDiatomics::DisplayResults() {
   cout << scientific << setprecision(10);
   cout << "Integral: " << result << '\n';
}
void Integrator3DDiatomics::DisplayProperties() {
   cout << "Using " << xrc.size() << " GLRadialCore base points." << '\n';
   cout << "Using " << xrv.size() << " GLRadialValence base points." << '\n';
   cout << "Using " << xpc.size() << " GLAngularCore base points." << '\n';
   cout << "Using " << xpv.size() << " GLAngularValence base points." << '\n';
   cout << "Total evaluation points: " << (2*(xt.size())) << '\n';
   cout << "x0: " << x0[0] << ' ' << x0[1] << ' ' << x0[2] << '\n';
   cout << "x1: " << x1[0] << ' ' << x1[1] << ' ' << x1[2] << '\n';
   cout << "xc: " << xc[0] << ' ' << xc[1] << ' ' << xc[2] << '\n';
   cout << "r0: " << r0 << '\n';
   cout << "r1: " << r1 << '\n';
   cout << " R: " << R << '\n';
}
void Integrator3DDiatomics::BuildHalfHemisphereCubature(const bool upper) {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("The Integrator3DDiatomics object is not setup!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return;
   }
   /* Core integral  */
   size_t nrc=wrc.size();
   size_t npc=wpc.size();
   double twi,apip1o2;
   double d,a;
   if ( upper ) {
      if ( firstisupper ) {
         d=fabs(x0[2]);
         a=r0;
      } else {
         d=fabs(x1[2]);
         a=r1;
      }
   } else {
      if ( firstisupper ) {
         d=fabs(x1[2]);
         a=r1;
      } else {
         d=fabs(x0[2]);
         a=r0;
      }
   }
   for ( size_t i=0 ; i<nrc ; ++i ) {
      apip1o2=0.5e0*a*(xrc[i]+1.0e0);
      twi=0.5e0*a*wrc[i]*apip1o2*apip1o2;
      for ( size_t j=0 ; j<npc ; ++j ) {
         wt[i*npc+j]=wpc[j]*twi;
         xt[i*npc+j][0]=apip1o2*sqrt(1.0e0-xpc[j]*xpc[j]);
         xt[i*npc+j][1]=0.0e0;
         xt[i*npc+j][2]=apip1o2*xpc[j]+d;
         if ( !upper ) { xt[i*npc+j][2]=-xt[i*npc+j][2]; }
      }
   }
   /* Valence M integral...  */
   size_t offset=nrc*npc;
   if ( upper ) {
      if ( firstisupper ) {
         d=x0[2]-xc[2];
         a=r0;
      } else {
         d=x1[2]-xc[2];
         a=r1;
      }
   } else {
      if ( firstisupper ) {
         d=xc[2]-x1[2];
         a=r1;
      } else {
         d=xc[2]-x0[2];
         a=r0;
      }
   }
   if ( d<0.0e0 ) {
      ScreenUtils::DisplayErrorMessage("d<0!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      d=fabs(d);
   }
   double phi1=M_PI-atan2(R,d);
   double nu1=cos(phi1);
   size_t nrv=wrv.size();
   size_t npv=wpv.size();
   double onemnu1o4wj,nuj,srterm,rhoij;
   for ( size_t j=0 ; j<npv ; ++j ) {
      onemnu1o4wj=wpv[j]*(1.0e0-nu1)*0.25e0;
      nuj=0.5e0*((1.0e0-nu1)*xpv[j]+nu1+1.0e0);
      srterm=sqrt(R*R-d*d*(1.0e0-nuj*nuj))-nuj*d;
      for ( size_t i=0 ; i<nrv ; ++i ) {
         rhoij=0.5e0*((srterm-a)*xrv[i]+srterm+a);
         wt[offset+j*nrv+i]=onemnu1o4wj*wrv[i]*(srterm-a)*rhoij*rhoij;
         xt[offset+j*nrv+i][0]=rhoij*sqrt(1.0e0-nuj*nuj);
         xt[offset+j*nrv+i][1]=0.0e0;
         if ( upper ) {
            xt[offset+j*nrv+i][2]=rhoij*nuj+d+xc[2];
         } else {
            xt[offset+j*nrv+i][2]=-rhoij*nuj-d+xc[2];
         }
      }
   }
   /* Valence N integral...  */
   offset=nrc*npc+nrv*npv;
   double twdoetc;
   double onepnu1o4wj;
   for ( size_t j=0 ; j<npv ; ++j ) {
      onepnu1o4wj=wpv[j]*(1.0e0+nu1)*0.25e0;
      nuj=0.5e0*((1.0e0+nu1)*xpv[j]+nu1-1.0e0);
      twdoetc=-a-2.0e0*d/((1.0e0+nu1)*xpv[j]+nu1-1.0e0);
      for ( size_t i=0 ; i<nrv ; ++i ) {
         rhoij=0.5e0*((-a-d/nuj)*xrv[i]+a-d/nuj);
         wt[offset+j*nrv+i]=onepnu1o4wj*twdoetc*wrv[i]*rhoij*rhoij;
         xt[offset+j*nrv+i][0]=rhoij*sqrt(1.0e0-nuj*nuj);
         xt[offset+j*nrv+i][1]=0.0e0;
         if ( upper ) {
            xt[offset+j*nrv+i][2]=rhoij*nuj+d+xc[2];
         } else {
            xt[offset+j*nrv+i][2]=-rhoij*nuj-d+xc[2];
         }
      }
   }
}
void Integrator3DDiatomics::BuildUpperCubature() {
   return BuildHalfHemisphereCubature(true);
}
void Integrator3DDiatomics::BuildLowerCubature() {
   return BuildHalfHemisphereCubature(false);
   /*
   size_t offset=nrc*rpc;
   size_t nrv=wrv.size();
   size_t npv=wpv.size();
   for ( size_t i=0 ; i<npv ; ++i ) {
      wt[offset+i*npv+j]=0.0e0;
   }
   // */
}
void Integrator3DDiatomics::GetWeightsAndAbscissas(vector<double> &uw,\
      vector<vector<double> > &ux) {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("Cannot provide weights and abscissas, since\n"
            "the integrator is not setup!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return;
   }
   BuildUpperCubature();
   size_t nn=xt.size();
   if ( uw.size()!=(2*nn) ) { uw.resize(2*nn,0.0e0); }
   if ( ux.size()!=(2*nn) ) {
      for ( size_t i=0 ; i<ux.size() ; ++i ) { ux[i].clear(); }
      ux.resize((2*nn));
      for ( size_t i=0 ; i<ux.size() ; ++i ) { ux[i].resize(3); }
   }
   for ( size_t i=0 ; i<nn ; ++i ) { uw[i]=wt[i]; }
   for ( size_t i=0 ; i<nn ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) { ux[i][j]=xt[i][j]; }
   }
   BuildLowerCubature();
   for ( size_t i=0 ; i<nn ; ++i ) { uw[nn+i]=wt[i]; }
   for ( size_t i=0 ; i<nn ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) { ux[nn+i][j]=xt[i][j]; }
   }
}

