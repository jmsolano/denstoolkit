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
#include "screenutils.h"
#include "integrator3d_legsphtd.h"
#include "basegausslegendre.h"
#include "basesphtdesign.h"

Integrator3DLegSphtDes::Integrator3DLegSphtDes() : Integrator3D() {
   for ( size_t i=0 ; i<xs.size() ; ++i ) { xs[i].clear(); }
   xs.clear();
   xt.clear();
   ws.clear();
   xl.clear();
   wl.clear();
   SetDomainType(DomainType::SPHERICAL_INT);
   SetIntegratorType(IntegratorType::CUBATURE);
   imsetup=false;
}
Integrator3DLegSphtDes::Integrator3DLegSphtDes(shared_ptr<Function3D> i)
   : Integrator3DLegSphtDes() {
   SetIntegrand(i);
}
Integrator3DLegSphtDes::~Integrator3DLegSphtDes() {
   for ( size_t i=0 ; i<xs.size() ; ++i ) { xs[i].clear(); }
   xs.clear();
   for ( size_t i=0 ; i<xt.size() ; ++i ) { xt[i].clear(); }
   xt.clear();
   ws.clear();
   xl.clear();
   wl.clear();
}
void Integrator3DLegSphtDes::ComputeIntegral() {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("Data is not setup!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   double sum,ri;
   size_t nt=xt.size();
   size_t nl=xl.size();
   result=0.0e0;
   for ( size_t i=0 ; i<nl ; ++i ) {
      ri=xl[i];
      ScaleAbscissas(ri);
      sum=0.0e0;
      for ( size_t j=0 ; j<nt ; ++j ) { sum+=(integrand->f(xt[j])); }
      result+=(sum*wl[i]*ri*ri);
   }
   result*=ws[0];
}
bool Integrator3DLegSphtDes::SetupCubature(const double a,const double b,\
      const int glord,const int sphtord) {
   bool havevecs=true;
   xs=BaseSphericalTDesign::GetAbscissas(sphtord);
   ws=BaseSphericalTDesign::GetWeights(sphtord);
   BaseGaussLegendre::GetWeightsAndAbscissas(wl,xl,a,b,glord);
   xt.resize(xs.size());
   for ( size_t i=0 ; i<xt.size() ; ++i ) { xt[i].resize(3); }
   imsetup=true;
   return havevecs;
}
void Integrator3DLegSphtDes::GetWeightsAndAbscissas(vector<double> &uw,\
      vector<vector<double> > &ux) {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("Data is not setup!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return;
   }
   size_t nt=xt.size();
   size_t nl=xl.size();
   if ( uw.size() != (nt*nl) || ux.size() != (nt*nl) ) {
      for ( size_t i=0 ; i<ux.size() ; ++i ) { ux[i].clear(); }
      ux.resize(nl*nt);
      for ( size_t i=0 ; i<ux.size() ; ++i ) { ux[i].resize(3); }
      uw.resize(nl*nt);
   }
   size_t offset;
   double ri,ws0i;
   for ( size_t i=0 ; i<nl ; ++i ) {
      offset=nt*i;
      ri=xl[i];
      ScaleAbscissas(ri);
      for ( size_t j=0 ; j<nt ; ++j ) {
         for ( size_t k=0 ; k<3 ; ++k ) { ux[offset+j][k]=xt[j][k]; }
      }
      ws0i=ws[0]*wl[i]*ri*ri;
      for ( size_t j=0 ; j<nt ; ++j ) { uw[offset+j]=ws0i; }
   }
}
void Integrator3DLegSphtDes::ScaleAbscissas(const double a) {
   size_t n=xt.size();
   for ( size_t i=0 ; i<n ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) { xt[i][j]=a*xs[i][j]; }
   }
}
void Integrator3DLegSphtDes::DisplayResults() {
   cout << scientific << setprecision(10);
   cout << "Integral: " << result << '\n';
}
void Integrator3DLegSphtDes::DisplayProperties() {
   cout << "Using " << wl.size() << " GLweights and " << xl.size() << " GLAbscissas" << '\n';
   cout << "Using " << ws.size() << " STDweights and " << xs.size() << " STAbscissas" << '\n';
   cout << "Total evaluation points: " << ((xl.size())*(xs.size())) << '\n';
   if ( verbosity>0 ) {
      cout << "Gauss-Legendre weights and abscissas:" << '\n';
      for ( size_t i=0 ; i<xl.size() ; ++i ) {
         cout << wl[i] << ' ' << xl[i] << '\n';
      }
      cout << "Spherical-t design weights and abscissas:" << '\n';
      for ( size_t i=0 ; i<xs.size() ; ++i ) {
         cout << ws[i] << ' ' << xs[i][0] << ' ' << xs[i][1] << ' ' << xs[i][2] << '\n';
      }
   }
}

