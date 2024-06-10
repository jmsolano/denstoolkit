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
#include <string>
using std::string;
#include "integrator3d.h"

Integrator3D::Integrator3D() {
   integrand=nullptr;
   havefunction=false;
   havedomaintype=false;
   haveintegratortype=false;
   domaintype=DomainType::UNDEFINED;
   integratortype=IntegratorType::UNDEFINED;
   result=0.0e0;
   verbosity=0;
}
Integrator3D::Integrator3D(shared_ptr<Function3D> uptr) : Integrator3D() {
   if ( uptr!=nullptr ) { SetIntegrand(uptr); }
}
void Integrator3D::SetIntegrand(shared_ptr<Function3D> uptr) {
   integrand=uptr;
   havefunction=true;
}
void Integrator3D::SetIntegratorType(Integrator3D::IntegratorType it) {
   integratortype=it;
   havedomaintype=true;
}
void Integrator3D::SetDomainType(Integrator3D::DomainType dt) {
   domaintype=dt;
}
Integrator3D::~Integrator3D() {
}
void Integrator3D::DisplayProperties() {
   BaseDisplayProperties();
}
void Integrator3D::DisplayResults() {
   BaseDisplayResults();
}
void Integrator3D::BaseDisplayProperties() {
   if ( verbosity>0 ) {
      string lbl;
      switch ( integratortype ) {
         case IntegratorType::MONTECARLO :
            lbl="Monte-Carlo";
            break;
         default :
            lbl="Unknown integrator type";
            break;
      }
      cout << "Integrator type: " << lbl << '\n';
      switch ( domaintype ) {
         case DomainType::CUBOIDBOX :
            lbl="Cuboid-box";
            break;
         default :
            lbl="Unknown integrator type";
            break;
      }
      cout << "Domain type: " << lbl << '\n';
   }
}
void Integrator3D::BaseDisplayResults() {
   cout << "Integral: " << result << '\n';
}
void Integrator3D::WriteProperties(ofstream &ofil) {
   ofil << "The function WriteProperties has not been overloaded...\n";
}
void Integrator3D::WriteResults(ofstream &ofil) {
   ofil << "Integral: " << result << '\n';
}

