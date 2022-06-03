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
void Integrator3D::DisplayResults() {
   BaseDisplayResults();
}
void Integrator3D::BaseDisplayResults() {
   cout << "Integral: " << result << '\n';
}
void Integrator3D::DisplayProperties() {
   BaseDisplayProperties();
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
void Integrator3D::WriteResults(ofstream &ofil) {
   ofil << "Integral: " << result << '\n';
}

