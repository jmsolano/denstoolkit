#ifndef _INTEGRATOR3D_H_
#define _INTEGRATOR3D_H_
#include <memory>
using std::shared_ptr;
#include <fstream>
using std::ofstream;
#include "function3d.h"

/* ************************************************************************** */
class Integrator3D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   enum class DomainType {UNDEFINED,CARTESIAN_INF,CUBOIDBOX,SPHERICAL_INT,\
      CYLINDRICAL_INT,INFINITE_3D_2CAVITIES,\
      SPHERE_1CAVITY_ZAL_AZSYM,SPHERE_2CAVITIES_ZAL_AZSYM};
   enum class IntegratorType {UNDEFINED,ANALYTICAL,MONTECARLO,\
      QUADRATURE};
   Integrator3D();
   Integrator3D(shared_ptr<Function3D> uptr);
   void SetIntegrand(shared_ptr<Function3D> uptr);
   void SetIntegratorType(Integrator3D::IntegratorType it);
   void SetDomainType(Integrator3D::DomainType dt);
   virtual void ComputeIntegral() = 0;
   virtual void DisplayProperties();
   virtual void DisplayResults();
   virtual void WriteResults(ofstream &ofil);
   double Result() {return result;}
   void SetVerbosityLevel(int vv) { verbosity=vv; }
   void BaseDisplayResults();
   void BaseDisplayProperties();
   virtual ~Integrator3D();
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   bool havefunction;
   bool havedomaintype;
   bool haveintegratortype;
   shared_ptr<Function3D> integrand;
   double result;
   Integrator3D::DomainType domaintype;
   Integrator3D::IntegratorType integratortype;
   int verbosity;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATOR3D_H_ */

