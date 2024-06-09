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
      DIATOMIC,CYLINDRICAL_INT,INFINITE_3D_2CAVITIES,\
      SPHERE_1CAVITY_ZAL_AZSYM,SPHERE_2CAVITIES_ZAL_AZSYM};
   enum class IntegratorType {UNDEFINED,ANALYTICAL,MONTECARLO,\
      QUADRATURE,CUBATURE};
   Integrator3D();
   Integrator3D(shared_ptr<Function3D> uptr);
   void SetIntegrand(shared_ptr<Function3D> uptr);
   void SetIntegratorType(Integrator3D::IntegratorType it);
   void SetDomainType(Integrator3D::DomainType dt);
   virtual void ComputeIntegral() = 0;
   virtual void DisplayProperties();
   virtual void DisplayResults();
   virtual void WriteResults(ofstream &ofil);
   virtual size_t NumberOfEvaluations() = 0;
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

