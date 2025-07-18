/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
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
   /** DomainType is an enum that will be useful to identify internally
    * domains of integration. */
   enum class DomainType {UNDEFINED,CARTESIAN_INF,CUBOIDBOX,SPHERICAL_INT,\
      DIATOMIC,CYLINDRICAL_INT,\
      INFINITE_3D_2CAV_ZAL,INFINITE_3D_2CAV_ZAL_AZSYM,\
      SPHERE_1CAV_ZAL,SPHERE_1CAV_ZAL_AZSYM,\
      SPHERE_2CAV_ZAL,SPHERE_2CAV_ZAL_AZSYM};
   /** IntegratorType will be useful to identify internally the type
    * of integrator. */
   enum class IntegratorType {UNDEFINED,ANALYTICAL,MONTECARLO,\
      QUADRATURE,CUBATURE};
   Integrator3D();
   /** This constructor is useful to setup the function to be
    * integrated. Especially, when the integrator is used
    * as object. */
   Integrator3D(shared_ptr<Function3D> uptr);
   /** Sets an internal pointer to Function3D; this is useful
    * to call the function during internal calculations. */
   void SetIntegrand(shared_ptr<Function3D> uptr);
   /** To manually set the integrator type. */
   void SetIntegratorType(Integrator3D::IntegratorType it);
   /** To manually set the integration domain. */
   void SetDomainType(Integrator3D::DomainType dt);
   /** Abstract function definition; this function <strong>must be
    * implemented</strong> in each derived class. */
   virtual void ComputeIntegral() = 0;
   /** An abstract method that returns the weights and abscissas after the object
    * has been correctly instantiated and setup. */
   virtual void GetWeightsAndAbscissas(vector<double> &uw,vector<vector<double> > &ux);
   /** Displays some internal parameters used by each derived
    * integrator. */
   virtual void DisplayProperties();
   /** Displays the computed integral. Each derived class can overload
    * this function so as to display extra information. */
   virtual void DisplayResults();
   /** In essence, this function outputs the same information, but
    * it does onto a file. */
   virtual void WriteProperties(ofstream &ofil);
   virtual void WriteResults(ofstream &ofil);
   virtual size_t NumberOfEvaluations() = 0;
   double Result() {return result;}
   /** Useful to display or hide debug or details of the integrator. */
   void SetVerbosityLevel(int vv) { verbosity=vv; }
   /** This is the base class DisplayResults. */
   void BaseDisplayResults();
   /** This is the base class DisplayProperties */
   void BaseDisplayProperties();
   virtual ~Integrator3D();
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   bool havefunction; /*!< True if the pointer integrand is not nullptr. */
   bool havedomaintype; /*!< True if the domain type has been setup. */
   bool haveintegratortype; /*!< True if the the integrator type has been given */
   shared_ptr<Function3D> integrand; /*!< Pointer to the integrand. */
   double result; /*!< Where the integral result is saved. */
   Integrator3D::DomainType domaintype; /*!< Internal variable to hold the DomainType */
   Integrator3D::IntegratorType integratortype; /*!< Internal variable to hold the IntegratorType */
   int verbosity; /*!< An integer to choose the level of verbosity. */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATOR3D_H_ */

