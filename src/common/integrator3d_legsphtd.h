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
#ifndef _INTEGRATOR3D_LEGSPHTD_H_
#define _INTEGRATOR3D_LEGSPHTD_H_
#include <vector>
using std::vector;
#include <memory>
using std::shared_ptr;
#include "integrator3d.h"

/* ************************************************************************** */
/** The class Integrator3DLegSphtDes is designed to  integrate
 * a Function3D over the sphere of over a hollowed sphere.
 * After the initial construction, use SetupCubature to
 * provide the domain.  */
class Integrator3DLegSphtDes : public Integrator3D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Integrator3DLegSphtDes();
   Integrator3DLegSphtDes(shared_ptr<Function3D> i);
   virtual ~Integrator3DLegSphtDes();
   void ComputeIntegral();
   void DisplayResults();
   void DisplayProperties();
   size_t NumberOfEvaluations() { return ((xl.size())*(xs.size())); }
   size_t NGaussLegendre() { return xl.size(); }
   size_t NSpherTDes() { return xs.size(); }
/* ************************************************************************** */
   /** This function prepares the internal weights and abscissas,
    * in order to perform the following integral:
    * \f$\int_0^{2\pi}d\varphi\int_0^{\pi}\sin\theta d\theta
    *    \int_a^bf(r,\theta\phi)\f$.
    * @param[in]: a       the lower radius bound
    * @param[in]: b       the upper radius bound
    * @param[in]: glord   the radial Gauss-Legendre order.
    * @param[in]: sphtord the spherical t design order. */
   bool SetupCubature(const double a,const double b,\
         const int glord,const int sphtord);
   void GetWeightsAndAbscissas(vector<double> &uw,vector<vector<double> > &ux);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   vector<vector<double> > xs;/*!< Spherical-t design abscissas.  */
   vector<vector<double> > xt;/*!< Internal temporary abscissas.  */
   vector<double> ws;/*!< Spherical-t design weights.  */
   vector<double> xl;/*!< Gauss-Legendre abscissas.  */
   vector<double> wl;/*!< Gauss-Legendre weights.  */
   void ScaleAbscissas(const double a);
   bool imsetup;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATOR3D_LEGSPHTD_H_ */

