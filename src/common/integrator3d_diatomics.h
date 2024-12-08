/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
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
#ifndef _INTEGRATOR3D_DIATOMICS_H
#define _INTEGRATOR3D_DIATOMICS_H
#include <vector>
using std::vector;
#include <memory>
using std::shared_ptr;
#include <fstream>
using std::ofstream;
#include "integrator3d.h"

/* ************************************************************************** */
class Integrator3DDiatomics : public Integrator3D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Integrator3DDiatomics();
   Integrator3DDiatomics(shared_ptr<Function3D> i);
   virtual ~Integrator3DDiatomics();
   void ComputeIntegral();
   void DisplayProperties();
   void DisplayResults();
   void WriteProperties(ofstream &ofil);
   void WriteResults(ofstream &ofil);
   size_t NumberOfEvaluations() { return (2*(xt.size()));}
   size_t NGaussLegR() { return xrc.size(); }
   size_t NGaussLegPhi() { return xpc.size(); }
   void GetWeightsAndAbscissas(vector<double> &uw,vector<vector<double> > &ux);
/* ************************************************************************** */
   /** This function sets up the member variables using the information
    * provided by the user. The user does not need to concern about
    * which atom is located on the upper/lower semiplane (XZ plane).
    * However, the atoms will be assumed to be located along the z-axis.
    * @param[in] xx0    The coordinates of the first atom.
    * @param[in] xx1    The coordinates of the second atom.
    * @param[in] xxc    The point at which the XY plane intersects the
    *                     z-axis.
    * @param[in] rr0    The radius of the first core.
    * @param[in] rr1    The radius of the second core.
    * @param[in] Rout   The outer integration radius; it includes
    *                      both atoms and delimits the integration
    *                      circle.
    * @param[in] radord The order of the radial quadrature scheme.
    * @param[in] phiord The order of the angular quadrature scheme.*/
   bool SetupCubature(const vector<double> &xx0,const vector<double> &xx1,\
         const vector<double> &xxc,const double rr0,const double rr1,\
         const double Rout,const int radord,const int phiord);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void BuildHalfHemisphereCubature(const bool upper);
   void BuildUpperCubature();
   void BuildLowerCubature();
/* ************************************************************************** */
   vector<vector<double> > xt;/*!< Internal temporary abscissas.  */
   vector<double> wt; /*!< Internal temporary weights.  */
   vector<double> xrc;/*!< Gauss-Legendre radial core abscissas.  */
   vector<double> xrv;/*!< Gauss-Legendre radial valence abscissas.  */
   vector<double> xpc;/*!< Gauss-Legendre phi core abscissas.  */
   vector<double> xpv;/*!< Gauss-Legendre phi valence abscissas.  */
   vector<double> wrc;/*!< Gauss-Legendre radial core weights.  */
   vector<double> wrv;/*!< Gauss-Legendre radial valence weights.  */
   vector<double> wpc;/*!< Gauss-Legendre phi core weights.  */
   vector<double> wpv;/*!< Gauss-Legendre phi valence weights.  */
   double xdir[3]; /*!< The direction of the line that joins the two atoms.
                        It must be (0,0,-1).  */
   double x0[3]; /*!< The coordinates of the upper atom.  */
   double x1[3]; /*!< The coordinates of the lower atom.  */
   double xc[3]; /*!< The coordinates of the intersection of the z-axis with the XY plane.  */
   double r0; /*!< The core-radius of the upper atom.  */
   double r1; /*!< The core-radius of the lower atom.  */
   double R; /*!< The outer radius of integration (the circle of radius
                  R includes both atoms.  */
   bool firstisupper;
   bool imsetup;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATOR3D_DIATOMICS_H */

