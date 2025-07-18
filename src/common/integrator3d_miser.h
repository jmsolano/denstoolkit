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
#ifndef _INTEGRATOR3D_MISER_H_
#define _INTEGRATOR3D_MISER_H_
#include <vector>
using std::vector;
#include <memory>
using std::shared_ptr;
#include "myrandom.h"
#include "fileutils.h"
#include "function3d.h"
#include "integrator3d.h"

#ifndef DEFAULTMISERINTEGNPTS
#define DEFAULTMISERINTEGNPTS 262144
#endif

/* ************************************************************************** */
/** The Miser Monte Carlo integral estimator. This implementation is based on
 * the original miser algorithm proposed in:
 * W. H. Press and G. R. Farrar, Recursive Stratified Sampling for
 * Multidimensional Monte Carlo Integration, Computers in Physics,
 * ?? (1990) 180 - 185.  */

class Integrator3DMiser : public Integrator3D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Integrator3DMiser();
   virtual ~Integrator3DMiser() {}
   Integrator3DMiser(shared_ptr<Function3D> i);
   void ComputeIntegral();
   void DisplayProperties();
   void DisplayResults();
   void WriteProperties(ofstream &ofil);
   void WriteResults(ofstream &ofil);
   size_t NumberOfEvaluations() { return neval; }
   void Miser(const vector<double> &xa,const vector<double> &xb,const size_t npts,double &mean,double &var);
   double Variance() {return variance;}
   void SetDith(double d) {dith=d;}
  /** Sets the number of points to be evaluated during the integration  */
   void SetNumPts(const size_t nnn) {totalnpts=nnn;} 
   void SetXMin(const vector<double> &xm);
   void SetXMax(const vector<double> &xm);
   size_t GetTotalEvaluations() {return neval;}
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   MyRandom rg;
   double variance;
   double dith;
   vector<double> xa,xb;
   size_t totalnpts;
   static size_t neval;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATOR3D_MISER_H_ */

