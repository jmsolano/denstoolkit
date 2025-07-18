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
#ifndef _HELPERSINTEGRATE_H_
#define _HELPERSINTEGRATE_H_
#include <memory>
using std::shared_ptr;
#include "optflags.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "integrator3d.h"

/* ************************************************************************** */
class FactoryIntegrator {
/* ************************************************************************** */
public:
   static shared_ptr<Integrator3D> CreateIntegrator(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw);
   static shared_ptr<Integrator3D> CreateIntegratorVegas(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw);
   static shared_ptr<Integrator3D> CreateIntegratorCubLegSphtDes(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw);
   static shared_ptr<Integrator3D> CreateIntegratorMiser(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw);
   static shared_ptr<Integrator3D> CreateIntegratorDiatomics(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw);
/* ************************************************************************** */
   /** Determines the integration limits (i.e., the bounding box), if the
    * system is a single atom.  */
   static void FindIntegralLimits(OptionFlags &options,char*argv[],\
         GaussWaveFunction &wf,BondNetWork &bn,char ft,vector<double> &rmin,vector<double> &rmax);
   /** Determines the integral limits to be used when the system is a diatomic
    * molecule.
    * @param [in]: ft is the field type
    * @param [out]: r0max is the maximum distance to be included in the
    *                     integral (away from the first atom).
    * @param [out]: r1max is the maximum distance to be included in the
    *                     integral (away from the second atom).
    * @param [out]: rmid  is the distance from the global coordinate
    *                     origin where the plane that divides the
    *                     integration domain is located.
    * r1max-r0max is directed along the line that joins the two atoms,
    * which must be along the z-axis. */
   static void DetermineDiatomicIntegralLimits(GaussWaveFunction &wf,\
         char ft,vector<double> &r0mx,vector<double> &r1mx,vector<double> &rmid);
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSINTEGRATE_H_ */


