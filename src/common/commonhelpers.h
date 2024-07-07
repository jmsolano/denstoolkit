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
#ifndef _COMMONHELPERS_H_
#define _COMMONHELPERS_H_
#include "bondnetwork.h"
#include "povraytools.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
/* ************************************************************************** */
class CommonHelpers {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   /** Adds the nuclei spheres to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmat is the string to be passed for the transparency of
    * the atoms. */
   static void PutNuclei(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmat,bool cpkview=false);
   /** Adds the bond cylinders to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmbnd is the string to be passed for the transparency of
    * the bonds. */
   static void PutBonds(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmbnd);
   /** Adds special spheres to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmbnd is the string to be passed for the transparency
    * of the spheres.
    * sp is a matrix. each row is a row-vector that has the following
    * structure:
    * sp[i][0], sp[i][1], and sp[i][2] are the coordinates of the i-th sphere centre.
    * sp[i][3] is the radius of the i-th sphere.
    * sp[i][4], sp[i][5], and sp[i][6] are the r, g, and b components of
    * the rgb colour that will be used to colour the ith sphere.  */
   static void PutSpecialSpheres(ofstream &ofil,const int ntbs,\
         const vector<vector<double> > &sp,const string trnsmbnd);
   static void WriteAngleDeclarations(ofstream &ofil,POVRayConfiguration &pvc);
   static void RenderPovfile(const string &povname,bool verbose=false,int w=1200);
   static void RotateCameraAroundLocCam(POVRayConfiguration &pvc, double angle);
   static void RotateCameraAroundUp(POVRayConfiguration &pvc, double angle);
   static void RotateCameraAroundRight(POVRayConfiguration &pvc, double angle);
   static double DetermineMonoatomicIntegralLimit(GaussWaveFunction &wf,
         const char ft);
   static void DetermineDiatomicIntegralLimits(GaussWaveFunction &wf,\
         const char ft, const vector<double> &r0, const vector<double> &r1,\
         const double a0,const double a1,vector<double> &r0mx,\
         vector<double> &r1mx,vector<double> &rmid);
   /** Checks whether two atoms are aligned along the z-axis. It only
    * works if wf has exactly two atoms. If wf has one or more than two atoms,
    * the function returns a false;  */
   static bool AtomsAreZAligned(GaussWaveFunction &wf);
   /** Checks whether a wavefunction is spherically symmetric. In this version,
    * the function calculates \f$\rho\f$ at \f$r=0.4R_{vdW}\f$ and over
    * compares the value along the x-, y- and z-axis.
    * possp is a boolean that indicates whether the symmetry must
    * be checked in position space; if possp==false, then the check
    * will be performed in momentum space.*/
   static bool IsSphericallySymmetric(GaussWaveFunction &wf,\
         bool possp=true,bool verbose=false);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */
#endif  /* _COMMONHELPERS_H_ */


