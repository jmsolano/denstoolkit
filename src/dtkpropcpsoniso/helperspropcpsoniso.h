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
#ifndef _HELPERSPROPCPSONISO_H_
#define _HELPERSPROPCPSONISO_H_
#include <vector>
using std::vector;
#include <memory>
using std::shared_ptr;
#include "optflags.h"
#include "bondnetwork.h"
#include "gausswavefunction.h"
#include "symmetricsurfacegrid.h"
#include "isosurface.h"

/* ************************************************************************** */
class HelpersPropCPsOnIso {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static bool HaveIncompatibleOptions(int argc,char *argv[],OptionFlags &opt);
   static void GetCenterIndexAndVectors(char *argv[],const OptionFlags &options,\
         const BondNetWork &bn,int &cat,vector<double> &xc,vector<double> &xd);
   static bool MakePovFile(const string &povname,OptionFlags &options,POVRayConfiguration &pvp,\
         BondNetWork &bn,shared_ptr<MeshGrid> grid,vector<vector<double> > &sp,\
         const string &palname="none");
   static void ProjectGridOntoIsosurface(GaussWaveFunction &wf,shared_ptr<SymmetricSurfaceGrid> g,\
         const char prop,const double iso);
   /** Looks for the isovalue along a line. The line passes through c and r0. c is
    * the reference point, and r0 is the first point at which the field is evaluated.
    * Subsequent points are ri=r0+ih(r0-c) (i=1,2,...).
    * This function assumes that the field values decrease as the point
    * moves away from r0 (along the direction r0-c).  */
   static double SearchValueAlongLineDescending(const vector<double> &c,const vector<double> &r0,\
         GaussWaveFunction &wf,vector<double> &ri,const char prop,const double iso);
   /** This function searches the critical points (local minima and local maxima)
    * using the centroid of each face (projected onto the isosurface) and
    * the vertices of each face to determine whether v(rctd) is a critical point.
    * This function must be called after calling ProjectGridOntoIsosurface  */
   static bool SearchCPsCap(shared_ptr<MeshGrid> g,GaussWaveFunction &wf,\
         vector<vector<double> > &rcp,vector<size_t> &poscp,vector<int> &sigcp,\
         vector<double> &valcp,const char prop='V');
   static bool SearchCPsIso(shared_ptr<MeshGrid> g,GaussWaveFunction &wf,\
         vector<vector<double> > &rcp,vector<size_t> &poscp,vector<int> &sigcp,\
         vector<double> &valcp,const char prop='V');
   static bool ComputeNormalsAtVertices(shared_ptr<MeshGrid> g,GaussWaveFunction &wf,const char prop='d');
   static vector<vector<double> > ComputeTextures(shared_ptr<MeshGrid> g,\
         const double valmin,const double valmax,const string &palname="blues");
   static shared_ptr<MeshGrid> BuildCapMesh(int argc,char *argv[],const OptionFlags &opt,\
         GaussWaveFunction &wf,BondNetWork &bn,const double isovalue);
   static shared_ptr<MeshGrid> BuildMeshFromCube(int argc,char *argv[],const OptionFlags &opt,\
         GaussWaveFunction &wf,BondNetWork &bn,const double isovalue);
   static shared_ptr<MeshGrid> BuildMesh(int argc,char *argv[],const OptionFlags &opt,\
         GaussWaveFunction &wf,BondNetWork &bn,const double isovalue);
   static int FindClosestAtom(const BondNetWork &bn,const vector<double> &usrr);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSPROPCPSONISO_H_ */

