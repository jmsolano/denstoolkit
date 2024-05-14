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
/* wfgrid3d.h
   Created by Juan Manuel Solano Altamirano on 06/05/13.
   Copyright 2013. All rights reserved.
*/
#ifndef _WFGRID3D_H_
#define _WFGRID3D_H_

#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "mymemory.h"
#include "fldtypesdef.h"

#ifndef DEFAULTPOINTSPERDIRECTION
#define DEFAULTPOINTSPERDIRECTION (80)
#endif
#ifndef USEPROGRESSBAR
#define USEPROGRESSBAR 0
#endif

#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>
using std::string;

/* ********************************************************************************** */
class WaveFunctionGrid3D {
/* ********************************************************************************** */
public:
/* ********************************************************************************** */
   WaveFunctionGrid3D();
   ~WaveFunctionGrid3D();
   /* ******************************************************************************* */
   double dx[3][3];
   double xin[3];
   string comments;
   double *prop1d;
   ScalarFieldType prop2plot;
   void SetUpSimpleGrid(GaussWaveFunction &wf,BondNetWork &bn);
   void SetUpSmartCuboidGrid(GaussWaveFunction &wf,BondNetWork &bn,const int nmx);
   void SetUpCenteredGrid(GaussWaveFunction &wf,BondNetWork &bn,\
         const int at1,const int at2,const double len,const int nmx);
   void SetNPts(int nx,int ny,int nz);
   void SetNPts(int nn);
   void SetExtraSpace(const double ll);
   int GetNPts(int ii);
   void WriteCubeRho(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeLapRho(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeELF(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeShannonEntropy(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeMagGradRho(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeLOL(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf);
   void MakeCube(string &onam,GaussWaveFunction &wf,ScalarFieldType ft);
   void MakeCube(string &onam,GaussWaveFunction &wf,char prop);
   void WriteCubeMagGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeMolElecPot(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeMagLED(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeRedDensGrad(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeRoSE(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeScalarCustFld(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeEllipticity(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeNCIRedDensGrad(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeNCIRho(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeDORI(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeSpinDensity(ofstream &ofil,GaussWaveFunction &wf);
   void WriteCubeOneElecDiseq(ofstream &ofil,GaussWaveFunction &wf);
/* ********************************************************************************** */
private:
   bool imsetup;
   int npts[3];
   double extraLen;
/* ********************************************************************************** */
};
/* ********************************************************************************** */
#endif//_WFGRID3D_H_


