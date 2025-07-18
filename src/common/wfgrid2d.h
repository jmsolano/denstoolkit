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

/* wfgrid2d.h
    Created by Juan Manuel Solano Altamirano on 10/05/13.
    Copyright 2013. All rights reserved. */


#ifndef _WFGRID2D_H_
#define _WFGRID2D_H_

#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "fldtypesdef.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef DEFAULTPOINTSPERDIRECTION
#define DEFAULTPOINTSPERDIRECTION (200)
#endif

#define EXTRASPACEPLANEFACTOR (1.0e0)

#ifndef USEPROGRESSBAR
#define USEPROGRESSBAR 0
#endif

#ifndef COLLINEAREPS
#define COLLINEAREPS (1.0e-02)
#endif

#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>
using std::string;

/* ****************************************************************************** */
class WaveFunctionGrid2D {
/* ****************************************************************************** */
public:
   /* *************************************************************************** */
   WaveFunctionGrid2D();
   ~WaveFunctionGrid2D();
   /* *************************************************************************** */
   double dircos1[3],dircos2[3],orig[3],dx[2],maxdim;
   double Ca[3],Cb[3],Cc[3],Cd[3];
   string comments;
   double *prop1d,**prop2d;
   ScalarFieldType prop2plot;
   /* *************************************************************************** */
   int GetNPts(int ii);   
   void SetNPts(int nx,int ny);
   void SetNPts(int nn);
   void SetUpSimplePlane(BondNetWork &bn,int na,int nb,int nc);
   void SetUpSimplePlane(BondNetWork &bn,int na,int nb);
   void SetUpSimplePlane(BondNetWork &bn,int na);
   void SetUpSimplePlane(BondNetWork &bn,double (&ta)[3],double (&tb)[3],double (&tc)[3]);
   void SetUpSimplePlane(BondNetWork &bn,double (&ta)[3],double (&tb)[3]);
   void SetUpSimplePlane(BondNetWork &bn,double (&ta)[3]);
   bool WritePlaneTsvRho(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvLapRho(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvELF(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvShannonEntropy(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvMagGradRho(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvLOL(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvMagGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvMolElecPot(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvLED(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvMagLED(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvRedDensMag(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvRoSE(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvDORI(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvSpinDensity(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvOneElecDiseq(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvScalarCustFld(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvVectorCustFld(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf);
   bool WritePlaneTsvEllipticity(ofstream &ofil,GaussWaveFunction &wf);
   void MakeTsv(string &onam,GaussWaveFunction &wf,ScalarFieldType ft);
/* *************************************************************************** */
private:
/* *************************************************************************** */
   bool imsetup;
   int npts[2];
/* *************************************************************************** */
};
/* *************************************************************************** */


#endif//_WFGRID2D_H_

