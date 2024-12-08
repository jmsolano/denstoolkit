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

/** wfgrid1d.h
  Created by Juan Manuel Solano on 2013-06-03. */

#ifndef _WFGRID1D_H_
#define _WFGRID1D_H_

#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "fldtypesdef.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef DEFAULTPOINTSPERDIRECTION
#define DEFAULTPOINTSPERDIRECTION (200)
#endif

#define EXTRASPACELINEFACTOR (1.0e0)

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

/* ********************************************************************************** */
class WaveFunctionGrid1D {
/* ********************************************************************************** */
public:
   /* ******************************************************************************* */
   WaveFunctionGrid1D();
   ~WaveFunctionGrid1D();
   /* ******************************************************************************* */
   double dx,maxdim;
   double Ca[3],Cb[3];
   string comments;
   double *prop1d;
   ScalarFieldType prop2plot;
   /* ******************************************************************************* */
   void SetNPts(int nx);
   int GetNPts(void);
   void SetUpSimpleLine(BondNetWork &bn,int na,int nb);
   void SetUpSimpleLine(BondNetWork &bn,int na);
   void SetUpSimpleLine(BondNetWork &bn,double (&ta)[3],double (&tb)[3]);
   //void SetUpSimpleLine(BondNetWork &bn,double (&ta)[3]);
   bool WriteLineDatRho(ofstream &ofil,GaussWaveFunction &wf);
   void MakeDat(string &onam,GaussWaveFunction &wf,char prop);
   void MakeDat(string &onam,GaussWaveFunction &wf,ScalarFieldType ft);
   bool WriteLineDatLapRho(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatELF(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatLOL(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatMagGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatShannonEntropy(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatMagGradRho(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatMolElecPot(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatMagLED(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatRedDensGrad(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatRoSE(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatDORI(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatSpinDensity(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatOneElecDiseq(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatScalarCustFld(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf);
   bool WriteLineDatEllipticity(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
private:
   bool imsetup;
   int npts;
/* ******************************************************************************* */
};
/* ********************************************************************************** */
#endif //_WFGRID1D_H_
