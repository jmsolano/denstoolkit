/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.3.0
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */

/*
 *  wfgrid3d.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 06/05/13.
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef _WFGRID3D_H_
#define _WFGRID3D_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_ 
typedef double solreal;
//typedef float solreal;
#endif

#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "solmemhand.h"
#include "fldtypesdef.h"

#ifndef DEFAULTPOINTSPERDIRECTION
#define DEFAULTPOINTSPERDIRECTION (80)
#endif
#define EXTRASPACECUBEFACTOR (1.0e0)
#ifndef USEPROGRESSBAR
#define USEPROGRESSBAR 0
#endif

#include <iostream>
using std::cout;
using std::cin;
using std::endl;
using std::ios;
#include <fstream>
using std::fstream;
using std::ifstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <math.h>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;
using std::scientific;

/* ********************************************************************************** */
class waveFunctionGrid3D
/* ********************************************************************************** */
{
public:
   /* ******************************************************************************* */
   waveFunctionGrid3D();
   ~waveFunctionGrid3D();
   /* ******************************************************************************* */
   solreal dx[3][3];
   solreal xin[3];
   string comments;
   solreal *prop1d;
   ScalarFieldType prop2plot;
   /* ******************************************************************************* */
   void setUpSimpleGrid(GaussWaveFunction &wf,bondNetWork &bn);
   /* ******************************************************************************* */
   void setUpSmartCuboidGrid(GaussWaveFunction &wf,bondNetWork &bn,const int nmx);
   /* ******************************************************************************* */
   void setNPts(int nx,int ny,int nz);
   /* ******************************************************************************* */
   void setNPts(int nn);
   /* ******************************************************************************* */
   int getNPts(int ii);
   /* ******************************************************************************* */
   void writeCubeRho(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeLapRho(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeELF(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeShannonEntropy(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeMagGradRho(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeLOL(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void makeCube(string &onam,GaussWaveFunction &wf,ScalarFieldType ft);
   /* ******************************************************************************* */
   void writeCubeMagGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeMolElecPot(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeMagLED(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeRedDensGrad(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeRoSE(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeScalarCustFld(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeVirialPotentialEnergyDensity(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeNCIRedDensGrad(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
   void writeCubeNCIRho(ofstream &ofil,GaussWaveFunction &wf);
   /* ******************************************************************************* */
private:
   bool imsetup;
   int npts[3];
};
/* ********************************************************************************** */
/* ********************************************************************************** */
#endif//_WFGRID3D_H_


