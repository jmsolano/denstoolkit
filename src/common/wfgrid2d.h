/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
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
 *  wfgrid2d.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 10/05/13.
 *  Copyright 2013. All rights reserved.
 *
 */


#ifndef _WFGRID2D_H_
#define _WFGRID2D_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

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


/* ****************************************************************************** */
class waveFunctionGrid2D
/* ****************************************************************************** */
{
public:
   /* *************************************************************************** */
   waveFunctionGrid2D();
   ~waveFunctionGrid2D();
   /* *************************************************************************** */
   solreal dircos1[3],dircos2[3],orig[3],dx[2],maxdim;
   solreal Ca[3],Cb[3],Cc[3],Cd[3];
   string comments;
   solreal *prop1d,**prop2d;
   ScalarFieldType prop2plot;
   /* *************************************************************************** */
   void setNPts(int nx,int ny);
   /* *************************************************************************** */
   void setNPts(int nn);
   /* *************************************************************************** */
   int getNPts(int ii);   
   /* *************************************************************************** */
   void setUpSimplePlane(bondNetWork &bn,int na,int nb,int nc);
   /* *************************************************************************** */
   void setUpSimplePlane(bondNetWork &bn,int na,int nb);
   /* *************************************************************************** */
   void setUpSimplePlane(bondNetWork &bn,int na);
   /* *************************************************************************** */
   void setUpSimplePlane(bondNetWork &bn,solreal (&ta)[3],solreal (&tb)[3],solreal (&tc)[3]);
   /* *************************************************************************** */
   void setUpSimplePlane(bondNetWork &bn,solreal (&ta)[3],solreal (&tb)[3]);
   /* *************************************************************************** */
   void setUpSimplePlane(bondNetWork &bn,solreal (&ta)[3]);
   /* *************************************************************************** */
   bool writePlaneTsvRho(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvLapRho(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvELF(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvShannonEntropy(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMagGradRho(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvLOL(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMagGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvGradLOL(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvKinetEnerDensG(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvKinetEnerDensK(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMolElecPot(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvLED(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMagLED(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvRedDensMag(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvRoSE(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvScalarCustFld(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   bool writePlaneTsvVectorCustFld(ofstream &ofil,GaussWaveFunction &wf);
   /* *************************************************************************** */
   /* *************************************************************************** */
   void makeTsv(string &onam,GaussWaveFunction &wf,ScalarFieldType ft);
   /* *************************************************************************** */
private:
   bool imsetup;
   int npts[2];
};




#endif//_WFGRID2D_H_

