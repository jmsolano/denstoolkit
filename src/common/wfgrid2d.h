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

#include "wavefunctionclass.h"
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
   bool writePlaneTsvRho(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvLapRho(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvELF(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvShannonEntropy(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMagGradRho(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvLOL(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMagGradLOL(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvGradLOL(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvKinetEnerDensG(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvKinetEnerDensK(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMolElecPot(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvLED(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvMagLED(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvRedDensMag(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   bool writePlaneTsvRoSE(ofstream &ofil,gaussWaveFunc &wf);
   /* *************************************************************************** */
   /* *************************************************************************** */
   void makeTsv(string &onam,gaussWaveFunc &wf,ScalarFieldType ft);
   /* *************************************************************************** */
private:
   bool imsetup;
   int npts[2];
};




#endif//_WFGRID2D_H_

