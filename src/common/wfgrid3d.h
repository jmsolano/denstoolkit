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

#include "wavefunctionclass.h"
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

//**********************************************************************************************
class waveFunctionGrid3D
//**********************************************************************************************
{
public:
   //*******************************************************************************************
   waveFunctionGrid3D();
   ~waveFunctionGrid3D();
   //*******************************************************************************************
   solreal dx[3][3];
   solreal xin[3];
   string comments;
   solreal *prop1d;
   ScalarFieldType prop2plot;
   //*******************************************************************************************
   void setUpSimpleGrid(gaussWaveFunc &wf,bondNetWork &bn);
   //*******************************************************************************************
   void setUpSmartCuboidGrid(gaussWaveFunc &wf,bondNetWork &bn,const int nmx);
   //*******************************************************************************************
   void setNPts(int nx,int ny,int nz);
   //*******************************************************************************************
   void setNPts(int nn);
   //*******************************************************************************************
   int getNPts(int ii);
   //*******************************************************************************************
   void writeCubeRho(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeLapRho(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeELF(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeShannonEntropy(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeMagGradRho(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeLOL(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeKinetEnerDensG(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeKinetEnerDensK(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void makeCube(string &onam,gaussWaveFunc &wf,ScalarFieldType ft);
   //*******************************************************************************************
   void writeCubeMagGradLOL(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
   void writeCubeMolElecPot(ofstream &ofil,gaussWaveFunc &wf);
   //*******************************************************************************************
private:
   bool imsetup;
   int npts[3];
};
//**********************************************************************************************
//**********************************************************************************************
#endif//_WFGRID3D_H_


