//
//  wfgrid1d.h
//  
//
//  Created by Juan Manuel Solano on 2013-06-03.
//
//

#ifndef _WFGRID1D_H_
#define _WFGRID1D_H_

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

#define EXTRASPACELINEFACTOR (1.0e0)

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

//**************************************************************************************************
class waveFunctionGrid1D
//**************************************************************************************************
{
public:
   //***********************************************************************************************
   waveFunctionGrid1D();
   ~waveFunctionGrid1D();
   //***********************************************************************************************
   solreal dx,maxdim;
   solreal Ca[3],Cb[3];
   string comments;
   solreal *prop1d;
   ScalarFieldType prop2plot;
   //***********************************************************************************************
   void setNPts(int nx);
   //***********************************************************************************************
   int getNPts(void);
   //***********************************************************************************************
   void setUpSimpleLine(bondNetWork &bn,int na,int nb);
   //***********************************************************************************************
   void setUpSimpleLine(bondNetWork &bn,int na);
   //***********************************************************************************************
   void setUpSimpleLine(bondNetWork &bn,solreal (&ta)[3],solreal (&tb)[3]);
   //***********************************************************************************************
   //void setUpSimpleLine(bondNetWork &bn,solreal (&ta)[3]);
   //***********************************************************************************************
   bool writeLineDatRho(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   void makeDat(string &onam,gaussWaveFunc &wf,ScalarFieldType ft);
   //***********************************************************************************************
   bool writeLineDatLapRho(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatELF(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatLOL(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatMagGradLOL(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatShannonEntropy(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatMagGradRho(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatKinetEnerDensG(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatKinetEnerDensK(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
   bool writeLineDatMolElecPot(ofstream &ofil,gaussWaveFunc &wf);
   //***********************************************************************************************
private:
   bool imsetup;
   int npts;
   //***********************************************************************************************
};
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
#endif //_WFGRID1D_H_
