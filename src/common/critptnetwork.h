/*
 *  critptnetwork.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 05/06/13.
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef _CRITPTNETWORK_H_
#define _CRITPTNETWORK_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#ifndef DISPLAYDEBUGINFOFILELINE
#define DISPLAYDEBUGINFOFILELINE (std::cout << __FILE__ << ", line: " << __LINE__ << std::endl)
#endif

#ifndef SIGNF(a)
#define SIGNF(a) ((a)>=0?(1):(-1))
#endif

#ifndef CHOOSEPOVVERSION36
#define CHOOSEPOVVERSION36 1
#endif

#ifndef MAXLOLACPSPERATOM
#define MAXLOLACPSPERATOM (8)
#endif

#ifndef MAXRHOACPSPERATOM
#define MAXRHOACPSPERATOM (2)
#endif

#ifndef EPSGRADMAG
#define EPSGRADMAG 1.00000e-14
#endif

#ifndef MAXGRADMAG
#define MAXGRADMAG 2.0e02
#endif

#ifndef MAXITERATIONACPSEARCH
#define MAXITERATIONACPSEARCH 60
#endif

#ifndef MAXITERATIONBCPSEARCH
#define MAXITERATIONBCPSEARCH 40
#endif

#ifndef MAXITERATIONRCPSEARCH
#define MAXITERATIONRCPSEARCH 60
#endif

#ifndef MAXITERATIONCCPSEARCH
#define MAXITERATIONCCPSEARCH 240
#endif

#ifndef MAXSTEPSIZEACPSEARCH
#define MAXSTEPSIZEACPSEARCH 0.3
#endif

#ifndef MAXSTEPSIZEACPRHOSEARCH
#define MAXSTEPSIZEACPRHOSEARCH 0.1
#endif

#ifndef MAXSTEPSIZEBCPSEARCH
#define MAXSTEPSIZEBCPSEARCH 0.4
#endif

#ifndef MAXSTEPSIZERCPSEARCH
#define MAXSTEPSIZERCPSEARCH 0.35
#endif

#ifndef MAXSTEPSIZECCPSEARCH
#define MAXSTEPSIZECCPSEARCH 0.3
#endif

#ifndef MAXSTEPSIZEACPLOLSEARCH
#define MAXSTEPSIZEACPLOLSEARCH 0.01
#endif

#ifndef MAXSTEPSIZEBCPLOLSEARCH
#define MAXSTEPSIZEBCPLOLSEARCH 0.1
#endif

#ifndef MAXSTEPSIZERCPLOLSEARCH
#define MAXSTEPSIZERCPLOLSEARCH 0.3
#endif

#ifndef MAXSTEPSIZECCPLOLSEARCH
#define MAXSTEPSIZECCPLOLSEARCH 0.15
#endif

#ifndef DEFAULTHGRADIENTPATHS
#define DEFAULTHGRADIENTPATHS 0.1
#endif

#ifndef ARRAYSIZEGRADPATH
#define ARRAYSIZEGRADPATH 100
#endif

#include <iostream>
using std::cout;
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

#include "fldtypesdef.h"
#include "solmath.h"
#include "solmemhand.h"
#include "bondnetwork.h"
#include "wavefunctionclass.h"

#ifndef MAXBONDINGATOMS
#define MAXBONDINGATOMS 8
#endif

#define MAXRINGATOMS 12
#define EPSRHOMAG (1.0e-10)
#define EPSFABSDIFFCOORD (1.0e-06)
#define EXTENDEDBONDDISTFACTOR (3.5e0)
#define MINRHOSIGNIFICATIVEVAL (9.0e-04)
#define BONDISTEXTCCPSEARCHFACTOR (2.0e0)
#define ATOMCRITICALPOINTSIZEFACTOR (0.2e0)

/* ************************************************************************************* */
class critPtNetWork
/* ************************************************************************************* */
{
   /* ********************************************************************************** */
public:
   /* ********************************************************************************** */
   //            Variables
   /* ********************************************************************************** */
   int nACP,nBCP,nRCP,nCCP,nBGP;
   int **atBCP; //This array contains the atoms associated with a BCP, the third number is
                // reserver to store the number of points for the gradient path (associated also
                // to the BCP) if the bond gradient paths are requested.
   solreal **RACP,**RBCP,**RRCP,**RCCP;
   solreal ***RBGP;
   solreal centMolecVec[3],RGP[100][3];
   string *lblACP,*lblBCP,*lblRCP,*lblCCP;
   /* ********************************************************************************** */
   //            Functions
   /* ********************************************************************************** */
   critPtNetWork();
   ~critPtNetWork();
   /* ********************************************************************************** */
   void setCriticalPoints(bondNetWork &bn,gaussWaveFunc &wf,ScalarFieldType ft);
   /* ********************************************************************************** */
   void removeRedundInLabel(string &lbl);
   /* ********************************************************************************** */
   string getFirstChunkOfLabel(string &lbl);
   /* ********************************************************************************** */
   int addRhoBCP(solreal (&x)[3],string &lbl,bondNetWork &bn);
   /* ********************************************************************************** */
   void addRhoRCP(solreal (&x)[3],string &lbl,bondNetWork &bn);
   /* ********************************************************************************** */
   void addRhoCCP(solreal (&x)[3],string &lbl,bondNetWork &bn);
   /* ********************************************************************************** */
   void addLOLACP(solreal (&x)[3],string &lbl,bondNetWork &bn);
   /* ********************************************************************************** */
   void addLOLBCP(solreal (&x)[3],string &lbl,bondNetWork &bn);
   /* ********************************************************************************** */
   void printAllFieldProperties(solreal &x,solreal &y,solreal &z,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void writeAllFieldProperties(ofstream &ofil,solreal &x,solreal &y,solreal &z,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void printCPProps(gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void writeCPProps(string &ofnam,string &wfnnam,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void putNuclei(ofstream &pof,bondNetWork &bn);
   /* ********************************************************************************** */
   void putBonds(ofstream &pof,bondNetWork &bn);
   /* ********************************************************************************** */
   void drawNuclei(bool dn);
   /* ********************************************************************************** */
   void drawBonds(bool db);
   /* ********************************************************************************** */
   void drawBondGradPaths(bool dbg);
   /* ********************************************************************************** */
   void tubeStyleBGP(bool stl);
   /* ********************************************************************************** */
   bool makePOVFile(string pnam,bondNetWork &bn,povRayConfProp &pvp,
                    int campos);
   /* ********************************************************************************** */
   void displayACPCoords(void);
   /* ********************************************************************************** */
   void displayBCPCoords(void);
   /* ********************************************************************************** */
   void displayRCPCoords(void);
   /* ********************************************************************************** */
   void displayCCPCoords(void);
   /* ********************************************************************************** */
   void displayIHVCoords(void);
   /* ********************************************************************************** */
   void getNextPointInGradientPathSimple(solreal (&xn)[3],solreal &stepsize,solreal &mgg,
                                      gaussWaveFunc &wf);
   /* ********************************************************************************** */
   int findGradientPathSimple(int iacp1,int iacp2,int ibcp,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   int findGradientPathRK5(int iacp1,int iacp2,int ibcp,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   int findSingleRhoGradientPathRK5(int at1,int at2,solreal hstep,
                                    int dima,solreal** (&arbgp),solreal (&ro)[3],
                                    gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void getNextPointInGradientPathRK5(solreal (&xn)[3],solreal &stepsize,solreal &mgg,
                                      gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void addGradientPathToPOVFile(int npts,ofstream &ofil);
   /* ********************************************************************************** */
   void setBondPaths(gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool seekSingleRhoBCP(int ata,int atb,gaussWaveFunc &wf,solreal (&x)[3]);
   /* ********************************************************************************** */
   void invertOrderBGPPoints(int dim,solreal** (&arr));
   /* ********************************************************************************** */
   
   /* ********************************************************************************** */
   bool iKnowACPs(void);
   /* ********************************************************************************** */
   bool iKnowBCPs(void);
   /* ********************************************************************************** */
   bool iKnowRCPs(void);
   /* ********************************************************************************** */
   bool iKnowCCPs(void);
   /* ********************************************************************************** */
   ScalarFieldType myCPType(void);
   /* ********************************************************************************** */
   bool iKnowBGPs(void);
   /* ********************************************************************************** */
   bool readFromFile(string inname);
   /* ********************************************************************************** */
   void displayStatus(bool lngdesc = false);
   /* ********************************************************************************** */
private:
   /* ********************************************************************************** */
   int dACP,dBCP,dRCP,dCCP;
   bool iknowacps,iknowbcps,iknowrcps,iknowccps, iknowallcps;
   bool iknowbgps;
   bool drawNuc,drawBnd,drawBGPs;
   bool tubeBGPStyle;
   ScalarFieldType mycptype;
   static const int nIHV=14; //It is actually the vertices of an icosahedron plus the origin
                                // (0,0,0)
   static solreal V0,V5,V8,IHV[nIHV][3];
   /* ********************************************************************************** */
   void getACPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
   /* ********************************************************************************** */
   void getBCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
   /* ********************************************************************************** */
   void getRCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
   /* ********************************************************************************** */
   void getCCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
   /* ********************************************************************************** */
   void seekRhoACP(solreal (&x)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setRhoACPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekRhoBCP(solreal (&x)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setRhoBCPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekRhoRCP(solreal (&x)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setRhoRCPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekRhoCCP(solreal (&x)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setRhoCCPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekLOLACP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setLOLACPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekLOLBCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setLOLBCPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekLOLRCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setLOLRCPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void seekLOLCCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],gaussWaveFunc &wf);
   /* ********************************************************************************** */
   bool setLOLCCPs(bondNetWork &bn,gaussWaveFunc &wf);
   /* ********************************************************************************** */
   void addArbCP(solreal (&x)[3],solreal* &arr,int dim);
   /* ********************************************************************************** */
   void putACPs(ofstream &pof);
   /* ********************************************************************************** */
   void centerMolecule(bondNetWork &bn);
   /* ********************************************************************************** */
   void invertOrderBGPPoints(int dim);
   /* ********************************************************************************** */
   void findTwoClosestAtoms(solreal (&xo)[3],gaussWaveFunc &wf,int &idx1st,\
         int &idx2nd);
   /* ********************************************************************************** */
   /* ********************************************************************************** */
   /* ********************************************************************************** */
   bool imNew(solreal (&x)[3],int dim,solreal ** (&arr),size_t &pos);
   //size_t maxiterracp,maxiterrbcp,maxiterrrcp,maxiterrccp;
   //size_t maxiterlacp,maxiterlbcp,maxiterlrcp,maxiterlccp;
   /* ********************************************************************************** */
};
/* ************************************************************************************* */
/* ************************************************************************************* */
#endif//_CRITPTNETWORK_H_

