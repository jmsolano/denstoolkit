/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.0
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



#ifndef _CRITPTNETWORK_H_
#define _CRITPTNETWORK_H_

#include <cstdlib>
#include <fstream>
using std::ofstream;
#include <cmath>

#include "fldtypesdef.h"

#ifndef CPNW_ARRAYSIZEGRADPATH
#define CPNW_ARRAYSIZEGRADPATH 100
#endif

/* ************************************************************************************ */
class critPtNetWork {
/* ************************************************************************************ */
public:
/* ************************************************************************************ */
//                    Variables
/* ************************************************************************************ */
   int nACP,nBCP,nRCP,nCCP,nBGP;
   int **atBCP; //This array contains the atoms (acps) associated with a BCP, the third number is
   // reserved to store the number of points for the gradient path (associated also
   // to the BCP) if the bond gradient paths are requested.
   // In the old version, atBCP only included atoms which were part of the wf.
   // In contrast, in this version atBCP contains indices to actual ACPs. This is
   // needed in order to correctly search BCPs and Bond paths between atoms and 
   // non-nuclear ACPs.
   solreal **RACP,**RBCP,**RRCP,**RCCP;
   solreal ***RBGP;
   solreal centMolecVec[3],RGP[CPNW_ARRAYSIZEGRADPATH][3];
   string *lblACP,*lblBCP,*lblRCP,*lblCCP;
/* ************************************************************************************ */
//                   Functions
/* ************************************************************************************ */
   critPtNetWork(class gaussWaveFunc &uwf,class bondNetWork &ubn);
   ~critPtNetWork();
/* ************************************************************************************ */
   void setMaxIterationsACP(int ii) {maxItACP=ii;}
   void setMaxIterationsBCP(int ii) {maxItBCP=ii;}
   void setMaxIterationsRCP(int ii) {maxItRCP=ii;}
   void setMaxIterationsCCP(int ii) {maxItCCP=ii;}
/* ************************************************************************************ */
   void setCriticalPoints(ScalarFieldType ft);
/* ************************************************************************************ */
   void displayXCPCoords(char cpt);
/* ************************************************************************************ */
   void displayAllCPCoords(void);
/* ************************************************************************************ */
   void displayACPCoords(void) {displayXCPCoords('a');}
/* ************************************************************************************ */
   void displayBCPCoords(void) {displayXCPCoords('b');}
/* ************************************************************************************ */
   void displayRCPCoords(void) {displayXCPCoords('r');}
/* ************************************************************************************ */
   void displayCCPCoords(void) {displayXCPCoords('c');}
/* ************************************************************************************ */
   void displayIHVCoords(void);
/* ************************************************************************************ */
   void seekRhoACPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekRhoBCPWithExtraACP(int acppos,solreal maxrad);
/* ************************************************************************************ */
   void seekRhoBCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekRhoRCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekRhoCCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekLOLACPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekLOLACPsOnASphere(int atIdx,int nDivR,int nDivT,int nDivP,solreal radmin,\
         solreal radmax);
/* ************************************************************************************ */
   void seekLOLBCPWithExtraACP(int acppos,solreal maxrad);
/* ************************************************************************************ */
   void extendedSearchCPs();
/* ************************************************************************************ */
   bool readFromFile(string inname);
/* ************************************************************************************ */
   void displayStatus(bool lngdesc = false);
/* ************************************************************************************ */
   void setBondPaths(void);
/* ************************************************************************************ */
   bool seekSingleRhoBCP(int ata,int atb,solreal (&x)[3]);
/* ************************************************************************************ */
   int findSingleRhoGradientPathRK5(int at1,int at2,solreal hstep,\
         int dima,solreal** (&arbgp),solreal (&ro)[3]);
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
   bool makePOVFile(string pnam,class povRayConfProp &pvp,int campos);
/* ************************************************************************************ */
   void drawNuclei(bool dn) {drawNuc=dn;}
/* ************************************************************************************ */
   void drawBonds(bool db) {drawBnd=db;}
/* ************************************************************************************ */
   void drawBondGradPaths(bool dbg) {drawBGPs=dbg;}
/* ************************************************************************************ */
   void tubeStyleBGP(bool stl) {tubeBGPStyle=stl;}
/* ************************************************************************************ */
   void setExtendedSearch(bool ss) {mkextsearch=ss;}
/* ************************************************************************************ */
   void writeCPProps(string &ofnam,string &wfnam);
/* ************************************************************************************ */
   void printAllFieldProperties(solreal x,solreal y,solreal z);
/* ************************************************************************************ */
   void writeAllFieldProperties(solreal x,solreal y,solreal z,ofstream &ofil);
/* ************************************************************************************ */   
   bool iKnowACPs(void) {return iknowacps;}
/* ************************************************************************************ */
   bool iKnowBCPs(void) {return iknowbcps;}
/* ************************************************************************************ */
   bool iKnowRCPs(void) {return iknowrcps;}
/* ************************************************************************************ */
   bool iKnowCCPs(void) {return iknowccps;}
/* ************************************************************************************ */
   bool iKnowBGPs(void) {return iknowbgps;}
/* ************************************************************************************ */
   ScalarFieldType myCPType(void) {return mycptype;}
/* ************************************************************************************ */
protected:
/* ************************************************************************************ */
   class gaussWaveFunc *wf;
   class bondNetWork *bn;
   int dACP,dBCP,dRCP,dCCP;
   int maxItACP,maxItBCP,maxItRCP,maxItCCP;
   int normalbcp;
   bool iknowacps,iknowbcps,iknowrcps,iknowccps, iknowallcps;
   bool iknowbgps;
   bool drawNuc,drawBnd,drawBGPs;
   bool tubeBGPStyle;
   bool mkextsearch;
   solreal stepSizeACP,stepSizeBCP,stepSizeRCP,stepSizeCCP;
   ScalarFieldType mycptype;
   static const int nIHV=16; //It is actually the vertices of an icosahedron plus the origin
   // (0,0,0)
   static solreal V0,V5,V8,IHV[nIHV][3];
/* ************************************************************************************ */
   void init();
   /** The constructor without arguments is not public. This will enforce the use
    * of the constructor for assigning the pointers of gaussWaveFunc and 
    * bondNetWork.  */
   critPtNetWork();
/* ************************************************************************************ */
   void removeRedundInLabel(string &lbl);
/* ************************************************************************************ */
   string getFirstChunkOfLabel(string &lbl);
/* ************************************************************************************ */
   inline solreal computeMagnitudeV3(solreal (&v)[3])
          {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);}
/* ************************************************************************************ */
   bool imNew(solreal (&x)[3],int dim,solreal ** (&arr),size_t &pos);
/* ************************************************************************************ */
   void getACPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
/* ************************************************************************************ */
   void getBCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
/* ************************************************************************************ */
   void getRCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
/* ************************************************************************************ */
   void getCCPStep(solreal (&g)[3],solreal (&hess)[3][3],solreal (&hh)[3],int &sig);
/* ************************************************************************************ */
   int computeSignature(solreal (&ev)[3]);
/* ************************************************************************************ */
   int computeSignature(solreal (&hh)[3][3]);
/* ************************************************************************************ */
   void seekRhoACP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekRhoBCP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekRhoRCP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekRhoCCP(solreal (&x)[3],solreal &rho2ret,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekLOLACP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekLOLBCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekLOLRCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   void seekLOLCCP(solreal (&x)[3],solreal &ll,solreal (&g)[3],int &sig);
/* ************************************************************************************ */
   bool setRhoACPs(void);
/* ************************************************************************************ */
   bool setRhoBCPs(void);
/* ************************************************************************************ */
   bool setRhoRCPs(void);
/* ************************************************************************************ */
   bool setRhoCCPs(void);
/* ************************************************************************************ */
   bool setLOLACPs(void);
/* ************************************************************************************ */
   bool setLOLBCPs(void);
/* ************************************************************************************ */
   bool setLOLRCPs(void);
/* ************************************************************************************ */
   bool setLOLCCPs(void);
/* ************************************************************************************ */
   bool addRhoACP(solreal (&x)[3],int sig,string &lbl);
/* ************************************************************************************ */
   bool addRhoBCP(solreal (&x)[3],int sig,string &lbl,int &pos);
/* ************************************************************************************ */
   bool addRhoRCP(solreal (&x)[3],int sig,string &lbl,int &pos);
/* ************************************************************************************ */
   bool addRhoCCP(solreal (&x)[3],int sig,string &lbl,int &pos);
/* ************************************************************************************ */
   void findTwoClosestAtoms(solreal (&xo)[3],int &idx1st,int &idx2nd);
/* ************************************************************************************ */
   void findTwoClosestACPs(solreal (&xo)[3],int &idx1st,int &idx2nd);
/* ************************************************************************************ */
   void invertOrderBGPPoints(int dim);
/* ************************************************************************************ */
   void invertOrderBGPPoints(int dim,solreal** (&arr));
/* ************************************************************************************ */
   void getNextPointInGradientPathRK5(solreal (&xn)[3],solreal &stepsize,solreal &mgg);
/* ************************************************************************************ */
/* ************************************************************************************ */
   void centerMolecule(void);
/* ************************************************************************************ */
   void putNuclei(ofstream &pof);
/* ************************************************************************************ */
   void putBonds(ofstream &pof);
/* ************************************************************************************ */

/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
};
/* ************************************************************************************ */


#endif  /* _CRITPTNETWORK_H_ */

