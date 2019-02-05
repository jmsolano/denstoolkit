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


#ifndef CPNW_MAXBCPSCONNECTEDTORCP
#define CPNW_MAXBCPSCONNECTEDTORCP (18)
#endif

#ifndef CPNW_MAXRCPSCONNECTEDTOCCP
#define CPNW_MAXRCPSCONNECTEDTOCCP (32)
#endif

/* ************************************************************************************ */
class critPtNetWork {
/* ************************************************************************************ */
public:
/* ************************************************************************************ */
//                    Variables
/* ************************************************************************************ */
   int nACP, /*!< Number of ACPs  */
       nBCP, /*!< Number of BCPs  */
       nRCP, /*!< Number of RCPs  */
       nCCP, /*!< Number of CCPs  */
       nBGP, /*!< Number of Bond Gradient Paths.  */
       nRGP, /*!< Number of Ring Gradient Paths.  */
       nCGP; /*!< Number of Cage Gradient Paths.  */
   /** This array (connectivity of a BCP) contains the acps associated with a BCP.
    * In conBCP[i][j], i refers to the i-th BCP in the list.
    * j=0 (j=1) contains the first (second) ACP connected to the BCP.
    * j=2 is reserved to store the number of points for the gradient path
    * (associated also to the BCP) if the bond gradient paths are requested.
    * In the old version, conBCP only included atoms which were part of the wf.
    * In contrast, in this version conBCP contains indices to actual ACPs. This is
    * needed in order to correctly search BCPs and Bond paths between atoms and 
    * non-nuclear ACPs.
    */
   int **conBCP;
   /** This array (connectivity of a RCP) contains the bcps associated with a RCP.
    * In conRCP[i][j][k], i refers to the i-th RCP in the list.
    * [j=0][k] contains the list of kth-bcps possibly connected to the rcp.
    * [j=1][k] contains the list of the number of points in the Ring Path
    * that connects the i-th rcp with the k-th bcp.
    * Here Ring Path is equivalent to Ring Grad Path, and are the Gradient
    * paths that connects RCPs with BCPs (BGP connects BCPs with ACPs).
    */
   int ***conRCP;
   /** This array (connectivity of a CCP) contains the rcps associated with a CCP.
    * In conCCP[i][j][k], i refers to the i-th CCP in the list.
    * [j=0][k] contains the list of kth-rcps possibly connected to the ccp.
    * [j=1][k] contains the list of the number of points in the Cage Path
    * that connects the i-th ccp with the k-th rcp.
    * Here Cage Path is equivalent to Cage Grad Path [CGP], and are the Gradient
    * paths that connects CCPs with RCPs (BGP connects BCPs with ACPs).
    */
   int ***conCCP;
   /** This array contains the coordinates of ACPs. in RACP[i][j], the j-th cartesian
    * coordinates of the i-th ACP are stored.  */
   solreal **RACP;
   /** This array contains the coordinates of BCPs. in RBCP[i][j], the j-th cartesian
    * coordinates of the i-th ACP are stored.  */
   solreal **RBCP;
   /** This array contains the coordinates of RCPs. in RRCP[i][j], the j-th cartesian
    * coordinates of the i-th RCP are stored.  */
   solreal **RRCP;
   /** This array contains the coordinates of CCPs. in RCCP[i][j], the j-th cartesian
    * coordinates of the i-th CCP are stored.  */
   solreal **RCCP;
   /** RBGP contains the coordinates of the <b>B</b>ond <b>G</b>radient
    * <b>P</b>aths. In this implementation the bond path consists of a set
    * of points along the real bond path (which is a continuos curve in 
    * space). Each bond path is uniquely associated to a BCP.
    * In RGBP[i][n][j], the first index indicates the index of the 
    * i-th BCP, the second indicates
    * the n-th point in the bond path, and the third index is the cartesian
    * coordinate of the n-th point.  */
   solreal ***RBGP;
   /** RRGP contains the coordinates of the <b>R</b>ing <b>G</b>radient
    * <b>P</b>aths. In this implementation the ring path consists of a set
    * of BCPs associated with the RCP. Each RCP-BCP pair consists, in turn,
    * of a set of points along the real ring path (which is a continuos curve in 
    * space). Several ring paths are associated to an RCP (each one of the
    * paths associated uniquely to a BCP-RCP pair.
    * In RRGP[i][j][n][k], the first index indicates the index of the 
    * i-th RCP, the second indicates
    * the j-th BCP connected to the i-th RCP,
    * n refers to the n-tn point in the ring path, and the fourth index is the cartesian
    * coordinate of the n-th point.  */
   solreal ****RRGP;
   /** RCGP contains the coordinates of the <b>C</b>age <b>G</b>radient
    * <b>P</b>aths. In this implementation the cage path consists of a set
    * of RCPs associated with the CCP. Each CCP-RCP pair consists, in turn,
    * of a set of points along the real cage path (which is a continuos curve in 
    * space). Several cage paths are associated to an CCP (each one of the
    * paths associated uniquely to a RCP-CCP pair.
    * In RCGP[i][j][n][k], the first index indicates the index of the 
    * i-th CCP, the second indicates
    * the j-th RCP connected to the i-th CCP,
    * n refers to the n-tn point in the cage path, and the fourth index is the k-th cartesian
    * coordinate of the n-th point.  */
   solreal ****RCGP;
   /** centMolecVec is a vector used for centering the molecule. This array
    * is used for producing POV files. <b>Warning: In the current version,
    * after the POV has been recorded, the coordinates of the critical
    * points will no longer correspond to the coordinate system of the
    * wavefunction.</b>  */
   solreal centMolecVec[3];
   /** The array RGP is an auxiliary array, used for manipulating the order
    * of the points in a bond path searching procedure.  */
   solreal **RGP;
   string *lblACP, /*!< An array to store the labels of the ACPs.  */\
      *lblBCP, /*!< An array to store the labels of the BCPs. */\
      *lblRCP, /*!< An array to store the labels of the RCPs.  */\
      *lblCCP; /*!< An array to store the labels of the CCPs.  */
/* ************************************************************************************ */
//                   Functions
/* ************************************************************************************ */
   /** The only allowed constructor. It requires a GaussWaveFunction and
    * a bondNetWork objects. <b>Warning: these objects must be
    * properly initialized before passing to the critPtNetWork object.</b>  */
   critPtNetWork(class GaussWaveFunction &uwf,class bondNetWork &ubn);
   ~critPtNetWork();
/* ************************************************************************************ */
   /** Self descriptive.  */
   void setMaxIterationsACP(int ii) {maxItACP=ii;}
   /** Self descriptive.  */
   void setMaxIterationsBCP(int ii) {maxItBCP=ii;}
   /** Self descriptive.  */
   void setMaxIterationsRCP(int ii) {maxItRCP=ii;}
   /** Self descriptive.  */
   void setMaxIterationsCCP(int ii) {maxItCCP=ii;}
   /** Self descriptive  */
   void setMaxGradPathNPts(int ii) {maxGradPathNPts=ii;}
   /** Self descriptive  */
   void setStepSizeBGP(solreal ss) {stepSizeBGP=ss;}
/* ************************************************************************************ */
   /** The main public function for searching all critical points.
    * Configuration, such as requesting extended search should be
    * done before calling this function.
    * ft is used for requesting the type of critical points.
    * In current implementation, only the electron density
    * (DENS) and some of the LOL (LOLD) CPs are found.
    */
   void setCriticalPoints(ScalarFieldType ft);
/* ************************************************************************************ */
   /** It allocates the arrays to work with ACPs. It is public because for non
    * regular applications, only selected parts of the whole critical point search
    * are needed. For instance, one would be interested in searching ONLY ACPs.  */
   void setupACPs(ScalarFieldType ft);
   /** This function searchs the ACPs. For regular calculations, use setCriticalPoints()  */
   void setACPs(ScalarFieldType ft);
   /** Basically, it allocates the arrays for working with BCPs. It needs to know
    * ACPs.  */
   void setupBCPs(ScalarFieldType ft);
   /** This function search the BCPs. For regular calculations, use setCriticalPoints().  */
   void setBCPs(ScalarFieldType ft);
/* ************************************************************************************ */
   /** Displays the coordinates of the critical points of type 'cpt'.
    * @param cpt: This char parameter is used to request the critical
    * point type. It can take the values: 'a', 'b', 'r', and 'c'.
    */
   void displayXCPCoords(char cpt);
/* ************************************************************************************ */
   /** Prints to the std::cout the coordinates of all CPs found.  */
   void displayAllCPCoords(void);
/* ************************************************************************************ */
   /** Self descriptive.  */
   void displayACPCoords(void) {displayXCPCoords('a');}
/* ************************************************************************************ */
   /** Self descriptive.  */
   void displayBCPCoords(void) {displayXCPCoords('b');}
/* ************************************************************************************ */
   /** Self descriptive.  */
   void displayRCPCoords(void) {displayXCPCoords('r');}
/* ************************************************************************************ */
   /** Self descriptive.  */
   void displayCCPCoords(void) {displayXCPCoords('c');}
/* ************************************************************************************ */
   /** Prints to the std::cout the coordinates used for seeding around a point.
    * Here seeding means to set starting points for a cp search. The letters
    * IHV are the acronynm of <b>I</b>cosa<b>H</b>edron <b>V</b>ertices.
    */
   void displayIHVCoords(void);
/* ************************************************************************************ */
   /** Self descriptive. 
    * This function uses a point as a spatial reference. A icosahedron, whose vertices
    * are considered to coincide with a sphere centered at the reference point, is
    * drawn, and its vertices are used as starting points for searching critical points.
    * @param oo: The coordinates of the point around of which ACPs will be looked for.
    * @param ddxx: The radius of the sphere around the point.
    * @param blbl: The base label used for naming the critical points found within
    * this function.
    * @param nvrt: Number of vertices. If only nvrt vertices must be used.
    * if nvrt=-1, the total number of vertices will be used as seeds. This is
    * useful for partial searches. For instance, the first four 'vertices' are not
    * really vertices, but the origin, and three points along x, y, and z axis,
    * displaced a small distance away from the center.
    */
   void seekRhoACPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekRhoBCPWithExtraACP(int acppos,solreal maxrad);
/* ************************************************************************************ */
   /** Same functionality as critPtNetWork::seekRhoACPsAroundAPoint, but for BCPs.  */
   void seekRhoBCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   /** Same functionality as critPtNetWork::seekRhoACPsAroundAPoint, but for RCPs.  */
   void seekRhoRCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   /** Same functionality as critPtNetWork::seekRhoACPsAroundAPoint, but for CCPs.  */
   void seekRhoCCPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   /** Same functionality as critPtNetWork::seekRhoACPsAroundAPoint, but for LOL ACPs.  */
   void seekLOLACPsAroundAPoint(solreal const (&oo)[3],solreal const ddxx,\
         string const &blbl,int nvrt=-1);
/* ************************************************************************************ */
   void seekLOLACPsOnASphere(int atIdx,int nDivR,int nDivT,int nDivP,solreal radmin,\
         solreal radmax);
/* ************************************************************************************ */
   void seekLOLBCPWithExtraACP(int acppos,solreal maxrad);
/* ************************************************************************************ */
   void customSearchTwoACPs(int acpIdx1,int acpIdx2);
/* ************************************************************************************ */
   void customSearchCPs(solreal (&xs)[3]);
/* ************************************************************************************ */
   void extendedSearchCPs();
/* ************************************************************************************ */
   bool readFromFile(string inname);
/* ************************************************************************************ */
   void displayStatus(bool lngdesc = false);
/* ************************************************************************************ */
   /** Basically, it allocates the necessary arrays for working with bond paths.
    * It requires to know ACPs.  */
   void setupBondPaths(void);
   /** The main function to compute bond paths. It requires to know BPCs (and ACPs).
    * This function calls setupBondPaths internally. Within regular calculations
    * this function should be used.  */
   void setBondPaths(void);
/* ************************************************************************************ */
   void setRingPaths(void);
/* ************************************************************************************ */
   void setCagePaths(void);
/* ************************************************************************************ */
   bool seekSingleRhoBCP(int ata,int atb,solreal (&x)[3]);
/* ************************************************************************************ */
   int findSingleRhoBondGradientPathRK5(int at1,int at2,solreal hstep,\
         int dima,solreal** (&arbgp),solreal (&ro)[3]);
/* ************************************************************************************ */
   int findSingleRhoRingGradientPathRK5(int rcpIdx,int bcpIdxInRRGP,\
         solreal hstep,int dima,solreal** (&arrgp));
/* ************************************************************************************ */
   int findSingleRhoCageGradientPathRK5(int ccpIdx,int rcpIdxInRCGP,\
         solreal hstep,int dima,solreal** (&arrgp));
/* ************************************************************************************ */
   /** This function will follow the grandient path that starts at x_1.
    *  It will use the point x_1 as the actual starting point. This is so, because
    *  usually paths starts at some CP, thus it is not well defined the direction
    *  to follow. Once x_1 is defined, the path has a starting point.
    *  x_e is used to track the case when the path is supposed to end at a 
    *  specific point. For instance, when looking for a path that connects
    *  a BCP and an RCP; if this is not required, then a dummy array must be
    *  passed to this function. x_m will save the closest point to x_e in
    *  the path, and dm the distance from x_e to x_m.
    *  h_{step} is the maximum distance between two consecutive
    *  points in the path. dima is the dimension of the array passed to 
    *  store the coordinates of the gradient path (arrgp). maxlen is the 
    *  maximum length the path should be. uphilldir is a bool to indicate
    *  that the gradient path must be uphill; if false, then it is downhill.
    *  The function returns true when the path ends at x_e. The number of
    *  points in the path is saved in npia. If the path does not end
    *  at x_e, the function saves the closest point in the path with 
    *  respect to x_e in x_m.
    * */
   bool walkGradientPathRK5ToEndPoint(solreal (&xi)[3],solreal (&x1)[3],\
         solreal (&xe)[3],solreal (&xm)[3],solreal &dm,solreal hstep,int dima,\
         solreal** (&arrgp),int &npia,solreal maxlen,bool uphilldir);
/* ************************************************************************************ */
   void forceBCPConnectivity(int bcpIdx,int acpIdx1,int acpIdx2);
/* ************************************************************************************ */
   void correctRCPConnectivity(void);
/* ************************************************************************************ */
   void removeFromConRCP(const int rcpIdx,const int pos2rem);
/* ************************************************************************************ */
   void addToConRCP(const int rcpIdx,const int bcpIdx);
/* ************************************************************************************ */
   void addToConCCP(const int ccpIdx,const int rcpIdx);
/* ************************************************************************************ */
   void findTwoClosestAtomsToBCP(int bcpIdx,int &at1Idx,int&at2Idx);
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
   void drawRingGradPaths(bool dbg) {drawRGPs=dbg;}
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
   bool iKnowRGPs(void) {return iknowrgps;}
/* ************************************************************************************ */
   bool iKnowCGPs(void) {return iknowcgps;}
/* ************************************************************************************ */
   ScalarFieldType myCPType(void) {return mycptype;}
/* ************************************************************************************ */
   solreal getMaxBondDist() {return maxBondDist;}
/* ************************************************************************************ */
   solreal getMaxDistBetwBCPAndBCP() {return maxBCPACPDist;}
/* ************************************************************************************ */
   int getNofRingPathsOfRCP(int rcpIdx);
/* ************************************************************************************ */
   int getNofCagePathsOfCCP(int ccpIdx);
/* ************************************************************************************ */
   int getTotalNofRingPaths(void);
/* ************************************************************************************ */
   int getTotalNofCagePaths(void);
/* ************************************************************************************ */
protected:
/* ************************************************************************************ */
   class GaussWaveFunction *wf;
   class bondNetWork *bn;
   int dACP,dBCP,dRCP,dCCP;
   int maxItACP,maxItBCP,maxItRCP,maxItCCP;
   int maxGradPathNPts;
   int normalbcp;
   bool iknowacps,iknowbcps,iknowrcps,iknowccps, iknowallcps;
   bool iknowbgps,iknowrgps,iknowcgps,iknowallgps;
   bool drawNuc,drawBnd,drawBGPs,drawRGPs,drawCGPs;
   bool tubeBGPStyle;
   bool mkextsearch;
   solreal stepSizeACP,stepSizeBCP,stepSizeRCP,stepSizeCCP;
   solreal stepSizeBGP;
   ScalarFieldType mycptype;
   solreal maxBondDist; /*!< The maximum distance between two ACPs related by a BCP  */
   solreal maxBCPACPDist; /*!< The maximum distance between a BCP and associated ACPs  */
   static const int nIHV=16; //It is actually the vertices of an icosahedron plus the origin
   // (0,0,0)
   static solreal V0,V5,V8,IHV[nIHV][3];
/* ************************************************************************************ */
   void init();
   /** The constructor without arguments is not public. This will enforce the use
    * of the constructor for assigning the pointers of GaussWaveFunction and 
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
   void addBCP2ConRCP(const int rcpIdx,const int bcpIdx);
/* ************************************************************************************ */
   void addRCP2ConCCP(const int ccpIdx,const int rcpIdx);
/* ************************************************************************************ */
   void invertOrderBGPPoints(int dim);
/* ************************************************************************************ */
   void invertOrderBGPPoints(int dim,solreal** (&arr));
/* ************************************************************************************ */
   void getNextPointInGradientPathRK5UpHill(solreal (&xn)[3],\
         solreal &stepsize,solreal &mgg);
/* ************************************************************************************ */
   void getNextPointInGradientPathRK5DownHill(solreal (&xn)[3],\
         solreal &stepsize,solreal &mgg);
/* ************************************************************************************ */
   void centerMolecule(void);
/* ************************************************************************************ */
   void putNuclei(ofstream &pof);
/* ************************************************************************************ */
   void putBonds(ofstream &pof);
/* ************************************************************************************ */
   void findMaxBondDist();
/* ************************************************************************************ */
   void copyRGP2Array(solreal** (&thearr),int nn);
/* ************************************************************************************ */
/* ************************************************************************************ */
};
/* ************************************************************************************ */


#endif  /* _CRITPTNETWORK_H_ */

