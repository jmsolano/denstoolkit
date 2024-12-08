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
class CritPtNetWork {
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
   double **RACP;
   /** This array contains the coordinates of BCPs. in RBCP[i][j], the j-th cartesian
    * coordinates of the i-th ACP are stored.  */
   double **RBCP;
   /** This array contains the coordinates of RCPs. in RRCP[i][j], the j-th cartesian
    * coordinates of the i-th RCP are stored.  */
   double **RRCP;
   /** This array contains the coordinates of CCPs. in RCCP[i][j], the j-th cartesian
    * coordinates of the i-th CCP are stored.  */
   double **RCCP;
   /** RBGP contains the coordinates of the <b>B</b>ond <b>G</b>radient
    * <b>P</b>aths. In this implementation the bond path consists of a set
    * of points along the real bond path (which is a continuos curve in 
    * space). Each bond path is uniquely associated to a BCP.
    * In RGBP[i][n][j], the first index indicates the index of the 
    * i-th BCP, the second indicates
    * the n-th point in the bond path, and the third index is the cartesian
    * coordinate of the n-th point.  */
   double ***RBGP;
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
   double ****RRGP;
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
   double ****RCGP;
   /** centMolecVec is a vector used for centering the molecule. This array
    * is used for producing POV files. <b>Warning: In the current version,
    * after the POV has been recorded, the coordinates of the critical
    * points will no longer correspond to the coordinate system of the
    * wavefunction.</b>  */
   double centMolecVec[3];
   /** The array RGP is an auxiliary array, used for manipulating the order
    * of the points in a bond path searching procedure.  */
   double **RGP;
   string *lblACP, /*!< An array to store the labels of the ACPs.  */\
      *lblBCP, /*!< An array to store the labels of the BCPs. */\
      *lblRCP, /*!< An array to store the labels of the RCPs.  */\
      *lblCCP; /*!< An array to store the labels of the CCPs.  */
/* ************************************************************************************ */
//                   Functions
/* ************************************************************************************ */
   /** The only allowed constructor. It requires a GaussWaveFunction and
    * a BondNetWork objects. <b>Warning: these objects must be
    * properly initialized before passing to the CritPtNetWork object.</b>  */
   CritPtNetWork(class GaussWaveFunction &uwf,class BondNetWork &ubn);
   ~CritPtNetWork();
/* ************************************************************************************ */
   /** Self descriptive.  */
   void SetMaxIterationsACP(int ii) {maxItACP=ii;}
   /** Self descriptive.  */
   void SetMaxIterationsBCP(int ii) {maxItBCP=ii;}
   /** Self descriptive.  */
   void SetMaxIterationsRCP(int ii) {maxItRCP=ii;}
   /** Self descriptive.  */
   void SetMaxIterationsCCP(int ii) {maxItCCP=ii;}
   /** Self descriptive  */
   void SetMaxGradPathNPts(int ii) {maxGradPathNPts=ii;}
   /** Self descriptive  */
   void SetStepSizeBGP(double ss) {stepSizeBGP=ss;}
   /** The main public function for searching all critical points.
    * Configuration, such as requesting extended search should be
    * done before calling this function.
    * ft is used for requesting the type of critical points.
    * In current implementation, only the electron density
    * (DENS) and some of the LOL (LOLD) CPs are found.
    */
   void SetCriticalPoints(ScalarFieldType ft);
   /** It allocates the arrays to work with ACPs. It is public because for non
    * regular applications, only selected parts of the whole critical point search
    * are needed. For instance, one would be interested in searching ONLY ACPs.  */
   void SetupACPs(ScalarFieldType ft);
   /** This function searchs the ACPs. For regular calculations, use SetCriticalPoints()  */
   void SetACPs(ScalarFieldType ft);
   /** Basically, it allocates the arrays for working with BCPs. It needs to know
    * ACPs.  */
   void SetupBCPs(ScalarFieldType ft);
   /** This function search the BCPs. For regular calculations, use SetCriticalPoints().  */
   void SetBCPs(ScalarFieldType ft);
   /** Displays the coordinates of the critical points of type 'cpt'.
    * @param cpt: This char parameter is used to request the critical
    * point type. It can take the values: 'a', 'b', 'r', and 'c'.
    */
   void DisplayXCPCoords(char cpt);
   /** Prints to the std::cout the coordinates of all CPs found.  */
   void DisplayAllCPCoords(void);
   /** Self descriptive.  */
   void DisplayACPCoords(void) {DisplayXCPCoords('a');}
   /** Self descriptive.  */
   void DisplayBCPCoords(void) {DisplayXCPCoords('b');}
   /** Self descriptive.  */
   void DisplayRCPCoords(void) {DisplayXCPCoords('r');}
   /** Self descriptive.  */
   void DisplayCCPCoords(void) {DisplayXCPCoords('c');}
   /** Prints to the std::cout the coordinates used for seeding around a point.
    * Here seeding means to set starting points for a cp search. The letters
    * IHV are the acronynm of <b>I</b>cosa<b>H</b>edron <b>V</b>ertices.
    */
   void DisplayIHVCoords(void);
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
   void SeekRhoACPsAroundAPoint(double const (&oo)[3],double const ddxx,\
         string const &blbl,int nvrt=-1);
   void SeekRhoBCPWithExtraACP(int acppos,double maxrad);
   /** Same functionality as CritPtNetWork::SeekRhoACPsAroundAPoint, but for BCPs.  */
   void SeekRhoBCPsAroundAPoint(double const (&oo)[3],double const ddxx,\
         string const &blbl,int nvrt=-1);
   /** Same functionality as CritPtNetWork::SeekRhoACPsAroundAPoint, but for RCPs.  */
   void SeekRhoRCPsAroundAPoint(double const (&oo)[3],double const ddxx,\
         string const &blbl,int nvrt=-1);
   /** Same functionality as CritPtNetWork::SeekRhoACPsAroundAPoint, but for CCPs.  */
   void SeekRhoCCPsAroundAPoint(double const (&oo)[3],double const ddxx,\
         string const &blbl,int nvrt=-1);
   /** Same functionality as CritPtNetWork::SeekRhoACPsAroundAPoint, but for LOL ACPs.  */
   void SeekLOLACPsAroundAPoint(double const (&oo)[3],double const ddxx,\
         string const &blbl,int nvrt=-1);
   void SeekLOLACPsOnASphere(int atIdx,int nDivR,int nDivT,int nDivP,double radmin,\
         double radmax);
   void SeekLOLBCPWithExtraACP(int acppos,double maxrad);
   void CustomSearchTwoACPs(int acpIdx1,int acpIdx2);
   void CustomSearchCPs(double (&xs)[3]);
   void ExtendedSearchCPs();
   bool ReadFromFile(string inname);
   void DisplayStatus(bool lngdesc = false);
   /** Basically, it allocates the necessary arrays for working with bond paths.
    * It requires to know ACPs.  */
   void SetupBondPaths(void);
   /** The main function to compute bond paths. It requires to know BPCs (and ACPs).
    * This function calls SetupBondPaths internally. Within regular calculations
    * this function should be used.  */
   void SetBondPaths(void);
   void SetRingPaths(void);
   void SetCagePaths(void);
   bool SeekSingleRhoBCP(int ata,int atb,double (&x)[3]);
   int FindSingleRhoBondGradientPathRK5(int at1,int at2,double hstep,\
         int dima,double** (&arbgp),double (&ro)[3]);
   int FindSingleRhoRingGradientPathRK5(int rcpIdx,int bcpIdxInRRGP,\
         double hstep,int dima,double** (&arrgp));
   int FindSingleRhoCageGradientPathRK5(int ccpIdx,int rcpIdxInRCGP,\
         double hstep,int dima,double** (&arrgp));
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
   bool WalkGradientPathRK5ToEndPoint(double (&xi)[3],double (&x1)[3],\
         double (&xe)[3],double (&xm)[3],double &dm,double hstep,int dima,\
         double** (&arrgp),int &npia,double maxlen,bool uphilldir);
   void ForceBCPConnectivity(int bcpIdx,int acpIdx1,int acpIdx2);
   void CorrectRCPConnectivity(void);
   void RemoveFromConRCP(const int rcpIdx,const int pos2rem);
   void AddToConRCP(const int rcpIdx,const int bcpIdx);
   void AddToConCCP(const int ccpIdx,const int rcpIdx);
   void FindTwoClosestAtomsToBCP(int bcpIdx,int &at1Idx,int&at2Idx);
   /* ************************************************************************** */
   /* ************************************************************************** */
   bool MakePOVFile(string pnam,class POVRayConfiguration &pvp,int campos);
   void DrawNuclei(bool dn) {drawNuc=dn;}
   void DrawBonds(bool db) {drawBnd=db;}
   void DrawBondGradPaths(bool dbg) {drawBGPs=dbg;}
   void DrawRingGradPaths(bool dbg) {drawRGPs=dbg;}
   void TubeStyleBGP(bool stl) {tubeBGPStyle=stl;}
   void SetExtendedSearch(bool ss) {mkextsearch=ss;}
   void WriteCPProps(string &ofnam,string &wfnam);
   void PrintAllFieldProperties(double x,double y,double z);
   void WriteAllFieldProperties(double x,double y,double z,ofstream &ofil);
   bool IKnowACPs(void) {return iknowacps;}
   bool IKnowBCPs(void) {return iknowbcps;}
   bool IKnowRCPs(void) {return iknowrcps;}
   bool IKnowCCPs(void) {return iknowccps;}
   bool IKnowBGPs(void) {return iknowbgps;}
   bool IKnowRGPs(void) {return iknowrgps;}
   bool IKnowCGPs(void) {return iknowcgps;}
   ScalarFieldType MyCPType(void) {return mycptype;}
   double GetMaxBondDist() {return maxBondDist;}
   double GetMaxDistBetwBCPAndBCP() {return maxBCPACPDist;}
   int GetNofRingPathsOfRCP(int rcpIdx);
   int GetNofCagePathsOfCCP(int ccpIdx);
   int GetTotalNofRingPaths(void);
   int GetTotalNofCagePaths(void);
/* ************************************************************************************ */
protected:
/* ************************************************************************************ */
   class GaussWaveFunction *wf;
   class BondNetWork *bn;
   int dACP,dBCP,dRCP,dCCP;
   int maxItACP,maxItBCP,maxItRCP,maxItCCP;
   int maxGradPathNPts;
   int normalbcp;
   bool iknowacps,iknowbcps,iknowrcps,iknowccps, iknowallcps;
   bool iknowbgps,iknowrgps,iknowcgps,iknowallgps;
   bool drawNuc,drawBnd,drawBGPs,drawRGPs,drawCGPs;
   bool tubeBGPStyle;
   bool mkextsearch;
   double stepSizeACP,stepSizeBCP,stepSizeRCP,stepSizeCCP;
   double stepSizeBGP;
   ScalarFieldType mycptype;
   double maxBondDist; /*!< The maximum distance between two ACPs related by a BCP  */
   double maxBCPACPDist; /*!< The maximum distance between a BCP and associated ACPs  */
   static const int nIHV=16; //It is actually the vertices of an icosahedron plus the origin
   // (0,0,0)
   static double V0,V5,V8,IHV[nIHV][3];
   void Init();
   /** The constructor without arguments is not public. This will enforce the use
    * of the constructor for assigning the pointers of GaussWaveFunction and 
    * BondNetWork.  */
   CritPtNetWork();
   void RemoveRedundInLabel(string &lbl);
   string GetFirstChunkOfLabel(string &lbl);
   inline double ComputeMagnitudeV3(double (&v)[3])
          {return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);}
   bool ImNew(double (&x)[3],int dim,double ** (&arr),size_t &pos);
   void GetACPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig);
   void GetBCPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig);
   void GetRCPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig);
   void GetCCPStep(double (&g)[3],double (&hess)[3][3],double (&hh)[3],int &sig);
   int ComputeSignature(double (&ev)[3]);
   int ComputeSignature(double (&hh)[3][3]);
   void SeekRhoACP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig);
   void SeekRhoBCP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig);
   void SeekRhoRCP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig);
   void SeekRhoCCP(double (&x)[3],double &rho2ret,double (&g)[3],int &sig);
   void SeekLOLACP(double (&x)[3],double &ll,double (&g)[3],int &sig);
   void SeekLOLBCP(double (&x)[3],double &ll,double (&g)[3],int &sig);
   void SeekLOLRCP(double (&x)[3],double &ll,double (&g)[3],int &sig);
   void SeekLOLCCP(double (&x)[3],double &ll,double (&g)[3],int &sig);
   bool SetRhoACPs(void);
   bool SetRhoBCPs(void);
   bool SetRhoRCPs(void);
   bool SetRhoCCPs(void);
   bool SetLOLACPs(void);
   bool SetLOLBCPs(void);
   bool SetLOLRCPs(void);
   bool SetLOLCCPs(void);
   bool AddRhoACP(double (&x)[3],int sig,string &lbl);
   bool AddRhoBCP(double (&x)[3],int sig,string &lbl,int &pos);
   bool AddRhoRCP(double (&x)[3],int sig,string &lbl,int &pos);
   bool AddRhoCCP(double (&x)[3],int sig,string &lbl,int &pos);
   void FindTwoClosestAtoms(double (&xo)[3],int &idx1st,int &idx2nd);
   void FindTwoClosestACPs(double (&xo)[3],int &idx1st,int &idx2nd);
   void AddBCP2ConRCP(const int rcpIdx,const int bcpIdx);
   void AddRCP2ConCCP(const int ccpIdx,const int rcpIdx);
   void InvertOrderBGPPoints(int dim);
   void InvertOrderBGPPoints(int dim,double** (&arr));
   void GetNextPointInGradientPathRK5UpHill(double (&xn)[3],\
         double &stepsize,double &mgg);
   void GetNextPointInGradientPathRK5DownHill(double (&xn)[3],\
         double &stepsize,double &mgg);
   void CenterMolecule(void);
   void PutNuclei(ofstream &pof);
   void PutBonds(ofstream &pof);
   void FindMaxBondDist();
   void CopyRGP2Array(double** (&thearr),int nn);
/* ************************************************************************************ */
};
/* ************************************************************************************ */

#endif  /* _CRITPTNETWORK_H_ */

