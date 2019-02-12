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

#ifndef _DEMAT1CRITPTNETWORKSL_H_
#define _DEMAT1CRITPTNETWORKSL_H_

#include <cstdlib>
#include <cstring>
using std::string;
#include <fstream>
using std::ofstream;

#ifndef DEMAT1MAXITERATIONACPSEARCH
#define DEMAT1MAXITERATIONACPSEARCH 30
#endif

#ifndef DEMAT1MAXITERATIONSCPSEARCH
#define DEMAT1MAXITERATIONSCPSEARCH 30
#endif

#ifndef DEMAT1MAXITERATIONRCPSEARCH
#define DEMAT1MAXITERATIONRCPSEARCH 30
#endif

#ifndef DEFAULTNDIVSACPSEARCH
#define DEFAULTNDIVSACPSEARCH 8
#endif

#ifndef DEFAULTNDIVSSCPSEARCH
#define DEFAULTNDIVSSCPSEARCH 10
#endif

#ifndef DEFAULTNDIVSRCPSEARCH
#define DEFAULTNDIVSRCPSEARCH 12
#endif

/* ************************************************************************************ */
/** 
 * This class will find, store and manipulate the critical points of the
 * density matrix of order 1.
 * For this class, it is assumed that the density matrix is studied in a plane
 * whose coordinates are given by pairs (u,v) (we denote this plane by
 * 'the plane UV' or the 'UV plane'). Therefore, there are only 
 * three types of critical points, which we will call as follows.
 * Also, it is assumed that only two nuclei are contained within the
 * plane UV.
 *
 * Abbrev.:         Name                  (range,signature)
 *
 *   ACP:   Attractive Critical Point   (2,-2)
 *   SCP:   Saddle Critical Point       (2, 0)
 *   RCP:   Repulsive CP                (2,+2)
 *
 *   The name is just for keeping some compatibility with the
 *   naming of three dimensional critical points (used in
 *   Bader's QTAIM).
 * 
 * When this class is instantiated, two atom indices (i,j) must be given.
 * In the first implementation, this will define the direction
 * of the straight line in real space that joins the atoms i and j,
 * and whose parametrization is (u,v). This indices are assumed to
 * be from zero to the number of nuclei contained in the 
 * wave function object, which is also required for instantiating 
 * the object DeMat1CriticalPointNetworkSL.
 * */
class DeMat1CriticalPointNetworkSL {
/* ************************************************************************************ */
public:
/* ************************************************************************************ */
   DeMat1CriticalPointNetworkSL(class GaussWaveFunction *usrwf,int at1,int at2);
   ~DeMat1CriticalPointNetworkSL();
/* ************************************************************************************ */
   int nACP,nSCP,nRCP;
   int asACP,asSCP,asRCP;
/* ************************************************************************************ */
   solreal **RACP,**RSCP,**RRCP;
/* ************************************************************************************ */
   string *lblACP,*lblSCP,*lblRCP;
/* ************************************************************************************ */
   solreal x1[3],x2[3],x2mx1[3],x2pmx1p[3];
   solreal lenline,oolenline;
/* ************************************************************************************ */
   void getXCoordinatesFromUV(solreal uu,solreal vv,solreal (&x)[3],solreal (&xp)[3]);
/* ************************************************************************************ */
   void evalUVGrad(solreal uu,solreal vv,solreal &gamm, solreal (&uvg)[2]);
/* ************************************************************************************ */
   void evalUVHessian(solreal uu,solreal vv,solreal &gamm,\
         solreal (&uvg)[2],solreal (&uvh)[2][2]);
/* ************************************************************************************ */
   void getACPStep(solreal (&g)[2],solreal (&hess)[2][2],solreal (&hh)[2],int &sig);
/* ************************************************************************************ */
   void getSCPStep(solreal (&g)[2],solreal (&hess)[2][2],solreal (&hh)[2],int &sig);
/* ************************************************************************************ */
   void getRCPStep(solreal (&g)[2],solreal (&hess)[2][2],solreal (&hh)[2],int &sig);
/* ************************************************************************************ */
   void seekGammaACP(solreal (&x)[2],solreal &gamm2ret,solreal (&g)[2],\
         int &sig,int maxit=DEMAT1MAXITERATIONACPSEARCH);
/* ************************************************************************************ */
   void seekGammaSCP(solreal (&x)[2],solreal &gamm2ret,solreal (&g)[2],int &sig,\
         int maxit=DEMAT1MAXITERATIONSCPSEARCH);
/* ************************************************************************************ */
   void seekGammaRCP(solreal (&x)[2],solreal &gamm2ret,solreal (&g)[2], int &sig,\
         int maxit=DEMAT1MAXITERATIONRCPSEARCH);
/* ************************************************************************************ */
   void setGammaACP(int ndivs=DEFAULTNDIVSACPSEARCH);
/* ************************************************************************************ */
   void seekGammaACPsAroundAPoint(solreal (&oo)[2],solreal ddxx);
/* ************************************************************************************ */
   void seekSingleGammaACP(solreal (&xs)[2],solreal &gamm,solreal (&gg)[2],string &lbl);
/* ************************************************************************************ */
   void setGammaSCP(int ndivs=DEFAULTNDIVSSCPSEARCH);
/* ************************************************************************************ */
   void seekGammaSCPsAroundAPoint(solreal (&oo)[2],solreal ddxx);
/* ************************************************************************************ */
   void seekSingleGammaSCP(solreal (&xs)[2],solreal &gamm,solreal (&gg)[2],string &lbl);
/* ************************************************************************************ */
   void setGammaRCP(int ndivs=DEFAULTNDIVSRCPSEARCH);
/* ************************************************************************************ */
   void seekGammaRCPsAroundAPoint(solreal (&oo)[2],solreal ddxx);
/* ************************************************************************************ */
   void seekSingleGammaRCP(solreal (&xs)[2],solreal &gamm,solreal (&gg)[2],string &lbl);
/* ************************************************************************************ */
/* ************************************************************************************ */
   void addGammaACP(solreal (&x)[2],string lbl);
/* ************************************************************************************ */
   void addGammaSCP(solreal (&x)[2],string lbl);
/* ************************************************************************************ */
   void addGammaRCP(solreal (&x)[2],string lbl);
/* ************************************************************************************ */
   bool imNew(solreal (&x)[2],int dim,solreal ** (&arr),size_t &pos);
/* ************************************************************************************ */
   void displayACPsInfo(void);
   void displaySCPsInfo(void);
   void displayRCPsInfo(void);
   void displayCPsInfo(void);
/* ************************************************************************************ */
   void writeACPsInfo(ofstream &ofil);
   void writeSCPsInfo(ofstream &ofil);
   void writeRCPsInfo(ofstream &ofil);
   void writeCPsInfo(ofstream &ofil);
/* ************************************************************************************ */
   void setGammaCriticalPoints(void);
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
   static const int nPolyV=13; //It is actually the vertices of an icosahedron plus the origin 
                             // (0,0,0)
   static solreal PolyV[nPolyV][2];
/* ************************************************************************************ */
private:
/* ************************************************************************************ */
   int ata,atb;
/* ************************************************************************************ */
   DeMat1CriticalPointNetworkSL();
/* ************************************************************************************ */
   class GaussWaveFunction* wf;
/* ************************************************************************************ */
   void init(void);
/* ************************************************************************************ */
   void computePolygonVertices(void);
/* ************************************************************************************ */
   inline solreal getV2Norm(solreal (&vv)[2]) {return sqrt(vv[0]*vv[0]+vv[1]*vv[1]);}
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
};
/* ************************************************************************************ */




#endif  /* _DEMAT1CRITPTNETWORKSL_H_ */

