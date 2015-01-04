#ifndef _DEMAT1CRITPTNETWORK_H_
#define _DEMAT1CRITPTNETWORK_H_

#include <cstdlib>
#include <cstring>
using std::string;

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
 * the object DeMat1CriticalPointNetwork.
 * */
class DeMat1CriticalPointNetwork {
/* ************************************************************************************ */
public:
/* ************************************************************************************ */
   DeMat1CriticalPointNetwork(class gaussWaveFunc *usrwf,int at1,int at2);
   ~DeMat1CriticalPointNetwork();
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
/* ************************************************************************************ */
   void displaySCPsInfo(void);
/* ************************************************************************************ */
   void displayRCPsInfo(void);
/* ************************************************************************************ */
   void displayCPsInfo(void);
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
   DeMat1CriticalPointNetwork();
/* ************************************************************************************ */
   class gaussWaveFunc* wf;
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




#endif  /* _DEMAT1CRITPTNETWORK_H_ */

