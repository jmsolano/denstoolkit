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
    gausswavefunction.h
 
   This class is for easy storaging/handling of a gaussian wave function. 
   The main class, namely GaussWaveFunction, is a class that storages all the information
   of a Gauss-Type wave function. It may contain also other common information contained in 
   *.wfn and *.wfx files, such as some title(s), some keywords/descriptor, atom labels, etc. 
   Below you can find the complete list of variables. For all arrays, you should remember
   that in c, arrays starts with index 0, rather than the standard 1 in fortran.
 
      int nTit          --> The number of titles (1 in wfn, arbitrary in wfx)
      string *title     --> The title(s) used in the *.wfn(*.wfx) file.
      string orbDesc    --> The orbital description (GAUSSIAN in wfn, GTO in wfx)
      int nNuc          --> The number of nuclei
      int nMOr          --> The number of Molecular Orbitals
      int nPri          --> The number of Primitives
      int EDFPri        --> The number of EDF primitives
      string *atLbl     --> The atom labels (for the wfn, the atom labels will be created as in the
                            following examples)
                              O    1    (CENTRE  1)  ----->   O1
                              C    10   (CENTRE 10)  ----->   C10
      int *primType     --> The primitive type (Type Assignments) 
      int *primCent     --> The primitive center (Centre assignments)
      int *myPN         --> The number of primitives associated with each nuclear center.
      solreal *R           --> The coordinates of the nuclei. It is a 1-dimensional array. It is 
                               this way to facilitate
                               the implementation of cuda code. Later, a function getR(i,j) will be 
                               provided to access easily the j-th coordinate of the i-th nucleus.
                               (Beware the 0 index in c arrays)
      solreal *atCharge    --> The atomic charge
      solreal *primExp     --> The primitive exponent
      solreal *MOCoeff     --> The Molecular Orbital Coefficients. Also a 1D array. The function getCoeff(i,j) is provided
                              to get easy access to the j-th coefficient of the i-th Molecular Orbital. 
                              (Beware the 0 index in c arrays)
      solreal *occN        --> The occupation number of each molecular orbital.
      solreal *MOEner      --> The energy of each Molecular Orbital.
      solreal *cab         --> The matrix whose coefficients are $cab[a][b]=\sum_{i=1}^{nMOr}occN[i]*c_{ai}*c_{bi}$
      solreal *chi         --> Auxiliar array used in some functions
      solreal *gx,*gy,*gz  --> Auxiliar arrays used in some functions
      solreal *hxx,*hyy
           *hzz,*hxy
           *hxz,*hyz    --> Auxiliar arrays used in some functions
      solreal totener      --> The total energy of the system.
      solreal virial       --> The virial ratio (-V/T)
      bool imldd        --> true if the wave function have been loaded (from a wf? file), false otherwise.
      bool ihaveEDF;    --> true if the wave function includes EDF information.
                               (this only works with wfx)
 
 The description of the functions is given below, right after every single function. 
 
 ------------------------
 
 Juan Manuel Solano Altamirano
 Adscription at the moment this project is initiated:
 Department of Chemistry, University of Guelph,
 Guelph, Ontario, Canada.
 e-mail: jmsolanoalt@gmail.com
 
 ------------------------
 
 This code is free code; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software 
 Foundation, Inc., 59 Temple Place - Suite 330, 
 Boston, MA  02111-1307, USA.
 
 WWW:  http://www.gnu.org/copyleft/gpl.html
 
 ----------------------
 
 */

#ifndef _GAUSSWAVEFUNCTION_H_
#define _GAUSSWAVEFUNCTION_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#ifndef PARALLELISEDTK
#define PARALLELISEDTK 0
#endif

#define MAXPRIMTYPEDEFINED 56
#define MAXSTEPSIZEBCPSEARCH 0.4
#define MAXSTEPSIZERCPSEARCH 0.35
#define MAXSTEPSIZECCPSEARCH 0.3
#define EPSGRADMAG 1.00000e-14

#ifndef MAXITERATIONBCPSEARCH
#define MAXITERATIONBCPSEARCH 80
#endif
#ifndef MAXITERATIONRCPSEARCH
#define MAXITERATIONRCPSEARCH 100
#endif
#define MAXITERATIONCCPSEARCH 240
#define MAXNUMBEROFPRIMITIVESFORMEMALLOC 10000

#ifndef SIGNF
#define SIGNF(a) ((a)>=0?(1):(-1))
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
#include <cmath>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;
#include <complex>
using std::complex;
#if PARALLELISEDTK
#include <omp.h>
#endif

class GaussWaveFunction {
public:
   /* *********************************************************************************** */
   GaussWaveFunction(); //Default constructor
   ~GaussWaveFunction(); //Destructor
   /* *********************************************************************************** */
   string *title,orbDesc; /* title */
   int nTit,nNuc,nMOr,nPri,EDFPri,totPri,coreElec;
   string *atLbl;
   int *primType, *primCent,*myPN;
   solreal *R, *atCharge, *primExp, *MOCoeff, *occN, *MOEner,*EDFCoeff;
   solreal *cab,*chi,*gx,*gy,*gz,*hxx,*hyy,*hzz,*hxy,*hxz,*hyz;
   solreal totener,virial;
   bool imldd,ihaveEDF;
   /* *********************************************************************************** */
   /**
      This function returns the value of the cart-th Cartesian coordinate of the nucnum-th nucleus.
      (0 for x, 1 for y and 1 for z coordinates.)
    */
   solreal getR(const int nucnum,const int cart);
   /* *********************************************************************************** */
   /** This function returns three angular exponents of the type pt (primitive type). */
   void getAng(int pt,int (&tt)[3]);
   /* *********************************************************************************** */
   /**
      This function returns the value of the primitive coefficient of the orbn-th orbital
      and primn-th primitive. As any function in c, the indices run from 0 to some number.
      You should be careful about this.
    */
   solreal getCoef(const int orbn,const int primn);
   /* *********************************************************************************** */
   /** 
      This function returns the value of the density (\f$\rho\f$) at the point \f$\vec{r}=(x,y,z)\f$ using
      a non optimized algorithm. It remains in this class for having a reference for the
      subsequent optimization of the function evalDensity. Unless you change evalDensity, use the 
      latter to obtain the density (it is the faster function).
    */
   solreal evalDensityArray(solreal x,solreal y, solreal z);
   /* *********************************************************************************** */
   /** This function calculates and prints to screen all the programmed field properties.
    */
   void displayAllFieldProperties(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /** This function calculates and writes to ofil all the programmed field properties.
    */
   void writeAllFieldProperties(solreal x,solreal y,solreal z,ofstream &ofil);
   /* *********************************************************************************** */
   /**
      This function returns the value of the density (\f$\rho\f$) at the point \f$\vec{r}=(x,y,z)\f$ using
      the most optimized algorithm.
    */
   solreal evalDensity(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   solreal evalOptimizedScalar(solreal px,solreal py,solreal pz);
   /* *********************************************************************************** */
   /**
      This function returns true if all the Molecular Orbital Occupation Numers have the same
      value. It only takes one occupation number to be different in order to this
      function returning false.
    */
   bool sameMolOrbOccNums(void);
   /* ************************************************************************************** */
   /** This function returns the value of the density (\f$\rho\f$) and the gradient of it 
      (\f$\nabla\rho=dx\hat{\imath}+dy\hat{\jmath}+dz\hat{k}\f$) at the point 
      \f$\vec{r}=(x,y,z)\f$ using the most optimized algorithm.
    */
   void evalRhoGradRho(solreal x, solreal y, solreal z,\
         solreal &rho, solreal &dx, solreal &dy, solreal &dz);
   /* ************************************************************************************** */
   void evalOptimizedVectorScalar(solreal x, solreal y, solreal z,\
         solreal &rho, solreal &dx, solreal &dy, solreal &dz);
   /* ************************************************************************************** */
   /**
      This function is the same as evalRhoGradRho, but using an array for the gradient instead of
      individual components of the gradient. i.e.
      evalRhoGradRho(x,y,z,rho,g[3])=evalRhoGradRho(x,y,z,rho,g[0],g[1],g[2])
    */
   void evalRhoGradRho(solreal x, solreal y, solreal z,solreal &rho, solreal (&grd)[3]);
   /* ************************************************************************************** */
   /**
      This function will load all the values of the wave function (title, orbDesc, etc.) from a 
      file, which name is inname. As the name of the function suggest, the file must be *.wfn The function
      will automatically allocate the corresponding arrays. And since a destructor is given, you
      one does not need to deallocate the arrays of the wave function.
    */
   bool readFromFileWFN(string inname);
   /* *********************************************************************************** */
   /**
      This function is essentially the same as readFromFileWFN, but using a *wfx file. 
      In the future, it is expected that both functions differ from each other, 
      since the files wfx are/will be capable of containing a larger amount of information.
    */
   bool readFromFileWFX(string inname);
   /* *********************************************************************************** */
   /**
      This function just look for the extension of the inname, if it is wfn(wfx), then calls
      readFromFileWFN(readFromFileWFX)
    */
   bool readFromFile(string inname);
   /* *********************************************************************************** */
   /** As the name suggests, the function performs a series of test to verify the sanity
    * (or suitability) of the loaded wavefunction. In the first implementation,
    * for instance, if the primitive type of any primitive center is higher than
    * 20, a warning message will be printed about the incompleteness of the implemented
    * fiedls.  */
   bool sanityChecks(void);
   /* *********************************************************************************** */
   /** This function allocates memory space for the auxiliar arrays the GaussWaveFunctiontion object
      uses for calculating numerical properties (rho, grad(rho), hess(rho), etc.).
    */
   bool allocAuxArrays(void);
   /* *********************************************************************************** */
   /** This function counts the number of primitives associated with each one
      of the nuclear centers.
    */
   void countPrimsPerCenter(void);
   /* *********************************************************************************** */
   /** This function calculates the values of the matrix coefficients array \f$C_{\dot{A}\dot{B}}\f$
      (see notes ******* for notation details.)
    */
   void calcCab(void);
   /* *********************************************************************************** */
   /** This function will write the wave function into a wfx file which name is outname. */
   bool writeToFileWFX(string outname);
   /* *********************************************************************************** */
   /** This function returns true if the *.wfn or *.wfx is supported. At the current version,
      only gaussian wave functions are handled.
    */
   bool testSupport(void);
      /* *********************************************************************************** */
   /** This funtion returns the value of the angular part of the primitive, i.e., it returns the
      value of \f$x^{a_1}y^{a_2}z^{a^3}\f$, where a1,a2,a3 are the values of the angular exponents.
      The values of a_i are coded in the value of pty (the type of primitive). For correct results,
      x must be the difference between the field point and the primitive center.
    */
   solreal evalAngACases(int &pty, solreal x, solreal y, solreal z);
   /* *********************************************************************************** */
   /** Let the field point be \f$\vec{\xi}\f$, and a primitive \f$\phi_{\dot{A}}\f$, with center 
      \f$R_{\dot{A}}\f$, type pty, and primitive exponent alp. The function
      evalDkAngCases returns the value of
      (\f$\nabla\phi_{\dot{A}}(x,y,z))\exp(-2\alpha r^2)\f$, where \f$x=\xi_x-R^x_{\dot{A}}\f$..., and
      \f$r^2=x^2+y^2+z^2\f$. The individual components are anx, any, anz, respectively
    */
   void evalDkAngCases(int &pty,solreal alp,solreal x, solreal y, solreal z,\
         solreal &anx, solreal &any, solreal &anz);
   /* *********************************************************************************** */
   /**
    Let the field point be \f$\vec{\xi}\f$, and a primitive \f$\phi_{\dot{A}}\f$, with center 
    \f$R_{\dot{A}}\f$, type pty, and primitive exponent alp. The function
    evalDkDlAngCases returns the value of
    (\f$\partial_k\partial_l\phi_{\dot{A}}(x,y,z))/exp(-2 alp r^2)\f$, where x=\xi_x-R^x_{\dot{A}}..., and
    \f$r^2=x^2+y^2+z^2\f$. The individual components are axx,ayy,azz,axy,axz,ayz, respectively.
    */
   void evalDkDlAngCases(int &pty,solreal alp,solreal x,solreal y,solreal z,\
         solreal &axx,solreal &ayy,solreal &azz,solreal &axy,solreal &axz,solreal &ayz);
   /* *********************************************************************************** */
   /** Let the field point be \f$\vec{\xi}\f$, and a primitive \f$\phi_{\dot{A}}\f$, with center 
    \f$R_{\dot{A}}\f$, type pty, and primitive exponent alp. The function
    evalLapAngCases returns the value of
    (\f$\nabla^2\phi_{\dot{A}}(x,y,z))/exp(-2 alp rr)\f$, where x=\xi_x-R^x_{\dot{A}}..., and
    \f$rr=x^2+y^2+z^2\f$.
    */
   solreal evalLapAngCases(int &pty,solreal alp,solreal x,solreal y,solreal z,solreal rr);
   /* *********************************************************************************** */
   /** This function evaluates the six independent components of the Hessian of \f$\rho\f$ 
      (dxx,dyy,dzz,dxy,dxz,dyz) at the point \f$(x,y,z)\f$
    */
   void evalHessian(solreal x, solreal y, solreal z,solreal &dxx, solreal &dyy, solreal &dzz,\
         solreal &dxy, solreal &dxz, solreal &dyz);
   /* *********************************************************************************** */
   /** This function evaluates the Hessian of \f$\rho\f$ at the point \f$(x,y,z)\f$ and 
    * store them in the array h.
    */
   void evalHessian(solreal x, solreal y, solreal z,solreal (&h)[3][3]);
   /* *********************************************************************************** */
   /** This function evaluates the gradient and Hessian of \f$\rho\f$ at the point 
    \f$(x,y,z)\f$ and store them in the arrays g and h.
    */
   void evalHessian(solreal x, solreal y, solreal z,solreal &dens,solreal (&g)[3],solreal (&h)[3][3]);
   /* *********************************************************************************** */
   void evalOptimizedScalVecHess(solreal x, solreal y, solreal z,solreal &dens,solreal (&g)[3],solreal (&h)[3][3]);
   /* *********************************************************************************** */
   /**
      This function returns the value of the Laplacian of the density
      (\f$\nabla^2\rho\f$) at the point \f$(x,y,z)\f$
    */
   solreal evalLapRho(solreal x, solreal y, solreal z);
   /* *********************************************************************************** */
   /**
      This function returns the value of the Electron Localization Function, EFL 
      (\f$\eta(x,y,z)=\frac{1}{1+[D(x,y,z)/D_h(x,y,z)]^2}\f$, where \f$D(x,y,z)=\frac{1}{2}
        \sum_i|\nabla\phi_i(x,y,z)|^2-\frac{1}{8}\frac{|\nabla\rho(x,y,z)|^2}{\rho(x,y,z)}\f$,
      and \f$D_h(x,y,z)=(3/10)(3\phi^2)^{2/3}\rho(x,y,z)^{5/3}\f$) at the
      point (x,y,z).
    */
   solreal evalELF(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /**
      This function returns the value of the Shannon entropy density (\f$-\rho\ln\rho\f$) at the
      point (x,y,z).
    */
   solreal evalShannonEntropy(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /**
    This function returns the value of the Shannon entropy density in the momentum space 
    (\f$-\pi\ln\pi\f$) at the point (px,py,pz).
    */
   solreal evalMomentumShannonEntropy(solreal px,solreal py,solreal pz);
   /* *********************************************************************************** */
   /**
    This function returns the value of the Magnitude of the Density Gradient 
    (\f$|\nabla\rho|\f$) at the point (x,y,z).
    */
   solreal evalMagGradRho(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /**
    This function returns the value of the Localized Orbital Locator, LOL, 
    (\f$\gamma(x,y,z)=\frac{\tau(x,y,z)}{1+\tau(x,y,z)}\f$, where
     \tau=2D_h(x,y,z)/(\sum_i|\nabla\chi_i|^2)) at the field point (x,y,z)
    */
   solreal evalLOL(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /**
    This function returns the value of the Kinetic Energy Density G, defined through
    \f$G(\vec{x})=\frac{1}{2}\sum_{\dot{A}\sum_{\dot{B}}}C_{\dot{A}\dot{B}}\nabla\phi_{\dot{A}}
     \cdot\nabla\phi_{\dot{B}}\f$
    */
   solreal evalKineticEnergyG(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /**
    This function returns the value of
    \f$\sum_{\dot{A}\sum_{\dot{B}}}C_{\dot{A}\dot{B}}\nabla\phi_{\dot{A}}
     \cdot\nabla\phi_{\dot{B}}\f$, and rho. It severs for evaluation of
     indices such as the 'reduced density gradient' or 'region of
     slow electrons'.
    */
   void evalNabPhi2(solreal const x,solreal const y,solreal const z,\
         solreal &rho2ret,solreal &twoG);
   /* *********************************************************************************** */
   /**
    This function returns the value of the Kinetic Energy Density K, defined through
    \f$K(\vec{x})=-\frac{1}{4}\sum_{\dot{A}\sum_{\dot{B}}}C_{\dot{A}\dot{B}}
    (\phi_{\dot{A}}\nabla^2\phi_{\dot{B}}+\phi_{\dot{B}}\nabla^2\phi_{\dot{A}})\f$
    */
   solreal evalKineticEnergyK(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /**
      This function seeks for a Bond Critical Point. The integers ii and jj are used to set
      \f$\vec{x}_0=\frac{1}{2}(\vec{R}_i+\vec{R}_j)\f$. The final values of the search are 
      stored in rx,ry,rz; in addition the values of the gradient at the point rx,ry,rz are 
      saved in gx,gy,gz. If the maximum number of iterations is reached and the critical
      was not found, then rx,ry,rz are the last values obtained from the search (and
      the gradient at such a point).
      <b>Important: For searching critical points, please use the class critPtNetWork instead of
      this function. This function will be deprecated within the next few revisions.</b>
    */
   void seekBondCP(int ii,int jj,solreal &rx,solreal &ry,solreal &rz,solreal &gx,solreal &gy,solreal &gz);
   /* *********************************************************************************** */
   /*
      This function uses the vector \f$\vec{x}\f$ as the original position, then it calculates 
      the step hh using the eigen-vector following algorithm. This algorithm is described with 
      detail for  this particular problem in 
         Chem. Phys. Lett. 228 (1994) 160--164, "A robust algorithm to locate automatically
            all types of critical points in the charge density and its Laplacian".
      More details (and ) can be found in references [9-11] of the above article.
      This particular function aims to locate Bond Critical Points.
   */
   void getBondCPStep(solreal (&x)[3],solreal (&hh)[3],solreal (&gg)[3]);
   /* *********************************************************************************** */
   /*
    This function seeks for a Ring Critical Point. The vector \f$\vec{x}_0=(r1,r2,r3)\f$ is used as
    the starting point for the search. The final values of the search are 
    stored in rx,ry,rz; in addition the values of the gradient at the point rx,ry,rz are 
    saved in gx,gy,gz. If the maximum number of iterations is reached and the critical
    was not found, then rx,ry,rz are the last values obtained from the search (and
    the gradient at such a point).
    */
   void seekRingCP(solreal &r1,solreal &r2,solreal &r3,solreal &gx,solreal &gy,solreal &gz);
   /* *********************************************************************************** */
   /*
    This function uses the vector \f$\vec{x}\f$ as the original position, then it calculates 
    the step hh using the eigen-vector following algorithm. This algorithm is described with 
    detail for  this particular problem in 
    Chem. Phys. Lett. 228 (1994) 160--164, "A robust algorithm to locate automatically
    all types of critical points in the charge density and its Laplacian".
    More details (and ) can be found in references [9-11] of the above article.
    This particular function aims to locate Ring Critical Points.
    */
   void getRingCPStep(solreal (&x)[3],solreal (&hh)[3],solreal (&g)[3]);
   /* *********************************************************************************** */
   /*
    This function seeks for a Cage Critical Point. The vector \f$\vec{x}_0=(r1,r2,r3)\f$ is used as
    the starting point for the search. The final values of the search are 
    stored in rx,ry,rz; in addition the values of the gradient at the point r1,r2,r3 are 
    saved in gx,gy,gz. If the maximum number of iterations is reached and the critical
    was not found, then rx,ry,rz are the last values obtained from the search (and
    the gradient at such a point).
    */
   void seekCageCP(solreal &r1,solreal &r2,solreal &r3,solreal &gx,solreal &gy,solreal &gz);
   /* *********************************************************************************** */
   /*
    This function uses the vector \f$\vec{x}\f$ as the original position, then it calculates 
    the step hh using the eigen-vector following algorithm. This algorithm is described with 
    detail for  this particular problem in 
    Chem. Phys. Lett. 228 (1994) 160--164, "A robust algorithm to locate automatically
    all types of critical points in the charge density and its Laplacian".
    More details (and ) can be found in references [9-11] of the above article.
    This particular function aims to locate Cage Critical Points.
    */
   void getCageCPStep(solreal (&x)[3],solreal (&hh)[3],solreal (&g)[3]);
   /* *********************************************************************************** */
   void evald1SingCartA(int &ang,solreal &t,solreal x,solreal &d0,solreal &d1);
   /* *********************************************************************************** */
   void evald2SingCartA(int &ang,solreal &t,solreal x,\
                      solreal &d0,solreal &d1,solreal &d2);
   /* *********************************************************************************** */
   void evald3SingCartA(int &ang,solreal &t,solreal &f,solreal &x,solreal &x2,
                      solreal &d0,solreal &d1,solreal &d2,solreal &d3);
   /* *********************************************************************************** */
   /** This function returns the third derivatives of the angular factors for evaluating
    * the third derivatives of the primitives. Since these derivatives are usually
    * needed when computing the derivatives of the electron density, the function
    * also evaluates the zero-th (angular factor), first, second, and third derivative
    * factors of the primitives. The final derivatives of the primitives are
    * obtained by multiplying the components of di times \f$\exp(-r^2\alpha)\f$.  */
   void evald3Ang(int (&a)[3],solreal &alp,solreal (&x)[3],solreal (&x2)[3],
                  solreal (&d0)[3],solreal (&d1)[3],solreal (&d2)[3],solreal (&d3)[3]);
   /* *********************************************************************************** */
   void evald4SingCartA(int &ang,solreal &t,solreal &f,solreal &x,solreal &x2,
                        solreal &d0,solreal &d1,solreal &d2,solreal &d3,solreal &d4);
   /* *********************************************************************************** */
   /** This function returns the fourth derivatives of the angular factors for evaluating
    * the fourth derivatives of the primitives. Since these derivatives are usually
    * needed when computing the derivatives of the electron density, the function
    * also evaluates the zero-th (angular factor), first, second, third, and fourth derivative
    * factors of the primitives. The final derivatives of the primitives are
    * obtained by multiplying the components of di times \f$\exp(-r^2\alpha)\f$.  */
   void evald4Ang(int (&a)[3],solreal &alp,solreal (&x)[3],solreal (&x2)[3],
                  solreal (&d0)[3],solreal (&d1)[3],solreal (&d2)[3],solreal (&d3)[3],
                  solreal (&d4)[3]);
   /* *********************************************************************************** */
   /* This function evaluates \f$\nabla^2\rho(x,y,z)\f$. It is implemented in order to
    * test evalDkDlAngCases(...)  */
   //solreal evalLapRhoUsingd2(solreal x,solreal y,solreal z);
   /* *********************************************************************************** */
   /* void evalDkDlAngCases(int &pty,solreal alp,solreal x,solreal y,solreal z,
    solreal &axx,solreal &ayy,solreal &azz,solreal &axy,solreal &axz,solreal &ayz); */
   void evalDiDjDkChi(int &pty,solreal &alp,solreal x,solreal y,solreal z,\
                      solreal (&dlm)[3][3],solreal (&dijk)[3][3][3]);
   /* *********************************************************************************** */
   /** Evaluates the Hessian of LOL (dxx,dyy,dzz,dxy,dxz,dyz). 
    * On the fly, it evaluates the electron density (dens),
    * the kinetic energy G (keG, LOL (lol), and the gradient of LOL (dx,dy,dz).  */
   void evalHessLOL(solreal x, solreal y, solreal z, solreal &dens,solreal &keG, solreal &lol,
         solreal &ddx, solreal &ddy, solreal &ddz,
         solreal &dxx, solreal &dyy, solreal &dzz,
         solreal &dxy, solreal &dxz, solreal &dyz);
   /* *********************************************************************************** */
   /** Evaluates the Hessian of LOL. It is an overloaded function.  */
   void evalHessLOL(solreal (&x)[3],solreal &lol,solreal (&glol)[3],solreal (&hlol)[3][3]);
   /* ************************************************************************************ */
   /**
     This function evaluates the "angular" part of the Fourier transform of every primitive.
       ang is the angular exponent for the component, ooa is 1/alpha, osra=1/sqrt(a),
       Rx is the component of the primitive center, px is the component of the momentum, and
       px2=px*px.
    */
   void evalFTASingCartA(int &ang,solreal &a,solreal &ooa,solreal &osra,\
         solreal &px,solreal &px2,solreal &Rx,\
         solreal &RePhi,solreal &ImPhi);
   /* ************************************************************************************ */
   void evalFTAng(int (&a)[3],solreal &alp,solreal &ooalp,solreal (&p)[3],solreal (&p2)[3],\
                  solreal (&Rx)[3],complex<solreal> &pang);
   /* ************************************************************************************ */
   void evalFTChi(int &pty,solreal &alp,solreal (&Rx)[3],solreal px,solreal py,solreal pz,
                      complex<solreal> &phi);
   /* ************************************************************************************ */
   /** This function computes the Fourier Transform of the electron density. <em>i.e.</em>
    * the momentum-space electron density.  */
   solreal evalFTDensity(solreal px,solreal py,solreal pz);
   /* ************************************************************************************ */
   /** This function computes the Kinetic energy density K in momentum space.  */
   solreal evalFTKineticEnergy(solreal px,solreal py,solreal pz);
   /* ************************************************************************************ */
   /** This function evaluates the Density Matrix of Order 1 at the points
    * \f$(x,y,z)\f$ and \f$(xp,yp,zp)\f$  */
   solreal evalDensityMatrix1(solreal x,solreal y,solreal z,solreal xp,solreal yp,solreal zp);
   /* ************************************************************************************ */
   /** This function evaluates the gradients of the Density Matrix of order 1, with
    * respect to the primed and non-primed variables.  */
   void evalGradDensityMatrix1(solreal x,solreal y,solreal z,\
         solreal xp,solreal yp,solreal zp,\
         solreal &gamm,solreal (&gg)[3],solreal (&gp)[3]);
   /* ************************************************************************************ */
   /** This function computes the Hessian of the Density Matrix of Order 1 (second
    * derivatives), with respect to the primed, non-primed and combined cases.  */
   void evalHessDensityMatrix1(solreal (&xx)[3],solreal (&xxp)[3],\
         solreal &gamm,solreal (&gg)[3],solreal (&gp)[3],\
         solreal (&hh)[3][3],solreal (&hph)[3][3],solreal (&hp)[3][3]);
   /* ************************************************************************************ */
   /** This function evaluates the magnitude of the grandient of LOL.  */
   solreal evalMagGradLOL(solreal x,solreal y,solreal z);
   /* ************************************************************************************ */
   /** This function returns the Molecular Electrostatic Potential (MEP) at the
    * point (x,y,z).  */
   solreal evalMolElecPot(solreal x,solreal y,solreal z);
   /* ************************************************************************************ */
   void evalHermiteCoefs(int (&aia)[3],int (&aib)[3],solreal &alpab,
                         solreal (&ra)[3],solreal (&rb)[3],
                         solreal (&rp)[3],
                         int (&maxl)[3],solreal (&Eijl)[3][7]);//tested on Dec 28th, 2013
   /* ************************************************************************************ */
   void evalRlmnIntegs(const int (&lmn)[3],const solreal &alpp,const solreal (&cp)[3],
                       solreal (&Rlmn)[7][7][7]);
   /* ************************************************************************************ */
   solreal evalVAB(solreal (&xx)[3],int (&aa)[3],int (&ab)[3],solreal &alpa,solreal &alpb,
                          solreal (&ra)[3],solreal (&rb)[3]);
   /* ************************************************************************************ */
   solreal evalOverlapIntegralAB(int (&aa)[3],int (&ab)[3],solreal &alpa,solreal &alpb,
                         solreal (&ra)[3],solreal (&rb)[3]);
   /* ************************************************************************************ */
   /** This function returns the integral of the electron density. It uses the analytical
    * properties of the Gauss-type orbital basis.  */
   solreal integralRho(void);
   /* ************************************************************************************ */
   /** This function returns the sum of the total number of protons of the molecule.  */
   solreal totalNuclearCharge(void);
   /* ************************************************************************************ */
   /** This function evaluates the Localized Electrons Detector vector (led) at the
    * point x.  */
   void evalLED(solreal const (&x)[3],solreal (&led)[3]);
   /* ************************************************************************************ */
   /** This function returns the magnitude of the vector LED.  */
   solreal evalMagLED(solreal x,solreal y,solreal z);
   /* ************************************************************************************ */
   /** This function returns the Reduced Density Gradient at the point (x,y,z)  */
   solreal evalReducedDensityGradient(solreal x,solreal y,solreal z);
   /* ************************************************************************************ */
   /** This function computes the Region of Slow Electron index at the point
    * (x,y,z).  */
   solreal evalRoSE(solreal x,solreal y,solreal z);
   /* ************************************************************************************ */
   /** This function is left to the final user for implementing its own custom 
    * scalar field.  */
   solreal evalCustomScalarField(solreal x,solreal y,solreal z);
   /* ************************************************************************************ */
   /** This function is left to the final user for implementing its own custom 
    * vector field.  */
   void evalCustomVectorField(solreal x,solreal y,solreal z,solreal (&v)[3]);
   /* ************************************************************************************ */
   void useScalarCustomField(bool ucf) {usescustfld=ucf;}
   /* ************************************************************************************ */
   void useVectorCustomField(bool ucf) {usevcustfld=ucf;}
   /* ************************************************************************************ */
   /** This function returns the Potential Energy Density at the point (x,y,z)
    * This field is taken from: "Hydrogen bond strengths revealed by topological
    * analyses of experimentally observed electron densities",
    * E. Espinosa, E. Molins, C. Lecomte, Chemical Physics Letters, 285 (1998), 170-173. */
   solreal evalVirialPotentialEnergyDensity(solreal x, solreal y, solreal z);
   /* ************************************************************************************ */
   /**The funtion evalNCIs returning values of Reduced Density Gradient 
    * applying NCI conditions.    
    */
   solreal evalNCIs(solreal x, solreal y, solreal z,solreal cutoff=2.0);
   /* ************************************************************************************ */
   /**The funtion evalNCILambda returning a value of the second eingevalues of the Hessian  
    * applying NCI conditions.   
    */
   solreal evalNCILambda(solreal x, solreal y, solreal z);
   /* ************************************************************************************ */
protected:
   static int prTy[MAXPRIMTYPEDEFINED*3];
   /* ************************************************************************************ */
   bool usescustfld,usevcustfld;
   int maxPrimType;
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
   /* ************************************************************************************ */
};
/* *********************************************************************************** */
/* *********************************************************************************** */
#endif//_GAUSSWAVEFUNCTION_H_

