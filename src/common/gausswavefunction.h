/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
      double *R           --> The coordinates of the nuclei. It is a 1-dimensional array. It is 
                               this way to facilitate
                               the implementation of cuda code. Later, a function GetR(i,j) will be 
                               provided to access easily the j-th coordinate of the i-th nucleus.
                               (Beware the 0 index in c arrays)
      double *atCharge    --> The atomic charge
      double *primExp     --> The primitive exponent
      double *MOCoeff     --> The Molecular Orbital Coefficients. Also a 1D array. The function getCoeff(i,j) is provided
                              to get easy access to the j-th coefficient of the i-th Molecular Orbital. 
                              (Beware the 0 index in c arrays)
      double *occN        --> The occupation number of each molecular orbital.
      double *MOEner      --> The energy of each Molecular Orbital.
      int *MOsptp         --> The spin type of each Molecular Orbital.
      double *cab         --> The matrix whose coefficients are $cab[a][b]=\sum_{i=1}^{nMOr}occN[i]*c_{ai}*c_{bi}$
      double *chi         --> Auxiliar array used in some functions
      double *gx,*gy,*gz  --> Auxiliar arrays used in some functions
      double *hxx,*hyy
           *hzz,*hxy
           *hxz,*hyz    --> Auxiliar arrays used in some functions
      double totener      --> The total energy of the system.
      double virial       --> The virial ratio (-V/T)
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

#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>
using std::string;
#include <complex>
using std::complex;

#if PARALLELISEDTK
#include <omp.h>
#endif

#ifndef NCIRHOMAX
#define NCIRHOMAX 6.5e-2
#endif

#ifndef NCIRHOMIN
#define NCIRHOMIN 5.0e-4
#endif

#ifndef NCISMAX
#define NCISMAX 2.0e0
#endif


/* *********************************************************************************** */
class GaussWaveFunction {
/* *********************************************************************************** */
public:
/* *********************************************************************************** */
   GaussWaveFunction(); //Default constructor
   ~GaussWaveFunction(); //Destructor
   /* *********************************************************************************** */
   string *title,orbDesc; /* title */
   int nTit,nNuc,nMOr,nPri,EDFPri,totPri,coreElec;
   string *atLbl;
   int *primType, *primCent,*myPN,*MOsptp;
   double *R, *atCharge, *primExp, *MOCoeff, *occN, *MOEner,*EDFCoeff;
   double *cab,*cabA,*cabB,*chi,*gx,*gy,*gz,*hxx,*hyy,*hzz,*hxy,*hxz,*hyz;
   double *prefMEP;
   double totener,virial;
   double nciRhoMin,nciRhoMax,nciSMax;
   bool imldd,ihaveEDF,ihaveSingleSpinOrbs,ihaveCABSingleSpin;
   /* *********************************************************************************** */
   /**
      This function returns the value of the cart-th Cartesian coordinate of the nucnum-th nucleus.
      (0 for x, 1 for y and 1 for z coordinates.)
    */
   double GetR(const int nucnum,const int cart);
   /** Sets the s cutoff for NCI  */
   void SetNCISMax(double sm) {nciSMax=sm;}
   /** Sets the rho lower cutoff for NCI  */
   void SetNCIRhoMin(double rr) {nciRhoMin=rr;}
   /** Sets the rho upper cutoff for NCI  */
   void SetNCIRhoMax(double rr) {nciRhoMax=rr;}
   /** This function returns three angular exponents of the type pt (primitive type). */
   void GetAng(int pt,int (&tt)[3]);
   /**
      This function returns the value of the primitive coefficient of the orbn-th orbital
      and primn-th primitive. As any function in c, the indices run from 0 to some number.
      You should be careful about this.
    */
   double GetCoef(const int orbn,const int primn);
   /** 
      This function returns the value of the density (\f$\rho\f$) at the point \f$\vec{r}=(x,y,z)\f$ using
      a non optimized algorithm. It remains in this class for having a reference for the
      subsequent optimization of the function EvalDensity. Unless you change EvalDensity, use the 
      latter to obtain the density (it is the faster function).
    */
   double EvalDensityArray(double x,double y, double z);
   /** This function calculates and prints to screen all the programmed field properties.
    */
   void DisplayAllFieldProperties(double x,double y,double z);
   /* *********************************************************************************** */
   /** This function calculates and writes to ofil all the programmed field properties.
    */
   void WriteAllFieldProperties(double x,double y,double z,ofstream &ofil);
   /**
      This function returns the value of the density (\f$\rho\f$) at the point \f$\vec{r}=(x,y,z)\f$ using
      the most optimized algorithm.
    */
   double EvalDensity(double x,double y,double z);
   double EvalSpinDensity(double x,double y,double z);
   double EvalOptimizedScalar(double px,double py,double pz);
   /* *********************************************************************************** */
   /**
      This function returns true if all the Molecular Orbital Occupation Numers have the same
      value. It only takes one occupation number to be different in order to this
      function returning false.
    */
   bool SameMolOrbOccNums(void);
   /** This function returns the value of the density (\f$\rho\f$) and the gradient of it 
      (\f$\nabla\rho=dx\hat{\imath}+dy\hat{\jmath}+dz\hat{k}\f$) at the point 
      \f$\vec{r}=(x,y,z)\f$ using the most optimized algorithm. */
   void EvalRhoGradRho(double x, double y, double z,\
         double &rho, double &dx, double &dy, double &dz);
   void EvalOptimizedVectorScalar(double x, double y, double z,\
         double &rho, double &dx, double &dy, double &dz);
   /** This function is the same as EvalRhoGradRho, but using an array for the gradient instead of
      individual components of the gradient. i.e.
      EvalRhoGradRho(x,y,z,rho,g[3])=EvalRhoGradRho(x,y,z,rho,g[0],g[1],g[2]) */
   void EvalRhoGradRho(double x, double y, double z,double &rho, double (&grd)[3]);
   /**
      This function will load all the values of the wave function (title, orbDesc, etc.) from a 
      file, which name is inname. As the name of the function suggest, the file must be *.wfn The function
      will automatically allocate the corresponding arrays. And since a destructor is given, you
      one does not need to deallocate the arrays of the wave function.
    */
   bool ReadFromFileWFN(string inname);
   /**
      This function is essentially the same as ReadFromFileWFN, but using a *wfx file. 
      So far, two additional capabilities are implemented:
      a) If the wfx has EDF obitals, the wavefunction will function accordingly.
    */
   bool ReadFromFileWFX(string inname);
   /**
      This function just look for the extension of the inname, if it is wfn(wfx), then calls
      ReadFromFileWFN(ReadFromFileWFX)
      */
   bool ReadFromFile(string inname);
   /** As the name suggests, the function performs a series of test to verify the sanity
      (or suitability) of the loaded wavefunction. In the first implementation,
      for instance, if the primitive type of any primitive center is higher than
      20, a warning message will be printed about the incompleteness of the implemented
      fields. 
   */
   bool SanityChecks(void);
   /** This function allocates memory space for the auxiliar arrays the GaussWaveFunctiontion object
      uses for calculating numerical properties (rho, grad(rho), hess(rho), etc.). */
   bool AllocAuxArrays(void);
   bool AllocAuxMEPArray(void);
   /** This function counts the number of primitives associated with each one
      of the nuclear centers. */
   void CountPrimsPerCenter(void);
   /** This function calculates the values of the matrix coefficients array \f$C_{\dot{A}\dot{B}}\f$
     (see notes ******* for notation details.) */
   void CalcCab(void);
   /** This function calculates the values of the single-spn-density matrices
     \f$C_{\dot{A}\dot{B}}^{\alpha}\f$ and \f$C_{\dot{A}\dot{B}}^{\beta}\f$
     (see notes ******* for notation details.) */
   void CalcCabAAndCabB(void);
   /** This function will write the wave function into a wfx file which name is outname. */
   bool WriteToFileWFX(string outname);
   /** This function returns true if the *.wfn or *.wfx is supported. At the current version,
      only gaussian wave functions are handled. */
   bool TestSupport(void);
   /** This funtion returns the value of the angular part of the primitive, i.e., it returns the
      value of \f$x^{a_1}y^{a_2}z^{a^3}\f$, where a1,a2,a3 are the values of the angular exponents.
      The values of a_i are coded in the value of pty (the type of primitive). For correct results,
      x must be the difference between the field point and the primitive center.*/
   double EvalAngACases(int &pty, double x, double y, double z);
   /** Let the field point be \f$\vec{\xi}\f$, and a primitive \f$\phi_{\dot{A}}\f$, with center 
      \f$R_{\dot{A}}\f$, type pty, and primitive exponent alp. The function
      EvalDkAngCases returns the value of
      (\f$\nabla\phi_{\dot{A}}(x,y,z))\exp(-2\alpha r^2)\f$, where \f$x=\xi_x-R^x_{\dot{A}}\f$..., and
      \f$r^2=x^2+y^2+z^2\f$. The individual components are anx, any, anz, respectively*/
   void EvalDkAngCases(int &pty,double alp,double x, double y, double z,\
         double &anx, double &any, double &anz);
   /** Let the field point be \f$\vec{\xi}\f$, and a primitive \f$\phi_{\dot{A}}\f$, with center 
     \f$R_{\dot{A}}\f$, type pty, and primitive exponent alp. The function
     EvalDkDlAngCases returns the value of
     (\f$\partial_k\partial_l\phi_{\dot{A}}(x,y,z))/exp(-2 alp r^2)\f$, where x=\xi_x-R^x_{\dot{A}}..., and
     \f$r^2=x^2+y^2+z^2\f$. The individual components are axx,ayy,azz,axy,axz,ayz, respectively. */
   void EvalDkDlAngCases(int &pty,double alp,double x,double y,double z,\
         double &axx,double &ayy,double &azz,double &axy,double &axz,double &ayz);
   /** Let the field point be \f$\vec{\xi}\f$, and a primitive \f$\phi_{\dot{A}}\f$, with center 
     \f$R_{\dot{A}}\f$, type pty, and primitive exponent alp. The function
     EvalLapAngCases returns the value of
     (\f$\nabla^2\phi_{\dot{A}}(x,y,z))/exp(-2 alp rr)\f$, where x=\xi_x-R^x_{\dot{A}}..., and
     \f$rr=x^2+y^2+z^2\f$.*/
   double EvalLapAngCases(int &pty,double alp,double x,double y,double z,double rr);
   /** This function evaluates the six independent components of the Hessian of \f$\rho\f$ 
      (dxx,dyy,dzz,dxy,dxz,dyz) at the point \f$(x,y,z)\f$ */
   void EvalHessian(double x, double y, double z,double &dxx, double &dyy, double &dzz,\
         double &dxy, double &dxz, double &dyz);
   /** This function evaluates the Hessian of \f$\rho\f$ at the point \f$(x,y,z)\f$ and 
     store them in the array h. */
   void EvalHessian(double x, double y, double z,double (&h)[3][3]);
   /** This function evaluates the gradient and Hessian of \f$\rho\f$ at the point 
     \f$(x,y,z)\f$ and store them in the arrays g and h. */
   void EvalHessian(double x, double y, double z,double &dens,double (&g)[3],double (&h)[3][3]);
   void EvalOptimizedScalVecHess(double x, double y, double z,double &dens,double (&g)[3],double (&h)[3][3]);
   /**
      This function returns the value of the Laplacian of the density
      (\f$\nabla^2\rho\f$) at the point \f$(x,y,z)\f$*/
   double EvalLapRho(double x, double y, double z);
   /**
      This function returns the value of the Electron Localization Function, EFL 
      (\f$\eta(x,y,z)=\frac{1}{1+[D(x,y,z)/D_h(x,y,z)]^2}\f$, where \f$D(x,y,z)=\frac{1}{2}
        \sum_i|\nabla\phi_i(x,y,z)|^2-\frac{1}{8}\frac{|\nabla\rho(x,y,z)|^2}{\rho(x,y,z)}\f$,
      and \f$D_h(x,y,z)=(3/10)(3\phi^2)^{2/3}\rho(x,y,z)^{5/3}\f$) at the
      point (x,y,z).*/
   double EvalELF(double x,double y,double z);
   /**
      This function returns the value of the Shannon entropy density (\f$-\rho\ln\rho\f$) at the
      point (x,y,z).*/
   double EvalShannonEntropy(double x,double y,double z);
   /**
    This function returns the value of the Shannon entropy density in the momentum space 
    (\f$-\pi\ln\pi\f$) at the point (px,py,pz). */
   double EvalMomentumShannonEntropy(double px,double py,double pz);
   /**
    This function returns the value of the Magnitude of the Density Gradient 
    (\f$|\nabla\rho|\f$) at the point (x,y,z). */
   double EvalMagGradRho(double x,double y,double z);
   /**
    This function returns the value of the Localized Orbital Locator, LOL, 
    (\f$\gamma(x,y,z)=\frac{\tau(x,y,z)}{1+\tau(x,y,z)}\f$, where
     \tau=2D_h(x,y,z)/(\sum_i|\nabla\chi_i|^2)) at the field point (x,y,z) */
   double EvalLOL(double x,double y,double z);
   /**
    This function returns the value of the Kinetic Energy Density G, defined through
    \f$G(\vec{x})=\frac{1}{2}\sum_{\dot{A}\sum_{\dot{B}}}C_{\dot{A}\dot{B}}\nabla\phi_{\dot{A}}
     \cdot\nabla\phi_{\dot{B}}\f$ */
   double EvalKineticEnergyG(double x,double y,double z);
   /** This function returns the value of
    \f$\sum_{\dot{A}\sum_{\dot{B}}}C_{\dot{A}\dot{B}}\nabla\phi_{\dot{A}}
     \cdot\nabla\phi_{\dot{B}}\f$, and rho. It severs for evaluation of
     indices such as the 'reduced density gradient' or 'region of
     slow electrons'. */
   void EvalNabPhi2(double const x,double const y,double const z,\
         double &rho2ret,double &twoG);
   /** This function returns the value of the Kinetic Energy Density K, defined through
     \f$K(\vec{x})=-\frac{1}{4}\sum_{\dot{A}\sum_{\dot{B}}}C_{\dot{A}\dot{B}}
     (\phi_{\dot{A}}\nabla^2\phi_{\dot{B}}+\phi_{\dot{B}}\nabla^2\phi_{\dot{A}})\f$ */
   double EvalKineticEnergyK(double x,double y,double z);
   /** Returns the ellipticity at the point (x,y,z). The ellipticity is defined as
      \f$\varepsilon(\vec x)=\frac{\lambda_1}{\lambda_2}-1\f$. Here \f$\lambda_1\ (\lambda_2)\f$
      is the first (second) Hessian eigenvalue. See 
      Chem. Phys. Lett., 143 (1988) 450 - 458 for an example of how to use
      ellipticyt profiles. */
   double EvalEllipticity(double x,double y,double z);
   /** This function seeks for a Bond Critical Point. The integers ii and jj are used to set
      \f$\vec{x}_0=\frac{1}{2}(\vec{R}_i+\vec{R}_j)\f$. The final values of the search are 
      stored in rx,ry,rz; in addition the values of the gradient at the point rx,ry,rz are 
      saved in gx,gy,gz. If the maximum number of iterations is reached and the critical
      was not found, then rx,ry,rz are the last values obtained from the search (and
      the gradient at such a point).
      <b>Important: For searching critical points, please use the class CritPtNetWork instead of
      this function. This function will be deprecated within the next few revisions.</b> */
   void SeekBondCP(int ii,int jj,double &rx,double &ry,double &rz,double &gx,double &gy,double &gz);
   /* This function uses the vector \f$\vec{x}\f$ as the original position, then it calculates 
      the step hh using the eigen-vector following algorithm. This algorithm is described with 
      detail for  this particular problem in 
         Chem. Phys. Lett. 228 (1994) 160--164, "A robust algorithm to locate automatically
            all types of critical points in the charge density and its Laplacian".
      More details (and ) can be found in references [9-11] of the above article.
      This particular function aims to locate Bond Critical Points. */
   void GetBondCPStep(double (&x)[3],double (&hh)[3],double (&gg)[3]);
   /* This function seeks for a Ring Critical Point. The vector \f$\vec{x}_0=(r1,r2,r3)\f$ is used as
    the starting point for the search. The final values of the search are 
    stored in rx,ry,rz; in addition the values of the gradient at the point rx,ry,rz are 
    saved in gx,gy,gz. If the maximum number of iterations is reached and the critical
    was not found, then rx,ry,rz are the last values obtained from the search (and
    the gradient at such a point). */
   void SeekRingCP(double &r1,double &r2,double &r3,double &gx,double &gy,double &gz);
   /* This function uses the vector \f$\vec{x}\f$ as the original position, then it calculates 
    the step hh using the eigen-vector following algorithm. This algorithm is described with 
    detail for  this particular problem in 
    Chem. Phys. Lett. 228 (1994) 160--164, "A robust algorithm to locate automatically
    all types of critical points in the charge density and its Laplacian".
    More details (and ) can be found in references [9-11] of the above article.
    This particular function aims to locate Ring Critical Points. */
   void GetRingCPStep(double (&x)[3],double (&hh)[3],double (&g)[3]);
   /* This function seeks for a Cage Critical Point. The vector \f$\vec{x}_0=(r1,r2,r3)\f$ is used as
    the starting point for the search. The final values of the search are 
    stored in rx,ry,rz; in addition the values of the gradient at the point r1,r2,r3 are 
    saved in gx,gy,gz. If the maximum number of iterations is reached and the critical
    was not found, then rx,ry,rz are the last values obtained from the search (and
    the gradient at such a point). */
   void SeekCageCP(double &r1,double &r2,double &r3,double &gx,double &gy,double &gz);
   /* This function uses the vector \f$\vec{x}\f$ as the original position, then it calculates 
      the step hh using the eigen-vector following algorithm. This algorithm is described with 
      detail for  this particular problem in 
      Chem. Phys. Lett. 228 (1994) 160--164, "A robust algorithm to locate automatically
      all types of critical points in the charge density and its Laplacian".
      More details (and ) can be found in references [9-11] of the above article.
      This particular function aims to locate Cage Critical Points. */
   void GetCageCPStep(double (&x)[3],double (&hh)[3],double (&g)[3]);
   void Evald1SingCartA(int &ang,double &t,double x,double &d0,double &d1);
   void Evald2SingCartA(int &ang,double &t,double x,\
                      double &d0,double &d1,double &d2);
   void Evald3SingCartA(int &ang,double &t,double &f,double &x,double &x2,
                      double &d0,double &d1,double &d2,double &d3);
   /** This function returns the third derivatives of the angular factors for evaluating
      the third derivatives of the primitives. Since these derivatives are usually
      needed when computing the derivatives of the electron density, the function
      also evaluates the zero-th (angular factor), first, second, and third derivative
      factors of the primitives. The final derivatives of the primitives are
      obtained by multiplying the components of di times \f$\exp(-r^2\alpha)\f$.  */
   void Evald3Ang(int (&a)[3],double &alp,double (&x)[3],double (&x2)[3],
                  double (&d0)[3],double (&d1)[3],double (&d2)[3],double (&d3)[3]);
   void Evald4SingCartA(int &ang,double &t,double &f,double &x,double &x2,
                        double &d0,double &d1,double &d2,double &d3,double &d4);
   /** This function returns the fourth derivatives of the angular factors for evaluating
      the fourth derivatives of the primitives. Since these derivatives are usually
      needed when computing the derivatives of the electron density, the function
      also evaluates the zero-th (angular factor), first, second, third, and fourth derivative
      factors of the primitives. The final derivatives of the primitives are
      obtained by multiplying the components of di times \f$\exp(-r^2\alpha)\f$.  */
   void Evald4Ang(int (&a)[3],double &alp,double (&x)[3],double (&x2)[3],
                  double (&d0)[3],double (&d1)[3],double (&d2)[3],double (&d3)[3],
                  double (&d4)[3]);
   /* This function evaluates \f$\nabla^2\rho(x,y,z)\f$. It is implemented in order to
    test EvalDkDlAngCases(...)  */
    //double EvalLapRhoUsingd2(double x,double y,double z);
   /* void EvalDkDlAngCases(int &pty,double alp,double x,double y,double z,
    double &axx,double &ayy,double &azz,double &axy,double &axz,double &ayz); */
   void EvalDiDjDkChi(int &pty,double &alp,double x,double y,double z,\
                      double (&dlm)[3][3],double (&dijk)[3][3][3]);
   /** Evaluates the Hessian of LOL (dxx,dyy,dzz,dxy,dxz,dyz). 
     On the fly, it evaluates the electron density (dens),
     the kinetic energy G (keG, LOL (lol), and the gradient of LOL (dx,dy,dz).  */
   void EvalHessLOL(double x, double y, double z, double &dens,double &keG, double &lol,
         double &ddx, double &ddy, double &ddz,
         double &dxx, double &dyy, double &dzz,
         double &dxy, double &dxz, double &dyz);
   /** Evaluates the Hessian of LOL. It is an overloaded function.  */
   void EvalHessLOL(double (&x)[3],double &lol,double (&glol)[3],double (&hlol)[3][3]);
   /**
     This function evaluates the "angular" part of the Fourier transform of every primitive.
       ang is the angular exponent for the component, ooa is 1/alpha, osra=1/sqrt(a),
       Rx is the component of the primitive center, px is the component of the momentum, and
       px2=px*px. */
   void EvalFTASingCartA(int &ang,double &a,double &ooa,double &osra,\
         double &px,double &px2,double &Rx,\
         double &RePhi,double &ImPhi);
   void EvalFTAng(int (&a)[3],double &alp,double &ooalp,double (&p)[3],double (&p2)[3],\
                  double (&Rx)[3],complex<double> &pang);
   void EvalFTChi(int &pty,double &alp,double (&Rx)[3],double px,double py,double pz,
                      complex<double> &phi);
   /** This function computes the Fourier Transform of the electron density. <em>i.e.</em>
      the momentum-space electron density.  */
   double EvalFTDensity(double px,double py,double pz);
   /** This function computes the Kinetic energy density K in momentum space.  */
   double EvalFTKineticEnergy(double px,double py,double pz);
   /** This function evaluates the Density Matrix of Order 1 at the points
      \f$(x,y,z)\f$ and \f$(xp,yp,zp)\f$  */
   double EvalDensityMatrix1(double x,double y,double z,double xp,double yp,double zp);
   /** This function evaluates the Density Matrix of Order 1 in momentum space at
       the points \f$(x,y,z)\f$ and \f$(xp,yp,zp)\f$, using coefficients
       \f$c_{\dot A\dot B}\f$ as provided through cabs. singlespin tells the
       function whether the cabs coefficients are for single spin orbitals or
       alpha-beta orbitals.  */
   complex<double> EvalGeneralFTDensityMatrix1(double p1x,double p1y,double p1z,\
         double p2x,double p2y,double p2z,const bool singlespin,double *cabs);
   inline complex<double> EvalFTDensityMatrix1(double p1x,double p1y,double p1z,\
         double p2x,double p2y,double p2z) {
      return EvalGeneralFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z,false,cab); }
   inline complex<double> EvalFTDensityMatrix1Alpha(double p1x,double p1y,double p1z,\
         double p2x,double p2y,double p2z) {
      return EvalGeneralFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z,true,cabA); }
   inline complex<double> EvalFTDensityMatrix1Beta(double p1x,double p1y,double p1z,\
         double p2x,double p2y,double p2z) {
      return EvalGeneralFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z,true,cabB); }
   /** This function evaluates the gradients of the Density Matrix of order 1, with
      respect to the primed and non-primed variables.  */
   void EvalGradDensityMatrix1(double x,double y,double z,\
         double xp,double yp,double zp,\
         double &gamm,double (&gg)[3],double (&gp)[3]);
   /** This function computes the Hessian of the Density Matrix of Order 1 (second
     derivatives), with respect to the primed, non-primed and combined cases.
     xx is the non primed evaluation point, xxp is the primed evaluation point,
     gamm is Gamma, gg is GradGamma, gp is Grad'gamma, hh is \$f\partial_i\partial_jGamma\$f,
     hph is \$f\partial_i\{partial'}_jGamma\$f, and hp is
     \$f{\partial'}_i{\partial'}_jGamma\$f.
     Check DeMat1CriticalPointNetworkBP::AssignHessian6D for more details. */
   void EvalHessDensityMatrix1(double (&xx)[3],double (&xxp)[3],\
         double &gamm,double (&gg)[3],double (&gp)[3],\
         double (&hh)[3][3],double (&hph)[3][3],double (&hp)[3][3]);
   double EvalLapDensityMatrix1(double (&xx)[3],double (&xxp)[3]);
   /** This function evaluates the generalized Density Matrix of Order 1 at the points
     \f$(x,y,z)\f$ and \f$(xp,yp,zp)\f$. Here by generalized it is to be understood that
     the wavefunction may or not be open-shell. Therefore, the function requires
     the matrix \f$c_{\dot{A}\dot{B}}\f$ to be passed as an argument as well as
     a bool singlespin. If singlespin==true, then the correction stemming from
     pseudopotentials will be applied (for single-spin \f$\Gamma_1^{\sigma}\f$,
     the pseudo-potential is half the contribution for closed-shell systems).
     Notice that this function does not check whether the cabs pointer is valid,
     this check must be done at a far higher level
     (when loading the wavefunction). */
   double EvalGeneralDensityMatrix1(const double x,const double y,const double z,\
         const double xp,const double yp,const double zp,const bool singlespin,double *cabs);
   /** Returns the \f$\Gamma_1^{\alpha}(x,y,z,x',y',z')\f$  */
   inline double EvalDensityMatrix1Alpha(const double x,const double y,const double z,\
         const double xp,const double yp,const double zp) {
      return EvalGeneralDensityMatrix1(x,y,z,xp,yp,zp,true,cabA);
   }
   /** Returns the \f$\Gamma_1^{\beta}(x,y,z,x',y',z')\f$  */
   inline double EvalDensityMatrix1Beta(const double x,const double y,const double z,\
         const double xp,const double yp,const double zp) {
      return EvalGeneralDensityMatrix1(x,y,z,xp,yp,zp,true,cabB);
   }
   /** Returns the spinless pair density function,
     \f$\rho_2(\vector{r}_1,\vector{r}_2)\f$,
     evaluated at the points \f$\vector{r}_1=(x_1,y_1,z_1)\f$ and
     \f$\vector{r}_2=(x_2,y_2,z_2)\f$, of a closed-shell
     system.   */
   double EvalRho2ClosedShell(const double x1,const double y1,const double z1,\
         const double x2,const double y2,const double z2);
   /** Returns the spinless pair density function,
     \f$\rho_2(\vector{r}_1,\vector{r}_2)\f$,
     evaluated at the points \f$\vector{r}_1=(x_1,y_1,z_1)\f$ and
     \f$\vector{r}_2=(x_2,y_2,z_2)\f$, of an open-shell
     system.   */
   double EvalRho2OpenShell(const double x1,const double y1,const double z1,\
         const double x2,const double y2,const double z2);
   /** This function evaluates the magnitude of the grandient of LOL.  */
   double EvalMagGradLOL(double x,double y,double z);
   /** This function returns the Molecular Electrostatic Potential (MEP) at the
      point (x,y,z).  */
   double EvalMolElecPot(double x,double y,double z);
   void EvalHermiteCoefs(int (&aia)[3],int (&aib)[3],double &alpab,
                         double (&ra)[3],double (&rb)[3],
                         double (&rp)[3],
                         int (&maxl)[3],double (&Eijl)[3][7]);//tested on Dec 28th, 2013
   void EvalRlmnIntegs(const int (&lmn)[3],const double &alpp,const double (&cp)[3],
                       double (&Rlmn)[7][7][7]);
   double EvalVAB(double (&xx)[3],int (&aa)[3],int (&ab)[3],double &alpa,double &alpb,
                          double (&ra)[3],double (&rb)[3]);
   double EvalVABCore(double &s00,double (&xx)[3],int idxA,int idxB,int (&aa)[3],int (&ab)[3],double &alpa,double &alpb,
                          double (&ra)[3],double (&rb)[3]);
   double EvalOverlapIntegralAB(int (&aa)[3],int (&ab)[3],double &alpa,double &alpb,
                         double (&ra)[3],double (&rb)[3]);
   /** This function returns the integral of the electron density. It uses the analytical
      properties of the Gauss-type orbital basis.  */
   double IntegralRho(void);
   /** This function returns the sum of the total number of protons of the molecule.  */
   double TotalNuclearCharge(void);
   /** This function evaluates the Localized Electrons Detector vector (led) at the
      point x. For the definition, see H. J. Boh\'orquez, C. F. Matta, and R. J. Boyd,
      Int. Jour. Quantum Chem., 110 (2010) 2418--2425. The field is defined as
      \f$\tilde{\boldsymbol P}=\frac{\nabla\rho}{2\rho}\f$.  */
   void EvalLED(double const (&x)[3],double (&led)[3]);
   /** This function returns the magnitude of the vector LED.
      See the dtk-manual for more details.  */
   double EvalMagLED(double x,double y,double z);
   /** This function returns the Reduced Density Gradient at the point (x,y,z).
      See the dtk-manual for more details. */
   double EvalReducedDensityGradient(double x,double y,double z);
   void EvalGradReducedDensityGradient(double x,double y,double z,double (&gs)[3]);
   /** This function computes the Region of Slow Electron index at the point
    * (x,y,z). See the dtk-manual for more details. */
   double EvalRoSE(double x,double y,double z);
   /** This function is left to the final user for implementing its own custom 
      scalar field.  */
   double EvalCustomScalarField(double x,double y,double z);
   /** This function is left to the final user for implementing its own custom 
      vector field.  */
   void EvalCustomVectorField(double x,double y,double z,double (&v)[3]);
   void UseScalarCustomField(bool ucf) {usescustfld=ucf;}
   void UseVectorCustomField(bool ucf) {usevcustfld=ucf;}
   /** This function returns the Potential Energy Density at the point (x,y,z)
      This field is taken from: "Hydrogen bond strengths revealed by topological
      analyses of experimentally observed electron densities",
      E. Espinosa, E. Molins, C. Lecomte, Chemical Physics Letters, 285 (1998), 170-173. */
   double EvalVirialPotentialEnergyDensity(double x,double y,double z);
   /**The funtion EvalNCIs returns values of Reduced Density Gradient 
      applying NCI conditions.
      To set the cutoffs, use setNCISMax, setNCIRhoMin, and setNCIRhoMax
      See J. Contreras-Garcia, E. R. Johnson, S. Keinan, R. Chaudret, J-P. Piquemal,
      D. N. Beratan, and W. Yang. J. Chem. Theory Comput. 2011, 7, pp 625-632
      for more details */
   double EvalNCIs(double x,double y,double z);
   /**Returns rho*sigma (if nciRhoMin<rho<nciRhoMax), where sigma is the second Hessian eigenvalue sign
      To set the cutoffs (nciRhoMin and nciRhoMax)), use setNCISMax, setNCIRhoMin, and setNCIRhoMax
      See J. Contreras-Garcia, E. R. Johnson, S. Keinan, R. Chaudret, J-P. Piquemal,
      D. N. Beratan, and W. Yang. J. Chem. Theory Comput. 2011, 7, pp 625-632
      for more details */
   double EvalNCILambda(double x,double y,double z);
   double EvalDORI(double x,double y,double z);
   int GetPrTy(int idx) {return prTy[idx];}
/* ************************************************************************************ */
protected:
/* ************************************************************************************ */
   static int prTy[MAXPRIMTYPEDEFINED*3];
   static bool prntHermiCoeffWrn;
   int maxPrimType;
   bool usescustfld,usevcustfld;
/* ************************************************************************************ */
};
/* *********************************************************************************** */
#endif//_GAUSSWAVEFUNCTION_H_

