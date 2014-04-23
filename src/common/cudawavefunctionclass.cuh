/*
    wavefunctionclass.h
 
   This class is for easy storaging/handling of a gaussian wave function. 
   The main class, namely gaussWaveFunc, is a class that storages all the information
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
      string *atLbl     --> The atom labels (for the wfn, the atom labels will be created as in the following examples)
                              O    1    (CENTRE  1)  ----->   O1
                              C    10   (CENTRE 10)  ----->   C10
      int *primType     --> The primitive type (Type Assignments) 
      int *primCent     --> The primitive center (Centre assignments)
      int *myPN         --> The number of primitives associated with each nuclear center.
      real *R           --> The coordinates of the nuclei. It is a 1-dimensional array. It is this way to facilitate
                              the implementation of cuda code. Later, a function getR(i,j) will be provided to access
                              easily the j-th coordinate of the i-th nucleus. (Beware the 0 index in c arrays)
      real *atCharge    --> The atomic charge
      real *primExp     --> The primitive exponent
      real *MOCoeff     --> The Molecular Orbital Coefficients. Also a 1D array. The function getCoeff(i,j) is provided
                              to get easy access to the j-th coefficient of the i-th Molecular Orbital. 
                              (Beware the 0 index in c arrays)
      real *occN        --> The occupation number of each molecular orbital.
      real *MOEner      --> The energy of each Molecular Orbital.
      real *cab         --> The matrix whose coefficients are $cab[a][b]=\sum_{i=1}^{nMOr}occN[i]*c_{ai}*c_{bi}$
      real *chi         --> Auxiliar array used in some functions
      real *gx,*gy,*gz  --> Auxiliar arrays used in some functions
      real *hxx,*hyy
           *hzz,*hxy
           *hxz,*hyz    --> Auxiliar arrays used in some functions
      real totener      --> The total energy of the system.
      real virial       --> The virial ratio (-V/T)
      bool imldd        --> true if the wave function have been loaded (from a wf? file), false otherwise.
 
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

#ifndef _SOLWAVEFUNCTIONCLASS_H_
#define _SOLWAVEFUNCTIONCLASS_H_

#ifndef _HAVE_DEF_REAL_TYPE_
#define _HAVE_DEF_REAL_TYPE_ 
typedef double real;
//typedef float real;
#endif

//#include <cuda.h>

#define MAXPRIMTYPEDEFINED 20
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

#ifndef SIGNF(a)
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
#include <math.h>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

extern "C" void setupWaveFunctionInGPU(int npr, \
                                       real** host_R,real** dev_R, \
                                       int** host_a,int** dev_a, \
                                                   real** dev_e, \
                                                   real** dev_c, \
                                       real** dev_chi,real** dev_rho, \
                                       real** dev_aux);
extern "C" void cleanWaveFunctionInGPU(real** host_R,real** dev_R, \
                                       int ** host_a,int ** dev_a, \
                                                     real** dev_e, \
                                                     real** dev_c, \
                                       real** dev_chi,real** dev_rho, \
                                       real** dev_aux);
extern "C" void copyDataToGPU(int npr,
                              real** host_R,real** dev_R, \
                              int ** host_a,int ** dev_a, \
                              real** host_e,real** dev_e, \
                              real** host_c,real** dev_c);
extern "C" void copyDataFromGPU(int npr,
                                real**dev_R,real**host_R, \
                                int **dev_a,int **host_a, \
                                real**dev_e,real**host_e, \
                                real**dev_c,real**host_c);
extern "C" void dosomethinginGPU(int npr, \
                                 real **dev_R,int **dev_a, \
                                 real **dev_e,real**dev_c);
extern "C" real evalDensityInGPU(real x,real y,real z,
                                 int npr,real **dev_R,int **dev_a,real **dev_e,real **dev_c, \
                                 real **host_chi,real **dev_chi,real **dev_rho);
extern "C" real evalDensityTileInGPU(real x,real y,real z,
                                     int npr,real **dev_R,int **dev_a,real **dev_e,real **dev_c, \
                                     real **host_chi,real **dev_chi,real **dev_rho,real **dev_aux);
__global__ void kernelR(int rdim,real *dev_R);
__global__ void kernela(int rdim,int *dev_a);
__global__ void kernele(int rdim,real *dev_e);
__global__ void kernelc(int rdim,real *dev_c);

__global__ void krnlCalcChi(int rdim,real *dev_R,int *dev_a,real *dev_e,real *dev_chi, \
                            real x,real y,real z);
__global__ void krnlCalcRho(int rdim,int offset,real *dev_c,real *dev_chi,real *dev_rho);


__device__ double atomicAdd(double* address, double val);
class gaussWaveFunc
{
public:
   //**********************************************************************************************
   gaussWaveFunc(); //Default constructor
   ~gaussWaveFunc(); //Destructor
   //**********************************************************************************************
   string *title,orbDesc; /* title */
   int nTit,nNuc,nMOr,nPri;
   string *atLbl;
   int *primType, *primCent,*myPN;
   real *R, *atCharge, *primExp, *MOCoeff, *occN, *MOEner;
   real *cab,*chi,*gx,*gy,*gz,*hxx,*hyy,*hzz,*hxy,*hxz,*hyz;
   real totener,virial;
   bool imldd;
   //*************************************************************************************************
   real *h_R,*d_R,*d_e,*d_c;
   int *h_a,*d_a;
   real *d_chi,*d_rho,*d_aux;
   //*************************************************************************************************
   bool readFromFileWFN(string inname);
   /*
      This function will load all the values of the wave function (title, orbDesc, etc.) from a 
      file, which name is inname. As the name of the function suggest, the file must be *.wfn The function
      will automatically allocate the corresponding arrays. And since a destructor is given, you
      one does not need to deallocate the arrays of the wave function.
    */
   //**********************************************************************************************
   bool readFromFileWFX(string inname);
   /*
      This function is essentially the same as readFromFileWFN, but using a *wfx file. In the future, it is
      expected that both functions differ from each other, since the files wfx are/will be capable of
      containing a larger amount of information.
    */
   //**********************************************************************************************
   bool readFromFile(string inname);
   /*
      This function just look for the extension of the inname, if it is wfn(wfx), then calls
      readFromFileWFN(readFromFileWFX)
    */
   //**********************************************************************************************
   bool allocAuxArrays(void);
   /*
      This function allocates memory space for the auxiliar arrays the gaussWaveFunction object
      uses for calculating numerical properties (rho, grad(rho), hess(rho), etc.).
    */
   //**********************************************************************************************
   bool setupGPU(void);
   //**********************************************************************************************
   void countPrimsPerCenter(void);
   /*
    This function counts the number of primitives associated with each one
    of the nuclear centers.
    */
   //**********************************************************************************************
   void calcCab(void);
   /*
      This function calculates the values of the matrix coefficients array $C_{\dot{A}\dot{B}}$
      (see notes ******* for notation details.)
    */
   //**********************************************************************************************
   bool testSupport(void);
   /*
      This function returns true if the *.wfn or *.wfx is supported. At the current version,
      only gaussian wave functions are handled.
    */
   //**********************************************************************************************
   real evalDensity(real x,real y,real z);
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
private:
   static int prTy[MAXPRIMTYPEDEFINED*3];
   //**********************************************************************************************
};
//**********************************************************************************************
//**********************************************************************************************
#endif//_SOLWAVEFUNCTIONCLASS_H_

