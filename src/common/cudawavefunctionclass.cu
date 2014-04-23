/*
 *  wavefunctionclass.cpp

 
 
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

#ifndef _SOLWAVEFUNCTIONCLASS_CPP_
#define _SOLWAVEFUNCTIONCLASS_CPP_

#include "cudawavefunctionclass.cuh"
#include "solmemhand.cpp"
#include "iofuncts-wfn.cpp"
#include "iofuncts-wfx.cpp"
#include "eig2-4.cpp"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef EPSFORELFVALUE
#define EPSFORELFVALUE (2.871e-05)
#endif

#ifndef EPSFORLOLVALUE
#define EPSFORLOLVALUE (2.871e-05)
#endif

#ifndef BASETHREADSPERBLOCK
#define BASETHREADSPERBLOCK 16
#endif

//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
gaussWaveFunc::gaussWaveFunc()
//*************************************************************************************************
{
   title=NULL;
   orbDesc=string("");
   nTit=0;
   nNuc=0;
   nMOr=0;
   nPri=0;
   atLbl=NULL;
   primType=NULL;
   primCent=NULL;
   myPN=NULL;
   R=NULL;
   atCharge=NULL;
   primExp=NULL;
   MOCoeff=NULL;
   occN=NULL;
   MOEner=NULL;
   cab=NULL;
   chi=NULL;
   gx=gy=gz=NULL;
   hxx=hyy=hzz=NULL;
   hxy=hxz=hyz=NULL;
   totener=0.00e0;
   virial=0.0e0;
   imldd=false;
   h_R=d_R=d_e=d_c=NULL;
   h_a=d_a=NULL;
   d_chi=d_rho=d_aux=NULL;
}
//*************************************************************************************************
int gaussWaveFunc::prTy[]={
   0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 0, 1,   2, 0, 0, 
   0, 2, 0,   0, 0, 2,   1, 1, 0,   1, 0, 1,   0, 1, 1, 
   3, 0, 0,   0, 3, 0,   0, 0, 3,   1, 2, 0,   2, 1, 0, 
   2, 0, 1,   1, 0, 2,   0, 1, 2,   0, 2, 1,   1, 1, 1
};
//*************************************************************************************************
gaussWaveFunc::~gaussWaveFunc()
//*************************************************************************************************
{
   dealloc1DStringArray(title);
   dealloc1DRealArray(R);
   dealloc1DStringArray(atLbl);
   dealloc1DRealArray(atCharge);
   dealloc1DIntArray(primCent);
   dealloc1DIntArray(primType);
   dealloc1DRealArray(primExp);
   dealloc1DRealArray(chi);
   dealloc1DRealArray(cab);
   dealloc1DRealArray(MOCoeff);
   dealloc1DRealArray(occN);
   dealloc1DRealArray(MOEner);
   dealloc1DRealArray(gx);
   dealloc1DRealArray(gy);
   dealloc1DRealArray(gz);
   dealloc1DRealArray(hxx);
   dealloc1DRealArray(hyy);
   dealloc1DRealArray(hzz);
   dealloc1DRealArray(hxy);
   dealloc1DRealArray(hxz);
   dealloc1DRealArray(hyz);
   cleanWaveFunctionInGPU(&h_R,&d_R,&h_a,&d_a,&d_e,&d_c,&d_chi,&d_rho,&d_aux);
   imldd=false;
}
//*************************************************************************************************
bool gaussWaveFunc::readFromFileWFN(string inname)
//*************************************************************************************************
{
   ifstream tif;
   tif.open(inname.c_str(),ios::in);
   if (!(tif.good())) {
      cout << "Error: File " << inname << "could not be opened...\n";
#if DEBUG
      cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
      return false;
   }
   tif.seekg(tif.beg);
   nTit=1;
   processFirstDataStringinWFNFile(tif,title,orbDesc,nMOr,nPri,nNuc);
   processCentersWFN(tif,nNuc,atLbl,R,atCharge);
   processPrimitivesWFN(tif,nPri,primCent,primType,primExp);
   processMolecularOrbitalPropsAndCoefs(tif,nMOr,nPri,occN,MOEner,MOCoeff);
   string liend;
   getline(tif,liend);
   //cout << "nPri%5: " << (nPri%5) << " len: " << liend.length() << endl;
   if (((nPri%5)==0)&&(liend.length()==0)) {
      getline(tif,liend);
      cout << liend << endl;
   }
   if (liend.substr(0,8)!="END DATA") {
      cout << "Error, expecting \"END DATA\" in file " << inname << endl;
      cout << "Line: " << liend << endl;
      return false;
   }
   getEnergyAndVirial(tif,totener,virial);
   allocAuxArrays();
   countPrimsPerCenter();
   calcCab();
   tif.close();
   imldd=testSupport();
   return true;
}
//*************************************************************************************************
bool gaussWaveFunc::readFromFileWFX(string inname)
//*************************************************************************************************
{
   ifstream tif;
   tif.open(inname.c_str(),ios::in);
   if (!(tif.good())) {
      cout << "Error: File " << inname << "could not be opened...\n";
#if DEBUG
      cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
      return false;
   }
   getTitleFromFileWFX(tif,nTit,title);
   getKeyWordsFromFileWFX(tif,orbDesc);
   if (orbDesc.substr(0,3)!="GTO") {
      cout << "Error: not supported wave function. Keyword: " << orbDesc << endl;
   }
   getNofNucleiFromFileWFX(tif,nNuc);
   getNofMolOrbFromFileWFX(tif,nMOr);
   getNofPrimFromFileWFX(tif,nPri);
   alloc1DStringArray("atLbl",nNuc,atLbl);
   alloc1DIntArray("primType",nPri,primType);
   alloc1DIntArray("primCent",nPri,primCent);
   alloc1DRealArray("R",(3*nNuc),R);
   alloc1DRealArray("atCharge",nNuc,atCharge);
   alloc1DRealArray("primExp",nPri,primExp);
   alloc1DRealArray("MOCoeff",(nMOr*nPri),MOCoeff);
   alloc1DRealArray("occN",nMOr,occN);
   alloc1DRealArray("MOEner",nMOr,MOEner);
   allocAuxArrays();
   getAtLabelsFromFileWFX(tif,nNuc,atLbl);
   getNucCartCoordsFromFileWFX(tif,nNuc,R);
   getAtChargesFromFileWFX(tif,nNuc,atCharge);
   getPrimCentersFromFileWFX(tif,nPri,primCent);
   getPrimTypesFromFileWFX(tif,nPri,primType);
   getPrimExponentsFromFileWFX(tif,nPri,primExp);
   getMolecOrbOccNumsFromFileWFX(tif,nMOr,occN);
   getMolecOrbEnergiesFromFileWFX(tif,nMOr,MOEner);
   getMolecOrbCoefficientsFromFileWFX(tif,nMOr,nPri,MOCoeff);
   getTotEnerAndVirialFromFileWFX(tif,totener,virial);
   countPrimsPerCenter();
   calcCab();
   tif.close();
   imldd=testSupport();
   return true;
}
//*************************************************************************************************
bool gaussWaveFunc::readFromFile(string inname)
{
   string extension;
   extension=inname.substr(inname.length()-3,3);
   if ((extension=="wfn")||(extension=="WFN")) {
      return readFromFileWFN(inname);
   } else if ((extension=="wfx")||(extension=="WFX")) {
      return readFromFileWFX(inname);
   } else {
      cout << "Error: unknown extension ("  << inname << ")!\nNothig to do, returning false...\n";
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
}
//*************************************************************************************************
bool gaussWaveFunc::testSupport()
{
   for (int i=0; i<nMOr; i++) {
      if (primType[i]>=MAXPRIMTYPEDEFINED) {
         cout << "Only " << MAXPRIMTYPEDEFINED << " types have been implemented in this version\n";
#if DEBUG
         cout << __FILE__ << "line " << __LINE__ << endl;
#endif
         return false;
      }
   }
   return true;
}
//*************************************************************************************************
void gaussWaveFunc::calcCab(void)
{
   int idx,indc;
   idx=0;
   if (nPri>MAXNUMBEROFPRIMITIVESFORMEMALLOC) {
      real memest=real(nPri*(nPri+12)*8)/real(1024*1024);
      char goon='n';
      cout << "The number of primitives is " << nPri <<". This will use approximatedly" << endl;
      cout << memest << "MB of RAM memory. Continue anyway (y/n)?" << endl;
      cin >> goon;
      if ((goon=='n')||(goon=='N')) {
         cout << "Perhaps you may want to recompile this program increasing the maximum number " << endl
              << "  of primitives. " << endl;
         exit(1);
      }
   }
   alloc1DRealArray(string("cab"),(nPri*nPri),cab);
   for (int i=0; i<nPri; i++) {
      for (int j=0; j<nPri; j++) {
         cab[idx]=0.0000000e0;
         for (int oi=0; oi<nMOr; oi++) {
            indc=oi*nPri;
            cab[idx]+=(occN[oi]*MOCoeff[indc+i]*MOCoeff[indc+j]);
         }
         idx++;
      }
   }
   return;
}
//*************************************************************************************************
void gaussWaveFunc::countPrimsPerCenter(void)
{
   alloc1DIntArray(string("myPN"),nNuc,myPN);
   for (int i=0; i<nPri; i++) {myPN[primCent[i]]++;}
   return;
}
//*************************************************************************************************
bool gaussWaveFunc::allocAuxArrays(void)
{
   bool allgood;
   allgood=alloc1DRealArray("chi",nPri,chi);
   if (!allgood) {
      cout << "Something wrong in allocating chi..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("gx"),nPri,gx);
   if (!allgood) {
      cout << "Something wrong in allocating gx..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("gy"),nPri,gy);
   if (!allgood) {
      cout << "Something wrong in allocating gy..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("gz"),nPri,gz);
   if (!allgood) {
      cout << "Something wrong in allocating gz..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hxx"),nPri,hxx);
   if (!allgood) {
      cout << "Something wrong in allocating hxx..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hyy"),nPri,hyy);
   if (!allgood) {
      cout << "Something wrong in allocating hyy..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hzz"),nPri,hzz);
   if (!allgood) {
      cout << "Something wrong in allocating hzz..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hxy"),nPri,hxy);
   if (!allgood) {
      cout << "Something wrong in allocating hxy..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hxz"),nPri,hxz);
   if (!allgood) {
      cout << "Something wrong in allocating hxz..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hyz"),nPri,hyz);
   if (!allgood) {
      cout << "Something wrong in allocating hyz..." << endl;
      return allgood;
   }
   return allgood;
}
//*************************************************************************************************
extern "C" void setupWaveFunctionInGPU(int npr, \
                                       real** host_R,real** dev_R, \
                                       int ** host_a,int ** dev_a, \
                                                     real** dev_e, \
                                                     real** dev_c, \
                                       real** dev_chi,real** dev_rho,real** dev_aux)
{
   //printf("Entering allocation function...\n");
   unsigned int mem_size=sizeof(real)*3*npr;
   *host_R=(real*)malloc(mem_size);
   checkCudaErrors(cudaMalloc((void **) &(*dev_R), mem_size));
   mem_size=sizeof(int)*3*npr;
   *host_a=(int*)malloc(mem_size);
   checkCudaErrors(cudaMalloc((void **) &(*dev_a), mem_size));
   mem_size=sizeof(real)*npr;
   checkCudaErrors(cudaMalloc((void **) &(*dev_e), mem_size));
   mem_size=sizeof(real)*npr*npr;
   checkCudaErrors(cudaMalloc((void **) &(*dev_c), mem_size));
   mem_size=sizeof(real)*npr;
   checkCudaErrors(cudaMalloc((void **) &(*dev_chi), mem_size));
   checkCudaErrors(cudaMalloc((void **) &(*dev_aux), mem_size));
   checkCudaErrors(cudaMalloc((void **) &(*dev_rho), sizeof(real)));
   //printf("Done\n");
   return;
}
//**************************************************************************************************
extern "C" void cleanWaveFunctionInGPU(real **host_R,real** dev_R,\
                                       int **host_a,int **dev_a, \
                                                    real**dev_e, \
                                                    real**dev_c, \
                                       real**dev_chi,real**dev_rho,real**dev_aux)
{
   if (*host_R!=NULL) {
      free(*host_R);
      *host_R=NULL;
   }
   checkCudaErrors(cudaFree(*dev_R));
   if (*host_a!=NULL) {
      free(*host_a);
      *host_a=NULL;
   }
   checkCudaErrors(cudaFree(*dev_a));
   checkCudaErrors(cudaFree(*dev_e));
   checkCudaErrors(cudaFree(*dev_c));
   checkCudaErrors(cudaFree(*dev_chi));
   checkCudaErrors(cudaFree(*dev_aux));
   checkCudaErrors(cudaFree(*dev_rho));
   return;
}
//**************************************************************************************************
bool gaussWaveFunc::setupGPU(void)
{
   //cout << "Checkpoint..." << endl;
   //cout << "h_R: " << h_R << ", d_R: " << d_R << endl;
   //cout << "h_a: " << h_a << ", d_a: " << d_a << endl;
   setupWaveFunctionInGPU(nPri,&h_R,&d_R,&h_a,&d_a,&d_e,&d_c,&d_chi,&d_rho,&d_aux);
   //cout << "h_R: " << h_R << ", d_R: " << d_R << endl;
   //cout << "h_a: " << h_a << ", d_a: " << d_a << endl;
   size_t idc,idr,idp;
   for (size_t i=0; i<nPri; i++) {
      idc=primCent[i]*3;
      idr=i*3;
      for (size_t j=0; j<3; j++) {
         h_R[idr+j]=R[idc+j];
      }
      idc=primType[i]*3;
      idp=i*3;
      for (size_t j=0; j<3; j++) {
         h_a[idp+j]=prTy[idc+j];
      }
   }
   copyDataToGPU(nPri,&h_R,&d_R,&h_a,&d_a,&primExp,&d_e,&cab,&d_c);
   //dosomethinginGPU(nPri,&d_R,&d_a,&d_e,&d_c);
   //copyDataFromGPU(nPri,&d_R,&h_R,&d_a,&h_a,&d_e,&primExp,&d_c,&cab);
   return true;
}
//**************************************************************************************************
real gaussWaveFunc::evalDensity(real x,real y,real z)
{
   return evalDensityInGPU(x,y,z,nPri,&d_R,&d_a,&d_e,&d_c,&chi,&d_chi,&d_rho);
}
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
extern "C" void copyDataToGPU(int npr,real **host_R,real **dev_R,int **host_a,int **dev_a, \
                                      real **host_e,real **dev_e,real**host_c,real**dev_c)
{
   unsigned int mem_size=sizeof(real)*3*npr;
   checkCudaErrors(cudaMemcpy(*dev_R, *host_R, mem_size, cudaMemcpyHostToDevice));
   mem_size=sizeof(int)*3*npr;
   checkCudaErrors(cudaMemcpy(*dev_a, *host_a, mem_size, cudaMemcpyHostToDevice));
   mem_size=sizeof(real)*npr;
   checkCudaErrors(cudaMemcpy(*dev_e, *host_e, mem_size, cudaMemcpyHostToDevice));
   mem_size=sizeof(real)*npr*npr;
   checkCudaErrors(cudaMemcpy(*dev_c, *host_c, mem_size, cudaMemcpyHostToDevice));
   return;
}
//**************************************************************************************************
extern "C" void copyDataFromGPU(int npr,real **dev_R,real **host_R,int **dev_a,int **host_a, \
                                        real **dev_e,real **host_e,real**dev_c,real**host_c)
{
   unsigned int mem_size=sizeof(real)*3*npr;
   checkCudaErrors(cudaMemcpy(*host_R, *dev_R, mem_size,cudaMemcpyDeviceToHost));
   mem_size=sizeof(int)*3*npr;
   checkCudaErrors(cudaMemcpy(*host_a, *dev_a, mem_size,cudaMemcpyDeviceToHost));
   mem_size=sizeof(real)*npr;
   checkCudaErrors(cudaMemcpy(*host_e, *dev_e, mem_size,cudaMemcpyDeviceToHost));
   mem_size=sizeof(real)*npr*npr;
   checkCudaErrors(cudaMemcpy(*host_c, *dev_c, mem_size,cudaMemcpyDeviceToHost));
   return;
}
//**************************************************************************************************
__global__ void kernelR(int rdim,real *dev_R)
{
   unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;;
   if (tid<rdim) {
      dev_R[tid]=1.0e0+0.1*real(tid);
   }
   return;
}
//**************************************************************************************************
__global__ void kernela(int rdim,int *dev_a)
{
   unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
   if (tid<rdim) {
      //dev_a[tid]=tid;
      dev_a[tid]+=10;
   }
   return;
}
//**************************************************************************************************
__global__ void kernele(int rdim,real *dev_e)
{
   unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
   if (tid<rdim) {
      dev_e[tid]*=(-1);
   }
   return;
}
//**************************************************************************************************
__global__ void kernelc(int rdim,real *dev_c)
{
   unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
   if (tid<rdim) {
      dev_c[tid]*=(-1);
   }
   return;
}
//**************************************************************************************************
__global__ void krnlCalcChi(int rdim,real *dev_R,int *dev_a,real *dev_e,real *dev_chi, \
                            real x,real y,real z)
{
   unsigned int btid = blockDim.x * blockIdx.x + threadIdx.x;
   unsigned int ttid = 3*btid;
   real r2,xmr,ymr,zmr,ang;
   if (btid<rdim) {
      xmr=x-dev_R[ttid];
      ymr=y-dev_R[ttid+1];
      zmr=z-dev_R[ttid+2];
      r2=(xmr*xmr)+(ymr*ymr)+(zmr*zmr);
      ang=1.0e0;
      for (int i=0; i<dev_a[ttid]; i++) {
         ang*=xmr;
      }
      for (int i=0; i<dev_a[ttid+1]; i++) {
         ang*=ymr;
      }
      for (int i=0; i<dev_a[ttid+2]; i++) {
         ang*=zmr;
      }
      dev_chi[btid]=ang*__expf(-r2*dev_e[btid]);
   }
   return;
}
//**************************************************************************************************
__global__ void krnlCalcRho(int rdim,int offset,real *dev_c,real *dev_chi,real *dev_rho)
{
   __shared__ real tmp[BASETHREADSPERBLOCK];
   unsigned int ptid = blockDim.x * blockIdx.x + threadIdx.x;
   
   if (ptid<rdim) {
      tmp[threadIdx.x]=dev_chi[ptid]*dev_c[offset + ptid];
   } else {
      tmp[threadIdx.x]=0.0e0;
   }
   __syncthreads();
   //if ( ptid == 0 ) {*dev_rho=0.0e0;}
   if (threadIdx.x == 0) {
      real sum=0.0e0;
      for (unsigned int j=0; j<BASETHREADSPERBLOCK; j++) {
         sum+=tmp[j];
      }
      //*dev_rho+=sum;
      atomicAdd(dev_rho,sum);
   }
}
//**************************************************************************************************
__device__ double atomicAdd(double* address, double val)
{
   unsigned long long int* address_as_ull =
   (unsigned long long int*)address;
   unsigned long long int old = *address_as_ull, assumed;
   do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed,
                      __double_as_longlong(val +
                                           __longlong_as_double(assumed)));
   } while (assumed != old);
   return __longlong_as_double(old);
}
//**************************************************************************************************
//**************************************************************************************************
extern "C" void dosomethinginGPU(int npr,real **dev_R,int **dev_a,real **dev_e,real **dev_c)
{
   int numElements=3*npr;
   int threadsPerBlock = BASETHREADSPERBLOCK;
   int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
   //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
   //printf("numElements= %d\n",numElements);
   kernelR<<<blocksPerGrid, threadsPerBlock>>>(numElements,*dev_R);
   kernela<<<blocksPerGrid, threadsPerBlock>>>(numElements,*dev_a);
   numElements=npr;
   blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
   kernele<<<blocksPerGrid, threadsPerBlock>>>(numElements,*dev_e);
   numElements=npr*npr;
   blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
   kernelc<<<blocksPerGrid, threadsPerBlock>>>(numElements,*dev_c);
   return;
}
//**************************************************************************************************
extern "C" real evalDensityInGPU(real x,real y,real z,
                                 int npr,real **dev_R,int **dev_a,real **dev_e,real **dev_c, \
                                 real **host_chi,real **dev_chi,real **dev_rho)
{
   int numElements=npr;
   int threadsPerBlock = BASETHREADSPERBLOCK;
   int blocksPerGrid =(numElements + threadsPerBlock - 1) / threadsPerBlock;
   krnlCalcChi<<<blocksPerGrid, threadsPerBlock>>>(numElements,*dev_R,*dev_a,*dev_e,*dev_chi, \
                                                   x,y,z);
   unsigned int mem_size=sizeof(real)*npr;
   checkCudaErrors(cudaMemcpy(*host_chi, *dev_chi, mem_size, cudaMemcpyDeviceToHost));
   real rho=0.0e0,trho;
   int offst;
   for (int i=0; i<npr; i++) {
      offst=i*npr;
      trho=0.0e0;
      checkCudaErrors(cudaMemcpy(*dev_rho, &trho, sizeof(real), cudaMemcpyHostToDevice));
      krnlCalcRho<<<blocksPerGrid, threadsPerBlock>>>(numElements,offst,*dev_c,*dev_chi,*dev_rho);
      checkCudaErrors(cudaMemcpy(&trho, *dev_rho, sizeof(real), cudaMemcpyDeviceToHost));
      rho+=trho*((*host_chi)[i]);
   }
   return rho;
   //unsigned int mem_size=sizeof(real)*npr;
   //checkCudaErrors(cudaMemcpy(*host_chi, *dev_chi, mem_size,cudaMemcpyDeviceToHost));
}
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
#endif//_SOLWAVEFUNCTIONCLASS_CPP_

