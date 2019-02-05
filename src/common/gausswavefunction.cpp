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

#ifndef _GAUSSWAVEFUNCTION_CPP_
#define _GAUSSWAVEFUNCTION_CPP_

#include <cstdlib>
using std::exit;
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
using std::ios;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <cmath>

#include "gausswavefunction.h"
#include "solscrutils.h"
#include "solmemhand.h"
#include "iofuncts-wfn.h"
#include "iofuncts-wfx.h"
#include "eig2-4.h"
#include "solmath.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#if DEBUG
#include "solstringtools.h"
#endif

#ifndef EPSFORELFVALUE
#define EPSFORELFVALUE (2.871e-05)
#endif

#ifndef EPSFORLOLVALUE
#define EPSFORLOLVALUE (2.871e-05)
#endif

#ifndef EPSFORMEPVALUE
#define EPSFORMEPVALUE (1.0e-8)
#endif


/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
GaussWaveFunction::GaussWaveFunction()
/* ************************************************************************************** */
{
   title=NULL;
   orbDesc=string("");
   nTit=0;
   nNuc=0;
   nMOr=0;
   nPri=0;
   EDFPri=0;
   totPri=0;
   coreElec=0;
   atLbl=NULL;
   primType=NULL;
   primCent=NULL;
   myPN=NULL;
   R=NULL;
   atCharge=NULL;
   primExp=NULL;
   MOCoeff=NULL;
   EDFCoeff=NULL;
   occN=NULL;
   MOEner=NULL;
   cab=NULL;
   chi=NULL;
   gx=gy=gz=NULL;
   hxx=hyy=hzz=NULL;
   hxy=hxz=hyz=NULL;
   totener=0.00e0;
   virial=0.0e0;
   nciRhoMin=NCIRHOMIN;
   nciRhoMax=NCIRHOMAX;
   nciSMax=NCISMAX;
   imldd=ihaveEDF=false;
   usescustfld=usevcustfld=false;
}
/* ************************************************************************************** */
int GaussWaveFunction::prTy[]={
   0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 0, 1,   2, 0, 0, 
   0, 2, 0,   0, 0, 2,   1, 1, 0,   1, 0, 1,   0, 1, 1, 
   3, 0, 0,   0, 3, 0,   0, 0, 3,   1, 2, 0,   2, 1, 0, 
   2, 0, 1,   1, 0, 2,   0, 1, 2,   0, 2, 1,   1, 1, 1,
   4, 0, 0,   0, 4, 0,   0, 0, 4,   3, 1, 0,   3, 0, 1,
   1, 3, 0,   0, 3, 1,   1, 0, 3,   0, 1, 3,   2, 2, 0,
   2, 0, 2,   0, 2, 2,   2, 1, 1,   1, 2, 1,   1, 1, 2,
   0, 0, 5,   0, 1, 4,   0, 2, 3,   0, 3, 2,   0, 4, 1,
   0, 5, 0,   1, 0, 4,   1, 1, 3,   1, 2, 2,   1, 3, 1,
   1, 4, 0,   2, 0, 3,   2, 1, 2,   2, 2, 1,   2, 3, 0,
   3, 0, 2,   3, 1, 1,   3, 2, 0,   4, 0, 1,   4, 1, 0,
   5, 0, 0
};
/* ************************************************************************************** */
GaussWaveFunction::~GaussWaveFunction()
/* ************************************************************************************** */
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
   if ( ihaveEDF ) {
      dealloc1DRealArray(EDFCoeff);
   }
   imldd=false;
}
/* ************************************************************************************** */
solreal GaussWaveFunction::getCoef(const int orbn,const int primn)
{
#if DEBUG
   if (orbn>=nMOr) {
      cout << "Error: attempting to look for a non existent Molecular Orbital!\nReturning zero...\n";
      return 0.0000e0;
   }
#endif
   return MOCoeff[(orbn*nPri)+primn];
}
/* *********************************************************************************** */
solreal GaussWaveFunction::getR(const int nucnum,const int cart)
{
#if DEBUG
   if (nucnum>=nNuc) {
      cout << "Error: attempting to look for a non-existent Nucleus!\nReturning zero...\n";
      return 0.00000e0;
   }
#endif
   return R[3*nucnum+cart];
}
/* ************************************************************************************** */
void GaussWaveFunction::getAng(int pt,int (&tt)[3])
{
   int tty=pt*3;
   tt[0]=prTy[tty+0];
   tt[1]=prTy[tty+1];
   tt[2]=prTy[tty+2];
}
/* ************************************************************************************** */
bool GaussWaveFunction::sameMolOrbOccNums()
{
   if (!occN) {
      cout << "Error: the wave function has not been properly allocated!\nReturning false...\n";
      return false;
   }
   for (int i=1; i<nMOr; i++) {if (occN[0]!=occN[i]) {return false;}}
   return true;
}
/* ************************************************************************************** */
bool GaussWaveFunction::readFromFileWFN(string inname)
/* ************************************************************************************** */
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
   totPri=nPri;
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
/* ************************************************************************************** */
bool GaussWaveFunction::readFromFileWFX(string inname)
/* ************************************************************************************** */
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
   totPri=nPri;
   getEDFExistenceFromFileWFX(tif,ihaveEDF);
   if ( ihaveEDF ) {
      int kk;
      countEDFCentersFromFileWFX(tif,kk);
      if ( kk>0 ) {
         cerr << "Error: In this version"
            <<" only combined EDF wave functions are supported!"
            << endl;
         return false;
      }
      getNofEDFPrimFromFileWFX(tif,1,EDFPri); //1: number of times to seek for
                     //EDF primitives. 1 is used for combined edfs, otherwise,
                     //it should be the number obtained wiht 
                     //countEDFCentersFromFileWFX.
      totPri+=EDFPri;
      alloc1DRealArray("EDFCoeff",EDFPri,EDFCoeff);
   }
   alloc1DStringArray("atLbl",nNuc,atLbl);
   alloc1DIntArray("primType",totPri,primType);
   alloc1DIntArray("primCent",totPri,primCent);
   alloc1DRealArray("R",(3*nNuc),R);
   alloc1DRealArray("atCharge",nNuc,atCharge);
   alloc1DRealArray("primExp",totPri,primExp);
   alloc1DRealArray("MOCoeff",(nMOr*nPri),MOCoeff);
   alloc1DRealArray("occN",(nMOr+1),occN); //The entry occN[nMOr] will save the
                    //number of core electrons when EDF information is present
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
   if ( ihaveEDF ) {
      getEDFPrimCentersFromFileWFX(tif,nPri,totPri,primCent);
      getEDFPrimTypesFromFileWFX(tif,nPri,totPri,primType);
      getEDFPrimExponentsFromFileWFX(tif,nPri,totPri,primExp);
      getNofCoreElectronsFromFileWFX(tif,coreElec);
      occN[nMOr]=1.0e0;
      getEDFPrimCoefficientsFromFileWFX(tif,EDFPri,EDFCoeff);
   }
   countPrimsPerCenter();
   calcCab();
   tif.close();
   imldd=testSupport();
   return true;
}
/* ************************************************************************************** */
bool GaussWaveFunction::readFromFile(string inname)
{
   string extension;
   extension=inname.substr(inname.length()-3,3);
   bool res;
   if ((extension=="wfn")||(extension=="WFN")) {
      res=readFromFileWFN(inname);
   } else if ((extension=="wfx")||(extension=="WFX")) {
      res=readFromFileWFX(inname);
   } else {
      cout << "Error: unknown extension ("  << inname << ")!\nNothig to do, returning false...\n";
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   if ( !res ) { return res; }
   res=(res&&sanityChecks());
   return res;
}
/* ************************************************************************************** */
bool GaussWaveFunction::sanityChecks(void)
{
   bool ret=true;
   /* Warns about primitive types.  */
   maxPrimType=0;
   for ( int i=0 ; i<nPri ; ++i ) {
      if ( primType[i]>maxPrimType ) { maxPrimType=primType[i]; }
   }
   if ( maxPrimType>19 ) {
      displayWarningMessage("Unfortunately, in current version ("
            CURRENTVERSION
            "),\nfunctions requiring computation of Boys Function"
            "\nterms are not implemented for"
            "\nprimitives with high angular momenta (i.e., primitive type > 20)");
   }
   return ret;
}
/* ************************************************************************************** */
void GaussWaveFunction::countPrimsPerCenter(void)
{
   alloc1DIntArray(string("myPN"),nNuc,myPN);
   for (int i=0; i<nPri; i++) {myPN[primCent[i]]++;}
   return;
}
/* ************************************************************************************** */
bool GaussWaveFunction::writeToFileWFX(string outname)
/* ************************************************************************************** */
{
   cout << "The function writeToFileWFX is not yet implemented!\nNothing done to file "
        << outname << "... \n";
   return true;
}
/* ************************************************************************************** */
bool GaussWaveFunction::testSupport()
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
/* ************************************************************************************** */
solreal GaussWaveFunction::evalDensityArray(solreal x,solreal y, solreal z)
{
   int indr,indc,indp;
   solreal xmr,ymr,zmr,rho,chia,chib,chit;
   rho=0.000000e0;
   for (int oi=0; oi<nMOr; oi++) {
      indc=oi*nPri;
      solreal rr;
      indr=0;
      indp=0;
      for (int i=0; i<nNuc; i++) {
         xmr=x-R[indr++];
         ymr=y-R[indr++];
         zmr=z-R[indr++];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         for (int j=0; j<myPN[i]; j++) {
            chi[indp]=evalAngACases(primType[indp],xmr,ymr,zmr);
            chi[indp]*=exp(primExp[indp]*rr);
            chi[indp++]*=MOCoeff[indc++];
         }
      }
      chia=0.000000000e0;
      for (int i=0; i<nPri; i++) {
         chib=0.0000000e0;
         chit=chi[i];
         for (int j=i+1; j<nPri; j++) {
            chib+=chi[j];
         }
         chib*=2.000000e0;
         chia+=(chit*(chit+chib));
      }
      rho+=(occN[oi]*chia);
   }
   return rho;
}
/* ************************************************************************************** */
void GaussWaveFunction::displayAllFieldProperties(solreal x,solreal y,solreal z)
{
   static solreal rho,lol,xx[3],g[3],hess[3][3];
   static solreal eivec[3][3],eival[3];
   xx[0]=x;
   xx[1]=y;
   xx[2]=z;
   evalRhoGradRho(x,y,z,rho,g);
   evalHessian(xx[0],xx[1],xx[2],hess);
   eigen_decomposition3(hess, eivec, eival);
   cout << scientific << setprecision(12);
   cout << "            R: " << setw(20) << x << setw(20) << y << setw(20) << z << endl;
   cout << "          Rho: " << setw(20) << rho << endl;
   cout << "      GradRho: " << setw(20) << g[0] << setw(20) << g[1] << setw(20) << g[2] << endl;
   cout << "    |GradRho|: " << setw(20) << sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]) << endl;
   cout << "      HessRho: " << setw(20) << hess[0][0] << setw(20) << hess[0][1]
                             << setw(20) << hess[0][2] << endl;
   cout << "               " << setw(20) << hess[1][0] << setw(20) << hess[1][1]
                             << setw(20) << hess[1][2] << endl;
   cout << "               " << setw(20) << hess[2][0] << setw(20) << hess[2][1]
                             << setw(20) << hess[2][2] << endl;
   cout << "  EigVal Hess: " << setw(20) << eival[0] << setw(20) << eival[1]
                             << setw(20) << eival[2] << endl;
   cout << "  EigVec Hess: " << setw(20) << eivec[0][0] << setw(20) << eivec[1][0]
                             << setw(20) << eivec[2][0] << endl;
   cout << "               " << setw(20) << eivec[0][1] << setw(20) << eivec[1][1]
                             << setw(20) << eivec[2][1] << endl;
   cout << "               " << setw(20) << eivec[0][2] << setw(20) << eivec[1][2]
                             << setw(20) << eivec[2][2] << endl;
   cout << "       LapRho: " << setw(20) << (hess[0][0]+hess[1][1]+hess[2][2]) << endl;
   evalHessLOL(xx,lol,g,hess);//(x,lol,gl,hl)
   eigen_decomposition3(hess, eivec, eival);
   cout << "          LOL: " << setw(20) << lol << endl;
   cout << "      gradLOL: " << setw(20) << g[0] << setw(20) << g[1] << setw(20) << g[2] << endl;
   cout << "    |gradLOL|: " << setw(20) << sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]) << endl;
   cout << "      HessLOL: " << setw(20) << hess[0][0] << setw(20) << hess[0][1]
                             << setw(20) << hess[0][2] << endl;
   cout << "               " << setw(20) << hess[1][0] << setw(20) << hess[1][1]
                             << setw(20) << hess[1][2] << endl;
   cout << "               " << setw(20) << hess[2][0] << setw(20) << hess[2][1]
                             << setw(20) << hess[2][2] << endl;
   cout << "  EigVal HLOL: " << setw(20) << eival[0] << setw(20) << eival[1]
                             << setw(20) << eival[2] << endl;
   cout << "  EigVec HLOL: " << setw(20) << eivec[0][0] << setw(20) << eivec[1][0]
                             << setw(20) << eivec[2][0] << endl;
   cout << "               " << setw(20) << eivec[0][1] << setw(20) << eivec[1][1]
                             << setw(20) << eivec[2][1] << endl;
   cout << "               " << setw(20) << eivec[0][2] << setw(20) << eivec[1][2]
                             << setw(20) << eivec[2][2] << endl;
   cout << "          ELF: " << setw(20) << evalELF(x,y,z) << endl;
   cout << "      K.E. G.: " << setw(20) << evalKineticEnergyG(x,y,z) << endl;
   cout << "      K.E. K.: " << setw(20) << evalKineticEnergyK(x,y,z) << endl;
   cout << "  Shann. Ent.: " << setw(20) << evalShannonEntropy(x,y,z) << endl;
   cout << "  Elect. Pot.: " << setw(20) << evalMolElecPot(x,y,z) << endl;
   evalLED(xx,g);
   cout << "          LED: " << setw(20) << g[0] << setw(20) << g[1] << setw(20) << g[2] << endl;
   cout << "        |LED|: " << setw(20) << sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]) << endl;
   cout << "  RedDensGrad: " << setw(20) << evalReducedDensityGradient(xx[0],xx[1],xx[2]) << endl;
   cout << "         RoSE: " << setw(20) << evalRoSE(xx[0],xx[1],xx[2]) << endl;
   cout << "     V.P.E.D.: " << setw(20) << evalVirialPotentialEnergyDensity(x,y,z) << endl;
   if ( usescustfld ) {
      cout << "Cust. S. Field: " << setw(20) << evalCustomScalarField(xx[0],xx[1],xx[2]) << endl;
   }
   if ( usevcustfld ) {
      evalCustomVectorField(xx[0],xx[1],xx[2],g);
      cout << "Cust. V. Field: " <<  setw(20) << g[0] << setw(20) << g[1] 
                                 << setw(20) << g[2] << endl;
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::writeAllFieldProperties(solreal x,solreal y,solreal z,ofstream &ofil)
{
   static solreal rho,lol,xx[3],g[3],hess[3][3];
   static solreal eivec[3][3],eival[3];
   xx[0]=x;
   xx[1]=y;
   xx[2]=z;
   evalRhoGradRho(x,y,z,rho,g);
   evalHessian(xx[0],xx[1],xx[2],hess);
   eigen_decomposition3(hess, eivec, eival);
   ofil << scientific << setprecision(12);
   ofil << "  R:           " << setw(20) << x << setw(20) << y << setw(20) << z << endl;
   ofil << "  Rho:         " << setw(20) << rho << endl;
   ofil << "  GradRho:     " << setw(20) << g[0] << setw(20) << g[1] << setw(20) << g[2] << endl;
   ofil << "  |GradRho|:   " << setw(20) << sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]) << endl;
   ofil << "  HessRho:     " << setw(20) << hess[0][0] << setw(20) << hess[0][1]
                             << setw(20) << hess[0][2] << endl;
   ofil << "               " << setw(20) << hess[1][0] << setw(20) << hess[1][1]
                             << setw(20) << hess[1][2] << endl;
   ofil << "               " << setw(20) << hess[2][0] << setw(20) << hess[2][1]
                             << setw(20) << hess[2][2] << endl;
   ofil << "  EigVal Hess: " << setw(20) << eival[0] << setw(20) << eival[1]
                             << setw(20) << eival[2] << endl;
   ofil << "  EigVec Hess: " << setw(20) << eivec[0][0] << setw(20) << eivec[1][0]
                             << setw(20) << eivec[2][0] << endl;
   ofil << "               " << setw(20) << eivec[0][1] << setw(20) << eivec[1][1]
                             << setw(20) << eivec[2][1] << endl;
   ofil << "               " << setw(20) << eivec[0][2] << setw(20) << eivec[1][2]
                             << setw(20) << eivec[2][2] << endl;
   ofil << "  LapRho:      " << setw(20) << (hess[0][0]+hess[1][1]+hess[2][2]) << endl;
   evalHessLOL(xx,lol,g,hess);//(x,lol,gl,hl)
   eigen_decomposition3(hess, eivec, eival);
   ofil << "  LOL:         " << setw(20) << lol << endl;
   ofil << "  gradLOL:     " << setw(20) << g[0] << setw(20) << g[1] << setw(20) << g[2] << endl;
   ofil << "  |gradLOL|:   " << setw(20) << sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]) << endl;
   ofil << "  HessLOL:     " << setw(20) << hess[0][0] << setw(20) << hess[0][1]
                             << setw(20) << hess[0][2] << endl;
   ofil << "               " << setw(20) << hess[1][0] << setw(20) << hess[1][1]
                             << setw(20) << hess[1][2] << endl;
   ofil << "               " << setw(20) << hess[2][0] << setw(20) << hess[2][1]
                             << setw(20) << hess[2][2] << endl;
   ofil << "  EigVal HLOL: " << setw(20) << eival[0] << setw(20) << eival[1]
                             << setw(20) << eival[2] << endl;
   ofil << "  EigVec HLOL: " << setw(20) << eivec[0][0] << setw(20) << eivec[1][0]
                             << setw(20) << eivec[2][0] << endl;
   ofil << "               " << setw(20) << eivec[0][1] << setw(20) << eivec[1][1]
                             << setw(20) << eivec[2][1] << endl;
   ofil << "               " << setw(20) << eivec[0][2] << setw(20) << eivec[1][2]
                             << setw(20) << eivec[2][2] << endl;
   ofil << "  ELF:         " << setw(20) << evalELF(x,y,z) << endl;
   ofil << "  K.E. G.:     " << setw(20) << evalKineticEnergyG(x,y,z) << endl;
   ofil << "  K.E. K.:     " << setw(20) << evalKineticEnergyK(x,y,z) << endl;
   ofil << "  Shann. Ent.: " << setw(20) << evalShannonEntropy(x,y,z) << endl;
   ofil << "  Elect. Pot.: " << setw(20) << evalMolElecPot(x,y,z) << endl;
   evalLED(xx,g);
   ofil << "  LED:         " << setw(20) << g[0] << setw(20) << g[1] << setw(20) << g[2] << endl;
   ofil << "  |LED|:       " << setw(20) << sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]) << endl;
   ofil << "  RedDensGrad: " << setw(20) << evalReducedDensityGradient(x,y,z) << endl;
   ofil << "  RoSE:        " << setw(20) << evalRoSE(x,y,z) << endl;
   ofil << "  V.P.E.D.:    " << setw(20) << evalVirialPotentialEnergyDensity(x,y,z) << endl;
   if ( usescustfld ) {
      ofil << "Cust. S. Field: " << setw(20) << evalCustomScalarField(xx[0],xx[1],xx[2]) << endl;
   }
   if ( usevcustfld ) {
      evalCustomVectorField(xx[0],xx[1],xx[2],g);
      ofil << "Cust. V. Field: " <<  setw(20) << g[0] << setw(20) << g[1] 
                                 << setw(20) << g[2] << endl;
   }
return;
}
/* ************************************************************************************** */
/*
solreal GaussWaveFunction::evalPrimCases(int &pty, solreal &alp, solreal x, solreal y, solreal z)
{
   solreal xx,yy,zz,rr,pv;
   xx=x*x;
   yy=y*y;
   zz=z*z;
   rr=xx+yy+zz;
   pv=exp(-alp*rr);
   switch (pty) {
      case 0:
         break;
      case 1:
         pv*=x;
         break;
      case 2:
         pv*=y;
         break;
      case 3:
         pv*=z;
         break;
      case 4:
         pv*=xx;
         break;
      case 5:
         pv*=yy;
         break;
      case 6:
         pv*=zz;
         break;
      case 7:
         pv*=x;
         pv*=y;
         break;
      case 8:
         pv*=x;
         pv*=z;
         break;
      case 9:
         pv*=y;
         pv*=z;
         break;
      case 10:
         pv*=xx;
         pv*=x;
         break;
      case 11:
         pv*=yy;
         pv*=y;
         break;
      case 12:
         pv*=zz;
         pv*=z;
         break;
      case 13:
         pv*=x;
         pv*=yy;
         break;
      case 14:
         pv*=xx;
         pv*=y;
         break;
      case 15:
         pv*=xx;
         pv*=z;
         break;
      case 16:
         pv*=x;
         pv*=zz;
         break;
      case 17:
         pv*=y;
         pv*=zz;
         break;
      case 18:
         pv*=yy;
         pv*=z;
         break;
      case 19:
         pv*=x;
         pv*=y;
         pv*=z;
         break;
      default:
         cout << "Not implemented angular function! (a=" << pty << ")\n";
         break;
   }
   return pv;
}
// */
/* ************************************************************************************** */
/* Preliminary tests indicate that, surprisingly, the case choosing is slightly
 * slower (around 1-2%) than the brute for(...) {pv*=x_i;}. Tested with phenantrene and 
 * f2.g09.wfn  */
solreal GaussWaveFunction::evalAngACases(int &pty, solreal x, solreal y, solreal z)
{
   solreal pv=1.00000000e0;
   //*
   int a[3];
   getAng(pty,a);
   for ( int i=0 ; i<a[0] ; ++i ) { pv*=x; }
   for ( int i=0 ; i<a[1] ; ++i ) { pv*=y; }
   for ( int i=0 ; i<a[2] ; ++i ) { pv*=z; }
   // */
   /*
   switch (pty) {
      case 0:
         break;
      case 1:
         pv*=x;
         break;
      case 2:
         pv*=y;
         break;
      case 3:
         pv*=z;
         break;
      case 4:
         pv*=(x*x);
         break;
      case 5:
         pv*=(y*y);
         break;
      case 6:
         pv*=(z*z);
         break;
      case 7:
         pv*=x;
         pv*=y;
         break;
      case 8:
         pv*=x;
         pv*=z;
         break;
      case 9:
         pv*=y;
         pv*=z;
         break;
      case 10:
         pv*=(x*x*x);
         break;
      case 11:
         pv*=(y*y*y);
         break;
      case 12:
         pv*=(z*z*z);
         break;
      case 13:
         pv*=x;
         pv*=(y*y);
         break;
      case 14:
         pv*=(x*x);
         pv*=y;
         break;
      case 15:
         pv*=(x*x);
         pv*=z;
         break;
      case 16:
         pv*=x;
         pv*=(z*z);
         break;
      case 17:
         pv*=y;
         pv*=(z*z);
         break;
      case 18:
         pv*=(y*y);
         pv*=z;
         break;
      case 19:
         pv*=x;
         pv*=y;
         pv*=z;
         break;
      case 20 :
      case 21 :
      case 22 :
      case 23 :
      case 24 :
      case 25 :
      case 26 :
      case 27 :
      case 28 :
      case 29 :
      case 30 :
      case 31 :
      case 32 :
      case 33 :
      case 34 :
      case 35 :
         int a[3];
         getAng(pty,a);
         for ( int i=0 ; i<a[0] ; ++i ) { pv*=x; }
         for ( int i=0 ; i<a[1] ; ++i ) { pv*=y; }
         for ( int i=0 ; i<a[2] ; ++i ) { pv*=z; }
         break;
      default:
         cout << "Not implemented angular function! (t=" << (pty+1) << ")\n";
         break;
   }
   // */
   return pv;
}
/* ************************************************************************************** */
void GaussWaveFunction::calcCab(void)
{
   int idx,indc;
   idx=0;
   if (nPri>MAXNUMBEROFPRIMITIVESFORMEMALLOC) {
      solreal memest=solreal(nPri*(nPri+12)*8)/solreal(1024*1024);
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
/* ************************************************************************************** */
#if PARALLELISEDTK
solreal GaussWaveFunction::evalDensity(solreal x,solreal y,solreal z)
{
   int indr,indp;
   solreal xmr,ymr,zmr,rho,chib;
   rho=0.000000e0;
   solreal rr;
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         chi[indp]=evalAngACases(primType[indp],xmr,ymr,zmr);
         chi[indp]*=exp(primExp[indp]*rr);
         indp++;
      }
   }
   indr=0;
   rho=0.000000e0;
   solreal chii;
   int i,j;
#pragma omp parallel for private(chii,indr,chib) \
firstprivate(j) lastprivate(i) reduction(+: rho)
   for (i=0; i<nPri; i++) {
      indr=i*(nPri);
      chib=0.0000000e0;
      for (j=(i+1); j<nPri; j++) {
         chib+=(cab[indr+j]*chi[j]);
      }
      chii=chi[i];
      rho+=(cab[indr+i]*chii*chii+2.00000000e0*chib*chii);
   }
   return rho;
}
#else
solreal GaussWaveFunction::evalDensity(solreal x,solreal y,solreal z)
{
   int indr,indp;
   solreal xmr,ymr,zmr,rho,chib;
   rho=0.000000e0;
   solreal rr;
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         chi[indp]=evalAngACases(primType[indp],xmr,ymr,zmr);
         chi[indp]*=exp(primExp[indp]*rr);
         indp++;
      }
   }
   /*
   indr=0;
   rho=0.000000e0;
   for (int i=0; i<nPri; i++) {
      indr=i*(nPri+1);
      rho+=(cab[indr++]*chi[i]*chi[i]);
      chib=0.0000000e0;
      for (int j=(i+1); j<nPri; j++) {
         chib+=(cab[indr++]*chi[j]);
      }
      rho+=(2.00000000e0*chib*chi[i]);
   }
   // */
   int lowPri=nPri-(nPri%4);
   rho=0.000000e0;
   indr=-1;
   for ( int i=0 ; i<nPri ; ++i ) {
      chib=0.0e0;
      for ( int j=0 ; j<lowPri ; j+=4 ) {
         chib+=(cab[++indr]*chi[j  ]);
         chib+=(cab[++indr]*chi[j+1]);
         chib+=(cab[++indr]*chi[j+2]);
         chib+=(cab[++indr]*chi[j+3]);
      }
      for ( int j=lowPri ; j<nPri ; ++j ) {
         chib+=(cab[++indr]*chi[j]);
      }
      rho+=chib*chi[i];
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         //cout << "pc: " << primCent[i] << ", rr: " << rr << endl;
         chi[i]=evalAngACases(primType[i],xmr,ymr,zmr);
         chi[i]*=exp(primExp[i]*rr);
      }
      chib=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         chib+=(EDFCoeff[i-nPri]*chi[i]);
      }
      //rho+=occN[nMOr]*chib*chib;
      rho+=chib;
   }
   return rho;
}
#endif
/* ************************************************************************************** */
solreal GaussWaveFunction::evalOptimizedScalar(solreal px,solreal py,solreal pz)
{
   int indr,indp,ppt;
   solreal rhop,Rx[3],alp;
   complex<solreal> chit,chiu;
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      Rx[0]=R[indr++];
      Rx[1]=R[indr++];
      Rx[2]=R[indr++];
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         evalFTChi(ppt,alp,Rx,px,py,pz,chit);
         chi[indp]=chit.real();
         gx[indp]=chit.imag();
         indp++;
      }
   }
   rhop=0.000000e0;
   /*
   complex<solreal> chiv;
   indr=0;
   for (int i=0; i<nPri; i++) {
      indr=i*(nPri+1);
      chiu=complex<solreal>(chi[i],gx[i]);
      //chiu.real(chi[i]); chiu.imag(gx[i]);
      rhop+=(cab[indr++]*norm(chiu));
      for (int j=(i+1); j<nPri; j++) {
         chiv=complex<solreal>(chi[j],gx[j]);
         //chiv.real(chi[j]); chiv.imag(gx[j]);
         chit=(chiv*conj(chiu));
         chit+=(conj(chiv)*chiu);
         rhop+=(cab[indr++]*(chit.real()));
      }
   }
   // */
   solreal sumim,sumre,cc;
   indr=0;
   for ( int i=0 ; i<nPri ; ++i ) {
      sumim=sumre=0.0e0;
      for ( int j=0 ; j<nPri ; ++j ) {
         cc=cab[indr++];
         sumre+=cc*chi[j];
         sumim+=cc*gx[j];
      }
      rhop+=sumre*chi[i];
      rhop+=sumim*gx[i];
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         Rx[0]=R[indr];
         Rx[1]=R[indr+1];
         Rx[2]=R[indr+2];
         ppt=primType[i];
         alp=0.5e0*primExp[i];
         evalFTChi(ppt,alp,Rx,px,py,pz,chit);
         chi[i]=chit.real();
         gx[i]=chit.imag();
      }
      for (int i=nPri; i<totPri; ++i) {
         chiu=complex<solreal>(chi[i],gx[i]);
         //chiu.real(chi[i]); chiu.imag(gx[i]);
         rhop+=(EDFCoeff[i-nPri]*norm(chiu));
      }
   }
   return rhop;
}
/* ************************************************************************************** */
#if PARALLELISEDTK
void GaussWaveFunction::evalRhoGradRho(solreal x, solreal y, solreal z,solreal &rho, solreal &dx, solreal &dy, solreal &dz)
{
   solreal nabx,naby,nabz,xmr,ymr,zmr,trho,cc,rr,alp,chib;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   nabx=naby=nabz=0.000000000000000e0;
   indp=0;
   trho=0.0000000e0;
   int i,j;
#pragma omp parallel for private(indp,chib) \
firstprivate(j) lastprivate(i) reduction(+: trho,nabx,naby,nabz)
   for (i=0; i<nPri; i++) {
      indp=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      for (j=0; j<nPri; j++) {
         chib+=(chi[j]*cab[indp+j]);
      }
      trho+=(chib*chi[i]);
      nabx+=(chib*gx[i]);
      naby+=(chib*gy[i]);
      nabz+=(chib*gz[i]);
   }
   dx=2.00000e0*nabx;
   dy=2.00000e0*naby;
   dz=2.00000e0*nabz;
   rho=trho;
   return;
}
#else
/* ************************************************************************************** */
void GaussWaveFunction::evalRhoGradRho(solreal x, solreal y, solreal z,solreal &rho, solreal &dx, solreal &dy, solreal &dz)
{
   solreal nabx,naby,nabz,xmr,ymr,zmr,trho,cc,rr,alp,chib;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         ++indp;
      }
   }
   nabx=naby=nabz=0.000000000000000e0;
   indp=0;
   trho=0.0000000e0;
   int lowPri=nPri-(nPri%4);
   nabx=naby=nabz=0.000000000000000e0;
   indp=0;
   trho=0.0000000e0;
   //*
   for (int i=0; i<nPri; i++) {
      chib=0.0000000e0;
      for (int j=0; j<lowPri; j+=4) {
         chib+=(cab[indp++]*chi[j  ]);
         chib+=(cab[indp++]*chi[j+1]);
         chib+=(cab[indp++]*chi[j+2]);
         chib+=(cab[indp++]*chi[j+3]);
      }
      for (int j=lowPri; j<nPri; j++) {
         chib+=(cab[indp++]*chi[j]);
      }
      trho+=(chib*chi[i]);
      nabx+=(chib*gx[i]);
      naby+=(chib*gy[i]);
      nabz+=(chib*gz[i]);
   }
   // */
   /*
   for (int i=0; i<nPri; i++) {
      //indr=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      for (int j=0; j<nPri; j++) {
         //cc=cab[indp++];
         chib+=(chi[j]*cab[indp++]);
      }
      trho+=(chib*chi[i]);
      nabx+=(chib*gx[i]);
      naby+=(chib*gy[i]);
      nabz+=(chib*gz[i]);
   }
   // */
   dx=2.00000e0*nabx;
   dy=2.00000e0*naby;
   dz=2.00000e0*nabz;
   rho=trho;
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      chib=nabx=naby=nabz=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         nabx+=(cc*gx[i]);
         naby+=(cc*gy[i]);
         nabz+=(cc*gz[i]);
      }
      /*
      rho+=occN[nMOr]*chib*chib;
      dx+=2.0e0*occN[nMOr]*chib*nabx;
      dy+=2.0e0*occN[nMOr]*chib*naby;
      dz+=2.0e0*occN[nMOr]*chib*nabz;
      // */
      rho+=chib;
      dx+=nabx;
      dy+=naby;
      dz+=nabz;
   }
   return;
}
#endif
/* ************************************************************************************** */
void GaussWaveFunction::evalOptimizedVectorScalar(solreal x, solreal y, solreal z,solreal &rho, solreal &dx, solreal &dy, solreal &dz)
{
   solreal nabx,naby,nabz,xmr,ymr,zmr,trho,cc,rr,alp,chib;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   int lowPri=nPri-(nPri%4);
   nabx=naby=nabz=0.000000000000000e0;
   indp=0;
   trho=0.0000000e0;
   for (int i=0; i<nPri; i++) {
      //indr=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      for (int j=0; j<lowPri; j+=4) {
         chib+=(cab[indp++]*chi[j  ]);
         chib+=(cab[indp++]*chi[j+1]);
         chib+=(cab[indp++]*chi[j+2]);
         chib+=(cab[indp++]*chi[j+3]);
      }
      for (int j=lowPri; j<nPri; ++j) {
         chib+=(cab[indp++]*chi[j]);
      }
      trho+=(chib*chi[i]);
      nabx+=(chib*gx[i]);
      naby+=(chib*gy[i]);
      nabz+=(chib*gz[i]);
   }
   dx=2.00000e0*nabx;
   dy=2.00000e0*naby;
   dz=2.00000e0*nabz;
   rho=trho;
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      chib=nabx=naby=nabz=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         nabx+=(cc*gx[i]);
         naby+=(cc*gy[i]);
         nabz+=(cc*gz[i]);
      }
      /*
      rho+=occN[nMOr]*chib*chib;
      dx+=2.0e0*occN[nMOr]*chib*nabx;
      dy+=2.0e0*occN[nMOr]*chib*naby;
      dz+=2.0e0*occN[nMOr]*chib*nabz;
      // */
      rho+=chib;
      dx+=nabx;
      dy+=naby;
      dz+=nabz;
   }
   return;
}
/* ************************************************************************************** */
bool GaussWaveFunction::allocAuxArrays(void)
{
   bool allgood;
   allgood=alloc1DRealArray("chi",totPri,chi);
   if (!allgood) {
      cout << "Something wrong in allocating chi..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("gx"),totPri,gx);
   if (!allgood) {
      cout << "Something wrong in allocating gx..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("gy"),totPri,gy);
   if (!allgood) {
      cout << "Something wrong in allocating gy..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("gz"),totPri,gz);
   if (!allgood) {
      cout << "Something wrong in allocating gz..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hxx"),totPri,hxx);
   if (!allgood) {
      cout << "Something wrong in allocating hxx..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hyy"),totPri,hyy);
   if (!allgood) {
      cout << "Something wrong in allocating hyy..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hzz"),totPri,hzz);
   if (!allgood) {
      cout << "Something wrong in allocating hzz..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hxy"),totPri,hxy);
   if (!allgood) {
      cout << "Something wrong in allocating hxy..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hxz"),totPri,hxz);
   if (!allgood) {
      cout << "Something wrong in allocating hxz..." << endl;
      return allgood;
   }
   allgood=alloc1DRealArray(string("hyz"),totPri,hyz);
   if (!allgood) {
      cout << "Something wrong in allocating hyz..." << endl;
      return allgood;
   }
   return allgood;
}
/* ************************************************************************************** */
/* ************************************************************************************** */
void GaussWaveFunction::evalDkAngCases(int &pty,solreal alp,solreal x, solreal y, solreal z, solreal &anx, solreal &any, solreal &anz)
{
   solreal cc=(-2.0000000e0*alp);
   /* As opposed to the case of evalAngACases, tests indicate that a combination
    * of cases-choose for small angular momentum, and the generic algorithm for
    * large angular exponent is the fastets choice.  */
   switch (pty) {
      case 0:
         anx=cc*x;
         any=cc*y;
         anz=cc*z;
         break;
      case 1:
         cc*=x;
         anx=1.0000000e0+cc*x;
         any=cc*y;
         anz=cc*z;
         break;
      case 2:
         cc*=y;
         anx=cc*x;
         any=1.0000000e0+cc*y;
         anz=cc*z;
         break;
      case 3:
         cc*=z;
         anx=cc*x;
         any=cc*y;
         anz=1.0000000e0+cc*z;
         break;
      case 4:
         cc*=(x*x);
         anx=(2.000000e0+cc)*x;
         any=cc*y;
         anz=cc*z;
         break;
      case 5:
         cc*=(y*y);
         anx=cc*x;
         any=(2.00000e0+cc)*y;
         anz=cc*z;
         break;
      case 6:
         cc*=(z*z);
         anx=cc*x;
         any=cc*y;
         anz=(2.00000e0+cc)*z;
         break;
      case 7:
         anx=(1.00000e0+cc*x*x)*y;
         any=(1.00000e0+cc*y*y)*x;
         anz=cc*x*y*z;
         break;
      case 8:
         anx=(1.00000e0+cc*x*x)*z;
         any=cc*x*y*z;
         anz=(1.00000e0+cc*z*z)*x;
         break;
      case 9:
         anx=cc*x*y*z;
         any=(1.00000e0+cc*y*y)*z;
         anz=(1.00000e0+cc*z*z)*y;
         break;
      case 10:
         cc*=(x*x);
         anx=(3.00000e0+cc)*x*x;
         any=cc*x*y;
         anz=cc*x*z;
         break;
      case 11:
         cc*=(y*y);
         anx=cc*x*y;
         any=(3.00000e0+cc)*y*y;
         anz=cc*z*y;
         break;
      case 12:
         cc*=(z*z);
         anx=cc*x*z;
         any=cc*y*z;
         anz=(3.00000e0+cc)*z*z;
         break;
      case 13:
         cc*=(y*y);
         anx=(y*y+cc*x*x);
         any=(2.00000e0+cc)*x*y;
         anz=cc*x*z;
         break; 
      case 14:
         cc*=(x*x);
         anx=(2.00000e0+cc)*x*y;
         any=x*x+cc*y*y;
         anz=cc*y*z;
         break;
      case 15:
         cc*=(x*x);
         anx=(2.00000e0+cc)*x*z;
         any=cc*y*z;
         anz=(x*x+cc*z*z);
         break;
      case 16:
         cc*=(z*z);
         anx=z*z+cc*x*x;
         any=cc*x*y;
         anz=(2.00000e0+cc)*x*z;
         break;
      case 17:
         cc*=(z*z);
         anx=cc*x*y;
         any=z*z+cc*y*y;
         anz=(2.00000e0+cc)*y*z;
         break;
      case 18:
         cc*=(y*y);
         anx=cc*x*z;
         any=(2.00000e0+cc)*y*z;
         anz=y*y+cc*z*z;
         break;
      case 19:
         cc*=(x*y*z);
         anx=y*z+cc*x;
         any=x*z+cc*y;
         anz=x*y+cc*z;
         break;
      default:
         {
            int a[3];
            getAng(pty,a);
            solreal tmpx,tmpy,tmpz;
            evald1SingCartA(a[0],cc,x,tmpx,anx);
            evald1SingCartA(a[1],cc,y,tmpy,any);
            evald1SingCartA(a[2],cc,z,tmpz,anz);
            anx*=(tmpy*tmpz);
            any*=(tmpx*tmpz);
            anz*=(tmpx*tmpy);
         }
         break;
   }
   return;
}
/* *********************************************************************************** */
/* Some loop tests shows that using directly evald2SingCartA (as opposed to
 * use cases) for evalDkDlAngCases is slightly faster (~0.5-1.0%).
 * For the time being, JMSA does not advise to remove old code, since
 * this is thoroughly tested (numerically), thus it is reliable.
 * So far, evald2SingCartA reproduces without fault values at single
 * points.  */
void GaussWaveFunction::evalDkDlAngCases(int &pty,solreal alp,solreal x,solreal y,solreal z,
      solreal &axx,solreal &ayy,solreal &azz,solreal &axy,solreal &axz,solreal &ayz)
{
   solreal ta=(-2.00000e0*alp);
   //*
   int a[3];
   getAng(pty,a);
   solreal d0x,d0y,d0z,d1x,d2x,d1y,d2y,d1z,d2z;
   evald2SingCartA(a[0],ta,x,d0x,d1x,d2x);
   evald2SingCartA(a[1],ta,y,d0y,d1y,d2y);
   evald2SingCartA(a[2],ta,z,d0z,d1z,d2z);
   axx=d2x*d0y*d0z;
   ayy=d0x*d2y*d0z;
   azz=d0x*d0y*d2z;
   axy=d1x*d1y*d0z;
   axz=d1x*d0y*d1z;
   ayz=d0x*d1y*d1z;
   // */
   /*
   solreal fa2,x2,y2,z2;
   fa2=ta*ta;
   x2=x*x;
   y2=y*y;
   z2=z*z;
   switch (pty) {
      case 0:
         axx=ta+fa2*x2;
         ayy=ta+fa2*y2;
         azz=ta+fa2*z2;
         axy=fa2*x*y;
         axz=fa2*x*z;
         ayz=fa2*y*z;
         break;
      case 1:
         axx=(3.00000e0*ta+fa2*x2)*x;
         ayy=(ta+fa2*y2)*x;
         azz=(ta+fa2*z2)*x;
         axy=(ta+fa2*x2)*y;
         axz=(ta+fa2*x2)*z;
         ayz=fa2*x*y*z;
         break;
      case 2:
         axx=(ta+fa2*x2)*y;
         ayy=(3.00000e0*ta+fa2*y2)*y;
         azz=(ta+fa2*z2)*y;
         axy=(ta+fa2*y2)*x;
         axz=fa2*x*y*z;
         ayz=(ta+fa2*y2)*z;
         break;
      case 3:
         axx=(ta+fa2*x2)*z;
         ayy=(ta+fa2*y2)*z;
         azz=(3.00000e0*ta+fa2*z2)*z;
         axy=fa2*x*y*z;
         axz=(ta+fa2*z2)*x;
         ayz=(ta+fa2*z2)*y;
         break;
      case 4:
         axx=2.00000e0+(5.00000e0*ta+fa2*x2)*x2;
         ayy=(ta+fa2*y2)*x2;
         azz=(ta+fa2*z2)*x2;
         axy=(2.00000e0*ta+fa2*x2)*x*y;
         axz=(2.00000e0*ta+fa2*x2)*x*z;
         ayz=fa2*x2*y*z;
         break;
      case 5:
         axx=(ta+fa2*x2)*y2;
         ayy=2.00000e0+(5.00000e0*ta+fa2*y2)*y2;
         azz=(ta+fa2*z2)*y2;
         axy=(2.00000e0*ta+fa2*y2)*x*y;
         axz=fa2*x*y2*z;
         ayz=(2.00000e0*ta+fa2*y2)*y*z;
         break;
      case 6:
         axx=(ta+fa2*x2)*z2;
         ayy=(ta+fa2*y2)*z2;
         azz=2.00000e0+(5.00000e0*ta+fa2*z2)*z2;
         axy=fa2*x*y*z2;
         axz=(2.00000e0*ta+fa2*z2)*x*z;
         ayz=(2.00000e0*ta+fa2*z2)*y*z;
         break;
      case 7:
         axx=(3.00000e0*ta+fa2*x2)*x*y;
         ayy=(3.00000e0*ta+fa2*y2)*x*y;
         azz=(ta+fa2*z2)*x*y;
         axy=(1.00000e0+ta*x2)*(1.00000e0+ta*y2);
         axz=(ta+fa2*x2)*y*z;
         ayz=(ta+fa2*y2)*x*z;
         break;
      case 8:
         axx=(3.00000e0*ta+fa2*x2)*x*z;
         ayy=(ta+fa2*y2)*x*z;
         azz=(3.00000e0*ta+fa2*z2)*x*z;
         axy=(ta+fa2*x2)*y*z;
         axz=(1.00000e0+ta*x2)*(1.0000e0+ta*z2);
         ayz=(ta+fa2*z2)*x*y;
         break;
      case 9:
         axx=(ta+fa2*x2)*y*z;
         ayy=(3.00000e0*ta+fa2*y2)*y*z;
         azz=(3.00000e0*ta+fa2*z2)*y*z;
         axy=(ta+fa2*y2)*x*z;
         axz=(ta+fa2*z2)*x*y;
         ayz=(1.00000e0+ta*y2)*(1.00000e0+ta*z2);
         break;
      case 10:
         axx=(6.00000e0+(7.00000e0*ta+fa2*x2)*x2)*x;
         ayy=(ta+fa2*y2)*x2*x;
         azz=(ta+fa2*z2)*x2*x;
         axy=(3.00000e0*ta+fa2*x2)*x2*y;
         axz=(3.00000e0*ta+fa2*x2)*x2*z;
         ayz=fa2*x2*x*y*z;
         break;
      case 11:
         axx=(ta+fa2*x2)*y2*y;
         ayy=(6.00000e0+(7.00000e0*ta+fa2*y2)*y2)*y;
         azz=(ta+fa2*z2)*y2*y;
         axy=(3.00000e0*ta+fa2*y2)*y2*x;
         axz=fa2*x*y2*y*z;
         ayz=(3.00000e0*ta+fa2*y2)*y2*z;
         break;
      case 12:
         axx=(ta+fa2*x2)*z2*z;
         ayy=(ta+fa2*y2)*z2*z;
         azz=(6.00000e0+(7.000000e0*ta+fa2*z2)*z2)*z;
         axy=fa2*x*y*z2*z;
         axz=(3.00000e0*ta+fa2*z2)*z2*x;
         ayz=(3.00000e0*ta+fa2*z2)*z2*y;
         break;
      case 13:
         axx=(3.00000e0*ta+fa2*x2)*x*y2;
         ayy=(2.00000e0+(5.00000e0*ta+fa2*y2)*y2)*x;
         azz=(ta+fa2*z2)*x*y2;
         axy=(1.00000e0+ta*x2)*(2.00000e0+ta*y2)*y;
         axz=(ta+fa2*x2)*y2*z;
         ayz=(2.00000e0*ta+fa2*y2)*x*y*z;
         break;
      case 14:
         axx=(2.00000e0+(5.00000e0*ta+fa2*x2)*x2)*y;
         ayy=(3.00000e0*ta+fa2*y2)*x2*y;
         azz=(ta+fa2*z2)*x2*y;
         axy=(1.00000e0+ta*y2)*(2.00000e0+ta*x2)*x;
         axz=(2.00000e0*ta+fa2*x2)*x*y*z;
         ayz=(ta+fa2*y2)*x2*z;
         break;
      case 15:
         axx=(2.00000e0+(5.00000e0*ta+fa2*x2)*x2)*z;
         ayy=(ta+fa2*y2)*x2*z;
         azz=(3.00000e0*ta+fa2*z2)*x2*z;
         axy=(2.00000e0*ta+fa2*x2)*x*y*z;
         axz=(1.00000e0+ta*z2)*(2.00000e0+ta*x2)*x;
         ayz=(ta+fa2*z2)*x2*y;
         break;
      case 16:
         axx=(3.00000e0*ta+fa2*x2)*x*z2;
         ayy=(ta+fa2*y2)*x*z2;
         azz=(2.00000e0+(5.00000e0*ta+fa2*z2)*z2)*x;
         axy=(ta+fa2*x2)*y*z2;
         axz=(1.00000e0+ta*x2)*(2.00000e0+ta*z2)*z;
         ayz=(2.00000e0*ta+fa2*z2)*x*y*z;
         break;
      case 17:
         axx=(ta+fa2*x2)*y*z2;
         ayy=(3.00000e0*ta+fa2*y2)*y*z2;
         azz=(2.00000e0+(5.00000e0*ta+fa2*z2)*z2)*y;
         axy=(ta+fa2*y2)*x*z2;
         axz=(2.00000e0*ta+fa2*z2)*x*y*z;
         ayz=(1.00000e0+ta*y2)*(2.00000e0+ta*z2)*z;
         break;
      case 18:
         axx=(ta+fa2*x2)*y2*z;
         ayy=(2.00000e0+(5.00000e0*ta+fa2*y2)*y2)*z;
         azz=(3.00000e0*ta+fa2*z2)*y2*z;
         axy=(2.00000e0*ta+fa2*y2)*x*y*z;
         axz=(ta+fa2*z2)*x*y2;
         ayz=(1.00000e0+ta*z2)*(2.00000e0+ta*y2)*y;
         break;
      case 19:
         axx=(3.00000e0*ta+fa2*x2)*x*y*z;
         ayy=(3.00000e0*ta+fa2*y2)*x*y*z;
         azz=(3.00000e0*ta+fa2*z2)*x*y*z;
         axy=(1.00000e0+ta*x2)*(1.00000e0+ta*y2)*z;
         axz=(1.00000e0+ta*x2)*(1.00000e0+ta*z2)*y;
         ayz=(1.00000e0+ta*y2)*(1.00000e0+ta*z2)*x;
         break;
      default:
         {
            int a[3];
            getAng(pty,a);
            solreal d0x,d0y,d0z,d1x,d2x,d1y,d2y,d1z,d2z;
            evald2SingCartA(a[0],ta,x,d0x,d1x,d2x);
            evald2SingCartA(a[1],ta,y,d0y,d1y,d2y);
            evald2SingCartA(a[2],ta,z,d0z,d1z,d2z);
            axx=d2x*d0y*d0z;
            ayy=d0x*d2y*d0z;
            azz=d0x*d0y*d2z;
            axy=d1x*d1y*d0z;
            axz=d1x*d0y*d1z;
            ayz=d0x*d1y*d1z;
         }
         break;
   }
   // */
   return;
}
#if PARALLELISEDTK
/* ************************************************************************************** */
void GaussWaveFunction::evalHessian(solreal x, solreal y, solreal z,
                                solreal &dxx, solreal &dyy, solreal &dzz,
                                solreal &dxy, solreal &dxz, solreal &dyz)
{
   solreal nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,xmr,ymr,zmr,cc,rr,alp,chii,gxi,gyi,gzi;
   solreal sxx,syy,szz,sxy,sxz,syz,gxs,gys,gzs;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[indp],hyy[indp],hzz[indp],
                          hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   int i,j;
   nabxx=nabyy=nabzz=nabxy=nabxz=nabyz=0.000000000000000e0;
   //indr=0;
#pragma omp parallel for private(indr,cc,gxs,gys,gzs,sxx,syy,szz,sxy,sxz,syz,\
                                 gxi,gyi,gzi) \
firstprivate(j) lastprivate(i) reduction(+: nabxx,nabyy,nabzz,nabxy,nabxz,nabyz)
   for (i=0; i<nPri; i++) {
      indr=i*nPri;
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chii=0.00000e0;
      for (j=0; j<nPri; j++) {
         cc=cab[indr+j];
         chii+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
      }
      gxi=gx[i];
      gyi=gy[i];
      gzi=gz[i];
      nabxx+=chii*hxx[i]+gxi*gxs;
      nabyy+=chii*hyy[i]+gyi*gys;
      nabzz+=chii*hzz[i]+gzi*gzs;
      nabxy+=chii*hxy[i]+gxi*gys;
      nabxz+=chii*hxz[i]+gxi*gzs;
      nabyz+=chii*hyz[i]+gyi*gzs;
   }
   dxx=2.00000e0*nabxx;
   dyy=2.00000e0*nabyy;
   dzz=2.00000e0*nabzz;
   dxy=2.00000e0*nabxy;
   dxz=2.00000e0*nabxz;
   dyz=2.00000e0*nabyz;
   return;
}
#else
/* ************************************************************************************** */
void GaussWaveFunction::evalHessian(solreal x, solreal y, solreal z,
                                solreal &dxx, solreal &dyy, solreal &dzz,
                                solreal &dxy, solreal &dxz, solreal &dyz)
{
   solreal nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,xmr,ymr,zmr,cc,rr,alp,chii,gxi,gyi,gzi;
   solreal sxx,syy,szz,sxy,sxz,syz,gxs,gys,gzs;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[indp],hyy[indp],hzz[indp],
                          hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   nabxx=nabyy=nabzz=nabxy=nabxz=nabyz=0.000000000000000e0;
   indr=0;
   for (int i=0; i<nPri; i++) {
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chii=0.00000e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indr++];
         chii+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
      }
      gxi=gx[i];
      gyi=gy[i];
      gzi=gz[i];
      nabxx+=chii*hxx[i]+gxi*gxs;
      nabyy+=chii*hyy[i]+gyi*gys;
      nabzz+=chii*hzz[i]+gzi*gzs;
      nabxy+=chii*hxy[i]+gxi*gys;
      nabxz+=chii*hxz[i]+gxi*gzs;
      nabyz+=chii*hyz[i]+gyi*gzs;
   }
   dxx=2.00000e0*nabxx;
   dyy=2.00000e0*nabyy;
   dzz=2.00000e0*nabzz;
   dxy=2.00000e0*nabxy;
   dxz=2.00000e0*nabxz;
   dyz=2.00000e0*nabyz;
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[i],hyy[i],hzz[i],
                          hxy[i],hxz[i],hyz[i]);
         hxx[i]*=cc;
         hyy[i]*=cc;
         hzz[i]*=cc;
         hxy[i]*=cc;
         hxz[i]*=cc;
         hyz[i]*=cc;
      }
      chii=0.0e0;
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chii+=(cc*chi[i]);
         gxs+=(cc*gx[i]);
         gys+=(cc*gy[i]);
         gzs+=(cc*gz[i]);
         sxx+=(cc*hxx[i]);
         syy+=(cc*hyy[i]);
         szz+=(cc*hzz[i]);
         sxy+=(cc*hxy[i]);
         sxz+=(cc*hxz[i]);
         syz+=(cc*hyz[i]);
      }
      dxx+=(sxx);
      dyy+=(syy);
      dzz+=(szz);
      dxy+=(sxy);
      dxz+=(sxz);
      dyz+=(syz);
   }
return;
}
#endif
/* ************************************************************************************** */
solreal GaussWaveFunction::evalLapAngCases(int &pty,solreal alp,solreal x,solreal y,solreal z,solreal rr)
{
   solreal ta,fr;
   ta=(-2.00000e0*alp);
   fr=ta*ta;
   fr*=rr;
   switch (pty) {
      case 0:
         return (3.00000e0*ta+fr);
         break;
      case 1:
         return ((5.00000e0*ta+fr)*x);
         break;
      case 2:
         return ((5.00000e0*ta+fr)*y);
         break;
      case 3:
         return ((5.00000e0*ta+fr)*z);
         break;
      case 4:
         return (2.00000e0+(7.00000e0*ta+fr)*x*x);
         break;
      case 5:
         return (2.00000e0+(7.00000e0*ta+fr)*y*y);
         break;
      case 6:
         return (2.00000e0+(7.00000e0*ta+fr)*z*z);
         break;
      case 7:
         return ((7.00000e0*ta+fr)*x*y);
         break;
      case 8:
         return ((7.00000e0*ta+fr)*x*z);
         break;
      case 9:
         return ((7.00000e0*ta+fr)*y*z);
         break;
      case 10:
         return ((6.00000e0+(9.00000e0*ta+fr)*x*x)*x);
         break;
      case 11:
         return ((6.00000e0+(9.00000e0*ta+fr)*y*y)*y);
         break;
      case 12:
         return ((6.00000e0+(9.00000e0*ta+fr)*z*z)*z);
         break;
      case 13:
         return ((2.00000e0+(9.00000e0*ta+fr)*y*y)*x);
         break;
      case 14:
         return ((2.00000e0+(9.00000e0*ta+fr)*x*x)*y);
         break;
      case 15:
         return ((2.00000e0+(9.00000e0*ta+fr)*x*x)*z);
         break;
      case 16:
         return ((2.00000e0+(9.00000e0*ta+fr)*z*z)*x);
         break;
      case 17:
         return ((2.00000e0+(9.00000e0*ta+fr)*z*z)*y);
         break;
      case 18:
         return ((2.00000e0+(9.00000e0*ta+fr)*y*y)*z);
         break;
      case 19:
         return ((9.00000e0*ta+fr)*x*y*z);
         break;
      default:
         {
            int a[3];
            getAng(pty,a);
            solreal d0x,d0y,d0z,d1x,d2x,d1y,d2y,d1z,d2z;
            evald2SingCartA(a[0],ta,x,d0x,d1x,d2x);
            evald2SingCartA(a[1],ta,y,d0y,d1y,d2y);
            evald2SingCartA(a[2],ta,z,d0z,d1z,d2z);
            return (d2x*d0y*d0z+d0x*d2y*d0z+d0x*d0y*d2z);
         }
         break;
   }
   return 0.00000e0;
}
#if PARALLELISEDTK
/* ************************************************************************************** */
solreal GaussWaveFunction::evalLapRho(solreal x, solreal y, solreal z)
{
   solreal lap,xmr,ymr,zmr,cc,rr,alp;
   solreal sxx,gxs,gys,gzs;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(-alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         hxx[indp]=evalLapAngCases(ppt,alp,xmr,ymr,zmr,rr);
         hxx[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   lap=0.000000000000000e0;
   indr=0;
   int i,k;
#pragma omp parallel for private(sxx,cc,gxs,gys,gzs) \
firstprivate(k) lastprivate(i) reduction(+: lap)
   for (i=0; i<nPri; i++) {
      sxx=0.00000e0;
      indp=i*(nPri);
      for (k=0; k<nPri; k++) {
         cc=cab[indp+k];
         sxx+=cc*hxx[k];
      }
      lap+=chi[i]*sxx;
      cc=cab[indp+i];
      lap+=cc*gx[i]*gx[i];
      lap+=cc*gy[i]*gy[i];
      lap+=cc*gz[i]*gz[i];
      gxs=gys=gzs=0.00000e0;
      for (k=(i+1); k<nPri; k++) {
         cc=cab[indp+k];
         gxs+=cc*gx[k];
         gys+=cc*gy[k];
         gzs+=cc*gz[k];
      }
      lap+=gx[i]*gxs*2.00000e0;
      lap+=gy[i]*gys*2.00000e0;
      lap+=gz[i]*gzs*2.00000e0;
   }
   return (2.00000e0*lap);
}
#else
/* ************************************************************************************** */
solreal GaussWaveFunction::evalLapRho(solreal x, solreal y, solreal z)
{
   solreal lap,xmr,ymr,zmr,cc,rr,alp;
   solreal sxx,gxs,gys,gzs;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(-alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         hxx[indp]=evalLapAngCases(ppt,alp,xmr,ymr,zmr,rr);
         hxx[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   lap=0.000000000000000e0;
   indr=0;
   for (int i=0; i<nPri; i++) {
      sxx=0.00000e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indr++];
         sxx+=cc*hxx[j];
      }
      lap+=chi[i]*sxx;
      indp=i*(nPri+1);
      cc=cab[indp++];
      lap+=cc*gx[i]*gx[i];
      lap+=cc*gy[i]*gy[i];
      lap+=cc*gz[i]*gz[i];
      gxs=gys=gzs=0.00000e0;
      for (int k=(i+1); k<nPri; k++) {
         cc=cab[indp++];
         gxs+=cc*gx[k];
         gys+=cc*gy[k];
         gzs+=cc*gz[k];
      }
      lap+=gx[i]*gxs*2.00000e0;
      lap+=gy[i]*gys*2.00000e0;
      lap+=gz[i]*gzs*2.00000e0;
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(-alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         hxx[i]=evalLapAngCases(ppt,alp,xmr,ymr,zmr,rr);
         hxx[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      sxx=gxs=gys=gzs=0.0e0;
      solreal chib=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         sxx+=(cc*hxx[i]);
         gxs+=(cc*gx[i]);
         gys+=(cc*gy[i]);
         gzs+=(cc*gz[i]);
      }
      //lap+=(occN[nMOr]*(chib*sxx+gxs*gxs+gys*gys+gzs*gzs));
      lap+=(0.5e0*sxx);
   }
   return (2.00000e0*lap);
}
#endif
/* ************************************************************************************** */
void GaussWaveFunction::seekBondCP(int ii,int jj,solreal &r1,solreal &r2,solreal &r3,solreal &gx,solreal &gy,solreal &gz)
{
   solreal x[3],delta[3],g[3];
   for (int i=0; i<3; i++) {x[i]=0.500000e0*(R[3*ii+i]+R[3*jj+i]);}
   getBondCPStep(x,delta,g);
   solreal magd=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);;
   solreal magh=magd;
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&(count<MAXITERATIONBCPSEARCH)){
      for (int i=0; i<3; i++) {x[i]+=delta[i];}
      getBondCPStep(x,delta,g);
      magd=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      count++;
   }
   //cout << "nIt: " << count << endl;
#if DEBUG
   if (count==MAXITERATIONBCPSEARCH) {
      cout << "Maximum iteration reached (" << MAXITERATIONBCPSEARCH << ")!\n";
   }
#endif
   r1=x[0];
   r2=x[1];
   r3=x[2];
   gx=g[0];
   gy=g[1];
   gz=g[2];
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evalHessian(solreal x, solreal y, solreal z,solreal (&h)[3][3])
{
   evalHessian(x,y,z,h[0][0],h[1][1],h[2][2],h[0][1],h[0][2],h[1][2]);
   h[1][0]=h[0][1];
   h[2][0]=h[0][2];
   h[2][1]=h[1][2];
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evalHessian(solreal x, solreal y, solreal z,solreal &dens,solreal (&g)[3],solreal (&h)[3][3])
{
   solreal nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,xmr,ymr,zmr,cc,rr,alp,
   chii,gxi,gyi,gzi,rho,delx,dely,delz;
   solreal sxx,syy,szz,sxy,sxz,syz,gxs,gys,gzs;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[indp],hyy[indp],hzz[indp],
                          hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   rho=delx=dely=delz=0.00000e0;
   nabxx=nabyy=nabzz=nabxy=nabxz=nabyz=0.000000000000000e0;
   /*
   indr=0;
   for (int i=0; i<nPri; i++) {
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chii=0.00000e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indr++];
         chii+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
      }
      rho+=chii*chi[i];
      gxi=gx[i];
      gyi=gy[i];
      gzi=gz[i];
      delx+=(gxi*chii);
      dely+=(gyi*chii);
      delz+=(gzi*chii);
      nabxx+=chii*hxx[i]+gxi*gxs;
      nabyy+=chii*hyy[i]+gyi*gys;
      nabzz+=chii*hzz[i]+gzi*gzs;
      nabxy+=chii*hxy[i]+gxi*gys;
      nabxz+=chii*hxz[i]+gxi*gzs;
      nabyz+=chii*hyz[i]+gyi*gzs;
   }
   // */
   //*
   int lowPri=nPri-(nPri%4);
   indr=0;
   for (int i=0; i<nPri; i++) {
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chii=0.00000e0;
      for (int j=0; j<lowPri; j+=4) {
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j  ];
         gxs+=cc*gx[j  ];
         gys+=cc*gy[j  ];
         gzs+=cc*gz[j  ];
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j+1];
         gxs+=cc*gx[j+1];
         gys+=cc*gy[j+1];
         gzs+=cc*gz[j+1];
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j+2];
         gxs+=cc*gx[j+2];
         gys+=cc*gy[j+2];
         gzs+=cc*gz[j+2];
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j+3];
         gxs+=cc*gx[j+3];
         gys+=cc*gy[j+3];
         gzs+=cc*gz[j+3];
      }
      for (int j=lowPri; j<nPri; ++j) {
         cc=cab[indr++];
         chii+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
      }
      rho+=chii*chi[i];
      gxi=gx[i];
      gyi=gy[i];
      gzi=gz[i];
      delx+=(gxi*chii);
      dely+=(gyi*chii);
      delz+=(gzi*chii);
      nabxx+=chii*hxx[i]+gxi*gxs;
      nabyy+=chii*hyy[i]+gyi*gys;
      nabzz+=chii*hzz[i]+gzi*gzs;
      nabxy+=chii*hxy[i]+gxi*gys;
      nabxz+=chii*hxz[i]+gxi*gzs;
      nabyz+=chii*hyz[i]+gyi*gzs;
   }
   // */
   dens=rho;
   g[0]=2.00000e0*delx;
   g[1]=2.00000e0*dely;
   g[2]=2.00000e0*delz;
   h[0][0]=2.00000e0*nabxx;
   h[1][1]=2.00000e0*nabyy;
   h[2][2]=2.00000e0*nabzz;
   h[0][1]=2.00000e0*nabxy;
   h[0][2]=2.00000e0*nabxz;
   h[1][2]=2.00000e0*nabyz;
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[i],hyy[i],hzz[i],
                          hxy[i],hxz[i],hyz[i]);
         hxx[i]*=cc;
         hyy[i]*=cc;
         hzz[i]*=cc;
         hxy[i]*=cc;
         hxz[i]*=cc;
         hyz[i]*=cc;
      }
      chii=0.0e0;
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chii+=(cc*chi[i]);
         gxs+=(cc*gx[i]);
         gys+=(cc*gy[i]);
         gzs+=(cc*gz[i]);
         sxx+=(cc*hxx[i]);
         syy+=(cc*hyy[i]);
         szz+=(cc*hzz[i]);
         sxy+=(cc*hxy[i]);
         sxz+=(cc*hxz[i]);
         syz+=(cc*hyz[i]);
      }
      dens+=(chii);
      g[0]+=(gxs);
      g[1]+=(gys);
      g[2]+=(gzs);
      h[0][0]+=(sxx);
      h[1][1]+=(syy);
      h[2][2]+=(szz);
      h[0][1]+=(sxy);
      h[0][2]+=(sxz);
      h[1][2]+=(syz);
   }
   h[1][0]=h[0][1];
   h[2][0]=h[0][2];
   h[2][1]=h[1][2];
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evalOptimizedScalVecHess(solreal x, solreal y, solreal z,solreal &dens,solreal (&g)[3],solreal (&h)[3][3])
{
   solreal nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,xmr,ymr,zmr,cc,rr,alp,
   chii,gxi,gyi,gzi,rho,delx,dely,delz;
   solreal sxx,syy,szz,sxy,sxz,syz,gxs,gys,gzs;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[indp],hyy[indp],hzz[indp],
                          hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   rho=delx=dely=delz=0.00000e0;
   nabxx=nabyy=nabzz=nabxy=nabxz=nabyz=0.000000000000000e0;
   int lowPri=nPri-(nPri%4);
   indr=0;
   for (int i=0; i<nPri; i++) {
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chii=0.00000e0;
      for (int j=0; j<lowPri; j+=4) {
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j  ];
         gxs+=cc*gx[j  ];
         gys+=cc*gy[j  ];
         gzs+=cc*gz[j  ];
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j+1];
         gxs+=cc*gx[j+1];
         gys+=cc*gy[j+1];
         gzs+=cc*gz[j+1];
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j+2];
         gxs+=cc*gx[j+2];
         gys+=cc*gy[j+2];
         gzs+=cc*gz[j+2];
         //----------------
         cc=cab[indr++];
         chii+=cc*chi[j+3];
         gxs+=cc*gx[j+3];
         gys+=cc*gy[j+3];
         gzs+=cc*gz[j+3];
      }
      for (int j=lowPri; j<nPri; ++j) {
         cc=cab[indr++];
         chii+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
      }
      rho+=chii*chi[i];
      gxi=gx[i];
      gyi=gy[i];
      gzi=gz[i];
      delx+=(gxi*chii);
      dely+=(gyi*chii);
      delz+=(gzi*chii);
      nabxx+=chii*hxx[i]+gxi*gxs;
      nabyy+=chii*hyy[i]+gyi*gys;
      nabzz+=chii*hzz[i]+gzi*gzs;
      nabxy+=chii*hxy[i]+gxi*gys;
      nabxz+=chii*hxz[i]+gxi*gzs;
      nabyz+=chii*hyz[i]+gyi*gzs;
   }
   dens=rho;
   g[0]=2.00000e0*delx;
   g[1]=2.00000e0*dely;
   g[2]=2.00000e0*delz;
   h[0][0]=2.00000e0*nabxx;
   h[1][1]=2.00000e0*nabyy;
   h[2][2]=2.00000e0*nabzz;
   h[0][1]=2.00000e0*nabxy;
   h[0][2]=2.00000e0*nabxz;
   h[1][2]=2.00000e0*nabyz;
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[i],hyy[i],hzz[i],
                          hxy[i],hxz[i],hyz[i]);
         hxx[i]*=cc;
         hyy[i]*=cc;
         hzz[i]*=cc;
         hxy[i]*=cc;
         hxz[i]*=cc;
         hyz[i]*=cc;
      }
      chii=0.0e0;
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chii+=(cc*chi[i]);
         gxs+=(cc*gx[i]);
         gys+=(cc*gy[i]);
         gzs+=(cc*gz[i]);
         sxx+=(cc*hxx[i]);
         syy+=(cc*hyy[i]);
         szz+=(cc*hzz[i]);
         sxy+=(cc*hxy[i]);
         sxz+=(cc*hxz[i]);
         syz+=(cc*hyz[i]);
      }
      /*
      cc=occN[nMOr];
      dens+=(cc*chii*chii);
      cc*=2.0e0;
      g[0]+=(cc*chii*gxs);
      g[1]+=(cc*chii*gys);
      g[2]+=(cc*chii*gzs);
      h[0][0]+=(cc*(sxx*chii+gxs*gxs));
      h[1][1]+=(cc*(syy*chii+gys*gys));
      h[2][2]+=(cc*(szz*chii+gzs*gzs));
      h[0][1]+=(cc*(sxy*chii+gxs*gys));
      h[0][2]+=(cc*(sxz*chii+gxs*gzs));
      h[1][2]+=(cc*(syz*chii+gys*gzs));
      // */
      dens+=(chii);
      g[0]+=(gxs);
      g[1]+=(gys);
      g[2]+=(gzs);
      h[0][0]+=(sxx);
      h[1][1]+=(syy);
      h[2][2]+=(szz);
      h[0][1]+=(sxy);
      h[0][2]+=(sxz);
      h[1][2]+=(syz);
   }
   h[1][0]=h[0][1];
   h[2][0]=h[0][2];
   h[2][1]=h[1][2];
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evalRhoGradRho(solreal x, solreal y, solreal z,solreal &rho, solreal (&grd)[3])
{
   evalRhoGradRho(x,y,z,rho,grd[0],grd[1],grd[2]);
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::getBondCPStep(solreal (&x)[3],solreal (&hh)[3],solreal (&g)[3])
{
   static solreal hess[3][3],eive[3][3],b[3];
   solreal rho;
   evalRhoGradRho(x[0],x[1],x[2],rho,g[0],g[1],g[2]);
   evalHessian(x[0],x[1],x[2],hess);
   eigen_decomposition3(hess, eive, b);
   static solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {hess[i][j]=0.00000e0;}
   }
   solreal ln=0.5e0*(b[2]-sqrt(b[2]*b[2]+4.0e0*F[2]*F[2]));
   hess[0][0]=b[0];
   hess[1][1]=b[1];
   hess[2][0]=hess[0][2]=F[0];
   hess[2][1]=hess[1][2]=F[1];
   static solreal m3[3][3],vv[3];
   eigen_decomposition3(hess, m3, vv);
   solreal lp=vv[2];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp);
      hh[1]-=eive[1][j]*F[j]/(b[j]-lp);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZEBCPSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZEBCPSEARCH;
      }
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::getRingCPStep(solreal (&x)[3],solreal (&hh)[3],solreal (&g)[3])
{
   static solreal hess[3][3],eive[3][3],b[3];
   solreal rho;
   evalRhoGradRho(x[0],x[1],x[2],rho,g[0],g[1],g[2]);
   evalHessian(x[0],x[1],x[2],hess);
   eigen_decomposition3(hess, eive, b);
   static solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {hess[i][j]=0.00000e0;}
   }
   solreal lp=0.5e0*(b[0]+sqrt(b[0]*b[0]+4.0e0*F[0]*F[0]));
   hess[0][0]=b[1];
   hess[1][1]=b[2];
   hess[2][0]=hess[0][2]=F[1];
   hess[2][1]=hess[1][2]=F[2];
   static solreal m3[3][3],vv[3];
   eigen_decomposition3(hess, m3, vv);
   solreal ln=vv[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-lp);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZERCPSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZERCPSEARCH;
      }
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::seekRingCP(solreal &r1,solreal &r2,solreal &r3,solreal &gx,solreal &gy,solreal &gz)
{
   solreal x[3],delta[3],g[3];
   x[0]=r1;
   x[1]=r2;
   x[2]=r3;
   getRingCPStep(x,delta,g);
   solreal magd=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   solreal magh=magd;
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&(count<MAXITERATIONRCPSEARCH)) {
      for (int i=0; i<3; i++) {x[i]+=delta[i];}
      getRingCPStep(x,delta,g);
      magd=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      count++;
   }
#if DEBUG
   if (count==MAXITERATIONRCPSEARCH) {
      cout << "Maximum iteration reached (" << MAXITERATIONRCPSEARCH << ")!\n";
   }
#endif
   //cout << "nIt: " << count << endl;
   r1=x[0];
   r2=x[1];
   r3=x[2];
   gx=g[0];
   gy=g[1];
   gz=g[2];
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::seekCageCP(solreal &r1,solreal &r2,solreal &r3,solreal &gx,solreal &gy,solreal &gz)
{
   solreal x[3],delta[3],g[3];
   x[0]=r1;
   x[1]=r2;
   x[2]=r3;
   getCageCPStep(x,delta,g);
   solreal magd=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
   solreal magh=magd;
   int count=0;
   while (((magd>EPSGRADMAG)&&(magh>EPSGRADMAG))&&(count<MAXITERATIONCCPSEARCH)) {
      for (int i=0; i<3; i++) {x[i]+=delta[i];}
      getCageCPStep(x,delta,g);
      magd=sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
      count++;
   }
#if DEBUG
   if (count==MAXITERATIONCCPSEARCH) {
      cout << "Maximum iteration reached (" << MAXITERATIONCCPSEARCH << ")!\n";
   }
#endif
   //cout << "nIt: " << count << endl;
   r1=x[0];
   r2=x[1];
   r3=x[2];
   gx=g[0];
   gy=g[1];
   gz=g[2];
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::getCageCPStep(solreal (&x)[3],solreal (&hh)[3],solreal (&g)[3])
{
   static solreal hess[3][3],eive[3][3],b[3];
   solreal rho;
   evalRhoGradRho(x[0],x[1],x[2],rho,g[0],g[1],g[2]);
   evalHessian(x[0],x[1],x[2],hess);
   eigen_decomposition3(hess, eive, b);
   static solreal F[3];
   for (int i=0; i<3; i++) {
      F[i]=0.00000e0;
      for (int j=0; j<3; j++) {
         F[i]+=g[j]*eive[j][i];
      }
   }
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {hess[i][j]=0.00000e0;}
   }
   static solreal h4[4][4],m4[4][4],v4[4];
   for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
         h4[i][j]=0.0e0;
         m4[i][j]=0.0e0;
      }
      v4[i]=0.0e0;
   }
   h4[0][0]=b[0]; h4[1][1]=b[1]; h4[2][2]=b[2];
   h4[0][3]=h4[3][0]=F[0];
   h4[1][3]=h4[3][1]=F[1];
   h4[2][3]=h4[3][2]=F[2];
   eigen_decomposition4(h4, m4, v4);
   solreal ln=v4[0];
   hh[2]=hh[1]=hh[0]=0.00000e0;
   for (int j=0; j<3; j++) {
      hh[0]-=eive[0][j]*F[j]/(b[j]-ln);
      hh[1]-=eive[1][j]*F[j]/(b[j]-ln);
      hh[2]-=eive[2][j]*F[j]/(b[j]-ln);
   }
   for (int i=0; i<3; i++) {
      if (fabs(hh[i])>MAXSTEPSIZECCPSEARCH) {
         hh[i]=SIGNF(hh[i])*MAXSTEPSIZECCPSEARCH;
      }
   }
   return;
}
#if PARALLELISEDTK
/* ************************************************************************************** */
solreal GaussWaveFunction::evalELF(solreal x,solreal y,solreal z)
{
   static const solreal mto3=-10.0e0/3.0e0;
   static const solreal ooferm2=0.121300564999911e0;
   static const solreal eps=EPSFORELFVALUE;
   solreal nabx,naby,nabz,xmr,ymr,zmr,rho,cc,rr,alp,chib;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   int i,j;
   indp=0;
   rho=0.0000000e0;
   solreal tgx,tgy,tgz,kej;
   nabx=naby=nabz=kej=0.000000000000000e0;
#pragma omp parallel for private(indr,chib,tgx,tgy,tgz,cc) \
firstprivate(j) lastprivate(i) reduction(+: rho,kej,nabx,naby,nabz)
   for (i=0; i<nPri; i++) {
      indp=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      tgx=tgy=tgz=0.0e0;
      for (j=0; j<nPri; j++) {
         cc=cab[indp+j];
         tgx+=(gx[j]*cc);
         tgy+=(gy[j]*cc);
         tgz+=(gz[j]*cc);
         chib+=(chi[j]*cc);
      }
      rho+=(chib*chi[i]);
      kej+=(tgx*gx[i]);
      kej+=(tgy*gy[i]);
      kej+=(tgz*gz[i]);
      nabx+=(chib*gx[i]);
      naby+=(chib*gy[i]);
      nabz+=(chib*gz[i]);
   }
   nabx*=nabx;
   nabx*=4.0e0;
   nabx+=(4.00000e0*naby*naby);
   nabx+=(4.00000e0*nabz*nabz);
   naby=0.125e0*nabx/rho;
   nabz=(0.5e0*kej-naby+eps);
   nabx=ooferm2*nabz*nabz*pow(rho,mto3);
   return 1.0e0/(1.0e0+nabx);
}
#else
/* ************************************************************************************** */
solreal GaussWaveFunction::evalELF(solreal x,solreal y,solreal z)
{
   static const solreal mto3=-10.0e0/3.0e0;
   static const solreal ooferm2=0.121300564999911e0;
   static const solreal eps=EPSFORELFVALUE;
   solreal nabx,naby,nabz,xmr,ymr,zmr,rho,cc,rr,alp,chib;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   indp=0;
   rho=0.0000000e0;
   solreal tgx,tgy,tgz,kej;
   nabx=naby=nabz=kej=0.000000000000000e0;
   for (int i=0; i<nPri; i++) {
      //indr=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      tgx=tgy=tgz=0.0e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indp++];
         tgx+=(gx[j]*cc);
         tgy+=(gy[j]*cc);
         tgz+=(gz[j]*cc);
         chib+=(chi[j]*cc);
      }
      rho+=(chib*chi[i]);
      kej+=(tgx*gx[i]);
      kej+=(tgy*gy[i]);
      kej+=(tgz*gz[i]);
      nabx+=(chib*gx[i]);
      naby+=(chib*gy[i]);
      nabz+=(chib*gz[i]);
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      chib=tgx=tgy=tgz=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         tgx+=(cc*gx[i]);
         tgy+=(cc*gy[i]);
         tgz+=(cc*gz[i]);
      }
      /*
      cc=occN[nMOr];
      rho+=cc*chib*chib;
      kej+=(cc*tgx*tgx);
      kej+=(cc*tgy*tgy);
      kej+=(cc*tgz*tgz);
      nabx+=(cc*tgx*chib);
      naby+=(cc*tgy*chib);
      nabz+=(cc*tgz*chib);
      // */
      rho+=chib;
      kej+=(tgx);
      kej+=(tgy);
      kej+=(tgz);
      nabx+=(0.25e0*tgx);
      naby+=(0.25e0*tgy);
      nabz+=(0.25e0*tgz);
   }
   nabx*=nabx;
   nabx*=4.0e0;
   nabx+=(4.00000e0*naby*naby);
   nabx+=(4.00000e0*nabz*nabz);
   naby=0.125e0*nabx/rho;
   nabz=(0.5e0*kej-naby+eps);
   solreal dodh=ooferm2*nabz*nabz*pow(rho,mto3);
   return 1.0e0/(1.0e0+dodh);
}
#endif
/* ************************************************************************************** */
solreal GaussWaveFunction::evalShannonEntropy(solreal x,solreal y,solreal z)
{
   solreal rho=evalDensity(x,y,z);
   return (-rho*log(rho));
}
/* ************************************************************************************** */
solreal GaussWaveFunction::evalMomentumShannonEntropy(solreal px,solreal py,solreal pz)
{
   solreal ppi=evalFTDensity(px,py,pz);
   return (-ppi*log(ppi));
}
/* ************************************************************************************** */
solreal GaussWaveFunction::evalMagGradRho(solreal x,solreal y,solreal z)
{
   solreal gx,gy,gz,rho;
   evalRhoGradRho(x,y,z,rho,gx,gy,gz);
   rho=gx*gx+gy*gy+gz*gz;
   return sqrt(rho);
}
/* ************************************************************************************** */
solreal GaussWaveFunction::evalMagGradLOL(solreal x,solreal y,solreal z)
{
   solreal lol,xx[3],gl[3],hl[3][3];
   xx[0]=x; xx[1]=y; xx[2]=z;
   evalHessLOL(xx,lol,gl,hl);
   return sqrt(gl[0]*gl[0]+gl[1]*gl[1]+gl[2]*gl[2]);
}
/* ************************************************************************************** */
#if PARALLELISEDTK
solreal GaussWaveFunction::evalLOL(solreal x,solreal y,solreal z)
{
   static const solreal fo3=5.0e0/3.0e0;
   static const solreal tferm=5.742468000376382e0;
   static const solreal eps=EPSFORLOLVALUE;
   solreal xmr,ymr,zmr,rho,cc,rr,alp,chib;
   int indp,indr,ppt,i,j;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   indp=0;
   rho=0.0000000e0;
   solreal gxj,gyj,gzj,kej;
   omp_set_num_threads(PARALLELISEDTK);
   kej=0.000000000000000e0;
#pragma omp parallel for private(chib,gxj,gyj,gzj,indp,cc) \
firstprivate(j) lastprivate(i) reduction(+: kej,rho)
   for (i=0; i<nPri; i++) {
      //indr=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      gxj=gyj=gzj=0.0e0;
      indp=i*(nPri);
      for (j=0; j<nPri; j++) {
         cc=cab[indp+j];
         gxj+=(cc*gx[j]);
         gyj+=(cc*gy[j]);
         gzj+=(cc*gz[j]);
         chib+=(cc*chi[j]);
      }
      rho+=(chib*chi[i]);
      kej+=(gxj*gx[i]);
      kej+=(gyj*gy[i]);
      kej+=(gzj*gz[i]);
   }
   kej+=eps;
   solreal tau=tferm*pow(rho,fo3)/kej;
   return tau/(1.0e0+tau);
}
#else
/* ************************************************************************************** */
solreal GaussWaveFunction::evalLOL(solreal x,solreal y,solreal z)
{
   static const solreal fo3=5.0e0/3.0e0;
   static const solreal tferm=5.742468000376382e0;
   static const solreal eps=EPSFORLOLVALUE;
   solreal xmr,ymr,zmr,rho,cc,rr,alp,chib;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   indp=0;
   rho=0.0000000e0;
   solreal gxj,gyj,gzj,kej;
   kej=0.000000000000000e0;
   for (int i=0; i<nPri; i++) {
      //indr=i*(nPri);
      //cc=cab[indr];
      chib=0.0000000e0;
      gxj=gyj=gzj=0.0e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indp++];
         gxj+=(cc*gx[j]);
         gyj+=(cc*gy[j]);
         gzj+=(cc*gz[j]);
         chib+=(cc*chi[j]);
      }
      rho+=(chib*chi[i]);
      kej+=(gxj*gx[i]);
      kej+=(gyj*gy[i]);
      kej+=(gzj*gz[i]);
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      chib=gxj=gyj=gzj=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         gxj+=(cc*gx[i]);
         gyj+=(cc*gy[i]);
         gzj+=(cc*gz[i]);
      }
      /*
      cc=occN[nMOr];
      rho+=cc*chib*chib;
      kej+=(cc*gxj*gxj);
      kej+=(cc*gyj*gyj);
      kej+=(cc*gzj*gzj);
      // */
      rho+=chib;
      kej+=(gxj);
      kej+=(gyj);
      kej+=(gzj);
   }
   kej+=eps;
   solreal tau=tferm*pow(rho,fo3)/kej;
   return tau/(1.0e0+tau);
}
#endif
//end if PARALLELISEDTK
/* ************************************************************************************** */
#if PARALLELISEDTK
solreal GaussWaveFunction::evalKineticEnergyG(solreal x, solreal y, solreal z)
{
   solreal nabx,naby,nabz,xmr,ymr,zmr,cc,rr,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         //chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         //chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   solreal gxj,gyj,gzj;
   nabx=naby=nabz=0.000000000000000e0;
   indp=0;
   int i,j;
#pragma omp parallel for private(indp,cc,gxj,gyj,gzj) \
firstprivate(j) lastprivate(i) reduction(+: nabx,naby,nabz)
   for (i=0; i<nPri; i++) {
      indp=i*(nPri);
      cc=cab[indp+i];
      nabx+=(cc*gx[i]*gx[i]);
      naby+=(cc*gy[i]*gy[i]);
      nabz+=(cc*gz[i]*gz[i]);
      gxj=gyj=gzj=0.0000000e0;
      for (j=(i+1); j<nPri; j++) {
         cc=cab[indp+j];
         gxj+=(cc*gx[j]);
         gyj+=(cc*gy[j]);
         gzj+=(cc*gz[j]);
      }
      nabx+=(2.0e0*gx[i]*gxj);
      naby+=(2.0e0*gy[i]*gyj);
      nabz+=(2.0e0*gz[i]*gzj);
   }
   gxj=nabx;
   gxj+=naby;
   gxj+=nabz;
   gxj*=0.50e0;
   return gxj;
}
#else
/* ************************************************************************************** */
solreal GaussWaveFunction::evalKineticEnergyG(solreal x, solreal y, solreal z)
{
   solreal nabx,naby,nabz,xmr,ymr,zmr,cc,rr,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         //chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         //chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   solreal gxj,gyj,gzj;
   nabx=naby=nabz=0.000000000000000e0;
   indp=0;
   for (int i=0; i<nPri; i++) {
      indp=i*(nPri+1);
      cc=cab[indp++];
      nabx+=(cc*gx[i]*gx[i]);
      naby+=(cc*gy[i]*gy[i]);
      nabz+=(cc*gz[i]*gz[i]);
      gxj=gyj=gzj=0.0000000e0;
      for (int j=i+1; j<nPri; j++) {
         cc=cab[indp++];
         gxj+=(cc*gx[j]);
         gyj+=(cc*gy[j]);
         gzj+=(cc*gz[j]);
      }
      nabx+=(2.0e0*gx[i]*gxj);
      naby+=(2.0e0*gy[i]*gyj);
      nabz+=(2.0e0*gz[i]*gzj);
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      gxj=gyj=gzj=0.0e0;
      /*
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         gxj+=(cc*gx[i]);
         gyj+=(cc*gy[i]);
         gzj+=(cc*gz[i]);
      }
      cc=occN[nMOr];
      nabx+=(cc*gxj*gxj);
      naby+=(cc*gyj*gyj);
      nabz+=(cc*gzj*gzj);
      // */
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         gxj+=(cc*gx[i]);
         gyj+=(cc*gy[i]);
         gzj+=(cc*gz[i]);
      }
      nabx+=(gxj*gxj);
      naby+=(gyj*gyj);
      nabz+=(gzj*gzj);
   }
   gxj=nabx;
   gxj+=naby;
   gxj+=nabz;
   gxj*=0.50e0;
   return gxj;
}
#endif
/* ************************************************************************************** */
#if PARALLELISEDTK
void GaussWaveFunction::evalNabPhi2(solreal const x,solreal const y,solreal const z,\
      solreal &rho2ret,solreal &twoG)
{
#if DEBUG
   bool prntonce=true;
   if ( prntonce ) {
      cout << "Using parallel version (OpenMP)." << endl;
      prntonce=false;
   }
#endif
   solreal nabx,naby,nabz,xmr,ymr,zmr,cc,rr,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         ++indp;
      }
   }
   solreal gxj,gyj,gzj,rho,chib;
   nabx=naby=nabz=rho=0.000000000000000e0;
   indp=0;
   int i,j;
#pragma omp parallel for private(indp,cc,gxj,gyj,gzj,chib) \
firstprivate(j) lastprivate(i) reduction(+: nabx,naby,nabz,rho)
   for (i=0; i<nPri; i++) {
      indp=i*(nPri);
      cc=cab[indp+i];
      nabx+=(cc*gx[i]*gx[i]);
      naby+=(cc*gy[i]*gy[i]);
      nabz+=(cc*gz[i]*gz[i]);
      rho+=(cc*chi[i]*chi[i]);
      gxj=gyj=gzj=chib=0.0000000e0;
      for (j=(i+1); j<nPri; j++) {
         cc=cab[indp+j];
         chib+=cc*chi[j];
         gxj+=(cc*gx[j]);
         gyj+=(cc*gy[j]);
         gzj+=(cc*gz[j]);
      }
      rho+=(2.0e0*chi[i]*chib);
      nabx+=(2.0e0*gx[i]*gxj);
      naby+=(2.0e0*gy[i]*gyj);
      nabz+=(2.0e0*gz[i]*gzj);
   }
   twoG=nabx+naby+nabz;
   rho2ret=rho;
}
#else
/* ************************************************************************************** */
void GaussWaveFunction::evalNabPhi2(solreal const x,solreal const y,solreal const z,\
      solreal &rho2ret,solreal &twoG)
{
   solreal nabx,naby,nabz,xmr,ymr,zmr,cc,rr,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         ++indp;
      }
   }
   solreal gxj,gyj,gzj,chib,rho;
   nabx=naby=nabz=rho=0.000000000000000e0;
   indp=0;
   for (int i=0; i<nPri; i++) {
      indp=i*(nPri+1);
      cc=cab[indp++];
      nabx+=(cc*gx[i]*gx[i]);
      naby+=(cc*gy[i]*gy[i]);
      nabz+=(cc*gz[i]*gz[i]);
      rho+=(cc*chi[i]*chi[i]);
      gxj=gyj=gzj=chib=0.0000000e0;
      for (int j=i+1; j<nPri; j++) {
         cc=cab[indp++];
         chib+=(cc*chi[j]);
         gxj+=(cc*gx[j]);
         gyj+=(cc*gy[j]);
         gzj+=(cc*gz[j]);
      }
      nabx+=(2.0e0*gx[i]*gxj);
      naby+=(2.0e0*gy[i]*gyj);
      nabz+=(2.0e0*gz[i]*gzj);
      rho+=(2.0e0*chi[i]*chib);
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=primExp[i];
         cc=exp(alp*rr);
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[i],gy[i],gz[i]);
         gx[i]*=cc;
         gy[i]*=cc;
         gz[i]*=cc;
      }
      gxj=gyj=gzj=chib=0.0e0;
      /*
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         gxj+=(cc*gx[i]);
         gyj+=(cc*gy[i]);
         gzj+=(cc*gz[i]);
      }
      cc=occN[nMOr];
      rho+=(cc*chib*chib);
      nabx+=(cc*gxj*gxj);
      naby+=(cc*gyj*gyj);
      nabz+=(cc*gzj*gzj);
      // */
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         chib+=(cc*chi[i]);
         gxj+=(cc*gx[i]);
         gyj+=(cc*gy[i]);
         gzj+=(cc*gz[i]);
      }
      rho2ret+=chib;
      nabx+=(gxj);
      naby+=(gyj);
      nabz+=(gzj);
}
   rho2ret=rho;
   twoG=nabx+naby+nabz;
}
#endif
/* ************************************************************************************** */
#if PARALLELISEDTK
solreal GaussWaveFunction::evalKineticEnergyK(solreal x, solreal y, solreal z)
{
   solreal lap,xmr,ymr,zmr,cc,rr,alp;
   solreal sxx;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(-alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         hxx[indp]=evalLapAngCases(ppt,alp,xmr,ymr,zmr,rr);
         hxx[indp]*=cc;
         indp++;
      }
   }
   lap=0.000000000000000e0;
   indr=0;
   int i,j;
#pragma omp parallel for private(sxx,indr) \
firstprivate(j) lastprivate(i) reduction(+: lap)
   for (i=0; i<nPri; i++) {
      sxx=0.00000e0;
      indr=i*nPri;
      for (j=0; j<nPri; j++) {
         sxx+=cab[indr+j]*hxx[j];
      }
      lap+=chi[i]*sxx;
   }
   return (-0.50000e0*lap);
}
#else
/* ************************************************************************************** */
solreal GaussWaveFunction::evalKineticEnergyK(solreal x, solreal y, solreal z)
{
   solreal lap,xmr,ymr,zmr,cc,rr,alp;
   solreal sxx;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(-alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         hxx[indp]=evalLapAngCases(ppt,alp,xmr,ymr,zmr,rr);
         hxx[indp]*=cc;
         indp++;
      }
   }
   lap=0.000000000000000e0;
   indr=0;
   for (int i=0; i<nPri; i++) {
      sxx=0.00000e0;
      for (int j=0; j<nPri; j++) {
         //cc=cab[indr++];
         sxx+=cab[indr++]*hxx[j];
      }
      lap+=chi[i]*sxx;
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         ppt=primType[i];
         alp=0.5e0*primExp[i];
         cc=exp(-alp*rr);
         chi[i]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[i]*=cc;
         hxx[i]=evalLapAngCases(ppt,alp,xmr,ymr,zmr,rr);
         hxx[i]*=cc;
      }
      sxx=0.0e0;
      for (int i=nPri; i<totPri; ++i) {
         cc=EDFCoeff[i-nPri];
         sxx+=(cc*hxx[i]*chi[i]);
      }
      lap+=(sxx);
   }
   return (-0.50000e0*lap);
}
#endif
/* ************************************************************************************** */
void GaussWaveFunction::evald4SingCartA(int &ang,solreal &t,solreal &f,solreal &x,solreal &x2,
                                  solreal &d0,solreal &d1,solreal &d2,solreal &d3,solreal &d4)
{
   solreal fax2=f*x2;
   switch (ang) {
      case 0:
         d0=1.00000e0;
         d1=t*x;
         d2=t+fax2;
         d3=f*x*(3.00000e0+t*x2);
         d4=f*(3.00000e0+x2*(6.00000e0*t+fax2));
         break;
      case 1:
         d0=x;
         d1=1.00000e0+t*x2;
         d2=x*(3.00000e0*t+fax2);
         d3=3.00000e0*t+fax2*(6.00000e0+t*x2);
         d4=f*x*(15.0000e0+x2*(10.0000e0*t+fax2));
         break;
      case 2:
         d0=x2;
         d1=x*(2.00000e0+t*x2);
         d2=2.00000e0+x2*(5.00000e0*t+fax2);
         d3=x*(12.0000e0*t+fax2*(9.00000e0+t*x2));
         d4=12.0000e0*t+fax2*(39.0000e0+x2*(14.0000e0*t+fax2));
         break;
      case 3:
         d0=x*x2;
         d1=x2*(3.00000e0+t*x2);
         d2=x*(6.00000e0+x2*(7.00000e0*t+fax2));
         d3=6.00000e0-x2*((9.00000e0*t-fax2)*(3.00000e0+t*x2));
         d4=x*(60.0000e0*t+fax2*(75.0000e0+x2*(18.0000e0*t+fax2)));
         break;
      case 4 :
         d0=x2*x2;
         d1=x2*x*(4.00000e0+t*x2);
         d2=x2*(12.0000e0+x2*(9.00000e0*t+fax2));
         d3=x*(24.0000e0+t*x2*(48.0000e0+x2*(15.0000e0*t+fax2)));
         d4=24.0e0+x2*(168.0e0*t+fax2*(123.0e0+t*x2*(22.0e0+t*x2)));
         break;
      case 5 :
         d0=x2*x2*x;
         d1=x2*x2*(5.00000e0+t*x2);
         d2=x2*x*(20.0000e0+x2*(11.0000e0*t+fax2));
         d3=x2*(60.0000e0+x2*(75.0000e0*t+fax2*(18.0000e0+t*x2)));
         d4=x*(120.0e0+x2*(360.0e0*t+fax2*(183.0e0+x2*(26.0e0*t+fax2))));
         break;
      default:
         d0=0.00000e0;
         d1=0.00000e0;
         d2=0.00000e0;
         d3=0.00000e0;
         d4=0.00000e0;
         break;
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evald4Ang(int (&a)[3],solreal &alp,solreal (&x)[3],solreal (&x2)[3],
                              solreal (&d0)[3],solreal (&d1)[3],solreal (&d2)[3],solreal (&d3)[3],solreal (&d4)[3])
{
   solreal tma,fa2;
   tma=-2.00000e0*alp;
   fa2=tma*tma;
   for (int i=0; i<3; i++) {
      evald4SingCartA(a[i],tma,fa2,x[i],x2[i],d0[i],d1[i],d2[i],d3[i],d4[i]);
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evald1SingCartA(int &ang,solreal &t,solreal x,\
      solreal &d0,solreal &d1)
{
   solreal x2=x*x;
   switch ( ang ) {
      case 0 :
         d0=1.0e0;
         d1=t*x;
         break;
      case 1 :
         d0=x;
         d1=1.0e0+t*x2;
         break;
      case 2 :
         d0=x2;
         d1=x*(2.00000e0+t*x2);
         break;
      case 3 :
         d0=x*x2;
         d1=x2*(3.00000e0+t*x2);
         break;
      case 4 :
         d0=x2*x2;
         d1=x2*x*(4.00000e0+t*x2);
         break;
      case 5 :
         d0=x2*x2*x;
         d1=x2*x2*(5.00000e0+t*x2);
         break;
      default :
         displayWarningMessage("Non supported primitive type. Numerical errors"
               "are expected!");
         d0=d1=0.0e0;
         break;
   }
}
/* ************************************************************************************** */
void GaussWaveFunction::evald2SingCartA(int &ang,solreal &t,solreal x,\
                      solreal &d0,solreal &d1,solreal &d2)
{
   solreal x2=x*x;
   solreal f=t*t;
   solreal fax2=f*x2;
   switch (ang) {
      case 0:
         d0=1.00000e0;
         d1=t*x;
         d2=t+fax2;
         break;
      case 1:
         d0=x;
         d1=1.00000e0+t*x2;
         d2=x*(3.00000e0*t+fax2);
         break;
      case 2:
         d0=x2;
         d1=x*(2.00000e0+t*x2);
         d2=2.00000e0+x2*(5.00000e0*t+fax2);
         break;
      case 3:
         d0=x*x2;
         d1=x2*(3.00000e0+t*x2);
         d2=x*(6.00000e0+x2*(7.00000e0*t+fax2));
         break;
      case 4 :
         d0=x2*x2;
         d1=x2*x*(4.00000e0+t*x2);
         d2=x2*(12.0000e0+x2*(9.00000e0*t+fax2));
         break;
      case 5 :
         d0=x2*x2*x;
         d1=x2*x2*(5.00000e0+t*x2);
         d2=x2*x*(20.0000e0+x2*(11.0000e0*t+fax2));
         break;
      default:
         displayWarningMessage("Non supported primitive type. Numerical errors"
               "are expected!");
         d0=0.00000e0;
         d1=0.00000e0;
         d2=0.00000e0;
         break;
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evald3SingCartA(int &ang,solreal &t,solreal &f,solreal &x,solreal &x2,
                                    solreal &d0,solreal &d1,solreal &d2,solreal &d3)
{
   solreal fax2=f*x2;
   switch (ang) {
      case 0:
         d0=1.00000e0;
         d1=t*x;
         d2=t+fax2;
         d3=f*x*(3.00000e0+t*x2);
         break;
      case 1:
         d0=x;
         d1=1.00000e0+t*x2;
         d2=x*(3.00000e0*t+fax2);
         d3=3.00000e0*t+fax2*(6.00000e0+t*x2);
         break;
      case 2:
         d0=x2;
         d1=x*(2.00000e0+t*x2);
         d2=2.00000e0+x2*(5.00000e0*t+fax2);
         d3=x*(12.0000e0*t+fax2*(9.00000e0+t*x2));
         break;
      case 3:
         d0=x*x2;
         d1=x2*(3.00000e0+t*x2);
         d2=x*(6.00000e0+x2*(7.00000e0*t+fax2));
         d3=6.00000e0-x2*((9.00000e0*t-fax2)*(3.00000e0+t*x2));
         break;
      case 4 :
         d0=x2*x2;
         d1=x2*x*(4.00000e0+t*x2);
         d2=x2*(12.0000e0+x2*(9.00000e0*t+fax2));
         d3=x*(24.0000e0+t*x2*(48.0000e0+x2*(15.0000e0*t+fax2)));
         break;
      case 5 :
         d0=x2*x2*x;
         d1=x2*x2*(5.00000e0+t*x2);
         d2=x2*x*(20.0000e0+x2*(11.0000e0*t+fax2));
         d3=x2*(60.0000e0+x2*(75.0000e0*t+fax2*(18.0000e0+t*x2)));
         break;
      default:
         displayWarningMessage("Non supported primitive type. Numerical errors"
               "are expected!");
         d0=0.00000e0;
         d1=0.00000e0;
         d2=0.00000e0;
         d3=0.00000e0;
         break;
   }
   return;
}
/* ************************************************************************************** */
void GaussWaveFunction::evald3Ang(int (&a)[3],solreal &alp,solreal (&x)[3],solreal (&x2)[3],
                              solreal (&d0)[3],solreal (&d1)[3],solreal (&d2)[3],solreal (&d3)[3])
{
   solreal tma,fa2;
   tma=-2.00000e0*alp;
   fa2=tma*tma;
   for (int i=0; i<3; i++) {
      evald3SingCartA(a[i],tma,fa2,x[i],x2[i],d0[i],d1[i],d2[i],d3[i]);
   }
   return;
}
/* ************************************************************************************** */
/*
solreal GaussWaveFunction::evalLapRhoUsingd2(solreal x, solreal y, solreal z)
{
   solreal lap,cc,rr,alp;
   solreal sxx,gxs,gys,gzs;
   solreal r[3];
   static solreal r2[3],D0[3],D1[3],D2[3],D3[3],D4[3];
   static int aa[3];
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      rr=0.00000e0;
      r[0]=x;
      r[1]=y;
      r[2]=z;
      for (int m=0; m<3; m++) {
         r[m]-=R[indr++];
         r2[m]=(r[m]*r[m]);
         rr+=r2[m];
      }
      for (int j=0; j<myPN[i]; j++) {
         ppt=3*primType[indp];
         for (int m=0; m<3; m++) {aa[m]=prTy[ppt++];}
         alp=primExp[indp];
         cc=exp(-alp*rr);
         evald4Ang(aa,alp,r,r2,D0,D1,D2,D3,D4);
         chi[indp]=D0[0]*D0[1]*D0[2];
         chi[indp]*=cc;
         hxx[indp]=(D2[0]*D0[1]*D0[2])+(D0[0]*D2[1]*D0[2])+(D0[0]*D0[1]*D2[2]);
         hxx[indp]*=cc;
         gx[indp]=(D1[0]*D0[1]*D0[2]);
         gy[indp]=(D0[0]*D1[1]*D0[2]);
         gz[indp]=(D0[0]*D0[1]*D1[2]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         indp++;
      }
   }
   lap=0.000000000000000e0;
   indr=0;
   for (int i=0; i<nPri; i++) {
      sxx=0.00000e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indr++];
         sxx+=cc*hxx[j];
      }
      lap+=chi[i]*sxx;
      indp=i*(nPri+1);
      cc=cab[indp++];
      lap+=cc*gx[i]*gx[i];
      lap+=cc*gy[i]*gy[i];
      lap+=cc*gz[i]*gz[i];
      gxs=gys=gzs=0.00000e0;
      for (int k=(i+1); k<nPri; k++) {
         cc=cab[indp++];
         gxs+=cc*gx[k];
         gys+=cc*gy[k];
         gzs+=cc*gz[k];
      }
      lap+=gx[i]*gxs*2.00000e0;
      lap+=gy[i]*gys*2.00000e0;
      lap+=gz[i]*gzs*2.00000e0;
   }
   return (2.00000e0*lap);
}
// */
/* *************************************************************************************** */
void GaussWaveFunction::evalDiDjDkChi(int &pty,solreal &alp,solreal x,solreal y,solreal z,
                                  solreal (&dlm)[3][3],solreal (&dijk)[3][3][3])
{
   static int aa[3];
   int ppt=3*pty;
   for (int m=0; m<3; m++) {aa[m]=prTy[ppt++];}
   static solreal X[3],X2[3],D0[3],D1[3],D2[3],D3[3];
   X[0]=x;
   X[1]=y;
   X[2]=z;
   for (int i=0; i<3; i++) {X2[i]=X[i]*X[i];}
   evald3Ang(aa,alp,X,X2,D0,D1,D2,D3);
   //------------------------------
   dijk[0][0][0]=D3[0]*D0[1]*D0[2];
   dijk[0][0][1]=D2[0]*D1[1]*D0[2];
   dijk[0][0][2]=D2[0]*D0[1]*D1[2];
   //
   dijk[0][1][0]=dijk[0][0][1];
   dijk[0][1][1]=D1[0]*D2[1]*D0[2];
   dijk[0][1][2]=D1[0]*D1[1]*D1[2];
   //
   dijk[0][2][0]=dijk[0][0][2];
   dijk[0][2][1]=dijk[0][1][2];
   dijk[0][2][2]=D1[0]*D0[1]*D2[2];
   //------------------------------
   dijk[1][0][0]=dijk[0][0][1];
   dijk[1][0][1]=dijk[0][1][1];
   dijk[1][0][2]=dijk[0][1][2];
   //
   dijk[1][1][0]=dijk[0][1][1];
   dijk[1][1][1]=D0[0]*D3[1]*D0[2];
   dijk[1][1][2]=D0[0]*D2[1]*D1[2];
   //
   dijk[1][2][0]=dijk[0][1][2];
   dijk[1][2][1]=dijk[1][1][2];
   dijk[1][2][2]=D0[0]*D1[1]*D2[2];
   //------------------------------
   dijk[2][0][0]=dijk[0][0][2];
   dijk[2][0][1]=dijk[0][1][2];
   dijk[2][0][2]=dijk[0][2][2];
   //
   dijk[2][1][0]=dijk[0][1][2];
   dijk[2][1][1]=dijk[1][1][2];
   dijk[2][1][2]=dijk[1][2][2];
   //
   dijk[2][2][0]=dijk[0][2][2];
   dijk[2][2][1]=dijk[1][2][2];
   dijk[2][2][2]=D0[0]*D0[1]*D3[2];
   //------------------------------
   //++++++++++++++++++++++++++++++
   //------------------------------
   dlm[0][0]=D2[0]*D0[1]*D0[2];
   dlm[0][1]=D1[0]*D1[1]*D0[2];
   dlm[0][2]=D1[0]*D0[1]*D1[2];
   //------------------------------
   dlm[1][0]=dlm[0][1];
   dlm[1][1]=D0[0]*D2[1]*D0[2];
   dlm[1][2]=D0[0]*D1[1]*D1[2];
   //------------------------------
   dlm[2][0]=dlm[0][2];
   dlm[2][1]=dlm[1][2];
   dlm[2][2]=D0[0]*D0[1]*D2[2];
   //------------------------------
   solreal expterm=exp(-alp*(X2[0]+X2[1]+X2[2]));
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         dlm[i][j]*=expterm;
         for (int k=0; k<3; k++) {dijk[i][j][k]*=expterm;}
      }
   }
}
#if (PARALLELISEDTK && DEBUG)
/* ************************************************************************************** */
void GaussWaveFunction::evalHessLOL(solreal x, solreal y, solreal z, solreal &dens, solreal &keG, solreal &lol,
                                solreal &ddx, solreal &ddy, solreal &ddz,
                                solreal &dxx, solreal &dyy, solreal &dzz,
                                solreal &dxy, solreal &dxz, solreal &dyz)
{
   solreal xmr,ymr,zmr,cc,rr,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[indp],hyy[indp],hzz[indp],
                          hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   static solreal dG[3],dR[3],ds[3],ddG[3][3],ddR[3][3],dds[3][3],ddC[3][3],dddC[3][3][3];
   for (int i=0; i<3; i++) {
      dR[i]=0.0e0;
      ds[i]=0.0e0;
      for (int j=0; j<3; j++) {
         ddR[i][j]=0.0e0;
         dds[i][j]=0.0e0;
         ddC[i][j]=0.0e0;
         for (int k=0; k<3; k++) {dddC[i][j][k]=0.0e0;}
      }
   }
   solreal dGx,dGy,dGz,ddGxx,ddGyy,ddGzz,ddGxy,ddGxz,ddGyz;
   solreal rho,G,s,nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,delx,dely,delz;
   rho=G=s=dGx=dGy=dGz=0.0e0;
   ddGxx=ddGyy=ddGzz=ddGxy=ddGxz=ddGyz=0.0e0;
   nabxx=nabyy=nabzz=nabxy=nabxz=nabyz=0.000000000000000e0;
   delx=dely=delz=0.0e0;
   solreal gxs,gys,gzs,sxx,syy,szz,sxy,sxz,syz,chij;
   indp=0;
   int i,j;
#pragma omp parallel for private(ppt,alp,indr,indp,xmr,ymr,zmr,gxs,gys,gzs,\
                                 sxx,syy,szz,sxy,sxz,syz,chij,cc) \
firstprivate(j) lastprivate(i) reduction(+: nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,\
                                            delx,dely,delz,rho,G,dGx,dGy,dGz,\
                                            ddGxx,ddGyy,ddGzz,ddGxy,ddGxz,ddGyz)
   for (i=0; i<nPri; i++) {
      ppt=primType[i];
      alp=primExp[i];
      indr=3*(primCent[i]);
      xmr=x-R[indr];
      ymr=y-R[indr+1];
      zmr=z-R[indr+2];
      evalDiDjDkChi(ppt,alp,xmr,ymr,zmr,ddC,dddC);
      indp=i*nPri;
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chij=0.00000e0;
      for (j=0; j<nPri; j++) {
         cc=cab[indp+j];
         chij+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
         sxx+=cc*hxx[j];
         syy+=cc*hyy[j];
         szz+=cc*hzz[j];
         sxy+=cc*hxy[j];
         sxz+=cc*hxz[j];
         syz+=cc*hyz[j];
      }
      nabxx+=(chij*hxx[i]+gx[i]*gxs);
      nabyy+=(chij*hyy[i]+gy[i]*gys);
      nabzz+=(chij*hzz[i]+gz[i]*gzs);
      nabxy+=(chij*hxy[i]+gx[i]*gys);
      nabxz+=(chij*hxz[i]+gx[i]*gzs);
      nabyz+=(chij*hyz[i]+gy[i]*gzs);
      delx+=(chij*gx[i]);
      dely+=(chij*gy[i]);
      delz+=(chij*gz[i]);
      rho+=(chij*chi[i]);
      //
      G+=(gxs*gx[i]+gys*gy[i]+gzs*gz[i]);
      dGx+=(gxs*hxx[i]+gys*hxy[i]+gzs*hxz[i]);
      dGy+=(gxs*hxy[i]+gys*hyy[i]+gzs*hyz[i]);
      dGz+=(gxs*hxz[i]+gys*hyz[i]+gzs*hzz[i]);
      //
      ddGxx+=(hxx[i]*sxx+hxy[i]*sxy+hxz[i]*sxz);
      ddGxy+=(hxx[i]*sxy+hxy[i]*syy+hxz[i]*syz);
      ddGxz+=(hxx[i]*sxz+hxy[i]*syz+hxz[i]*szz);
      ddGyy+=(hxy[i]*sxy+hyy[i]*syy+hyz[i]*syz);
      ddGyz+=(hxy[i]*sxz+hyy[i]*syz+hyz[i]*szz);
      ddGzz+=(hxz[i]*sxz+hyz[i]*syz+hzz[i]*szz);
      //
      ddGxx+=(dddC[0][0][0]*gxs+dddC[0][0][1]*gys+dddC[0][0][2]*gzs);
      ddGxy+=(dddC[0][1][0]*gxs+dddC[0][1][1]*gys+dddC[0][1][2]*gzs);
      ddGxz+=(dddC[0][2][0]*gxs+dddC[0][2][1]*gys+dddC[0][2][2]*gzs);
      ddGyy+=(dddC[1][1][0]*gxs+dddC[1][1][1]*gys+dddC[1][1][2]*gzs);
      ddGyz+=(dddC[1][2][0]*gxs+dddC[1][2][1]*gys+dddC[1][2][2]*gzs);
      ddGzz+=(dddC[2][2][0]*gxs+dddC[2][2][1]*gys+dddC[2][2][2]*gzs);
   }
   //cout << "ddC test: " << ddC[0][0] << " " << ddC[1][1] << endl
   //     << "          " << ddC[2][2] << " " << ddC[0][1] << endl
   //     << "          " << ddC[0][2] << " " << ddC[1][2] << endl;
   ddG[0][0]=ddGxx;
   ddG[0][1]=ddGxy;
   ddG[0][2]=ddGxz;
   ddG[1][0]=ddGxy;
   ddG[1][1]=ddGyy;
   ddG[1][2]=ddGyz;
   ddG[2][0]=ddGxz;
   ddG[2][1]=ddGyz;
   ddG[2][2]=ddGzz;
   ddR[0][0]=2.00000e0*nabxx;
   ddR[1][1]=2.00000e0*nabyy;
   ddR[2][2]=2.00000e0*nabzz;
   ddR[0][1]=2.00000e0*nabxy;
   ddR[0][2]=2.00000e0*nabxz;
   ddR[1][2]=2.00000e0*nabyz;
   ddR[1][0]=ddR[0][1];
   ddR[2][0]=ddR[0][2];
   ddR[2][1]=ddR[1][2];
   dR[0]=2.0e0*delx;
   dR[1]=2.0e0*dely;
   dR[2]=2.0e0*delz;
   dG[0]=dGx; dG[1]=dGy; dG[2]=dGz;
   for (i=0; i<3; i++) {
      dG[i]*=2.0e0;
      dG[i]+=(SIGNF(dG[i])*EPSFORLOLVALUE*0.01);
      for (j=0; j<3; j++) {
         ddG[i][j]*=2.0e0;
         //cout << ddG[i][j] << " ";
         //cout << ddR[i][j] << " ";
         //ddG[i][j]+=(SIGNF(ddG[i][j])*EPSFORLOLVALUE*0.001);
      }
      //cout << endl;
   }
   //cout << endl;
   rho+=EPSFORLOLVALUE;
   G+=EPSFORLOLVALUE;
   //==============================================================
   solreal oo2ferm=0.17414115323489e0,fo3=5.0e0/3.0e0;
   solreal rhoto5o3=pow(rho,fo3),oorho=1.0e0/(rho);
   s=oo2ferm*(G)/rhoto5o3;
   for (i=0; i<3; i++) {
      ds[i]=oo2ferm*(dG[i]-fo3*oorho*G*dR[i])/rhoto5o3;
      for (j=0; j<3; j++) {
         dds[i][j]=oo2ferm*(ddG[i][j]
                            -fo3*oorho*(dG[i]*dR[j]+dG[j]*dR[i]+G*ddR[i][j])
                            +(40.0e0*G*oorho*oorho/9.0e0)*dR[i]*dR[j])/rhoto5o3;
         //cout << dds[i][j] << " ";
      }
      //cout << endl;
   }
   //cout <<endl;
   //cout << "rho: " << rho << endl;
   //cout << "oor: " << oorho << endl;
   //cout << "s:   " << s << endl;
   //cout << ds[0] << " " << ds[1] << " " << ds[2] << endl;
   //cout << "G:   " << G << endl;
   //cout << "dG:  " << dG[0] << " " << dG[1] << " " << dG[2] << endl;
   //==============================================================
   dens=rho;
   keG=0.5e0*G;
   lol=1.0e0/(1.0e0+s);
   solreal gam2=lol*lol,gam3=lol*gam2;
   ddx=(-ds[0]*gam2);
   ddy=(-ds[1]*gam2);
   ddz=(-ds[2]*gam2);
   dxx=2.0e0*gam3*ds[0]*ds[0]-gam2*dds[0][0];
   dyy=2.0e0*gam3*ds[1]*ds[1]-gam2*dds[1][1];
   dzz=2.0e0*gam3*ds[2]*ds[2]-gam2*dds[2][2];
   dxy=2.0e0*gam3*ds[0]*ds[1]-gam2*dds[0][1];
   dxz=2.0e0*gam3*ds[0]*ds[2]-gam2*dds[0][2];
   dyz=2.0e0*gam3*ds[1]*ds[2]-gam2*dds[1][2];
   //cout << ds[0] << " " << ds[1] << " " << ds[2] << endl;
   //==============================================================
   /*
    ddx=dR[0];
    ddy=dR[1];
    ddz=dR[2];
    dxx=ddR[0][0];
    dyy=ddR[1][1];
    dzz=ddR[2][2];
    dxy=ddR[0][1];
    dxz=ddR[0][2];
    dyz=ddR[1][2];
    // */
   //==============================================================
   return;
}
#else
/* ************************************************************************************** */
void GaussWaveFunction::evalHessLOL(solreal x, solreal y, solreal z, solreal &dens, solreal &keG, solreal &lol,
                                solreal &ddx, solreal &ddy, solreal &ddz,
                                solreal &dxx, solreal &dyy, solreal &dzz,
                                solreal &dxy, solreal &dxz, solreal &dyz)
{
   solreal xmr,ymr,zmr,cc,rr,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,
                        gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
                          hxx[indp],hyy[indp],hzz[indp],
                          hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   static solreal dG[3],dR[3],ds[3],ddG[3][3],ddR[3][3],dds[3][3],ddC[3][3],dddC[3][3][3];
   for (int i=0; i<3; i++) {
      dG[i]=0.0e0;
      dR[i]=0.0e0;
      ds[i]=0.0e0;
      for (int j=0; j<3; j++) {
         ddR[i][j]=0.0e0;
         ddG[i][j]=0.0e0;
         dds[i][j]=0.0e0;
         ddC[i][j]=0.0e0;
         for (int k=0; k<3; k++) {dddC[i][j][k]=0.0e0;}
      }
   }
   solreal rho,G,s,nabxx,nabyy,nabzz,nabxy,nabxz,nabyz,delx,dely,delz;
   rho=G=s=0.0e0;
   nabxx=nabyy=nabzz=nabxy=nabxz=nabyz=0.000000000000000e0;
   delx=dely=delz=0.0e0;
   solreal gxs,gys,gzs,sxx,syy,szz,sxy,sxz,syz,chij;
   indp=0;
   for (int i=0; i<nPri; i++) {
      ppt=primType[i];
      alp=primExp[i];
      indr=3*(primCent[i]);
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      evalDiDjDkChi(ppt,alp,xmr,ymr,zmr,ddC,dddC);
      gxs=gys=gzs=0.00000e0;
      sxx=syy=szz=sxy=sxz=syz=0.00000e0;
      chij=0.00000e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indp++];
         chij+=cc*chi[j];
         gxs+=cc*gx[j];
         gys+=cc*gy[j];
         gzs+=cc*gz[j];
         sxx+=cc*hxx[j];
         syy+=cc*hyy[j];
         szz+=cc*hzz[j];
         sxy+=cc*hxy[j];
         sxz+=cc*hxz[j];
         syz+=cc*hyz[j];
      }
      nabxx+=chij*hxx[i]+gx[i]*gxs;
      nabyy+=chij*hyy[i]+gy[i]*gys;
      nabzz+=chij*hzz[i]+gz[i]*gzs;
      nabxy+=chij*hxy[i]+gx[i]*gys;
      nabxz+=chij*hxz[i]+gx[i]*gzs;
      nabyz+=chij*hyz[i]+gy[i]*gzs;
      delx+=chij*gx[i];
      dely+=chij*gy[i];
      delz+=chij*gz[i];
      rho+=chij*chi[i];
      //
      G+=(gxs*gx[i]+gys*gy[i]+gzs*gz[i]);
      dG[0]+=(gxs*hxx[i]+gys*hxy[i]+gzs*hxz[i]);
      dG[1]+=(gxs*hxy[i]+gys*hyy[i]+gzs*hyz[i]);
      dG[2]+=(gxs*hxz[i]+gys*hyz[i]+gzs*hzz[i]);
      //
      ddG[0][0]+=(hxx[i]*sxx+hxy[i]*sxy+hxz[i]*sxz);
      ddG[0][1]+=(hxx[i]*sxy+hxy[i]*syy+hxz[i]*syz);
      ddG[0][2]+=(hxx[i]*sxz+hxy[i]*syz+hxz[i]*szz);
      ddG[1][1]+=(hxy[i]*sxy+hyy[i]*syy+hyz[i]*syz);
      ddG[1][2]+=(hxy[i]*sxz+hyy[i]*syz+hyz[i]*szz);
      ddG[2][2]+=(hxz[i]*sxz+hyz[i]*syz+hzz[i]*szz);
      //
      ddG[0][0]+=(dddC[0][0][0]*gxs+dddC[0][0][1]*gys+dddC[0][0][2]*gzs);
      ddG[0][1]+=(dddC[0][1][0]*gxs+dddC[0][1][1]*gys+dddC[0][1][2]*gzs);
      ddG[0][2]+=(dddC[0][2][0]*gxs+dddC[0][2][1]*gys+dddC[0][2][2]*gzs);
      ddG[1][1]+=(dddC[1][1][0]*gxs+dddC[1][1][1]*gys+dddC[1][1][2]*gzs);
      ddG[1][2]+=(dddC[1][2][0]*gxs+dddC[1][2][1]*gys+dddC[1][2][2]*gzs);
      ddG[2][2]+=(dddC[2][2][0]*gxs+dddC[2][2][1]*gys+dddC[2][2][2]*gzs);
   }
   //cout << "ddC test: " << ddC[0][0] << " " << ddC[1][1] << endl
   //     << "          " << ddC[2][2] << " " << ddC[0][1] << endl
   //     << "          " << ddC[0][2] << " " << ddC[1][2] << endl;
   ddG[1][0]=ddG[0][1];
   ddG[2][0]=ddG[0][2];
   ddG[2][1]=ddG[1][2];
   ddR[0][0]=2.00000e0*nabxx;
   ddR[1][1]=2.00000e0*nabyy;
   ddR[2][2]=2.00000e0*nabzz;
   ddR[0][1]=2.00000e0*nabxy;
   ddR[0][2]=2.00000e0*nabxz;
   ddR[1][2]=2.00000e0*nabyz;
   ddR[1][0]=ddR[0][1];
   ddR[2][0]=ddR[0][2];
   ddR[2][1]=ddR[1][2];
   dR[0]=2.0e0*delx;
   dR[1]=2.0e0*dely;
   dR[2]=2.0e0*delz;
   for (int i=0; i<3; i++) {
      dG[i]*=2.0e0;
      dG[i]+=(SIGNF(dG[i])*EPSFORLOLVALUE*0.01);
      for (int j=0; j<3; j++) {
         ddG[i][j]*=2.0e0;
         //cout << ddG[i][j] << " ";
         //cout << ddR[i][j] << " ";
         //ddG[i][j]+=(SIGNF(ddG[i][j])*EPSFORLOLVALUE*0.001);
      }
      //cout << endl;
   }
   //cout << endl;
   rho+=EPSFORLOLVALUE;
   G+=EPSFORLOLVALUE;
   //==============================================================
   solreal oo2ferm=0.17414115323489e0,fo3=5.0e0/3.0e0;
   solreal rhoto5o3=pow(rho,fo3),oorho=1.0e0/(rho);
   s=oo2ferm*(G)/rhoto5o3;
   for (int i=0; i<3; i++) {
      ds[i]=oo2ferm*(dG[i]-fo3*oorho*G*dR[i])/rhoto5o3;
      for (int j=0; j<3; j++) {
         dds[i][j]=oo2ferm*(ddG[i][j]
                            -fo3*oorho*(dG[i]*dR[j]+dG[j]*dR[i]+G*ddR[i][j])
                            +(40.0e0*G*oorho*oorho/9.0e0)*dR[i]*dR[j])/rhoto5o3;
      }
   }
   dens=rho;
   keG=0.5e0*G;
   lol=1.0e0/(1.0e0+s);
   solreal gam2=lol*lol,gam3=lol*gam2;
   ddx=(-ds[0]*gam2);
   ddy=(-ds[1]*gam2);
   ddz=(-ds[2]*gam2);
   dxx=2.0e0*gam3*ds[0]*ds[0]-gam2*dds[0][0];
   dyy=2.0e0*gam3*ds[1]*ds[1]-gam2*dds[1][1];
   dzz=2.0e0*gam3*ds[2]*ds[2]-gam2*dds[2][2];
   dxy=2.0e0*gam3*ds[0]*ds[1]-gam2*dds[0][1];
   dxz=2.0e0*gam3*ds[0]*ds[2]-gam2*dds[0][2];
   dyz=2.0e0*gam3*ds[1]*ds[2]-gam2*dds[1][2];
   if ( ihaveEDF ) {
      cout << "Warning: HessLOL does not include EDF contributions!" << endl;
   }
   return;
}
#endif
/* *************************************************************************************** */
void GaussWaveFunction::evalHessLOL(solreal (&x)[3],solreal &lol,solreal (&glol)[3],solreal (&hlol)[3][3])
{
   static solreal rho,ke;
   evalHessLOL(x[0],x[1],x[2],rho,ke,lol,glol[0],glol[1],glol[2],
               hlol[0][0],hlol[1][1],hlol[2][2],hlol[0][1],hlol[0][2],hlol[1][2]);
   hlol[1][0]=hlol[0][1];
   hlol[2][0]=hlol[0][2];
   hlol[2][1]=hlol[1][2];
   return;
}
/* *************************************************************************************** */
void GaussWaveFunction::evalFTASingCartA(int &ang,solreal &a,solreal &ooa,solreal &osra,
                                     solreal &px,solreal &px2,solreal &Rx,
                                     solreal &RePhi,solreal &ImPhi)
{
   //static const solreal srpi=1.77245385090551602729817; //sqrt(pi)
   static const solreal srpi=0.707106781186547524405; //sqrt(pi/(2pi)) --This includes the
                                                      //normalization constant for the FT.
   solreal srpoa=srpi*osra; //sqrt(pi/alpha)
   solreal pr=px*Rx;
   //cout << "pr: " << pr << ", sin(pr): " << sin(pr) << endl;
   switch (ang) {
      case 0:
         RePhi=srpoa*cos(pr);
         ImPhi=(-srpoa*sin(pr));
         break;
      case 1:
         srpoa*=(0.5e0*ooa*px);
         RePhi=-srpoa*sin(pr);
         ImPhi=-srpoa*cos(pr);
         break;
      case 2:
         srpoa*=(0.25e0*ooa*ooa*(2.0e0*a-px2));
         RePhi=srpoa*cos(pr);
         ImPhi=-srpoa*sin(pr);
         break;
      case 3:
         srpoa*=(0.125e0*ooa*ooa*ooa*px*(px2-6.0e0*a));
         RePhi=srpoa*sin(pr);
         ImPhi=srpoa*cos(pr);
         break;
      default:
         RePhi=0.0e0;
         ImPhi=0.0e0;
         break;
   }
   return;
}
/* *************************************************************************************** */
void GaussWaveFunction::evalFTAng(int (&a)[3],solreal &alp,solreal &ooalp,solreal (&p)[3],solreal (&p2)[3],
                              solreal (&Rx)[3],complex<solreal> &pang)
{
   solreal osa,rp,ip;
   //oa=1.0e0/alp;
   osa=sqrt(ooalp);
   evalFTASingCartA(a[0],alp,ooalp,osa,p[0],p2[0],Rx[0],rp,ip);
   complex<solreal> tmp(rp,ip);
   evalFTASingCartA(a[1],alp,ooalp,osa,p[1],p2[1],Rx[1],rp,ip);
   tmp*=(complex<solreal>(rp,ip));
   evalFTASingCartA(a[2],alp,ooalp,osa,p[2],p2[2],Rx[2],rp,ip);
   tmp*=(complex<solreal>(rp,ip));
   pang=tmp;
   return;
}
/* *************************************************************************************** */
void GaussWaveFunction::evalFTChi(int &pty,solreal &alp,solreal (&Rx)[3],
                              solreal px,solreal py,solreal pz,
                              complex<solreal> &phi)
{
   static int aa[3];
   int ppt=3*pty;
   for (int m=0; m<3; m++) {aa[m]=prTy[ppt++];}
   solreal P[3],P2[3],ooalp;
   ooalp=1.0e0/alp;
   P[0]=px; P[1]=py; P[2]=pz;
   for (int i=0; i<3; i++) {P2[i]=P[i]*P[i];}
   complex<solreal> chit;
   evalFTAng(aa,alp,ooalp,P,P2,Rx,chit);
   solreal expterm=exp(-0.25*ooalp*(P2[0]+P2[1]+P2[2]));
   chit*=expterm;
   phi=chit;
   return;
}
#if PARALLELISEDTK
/* *************************************************************************************** */
solreal GaussWaveFunction::evalFTDensity(solreal px,solreal py,solreal pz)
{
   int indr,indp,ppt;
   solreal rhop,Rx[3],alp;
   complex<solreal> chit;
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      Rx[0]=R[indr++];
      Rx[1]=R[indr++];
      Rx[2]=R[indr++];
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         evalFTChi(ppt,alp,Rx,px,py,pz,chit);
         chi[indp]=chit.real();
         gx[indp]=chit.imag();
         //cout << "re: " << chi[indp] << ", im: " << gx[indp] << endl;
         indp++;
      }
   }
   indr=0;
   rhop=0.000000e0;
   complex<solreal> chiu,chiv;
   int i,j;
#pragma omp parallel for private(indr,chiu,chiv,chit) \
firstprivate(j) lastprivate(i) reduction(+: rhop)
   for (i=0; i<nPri; i++) {
      indr=i*(nPri);
      chiu=complex<solreal>(chi[i],gx[i]);
      rhop+=(cab[indr+i]*norm(chiu));
      for (j=(i+1); j<nPri; j++) {
         chiv=complex<solreal>(chi[j],gx[j]);
         chit=((chiv*conj(chiu))+(conj(chiv)*chiu));
         rhop+=(cab[indr+j]*(chit.real()));
      }
   }
   return rhop;
}
#else
/* *************************************************************************************** */
solreal GaussWaveFunction::evalFTDensity(solreal px,solreal py,solreal pz)
{
   int indr,indp,ppt;
   solreal rhop,Rx[3],alp;
   complex<solreal> chit,chiu;
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      Rx[0]=R[indr++];
      Rx[1]=R[indr++];
      Rx[2]=R[indr++];
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         evalFTChi(ppt,alp,Rx,px,py,pz,chit);
         chi[indp]=chit.real();
         gx[indp]=chit.imag();
         indp++;
      }
   }
   indr=0;
   rhop=0.000000e0;
   /*
   complex<solreal> chiv;
   for (int i=0; i<nPri; i++) {
      indr=i*(nPri+1);
      chiu=complex<solreal>(chi[i],gx[i]);
      rhop+=(cab[indr++]*norm(chiu));
      for (int j=(i+1); j<nPri; j++) {
         chiv=complex<solreal>(chi[j],gx[j]);
         chit=((chiv*conj(chiu))+(conj(chiv)*chiu));
         rhop+=(cab[indr++]*(chit.real()));
      }
   }
   / */
   solreal sumim,sumre,cc;
   indr=0;
   for ( int i=0 ; i<nPri ; ++i ) {
      sumim=sumre=0.0e0;
      for ( int j=0 ; j<nPri ; ++j ) {
         cc=cab[indr++];
         sumre+=cc*chi[j];
         sumim+=cc*gx[j];
      }
      rhop+=sumre*chi[i];
      rhop+=sumim*gx[i];
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         Rx[0]=R[indr];
         Rx[1]=R[indr+1];
         Rx[2]=R[indr+2];
         ppt=primType[i];
         alp=0.5e0*primExp[i];
         evalFTChi(ppt,alp,Rx,px,py,pz,chit);
         chi[i]=chit.real();
         gx[i]=chit.imag();
      }
      for (int i=nPri; i<totPri; ++i) {
         chiu=complex<solreal>(chi[i],gx[i]);
         rhop+=(EDFCoeff[i-nPri]*norm(chiu));
      }
   }
   return rhop;
}
#endif
/* ************************************************************************************ */
solreal GaussWaveFunction::evalFTKineticEnergy(solreal px,solreal py,solreal pz)
{
   solreal pp=evalFTDensity(px,py,pz);
   pp*=(px*px+py*py+pz*pz);
   return pp;
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalDensityMatrix1(solreal x,solreal y,solreal z,
      solreal xp,solreal yp,solreal zp)
{
   int indr,indp;
   solreal xmr,ymr,zmr,gamm,chib;
   gamm=0.000000e0;
   solreal rr;
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         chi[indp]=evalAngACases(primType[indp],xmr,ymr,zmr);
         chi[indp]*=exp(primExp[indp]*rr);
         indp++;
      }
   }
   indr=0;
   indp=0;
   for (int i=0; i<nNuc; i++) {
      xmr=xp-R[indr++];
      ymr=yp-R[indr++];
      zmr=zp-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         gx[indp]=evalAngACases(primType[indp],xmr,ymr,zmr);
         gx[indp]*=exp(primExp[indp]*rr);
         indp++;
      }
   }
   indp=0;
   gamm=0.0000000e0;
   for (int i=0; i<nPri; i++) {
      chib=0.0000000e0;
      for (int j=0; j<nPri; j++) {
         chib+=(chi[j]*cab[indp++]);
      }
      gamm+=(chib*gx[i]);
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=x-R[indr];
         ymr=y-R[indr+1];
         zmr=z-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         chi[i]=evalAngACases(primType[i],xmr,ymr,zmr);
         chi[i]*=exp(0.5e0*primExp[i]*rr);
      }
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*(primCent[i]);
         xmr=xp-R[indr];
         ymr=yp-R[indr+1];
         zmr=zp-R[indr+2];
         rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
         gx[i]=evalAngACases(primType[i],xmr,ymr,zmr);
         gx[i]*=exp(0.5e0*primExp[i]*rr);
      }
      for (int i=nPri; i<totPri; ++i) {
         gamm+=(EDFCoeff[i-nPri]*chi[i]*gx[i]);
      }
   }
return gamm;
}
/* *************************************************************************************** */
void GaussWaveFunction::evalGradDensityMatrix1(solreal x,solreal y,solreal z,\
      solreal xp,solreal yp,solreal zp,\
      solreal &gamm,solreal (&gg)[3],solreal (&gp)[3])
{
   solreal xmr,xpmr,ymr,ypmr,zmr,zpmr,rr,rrp;
   solreal cc,ccp,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr];
      xpmr=xp-R[indr++];
      ymr=y-R[indr];
      ypmr=yp-R[indr++];
      zmr=z-R[indr];
      zpmr=zp-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      rrp=-((xpmr*xpmr)+(ypmr*ypmr)+(zpmr*zpmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         ccp=exp(alp*rrp);
         hyz[indp]=evalAngACases(ppt,xpmr,ypmr,zpmr);
         hyz[indp]*=ccp;
         evalDkAngCases(ppt,alp,xpmr,ypmr,zpmr,hxx[indp],hyy[indp],hzz[indp]);
         hxx[indp]*=ccp;
         hyy[indp]*=ccp;
         hzz[indp]*=ccp;
         indp++;
      }
   }
   solreal nabx,naby,nabz,nabxp,nabyp,nabzp;
   nabx=naby=nabz=0.000000000000000e0;
   nabxp=nabyp=nabzp=0.000000000000000e0;
   indp=0;
   solreal trho=0.0000000e0,chib,chibp;
   for (int i=0; i<nPri; i++) {
      chib=chibp=0.0000000e0;
      for (int j=0; j<nPri; j++) {
         //cc=cab[indp++];
         chibp+=(hyz[j]*cab[indp]);
         chib+=(chi[j]*cab[indp++]);
      }
      trho+=(chibp*chi[i]);
      nabx+=(chibp*gx[i]);
      naby+=(chibp*gy[i]);
      nabz+=(chibp*gz[i]);
      nabxp+=(chib*hxx[i]);
      nabyp+=(chib*hyy[i]);
      nabzp+=(chib*hzz[i]);
   }
   gg[0]=nabx;
   gg[1]=naby;
   gg[2]=nabz;
   gp[0]=nabxp;
   gp[1]=nabyp;
   gp[2]=nabzp;
   gamm=trho;
   if ( ihaveEDF ) {
      cout << "Warning: evalGradDensityMatrix1 does not include EDF " <<
        "contributions!" << endl;
   }
   return;

}
/* ************************************************************************************ */
void GaussWaveFunction::evalHessDensityMatrix1(solreal (&xx)[3],solreal (&xxp)[3],\
      solreal &gamm,solreal (&gg)[3],solreal (&gp)[3],\
      solreal (&hh)[3][3],solreal (&hph)[3][3],solreal (&hp)[3][3])
{
   solreal x=xx[0],y=xx[1],z=xx[2];
   solreal xp=xxp[0],yp=xxp[1],zp=xxp[2];
   solreal xmr,xpmr,ymr,ypmr,zmr,zpmr,rr,rrp;
   solreal cc,ccp,alp;
   int indp,indr,ppt;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr];
      xpmr=xp-R[indr++];
      ymr=y-R[indr];
      ypmr=yp-R[indr++];
      zmr=z-R[indr];
      zpmr=zp-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      rrp=-((xpmr*xpmr)+(ypmr*ypmr)+(zpmr*zpmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         chi[indp]=evalAngACases(ppt,xmr,ymr,zmr);
         chi[indp]*=cc;
         evalDkAngCases(ppt,alp,xmr,ymr,zmr,gx[indp],gy[indp],gz[indp]);
         gx[indp]*=cc;
         gy[indp]*=cc;
         gz[indp]*=cc;
         ccp=exp(alp*rrp);
         hyz[indp]=evalAngACases(ppt,xpmr,ypmr,zpmr);
         hyz[indp]*=ccp;
         evalDkAngCases(ppt,alp,xpmr,ypmr,zpmr,hxx[indp],hyy[indp],hzz[indp]);
         hxx[indp]*=ccp;
         hyy[indp]*=ccp;
         hzz[indp]*=ccp;
         indp++;
      }
   }
   solreal nabx,naby,nabz,nabxp,nabyp,nabzp;
   nabx=naby=nabz=0.000000000000000e0;
   nabxp=nabyp=nabzp=0.000000000000000e0;
   indp=0;
   solreal trho=0.0000000e0,chib,chibp;
   solreal sumdiphiax,sumdiphiay,sumdiphiaz;
   solreal sumhxx,sumhyy,sumhzz,sumhxy,sumhxz,sumhyz;
   solreal sumxpy,sumxpz,sumypz;
   sumhxx=sumhyy=sumhzz=sumhxy=sumhxz=sumhyz=0.00000e0;
   sumxpy=sumxpz=sumypz=0.0e0;
   for (int i=0; i<nPri; i++) {
      chib=chibp=0.0000000e0;
      sumdiphiax=sumdiphiay=sumdiphiaz=0.0000000e0;
      for (int j=0; j<nPri; j++) {
         cc=cab[indp++];
         chibp+=(hyz[j]*cc); //double-checked
         chib+=(chi[j]*cc); //double-checked
         sumdiphiax+=(gx[j]*cc); //double-checked
         sumdiphiay+=(gy[j]*cc); //double-checked
         sumdiphiaz+=(gz[j]*cc); //double-checked
      }
      trho+=(chibp*chi[i]); //double-checked
      nabx+=(chibp*gx[i]); //double-checked
      naby+=(chibp*gy[i]); //double-checked
      nabz+=(chibp*gz[i]); //double-checked
      nabxp+=(chib*hxx[i]); //double-checked
      nabyp+=(chib*hyy[i]); //double-checked
      nabzp+=(chib*hzz[i]); //double-checked
      sumhxx+=(sumdiphiax*hxx[i]); //double-checked
      sumhyy+=(sumdiphiay*hyy[i]); //double-checked
      sumhzz+=(sumdiphiaz*hzz[i]); //double-checked
      sumhxy+=(sumdiphiax*hyy[i]); //double-checked
      sumhxz+=(sumdiphiax*hzz[i]); //double-checked
      sumhyz+=(sumdiphiay*hzz[i]); //double-checked
      sumxpy+=(sumdiphiay*hxx[i]); //double-checked
      sumxpz+=(sumdiphiaz*hxx[i]); //double-checked
      sumypz+=(sumdiphiaz*hyy[i]); //double-checked
   }
   gg[0]=nabx;
   gg[1]=naby;
   gg[2]=nabz;
   gp[0]=nabxp;
   gp[1]=nabyp;
   gp[2]=nabzp;
   gamm=trho;
   hph[0][0]=sumhxx; hph[0][1]=sumhxy; hph[0][2]=sumhxz;
   hph[1][0]=sumxpy; hph[1][1]=sumhyy; hph[1][2]=sumhyz;
   hph[2][0]=sumxpz; hph[2][1]=sumypz; hph[2][2]=sumhzz;
   for ( int i=0 ; i<nPri ; i++ ) {gx[i]=hyz[i];} //gx is now chi(xp)
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=x-R[indr++];
      ymr=y-R[indr++];
      zmr=z-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
               hxx[indp],hyy[indp],hzz[indp],
               hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   indp=0;
   sumhxx=sumhyy=sumhzz=sumhxy=sumhxz=sumhyz=0.00000e0;
   for (int i=0; i<nPri; i++) {
      //indr=i*(nPri);
      chibp=0.0000000e0;
      for (int j=0; j<nPri; j++) {
         chibp+=(gx[j]*cab[indp++]);
      }
      sumhxx+=(chibp*hxx[i]);
      sumhyy+=(chibp*hyy[i]);
      sumhzz+=(chibp*hzz[i]);
      sumhxy+=(chibp*hxy[i]);
      sumhxz+=(chibp*hxz[i]);
      sumhyz+=(chibp*hyz[i]);
   }
   hh[0][0]=sumhxx; hh[0][1]=sumhxy; hh[0][2]=sumhxz;
   hh[1][0]=sumhxy; hh[1][1]=sumhyy; hh[1][2]=sumhyz;
   hh[2][0]=sumhxz; hh[2][1]=sumhyz; hh[2][2]=sumhzz;
   indp=0;
   indr=0;
   for (int i=0; i<nNuc; i++) {
      xmr=xp-R[indr++];
      ymr=yp-R[indr++];
      zmr=zp-R[indr++];
      rr=-((xmr*xmr)+(ymr*ymr)+(zmr*zmr));
      for (int j=0; j<myPN[i]; j++) {
         ppt=primType[indp];
         alp=primExp[indp];
         cc=exp(alp*rr);
         evalDkDlAngCases(ppt,alp,xmr,ymr,zmr,
               hxx[indp],hyy[indp],hzz[indp],
               hxy[indp],hxz[indp],hyz[indp]);
         hxx[indp]*=cc;
         hyy[indp]*=cc;
         hzz[indp]*=cc;
         hxy[indp]*=cc;
         hxz[indp]*=cc;
         hyz[indp]*=cc;
         indp++;
      }
   }
   indp=0;
   sumhxx=sumhyy=sumhzz=sumhxy=sumhxz=sumhyz=0.00000e0;
   for (int i=0; i<nPri; i++) {
      chib=0.0000000e0;
      for (int j=0; j<nPri; j++) {
         chib+=(chi[j]*cab[indp++]);
      }
      sumhxx+=(chib*hxx[i]);
      sumhyy+=(chib*hyy[i]);
      sumhzz+=(chib*hzz[i]);
      sumhxy+=(chib*hxy[i]);
      sumhxz+=(chib*hxz[i]);
      sumhyz+=(chib*hyz[i]);
   }
   hp[0][0]=sumhxx; hp[0][1]=sumhxy; hp[0][2]=sumhxz;
   hp[1][0]=sumhxy; hp[1][1]=sumhyy; hp[1][2]=sumhyz;
   hp[2][0]=sumhxz; hp[2][1]=sumhyz; hp[2][2]=sumhzz;
   if ( ihaveEDF ) {
      cout << "Warning: evalHessDensityMatrix1 does not include EDF " <<
        "contributions!" << endl;
   }
   return;
}
/* ************************************************************************************ */
void GaussWaveFunction::evalHermiteCoefs(int (&aia)[3],int (&aib)[3],solreal &alpab,
      solreal (&ra)[3],solreal (&rb)[3],
      solreal (&rp)[3],
      int (&maxl)[3],solreal (&Eijl)[3][7])
{
   //for (int i=0; i<3; i++) {for (int j=0; j<=6; j++) {Eijl[i][j]=0.0e0;}}
   int ij[3];
   for (int i=0; i<3; i++) {ij[i]=10*aia[i]+aib[i];}
   if (ij[0]==0&&ij[1]==0&&ij[2]==0) {
      for (int i=0; i<3; i++) {
         maxl[i]=0;
         Eijl[i][0]=1.0e0;
      }
      return;
   }
   solreal PA[3],PB[3];
   solreal eta=0.5e0/alpab;
   for (int i=0; i<3; i++) {
      maxl[i]=aia[i]+aib[i];
      PA[i]=rp[i]-ra[i];
      PB[i]=rp[i]-rb[i];
   }
   solreal eta2,eta3,pa2,pb2,pa3,pb3,pa2p3papb,papb;
   for (int i=0; i<3; i++) {
      switch (ij[i]) {
         case 0:
            Eijl[i][0]=1.0e0;
            break;
         case 1:
            Eijl[i][0]=PB[i];
            Eijl[i][1]=eta;
            break;
         case 2:
            Eijl[i][0]=PB[i]*PB[i]+eta;
            Eijl[i][1]=2.0e0*eta*PB[i];
            Eijl[i][2]=eta*eta;
            break;
         case 3:
            pb2=PB[i]*PB[i];
            Eijl[i][0]=PB[i]*(pb2+3.0e0*eta);
            Eijl[i][1]=3.0e0*eta*(pb2+eta);
            Eijl[i][2]=3.0e0*PB[i]*eta*eta;
            Eijl[i][3]=eta*eta*eta;
            break;
         case 10:
            Eijl[i][0]=PA[i];
            Eijl[i][1]=eta;
            break;
         case 11:
            Eijl[i][0]=eta+PA[i]*PB[i];
            Eijl[i][1]=eta*(PA[i]+PB[i]);
            Eijl[i][2]=eta*eta;
            break;
         case 12:
            pb2=PB[i]*PB[i]; eta2=eta*eta;
            Eijl[i][0]=2.0e0*PB[i]*eta+PA[i]*(pb2+eta);
            Eijl[i][1]=eta*(2.0e0*PA[i]*PB[i]+pb2+3.0e0*eta);
            Eijl[i][2]=eta2*(PA[i]+2.0e0*PB[i]);
            Eijl[i][3]=eta2*eta;
            break;
         case 13:
            pb2=PB[i]*PB[i]; eta2=eta*eta;
            Eijl[i][0]=PA[i]*PB[i]*pb2+3.0e0*PB[i]*(PA[i]+PB[i])*eta+3.0e0*eta2;
            Eijl[i][1]=eta*(PB[i]*pb2+9.0e0*PB[i]*eta+3.0e0*PA[i]*(pb2+eta));
            Eijl[i][2]=3.0e0*eta2*(PB[i]*(PA[i]+PB[i])+2.0e0*eta);
            Eijl[i][3]=eta*eta2*(PA[i]+3.0e0*PB[i]);
            Eijl[i][4]=eta2*eta2;
            break;
         case 20:
            Eijl[i][0]=PA[i]*PA[i]+eta;
            Eijl[i][1]=2.0e0*PA[i]*eta;
            Eijl[i][2]=eta*eta;
            break;
         case 21:
            pa2=PA[i]*PA[i];
            Eijl[i][0]=pa2*PB[i]+(2.0e0*PA[i]+PB[i])*eta;
            Eijl[i][1]=eta*(pa2+2.0e0*PA[i]*PB[i]+3.0e0*eta);
            Eijl[i][2]=eta*eta*(2.0e0*PA[i]+PB[i]);
            Eijl[i][3]=eta*eta*eta;
            break;
         case 22:
            eta2=eta*eta; pa2=PA[i]*PA[i]; pb2=PB[i]*PB[i]; papb=PA[i]*PB[i];
            Eijl[i][0]=4.0e0*papb*eta+pa2*(pb2+eta)+eta*(pb2+3.0e0*eta);
            Eijl[i][1]=2.0e0*(PA[i]+PB[i])*eta*(papb+3.0e0*eta);
            Eijl[i][2]=eta2*(pa2+4.0e0*papb+pb2+6.0e0*eta);
            Eijl[i][3]=2.0e0*(PA[i]+PB[i])*eta*eta2;
            Eijl[i][4]=eta2*eta2;
            break;
         case 23:
            eta2=eta*eta; pa2=PA[i]*PA[i]; pb2=PB[i]*PB[i]; pb3=PB[i]*pb2;
            Eijl[i][0]=6.0e0*PA[i]*eta*(pb2+eta)+pb3*eta+9.0e0*PB[i]*eta2+pa2*(pb3+3.0e0*PB[i]*eta);
            Eijl[i][1]=eta*(3.0e0*pa2*(pb2+eta)\
                            +eta*(9.0e0*pb2+15.0e0*eta)+2.0e0*PA[i]*(pb3+9.0e0*PB[i]*eta));
            Eijl[i][2]=eta2*(3.0e0*pa2*PB[i]+pb3+18.0e0*PB[i]*eta+6.0e0*PA[i]*(pb2+2.0e0*eta));
            Eijl[i][3]=eta*eta2*(pa2+6.0e0*PA[i]*PB[i]+3.0e0*pb2+10.0e0*eta);
            Eijl[i][4]=eta2*eta2*(2.0e0*PA[i]+3.0e0*PB[i]);
            Eijl[i][5]=eta2*eta2*eta;
            break;
         case 30:
            pa2=PA[i]*PA[i];
            Eijl[i][0]=PA[i]*(pa2+3.0e0*eta);
            Eijl[i][1]=3.0e0*eta*(pa2+eta);
            Eijl[i][2]=3.0e0*PA[i]*eta*eta;
            Eijl[i][3]=eta*eta*eta;
            break;
         case 31:
            eta2=eta*eta; pa2=PA[i]*PA[i]; pa3=pa2*PA[i];
            Eijl[i][0]=pa3*PB[i]+3.0e0*PA[i]*(PA[i]+PB[i])*eta+3.0e0*eta2;
            Eijl[i][1]=eta*(pa2*(PA[i]+3.0e0*PB[i]))+eta2*(9.0e0*PA[i]+3.0e0*PB[i]);
            Eijl[i][2]=3.0e0*eta2*(PA[i]*(PA[i]+PB[i])+2.0e0*eta);
            Eijl[i][3]=eta2*eta*(3.0e0*PA[i]+PB[i]);
            Eijl[i][4]=eta2*eta2;
            break;
         case 32:
            eta2=eta*eta; eta3=eta2*eta; pa2=PA[i]*PA[i]; pa3=pa2*PA[i]; pb2=PB[i]*PB[i];
            Eijl[i][0]=6.0e0*PB[i]*(pa2*eta+eta2)+pa3*(pb2+eta)+3.0e0*PA[i]*eta*(pb2+3.0e0*eta);
            Eijl[i][1]=eta*(2.0e0*pa3*PB[i]+18.0e0*PA[i]*PB[i]*eta\
                            +3.0e0*pa2*(pb2+3.0e0*eta)+3.0e0*eta*(pb2+5.0e0*eta));
            Eijl[i][2]=eta2*(pa3+PB[i]*(6.0e0*pa2+12.0e0*eta)+3.0e0*PA[i]*(pb2+6.0e0*eta));
            Eijl[i][3]=eta3*(3.0e0*pa2+6.0e0*PA[i]*PB[i]+pb2+10.0e0*eta);
            Eijl[i][4]=eta2*eta2*(3.0e0*PA[i]+2.0e0*PB[i]);
            Eijl[i][5]=eta3*eta2;
            break;
         case 33:
            eta2=eta*eta; eta3=eta2*eta;
            pa2=PA[i]*PA[i]; pa3=pa2*PA[i]; pb2=PB[i]*PB[i]; pb3=pb2*PB[i]; papb=PA[i]*PB[i];
            pa2p3papb=pa2+3.0e0*papb+pb2;
            Eijl[i][0]=pa3*pb3+pa2p3papb*(3.0e0*papb*eta+9.0e0*eta2)+15.0e0*eta3;
            Eijl[i][1]=3.0e0*(PA[i]+PB[i])*eta*(8.0e0*papb*eta+pa2*(pb2+eta)+eta*(pb2+15.0e0*eta));
            Eijl[i][2]=3.0e0*eta2*(pa2p3papb*(papb+6.0e0*eta)+15.0e0*eta2);
            Eijl[i][3]=eta3*(PA[i]+PB[i])*(pa2+8.0e0*papb+pb2+30.0e0*eta);
            Eijl[i][4]=3.0e0*eta2*eta2*(pa2p3papb+5.0e0*eta);
            Eijl[i][5]=3.0e0*eta2*eta3*(PA[i]+PB[i]);
            Eijl[i][6]=eta3*eta3;
            break;
         default:
            cout << "Underconstruction..." << endl;
            break;
      }
   }
}
/* *************************************************************************************** */
void GaussWaveFunction::evalRlmnIntegs(const int (&lmn)[3],const solreal &alpp,const solreal (&cp)[3],
                                   solreal (&Rlmn)[7][7][7])
{
   int maxj=lmn[0]+lmn[1]+lmn[2];
   solreal a=cp[0],b=cp[1],c=cp[2];
   solreal T=alpp*(a*a+b*b+c*c);
   solreal Fj[7];
   BoysFunction(T,maxj,Fj);
   solreal Rjlmn[7][7][7][7];
   solreal m2a=-2.0e0*alpp,m2aj=1.0e0;
   for (int j=0; j<=maxj; j++) {
      Rjlmn[j][0][0][0]=Fj[j]*m2aj;
      m2aj*=m2a;
   }
   solreal rel,rem,ren;
   for (int l=1; l<=lmn[0]; l++) {
      for (int j=maxj-l; j>=0; j--) {
         Rjlmn[j][l][0][0]=a*Rjlmn[j+1][l-1][0][0];
         if (l>1) {
            rel=solreal(l-1);
            Rjlmn[j][l][0][0]+=(rel*Rjlmn[j+1][l-2][0][0]);
         }
      }
   }
   for (int l=0; l<=lmn[0]; l++) {
      for (int m=1; m<=lmn[1]; m++) {
         for (int j=maxj-l-m; j>=0; j--) {
            Rjlmn[j][l][m][0]=b*Rjlmn[j+1][l][m-1][0];
            if (m>1) {
               rem=solreal(m-1);
               Rjlmn[j][l][m][0]+=(rem*Rjlmn[j+1][l][m-2][0]);
            }
         }
      }
   }
   for (int l=0; l<=lmn[0]; l++) {
      for (int m=0; m<=lmn[1]; m++) {
         for (int n=1; n<=lmn[2]; n++) {
            for (int j=maxj-l-m-n; j>=0; j--) {
               Rjlmn[j][l][m][n]=c*Rjlmn[j+1][l][m][n-1];
               if (n>1) {
                  ren=solreal(n-1);
                  Rjlmn[j][l][m][n]+=(ren*Rjlmn[j+1][l][m][n-2]);
               }
            }
         }
      }
   }
   for (int l=0; l<=lmn[0]; l++) {
      for (int m=0; m<=lmn[1]; m++) {
         for (int n=0; n<=lmn[2]; n++) {Rlmn[l][m][n]=Rjlmn[0][l][m][n];}
      }
   }
   return;
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalVAB(solreal (&xx)[3],int (&aa)[3],int (&ab)[3],solreal &alpa,solreal &alpb,
                               solreal (&xa)[3],solreal (&xb)[3])
{
   solreal alpp=alpa+alpb;
   solreal ooalpp=1.0e0/alpp,xp[3],cp[3],ctmp=0.0e0,S00;
   int maxl[3];
   for (int i=0; i<3; i++) {
      xp[i]=ooalpp*(alpa*xa[i]+alpb*xb[i]);
      cp[i]=xp[i]-xx[i];
      S00=xa[i]-xb[i];
      ctmp+=S00*S00;
   }
   ctmp*=(alpa*alpb*ooalpp);
   S00=exp(-ctmp);
   if (S00<EPSFORMEPVALUE) {return 0.0e0;}
   solreal Eabk[3][7],Rijk[7][7][7];
   evalHermiteCoefs(aa,ab,alpp,xa,xb,xp,maxl,Eabk);
   evalRlmnIntegs(maxl,alpp,cp,Rijk);
   ctmp=0.0e0;
   for (int i=0; i<=maxl[0]; i++) {
      for (int j=0; j<=maxl[1]; j++) {
         for (int k=0; k<=maxl[2]; k++) {
            ctmp+=(Eabk[0][i]*Eabk[1][j]*Eabk[2][k]*Rijk[i][j][k]);
         }
      }
   }
   static const solreal twopi=6.2831853071795864769;
   return (twopi*ooalpp*S00*ctmp);
}
/* *************************************************************************************** */
#if PARALLELISEDTK
solreal GaussWaveFunction::evalMolElecPot(solreal x,solreal y,solreal z)
{
#if DEBUG
   static bool showmsg=true;
#endif
   if ( maxPrimType>19 ) {
#if DEBUG
      if (showmsg) {
         displayErrorMessage(string("Non supported angular momenta of"
            "primitives, requesting type: ")\
            +getStringFromInt(maxPrimType));
         showmsg=false;

      }
#endif
      return 0.0e0;
   }
   int indr,inda,indp;
   solreal xx[3],ra[3],rb[3],alpa,alpb,mepaa,mepab,mepelec,cc;
   int aa[3],ab[3],i,j;
   xx[0]=x; xx[1]=y; xx[2]=z;
   mepaa=mepab=mepelec=0.000000e0;
   inda=indp=0;
#pragma omp parallel for private(indr,inda,ra,rb,aa,ab,alpa,alpb,indp,cc) \
firstprivate(j) lastprivate(i) shared(xx) reduction(+: mepaa,mepab)
   for (i=0; i<nPri; i++) {
      indr=3*primCent[i];
      ra[0]=R[indr];
      ra[1]=R[indr+1];
      ra[2]=R[indr+2];
      inda=3*primType[i];
      aa[0]=prTy[inda];
      aa[1]=prTy[inda+1];
      aa[2]=prTy[inda+2];
      alpa=primExp[i];
      indp=i*(nPri);
      for (j=(i+1); j<nPri; j++) {
         cc=cab[indp+j];
         indr=3*primCent[j];
         rb[0]=R[indr];
         rb[1]=R[indr+1];
         rb[2]=R[indr+2];
         inda=3*primType[j];
         ab[0]=prTy[inda];
         ab[1]=prTy[inda+1];
         ab[2]=prTy[inda+2];
         alpb=primExp[j];
         mepab+=(cc*evalVAB(xx,aa,ab,alpa,alpb,ra,rb));
      }
      mepaa+=(cab[indp+i]*evalVAB(xx,aa,aa,alpa,alpa,ra,ra));
   }
   mepelec=mepaa+2.0e0*mepab;
   mepab=0.0e0;
   for (i=0; i<nNuc; i++) {
      indr=3*i;
      alpb=0.0e0;
      for (int k=0; k<3; k++) {
         alpa=(xx[k]-R[indr+k]);
         alpb+=(alpa*alpa);
      }
      mepab+=(atCharge[i]/sqrt(alpb));
   }
   return (mepab-mepelec);
}
#else
/* *************************************************************************************** */
solreal GaussWaveFunction::evalMolElecPot(solreal x,solreal y,solreal z)
{
#if DEBUG
   static bool showmsg=true;
#endif
   if ( maxPrimType>19 ) {
#if DEBUG
      if (showmsg) {
         displayErrorMessage(string("Non supported angular momenta of"
            "primitives, requesting type: ")\
            +getStringFromInt(maxPrimType));
         showmsg=false;

      }
#endif
      return 0.0e0;
   }
   int indr,inda,indp;
   solreal xx[3],ra[3],rb[3],alpa,alpb,mepaa,mepab,mepelec,cc;
   int aa[3],ab[3];
   xx[0]=x; xx[1]=y; xx[2]=z;
   mepaa=mepab=mepelec=0.000000e0;
   inda=indp=0;
   for (int i=0; i<nPri; i++) {
      indr=3*primCent[i];
      ra[0]=R[indr];
      ra[1]=R[indr+1];
      ra[2]=R[indr+2];
      inda=3*primType[i];
      aa[0]=prTy[inda];
      aa[1]=prTy[inda+1];
      aa[2]=prTy[inda+2];
      alpa=primExp[i];
      indp=i*(nPri+1);
      cc=cab[indp++];
      mepaa+=(cc*evalVAB(xx,aa,aa,alpa,alpa,ra,ra));
      for (int j=(i+1); j<nPri; j++) {
         cc=cab[indp++];
         indr=3*primCent[j];
         rb[0]=R[indr];
         rb[1]=R[indr+1];
         rb[2]=R[indr+2];
         inda=3*primType[j];
         ab[0]=prTy[inda];
         ab[1]=prTy[inda+1];
         ab[2]=prTy[inda+2];
         alpb=primExp[j];
         mepab+=(cc*evalVAB(xx,aa,ab,alpa,alpb,ra,rb));
      }
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*primCent[i];
         ra[0]=R[indr];
         ra[1]=R[indr+1];
         ra[2]=R[indr+2];
         inda=3*primType[i];
         aa[0]=prTy[inda];
         aa[1]=prTy[inda+1];
         aa[2]=prTy[inda+2];
         alpa=0.5e0*primExp[i];
         cc=EDFCoeff[i-nPri];
         mepaa+=(cc*evalVAB(xx,aa,aa,alpa,alpa,ra,ra));
      }
   }
   mepelec=mepaa+2.0e0*mepab;
   mepab=0.0e0;
   for (int i=0; i<nNuc; i++) {
      indr=3*i;
      alpb=0.0e0;
      for (int k=0; k<3; k++) {
         alpa=(xx[k]-R[indr+k]);
         alpb+=(alpa*alpa);
      }
      mepab+=(atCharge[i]/sqrt(alpb));
   }
   return (mepab-mepelec);
}
#endif
/* *************************************************************************************** */
solreal GaussWaveFunction::integralRho(void)
{
   int indr,inda,indp;
   solreal ra[3],rb[3],alpa,alpb,rhoaa,rhoab,cc;
   int aa[3],ab[3];
   rhoaa=rhoab=0.000000e0;
   inda=indp=0;
   for (int i=0; i<nPri; i++) {
      indr=3*primCent[i];
      ra[0]=R[indr];
      ra[1]=R[indr+1];
      ra[2]=R[indr+2];
      inda=3*primType[i];
      aa[0]=prTy[inda];
      aa[1]=prTy[inda+1];
      aa[2]=prTy[inda+2];
      alpa=primExp[i];
      indp=i*(nPri+1);
      cc=cab[indp++];
      rhoaa+=(cc*evalOverlapIntegralAB(aa,aa,alpa,alpa,ra,ra));
      //mepab=0.0e0;
      for (int j=(i+1); j<nPri; j++) {
         cc=cab[indp++];
         indr=3*primCent[j];
         rb[0]=R[indr];
         rb[1]=R[indr+1];
         rb[2]=R[indr+2];
         inda=3*primType[j];
         ab[0]=prTy[inda];
         ab[1]=prTy[inda+1];
         ab[2]=prTy[inda+2];
         alpb=primExp[j];
         rhoab+=(cc*evalOverlapIntegralAB(aa,ab,alpa,alpb,ra,rb));
      }
   }
   if ( ihaveEDF ) {
      for ( int i=nPri ; i<totPri ; ++i ) {
         indr=3*primCent[i];
         ra[0]=R[indr];
         ra[1]=R[indr+1];
         ra[2]=R[indr+2];
         inda=3*primType[i];
         aa[0]=prTy[inda];
         aa[1]=prTy[inda+1];
         aa[2]=prTy[inda+2];
         alpa=0.5e0*primExp[i];
         cc=EDFCoeff[i-nPri];
         rhoaa+=(cc*evalOverlapIntegralAB(aa,aa,alpa,alpa,ra,ra));
      }
   }
   return (rhoaa+2.0e0*rhoab);
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalOverlapIntegralAB(int (&aa)[3],int (&ab)[3],solreal &alpa,solreal &alpb,
      solreal (&ra)[3],solreal (&rb)[3])
{
   solreal alpp=alpa+alpb;
   solreal ooalpp=1.0e0/alpp,xp[3],ctmp=0.0e0,S00;
   int maxl[3];
   for (int i=0; i<3; i++) {
      xp[i]=ooalpp*(alpa*ra[i]+alpb*rb[i]);
      S00=ra[i]-rb[i];
      ctmp+=S00*S00;
   }
   ctmp*=(alpa*alpb*ooalpp);
   S00=exp(-ctmp);
   solreal Eabk[3][7];
   evalHermiteCoefs(aa,ab,alpp,ra,rb,xp,maxl,Eabk);
   ctmp=Eabk[0][0]*Eabk[1][0]*Eabk[2][0];
   ooalpp*=3.1415926535897932385;
   return (ctmp*sqrt(ooalpp*ooalpp*ooalpp)*S00);
}
/* *************************************************************************************** */
solreal GaussWaveFunction::totalNuclearCharge(void)
{
   solreal nc=0.0e0;
   for (int i=0; i<nNuc; i++) {nc+=atCharge[i];}
   if ( ihaveEDF ) { nc+=coreElec; }
   return nc;
}
/* *************************************************************************************** */
void GaussWaveFunction::evalLED(solreal const (&x)[3],solreal (&led)[3])
{
   solreal rho,g[3];
   evalRhoGradRho(x[0],x[1],x[2],rho,g);
   if ( rho<1.0e-12 ) {rho=1.0e-12;}
   for ( int i=0 ; i<3 ; i++ ) {led[i]=(-0.5e0*g[i]/rho);}
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalMagLED(solreal x,solreal y,solreal z)
{
   solreal rho,g[3];
   evalRhoGradRho(x,y,z,rho,g);
   if ( rho<1.0e-12 ) {rho=1.0e-12;}
   solreal led=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {led+=(g[i]*g[i]);}
   led=sqrt(led);
   led*=0.5e0;
   led/=rho;
   return led;
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalReducedDensityGradient(solreal x,solreal y,solreal z)
{
   static const solreal cc=0.161620459673995481331661e0; /* $(2(3\pi^2)^{1/3})^{-1}$  */
   static const solreal fouo3=4.0e0/3.0e0;
   solreal rho,g[3];
   evalRhoGradRho(x,y,z,rho,g);
   if ( rho<1.0e-10 ) {rho=1.0e-10;}
   return (cc*sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2])/pow(rho,fouo3));
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalRoSE(solreal x,solreal y,solreal z)
{
   static const solreal cc=2.87123400018819181594250e0; /* (3/10)(3\pi^2)^{2/3}  */
   static const solreal fivo3=5.0e0/3.0e0;
   solreal rho,tau;
   evalNabPhi2(x,y,z,rho,tau); //Here tau=2G
   tau*=0.5e0; //tau=G
   solreal tau0=cc*pow(rho,fivo3);
   return ((tau0-tau)/(tau0+tau));
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalVirialPotentialEnergyDensity(solreal x, solreal y, solreal z)
{
   solreal lapRho,G;
   lapRho=evalLapRho(x,y,z);
   G=evalKineticEnergyG(x,y,z);
   return (0.25e0*lapRho-2.0e0*G);
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalNCIs(solreal x,solreal y,solreal z)
{
   solreal s,rho=evalDensity(x,y,z);
   if (nciRhoMin<rho && rho <= nciRhoMax) {
      s=evalReducedDensityGradient(x,y,z);
      s=(s<=nciSMax ?  s : 100.0e0);
   } else {
      s=100.0e0;
   }
   return s;
}
/* *************************************************************************************** */
solreal GaussWaveFunction::evalNCILambda(solreal x,solreal y,solreal z)
{
   solreal rho=evalDensity(x,y,z);
   if (nciRhoMin<rho && rho <= nciRhoMax) {
      solreal hess[3][3];
      solreal eingvectors[3][3],eingvalues[3];
      evalHessian(x,y,z,hess);
      eigen_decomposition3(hess,eingvectors,eingvalues);
      rho=(eingvalues[1]<0.0e0 ? (-rho) : rho);
   } else{
      rho=100.0e0;
   }
   return rho;
}
/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */
/* The following line includes the implementation of the custom scalar and vector fields.  */
#include "custfld-wfnclass.cpp"
/* *************************************************************************************** */
#endif//_GAUSSWAVEFUNCTION_CPP_

