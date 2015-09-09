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



#ifndef _BONDNETWORK_CPP_
#define _BONDNETWORK_CPP_

#ifndef CHOOSEPOVVERSION36
#define CHOOSEPOVVERSION36 0
#endif

#include "bondnetwork.h"
#include "solpovtools.h"
#include "iofuncts-wfn.h"
#include "iofuncts-wfx.h"
#include "atomradiicust.h"
// The first 94 atomic radii are given,
//  the rest are set to be 0.80e0
//
//#include "atomcolschcust.h" //Choose this for the palette defined by JMHP
#include "atomcolschjmol.h" //Choose this for the palette used in JMol

//**********************************************************************************************
bondNetWork::bondNetWork()
//**********************************************************************************************
{
   R=NULL;
   bondDist=NULL;
   nNuc=0;
   atNum=NULL;
   atLbl=NULL;
   bNet=NULL;
   nBonds=0;
   nTit=0;
   title=NULL;
   drawAtSize=AUTOMATICDRAWATOMSIZE;
   drawStickSize=AUTOMATICDRAWATOMSIZE*AUTOMATICBALLANDSTICKRATIO;
   for (int i=0; i<3; i++) {rmax[i]=-1.0e50; rmin[i]=1.0e+50;}
   for (int i=0; i<3; i++) {bbmax[i]=-1.0e50; bbmin[i]=1.0e+50;}
   maxBondDist=-1.0e50;
   ballAndStickMode=true;
   spaceFillingMode=false;
   wireMode=false;
   isSTP=false;
}
//**********************************************************************************************
bondNetWork::~bondNetWork()
//**********************************************************************************************
{
   dealloc2DRealArray(R,nNuc);
   dealloc2DIntArray(bNet,nNuc);
   dealloc1DStringArray(atLbl);
   dealloc1DStringArray(title);
   dealloc1DIntArray(atNum);
   dealloc2DRealArray(bondDist,nNuc);
   isSTP=false;
}
//**********************************************************************************************
bool bondNetWork::imstp()
{
   return isSTP;
}
//**********************************************************************************************
bool bondNetWork::readFromFileWFX(string inname)
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
   getNofNucleiFromFileWFX(tif,nNuc);
   alloc1DStringArray("atLbl",nNuc,atLbl);
   alloc2DRealArray(string("R"),nNuc,3,R);
   alloc1DIntArray(string("atNum"),nNuc,atNum);
   getAtLabelsFromFileWFX(tif,nNuc,atLbl);
   getNucCartCoordsFromFileWFX(tif,nNuc,R);
   getAtNumbersFromFileWFX(tif,nNuc,atNum);
   tif.close();
   return true;
}
//**********************************************************************************************
bool bondNetWork::readFromFileWFN(string inname)
{
   ifstream tif;
   string orDe;
   int nmo,npr;
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
   processFirstDataStringinWFNFile(tif,title,orDe,nmo,npr,nNuc);
   solreal *atch,*tmprad;
   processCentersWFN(tif,nNuc,atLbl,tmprad,atch);
   alloc2DRealArray(string("R-readwfn-"),nNuc,3,R);
   alloc1DIntArray(string("atNum"),nNuc,atNum);
   for (int i=0; i<nNuc; i++) {
      atNum[i]=floor(atch[i])-1;
      R[i][0]=tmprad[3*i];
      R[i][1]=tmprad[3*i+1];
      R[i][2]=tmprad[3*i+2];
   }
   dealloc1DRealArray(atch);
   dealloc1DRealArray(tmprad);
   tif.close();
   return true;
}
//*************************************************************************************************
bool bondNetWork::readFromFile(string inname)
{
   string extension;
   extension=inname.substr(inname.length()-3,3);
   if ((extension=="wfn")||(extension=="WFN")) {
      return readFromFileWFN(inname);
   } else if ((extension=="wfx")||(extension=="WFX")) {
      return readFromFileWFX(inname);
   } else {
      cout << "Error: unknown extension!(" << inname << ")\nNothig to do, returning false...\n";
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
}
//**********************************************************************************************
solreal bondNetWork::dist(int i, int j)
{
#if DEBUG
   if ((i>=nNuc)||(j>=nNuc)) {
      cout << "Error: Attempting to obtain information from a non-existent nucleus...\n";
   }
#endif
   solreal d;
   d=(R[i][0]-R[j][0])*(R[i][0]-R[j][0]);
   d+=(R[i][1]-R[j][1])*(R[i][1]-R[j][1]);
   d+=(R[i][2]-R[j][2])*(R[i][2]-R[j][2]);
   return sqrt(d);
}
//**********************************************************************************************
bool bondNetWork::setUpBNW(void)
{
   lookForBonds();
   seekRMaxMin();
   isSTP=true;
   setBoundingBox();
   return isSTP;
}
//**********************************************************************************************
bool bondNetWork::lookForBonds(void)
{
   if (!(alloc2DRealArray(string("bondDist"),nNuc,MAXBONDINGATOMS,bondDist))) {
      cout << "Unknonw error from " << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   if (!(alloc2DIntArray(string("bNet"),nNuc,MAXBONDINGATOMS,bNet))) {
      cout << "Unknonw error from " << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   solreal d,vdwd,clsstatd=1.0e+50;
   for (int i=0; i<nNuc; i++) {
      for (int j=i+1; j<nNuc; j++) {
         d=dist(i,j);
         //vdwd=(atomicRadius[atNum[i]]+atomicRadius[atNum[j]])/BOHRTOANGSTROM;
         vdwd=(getAtomicVDWRadius(atNum[i])+getAtomicVDWRadius(atNum[j]))/BOHRTOANGSTROM;
         if (d<=vdwd) {
            addBond(i,j,d);
            if (d>maxBondDist) {
               maxBondDist=d;
            }
            nBonds++;
         }
         if ( d<clsstatd ) {clsstatd=d;}
      }
   }
   if (nNuc==1) {maxBondDist=AUTOMATICMAXBONDDIST;}
   if ((maxBondDist<0.0e0)&&(nNuc==2)) {maxBondDist=dist(0,1);}
   if (maxBondDist<0.0e0) {maxBondDist=clsstatd;}
   return true;
}
//**********************************************************************************************
void bondNetWork::addBond(int m, int n,solreal dd)
{
   int k=0;
   while (k<MAXBONDINGATOMS) {
      if (bNet[m][k]>0) {
         k++;
      } else {
         bNet[m][k]=n;
         bondDist[m][k]=dd;
         break;
      }
   }
}
//**********************************************************************************************
void bondNetWork::seekRMaxMin(void)
{
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<3; j++) {
         if (R[i][j]<rmin[j]) {rmin[j]=R[i][j];}
         if (R[i][j]>rmax[j]) {rmax[j]=R[i][j];}
      }
   }
   return;
}
//**********************************************************************************************
bool bondNetWork::makePOVFile(string pnam,povRayConfProp &pvp)
{
   ofstream pof;
   pof.open(pnam.c_str(),ios::out);
   if (!(pof.good())) {
      cout << "Error: File " << pnam << "could not be opened...\n";
#if DEBUG
      cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
      return false;
   }
   //povRayConfProp pvp;
   if (!isSTP) {setUpBNW();}
   if (!ballAndStickMode) {drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
   centerMolecule();
   calcViewRadius();
   cout << "rView: " << rView << endl;
   for (int i=0; i<3; i++) {
      pvp.locCam[i]*=rView;
      for (int j=0; j<2; j++) {
         pvp.lightSource[j][i]*=(rView*2.0e0);
      }
   }
#if CHOOSEPOVVERSION36
   pof << "#version 3.6;" << endl;
#endif
   pvp.writeHeader(pof);
   pof << endl;
   if (spaceFillingMode) {
      pof << indTabsStr(pvp.currIndLev) << "merge{" << endl;
      pvp.currIndLev++;
   }
   putNuclei(pof);
   if (spaceFillingMode) {
      pvp.currIndLev--;
      pof << indTabsStr(pvp.currIndLev) << "}" << endl;
   }
   if (ballAndStickMode||wireMode) {
      putBonds(pof);
   }
   pof.close();
   return true;
}
//**********************************************************************************************
void bondNetWork::putNuclei(ofstream & pof)
{
   int atomn;
   solreal atrad;
   for (int i=0; i<nNuc; i++) {
      atomn=atNum[i];
      if (spaceFillingMode) {
         //atrad=atomicRadius[atomn]*AUTOMATICSPACEFILLINGRATIO;
         atrad=getAtomicVDWRadius(atomn)*AUTOMATICSPACEFILLINGRATIO;
      } else {
         atrad=drawAtSize;
      }
      //writePOVSphere(pof,0,R[i][0],R[i][1],R[i][2],atrad,
      //               atomicColor[atomn][0],atomicColor[atomn][1],atomicColor[atomn][2]);
      writePOVSphere(pof,0,R[i][0],R[i][1],R[i][2],atrad,
                     getAtomicRColorReal(atomn),getAtomicGColorReal(atomn),
                     getAtomicBColorReal(atomn));
   }
   return;
}
//**********************************************************************************************
void bondNetWork::putBonds(ofstream &pof)
{
   pof << "union{" << endl;
   int k=0,atni,atnk;
   solreal startpt[3],frak1;
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bNet[i][j];
         atni=atNum[i];
         atnk=atNum[k];
         //frak1=atomicRadius[atni]/(atomicRadius[atni]+atomicRadius[atnk]);
         frak1=getAtomicVDWRadius(atni)/(getAtomicVDWRadius(atni)+getAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=R[i][l]*(1.0e0-frak1)+R[k][l]*frak1;
         }
         if (k>0) {
            writePOVCylinder(pof,1,
                             R[i][0],R[i][1],R[i][2],
                             startpt[0],startpt[1],startpt[2],drawStickSize,
                             getAtomicRColorReal(atni),getAtomicGColorReal(atni),
                             getAtomicBColorReal(atni));
            writePOVCylinder(pof,1,
                             startpt[0],startpt[1],startpt[2],
                             R[k][0],R[k][1],R[k][2],drawStickSize,
                             getAtomicRColorReal(atnk),getAtomicGColorReal(atnk),
                             getAtomicBColorReal(atnk));
         }
      }
   }
   pof << "}" << endl;
   return;
}
//**********************************************************************************************
void bondNetWork::centerMolecule()
{
   solreal trn[3];
   for (int i=0; i<3; i++) {
      trn[i]=0.5e0*(rmax[i]+rmin[i]);
   }
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<3; j++) {
         R[i][j]-=trn[j];
      }
   }
   return;
}
//**********************************************************************************************
void bondNetWork::calcViewRadius(void)
{
   seekRMaxMin();
   solreal rmagmax=0.0e0,rmagmin=0.0e0;
   for (int i=0; i<3; i++) {
      rmagmax+=(rmax[i]*rmax[i]);
      rmagmin+=(rmin[i]*rmin[i]);
      if (rmagmax>rmagmin) {
         rView=sqrt(rmagmax);
      } else {
         rView=sqrt(rmagmin);
      }
   }
}
//**********************************************************************************************
void bondNetWork::setBoundingBox(void)
{
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<3; j++) {
         if (R[i][j]<=bbmin[j]) {bbmin[j]=R[i][j];}
         if (R[i][j]>=bbmax[j]) {bbmax[j]=R[i][j];}
      }
   }
   if (!isSTP) {
      cout << "Error: The bondNetWork object is not set up!\n";
   } else {
      for (int i=0; i<3; i++) {
         bbmax[i]+=(0.5e0*maxBondDist);
         bbmin[i]-=(0.5e0*maxBondDist);
      }
   }
   return;
}
//**********************************************************************************************
int bondNetWork::countAtomsOfAtomicNumber(int nat)
{
#if DEBUG
   if (nat<1) {
      displayErrorMessage("Invalid atomic number!\nReturning zero.");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return 0;
   }
#endif
   int n=nat-1,k=0;
   for (int i=0; i<nNuc; i++) {if (atNum[i]==n) {k++;}}
   return k;
}
//**********************************************************************************************
//**********************************************************************************************
#endif//_BONDNETWORK_CPP_

