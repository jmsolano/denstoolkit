/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
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
#ifndef _BONDNETWORK_CPP_
#define _BONDNETWORK_CPP_

#ifndef CHOOSEPOVVERSION36
#define CHOOSEPOVVERSION36 0
#endif
#include <cstdlib>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include "bondnetwork.h"
#include "povraytools.h"
#include "iofuncts-wfn.h"
#include "iofuncts-wfx.h"
#include "iofuncts-cpx.h"
#include "atomradiicust.h"
#include "inputmolecule_xyz.h"
// The first 94 atomic radii are given,
//  the rest are set to be 0.80e0
//
//#include "atomcolschcust.h" //Choose this for the palette defined by JMHP
#include "atomcolschjmol.h" //Choose this for the palette used in JMol

BondNetWork::BondNetWork() {
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
BondNetWork::~BondNetWork() {
   MyMemory::Dealloc2DRealArray(R,nNuc);
   MyMemory::Dealloc2DIntArray(bNet,nNuc);
   MyMemory::Dealloc1DStringArray(atLbl);
   MyMemory::Dealloc1DStringArray(title);
   MyMemory::Dealloc1DIntArray(atNum);
   MyMemory::Dealloc2DRealArray(bondDist,nNuc);
   isSTP=false;
}
bool BondNetWork::ImStp() {
   return isSTP;
}
bool BondNetWork::ReadFromFileWFX(string inname) {
   ifstream tif;
   tif.open(inname.c_str(),std::ios::in);
   if (!(tif.good())) {
      ScreenUtils::DisplayErrorFileNotOpen(inname);
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   GetTitleFromFileWFX(tif,nTit,title);
   GetNofNucleiFromFileWFX(tif,nNuc);
   MyMemory::Alloc1DStringArray("atLbl",nNuc,atLbl);
   MyMemory::Alloc2DRealArray(string("R"),nNuc,3,R);
   MyMemory::Alloc1DIntArray(string("atNum"),nNuc,atNum);
   GetAtLabelsFromFileWFX(tif,nNuc,atLbl);
   GetNucCartCoordsFromFileWFX(tif,nNuc,R);
   GetAtNumbersFromFileWFX(tif,nNuc,atNum);
   tif.close();
   return true;
}
bool BondNetWork::ReadFromFileWFN(string inname) {
   ifstream tif;
   string orDe;
   int nmo,npr;
   tif.open(inname.c_str(),std::ios::in);
   if (!(tif.good())) {
      ScreenUtils::DisplayErrorFileNotOpen(inname);
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   tif.seekg(tif.beg);
   nTit=1;
   ProcessFirstDataStringinWFNFile(tif,title,orDe,nmo,npr,nNuc);
   double *atch,*tmprad;
   ProcessCentersWFN(tif,nNuc,atLbl,tmprad,atch);
   MyMemory::Alloc2DRealArray(string("R-readwfn-"),nNuc,3,R);
   MyMemory::Alloc1DIntArray(string("atNum"),nNuc,atNum);
   for (int i=0; i<nNuc; i++) {
      atNum[i]=floor(atch[i])-1;
      R[i][0]=tmprad[3*i];
      R[i][1]=tmprad[3*i+1];
      R[i][2]=tmprad[3*i+2];
   }
   MyMemory::Dealloc1DRealArray(atch);
   MyMemory::Dealloc1DRealArray(tmprad);
   tif.close();
   return true;
}
bool BondNetWork::ReadFromFileXYZ(string inname) {
   InputMoleculeXYZ mol(inname);
   nTit=1;
   MyMemory::Alloc1DStringArray("title",nTit,title);
   title[0]=mol.title;
   nNuc=int(mol.Size());
   MyMemory::Alloc1DStringArray("atLbl",nNuc,atLbl);
   MyMemory::Alloc2DRealArray(string("R"),nNuc,3,R);
   MyMemory::Alloc1DIntArray(string("atNum"),nNuc,atNum);
   for ( size_t i=0 ; i<mol.Size() ; ++i ) {
      atLbl[i]=mol.atom[i].symbol+std::to_string(i+1);
   }
   for ( size_t i=0 ; i<mol.Size() ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) {
         R[i][j]=mol.atom[i].x[j];
         R[i][j]/=(BOHRTOANGSTROM);
      }
   }
   for ( size_t i=0 ; i<mol.Size() ; ++i ) { atNum[i]=mol.atom[i].num-1; }
   return true;
}
bool BondNetWork::ReadFromFileCPX(string inname) {
   ifstream ifil(inname);
   string wfname=cpxGetWFXFileName(ifil);
   ifil.close();
   cout << "wfname: '" << wfname << "'" << '\n';
   string extension=wfname.substr(wfname.length()-3,3);
   if ((extension=="wfn")||(extension=="WFN")) {
      return ReadFromFileWFN(wfname);
   } else if ((extension=="wfx")||(extension=="WFX")) {
      return ReadFromFileWFX(wfname);
   }
   return false;
}
bool BondNetWork::ReadFromFile(string inname) {
   string extension;
   extension=inname.substr(inname.length()-3,3);
   if ((extension=="wfn")||(extension=="WFN")) {
      return ReadFromFileWFN(inname);
   } else if ((extension=="wfx")||(extension=="WFX")) {
      return ReadFromFileWFX(inname);
   } else if ((extension=="xyz")||(extension=="XYZ")) {
     return ReadFromFileXYZ(inname);
   } else if ((extension=="cpx")||(extension=="CPX")) {
     return ReadFromFileCPX(inname);
   } else {
      cout << "Error: unknown extension!(" << inname << ")\nNothig to do, returning false...\n";
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
}
double BondNetWork::Dist(int i, int j) {
#if DEBUG
   if ((i>=nNuc)||(j>=nNuc)) {
      cout << "Error: Attempting to obtain information from a non-existent nucleus...\n";
   }
#endif
   double d;
   d=(R[i][0]-R[j][0])*(R[i][0]-R[j][0]);
   d+=(R[i][1]-R[j][1])*(R[i][1]-R[j][1]);
   d+=(R[i][2]-R[j][2])*(R[i][2]-R[j][2]);
   return sqrt(d);
}
bool BondNetWork::SetUpBNW(void) {
   LookForBonds();
   SeekRMaxMin();
   isSTP=true;
   SetBoundingBox();
   return isSTP;
}
bool BondNetWork::LookForBonds(void) {
   if (!(MyMemory::Alloc2DRealArray(string("bondDist"),nNuc,MAXBONDINGATOMS,bondDist))) {
      cout << "Unknonw error from " << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   if (!(MyMemory::Alloc2DIntArray(string("bNet"),nNuc,MAXBONDINGATOMS,bNet))) {
      cout << "Unknonw error from " << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   double d,vdwd,clsstatd=1.0e+50;
   for (int i=0; i<nNuc; i++) {
      for (int j=i+1; j<nNuc; j++) {
         d=Dist(i,j);
         //vdwd=(atomicRadius[atNum[i]]+atomicRadius[atNum[j]])/BOHRTOANGSTROM;
         vdwd=(GetAtomicVDWRadius(atNum[i])+GetAtomicVDWRadius(atNum[j]))/BOHRTOANGSTROM;
         if (d<=vdwd) {
            AddBond(i,j,d);
            if (d>maxBondDist) {
               maxBondDist=d;
            }
            nBonds++;
         }
         if ( d<clsstatd ) {clsstatd=d;}
      }
   }
   if (nNuc==1) {maxBondDist=AUTOMATICMAXBONDDIST;}
   if ((maxBondDist<0.0e0)&&(nNuc==2)) {maxBondDist=Dist(0,1);}
   if (maxBondDist<0.0e0) {maxBondDist=clsstatd;}
   return true;
}
void BondNetWork::AddBond(int m, int n,double dd) {
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
void BondNetWork::SeekRMaxMin(void) {
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<3; j++) {
         if (R[i][j]<rmin[j]) {rmin[j]=R[i][j];}
         if (R[i][j]>rmax[j]) {rmax[j]=R[i][j];}
      }
   }
   return;
}
bool BondNetWork::MakePOVFile(string pnam,POVRayConfiguration &pvp) {
   ofstream pof;
   pof.open(pnam.c_str(),std::ios::out);
   if (!(pof.good())) {
      ScreenUtils::DisplayErrorFileNotOpen(pnam);
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return false;
   }
   //POVRayConfiguration pvp;
   if (!isSTP) {SetUpBNW();}
   if (!ballAndStickMode) {drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
   CenterMolecule();
   CalcViewRadius();
   cout << "rView: " << rView << endl;
   for (int i=0; i<3; i++) {
      pvp.locCam[i]*=rView;
      for (int j=0; j<2; j++) {
         pvp.lightSource[j][i]*=(rView*2.0e0);
      }
   }
   pvp.WriteHeader(pof);
   pof << endl;
   if (spaceFillingMode) {
      pof << HelpersPOVRay::IndTabsStr(pvp.currIndLev) << "merge{" << endl;
      pvp.currIndLev++;
   }
   PutNuclei(pof);
   if (spaceFillingMode) {
      pvp.currIndLev--;
      pof << HelpersPOVRay::IndTabsStr(pvp.currIndLev) << "}" << endl;
   }
   if (ballAndStickMode||wireMode) {
      PutBonds(pof);
   }
   pof.close();
   return true;
}
void BondNetWork::PutNuclei(ofstream & pof) {
   int atomn;
   double atrad;
   for (int i=0; i<nNuc; i++) {
      atomn=atNum[i];
      if (spaceFillingMode) {
         //atrad=atomicRadius[atomn]*AUTOMATICSPACEFILLINGRATIO;
         atrad=GetAtomicVDWRadius(atomn)*AUTOMATICSPACEFILLINGRATIO;
      } else {
         atrad=drawAtSize;
      }
      //HelpersPOVRay::WriteSphere(pof,0,R[i][0],R[i][1],R[i][2],atrad,
      //               atomicColor[atomn][0],atomicColor[atomn][1],atomicColor[atomn][2]);
      HelpersPOVRay::WriteSphere(pof,0,R[i][0],R[i][1],R[i][2],atrad,
                     GetAtomicRColorReal(atomn),GetAtomicGColorReal(atomn),
                     GetAtomicBColorReal(atomn));
   }
   return;
}
void BondNetWork::PutBonds(ofstream &pof) {
   pof << "union{" << endl;
   int k=0,atni,atnk;
   double startpt[3],frak1;
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bNet[i][j];
         atni=atNum[i];
         atnk=atNum[k];
         //frak1=atomicRadius[atni]/(atomicRadius[atni]+atomicRadius[atnk]);
         frak1=GetAtomicVDWRadius(atni)/(GetAtomicVDWRadius(atni)+GetAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=R[i][l]*(1.0e0-frak1)+R[k][l]*frak1;
         }
         if (k>0) {
            HelpersPOVRay::WriteCylinder(pof,1,
                             R[i][0],R[i][1],R[i][2],
                             startpt[0],startpt[1],startpt[2],drawStickSize,
                             GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
                             GetAtomicBColorReal(atni));
            HelpersPOVRay::WriteCylinder(pof,1,
                             startpt[0],startpt[1],startpt[2],
                             R[k][0],R[k][1],R[k][2],drawStickSize,
                             GetAtomicRColorReal(atnk),GetAtomicGColorReal(atnk),
                             GetAtomicBColorReal(atnk));
         }
      }
   }
   pof << "}" << endl;
   return;
}
void BondNetWork::CenterMolecule() {
   double trn[3];
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
void BondNetWork::CalcViewRadius(void) {
   SeekRMaxMin();
   double rmagmax=0.0e0,rmagmin=0.0e0;
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
void BondNetWork::SetBoundingBox(void) {
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<3; j++) {
         if (R[i][j]<=bbmin[j]) {bbmin[j]=R[i][j];}
         if (R[i][j]>=bbmax[j]) {bbmax[j]=R[i][j];}
      }
   }
   if (!isSTP) {
      cout << "Error: The BondNetWork object is not set up!\n";
   } else {
      for (int i=0; i<3; i++) {
         bbmax[i]+=(0.5e0*maxBondDist);
         bbmin[i]-=(0.5e0*maxBondDist);
      }
   }
   return;
}
int BondNetWork::CountAtomsOfAtomicNumber(int nat) {
#if DEBUG
   if (nat<1) {
      ScreenUtils::DisplayErrorMessage("Invalid atomic number!\nReturning zero.");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return 0;
   }
#endif
   int n=nat-1,k=0;
   for (int i=0; i<nNuc; i++) {if (atNum[i]==n) {k++;}}
   return k;
}

#endif//_BONDNETWORK_CPP_

