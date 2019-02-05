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



#ifndef _IOFUNCTS_WFX_CPP_
#define _IOFUNCTS_WFX_CPP_

#include "iofuncts-wfx.h"
#include <iostream>
using std::cout;
using std::endl;
#include "solmemhand.h"
#include "solstringtools.h"

#define MAXWFXKEYSDEF 36
static const string wfxKeysTab[MAXWFXKEYSDEF]=
{"Title","Keywords"
   ,"Number of Nuclei","Number of Primitives"
   ,"Number of Occupied Molecular Orbitals","Number of Perturbations"
   ,"Nuclear Names","Atomic Numbers"
   ,"Nuclear Charges","Nuclear Cartesian Coordinates"
   ,"Net Charge","Number of Electrons"
   ,"Number of Alpha Electrons","Number of Beta Electrons"
   ,"Electronic Spin Multiplicity","Model"
   ,"Primitive Centers","Primitive Types"
   ,"Primitive Exponents","Molecular Orbital Occupation Numbers"
   ,"Molecular Orbital Energies","Molecular Orbital Spin Types"
   ,"Molecular Orbital Primitive Coefficients","MO Number"
   ,"Energy = T + Vne + Vee + Vnn","Virial Ratio (-V/T)"
   ,"Nuclear Cartesian Energy Gradients","Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W"
   ,"Full Virial Ratio, -(V - W)/T","Number of Core Electrons"
   ,"Additional Electron Density Function (EDF)","Number of EDF Primitives"
   ,"EDF Primitive Centers","EDF Primitive Types"
   ,"EDF Primitive Exponents","EDF Primitive Coefficients"
};

/* ************************************************************************** */
string getWfxKey(const int ii)
{
#if DEBUG
   if (ii>=MAXWFXKEYSDEF) {
      cout << "The maximum number of defined keys is " << MAXWFXKEYSDEF << endl
           << "Warning from: " << __FILE__ << ", line: " << __LINE__ << endl;
      return string("Unknown");
   }
#endif
   return wfxKeysTab[ii];
}
/* ************************************************************************** */
void getPosInFile(ifstream &ifil,bool frombeg,string idkey,int & ist,int &ien)
{
   int orpos,pos;
   string line,tk;
   orpos=ifil.tellg();
#if DEBUG
   bool havep1;
   havep1=false;
   for (int i=0; i<MAXWFXKEYSDEF; i++) {
      if (idkey==wfxKeysTab[i]) {
         havep1=true;
         break;
      }
   }
   if (!havep1) {
      cout << "Not a valid key!\n" << __FILE__ << ", " << __LINE__ << endl;
   }
   havep1=false;
#endif
   if (frombeg) {ifil.seekg(ifil.beg);}
   tk=idkey;
   tk.append(">"); tk.insert(0,"<");
   //cout <<  "line= " << line << endl;
   //*
   while (!ifil.eof()) {
      getline(ifil,line);
      pos=ifil.tellg();
      //cout << line << endl;
      if (tk==line) {
         //cout << "key found: " << line << endl;
         ist=pos;
         break;
      }
   }
   tk.insert(1,"/");
   while (!ifil.eof()) {
      pos=ifil.tellg();
      getline(ifil,line);
      //cout << line << endl;
      if (tk==line) {
         //havep1=true;
         //cout << "key found: " << line << endl;
         ien=pos;
         break;
      }
   }
   if (ifil.eof()) {
      cout << "End of file and no string found\n";
   } else {
      ifil.seekg(orpos);
   }
   // */
}
/* ************************************************************************** */
int getInitPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey)
{
   int orpos,pos;
   string line,tk;
   orpos=ifil.tellg();
#if DEBUG
   bool havep1;
   havep1=false;
   for (int i=0; i<MAXWFXKEYSDEF; i++) {
      if (idkey==wfxKeysTab[i]) {
         havep1=true;
         break;
      }
   }
   if (!havep1) {
      cout << "Not a valid key!\n" << __FILE__ << ", " << __LINE__ << endl;
   }
   havep1=false;
#endif
   if (frombeg) {ifil.seekg(ifil.beg);}
   tk=idkey;
   tk.append(">"); tk.insert(0,"<");
   //*
   while (!ifil.eof()) {
      getline(ifil,line);
      pos=ifil.tellg();
      while ((line[0]==' ')||(line[0]=='\t')) {line.erase(0,1);}
      //cout << line << endl;
      if (tk==line) {
         //cout << "key found: " << line << endl;
         ifil.seekg(orpos);
         return pos;
      }
   }
   //cout << "Not found getInitPosOfKeyInFile..." << endl;
   return -1; //if it does not return within the loop, there is an error...
}
/* ************************************************************************** */
int getFinPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey)
{
   int orpos,pos;
   string line,tk;
   orpos=ifil.tellg();
#if DEBUG
   bool havep1;
   havep1=false;
   for (int i=0; i<MAXWFXKEYSDEF; i++) {
      if (idkey==wfxKeysTab[i]) {
         havep1=true;
         break;
      }
   }
   if (!havep1) {
      cout << "Not a valid key!\n" << __FILE__ << ", " << __LINE__ << endl;
   }
   havep1=false;
#endif
   if (frombeg) {ifil.seekg(ifil.beg);}
   tk=idkey;
   tk.append(">"); tk.insert(0,"</");
   //*
   while (!ifil.eof()) {
      getline(ifil,line);
      pos=ifil.tellg();
      while ((line[0]==' ')||(line[0]=='\t')) {line.erase(0,1);}
      //cout << line << endl;
      if (tk==line) {
         //cout << "key found: " << line << endl;
         ifil.seekg(orpos);
         return pos;
      }
   }
   return -1; //if the function does not return before, there is an error...
}
/* ************************************************************************** */
bool getTitleFromFileWFX(ifstream &ifil,int &nt,string* &tit)
{
   ifil.seekg(ifil.beg);
   int ipos,fpos;
   ipos=getInitPosOfKeyInFile(ifil,true,string("Title"));
   fpos=getFinPosOfKeyInFile(ifil,false,string("Title"));
   ifil.seekg(ipos);
   nt=0;
   string line;
   while (ifil.tellg()<(fpos-9)) {//9 is the length of the string "\n</Title>"
      getline(ifil,line);
      nt++;
   }
   alloc1DStringArray("tit",nt,tit);
   ifil.seekg(ipos);
   for (int i=0; i<nt; i++) {
      getline(ifil,tit[i]);
      removeSpacesLeftAndRight(tit[i]);
      removeRedundantSpaces(tit[i]);
   }
   return true;
}
/* ************************************************************************** */
bool getKeyWordsFromFileWFX(ifstream &ifil,string &kw)
{
   ifil.seekg(ifil.beg);
   int ipos,fpos;
   ipos=getInitPosOfKeyInFile(ifil,true,string("Keywords"));
   fpos=getFinPosOfKeyInFile(ifil,false,string("Keywords"));
   ifil.seekg(ipos);
   string line;
   while (ifil.tellg()<(fpos-12)) {//12 is the length of the string "\n</Keywords>"
      getline(ifil,line);
      kw+=line;
      kw+=" ";
   }
   while (kw[0]==' ') {kw.erase(0,1);}
   //cout << kw;
   return true;
}
/* ************************************************************************** */
void getNetCharge(ifstream &ifil,int &nc)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Net Charge")));
   ifil >> nc;
   return;
}
/* ************************************************************************** */
void getNofNucleiFromFileWFX(ifstream &ifil, int & nnuc)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Number of Nuclei")));
   ifil >> nnuc;
   return;
}
/* ************************************************************************** */
void getNofElectronsFromFileWFX(ifstream &ifil, int & nel)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Number of Electrons")));
   ifil >> nel;
   return;
}
/* ************************************************************************** */
void getNofCoreElectronsFromFileWFX(ifstream &ifil, int & ncel)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Number of Core Electrons")));
   ifil >> ncel;
   return;
}
/* ************************************************************************** */
void getNofMolOrbFromFileWFX(ifstream &ifil,int &nmo)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Number of Occupied Molecular Orbitals")));
   ifil >> nmo;
}
/* ************************************************************************** */
void getNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,solreal* &rr)
{
   if (!rr) {alloc1DRealArray(string("NuclCartCoords"),3*nn,rr);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Nuclear Cartesian Coordinates")));
   int ind=0;
   for (int i=0; i<nn; i++) {
      ind=3*i;
      ifil >> rr[ind];
      ifil >> rr[ind+1];
      ifil >> rr[ind+2];
   }
   return;
}
/* ************************************************************************** */
void getNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,solreal** &rr3)
{
   if (!rr3) {alloc2DRealArray(string("NuclCartCoords"),nn,3,rr3);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Nuclear Cartesian Coordinates")));
   for (int i=0; i<nn; i++) {
      ifil >> rr3[i][0];
      ifil >> rr3[i][1];
      ifil >> rr3[i][2];
   }
   return;
}
/* ************************************************************************** */
void getNofPrimFromFileWFX(ifstream &ifil,int &npri)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Number of Primitives")));
   ifil >> npri;
   return;
}
/* ************************************************************************** */
void getNofEDFPrimFromFileWFX(ifstream &ifil,const int nedfc,int &nedfprim)
{
   int nc=0;
   int totp=0,tmp;
   ifil.seekg(ifil.beg);
   size_t tt;
   while ( nedfc>nc ) {
      if ( (tt=getInitPosOfKeyInFile(ifil,false,\
                  string("Number of EDF Primitives"))) == string::npos ) {
         cout << "Error: key not found!" << endl;
         break;
      }
      ifil.seekg(tt);
      ifil >> tmp;
      totp+=tmp;
      nc++;
   }
   nedfprim=totp;
   ifil.clear();
   return;
}
/* ************************************************************************** */
void getAtLabelsFromFileWFX(ifstream &ifil,const int nn,string* &al)
{
   if (!al) {alloc1DStringArray(string("al"),nn,al);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Nuclear Names")));
   for (int i=0; i<nn; i++) {
      getline(ifil,al[i]);
      removeSpacesLeftAndRight(al[i]);
   }
   return;
}
/* ************************************************************************** */
void getAtChargesFromFileWFX(ifstream &ifil,const int nn,solreal* &ach)
{
   if (!ach) {alloc1DRealArray(string("ach"),nn,ach);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Nuclear Charges")));
   for (int i=0; i<nn; i++) {ifil >> ach[i];}
   return;
}
/* ************************************************************************** */
int getAtNumbersFromFileWFX(ifstream &ifil,const int nn,int* &anu)
{
   if(!anu){alloc1DIntArray(string("anu"),nn,anu);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Atomic Numbers")));
   int nprot=0;
   for (int i=0; i<nn; i++) {
      ifil >> anu[i];
      nprot+=anu[i];
      anu[i]--;
   }
   return nprot;
}
/* ************************************************************************** */
void getPrimCentersFromFileWFX(ifstream &ifil,const int npr,int* &pc)
{
   if (!pc) {alloc1DIntArray(string("pc"),npr,pc);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Primitive Centers")));
   for (int i=0; i<npr; i++) {
      ifil >> pc[i];
      pc[i]--;
   }
   return;
}
/* ************************************************************************** */
void getEDFPrimCentersFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pc)
{
   if (pc==NULL) {alloc1DIntArray(string("pc"),ntot,pc);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("EDF Primitive Centers")));
   for (int i=npr; i<ntot; ++i) {
      ifil >> pc[i];
      pc[i]--;
   }
   return;
}
/* ************************************************************************** */
void getPrimTypesFromFileWFX(ifstream &ifil,const int npr,int* &pt)
{
   if (!pt) {alloc1DIntArray(string("pt"),npr,pt);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Primitive Types")));
   for (int i=0; i<npr; i++) {
      ifil >> pt[i];
      pt[i]--;
   }
   return;
}
/* ************************************************************************** */
void getEDFPrimTypesFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pt)
{
   if (pt==NULL) {alloc1DIntArray(string("pt"),ntot,pt);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("EDF Primitive Types")));
   for (int i=npr; i<ntot; ++i) {
      ifil >> pt[i];
      pt[i]--;
   }
   return;
}
/* ************************************************************************** */
void getPrimExponentsFromFileWFX(ifstream &ifil,const int npr,solreal* &pex)
{
   if (!pex) {alloc1DRealArray(string("pex"),npr,pex);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Primitive Exponents")));
   for (int i=0; i<npr; i++) {ifil >> pex[i];}
   return;
}
/* ************************************************************************** */
void getEDFPrimExponentsFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,solreal* &pex)
{
   if (pex==NULL) {alloc1DRealArray(string("pex"),ntot,pex);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("EDF Primitive Exponents")));
   for (int i=npr; i<ntot; ++i) {ifil >> pex[i];}
   return;
}
/* ************************************************************************** */
void getMolecOrbOccNumsFromFileWFX(ifstream &ifil,const int nmo,solreal* &ocnu)
{
   if (!ocnu) {alloc1DRealArray(string("ocnu"),nmo,ocnu);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Occupation Numbers")));
   for (int i=0; i<nmo; i++) {ifil >> ocnu[i];}
   return;
}
/* ************************************************************************** */
void getMolecOrbEnergiesFromFileWFX(ifstream &ifil,const int nmo,solreal* &orben)
{
   if (!orben) {alloc1DRealArray(string("orben"),nmo,orben);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Energies")));
   for (int i=0; i<nmo; i++) {ifil >> orben[i];}
   return;
}
/* ************************************************************************** */
void getMolecOrbCoefficientsFromFileWFX(ifstream &ifil,const int nmo,const int npr,solreal* &tcf)
{
   if (!tcf) {alloc1DRealArray(string("tcf"),(nmo*npr),tcf);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Primitive Coefficients")));
   for (int i=0; i<nmo; i++) {
      ifil.seekg(getFinPosOfKeyInFile(ifil,false,string("MO Number")));
      for (int j=0; j<npr; j++) {
         ifil >> tcf[npr*i+j];
      }
   }
   return;
}
/* ************************************************************************** */
void getEDFPrimCoefficientsFromFileWFX(ifstream &ifil,const int nedfp,\
      solreal* &edfc)
{
   if (edfc==NULL) {alloc1DRealArray(string("edfc"),nedfp,edfc);}
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,\
            string("EDF Primitive Coefficients")));
   for (int i=0; i<nedfp; ++i) {ifil >> edfc[i];}
   return;
}
/* ************************************************************************** */
void getTotEnerAndVirialFromFileWFX(ifstream &ifil,solreal &tote,solreal &vir)
{
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Energy = T + Vne + Vee + Vnn")));
   ifil >> tote;
   ifil.seekg(getInitPosOfKeyInFile(ifil,true,string("Virial Ratio (-V/T)")));
   ifil >> vir;
   return;
}
/* ************************************************************************** */
void countEDFCentersFromFileWFX(ifstream &ifil,int &nedfc)
{
   string line;
   int nn=0;
   size_t orpos=ifil.tellg();
   ifil.seekg(ifil.beg);
   while ( !ifil.eof() ) {
      getline(ifil,line);
      removeSpacesLeftAndRight(line);
      if ( line==string("<EDF Name>") ) {++nn;}
   }
   nedfc=nn;
   ifil.seekg(orpos);
   ifil.clear();
}
/* ************************************************************************** */
void getEDFExistenceFromFileWFX(ifstream &ifil,bool &ihaveEDF)
{
   string line;
   size_t orpos=ifil.tellg();
   ifil.seekg(ifil.beg);
   while ( !ifil.eof() ) {
      getline(ifil,line);
      removeSpacesLeftAndRight(line);
      if ( line==string("<Number of EDF Primitives>") ) {
         ihaveEDF=true;
         break;
      }
   }
   if ( ifil.eof() ) { ifil.clear(); }
   ifil.seekg(orpos);
}
/* ************************************************************************** */
/* ************************************************************************** */
#endif//_IOFUNCTS_WFX_CPP_

