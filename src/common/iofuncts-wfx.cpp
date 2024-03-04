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
#ifndef _IOFUNCTS_WFX_CPP_
#define _IOFUNCTS_WFX_CPP_

#include "iofuncts-wfx.h"
#include <iostream>
using std::cout;
using std::endl;
#include "mymemory.h"
#include "stringtools.h"
#include "screenutils.h"

#define MAXWFXKEYSDEF 36
static const string wfxKeysTab[MAXWFXKEYSDEF]= {
   "Title","Keywords"
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

string GetWfxKey(const int ii) {
#if DEBUG
   if (ii>=MAXWFXKEYSDEF) {
      cout << "The maximum number of defined keys is " << MAXWFXKEYSDEF << endl
           << "Warning from: " << __FILE__ << ", line: " << __LINE__ << endl;
      return string("Unknown");
   }
#endif
   return wfxKeysTab[ii];
}
void GetPosInFile(ifstream &ifil,bool frombeg,string idkey,int & ist,int &ien) {
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
      //cout << '"' << line << '"' << std::endl;
      while ((line.size()>0)&&((line[0]==' ')||(line[0]=='\t'))) {line.erase(0,1);}
      while ((line.size()>0)&&((line.back()==' ')||(line.back()=='\t'))) {line.pop_back();}
      if ( line.size()==0 ) { continue; }
      if ( line.back()=='>' && line[0]!='<' ) {
         pos=line.find_first_of('<');
         line=line.substr(pos);
      }
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
      //cout << '"' << line << '"' << std::endl;
      while ((line.size()>0)&&((line[0]==' ')||(line[0]=='\t'))) {line.erase(0,1);}
      while ((line.size()>0)&&((line.back()==' ')||(line.back()=='\t'))) {line.pop_back();}
      if ( line.size()==0 ) { cout << "c/\n"; continue; }
      if ( line.back()=='>' && line[0]!='<' ) {
         pos=line.find_first_of('<');
         line=line.substr(pos);
      }
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
int GetInitPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey) {
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
      while ((line.size()>0)&&((line[0]==' ')||(line[0]=='\t'))) {line.erase(0,1);}
      while ((line.size()>0)&&((line.back()==' ')||(line.back()=='\t'))) {line.pop_back();}
      if ( line.size()==0 ) { continue; }
      if ( line.back()=='>' && line[0]!='<' ) {
         pos=line.find_first_of('<');
         line=line.substr(pos);
      }
      //cout << '"' << line << '"' << endl;
      if (tk==line) {
         //cout << "key found: " << line << endl;
         ifil.seekg(orpos);
         return pos;
      }
   }
   //cout << "Not found GetInitPosOfKeyInFile..." << endl;
   return -1; //if it does not return within the loop, there is an error...
}
int GetFinPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey) {
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
      while ((line.size()>0)&&((line[0]==' ')||(line[0]=='\t'))) {line.erase(0,1);}
      while ((line.size()>0)&&((line.back()==' ')||(line.back()=='\t'))) {line.pop_back();}
      if ( line.size()==0 ) { continue; }
      if ( line.back()=='>' && line[0]!='<' ) {
         pos=line.find_first_of('<');
         line=line.substr(pos);
      }
      //cout << line << endl;
      if (tk==line) {
         //cout << "key found: " << line << endl;
         ifil.seekg(orpos);
         return pos;
      }
   }
   return -1; //if the function does not return before, there is an error...
}
bool GetTitleFromFileWFX(ifstream &ifil,int &nt,string* &tit) {
   ifil.seekg(ifil.beg);
   int ipos,fpos;
   ipos=GetInitPosOfKeyInFile(ifil,true,string("Title"));
   fpos=GetFinPosOfKeyInFile(ifil,false,string("Title"));
   ifil.seekg(ipos);
   nt=0;
   string line;
   while (ifil.tellg()<(fpos-9)) {//9 is the length of the string "\n</Title>"
      getline(ifil,line);
      nt++;
   }
   MyMemory::Alloc1DStringArray("tit",nt,tit);
   ifil.seekg(ipos);
   for (int i=0; i<nt; i++) {
      getline(ifil,tit[i]);
      StringTools::RemoveSpacesLeftAndRight(tit[i]);
      StringTools::RemoveRedundantSpaces(tit[i]);
   }
   return true;
}
bool GetKeyWordsFromFileWFX(ifstream &ifil,string &kw) {
   ifil.seekg(ifil.beg);
   int ipos,fpos;
   ipos=GetInitPosOfKeyInFile(ifil,true,string("Keywords"));
   fpos=GetFinPosOfKeyInFile(ifil,false,string("Keywords"));
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
void GetNetCharge(ifstream &ifil,int &nc) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Net Charge")));
   ifil >> nc;
   return;
}
void GetNofNucleiFromFileWFX(ifstream &ifil, int & nnuc) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Number of Nuclei")));
   ifil >> nnuc;
   return;
}
void GetNofElectronsFromFileWFX(ifstream &ifil, int & nel) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Number of Electrons")));
   ifil >> nel;
   return;
}
void GetNofCoreElectronsFromFileWFX(ifstream &ifil, int & ncel) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Number of Core Electrons")));
   ifil >> ncel;
   return;
}
void GetNofMolOrbFromFileWFX(ifstream &ifil,int &nmo) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Number of Occupied Molecular Orbitals")));
   ifil >> nmo;
}
void GetNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,double* &rr) {
   if (!rr) {MyMemory::Alloc1DRealArray(string("NuclCartCoords"),3*nn,rr);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Nuclear Cartesian Coordinates")));
   int ind=0;
   for (int i=0; i<nn; i++) {
      ind=3*i;
      ifil >> rr[ind];
      ifil >> rr[ind+1];
      ifil >> rr[ind+2];
   }
   return;
}
void GetNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,double** &rr3) {
   if (!rr3) {MyMemory::Alloc2DRealArray(string("NuclCartCoords"),nn,3,rr3);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Nuclear Cartesian Coordinates")));
   for (int i=0; i<nn; i++) {
      ifil >> rr3[i][0];
      ifil >> rr3[i][1];
      ifil >> rr3[i][2];
   }
   return;
}
void GetNofPrimFromFileWFX(ifstream &ifil,int &npri) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Number of Primitives")));
   ifil >> npri;
   return;
}
void GetNofEDFPrimFromFileWFX(ifstream &ifil,const int nedfc,int &nedfprim) {
   int nc=0;
   int totp=0,tmp;
   ifil.seekg(ifil.beg);
   size_t tt;
   while ( nedfc>nc ) {
      if ( (tt=GetInitPosOfKeyInFile(ifil,false,\
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
void GetAtLabelsFromFileWFX(ifstream &ifil,const int nn,string* &al) {
   if (!al) {MyMemory::Alloc1DStringArray(string("al"),nn,al);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Nuclear Names")));
   for (int i=0; i<nn; i++) {
      getline(ifil,al[i]);
      StringTools::RemoveSpacesLeftAndRight(al[i]);
   }
   return;
}
void GetAtChargesFromFileWFX(ifstream &ifil,const int nn,double* &ach) {
   if (!ach) {MyMemory::Alloc1DRealArray(string("ach"),nn,ach);}
   //Reading atomic charges has a conflict with wfx having pseudopotentials.
   //ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Nuclear Charges")));
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Atomic Numbers")));
   for (int i=0; i<nn; i++) {ifil >> ach[i];}
   return;
}
int GetAtNumbersFromFileWFX(ifstream &ifil,const int nn,int* &anu) {
   if(!anu){MyMemory::Alloc1DIntArray(string("anu"),nn,anu);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Atomic Numbers")));
   int nprot=0;
   for (int i=0; i<nn; i++) {
      ifil >> anu[i];
      nprot+=anu[i];
      anu[i]--;
   }
   return nprot;
}
void GetPrimCentersFromFileWFX(ifstream &ifil,const int npr,int* &pc) {
   if (!pc) {MyMemory::Alloc1DIntArray(string("pc"),npr,pc);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Primitive Centers")));
   for (int i=0; i<npr; i++) {
      ifil >> pc[i];
      pc[i]--;
   }
   return;
}
void GetEDFPrimCentersFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pc) {
   if (pc==NULL) {MyMemory::Alloc1DIntArray(string("pc"),ntot,pc);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("EDF Primitive Centers")));
   for (int i=npr; i<ntot; ++i) {
      ifil >> pc[i];
      pc[i]--;
   }
   return;
}
void GetPrimTypesFromFileWFX(ifstream &ifil,const int npr,int* &pt) {
   if (!pt) {MyMemory::Alloc1DIntArray(string("pt"),npr,pt);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Primitive Types")));
   for (int i=0; i<npr; i++) {
      ifil >> pt[i];
      pt[i]--;
   }
   return;
}
void GetEDFPrimTypesFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pt) {
   if (pt==NULL) {MyMemory::Alloc1DIntArray(string("pt"),ntot,pt);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("EDF Primitive Types")));
   for (int i=npr; i<ntot; ++i) {
      ifil >> pt[i];
      pt[i]--;
   }
   return;
}
void GetPrimExponentsFromFileWFX(ifstream &ifil,const int npr,double* &pex) {
   if (!pex) {MyMemory::Alloc1DRealArray(string("pex"),npr,pex);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Primitive Exponents")));
   for (int i=0; i<npr; i++) {ifil >> pex[i];}
   return;
}
void GetEDFPrimExponentsFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,double* &pex) {
   if (pex==NULL) {MyMemory::Alloc1DRealArray(string("pex"),ntot,pex);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("EDF Primitive Exponents")));
   for (int i=npr; i<ntot; ++i) {ifil >> pex[i];}
   return;
}
void GetMolecOrbOccNumsFromFileWFX(ifstream &ifil,const int nmo,double* &ocnu) {
   if (!ocnu) {MyMemory::Alloc1DRealArray(string("ocnu"),nmo,ocnu);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Occupation Numbers")));
   for (int i=0; i<nmo; i++) {ifil >> ocnu[i];}
   return;
}
void GetMolecOrbEnergiesFromFileWFX(ifstream &ifil,const int nmo,double* &orben) {
   if (!orben) {MyMemory::Alloc1DRealArray(string("orben"),nmo,orben);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Energies")));
   for (int i=0; i<nmo; i++) {ifil >> orben[i];}
   return;
}
bool GetMolecOrbSpinTypesFromFileWFX(ifstream &ifil,const int nmo,int* &orbsptp) {
   if (!orbsptp) {MyMemory::Alloc1DIntArray(string("orbsptp"),nmo,orbsptp);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Spin Types")));
   bool res=false;
   bool dispmsg=true;
   string spstr;
   for (int i=0; i<nmo; i++) {
      getline(ifil,spstr);
      while ((spstr.size()>0)&&((spstr[0]==' ')||(spstr[0]=='\t'))) {spstr.erase(0,1);}
      while ((spstr.size()>0)&&((spstr.back()==' ')||(spstr.back()=='\t'))) {spstr.pop_back();}
      if ( spstr==string("Alpha and Beta") ) {
         orbsptp[i]=0;
      } else if ( spstr==string("Alpha") ) {
         orbsptp[i]=1;
         res=true;
      } else if ( spstr==string("Beta") ) {
         orbsptp[i]=2;
         res=true;
      } else {
         orbsptp[i]=0;
         if ( dispmsg ) {
            ScreenUtils::DisplayWarningMessage(string("Spin type '")+spstr
                  +string("' not recognized. Setting MO ")+std::to_string(i)
                  +string(" as 'Alpha and Beta'."));
            cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
            dispmsg=false;
         }
      }
   }
   return res;
}
void GetMolecOrbCoefficientsFromFileWFX(ifstream &ifil,const int nmo,const int npr,double* &tcf) {
   if (!tcf) {MyMemory::Alloc1DRealArray(string("tcf"),(nmo*npr),tcf);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Molecular Orbital Primitive Coefficients")));
   for (int i=0; i<nmo; i++) {
      ifil.seekg(GetFinPosOfKeyInFile(ifil,false,string("MO Number")));
      for (int j=0; j<npr; j++) {
         ifil >> tcf[npr*i+j];
      }
   }
   return;
}
void GetEDFPrimCoefficientsFromFileWFX(ifstream &ifil,const int nedfp,\
      double* &edfc) {
   if (edfc==NULL) {MyMemory::Alloc1DRealArray(string("edfc"),nedfp,edfc);}
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,\
            string("EDF Primitive Coefficients")));
   for (int i=0; i<nedfp; ++i) {ifil >> edfc[i];}
   return;
}
void GetTotEnerAndVirialFromFileWFX(ifstream &ifil,double &tote,double &vir) {
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Energy = T + Vne + Vee + Vnn")));
   ifil >> tote;
   ifil.seekg(GetInitPosOfKeyInFile(ifil,true,string("Virial Ratio (-V/T)")));
   ifil >> vir;
   return;
}
void countEDFCentersFromFileWFX(ifstream &ifil,int &nedfc) {
   string line;
   int nn=0;
   size_t orpos=ifil.tellg();
   ifil.seekg(ifil.beg);
   while ( !ifil.eof() ) {
      getline(ifil,line);
      StringTools::RemoveSpacesLeftAndRight(line);
      if ( line==string("<EDF Name>") ) {++nn;}
   }
   nedfc=nn;
   ifil.seekg(orpos);
   ifil.clear();
}
void GetEDFExistenceFromFileWFX(ifstream &ifil,bool &ihaveEDF) {
   string line;
   size_t orpos=ifil.tellg();
   ifil.seekg(ifil.beg);
   while ( !ifil.eof() ) {
      getline(ifil,line);
      StringTools::RemoveSpacesLeftAndRight(line);
      if ( line==string("<Number of EDF Primitives>") ) {
         ihaveEDF=true;
         break;
      }
   }
   if ( ifil.eof() ) { ifil.clear(); }
   ifil.seekg(orpos);
}
#endif//_IOFUNCTS_WFX_CPP_

