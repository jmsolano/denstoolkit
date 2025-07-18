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



#ifndef _IOFUNCTS_WFN_CPP_
#define _IOFUNCTS_WFN_CPP_
#include <cstdlib>
using std::atof;
#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
#include "mymemory.h"
#include "iofuncts-wfn.h"

string GetTitleFromFileWFN(ifstream &ifil) {
   string line;
   ifil.seekg(ifil.beg);
   getline(ifil,line);
   return line;
}
void ProcessFirstDataStringinWFNFile(ifstream &ifil,string* &tit,string &orbdesc,int &nmo,int &npr,int &nnu) {
   string line;
   ifil.seekg(ifil.beg);
   MyMemory::Alloc1DStringArray(string("tit"),1,tit);
   getline(ifil,tit[0]);
   getline(ifil,line);
   while (line[0]==' ') {line.erase(0,1);}
   int pos;
   pos=line.find_first_of(' ');
   orbdesc=line.substr(0,pos);
   if (line.substr(0,pos)!="GAUSSIAN") {
      cout << "Error: only gaussian orbitals are implemented in this version...\n";
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      ifil.close();
      exit(1);
   }
   line.erase(0,pos);
   while (line[0]==' ') {line.erase(0,1);}
   sscanf(line.c_str(),"%d",&nmo);
   pos=line.find_first_of(' ');
   line.erase(0,pos);
   while (line[0]==' ') {line.erase(0,1);}
   if (line.substr(0,12)!="MOL ORBITALS") {
      cout << "Error: Expecting \"MOL ORBITALS\" in file...\n";
      ifil.close();
      exit(1);
   }
   line.erase(0,12);
   while (line[0]==' ') {line.erase(0,1);}
   sscanf(line.c_str(),"%d",&npr);
   pos=line.find_first_of(' ');
   line.erase(0,pos);
   while (line[0]==' ') {line.erase(0,1);}
   if (line.substr(0,10)!="PRIMITIVES") {
      cout << "Error: Expecting \"PRIMITIVES\" in file...\n";
      ifil.close();
      exit(1);
   }
   line.erase(0,10);
   while (line[0]==' ') {line.erase(0,1);}
   sscanf(line.c_str(),"%d",&nnu);
   pos=line.find_first_of(' ');
   line.erase(0,pos);
   while (line[0]==' ') {line.erase(0,1);}
   if (line.substr(0,6)!="NUCLEI") {
      cout << "Error: Expecting \"NUCLEI\" in file...\n";
      ifil.close();
      exit(1);
   }
}
void ProcessCentersWFN(ifstream &ifil,const int nnu,string* &atlbl,double* &rr,double* &atch) {
   MyMemory::Alloc1DStringArray("atlbl",nnu,atlbl);
   MyMemory::Alloc1DRealArray("rr",3*nnu,rr);
   MyMemory::Alloc1DRealArray("atch",nnu,atch);
   string line;
   int fipos,fepos;
   size_t pos;
   //bool isspace=true;
   for (int i=0; i<nnu; i++) {
      fipos=ifil.tellg();
      getline(ifil,line);
      fepos=ifil.tellg();
      atlbl[i]=line.substr(0,12);
      while (atlbl[i][0]==' ') {atlbl[i].erase(0,1);}
      while (atlbl[i][atlbl[i].length()-1]==' ') {
         atlbl[i].erase(atlbl[i].length()-1);
      }
      pos=atlbl[i].find_first_of(' ');
      while (pos!=string::npos) {
         atlbl[i].erase(pos,1);
         pos=atlbl[i].find_first_of(' ');
      }
      line.erase(0,23);
      //sscanf(line.c_str(),"%g %g %g",&rr[3*i],&rr[3*i+1],&rr[3*i+2]);
      ifil.seekg(fipos+23);
      ifil >> rr[3*i] >> rr[3*i+1] >> rr[3*i+2];
      //cout << atlbl[i] << endl;
      //cout << "R[" << i << "]: " << rr[3*i] << " " << rr[3*i+1] << " " << rr[3*i+2] << endl;
      ifil.seekg(fipos+70);
      ifil >> atch[i];
      //cout << "Charge[" << i << "]: " << atch[i] << endl;
      ifil.seekg(fepos);
   }
}
void ProcessPrimitivesWFN(ifstream &ifil,const int npr,int* &pricen,int* &primty,double* &prexp) {
   MyMemory::Alloc1DIntArray("pricen",npr,pricen);
   MyMemory::Alloc1DIntArray("primty",npr,primty);
   MyMemory::Alloc1DRealArray("prexp",npr,prexp);
   int count=0;
   int pos;
   /* Getting the primitive centers */
   pos=ifil.tellg();
   ifil.seekg(pos+20);
   for (int i=0; i<npr; i++) {
      ifil >> pricen[i];
      pricen[i]--;
      count++;
      if (count==20) {
         pos=ifil.tellg();
         ifil.seekg(pos+20);
         count=0;
      }
      //cout << pricen[i] << " ";
   }
   if ((npr%20)==0) {
      pos=ifil.tellg();
      ifil.seekg(pos-20);
   }
   /* Getting the primitive types */
   pos=ifil.tellg();
   ifil.seekg(pos+20);
   count=0;
   for (int i=0; i<npr; i++) {
      ifil >> primty[i];
      primty[i]--;
      count++;
      if (count==20) {
         pos=ifil.tellg();
         ifil.seekg(pos+20);
         count=0;
      }
      //cout << primty[i] << " ";
   }
   if ((npr%20)==0) {
      pos=ifil.tellg();
      ifil.seekg(pos-20);
   }
   /* Getting the primitive exponents */
   string line;
   int nolines;
   nolines=floor(npr/5);
   if ((npr%5)>0) {nolines++;}
   count=0;
   for (int i=0; i<nolines; i++) {
      pos=ifil.tellg();
      getline(ifil,line); 
      if (line.length()==0) {getline(ifil,line); pos++;}
      line.erase(0,10);
      for (int j=0; j<5; j++) {
         if (count<npr) {
            if ((line[10]=='D')||(line[10]=='d')) {
               line[10]='e';
            }
            prexp[count]=(double)atof((line.substr(0,14)).c_str());
            count++;
            line.erase(0,14);
         }
      }
   }
}
void ProcessMolecularOrbitalPropsAndCoefs(ifstream &ifil,const int norb,const int npr
                                 ,double* &ocn,double* &moe,double* &moc) {
   MyMemory::Alloc1DRealArray("moc",(norb*npr),moc);
   MyMemory::Alloc1DRealArray("ocn",norb,ocn);
   MyMemory::Alloc1DRealArray("moe",norb,moe);
   int pos;
   /* Getting the molecular orbital occupation number, energies, and coefficients */
   string line;
   int count=0;
   int nolines;
   nolines=floor(npr/5);
   //cout << "npr: " << npr << " nolines: " << npr/5 << endl;
   if ((npr%5)>0) {nolines++;}
   for (int k=0; k<norb; k++) {
      pos=ifil.tellg();
      ifil.seekg(pos+35);
      ifil >> ocn[k];
      ifil.seekg(pos+62);
      ifil >> moe[k];
      ifil.seekg(pos+75);
      count=0;
      for (int i=0; i<nolines; i++) {
         pos=ifil.tellg();
         getline(ifil,line);
         if (line.length()==0) {getline(ifil,line); pos++;}
         for (int j=0; j<5; j++) {
            if (count<npr) {
               if ((line[12]=='D')||(line[12]=='d')) {
                  line[12]='e';
               }
               moc[(k*npr)+count]=(double)atof((line.substr(0,16)).c_str());
               count++;
               line.erase(0,16);
            }
         }
      }
   }
}
void GetEnergyAndVirial(ifstream &ifil,double &theener,double &thevir) {
   string line;
   getline(ifil,line);
   int pos;
   pos=line.find_first_of('=');
   line.erase(0,pos+1);
   while (line[0]==' ') {line.erase(0,1);}
   pos=line.find_first_of(' ');
   theener=(double)atof((line.substr(0,pos)).c_str());
   pos=line.find_first_of('=');
   line.erase(0,pos+1);
   while (line[0]==' ') {line.erase(0,1);}
   thevir=(double)atof(line.c_str());
}
#endif//_IOFUNCTS_WFN_CPP_

