

#ifndef _IOFUNCTS_WFN_CPP_
#define _IOFUNCTS_WFN_CPP_

#include "iofuncts-wfn.h"
#include <cstdlib>
using std::atof;
#include <iostream>
using std::cout;
using std::endl;
#include "solmemhand.h"
#include <cmath>
/* ********************************************************************************************* */
/* ********************************************************************************************* */
string getTitleFromFileWFN(ifstream &ifil)
{
   string line;
   ifil.seekg(ifil.beg);
   getline(ifil,line);
   return line;
}
/* ********************************************************************************************* */
void processFirstDataStringinWFNFile(ifstream &ifil,string* &tit,string &orbdesc,int &nmo,int &npr,int &nnu)
{
   string line;
   ifil.seekg(ifil.beg);
   alloc1DStringArray(string("tit"),1,tit);
   getline(ifil,tit[0]);
   getline(ifil,line);
   while (line[0]==' ') {line.erase(0,1);}
   int pos;
   pos=line.find_first_of(' ');
   orbdesc=line.substr(0,pos);
   if (line.substr(0,pos)!="GAUSSIAN") {
      cout << "Error: only gaussian orbitals are implemented in this version...\n";
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
/* ********************************************************************************************* */
void processCentersWFN(ifstream &ifil,const int nnu,string* &atlbl,solreal* &rr,solreal* &atch)
{
   alloc1DStringArray("atlbl",nnu,atlbl);
   alloc1DRealArray("rr",3*nnu,rr);
   alloc1DRealArray("atch",nnu,atch);
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
/* ********************************************************************************************* */
void processPrimitivesWFN(ifstream &ifil,const int npr,int* &pricen,int* &primty,solreal* &prexp)
{
   alloc1DIntArray("pricen",npr,pricen);
   alloc1DIntArray("primty",npr,primty);
   alloc1DRealArray("prexp",npr,prexp);
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
            prexp[count]=(solreal)atof((line.substr(0,14)).c_str());
            count++;
            line.erase(0,14);
         }
      }
   }
}
/* ********************************************************************************************* */
void processMolecularOrbitalPropsAndCoefs(ifstream &ifil,const int norb,const int npr
                                 ,solreal* &ocn,solreal* &moe,solreal* &moc)
{
   alloc1DRealArray("moc",(norb*npr),moc);
   alloc1DRealArray("ocn",norb,ocn);
   alloc1DRealArray("moe",norb,moe);
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
      ifil.seekg(pos+63);
      ifil >> moe[k];
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
               moc[(k*npr)+count]=(solreal)atof((line.substr(0,16)).c_str());
               count++;
               line.erase(0,16);
            }
         }
      }
   }
}
/* ********************************************************************************************* */
void getEnergyAndVirial(ifstream &ifil,solreal &theener,solreal &thevir)
{
   string line;
   getline(ifil,line);
   int pos;
   pos=line.find_first_of('=');
   line.erase(0,pos+1);
   while (line[0]==' ') {line.erase(0,1);}
   pos=line.find_first_of(' ');
   theener=(solreal)atof((line.substr(0,pos)).c_str());
   pos=line.find_first_of('=');
   line.erase(0,pos+1);
   while (line[0]==' ') {line.erase(0,1);}
   thevir=(solreal)atof(line.c_str());
}
/* ********************************************************************************************* */
/* ********************************************************************************************* */
/* ********************************************************************************************* */
#endif//_IOFUNCTS_WFN_CPP_

