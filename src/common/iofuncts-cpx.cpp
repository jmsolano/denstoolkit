//
//  iofuncts-cpx.cpp
//  
//
//  Created by Juan Manuel Solano on 2013-10-28.
//
//

#ifndef _IOFUNCTS_CPX_CPP_
#define _IOFUNCTS_CPX_CPP_

#include "iofuncts-cpx.h"
#include "solscrutils.h"
#include "solstringtools.h"

//**************************************************************************************************
#define MAXCPXKEYSDEFINED 24
static const string cpxKeysTab[MAXCPXKEYSDEFINED]={
   "WaveFunctionFileName","CriticalPointType",\
   "NumberOfCriticalPoints","NumberOfACPs",\
   "NumberOfBCPs","NumberOfRCPs",\
   "NumberOfCCPs","ACPCartesianCoordinates",\
   "BCPCartesianCoordinates","RCPCartesianCoordinates",\
   "CCPCartesianCoordinates","ACPLabels",\
   "BCPLabels","RCPLabels",\
   "CCPLabels","NumberOfBondPaths",\
   "NumbersOfPointsPerBondPath","BondPathsData",\
   "BondPathIndex","CoordinatesOfBondPathPoints",\
   "ACPConnectivity","BCPConnectivity",\
   "RCPConnectivity","CCPConnectivity"
};
//**************************************************************************************************
void writeCPXFile(string cpxname,string wfname,critPtNetWork &cp)
{
   ofstream cpxfil;
   cpxfil.open(cpxname.c_str(),ios::out);
   writeWFFileName(cpxfil,wfname);
   writeTypeOfCriticalPoints(cpxfil,cp);
   writeNumberOfCriticalPoints(cpxfil,cp);
   writeCoordinatesACPs(cpxfil,cp);
   writeCoordinatesBCPs(cpxfil,cp);
   writeCoordinatesRCPs(cpxfil,cp);
   writeCoordinatesCCPs(cpxfil,cp);
   writeConnectivityBCPs(cpxfil,cp);
   writeLabelsACPs(cpxfil,cp);
   writeLabelsBCPs(cpxfil,cp);
   writeLabelsRCPs(cpxfil,cp);
   writeLabelsCCPs(cpxfil,cp);
   writeNumberOfBondPaths(cpxfil,cp);
   writeNumberOfPointsPerBondPath(cpxfil,cp);
   writeBondPathsCoordinates(cpxfil,cp);
   cpxfil.close();
   return;
}
//**************************************************************************************************
void writeOpenningAttribute(ofstream &ofil,const char* attrib,bool nl = true)
{
   ofil << "<" << attrib << ">";
   if (nl) {ofil << endl << " ";}
   return;
}
//**************************************************************************************************
void writeClosingAttribute(ofstream &ofil,const char* attrib,bool nl = true)
{
   ofil << "</" << attrib << ">";
   if (nl) {ofil << endl;}
   return;
}
//**************************************************************************************************
void writeWFFileName(ofstream &ofil,string wfname)
{
   string key="WaveFunctionFileName";
   writeOpenningAttribute(ofil,key.c_str());
   ofil << wfname << endl;
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeTypeOfCriticalPoints(ofstream &ofil,critPtNetWork &cp)
{
   string key="CriticalPointType";
   writeOpenningAttribute(ofil,key.c_str());
   switch (cp.myCPType()) {
      case DENS:
         ofil << "Electron Density";
         break;
      case LOLD:
         ofil << "LOL";
         break;
      default:
         ofil << "Unknown";
#if DEBUG
         displayErrorMessage("Unknown Field Type!");
#endif
         break;
   }
   ofil << endl;
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeNumberOfCriticalPoints(ofstream &ofil,critPtNetWork &cp)
{
   string mkey="NumberOfCriticalPoints";
   string akey="NumberOfACPs";
   string bkey="NumberOfBCPs";
   string rkey="NumberOfRCPs";
   string ckey="NumberOfCCPs";
   //
   writeOpenningAttribute(ofil,mkey.c_str(),false);
   ofil << endl;
   //
   writeOpenningAttribute(ofil,akey.c_str());
   if (cp.iKnowACPs()) {ofil << cp.nACP;} else {ofil << -1;}
   ofil << endl;
   writeClosingAttribute(ofil,akey.c_str());
   //
   writeOpenningAttribute(ofil,bkey.c_str());
   if (cp.iKnowBCPs()) {ofil << cp.nBCP;} else {ofil << -1;}
   ofil << endl;
   writeClosingAttribute(ofil,bkey.c_str());
   //
   writeOpenningAttribute(ofil,rkey.c_str());
   if (cp.iKnowRCPs()) {ofil << cp.nRCP;} else {ofil << -1;}
   ofil << endl;
   writeClosingAttribute(ofil,rkey.c_str());
   //
   writeOpenningAttribute(ofil,ckey.c_str());
   if (cp.iKnowCCPs()) {ofil << cp.nCCP;} else {ofil << -1;}
   ofil << endl;
   writeClosingAttribute(ofil,ckey.c_str());
   //
   writeClosingAttribute(ofil,mkey.c_str());
   return;
}
//**************************************************************************************************
void writeCoordinatesACPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowACPs())) {
      return;
   }
#endif
   string key="ACPCartesianCoordinates";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nACP; i++) {
      ofil << " " << setw(19) << cp.RACP[i][0];
      ofil << " " << setw(19) << cp.RACP[i][1];
      ofil << " " << setw(19) << cp.RACP[i][2] << endl;
   }
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeCoordinatesBCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowBCPs())) {
      return;
   }
#endif
   string key="BCPCartesianCoordinates";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nBCP; i++) {
      ofil << " " << setw(19) << cp.RBCP[i][0];
      ofil << " " << setw(19) << cp.RBCP[i][1];
      ofil << " " << setw(19) << cp.RBCP[i][2] << endl;
   }
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeCoordinatesRCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowRCPs())) {
      return;
   }
#endif
   string key="RCPCartesianCoordinates";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nRCP; i++) {
      ofil << " " << setw(19) << cp.RRCP[i][0];
      ofil << " " << setw(19) << cp.RRCP[i][1];
      ofil << " " << setw(19) << cp.RRCP[i][2] << endl;
   }
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeCoordinatesCCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowCCPs())) {
      return;
   }
#endif
   string key="CCPCartesianCoordinates";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nCCP; i++) {
      ofil << " " << setw(19) << cp.RCCP[i][0];
      ofil << " " << setw(19) << cp.RCCP[i][1];
      ofil << " " << setw(19) << cp.RCCP[i][2] << endl;
   }
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeConnectivityBCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowBCPs())) {
      return;
   }
#endif
   string key="BCPConnectivity";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   for (int i=0; i<cp.nBCP; i++) {
      ofil << (i+1) << " 2 " << (cp.atBCP[i][0]+1) << " " << (cp.atBCP[i][1]+1) << endl;
   }
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeLabelsACPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowACPs())) {
      return;
   }
#endif
   string key="ACPLabels";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nACP; i++) {
      ofil << " " << cp.lblACP[i];
      if ((++ninps)==5) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeLabelsBCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowBCPs())) {
      return;
   }
#endif
   string key="BCPLabels";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nBCP; i++) {
      ofil << " " << cp.lblBCP[i];
      if ((++ninps)==5) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeLabelsRCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowRCPs())) {
      return;
   }
#endif
   string key="RCPLabels";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nRCP; i++) {
      ofil << " " << cp.lblRCP[i];
      if ((++ninps)==4) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeLabelsCCPs(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowCCPs())) {
      return;
   }
#endif
   string key="CCPLabels";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   //int ninps=0;
   for (int i=0; i<cp.nCCP; i++) {
      ofil << " " << cp.lblCCP[i];
      //if ((++ninps)==5) {ofil << endl; ninps=0;}
      ofil << endl;
   }
   //if (ninps!=0) {ofil << endl;}
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeNumberOfBondPaths(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowBGPs())) {
      return;
   }
#endif
   string key="NumberOfBondPaths";
   writeOpenningAttribute(ofil,key.c_str());
   ofil << cp.nBCP << endl;
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeNumberOfPointsPerBondPath(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowBGPs())) {
      return;
   }
#endif
   string key="NumbersOfPointsPerBondPath";
   writeOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nBCP; i++) {
      ofil << " " << cp.atBCP[i][2];
      if ((++ninps)==5) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   writeClosingAttribute(ofil,key.c_str());
   return;
}
//**************************************************************************************************
void writeBondPathsCoordinates(ofstream &ofil,critPtNetWork &cp)
{
#if DEBUG
   if (!(cp.iKnowBGPs())) {
      return;
   }
#endif
   string mkey="BondPathsData";
   string ikey="BondPathIndex";
   string ckey="CoordinatesOfBondPathPoints";
   writeOpenningAttribute(ofil,mkey.c_str(),false);
   ofil << endl;
   ofil << scientific << setprecision(12);
   int nbps=cp.nBCP,npts;
   for (int i=0; i<nbps; i++) {
      writeOpenningAttribute(ofil,ikey.c_str());
      ofil << i << endl;
      writeClosingAttribute(ofil,ikey.c_str());
      writeOpenningAttribute(ofil,ckey.c_str(),false);
      ofil << endl;
      npts=cp.atBCP[i][2];
      for (int j=0; j<npts; j++) {
         for (int k=0; k<3; k++) {
            ofil << " " << setw(19) << cp.RBGP[i][j][k];
         }
         ofil << endl;
      }
      writeClosingAttribute(ofil,ckey.c_str());
   }
   writeClosingAttribute(ofil,mkey.c_str());
   return;
}
//**************************************************************************************************
bool cpxGetPosBetweenKeyInFile(ifstream &ifil,int &aftik,int &befek,string key,bool frombeg)
{
#if DEBUG
   bool iknowthiskey=false;
   for (int i=0; i<MAXCPXKEYSDEFINED; i++) {
      if (string(cpxKeysTab[i])==key) {
         iknowthiskey=true;
         break;
      }
   }
   if (!iknowthiskey) {
      displayErrorMessage("Not a CPX valid key!");
      DISPLAYDEBUGINFOFILELINE;
      return -1;
   }
#endif
   bool foundkey=false;
   int pos,orpos=ifil.tellg();
   if (frombeg) {ifil.seekg(ifil.beg);}
   string line,tk=key;
   tk.insert(0,"<"); tk.append(">");
   while (!ifil.eof()) {
      getline(ifil,line);
      pos=ifil.tellg();
      if (tk==line) {
         aftik=pos;
         foundkey=true;
         break;
      }
   }
   if (!foundkey) {return foundkey;}
   foundkey=false;
   tk.insert(1,"/");
   while (!ifil.eof()) {
      pos=ifil.tellg();
      getline(ifil,line);
      if (tk==line) {
         befek=pos-1;
         foundkey=true;
         break;
      }
   }
   if (ifil.eof()) {
      displayErrorMessage("End of file and no string found...");
   } else {
      ifil.seekg(orpos);
   }
   return foundkey;
}
//**************************************************************************************************
int cpxGetPosAfterOpenningKeyInFile(ifstream &ifil,string key,bool frombeg)
{
   int ip,fp;
   if (!(cpxGetPosBetweenKeyInFile(ifil,ip,fp,key,frombeg))) {
      displayErrorMessage("Key not found!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
      return -1;
   }
   return ip;
}
//**************************************************************************************************
int cpxGetPosBeforeClosingKeyInFile(ifstream &ifil,string key,bool frombeg)
{
   int ip,fp;
   if (!(cpxGetPosBetweenKeyInFile(ifil,ip,fp,key,frombeg))) {
      displayErrorMessage("Key not found!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
      return -1;
   }
   return fp;
}
//**************************************************************************************************
void cpxSetPosOfFileAfterOpenningKey(ifstream &ifil,string key,bool frombeg)
{
   int pos=cpxGetPosAfterOpenningKeyInFile(ifil,key,frombeg);
   ifil.seekg(pos);
   return;
}
//**************************************************************************************************
string cpxGetWFXFileName(ifstream &ifil)
{
   string wfxnam;
   int pos=cpxGetPosAfterOpenningKeyInFile(ifil,"WaveFunctionFileName",false);
   ifil.seekg(pos);
   getline(ifil,wfxnam);
   removeSpacesLeftAndRight(wfxnam);
   return wfxnam;
}
//**************************************************************************************************
ScalarFieldType cpxGetCriticalPointFieldType(ifstream &ifil)
{
   ScalarFieldType sft;
   int pos=cpxGetPosAfterOpenningKeyInFile(ifil,"CriticalPointType",false);
   ifil.seekg(pos);
   string line;
   getline(ifil,line);
   removeSpacesLeftAndRight(line);
   if (line=="Electron Density") {
      sft=DENS;
   } else if (line=="LOL") {
      sft=LOLD;
   } else {
      displayErrorMessage("Unknown Critical Points Field Type!");
      DISPLAYDEBUGINFOFILELINE;
      sft=NONE;
   }
   return sft;
}
//**************************************************************************************************
int cpxGetNOfACPs(ifstream &ifil)
{
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfACPs",true);
   ifil >> res;
   return res;
}
//**************************************************************************************************
int cpxGetNOfBCPs(ifstream &ifil)
{
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfBCPs",true);
   ifil >> res;
   return res;
}
//**************************************************************************************************
int cpxGetNOfRCPs(ifstream &ifil)
{
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfRCPs",true);
   ifil >> res;
   return res;
}
//**************************************************************************************************
int cpxGetNOfCCPs(ifstream &ifil)
{
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfCCPs",true);
   ifil >> res;
   return res;
}
//**************************************************************************************************
void cpxGetACPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr)
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"ACPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
//**************************************************************************************************
void cpxGetBCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr)
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"BCPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
//**************************************************************************************************
void cpxGetRCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr)
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"RCPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
//**************************************************************************************************
void cpxGetCCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr)
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"CCPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
//**************************************************************************************************
void cpxGetBCPConnectivityFromFile(ifstream &ifil,const int nn,int** (&ii))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"BCPConnectivity",true);
   int itmp;
   for (int i=0; i<nn; i++) {
      ifil >> itmp;
      if ((itmp-1)!=i) {displayWarningMessage("Disordered BCPs!");}
      ifil >> itmp;
      if (itmp!=2) {displayErrorMessage("Bad BCP connectivity!");}
      ifil >> ii[i][0]; ii[i][0]--;
      ifil >> ii[i][1]; ii[i][1]--;
   }
   return;
}
//**************************************************************************************************
void cpxGetACPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"ACPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
//**************************************************************************************************
void cpxGetBCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"BCPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
//**************************************************************************************************
void cpxGetRCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"RCPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
//**************************************************************************************************
void cpxGetCCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"CCPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
//**************************************************************************************************
int cpxGetNOfBondPaths(ifstream &ifil)
{
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfBondPaths",true);
   ifil >> res;
   return res;
}
//**************************************************************************************************
void cpxGetNOfPtsPerBondPath(ifstream &ifil,const int nn,int** (&ii))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumbersOfPointsPerBondPath",true);
   for (int k=0; k<nn; k++) {ifil >> ii[k][2];}
   return;
}
//**************************************************************************************************
void cpxGetBondPathData(ifstream &ifil,const int nn,int** (&ii),solreal*** (&rrr))
{
   cpxSetPosOfFileAfterOpenningKey(ifil,"BondPathsData",true);
   int ktmp;
   for (int k=0; k<nn; k++) {
      cpxSetPosOfFileAfterOpenningKey(ifil,"BondPathIndex",false);
      ifil >> ktmp;
      if (ktmp!=k) {displayWarningMessage("Disordered bond paths!");}
      cpxSetPosOfFileAfterOpenningKey(ifil,"CoordinatesOfBondPathPoints",false);
      ktmp=ii[k][2];
      for (int j=0; j<ktmp; j++) {for (int l=0; l<3; l++) {ifil >> rrr[k][j][l];}}
   }
   return;
}
//**************************************************************************************************
//**************************************************************************************************
#endif /* defined(_IOFUNCTS_CPX_CPP_) */


