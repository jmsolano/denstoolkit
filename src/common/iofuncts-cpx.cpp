/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
#ifndef _IOFUNCTS_CPX_CPP_
#define _IOFUNCTS_CPX_CPP_
#include <iostream>
using std::endl;
using std::ios;
#include <iomanip>
using std::setprecision;
using std::scientific;
using std::setw;

#include "iofuncts-cpx.h"
#include "screenutils.h"
#include "stringtools.h"

#define MAXCPXKEYSDEFINED 34
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
   "NumberOfRingPaths","NumbersOfPointsPerRingPath",\
   "RingPathsData","RingPathIndex",\
   "CoordinatesOfRingPathPoints",\
   "NumberOfCagePaths","NumbersOfPointsPerCagePath",\
   "CagePathsData","CagePathIndex",\
   "CoordinatesOfCagePathPoints","ACPConnectivity",\
   "BCPConnectivity","RCPConnectivity",\
   "CCPConnectivity"
};
void WriteCPXFile(string cpxname,string wfname,CritPtNetWork &cp) {
   ofstream cpxfil;
   cpxfil.open(cpxname.c_str(),ios::out);
   WriteWFFileName(cpxfil,wfname);
   WriteTypeOfCriticalPoints(cpxfil,cp);
   WriteNumberOfCriticalPoints(cpxfil,cp);
   WriteCoordinatesACPs(cpxfil,cp);
   WriteCoordinatesBCPs(cpxfil,cp);
   WriteCoordinatesRCPs(cpxfil,cp);
   WriteCoordinatesCCPs(cpxfil,cp);
   WriteConnectivityBCPs(cpxfil,cp);
   if (cp.IKnowRGPs()) {WriteConnectivityRCPs(cpxfil,cp);}
   if (cp.IKnowCGPs()) {WriteConnectivityCCPs(cpxfil,cp);}
   WriteLabelsACPs(cpxfil,cp);
   WriteLabelsBCPs(cpxfil,cp);
   WriteLabelsRCPs(cpxfil,cp);
   WriteLabelsCCPs(cpxfil,cp);
   WriteNumberOfBondPaths(cpxfil,cp);
   WriteNumberOfPointsPerBondPath(cpxfil,cp);
   WriteBondPathsCoordinates(cpxfil,cp);
   if ( cp.IKnowRGPs() ) {
      WriteNumberOfRingPaths(cpxfil,cp);
      WriteNumberOfPointsPerRingPath(cpxfil,cp);
      WriteRingPathsCoordinates(cpxfil,cp);
   }
   if ( cp.IKnowCGPs() ) {
      WriteNumberOfCagePaths(cpxfil,cp);
      WriteNumberOfPointsPerCagePath(cpxfil,cp);
      WriteCagePathsCoordinates(cpxfil,cp);
   }
   cpxfil.close();
   return;
}
void WriteOpenningAttribute(ofstream &ofil,const char* attrib,bool nl = true) {
   ofil << "<" << attrib << ">";
   if (nl) {ofil << endl << " ";}
   return;
}
void WriteClosingAttribute(ofstream &ofil,const char* attrib,bool nl = true) {
   ofil << "</" << attrib << ">";
   if (nl) {ofil << endl;}
   return;
}
void WriteWFFileName(ofstream &ofil,string wfname) {
   string key="WaveFunctionFileName";
   WriteOpenningAttribute(ofil,key.c_str());
   ofil << wfname << endl;
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteTypeOfCriticalPoints(ofstream &ofil,CritPtNetWork &cp) {
   string key="CriticalPointType";
   WriteOpenningAttribute(ofil,key.c_str());
   switch (cp.MyCPType()) {
      case DENS:
         ofil << "Electron Density";
         break;
      case LOLD:
         ofil << "LOL";
         break;
      default:
         ofil << "Unknown";
#if DEBUG
         ScreenUtils::DisplayErrorMessage("Unknown Field Type!");
#endif
         break;
   }
   ofil << endl;
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfCriticalPoints(ofstream &ofil,CritPtNetWork &cp) {
   string mkey="NumberOfCriticalPoints";
   string akey="NumberOfACPs";
   string bkey="NumberOfBCPs";
   string rkey="NumberOfRCPs";
   string ckey="NumberOfCCPs";
   //
   WriteOpenningAttribute(ofil,mkey.c_str(),false);
   ofil << endl;
   //
   WriteOpenningAttribute(ofil,akey.c_str());
   if (cp.IKnowACPs()) {ofil << cp.nACP;} else {ofil << -1;}
   ofil << endl;
   WriteClosingAttribute(ofil,akey.c_str());
   //
   WriteOpenningAttribute(ofil,bkey.c_str());
   if (cp.IKnowBCPs()) {ofil << cp.nBCP;} else {ofil << -1;}
   ofil << endl;
   WriteClosingAttribute(ofil,bkey.c_str());
   //
   WriteOpenningAttribute(ofil,rkey.c_str());
   if (cp.IKnowRCPs()) {ofil << cp.nRCP;} else {ofil << -1;}
   ofil << endl;
   WriteClosingAttribute(ofil,rkey.c_str());
   //
   WriteOpenningAttribute(ofil,ckey.c_str());
   if (cp.IKnowCCPs()) {ofil << cp.nCCP;} else {ofil << -1;}
   ofil << endl;
   WriteClosingAttribute(ofil,ckey.c_str());
   //
   WriteClosingAttribute(ofil,mkey.c_str());
   return;
}
void WriteCoordinatesACPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowACPs())) {
      return;
   }
#endif
   string key="ACPCartesianCoordinates";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nACP; i++) {
      ofil << " " << setw(19) << cp.RACP[i][0];
      ofil << " " << setw(19) << cp.RACP[i][1];
      ofil << " " << setw(19) << cp.RACP[i][2] << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteCoordinatesBCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBCPs())) {
      return;
   }
#endif
   string key="BCPCartesianCoordinates";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nBCP; i++) {
      ofil << " " << setw(19) << cp.RBCP[i][0];
      ofil << " " << setw(19) << cp.RBCP[i][1];
      ofil << " " << setw(19) << cp.RBCP[i][2] << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteCoordinatesRCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowRCPs())) {
      return;
   }
#endif
   string key="RCPCartesianCoordinates";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nRCP; i++) {
      ofil << " " << setw(19) << cp.RRCP[i][0];
      ofil << " " << setw(19) << cp.RRCP[i][1];
      ofil << " " << setw(19) << cp.RRCP[i][2] << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteCoordinatesCCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowCCPs())) {
      return;
   }
#endif
   string key="CCPCartesianCoordinates";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl << scientific << setprecision(12);
   for (int i=0; i<cp.nCCP; i++) {
      ofil << " " << setw(19) << cp.RCCP[i][0];
      ofil << " " << setw(19) << cp.RCCP[i][1];
      ofil << " " << setw(19) << cp.RCCP[i][2] << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteConnectivityBCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBCPs())) {
      return;
   }
#endif
   string key="BCPConnectivity";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   for (int i=0; i<cp.nBCP; i++) {
      ofil << (i+1) << " 2 " << (cp.conBCP[i][0]+1) << " " << (cp.conBCP[i][1]+1) << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteConnectivityRCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowRCPs())) {
      return;
   }
#endif
   ofil << "#RCPConnectivity format:" << endl
        << "#<rcpIdx> <n = Num Of BCPs associated with RCP rcpIdx>"
        << " <bcpIdx 1> ... <bcpIdx n>" << endl;
   string key="RCPConnectivity";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int nofrgps;
   for (int i=0; i<cp.nRCP; ++i) {
      nofrgps=cp.GetNofRingPathsOfRCP(i);
      ofil << (i+1) << " " << nofrgps;
      for ( int k=0 ; k<nofrgps ; ++k ) {ofil << " " << cp.conRCP[i][0][k];}
      ofil << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteConnectivityCCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowCCPs())) {
      return;
   }
#endif
   ofil << "#CCPConnectivity format:" << endl
        << "#<ccpIdx> <n = Num Of RCPs associated with CCP ccpIdx>"
        << " <rcpIdx 1> ... <rcpIdx n>" << endl;
   string key="CCPConnectivity";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int nofcgps;
   for (int i=0; i<cp.nCCP; ++i) {
      nofcgps=cp.GetNofCagePathsOfCCP(i);
      ofil << (i+1) << " " << nofcgps;
      for ( int k=0 ; k<nofcgps ; ++k ) {ofil << " " << cp.conCCP[i][0][k];}
      ofil << endl;
   }
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteLabelsACPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowACPs())) {
      return;
   }
#endif
   string key="ACPLabels";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nACP; i++) {
      ofil << " " << cp.lblACP[i];
      if ((++ninps)==5) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteLabelsBCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBCPs())) {
      return;
   }
#endif
   string key="BCPLabels";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nBCP; i++) {
      ofil << " " << cp.lblBCP[i];
      if ((++ninps)==5) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteLabelsRCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowRCPs())) {
      return;
   }
#endif
   string key="RCPLabels";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nRCP; i++) {
      ofil << " " << cp.lblRCP[i];
      if ((++ninps)==2) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteLabelsCCPs(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowCCPs())) {
      return;
   }
#endif
   string key="CCPLabels";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   //int ninps=0;
   for (int i=0; i<cp.nCCP; i++) {
      ofil << " " << cp.lblCCP[i];
      //if ((++ninps)==5) {ofil << endl; ninps=0;}
      ofil << endl;
   }
   //if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfBondPaths(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBGPs())) {
      return;
   }
#endif
   string key="NumberOfBondPaths";
   WriteOpenningAttribute(ofil,key.c_str());
   ofil << cp.nBCP << endl;
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfRingPaths(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowRGPs())) {
      return;
   }
#endif
   string key="NumberOfRingPaths";
   ofil << "#Reminder: NumberOfRingPaths is the total number." << endl;
   ofil << "#The connectivity is in RCPConnectivity." << endl;
   WriteOpenningAttribute(ofil,key.c_str());
   ofil << cp.GetTotalNofRingPaths() << endl;
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfCagePaths(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowCGPs())) {
      return;
   }
#endif
   string key="NumberOfCagePaths";
   ofil << "#Reminder: NumberOfCagePaths is the total number." << endl;
   ofil << "#The connectivity is in CCPConnectivity." << endl;
   WriteOpenningAttribute(ofil,key.c_str());
   ofil << cp.GetTotalNofCagePaths() << endl;
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfPointsPerBondPath(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBGPs())) {
      return;
   }
#endif
   string key="NumbersOfPointsPerBondPath";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int ninps=0;
   for (int i=0; i<cp.nBCP; i++) {
      ofil << " " << cp.conBCP[i][2];
      if ((++ninps)==5) {ofil << endl; ninps=0;}
   }
   if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfPointsPerRingPath(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBGPs())) {
      return;
   }
#endif
   string key="NumbersOfPointsPerRingPath";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int nofrgps,ninps=0;
   for (int i=0; i<cp.nRCP; ++i) {
      nofrgps=cp.GetNofRingPathsOfRCP(i);
      for ( int k=0 ; k<nofrgps ; ++k ) {
         ofil << " " << cp.conRCP[i][1][k];
         if ((++ninps)==5) {ofil << endl; ninps=0;}
      }
   }
   if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteNumberOfPointsPerCagePath(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowCGPs())) {
      return;
   }
#endif
   string key="NumbersOfPointsPerCagePath";
   WriteOpenningAttribute(ofil,key.c_str(),false);
   ofil << endl;
   int nofcgps,ninps=0;
   for (int i=0; i<cp.nCCP; ++i) {
      nofcgps=cp.GetNofCagePathsOfCCP(i);
      for ( int k=0 ; k<nofcgps ; ++k ) {
         ofil << " " << cp.conCCP[i][1][k];
         if ((++ninps)==5) {ofil << endl; ninps=0;}
      }
   }
   if (ninps!=0) {ofil << endl;}
   WriteClosingAttribute(ofil,key.c_str());
   return;
}
void WriteBondPathsCoordinates(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowBGPs())) {
      return;
   }
#endif
   string mkey="BondPathsData";
   string ikey="BondPathIndex";
   string ckey="CoordinatesOfBondPathPoints";
   WriteOpenningAttribute(ofil,mkey.c_str(),false);
   ofil << endl;
   ofil << scientific << setprecision(12);
   int nbps=cp.nBCP,npts;
   for (int i=0; i<nbps; i++) {
      WriteOpenningAttribute(ofil,ikey.c_str());
      ofil << i << endl;
      WriteClosingAttribute(ofil,ikey.c_str());
      WriteOpenningAttribute(ofil,ckey.c_str(),false);
      ofil << endl;
      npts=cp.conBCP[i][2];
      for (int j=0; j<npts; j++) {
         for (int k=0; k<3; k++) {
            ofil << " " << setw(19) << cp.RBGP[i][j][k];
         }
         ofil << endl;
      }
      WriteClosingAttribute(ofil,ckey.c_str());
   }
   WriteClosingAttribute(ofil,mkey.c_str());
   return;
}
void WriteRingPathsCoordinates(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowRGPs())) {
      return;
   }
#endif
   string mkey="RingPathsData";
   string ikey="RingPathIndex";
   string ckey="CoordinatesOfRingPathPoints";
   WriteOpenningAttribute(ofil,mkey.c_str(),false);
   ofil << endl;
   ofil << scientific << setprecision(12);
   int npts,nrps,count=0;
   for (int rcpIdx=0; rcpIdx<cp.nRCP; ++rcpIdx) {
      nrps=cp.GetNofRingPathsOfRCP(rcpIdx);
      for ( int bcpIdxInRRGP=0 ; bcpIdxInRRGP<nrps ; ++bcpIdxInRRGP ) {
         WriteOpenningAttribute(ofil,ikey.c_str());
         ofil << count << endl;
         WriteClosingAttribute(ofil,ikey.c_str());
         WriteOpenningAttribute(ofil,ckey.c_str(),false);
         ofil << endl;
         npts=cp.conRCP[rcpIdx][1][bcpIdxInRRGP];
         for (int j=0; j<npts; j++) {
            for (int k=0; k<3; k++) {
               ofil << " " << setw(19) << cp.RRGP[rcpIdx][bcpIdxInRRGP][j][k];
            }
            ofil << endl;
         }
         WriteClosingAttribute(ofil,ckey.c_str());
         ++count;
      }
   }
   WriteClosingAttribute(ofil,mkey.c_str());
   return;
}
void WriteCagePathsCoordinates(ofstream &ofil,CritPtNetWork &cp) {
#if DEBUG
   if (!(cp.IKnowCGPs())) {
      return;
   }
#endif
   string mkey="CagePathsData";
   string ikey="CagePathIndex";
   string ckey="CoordinatesOfCagePathPoints";
   WriteOpenningAttribute(ofil,mkey.c_str(),false);
   ofil << endl;
   ofil << scientific << setprecision(12);
   int npts,ncps,count=0;
   for (int ccpIdx=0; ccpIdx<cp.nCCP; ++ccpIdx) {
      ncps=cp.GetNofCagePathsOfCCP(ccpIdx);
      for ( int rcpIdxInRCGP=0 ; rcpIdxInRCGP<ncps ; ++rcpIdxInRCGP ) {
         WriteOpenningAttribute(ofil,ikey.c_str());
         ofil << count << endl;
         WriteClosingAttribute(ofil,ikey.c_str());
         WriteOpenningAttribute(ofil,ckey.c_str(),false);
         ofil << endl;
         npts=cp.conCCP[ccpIdx][1][rcpIdxInRCGP];
         for (int j=0; j<npts; j++) {
            for (int k=0; k<3; k++) {
               ofil << " " << setw(19) << cp.RCGP[ccpIdx][rcpIdxInRCGP][j][k];
            }
            ofil << endl;
         }
         WriteClosingAttribute(ofil,ckey.c_str());
         ++count;
      }
   }
   WriteClosingAttribute(ofil,mkey.c_str());
   return;
}
bool cpxGetPosBetweenKeyInFile(ifstream &ifil,int &aftik,int &befek,string key,bool frombeg) {
#if DEBUG
   bool iknowthiskey=false;
   for (int i=0; i<MAXCPXKEYSDEFINED; i++) {
      if (string(cpxKeysTab[i])==key) {
         iknowthiskey=true;
         break;
      }
   }
   if (!iknowthiskey) {
      ScreenUtils::DisplayErrorMessage("Not a CPX valid key!");
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
      ScreenUtils::DisplayErrorMessage("End of file and no string found...");
   } else {
      ifil.seekg(orpos);
   }
   return foundkey;
}
int cpxGetPosAfterOpenningKeyInFile(ifstream &ifil,string key,bool frombeg) {
   int ip,fp;
   if (!(cpxGetPosBetweenKeyInFile(ifil,ip,fp,key,frombeg))) {
      ScreenUtils::DisplayErrorMessage("Key not found!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
      return -1;
   }
   return ip;
}
int cpxGetPosBeforeClosingKeyInFile(ifstream &ifil,string key,bool frombeg) {
   int ip,fp;
   if (!(cpxGetPosBetweenKeyInFile(ifil,ip,fp,key,frombeg))) {
      ScreenUtils::DisplayErrorMessage("Key not found!");
#if DEBUG
      DISPLAYDEBUGINFOFILELINE;
#endif
      return -1;
   }
   return fp;
}
void cpxSetPosOfFileAfterOpenningKey(ifstream &ifil,string key,bool frombeg) {
   int pos=cpxGetPosAfterOpenningKeyInFile(ifil,key,frombeg);
   ifil.seekg(pos);
   return;
}
string cpxGetWFXFileName(const string cpxname) {
   ifstream ifil(cpxname.c_str());
   string res=cpxGetWFXFileName(ifil);
   ifil.close();
   ifil.clear();
   return res;
}
string cpxGetWFXFileName(ifstream &ifil) {
   string wfxnam;
   int pos=cpxGetPosAfterOpenningKeyInFile(ifil,"WaveFunctionFileName",false);
   ifil.seekg(pos);
   getline(ifil,wfxnam);
   StringTools::RemoveSpacesLeftAndRight(wfxnam);
   return wfxnam;
}
ScalarFieldType cpxGetCriticalPointFieldType(ifstream &ifil) {
   ScalarFieldType sft;
   int pos=cpxGetPosAfterOpenningKeyInFile(ifil,"CriticalPointType",false);
   ifil.seekg(pos);
   string line;
   getline(ifil,line);
   StringTools::RemoveSpacesLeftAndRight(line);
   if (line=="Electron Density") {
      sft=DENS;
   } else if (line=="LOL") {
      sft=LOLD;
   } else {
      ScreenUtils::DisplayErrorMessage("Unknown Critical Points Field Type!");
      DISPLAYDEBUGINFOFILELINE;
      sft=NONE;
   }
   return sft;
}
int cpxGetNOfACPs(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfACPs",true);
   ifil >> res;
   return res;
}
int cpxGetNOfBCPs(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfBCPs",true);
   ifil >> res;
   return res;
}
int cpxGetNOfRCPs(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfRCPs",true);
   ifil >> res;
   return res;
}
int cpxGetNOfCCPs(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfCCPs",true);
   ifil >> res;
   return res;
}
void cpxGetACPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"ACPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
void cpxGetBCPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"BCPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
void cpxGetRCPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"RCPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
void cpxGetCCPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"CCPCartesianCoordinates",true);
   for (int i=0; i<nn; i++) {for (int j=0; j<3; j++) {ifil >> rr[i][j];}}
   return;
}
void cpxGetBCPConnectivityFromFile(ifstream &ifil,const int nn,int** (&ii)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"BCPConnectivity",true);
   int itmp;
   for (int i=0; i<nn; i++) {
      ifil >> itmp;
      if ((itmp-1)!=i) {ScreenUtils::DisplayWarningMessage("Disordered BCPs!");}
      ifil >> itmp;
      if (itmp!=2) {ScreenUtils::DisplayErrorMessage("Bad BCP connectivity!");}
      ifil >> ii[i][0]; ii[i][0]--;
      ifil >> ii[i][1]; ii[i][1]--;
   }
   return;
}
void cpxGetACPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"ACPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
void cpxGetBCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"BCPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
void cpxGetRCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"RCPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
void cpxGetCCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"CCPLabels",true);
   for (int i=0; i<nn; i++) {ifil >> ss[i];}
   return;
}
int cpxGetNOfBondPaths(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfBondPaths",true);
   ifil >> res;
   return res;
}
void cpxGetNOfPtsPerBondPath(ifstream &ifil,const int nn,int** (&ii)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumbersOfPointsPerBondPath",true);
   for (int k=0; k<nn; k++) {ifil >> ii[k][2];}
   return;
}
void cpxGetBondPathData(ifstream &ifil,const int nn,int** (&ii),double*** (&rrr)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"BondPathsData",true);
   int ktmp;
   for (int k=0; k<nn; k++) {
      cpxSetPosOfFileAfterOpenningKey(ifil,"BondPathIndex",false);
      ifil >> ktmp;
      if (ktmp!=k) {ScreenUtils::DisplayWarningMessage("Disordered bond paths!");}
      cpxSetPosOfFileAfterOpenningKey(ifil,"CoordinatesOfBondPathPoints",false);
      ktmp=ii[k][2];
      for (int j=0; j<ktmp; j++) {for (int l=0; l<3; l++) {ifil >> rrr[k][j][l];}}
   }
   return;
}
void cpxGetRCPConnectivityFromFile(ifstream &ifil,const int nn,int*** (&cc)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"RCPConnectivity",true);
   int nrgps;
   for (int i=0; i<nn; i++) {
      ifil >> nrgps;
      if ((nrgps-1)!=i) {ScreenUtils::DisplayWarningMessage("Disordered RCPs!");}
      ifil >> nrgps;
      for ( int j=0 ; j<nrgps ; ++j ) {
         ifil >> cc[i][0][j];
      }
   }
}
int cpxGetNOfRingPaths(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfRingPaths",true);
   ifil >> res;
   return res;
}
void cpxGetNOfPtsPerRingPath(ifstream &ifil,const int nn,int*** (&ii)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumbersOfPointsPerRingPath",true);
   int k;
   for ( int i=0 ; i<nn ; ++i ) {
      k=0;
      while ( ii[i][0][k]>=0 ) {
         ifil >> ii[i][1][k++];
      }
   }
   return;
}
void cpxGetRingPathData(ifstream &ifil,const int nn,int*** (&ii),double**** (&rrr)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"RingPathsData",true);
   int nnpts,ktmp,glblIdx=0;
   for (int rcpIdx=0; rcpIdx<nn; ++rcpIdx) {
      cpxSetPosOfFileAfterOpenningKey(ifil,"RingPathIndex",false);
      ifil >> nnpts;
      if (nnpts!=glblIdx) {ScreenUtils::DisplayWarningMessage("Disordered bond paths!");}
      ktmp=0;
      while ( ii[rcpIdx][0][ktmp]>=0 ) {
         nnpts=ii[rcpIdx][1][ktmp];
         cpxSetPosOfFileAfterOpenningKey(ifil,"CoordinatesOfRingPathPoints",false);
         for ( int l=0 ; l<nnpts ; ++l ) {
            ifil >> rrr[rcpIdx][ktmp][l][0];
            ifil >> rrr[rcpIdx][ktmp][l][1];
            ifil >> rrr[rcpIdx][ktmp][l][2];
         }
         ++ktmp;
         ++glblIdx;
      }
   }
   return;
}
void cpxGetCCPConnectivityFromFile(ifstream &ifil,const int nn,int*** (&cc)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"CCPConnectivity",true);
   int ncgps;
   for (int i=0; i<nn; i++) {
      ifil >> ncgps;
      if ((ncgps-1)!=i) {ScreenUtils::DisplayWarningMessage("Disordered CCPs!");}
      ifil >> ncgps;
      for ( int j=0 ; j<ncgps ; ++j ) {
         ifil >> cc[i][0][j];
      }
   }
}
int cpxGetNOfCagePaths(ifstream &ifil) {
   int res;
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumberOfCagePaths",true);
   ifil >> res;
   return res;
}
void cpxGetNOfPtsPerCagePath(ifstream &ifil,const int nn,int*** (&ii)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"NumbersOfPointsPerCagePath",true);
   int k;
   for ( int i=0 ; i<nn ; ++i ) {
      k=0;
      while ( ii[i][0][k]>=0 ) {
         ifil >> ii[i][1][k++];
      }
   }
   return;
}
void cpxGetCagePathData(ifstream &ifil,const int nn,int*** (&ii),double**** (&rrr)) {
   cpxSetPosOfFileAfterOpenningKey(ifil,"CagePathsData",true);
   int nnpts,ktmp,glblIdx=0;
   for (int ccpIdx=0; ccpIdx<nn; ++ccpIdx) {
      cpxSetPosOfFileAfterOpenningKey(ifil,"CagePathIndex",false);
      ifil >> nnpts;
      if (nnpts!=glblIdx) {ScreenUtils::DisplayWarningMessage("Disordered bond paths!");}
      ktmp=0;
      while ( ii[ccpIdx][0][ktmp]>=0 ) {
         nnpts=ii[ccpIdx][1][ktmp];
         cpxSetPosOfFileAfterOpenningKey(ifil,"CoordinatesOfCagePathPoints",false);
         for ( int l=0 ; l<nnpts ; ++l ) {
            ifil >> rrr[ccpIdx][ktmp][l][0];
            ifil >> rrr[ccpIdx][ktmp][l][1];
            ifil >> rrr[ccpIdx][ktmp][l][2];
         }
         ++ktmp;
         ++glblIdx;
      }
   }
   return;
}
#endif /* defined(_IOFUNCTS_CPX_CPP_) */


