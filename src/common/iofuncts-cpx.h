/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
#ifndef _IOFUNCTS_CPX_H_
#define _IOFUNCTS_CPX_H_

#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;

#include "fldtypesdef.h"
#include "critptnetwork.h"

/* ************************************************************************** */
void WriteCPXFile(string cpxname,string wfname,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteOpenningAttribute(ofstream &ofil,const char* attrib,bool nl);
/* ************************************************************************** */
void WriteClosingAttribute(ofstream &ofil,const char* attrib,bool nl);
/* ************************************************************************** */
void WriteWFFileName(ofstream &ofil,string wfname);
/* ************************************************************************** */
void WriteTypeOfCriticalPoints(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfCriticalPoints(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteCoordinatesACPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteCoordinatesBCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteCoordinatesRCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteCoordinatesCCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteConnectivityBCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteConnectivityRCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteConnectivityCCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteLabelsACPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteLabelsBCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************* */
void WriteLabelsRCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteLabelsCCPs(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfBondPaths(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfRingPaths(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfCagePaths(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfPointsPerBondPath(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfPointsPerRingPath(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteNumberOfPointsPerCagePath(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteBondPathsCoordinates(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteRingPathsCoordinates(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
void WriteCagePathsCoordinates(ofstream &ofil,CritPtNetWork &cp);
/* ************************************************************************** */
bool cpxGetPosBetweenKeyInFile(ifstream &ifil,int &aftik,int &befek,string key,bool frombeg);
/* ************************************************************************** */
int cpxGetPosAfterOpenningKeyInFile(ifstream &ifil,string key,bool frombeg);
/* ************************************************************************** */
void cpxSetPosOfFileAfterOpenningKey(ifstream &ifil,string key,bool frombeg);
/* ************************************************************************** */
int cpxGetPosBeforeClosingKeyInFile(ifstream &ifil,string key,bool frombeg);
/* ************************************************************************** */
string cpxGetWFXFileName(const string cpxname);
/* ************************************************************************** */
string cpxGetWFXFileName(ifstream &ifil);
/* ************************************************************************** */
ScalarFieldType cpxGetCriticalPointFieldType(ifstream &ifil);
/* ************************************************************************** */
int cpxGetNOfACPs(ifstream &ifil);
/* ************************************************************************** */
int cpxGetNOfBCPs(ifstream &ifil);
/* ************************************************************************** */
int cpxGetNOfRCPs(ifstream &ifil);
/* ************************************************************************** */
int cpxGetNOfCCPs(ifstream &ifil);
/* ************************************************************************** */
void cpxGetACPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr);
/* ************************************************************************** */
void cpxGetBCPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr);
/* ************************************************************************** */
void cpxGetRCPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr);
/* ************************************************************************** */
void cpxGetCCPCartCoordFromFile(ifstream &ifil,const int nn,double** &rr);
/* ************************************************************************** */
void cpxGetBCPConnectivityFromFile(ifstream &ifil,const int nn,int** (&ii));
/* ************************************************************************** */
void cpxGetACPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
/* ************************************************************************** */
void cpxGetBCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
/* ************************************************************************** */
void cpxGetRCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
/* ************************************************************************** */
void cpxGetCCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
/* ************************************************************************** */
int cpxGetNOfBondPaths(ifstream &ifil);
/* ************************************************************************** */
void cpxGetNOfPtsPerBondPath(ifstream &ifil,const int nn,int** (&ii));
/* ************************************************************************** */
void cpxGetBondPathData(ifstream &ifil,const int nn,int** (&ii),double*** (&rrr));
/* ************************************************************************** */
void cpxGetRCPConnectivityFromFile(ifstream &ifil,const int nn,int*** (&cc));
/* ************************************************************************** */
int cpxGetNOfRingPaths(ifstream &ifil);
/* ************************************************************************** */
void cpxGetNOfPtsPerRingPath(ifstream &ifil,const int nn,int*** (&ii));
/* ************************************************************************** */
void cpxGetRingPathData(ifstream &ifil,const int nn,int*** (&ii),double**** (&rrr));
/* ************************************************************************** */
void cpxGetCCPConnectivityFromFile(ifstream &ifil,const int nn,int*** (&cc));
/* ************************************************************************** */
int cpxGetNOfCagePaths(ifstream &ifil);
/* ************************************************************************** */
void cpxGetNOfPtsPerCagePath(ifstream &ifil,const int nn,int*** (&ii));
/* ************************************************************************** */
void cpxGetCagePathData(ifstream &ifil,const int nn,int*** (&ii),double**** (&rrr));
/* ************************************************************************** */
#endif /* defined(_IOFUNCTS_CPX_H_) */
