/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.1
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
void writeCPXFile(string cpxname,string wfname,critPtNetWork &cp);
/* ************************************************************************** */
void writeOpenningAttribute(ofstream &ofil,const char* attrib,bool nl);
/* ************************************************************************** */
void writeClosingAttribute(ofstream &ofil,const char* attrib,bool nl);
/* ************************************************************************** */
void writeWFFileName(ofstream &ofil,string wfname);
/* ************************************************************************** */
void writeTypeOfCriticalPoints(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfCriticalPoints(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeCoordinatesACPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeCoordinatesBCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeCoordinatesRCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeCoordinatesCCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeConnectivityBCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeConnectivityRCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeConnectivityCCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeLabelsACPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeLabelsBCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************* */
void writeLabelsRCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeLabelsCCPs(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfBondPaths(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfRingPaths(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfCagePaths(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfPointsPerBondPath(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfPointsPerRingPath(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeNumberOfPointsPerCagePath(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeBondPathsCoordinates(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeRingPathsCoordinates(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
void writeCagePathsCoordinates(ofstream &ofil,critPtNetWork &cp);
/* ************************************************************************** */
bool cpxGetPosBetweenKeyInFile(ifstream &ifil,int &aftik,int &befek,string key,bool frombeg);
/* ************************************************************************** */
int cpxGetPosAfterOpenningKeyInFile(ifstream &ifil,string key,bool frombeg);
/* ************************************************************************** */
void cpxSetPosOfFileAfterOpenningKey(ifstream &ifil,string key,bool frombeg);
/* ************************************************************************** */
int cpxGetPosBeforeClosingKeyInFile(ifstream &ifil,string key,bool frombeg);
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
void cpxGetACPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
/* ************************************************************************** */
void cpxGetBCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
/* ************************************************************************** */
void cpxGetRCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
/* ************************************************************************** */
void cpxGetCCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
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
void cpxGetBondPathData(ifstream &ifil,const int nn,int** (&ii),solreal*** (&rrr));
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
#endif /* defined(_IOFUNCTS_CPX_H_) */
