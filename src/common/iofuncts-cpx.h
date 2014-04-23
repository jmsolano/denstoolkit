//
//  iofuncts-cpx.h
//  
//
//  Created by Juan Manuel Solano on 2013-10-28.
//
//

#ifndef _IOFUNCTS_CPX_H_
#define _IOFUNCTS_CPX_H_

#include <iostream>
using std::endl;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;

#include "fldtypesdef.h"
#include "critptnetwork.h"

//**************************************************************************************************
void writeCPXFile(string cpxname,string wfname,critPtNetWork &cp);
//**************************************************************************************************
void writeOpenningAttribute(ofstream &ofil,const char* attrib,bool nl);
//**************************************************************************************************
void writeClosingAttribute(ofstream &ofil,const char* attrib,bool nl);
//**************************************************************************************************
void writeWFFileName(ofstream &ofil,string wfname);
//**************************************************************************************************
void writeTypeOfCriticalPoints(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeNumberOfCriticalPoints(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeCoordinatesACPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeCoordinatesBCPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeCoordinatesRCPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeCoordinatesCCPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeConnectivityBCPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeLabelsACPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeLabelsBCPs(ofstream &ofil,critPtNetWork &cp);
//*************************************************************************************************
void writeLabelsRCPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeLabelsCCPs(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeNumberOfBondPaths(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeNumberOfPointsPerBondPath(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
void writeBondPathsCoordinates(ofstream &ofil,critPtNetWork &cp);
//**************************************************************************************************
bool cpxGetPosBetweenKeyInFile(ifstream &ifil,int &aftik,int &befek,string key,bool frombeg);
//**************************************************************************************************
int cpxGetPosAfterOpenningKeyInFile(ifstream &ifil,string key,bool frombeg);
//**************************************************************************************************
void cpxSetPosOfFileAfterOpenningKey(ifstream &ifil,string key,bool frombeg);
//**************************************************************************************************
int cpxGetPosBeforeClosingKeyInFile(ifstream &ifil,string key,bool frombeg);
//**************************************************************************************************
string cpxGetWFXFileName(ifstream &ifil);
//**************************************************************************************************
ScalarFieldType cpxGetCriticalPointFieldType(ifstream &ifil);
//**************************************************************************************************
int cpxGetNOfACPs(ifstream &ifil);
//**************************************************************************************************
int cpxGetNOfBCPs(ifstream &ifil);
//**************************************************************************************************
int cpxGetNOfRCPs(ifstream &ifil);
//**************************************************************************************************
int cpxGetNOfCCPs(ifstream &ifil);
//**************************************************************************************************
void cpxGetACPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
//**************************************************************************************************
void cpxGetBCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
//**************************************************************************************************
void cpxGetRCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
//**************************************************************************************************
void cpxGetCCPCartCoordFromFile(ifstream &ifil,const int nn,solreal** &rr);
//**************************************************************************************************
void cpxGetBCPConnectivityFromFile(ifstream &ifil,const int nn,int** (&ii));
//**************************************************************************************************
void cpxGetACPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
//**************************************************************************************************
void cpxGetBCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
//**************************************************************************************************
void cpxGetRCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
//**************************************************************************************************
void cpxGetCCPLabelsFromFile(ifstream &ifil,const int nn,string* (&ss));
//**************************************************************************************************
int cpxGetNOfBondPaths(ifstream &ifil);
//**************************************************************************************************
void cpxGetNOfPtsPerBondPath(ifstream &ifil,const int nn,int** (&ii));
//**************************************************************************************************
void cpxGetBondPathData(ifstream &ifil,const int nn,int** (&ii),solreal*** (&rrr));
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
#endif /* defined(_IOFUNCTS_CPX_H_) */
