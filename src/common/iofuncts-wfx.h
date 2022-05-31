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
#ifndef _IOFUNCTS_WFX_H_
#define _IOFUNCTS_WFX_H_
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;

#ifndef DEBUG
#define DEBUG 0
#endif

string GetWfxKey(const int ii);
void GetPosInFile(ifstream &ifil,bool frombeg,string idkey,int & ist,int &ien);
int GetInitPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey);
int GetFinPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey);
bool GetTitleFromFileWFX(ifstream &ifil,int &nt,string* &tit);
bool GetKeyWordsFromFileWFX(ifstream &ifil,string &kw);
void GetNetCharge(ifstream &ifil,int &nc);
void GetNofNucleiFromFileWFX(ifstream &ifil, int &nnuc);
void GetNofPrimFromFileWFX(ifstream &ifil,int &npri);
void GetNofEDFPrimFromFileWFX(ifstream &ifil,const int nedfc,int &nedfprim);
void GetNofMolOrbFromFileWFX(ifstream &ifil,int &nmo);
void GetNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,double* &rr);
void GetNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,double** &rr3);
void GetAtLabelsFromFileWFX(ifstream &ifil,const int nn,string* &al);
void GetAtChargesFromFileWFX(ifstream &ifil,const int nn,double* &ach);
int GetAtNumbersFromFileWFX(ifstream &ifil,const int nn,int* &anu);
void GetNofElectronsFromFileWFX(ifstream &ifil,int &nel);
void GetNofCoreElectronsFromFileWFX(ifstream &ifil,int &ncel);
void GetPrimCentersFromFileWFX(ifstream &ifil,const int npr,int* &pc);
void GetEDFPrimCentersFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pc);
void GetPrimTypesFromFileWFX(ifstream &ifil,const int npr,int* &pt);
void GetEDFPrimTypesFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pt);
void GetPrimExponentsFromFileWFX(ifstream &ifil,const int npr,double* &pex);
void GetEDFPrimExponentsFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,double* &pex);
void GetMolecOrbOccNumsFromFileWFX(ifstream &ifil,const int nmo,double* &ocnu);
void GetMolecOrbEnergiesFromFileWFX(ifstream &ifil,const int nmo,double* &orben);
void GetMolecOrbCoefficientsFromFileWFX(ifstream &ifil,const int nmo,const int npr,double* &tcf);
void GetTotEnerAndVirialFromFileWFX(ifstream &ifil,double &tote,double &vir);
void countEDFCentersFromFileWFX(ifstream &ifil,int &nedfc);
void GetEDFExistenceFromFileWFX(ifstream &ifil,bool &ihaveEDF);
void GetEDFPrimCoefficientsFromFileWFX(ifstream &ifil,\
      const int nedfp,double* &edfc);

#endif//_IOFUNCTS_WFX_H_

