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
#ifndef _IOFUNCTS_WFN_H_
#define _IOFUNCTS_WFN_H_
#include <cstdlib>
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;

#ifndef DEBUG
#define DEBUG 0
#endif

string GetTitleFromFileWFN(ifstream &ifil);
void ProcessFirstDataStringinWFNFile(ifstream &ifil,string* &tit,string & orbdesc,int &nmo,int &npr,int &nnu);
void ProcessCentersWFN(ifstream &ifil,const int nnu,string* &atlbl,double* &rr,double* &atch);
void ProcessPrimitivesWFN(ifstream &ifil,const int npr,int* &pricen,int* &primty,double* &prexp);
void ProcessMolecularOrbitalPropsAndCoefs(ifstream &ifil,const int norb,const int npr
                                 ,double* &ocn,double* &moe,double* &moc);
void GetEnergyAndVirial(ifstream &ifil,double &theener,double &thevir);
#endif//_IOFUNCTS_WFN_H_

