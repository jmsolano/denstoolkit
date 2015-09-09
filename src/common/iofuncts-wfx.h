/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
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



#ifndef _IOFUNCTS_WFX_H_
#define _IOFUNCTS_WFX_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
#endif

#include <string>
using namespace std;
#include <fstream>
using std::ifstream;
using std::ofstream;

#ifndef DEBUG
#define DEBUG 0
#endif

string getWfxKey(const int ii);
void getPosInFile(ifstream &ifil,bool frombeg,string idkey,int & ist,int &ien);
int getInitPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey);
int getFinPosOfKeyInFile(ifstream &ifil,bool frombeg,string idkey);
bool getTitleFromFileWFX(ifstream &ifil,int &nt,string* &tit);
bool getKeyWordsFromFileWFX(ifstream &ifil,string &kw);
void getNetCharge(ifstream &ifil,int &nc);
void getNofNucleiFromFileWFX(ifstream &ifil, int &nnuc);
void getNofPrimFromFileWFX(ifstream &ifil,int &npri);
void getNofEDFPrimFromFileWFX(ifstream &ifil,const int nedfc,int &nedfprim);
void getNofMolOrbFromFileWFX(ifstream &ifil,int &nmo);
void getNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,solreal* &rr);
void getNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,solreal** &rr3);
void getAtLabelsFromFileWFX(ifstream &ifil,const int nn,string* &al);
void getAtChargesFromFileWFX(ifstream &ifil,const int nn,solreal* &ach);
int getAtNumbersFromFileWFX(ifstream &ifil,const int nn,int* &anu);
void getNofElectronsFromFileWFX(ifstream &ifil,int &nel);
void getNofCoreElectronsFromFileWFX(ifstream &ifil,int &ncel);
void getPrimCentersFromFileWFX(ifstream &ifil,const int npr,int* &pc);
void getEDFPrimCentersFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pc);
void getPrimTypesFromFileWFX(ifstream &ifil,const int npr,int* &pt);
void getEDFPrimTypesFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,int* &pt);
void getPrimExponentsFromFileWFX(ifstream &ifil,const int npr,solreal* &pex);
void getEDFPrimExponentsFromFileWFX(ifstream &ifil,const int npr,\
      const int ntot,solreal* &pex);
void getMolecOrbOccNumsFromFileWFX(ifstream &ifil,const int nmo,solreal* &ocnu);
void getMolecOrbEnergiesFromFileWFX(ifstream &ifil,const int nmo,solreal* &orben);
void getMolecOrbCoefficientsFromFileWFX(ifstream &ifil,const int nmo,const int npr,solreal* &tcf);
void getTotEnerAndVirialFromFileWFX(ifstream &ifil,solreal &tote,solreal &vir);
void countEDFCentersFromFileWFX(ifstream &ifil,int &nedfc);
void getEDFExistenceFromFileWFX(ifstream &ifil,bool &ihaveEDF);
void getEDFPrimCoefficientsFromFileWFX(ifstream &ifil,\
      const int nedfp,solreal* &edfc);

#endif//_IOFUNCTS_WFX_H_

