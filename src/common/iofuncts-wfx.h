/*
 *  iofuncts-wfx.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 25/03/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
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
void getNofNucleiFromFileWFX(ifstream &ifil, int &nnuc);
void getNofPrimFromFileWFX(ifstream &ifil,int &npri);
void getNofMolOrbFromFileWFX(ifstream &ifil,int &nmo);
void getNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,solreal* &rr);
void getNucCartCoordsFromFileWFX(ifstream &ifil,const int nn,solreal** &rr3);
void getAtLabelsFromFileWFX(ifstream &ifil,const int nn,string* &al);
void getAtChargesFromFileWFX(ifstream &ifil,const int nn,solreal* &ach);
void getAtNumbersFromFileWFX(ifstream &ifil,const int nn,int* &anu);
void getPrimCentersFromFileWFX(ifstream &ifil,const int npr,int* &pc);
void getPrimTypesFromFileWFX(ifstream &ifil,const int npr,int* &pt);
void getPrimExponentsFromFileWFX(ifstream &ifil,const int npr,solreal* &pex);
void getMolecOrbOccNumsFromFileWFX(ifstream &ifil,const int nmo,solreal* &ocnu);
void getMolecOrbEnergiesFromFileWFX(ifstream &ifil,const int nmo,solreal* &orben);
void getMolecOrbCoefficientsFromFileWFX(ifstream &ifil,const int nmo,const int npr,solreal* &tcf);
void getTotEnerAndVirialFromFileWFX(ifstream &ifil,solreal &tote,solreal &vir);

#endif//_IOFUNCTS_WFX_H_

