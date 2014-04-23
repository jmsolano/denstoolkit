/*
 *  iofuncts-wfn.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 25/03/13.
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef _IOFUNCTS_WFN_H_
#define _IOFUNCTS_WFN_H_

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

string getTitleFromFileWFN(ifstream &ifil);
void processFirstDataStringinWFNFile(ifstream &ifil,string* &tit,string & orbdesc,int &nmo,int &npr,int &nnu);
void processCentersWFN(ifstream &ifil,const int nnu,string* &atlbl,solreal* &rr,solreal* &atch);
void processPrimitivesWFN(ifstream &ifil,const int npr,int* &pricen,int* &primty,solreal* &prexp);
void processMolecularOrbitalPropsAndCoefs(ifstream &ifil,const int norb,const int npr
                                 ,solreal* &ocn,solreal* &moe,solreal* &moc);
void getEnergyAndVirial(ifstream &ifil,solreal &theener,solreal &thevir);
#endif//_IOFUNCTS_WFN_H_

