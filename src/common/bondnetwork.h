/*
 *  bondnetwork.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 03/04/13.
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef _BONDNETWORK_H_
#define _BONDNETWORK_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#define BOHRTOANGSTROM (0.529177249e0)
#define VDWFACTOR (1.2e0)
#define AUTOMATICDRAWATOMSIZE (0.5e0)
#define AUTOMATICBALLANDSTICKRATIO (0.4e0)
#define AUTOMATICSPACEFILLINGRATIO (3.3e0)
#define BOUNDINGBOXSCALEFACTOR (1.2e0)
#define AUTOMATICMAXBONDDIST (3.0e0)

#include <iostream>
using std::cout;
using std::cin;
using std::endl;
using std::ios;
#include <fstream>
using std::fstream;
using std::ifstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <math.h>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;

#include "solscrutils.h"
#include "solmemhand.h"
#include "solpovtools.h"
#define MAXBONDINGATOMS 8

//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
class bondNetWork
//**********************************************************************************************
{
public:
   //**********************************************************************************************
   bondNetWork(); /* Default constructor */
   bondNetWork(int nn,solreal* &rin);
   ~bondNetWork(); /* Destructor */
   bool readFromFileWFX(string inname);
   bool readFromFileWFN(string inname);
   bool readFromFile(string inname);
   //**********************************************************************************************
   solreal **R;         //Contains the atoms' radius-vectors.
   solreal **bondDist;  //It contains the bond distances of the actual bonded atoms. It will
                     //  not be allocated by its own, but rather one needs to call some setup
                     //  function.
   int nNuc;         //The number of nuclei.
   int *atNum;       //The atomic number of each nuclei
   string *atLbl;    //The atom labels
   int **bNet;       //The bonding network, contains the labels of the
                     //   atoms that each atom is bonded to.
   int nBonds;
   int nTit;
   string *title;
   solreal drawAtSize,drawStickSize;
   solreal rmax[3],rmin[3],rView;
   solreal bbmax[3],bbmin[3],maxBondDist;
   bool ballAndStickMode;
   bool spaceFillingMode;
   bool wireMode;
   //**********************************************************************************************
   bool setUpBNW(void);
   bool imstp(void);
   bool lookForBonds(void);
   //**********************************************************************************************
   solreal dist(int i, int k);
   //**********************************************************************************************
   void addBond(int i,int j,solreal dd);
   //**********************************************************************************************
   bool makePOVFile(string pnam, povRayConfProp &pvp);
   //**********************************************************************************************
   void seekRMaxMin(void);
   //**********************************************************************************************
   void calcViewRadius(void);
   //**********************************************************************************************
   void putNuclei(ofstream &pof);
   //**********************************************************************************************
   void putBonds(ofstream &pof);
   //**********************************************************************************************
   void centerMolecule(void);
   //**********************************************************************************************
   void setBoundingBox(void);
   //**********************************************************************************************
   int countAtomsOfAtomicNumber(int nat);
   //**********************************************************************************************
private:
   bool isSTP;       //Just to ensure that the bond network has been set up.
   //**********************************************************************************************
};

#endif//_BONDNETWORK_H_



