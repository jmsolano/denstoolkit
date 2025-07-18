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
#ifndef _BONDNETWORK_H_
#define _BONDNETWORK_H_

#define BOHRTOANGSTROM (0.529177249e0)
#define VDWFACTOR (1.2e0)
#define AUTOMATICDRAWATOMSIZE (0.5e0)
#define AUTOMATICBALLANDSTICKRATIO (0.4e0)
#define AUTOMATICSPACEFILLINGRATIO (3.3e0)
#define BOUNDINGBOXSCALEFACTOR (1.2e0)
#define AUTOMATICMAXBONDDIST (3.0e0)

#include <string>
using std::string;
#include "screenutils.h"
#include "mymemory.h"
#include "povraytools.h"

#define MAXBONDINGATOMS 8

class BondNetWork {
   /* ************************************************************************** */
public:
/* ************************************************************************** */
   BondNetWork(); /* Default constructor */
   BondNetWork(int nn,double* &rin);
   ~BondNetWork(); /* Destructor */
   bool ReadFromFileWFX(string inname);
   bool ReadFromFileWFN(string inname);
   bool ReadFromFileXYZ(string inname);
   bool ReadFromFileCPX(string inname);
   bool ReadFromFile(string inname);
   /* ************************************************************************** */
   double **R;         //Contains the atoms' radius-vectors.
   double **bondDist;  //It contains the bond distances of the actual bonded atoms. It will
                     //  not be allocated by its own, but rather one needs to call some setup
                     //  function.
   int nNuc;         //The number of nuclei.
   int *atNum;       //The atomic number of each nuclei atNum(H)=0
   string *atLbl;    //The atom labels
   int **bNet;       //The bonding network, contains the labels of the
                     //   atoms that each atom is bonded to.
   int nBonds;
   int nTit;
   string *title;
   double drawAtSize,drawStickSize;
   double rmax[3],rmin[3],rView;
   double bbmax[3],bbmin[3],maxBondDist;
   bool ballAndStickMode;
   bool spaceFillingMode;
   bool wireMode;
   /* ************************************************************************** */
   bool SetUpBNW(void);
   bool ImStp(void);
   bool LookForBonds(void);
   /* ************************************************************************** */
   double Dist(int i, int k);
   void AddBond(int i,int j,double dd);
   bool MakePOVFile(string pnam, POVRayConfiguration &pvp);
   void SeekRMaxMin(void);
   void CalcViewRadius(void);
   void PutNuclei(ofstream &pof);
   void PutBonds(ofstream &pof);
   void CenterMolecule(void);
   void SetBoundingBox(void);
   int CountAtomsOfAtomicNumber(int nat);
/* ************************************************************************** */
private:
   bool isSTP;       //Just to ensure that the bond network has been set up.
};

#endif//_BONDNETWORK_H_

