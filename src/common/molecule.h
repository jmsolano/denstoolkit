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
#ifndef _MOLECULE_H_
#define _MOLECULE_H_
#include <vector>
using std::vector;
#include <string>
using std::string;
#include "atom.h"

#ifndef SINGLECOORDEPS
#define SINGLECOORDEPS 1.0e-04
#endif

/* ************************************************************************** */
class Molecule {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   /** The Molecule class contains the basic information of a chemical compound.
    * Basically, it is a set of atoms, each of which is located at a certain
    * 3D position. Internally, Molecule will assume that the coordinates
    * are expressed in Angstroms.  */
   Molecule();
   virtual ~Molecule();
/* ************************************************************************** */
   void AddAtom(vector<double> &ux,int an);
   void AddAtom(vector<double> &ux,string &usymb);
   size_t Size() const {return atom.size();}
   void DisplayAtomProperties();
   virtual void DisplayProperties();
   string EmpiricalFormula() const;
   int CountAtomsOfType(const char* cc) {return CountAtomsOfType(string(cc));}
   int CountAtomsOfType(string ss);
   int CountAtomsOfType(int nn);
   void SortCoordinates();
   /** Determine the bounding box of the molecule.  */
   void DetermineBoundingBox();
   void ComputeCenterOfMass();
   void ComputeCentroid();
   void CenterAtCentroid();
   void CenterAtCenterOfMass();
   /** Sets the center of the coordinate system as in the original reading/loading
    * (the molecule is translated to origCent).  */
   void ResetOriginOfCoordinates();
   void SetupBonds();
   void SetupCells();
   void DisplayBondProperties();
   /** Returns a vector of indices. The first element of the vector
    * is the observed atom, and the following elements are the indices of
    * the observed atom neighbours.  */
   vector<size_t> ListOfNeighbours(const size_t atpos);
/* ************************************************************************** */
   vector<Atom> atom;
   vector<double> cm; /*!< Center of mass  */
   vector<double> cd; /*!< Centroid  */
   vector<double> xmin; /** Lower, back, left bounding box.  */
   vector<double> xmax; /** Upper, front, right bounding box. */
   double rmax; /*!< Holds the maximum atom distance from the center of mass.  */
   vector<vector<int> > bond;
   vector<vector<double> > bndDist;
   double maxBondDist;
   vector<vector<vector<vector<int > > > > cell;
   static double constexpr cellLen=5.0e0;
   bool ImSetup() const;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void Init();
   void QuickSort(int srtIdx) {return QuickSort((atom.size()-1),0,srtIdx);}
   void QuickSort(int high, int low,int srtIdx=0);
   bool imsetup;
   vector<double> origCent; /*!< Saves the initial origin of the coordinate system.  */
/* ************************************************************************** */
   friend bool operator== (const Molecule &m1, const Molecule &m2);
/* ************************************************************************** */
};
/* ************************************************************************** */
std::ostream &operator<<(std::ostream &out,const Molecule (&mol));
std::ostream &operator<<(std::ostream &out,const Molecule* mol);
bool operator== (const Molecule &m1, const Molecule &m2);
void DuplicateAtomOrder(Molecule &m1,Molecule &m2);

#endif  /* _MOLECULE_H_ */

