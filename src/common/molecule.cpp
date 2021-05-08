/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
#ifndef _MOLECULE_CPP_
#define _MOLECULE_CPP_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <algorithm>
#include <cmath>
#include <cassert>
#include "molecule.h"
#include "matrixvectoroperations3d.h"

Molecule::Molecule() {
   Init();
}
void Molecule::Init() {
   atom.clear();
   imsetup=false;
}
Molecule::~Molecule() {
}
void Molecule::AddAtom(vector<double> &ux,int an) {
   atom.push_back(Atom(ux,an));
}
void Molecule::AddAtom(vector<double> &ux,string &usymb) {
   atom.push_back(Atom(ux,usymb));
}
void Molecule::DisplayAtomProperties() {
   size_t k=atom.size();
   for ( size_t i=0 ; i<k ; ++i ) { atom[i].DisplayProperties(); }
}
void Molecule::DisplayProperties() {
   DisplayAtomProperties();
}
void Molecule::SortCoordinates() {
   //Sorts first column
   QuickSort(0);
   //Sorts second column...
   size_t n=atom.size()-1;
   size_t count,currpos;
   for ( size_t i=0 ; i<n ; ++i ) {
      count=0;
      currpos=i;
      //cout << "currpos: " << currpos << "(" << i << ")";
      while ( currpos<n && (fabs(atom[currpos+1].x[0]-atom[currpos].x[0])<=SINGLECOORDEPS) ) { //see globaldefs.h
         ++count;
         ++currpos;
      }
      //cout << ", count: " << count << endl;
      if ( count!=0 ) {
         QuickSort(i+count,i,1);
         i+=count;
      }
   }
   //Sorts third column
   for ( size_t i=0 ; i<n ; ++i ) {
      count=0;
      currpos=i;
      //cout << "currpos: " << currpos << "(" << i << ")";
      while ( currpos<n &&\
              (fabs(atom[currpos+1].x[1]-atom[currpos].x[1])<=SINGLECOORDEPS) &&\
              (fabs(atom[currpos+1].x[0]-atom[currpos].x[0])<=SINGLECOORDEPS) ) { //see globaldefs.h
         ++count;
         ++currpos;
      }
      //cout << ", count: " << count << endl;
      if ( count!=0 ) {
         QuickSort(i+count,i,2);
         i+=count;
      }
   }
}
void Molecule::QuickSort(int high,int low,int srtIdx) {
   int i=low;
   int j=high;
   int hplo2=(high+low)/2;
   double mid[3];
   mid[0]=atom[hplo2].x[0];
   mid[1]=atom[hplo2].x[1];
   mid[2]=atom[hplo2].x[2];
   do {
      while (atom[i].x[srtIdx]<mid[srtIdx]) { ++i; }
      while (atom[j].x[srtIdx]>mid[srtIdx]) { --j; }
      if ( i<=j ) {
         std::swap(atom[i],atom[j]);
         ++i;
         --j;
      }
   } while ( i<j );
   if ( low<j ) { QuickSort(j,low,srtIdx); }
   if ( i<high ) { QuickSort(high,i,srtIdx); }
}
string Molecule::EmpiricalFormula() const {
   size_t n=MAXATNUMDEF; //see globaldefs.h
   vector<int> count(n);
   for ( size_t i=0 ; i<n ; ++i ) { count[i]=0; }
   for ( size_t i=0 ; i<atom.size() ; ++i ) { count[atom[i].num-1]++; }
   string res;
   for ( size_t i=0 ; i<n ; ++i ) {
      if ( count[i]>0 ) {
         res+=Atom::GetAtomicSymbol(i+1);
         res+=std::to_string(count[i]);
      }
   }
   return res;
}
int Molecule::CountAtomsOfType(string ss) {
   int atnum=Atom::GetAtomicNumberFromSymbol(ss);
   return CountAtomsOfType(atnum);
}
int Molecule::CountAtomsOfType(int nn) {
   int count=0;
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      if ( atom[i].num==nn ) { ++count; }
   }
   return count;
}
std::ostream &operator<<(std::ostream &out,const Molecule (&mol)) {
   out << "The molecule has " << mol.Size() << " atoms." << endl;
   for ( size_t i=0 ; i<(mol.Size()-1) ; ++i ) {
      out << mol.atom[i] << endl;
   }
   out << mol.atom[mol.Size()-1] << endl;
   out << "Empirical formula: " << mol.EmpiricalFormula();
   return out;
}
std::ostream &operator<<(std::ostream &out,const Molecule* mol) {
   out << *mol;
   return out;
}
bool operator== (const Molecule &m1, const Molecule &m2) {
   size_t n1=m1.Size();
   if ( n1!=m2.Size() ) { return false; }
   bool tmpbool;
   vector<double> tmpvec;
   double d;
   int nat;
   for ( size_t i=0 ; i<n1 ; ++i ) {
      tmpvec=m1.atom[i].x;
      tmpbool=false;
      nat=m1.atom[i].num;
      for ( size_t j=0 ; j<n1 ; ++j ) {
         d=MatrixVectorOperations3D::Distance(tmpvec,m2.atom[j].x);
         if ( d<COORDSEPSILON && (nat==m2.atom[j].num) ) {
            tmpbool=true;
            break;
         }
      }
      if ( !tmpbool ) { return false; }
   }
   return true;
}
void DuplicateAtomOrder(Molecule &m1,Molecule &m2) {
   size_t n1=m1.Size();
   vector<double> tmpvec;
   double d;
   int nat;
   for ( size_t i=0 ; i<n1 ; ++i ) {
      tmpvec=m1.atom[i].x;
      nat=m1.atom[i].num;
      for ( size_t j=0 ; j<n1 ; ++j ) {
         d=MatrixVectorOperations3D::Distance(tmpvec,m2.atom[j].x);
         if ( d<COORDSEPSILON && (nat==m2.atom[j].num) ) {
            std::swap(m2.atom[i],m2.atom[j]);
            break;
         }
      }
   }
}

#endif  /* _MOLECULE_CPP_ */

