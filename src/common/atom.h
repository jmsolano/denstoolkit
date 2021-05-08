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
#ifndef _ATOM_H_
#define _ATOM_H_
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <memory>

#define MAXATNUMDEF 109
#define COORDSEPSILON 5.0e-01

/* ************************************************************************** */
/** This class holds the properties of a single atom. The variables
 * x, symbol, weight, and name are public, in order to make them
 * accessible without Getters.  */
class Atom {
/* ************************************************************************** */
public:
   Atom(int an);
   Atom(vector<double> &ux,string &usymb);
   Atom(vector<double> &ux,int an);
   Atom(const Atom &p);
   Atom& operator=(const Atom& other);
   void SetupAtom(int an);
   /* ************************************************************************** */
   /** Returns the atomic weight of the atom whose atomic number is n
    * Notice: it must be the atomic number, not the index!  */
   static double GetAtomicWeight(int n);
   /** Returns the name of the atom whose atomic number is n
    * Notice: it must be the atomic number, not the index!  */
   static string GetName(int n);
   /** Returns the atomic symbol of the atom whose atomic number is n
    * Notice: it must be the atomic number, not the index!  */
   static string GetAtomicSymbol(int n);
   static int GetAtomicNumberFromSymbol(string smb);
   /** Returns the valence electrons of the atom whose atomic
    * number is n. In the current version, only neutral atoms
    * are treated, i.e. ions are not considered.  */
   static int GetValenceElectrons(int n);
   static int GetValenceElectrons(string s) { return GetValenceElectrons(GetAtomicNumberFromSymbol(s));}
   bool IsMyPosition(vector<double> &xp);
   void DisplayProperties();
   /* ************************************************************************** */
   vector<double> x;
   string symbol;
   string name;
   double weight;
   int num;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   Atom();
   void Init();
   static vector<string> tab_symbol;
};
/* ************************************************************************** */
std::ostream &operator<<(std::ostream &out,const Atom (&atom));

#endif  /* _ATOM_H_ */

