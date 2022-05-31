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
#ifndef _INPUTMOLECULE_CUB_H_
#define _INPUTMOLECULE_CUB_H_
#include <string>
using std::string;
#include "molecule.h"
#include <fstream>
using std::ifstream;

/* ************************************************************************** */
class InputMoleculeCub : public Molecule {
/* ************************************************************************** */
public:
   InputMoleculeCub();
   InputMoleculeCub(const string fname);
   ~InputMoleculeCub() {}
   /** Reads the molecule geometry from fname file. 
    * This opens and closes the file, internally. */
   bool ReadFromFile(const string fname);
   /** Reads the molecule geometry from the 
    * ifstream ifil. This will set the ifil buffer position
    * at 0 after loading the molecule geometry.*/
   void ReadFromFile(ifstream &ifil);
   void LoadCoordinatesNumbers(ifstream &ifil,int nat);
   /** As far as JMSA knows, cube files always contain
    * atomic numbers, however, this allows for cubes
    * wherein atomic symbols are used.  */
   void LoadCoordinatesSymbols(ifstream &ifil,int nat);
   void DisplayProperties();
   string title1,title2;
   vector<double> charge;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */
std::ostream &operator<<(std::ostream &out,const InputMoleculeCub (&mol));

#endif  /* _INPUTMOLECULE_CUB_H_ */

