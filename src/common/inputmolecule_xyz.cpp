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
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;
#include <iomanip>
#include "inputmolecule_xyz.h"
#include "screenutils.h"
#include "stringtools.h"

InputMoleculeXYZ::InputMoleculeXYZ() : Molecule() {
}
InputMoleculeXYZ::InputMoleculeXYZ(string fname) : InputMoleculeXYZ() {
   ReadFromFile(fname);
}
void InputMoleculeXYZ::ReadFromFile(string fname) {
   ifstream ifil(fname.c_str());
   if ( !ifil.good() ) {
      ScreenUtils::DisplayErrorMessage(string("Could not open the file \"")+fname+string("\"!"));
      ifil.close();
      imsetup=false;
      return;
   }
   ReadFromFile(ifil);
   ifil.close();
   imsetup=true;
}
void InputMoleculeXYZ::ReadFromFile(ifstream &ifil) {
   string symb;
   std::getline(ifil,symb);
   int nat=std::stoi(symb);
   std::getline(ifil,title);
   if ( title.size()==0 || title==string("\r") ) {
      title=string("No previous title. XYZ created by InputMoleculeXYZ class.");
   }
   size_t pos=ifil.tellg();
   std::getline(ifil,symb);
   ifil.seekg(pos);
   StringTools strtools;
   strtools.RemoveSpacesLeftAndRight(symb);
   if ( ScreenUtils::IsDigit(symb[0]) ) {
      LoadCoordinatesNumbers(ifil,nat);
   } else {
      LoadCoordinatesSymbols(ifil,nat);
   } 
}
void InputMoleculeXYZ::LoadCoordinatesNumbers(ifstream &ifil,int nat) {
   int kk;
   vector<double> xt(3);
   for ( int i=0 ; i<nat ; ++i ) {
      ifil >> kk;
      for ( int j=0 ; j<3 ; ++j ) { ifil >> xt[j]; }
      AddAtom(xt,kk);
   }
}
void InputMoleculeXYZ::LoadCoordinatesSymbols(ifstream &ifil, int nat) {
   string symb;
   vector<double> xt(3);
   for ( int i=0 ; i<nat ; ++i ) {
      ifil >> symb;
      for ( int j=0 ; j<3 ; ++j ) { ifil >> xt[j]; }
      AddAtom(xt,symb);
   }
}
void InputMoleculeXYZ::DisplayProperties() {
   cout << title << endl;
   DisplayAtomProperties();
}
void InputMoleculeXYZ::Save(const string &onam) const {
   ofstream ofil(onam.c_str());
   ofil << Size() << endl;
   ofil << title << endl;
   ofil << std::scientific << std::setprecision(10);
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      ofil << std::setw(2) << atom[i].symbol
           << std::setw(20) << atom[i].x[0]
           << std::setw(20) << atom[i].x[1]
           << std::setw(20) << atom[i].x[2]
           << endl;
   }
   ofil.close();
}
std::ostream &operator<<(std::ostream &out,const InputMoleculeXYZ (&mol)) {
   const Molecule* pmol=&mol;
   out << mol.title << endl;
   out << pmol;
   return out;
}

