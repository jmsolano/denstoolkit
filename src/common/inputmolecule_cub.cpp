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
#include "inputmolecule_cub.h"
#include "fileutils.h"
#include "screenutils.h"
#include "stringtools.h"
#include "unitconversion.h"

InputMoleculeCub::InputMoleculeCub() : Molecule() {
}
InputMoleculeCub::InputMoleculeCub(const string fname) : InputMoleculeCub() {
   ReadFromFile(fname);
}
void InputMoleculeCub::ReadFromFile(ifstream &ifil) {
   size_t initPos=ifil.tellg();
   std::getline(ifil,title1);
   std::getline(ifil,title2);
   int nat;
   ifil >> nat;
   charge.resize(nat);
   string symb;
   for ( size_t i=0 ; i<4 ; ++i ) {
      std::getline(ifil,symb);
   }
   StringTools::RemoveSpacesLeftAndRight(symb);
   string tmp=StringTools::GetFirstChunk(symb);
   int npts=std::stoi(tmp);
   bool inangstroms=false;
   if ( npts<0 ) { inangstroms=true; }
   size_t pos=ifil.tellg();
   ifil >> symb;
   ifil.seekg(pos);
   StringTools strtools;
   strtools.RemoveSpacesLeftAndRight(symb);
   if ( ScreenUtils::IsDigit(symb[0]) ) {
      LoadCoordinatesNumbers(ifil,nat,inangstroms);
   } else {
      LoadCoordinatesSymbols(ifil,nat,inangstroms);
   }
   ifil.seekg(initPos);
}
bool InputMoleculeCub::ReadFromFile(const string fname) {
   if ( !(FileUtils::ExtensionMatches(fname,"cub") ||
            FileUtils::ExtensionMatches(fname,"cube")) ) {
      return false;
   }
   ifstream ifil(fname.c_str());
   if ( !ifil.good() ) {
      ScreenUtils::DisplayErrorMessage(string("Could not open the file \"")+fname+string("\"!"));
      ifil.close();
      imsetup=false;
      return false;
   }
   ReadFromFile(ifil);
   ifil.close();
   imsetup=true;
   return true;
}
void InputMoleculeCub::LoadCoordinatesNumbers(ifstream &ifil,int nat,bool angs) {
   int kk;
   double convfact=1.0e0;
   if ( !angs ) { convfact=unitconv::bohr2angstrom; }
   vector<double> xt(3);
   for ( int i=0 ; i<nat ; ++i ) {
      ifil >> kk >> charge[i];//this line prevents reusing the code from xyz files.
      for ( int j=0 ; j<3 ; ++j ) { ifil >> xt[j]; xt[j]*=convfact; }
      AddAtom(xt,kk);
   }
}
void InputMoleculeCub::LoadCoordinatesSymbols(ifstream &ifil, int nat,bool angs) {
   string symb;
   double convfact=1.0e0;
   if ( !angs ) { convfact=unitconv::bohr2angstrom; }
   vector<double> xt(3);
   for ( int i=0 ; i<nat ; ++i ) {
      ifil >> symb >> charge[i];//this line prevents reusing the code from xyz files.
      for ( int j=0 ; j<3 ; ++j ) { ifil >> xt[j]; xt[j]*=convfact; }
      AddAtom(xt,symb);
   }
}
void InputMoleculeCub::DisplayProperties() {
   cout << title1 << endl;
   cout << title2 << endl;
   DisplayAtomProperties();
}
std::ostream &operator<<(std::ostream &out,const InputMoleculeCub (&mol)) {
   const Molecule* pmol=&mol;
   out << mol.title1 << endl;
   out << mol.title2 << endl;
   if ( mol.Size()>0 ) { out << pmol; } else {cout << "Empty molecule!" << endl;}
   return out;
}
