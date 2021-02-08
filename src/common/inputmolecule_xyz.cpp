#ifndef _INPUTMOLECULE_XYZ_CPP_
#define _INPUTMOLECULE_XYZ_CPP_

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
   ifil.close();
   imsetup=true;
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

#endif  /* _INPUTMOLECULE_XYZ_CPP_ */

