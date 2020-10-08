#ifndef _INPUTMOLECULE_CUB_CPP_
#define _INPUTMOLECULE_CUB_CPP_

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
#include "screenutils.h"
#include "stringtools.h"

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
   size_t pos=ifil.tellg();
   ifil >> symb;
   ifil.seekg(pos);
   StringTools strtools;
   strtools.RemoveSpacesLeftAndRight(symb);
   if ( ScreenUtils::IsDigit(symb[0]) ) {
      LoadCoordinatesNumbers(ifil,nat);
   } else {
      LoadCoordinatesSymbols(ifil,nat);
   }
   ifil.seekg(initPos);
}
bool InputMoleculeCub::ReadFromFile(const string fname) {
   size_t pos=fname.find_last_of(".");
   if ( fname.substr(pos,3)!=string("cub") ) {
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
void InputMoleculeCub::LoadCoordinatesNumbers(ifstream &ifil,int nat) {
   int kk;
   vector<double> xt(3);
   for ( int i=0 ; i<nat ; ++i ) {
      ifil >> kk >> charge[i];//this line prevents reusing the code from xyz files.
      for ( int j=0 ; j<3 ; ++j ) { ifil >> xt[j]; }
      AddAtom(xt,kk);
   }
}
void InputMoleculeCub::LoadCoordinatesSymbols(ifstream &ifil, int nat) {
   string symb;
   vector<double> xt(3);
   for ( int i=0 ; i<nat ; ++i ) {
      ifil >> symb >> charge[i];//this line prevents reusing the code from xyz files.
      for ( int j=0 ; j<3 ; ++j ) { ifil >> xt[j]; }
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
#endif  /* _INPUTMOLECULE_CUB_CPP_ */

