#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cmath>
#include <string>
using std::string;
#include "molecule.h"

int main (int argc, char *argv[]) {
   bool verbose=false;
   if ( argc==2 ) {
      if ( string(argv[1])==string("-v") ) {
         verbose=true;
      }
   }
   bool passed=true;
   if ( verbose ) {
      cout << "************************************************" << endl;
      cout << "************************************************" << endl;
   }
   Molecule molecule;
   vector<double> x(3); x[1]=1.0e0;
   molecule.AddAtom(x,3);
   passed=passed&&(molecule.Size()==1);
   //---------------------------------------
   cout << (passed? "PASSED" : "FAILED") << endl;
   //---------------------------------------
   if ( verbose ) {
      x[1]=-1.0e0; molecule.AddAtom(x,1);
      cout << molecule << endl;
   }
   return EXIT_SUCCESS;
}



