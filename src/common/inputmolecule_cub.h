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

