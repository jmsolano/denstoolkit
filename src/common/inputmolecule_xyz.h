#ifndef _INPUTMOLECULE_XYZ_H_
#define _INPUTMOLECULE_XYZ_H_
#include <string>
using std::string;
#include "molecule.h"
#include <fstream>
using std::ifstream;

/* ************************************************************************** */
class InputMoleculeXYZ : public Molecule {
/* ************************************************************************** */
public:
   InputMoleculeXYZ();
   InputMoleculeXYZ(string fname);
   void ReadFromFile(string fname);
   void LoadCoordinatesNumbers(ifstream &ifil,int nat);
   void LoadCoordinatesSymbols(ifstream &ifil,int nat);
   void Save(const string &onam) const;
   void DisplayProperties();
   string title;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */
std::ostream &operator<<(std::ostream &out,const InputMoleculeXYZ (&mol));

#endif  /* _INPUTMOLECULE_XYZ_H_ */

