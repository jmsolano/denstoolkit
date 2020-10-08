#ifndef _MOLECULE_H_
#define _MOLECULE_H_
#include <vector>
using std::vector;
#include <string>
using std::string;
#include "atom.h"

#ifndef SINGLECOORDEPS
#define SINGLECOORDEPS 1.0e-04
#endif

/* ************************************************************************** */
class Molecule {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Molecule();
   virtual ~Molecule();
/* ************************************************************************** */
   void AddAtom(vector<double> &ux,int an);
   void AddAtom(vector<double> &ux,string &usymb);
   size_t Size() const {return atom.size();}
   void DisplayAtomProperties();
   virtual void DisplayProperties();
   string EmpiricalFormula() const;
   int CountAtomsOfType(const char* cc) {return CountAtomsOfType(string(cc));}
   int CountAtomsOfType(string ss);
   int CountAtomsOfType(int nn);
   void SortCoordinates();
   bool ImSetup() {return imsetup;}
/* ************************************************************************** */
   vector<Atom> atom;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void Init();
   void QuickSort(int srtIdx) {return QuickSort((atom.size()-1),0,srtIdx);}
   void QuickSort(int high, int low,int srtIdx=0);
   bool imsetup;
/* ************************************************************************** */
   friend bool operator== (const Molecule &m1, const Molecule &m2);
/* ************************************************************************** */
};
/* ************************************************************************** */
std::ostream &operator<<(std::ostream &out,const Molecule (&mol));
std::ostream &operator<<(std::ostream &out,const Molecule* mol);
bool operator== (const Molecule &m1, const Molecule &m2);
void DuplicateAtomOrder(Molecule &m1,Molecule &m2);

#endif  /* _MOLECULE_H_ */

