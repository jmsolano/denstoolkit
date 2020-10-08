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

