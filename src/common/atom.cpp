#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <iomanip>
#include <cmath>
#include "atom.h"
#include "screenutils.h"

vector<string> Atom::tab_symbol={
      "H" , "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
      "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
      "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
      "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
      "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
      "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
      "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
      "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
      "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt" };
/* Most of the hard-coded data was taken from:
 * https://www.science.co.il/elements/
 * Last access: 2018-Apr-18  */
Atom::Atom() {
   Init();
}
Atom::Atom(vector<double> &ux,int an) : Atom(an) {
   for ( int i=0 ; i<3 ; ++i ) { x[i]=ux[i]; }
}
Atom::Atom(vector<double> &ux,string &usymb) : Atom() {
   int nn=GetAtomicNumberFromSymbol(usymb);
   for ( int i=0 ; i<3 ; ++i ) { x[i]=ux[i]; }
   SetupAtom(nn);
}
Atom::Atom(int an) : Atom() {
   for ( int i=0 ; i<3 ; ++i ) { x[i]=0.0e0; }
   SetupAtom(an);
}
Atom::Atom(const Atom &other) : Atom() {
      x=other.x;
      symbol=other.symbol;
      name=other.name;
      weight=other.weight;
      num=other.num;
}
void Atom::Init() {
   x.resize(3);
   symbol="";
   name="";
   weight=0.0e0;
   num=0;
}
Atom& Atom::operator=(const Atom& other) {
   if (this != &other) {
      x=other.x;
      symbol=other.symbol;
      name=other.name;
      weight=other.weight;
      num=other.num;
   }
   return *this;
}
void Atom::SetupAtom(int an) {
   num=an;
   weight=GetAtomicWeight(num);
   name=GetName(num);
   symbol=GetAtomicSymbol(num);
}
double Atom::GetAtomicWeight(int n) {
   static double tab_weight[MAXATNUMDEF] = {
      1.008e0, 4.003e0, 6.941e0, 9.012e0, 10.811e0,
      12.011e0, 14.007e0, 15.999e0, 18.998e0, 20.180e0,
      22.990e0, 24.305e0, 26.982e0, 28.086e0, 30.974e0,
      32.065e0, 35.453e0, 39.948e0, 39.098e0, 40.078e0,
      44.956e0, 47.867e0, 50.942e0, 51.996e0, 54.938e0,
      55.845e0, 58.933e0, 58.693e0, 63.546e0, 65.390e0,
      69.723e0, 72.640e0, 74.922e0, 78.960e0, 79.904e0,
      83.800e0, 85.468e0, 87.620e0, 88.906e0, 91.224e0,
      92.906e0, 95.940e0, 98.000e0, 101.070e0, 102.906e0,
      106.420e0, 107.868e0, 112.411e0, 114.818e0, 118.710e0,
      121.760e0, 127.600e0, 126.905e0, 131.293e0, 132.906e0,
      137.327e0, 138.906e0, 140.116e0, 140.908e0, 144.240e0,
      145.000e0, 150.360e0, 151.964e0, 157.250e0, 158.925e0,
      162.500e0, 164.930e0, 167.259e0, 168.934e0, 173.040e0,
      174.967e0, 178.490e0, 180.948e0, 183.840e0, 186.207e0,
      190.230e0, 192.217e0, 195.078e0, 196.967e0, 200.590e0,
      204.383e0, 207.200e0, 208.980e0, 209.000e0, 210.000e0,
      222.000e0, 223.000e0, 226.000e0, 227.000e0, 232.038e0,
      231.036e0, 238.029e0, 237.000e0, 244.000e0, 243.000e0,
      247.000e0, 247.000e0, 251.000e0, 252.000e0, 257.000e0,
      258.000e0, 259.000e0, 262.000e0, 261.000e0, 262.000e0,
      266.000e0, 264.000e0, 277.000e0, 268.000e0
   };
   if ( n<1 ) {
      ScreenUtils::DisplayErrorMessage("n<1!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return 0.0e0;
   }
   return tab_weight[n-1];
}
string Atom::GetName(int n) {
   static string tab_name[MAXATNUMDEF]={
      "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
      "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
      "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorus",
      "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
      "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese",
      "Iron", "Cobalt", "Nickel", "Copper", "Zinc",
      "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine",
      "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium",
      "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium",
      "Palladium", "Silver", "Cadmium", "Indium", "Tin",
      "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium",
      "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
      "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium",
      "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium",
      "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium",
      "Osmium", "Iridium", "Platinum", "Gold", "Mercury",
      "Thallium", "Lead", "Bismuth", "Polonium", "Astatine",
      "Radon", "Francium", "Radium", "Actinium", "Thorium",
      "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium",
      "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium",
      "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium",
      "Seaborgium", "Bohrium", "Hassium", "Meitnerium"
   };
   if ( n<1 ) {
      ScreenUtils::DisplayErrorMessage("n<1!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return "Error";
   }
   return tab_name[n-1];
}
string Atom::GetAtomicSymbol(int n) {
   if ( n<1 ) {
      ScreenUtils::DisplayErrorMessage("n<1!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return "Error";
   }
   return tab_symbol[n-1];
}
int Atom::GetValenceElectrons(int n) {
   if ( n<=2 ) {
      return n;
   } else if ( n<=10 ) {
      return n-2;
   } else if ( n<=18 ) {
      return n-10;
   } else {
      ScreenUtils::DisplayErrorMessage("Only atoms of the first three rows of "\
            "the periodic table are implemented");
   }
   return -1;
}
bool Atom::IsMyPosition(vector<double> &xp) {
   if ( xp.size()!=3 ) {
#if DEBUG
      ScreenUtils::DisplayErrorMessage("vector size is != 3 !");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
      return false;
   }
   bool res=true;
   for ( int i=0 ; i<3 ; ++i ) {
      res=res&&(fabs(xp[i]-x[i])<COORDSEPSILON);
   }
   return res;
}
int Atom::GetAtomicNumberFromSymbol(string smb) {
   int k=0;
   while(smb!=tab_symbol[k] && k<MAXATNUMDEF) { ++k; }
   return k+1;
}
void Atom::DisplayProperties() {
   cout << "Symbol: " << symbol << endl;
   cout << "Number: " << num << endl;
   cout << "  Name: " << name << endl;
   cout << "Weight: " << weight << endl;
   cout << "     x: " << x[0] << " " << x[1] << " " << x[2] << endl;
}
std::ostream &operator<<(std::ostream &out,const Atom (&atm)) {
   return out << std::setw(2) << atm.symbol << " (" << atm.num
              << ", '" << atm.name << "', " << atm.weight << "): "
              << atm.x[0] << " " << atm.x[1] << " " << atm.x[2];
}

