

#ifndef _COL_SCHEME_JMOL_H_
#define _COL_SCHEME_JMOL_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#define _HAVE_SELECTED_ATOM_PALETTE_ 1

solreal getAtomicRGBColorReal(int nat,int rgb);
solreal getAtomicRColorReal(int nat);
solreal getAtomicGColorReal(int nat);
solreal getAtomicBColorReal(int nat);

int getAtomicRGBColorInt(int nat,int rgb);
int getAtomicRColorInt(int nat);
int getAtomicGColorInt(int nat);
int getAtomicBColorInt(int nat);

void getAtomicRGBColorsReal(int nat,solreal &rr,solreal &gg,solreal &bb);
void getAtomicRGBColorsInt(int nat,int &rr,int &gg,int &bb);


#endif /* defined(_COL_SCHEME_JMOL_H_) */
