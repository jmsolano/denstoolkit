

#ifndef _RADII_SCHEME_JMHP_H_
#define _RADII_SCHEME_JMHP_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#define _HAVE_SELECTED_ATOM_RADII_ 1

#define MAXDEFINEDATOMICRADII 112

solreal getAtomicVDWRadius(int atn);

#endif /* defined(_RADII_SCHEME_JMHP_H_) */
