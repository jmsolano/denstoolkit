#ifndef _COMMONHELPERS_H_
#define _COMMONHELPERS_H_
#include "bondnetwork.h"
/* ************************************************************************** */
class CommonHelpers {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   /** Adds the nuclei spheres to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmat is the string to be passed for the transparency of
    * the atoms. */
   static void PutNuclei(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmat);
   /** Adds the bond cylinders to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmbnd is the string to be passed for the transparency of
    * the bonds. */
   static void PutBonds(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmbnd);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */
#endif  /* _COMMONHELPERS_H_ */


