#ifndef _HELPERSDRAWER_H_
#define _HELPERSDRAWER_H_
#include <string>
using std::string;
#include "optflags.h"
#include "bondnetwork.h"

/* ************************************************************************** */
class HelpersDrawer {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static bool MakePovFile(const string &basename,const OptionFlags &options,\
         char *argv[],BondNetWork &bn);
   static bool AlignMolecule(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]);
   static bool AlignMolecule3Atoms(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSDRAWER_H_ */

