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
   static bool MakeCriticalPointsPovFile(const string &basename,const OptionFlags &options,\
         char *argv[],BondNetWork &bn,const string cpname);
   static bool AlignMolecule(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]);
   static bool AlignMolecule3Atoms(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]);
   static void SetupPovConf(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char* argv[]);
   static void WriteCameraSettings(ofstream &ofil,POVRayConfiguration &pvc,\
         bool lookat=false);
   static void WriteLightSources(ofstream &ofil,POVRayConfiguration &pvc);
   static void WritePovHeader(ofstream &ofil);
   static void WriteBondNetworkElements(ofstream &ofil,BondNetWork &bn,\
         POVRayConfiguration &pvc,const OptionFlags &options,char *argv[]);
   static bool DrawCPofType(char t,const OptionFlags &options,char *argv[]);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSDRAWER_H_ */

