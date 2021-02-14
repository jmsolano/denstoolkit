#ifndef _HELPERSNCI_H_
#define _HELPERSNCI_H_
#include <string>
using std::string;
#include "povraytools.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "isosurface.h"
#include "optflags.h"

/* ************************************************************************** */
class HelpersNCI {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static bool ComputeLambdaOnCentroids(GaussWaveFunction &wf,Isosurface &iso);
   static bool ComputeLambdaOnVertices(GaussWaveFunction &wf,Isosurface &iso);
   static bool ComputeNormalsAtVertices(GaussWaveFunction &wf,Isosurface &iso);
   static bool ComputeTriangleNormals(GaussWaveFunction &wf,Isosurface &iso);
   static bool MakePovFile(const string &povname,POVRayConfiguration &pvp,BondNetWork &bn,Isosurface &iso,\
         OptionFlags &options,char *argv[]);
   static void CenterMolecule(BondNetWork &bn,Isosurface &iso);
   static bool AlignMolecule(POVRayConfiguration &pvp,BondNetWork &bn,const OptionFlags &options,char *argv[]);
   static bool AlignMolecule3Atoms(POVRayConfiguration &pvp,BondNetWork &bn,const OptionFlags &options,char *argv[]);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSNCI_H_ */

