#ifndef _HELPERSNCI_H_
#define _HELPERSNCI_H_
#include <string>
using std::string;
#include "povraytools.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "isosurface.h"

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
         const string &palname="greens",bool render=false);
   static void CenterMolecule(BondNetWork &bn,Isosurface &iso);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSNCI_H_ */

