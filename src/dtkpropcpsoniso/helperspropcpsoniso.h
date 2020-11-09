#ifndef _HELPERSPROPCPSONISO_H_
#define _HELPERSPROPCPSONISO_H_
#include <vector>
using std::vector;
#include "optflags.h"
#include "bondnetwork.h"
#include "symmetricsurfacegrid.h"
#include "gausswavefunction.h"

/* ************************************************************************** */
class HelpersPropCPsOnIso {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static bool HaveIncompatibleOptions(int argc,char *argv[],OptionFlags &opt);
   static void GetCenterIndexAndVectors(char *argv[],const OptionFlags &options,\
         const BondNetWork &bn,int &cat,vector<double> &xc,vector<double> &xd);
   static bool MakePovFile(const string &povname,OptionFlags &options,POVRayConfiguration &pvp,\
         BondNetWork &bn,SymmetricSurfaceGrid &grid,vector<vector<double> > &sp);
   static void ProjectGridOntoIsosurface(GaussWaveFunction &wf,SymmetricSurfaceGrid &g,\
         const char prop,const double iso);
   /** Looks for the isovalue along a line. The line passes through c and r0. c is
    * the reference point, and r0 is the first point at which the field is evaluated.
    * Subsequent points are ri=r0+ih(r0-c) (i=1,2,...).
    * This function assumes that the field values decrease as the point
    * moves away from r0 (along the direction r0-c).  */
   static double SearchValueAlongLineDescending(const vector<double> &c,const vector<double> &r0,\
         GaussWaveFunction &wf,vector<double> &ri,const char prop,const double iso);
   /** This function searches the critical points (local minima and local maxima)
    * using the centroid of each face (projected onto the isosurface) and
    * the vertices of each face to determine whether v(rctd) is a critical point.
    * This function must be called after calling ProjectGridOntoIsosurface  */
   static bool SearchCPs(SymmetricSurfaceGrid &g,GaussWaveFunction &wf,\
         vector<vector<double> > &rcp,vector<size_t> &poscp,vector<int> &sigcp,\
         vector<double> &valcp,const char prop='V');
   static bool ComputeNormalsAtVertices(SymmetricSurfaceGrid &g,GaussWaveFunction &wf,const char prop='d');
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSPROPCPSONISO_H_ */

