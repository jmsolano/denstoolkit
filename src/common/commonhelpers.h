#ifndef _COMMONHELPERS_H_
#define _COMMONHELPERS_H_
#include "bondnetwork.h"
#include "povraytools.h"
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
         const string trnsmat,bool cpkview=false);
   /** Adds the bond cylinders to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmbnd is the string to be passed for the transparency of
    * the bonds. */
   static void PutBonds(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmbnd);
   /** Adds special spheres to ofil. ntbs is the number of tabulators
    * that are added at the beginning of each line.
    * trnsmbnd is the string to be passed for the transparency
    * of the spheres.
    * sp is a matrix. each row is a row-vector that has the following
    * structure:
    * sp[i][0], sp[i][1], and sp[i][2] are the coordinates of the i-th sphere centre.
    * sp[i][3] is the radius of the i-th sphere.
    * sp[i][4], sp[i][5], and sp[i][6] are the r, g, and b components of
    * the rgb colour that will be used to colour the ith sphere.  */
   static void PutSpecialSpheres(ofstream &ofil,const int ntbs,\
         const vector<vector<double> > &sp,const string trnsmbnd);
   static void WriteAngleDeclarations(ofstream &ofil,POVRayConfiguration &pvc);
   static void RenderPovfile(const string &povname,bool verbose=false);
   static void RotateCameraAroundLocCam(POVRayConfiguration &pvc, double angle);
   static void RotateCameraAroundUp(POVRayConfiguration &pvc, double angle);
   static void RotateCameraAroundRight(POVRayConfiguration &pvc, double angle);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */
#endif  /* _COMMONHELPERS_H_ */


