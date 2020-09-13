#ifndef _INTERPOLATORS_H_
#define _INTERPOLATORS_H_
#include <vector>
using std::vector;
#include "spline.h"

/* ************************************************************************** */
class Interpolator {
/* ************************************************************************** */
public:
   Interpolator();
/* ************************************************************************** */
   static double Linear(const double x1,const double y1,const double x2,const double y2,const double xi);
   static double Linear(const vector<double> &x,const vector<double> &y,const double xi);
   double Spline(const vector<double> &x,const vector<double> &y,const double xi);
   inline double Spline(const double xi) {return sp(xi);} //Unsafe, but faster (in principle).
   /** Returns true if v[0]<=vi<=v[n-1]. It saves the left (l) and right (r) positions that encloses the value vi in v. If
    * there is a v[k] element that exactly matches vi (i.e. vi==v[k]), then it sets
    * l=k, and r=k+1.  */
   static bool BinarySearch(const vector<double> &v,size_t &ul,size_t &ur,const double vi);
   void SetPoints(const vector<double> &x,const vector<double> &y);
   /** This will prepare the interpolator object to be ready for a new interpolation
    * set. Usefor for reusing the same Interpolator instance with different
    * interpolation vectors.  */
   void Reset() { spIsSetup=false; }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   tk::spline sp;
   bool spIsSetup;
   vector<double> xx,yy;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTERPOLATORS_H_ */

