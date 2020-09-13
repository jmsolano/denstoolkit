#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <algorithm>
#include "interpolators.h"
#include "screenutils.h"

#ifndef DEBUG
#define DEBUG 0
#endif

Interpolator::Interpolator() : sp() {
   spIsSetup=false;
}
double Interpolator::Linear(const double x1,const double y1,\
      const double x2,const double y2,const double xi) {
#if DEBUG
   if ( x2==x1 ) {
      ScreenUtils::DisplayErrorMessage("x2==x1");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
#endif
   double m=(y2-y1)/(x2-x1);
   return (y1+m*(xi-x1));
}
double Interpolator::Linear(const vector<double> &x,const vector<double> &y,\
      const double xi) {
   /*
   auto x1=std::lower_bound(x.begin(), x.end(), xi);
   //cout << "x[pos]: " <<  *x1 << endl; // gives the actual value
   if ( *x1==xi ) { return y[pos]; }
   size_t pos=x1-x.begin();
   //cout << "pos: " << pos << ", x: " << x[pos] << endl;
   return Linear(x[pos-1],y[pos-1],x[pos],y[pos],xi);
   // */
   size_t l,r;
   bool res=BinarySearch(x,l,r,xi);
   if ( !res ) {
      ScreenUtils::DisplayErrorMessage("Not an element!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return y[res-1];
   }
   if ( x[l]==xi ) { return y[l]; }
   //cout << x[l] << "," << xi << "," << x[r] << endl;
   return Linear(x[l],y[l],x[r],y[r],xi);
}
bool Interpolator::BinarySearch(const vector<double> &v,size_t &ul,size_t &ur,const double val) {
   if ( val<v[0] || val>v.back() ) {
      ScreenUtils::DisplayErrorMessage("Out of range!");
      ul=0;
      ur=v.size()-1;
      return false;
   }
   size_t l=0,r=v.size()-1;;
   size_t m=l+(r-l)/2; 
   while ((r-l)>1) { 
      // Check if x is present at mid 
      if ( v[m]==val ) { ul=m; ur=m+1; return true; }
      // If x greater, ignore left half 
      if ( val>v[m] )  { l=m; }
      // If x is smaller, ignore right half 
      else { r=m; }
      m=l+(r-l)/2;
   } 
   ul=l; ur=r;
   //cout << "[l,r]: " << l << "," << r << endl;
   return true;
}
void Interpolator::SetPoints(const vector<double> &x,const vector<double> &y) {
   //vector<double> tx(x.begin(),x.end());
   //vector<double> ty(y.begin(),y.end());
   sp.set_points(x,y);
   spIsSetup=true;
}
double Interpolator::Spline(const vector<double> &x,const vector<double> &y,const double xi) {
   if ( spIsSetup ) { return sp(xi); }
   SetPoints(x,y);
   return sp(xi);
}


