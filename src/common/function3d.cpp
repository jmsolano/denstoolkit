#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cmath>
#include "function3d.h"

Function3D::Function3D() {
   x1.resize(3);
   x2.resize(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      x1[i]=-1.75e0;
      x2[i]=1.24e0;
   }
}
double Function3D::testF1(const vector<double> &x) {
   double res1=0.0e0,res2=0.0e0;;
   for ( size_t i=0 ; i<x.size() ; ++i ) {
      res1+=((x[i]-x1[i])*(x[i]-x1[i]));
      res2+=((x[i]-x2[i])*(x[i]-x2[i]));
   }
   return (exp(-res1)+exp(-res2));
}
double Function3D::analytF1Res() {
   return 1.74685581658e0; //int _{3/2}^{3/2}testF1
}


