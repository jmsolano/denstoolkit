#include <cstdlib>
#include <iostream>
using std::cout;
#include <iomanip>
using std::scientific;
using std::setprecision;
#include <memory>
using std::shared_ptr;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "../common/mytimer.h"
#include "../common/gausswavefunction.h"
#include "../common/fldtypesdef.h"
#include "../common/dtkscalarfunction6d.h"
#include "../common/screenutils.h"

int main (int argc, char *argv[]) {
   string fname="ch4.wfx";
   bool verbose=false;
   if ( argc==2 ) {
      if ( string(argv[1])==string("-h") ) {
         cout << "Usage:\n   " << argv[0] << " [-v,-h]" << '\n';
         return EXIT_SUCCESS;
      } else if ( string(argv[1])==string("-v") ) {
        verbose=true;
      } else {
         cout << "Unknown option..." << '\n';
         return EXIT_FAILURE;
      }
   }
   shared_ptr<GaussWaveFunction> wf=shared_ptr<GaussWaveFunction>(new GaussWaveFunction());
   wf->ReadFromFile(fname);
   bool passed=true;
   double gamma=wf->EvalDensityMatrix1(0.0e0,0.0e0,0.0e0,0.1,0.1,0.1);
   double refval1=35.258420520776;
   passed=passed&&(fabs(gamma-refval1)<1.0e-12);
   SET_MY_PRECISION;
   if ( verbose ) {
      cout << "gamma(raw): " << gamma << " (" << fabs(gamma-refval1) << ')' << '\n';
   }
   DTKScalarFunction6D f(*wf);
   f.SetScalarFunction('g');
   double gammac=f.f(0.0e0,0.0e0,0.0e0,0.1,0.1,0.1);
   passed=passed&&(fabs(gammac-refval1)<1.0e-12);
   if ( verbose ) {
      cout << "gamma(cas): " << gammac << " (" << fabs(gammac-refval1) << ')' << '\n';
   }

   f.SetScalarFunction('g');
   int iter=1000;
   int NN=16;
#if PARALLELISEDTK
   NN=4;
#endif
   vector<double> tmp(NN*NN*NN,0.0e0);
   MyTimer timer;
   timer.Start();
   int count=0;
   double x=0.0e0,y=0.0e0,z=0.0e0;
   double xp=0.1e0,yp=0.1e0,zp=0.1e0;
   for ( int l=0 ; l<iter ; ++l ) {
      x=0.0e0;
      count=0;
      for ( int i=0 ; i<NN ; ++i ) {
         y=0.0e0;
         for ( int j=0 ; j<NN ; ++j ) {
            z=0.0e0;
            for ( int k=0 ; k<NN ; ++k ) { tmp[count++]=f.f(x,y,z,xp,yp,zp); z+=0.1e0; }
            y+=0.11e0;
         }
         x+=0.052e0;
      }
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Function switch"); }

   //
   timer.Start();
   count=0;
   x=0.0e0; y=0.0e0; z=0.0e0;
   for ( int l=0 ; l<iter ; ++l ) {
      x=0.0e0;
      count=0;
      for ( int i=0 ; i<NN ; ++i ) {
         y=0.0e0;
         for ( int j=0 ; j<NN ; ++j ) {
            z=0.0e0;
            for ( int k=0 ; k<NN ; ++k ) { tmp[count++]=wf->EvalDensityMatrix1(x,y,z,xp,yp,zp); z+=0.1e0; }
            y+=0.11e0;
         }
         x+=0.052e0;
      }
   }
   timer.End();
   if ( verbose ) {
      timer.PrintElapsedTimeMilliSec("Raw function call");
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   cout << (passed? "PASSED" : "FAILED") << '\n';
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   return EXIT_SUCCESS;
}


