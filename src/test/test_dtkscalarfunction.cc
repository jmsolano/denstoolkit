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
#include "../common/dtkscalarfunction3d.h"
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
   SET_MY_PRECISION;
   bool passed=true;
   double refval1=77.2060292975;
   double rho=wf->EvalDensity(0.0e0,0.0e0,0.0e0);
   passed=passed&&(fabs(rho-refval1)<4.0e-11);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "rho(raw): " << rho << " (" << fabs(rho-refval1) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   DTKScalarFunction f(*wf);
   f.SetScalarFunction('d');
   double rhop=f.fptr(0.0e0,0.0e0,0.0e0);
   passed=passed&&(fabs(rhop-refval1)<4.0e-11);
   double rhoc=f.f(0.0e0,0.0e0,0.0e0);
   passed=passed&&(fabs(rhoc-refval1)<4.0e-11);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "rho(ptr): " << rhop << " (" << fabs(rhop-refval1) << ')' << '\n';
      cout << "rho(cas): " << rhoc << " (" << fabs(rhoc-refval1) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   f.SetScalarFunction('d');
   size_t NN=16;
#if PARALLELISEDTK
   NN=4;
#endif
   vector<double> tmp(NN*NN*NN);
   int iter=1000;
   MyTimer timer;
   timer.Start();
   size_t count=0;
   double x=0.0e0,y=0.0e0,z=0.0e0;
   for ( int l=0 ; l<iter ; ++l ) {
      x=0.0e0;
      count=0;
      for ( size_t i=0 ; i<NN ; ++i ) {
         y=0.0e0;
         for ( size_t j=0 ; j<NN ; ++j ) {
            z=0.0e0;
            for ( size_t k=0 ; k<NN ; ++k ) { tmp[count++]=f.fptr(x,y,z); z+=0.1e0; }
            y+=0.11e0;
         }
         x+=0.052e0;
      }
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Pointer to function"); }
   //
   timer.Start();
   x=0.0e0; y=0.0e0; z=0.0e0;
   for ( int l=0 ; l<iter ; ++l ) {
      count=0;
      x=0.0e0;
      for ( size_t i=0 ; i<NN ; ++i ) {
         y=0.0e0;
         for ( size_t j=0 ; j<NN ; ++j ) {
            z=0.0e0;
            for ( size_t k=0 ; k<NN ; ++k ) { tmp[count++]=f.f(x,y,z); z+=0.1e0; }
            y+=0.11e0;
         }
         x+=0.052e0;
      }
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Function switch"); }

   //
   timer.Start();
   x=0.0e0; y=0.0e0; z=0.0e0;
   for ( int l=0 ; l<iter ; ++l ) {
      count=0;
      x=0.0e0;
      for ( size_t i=0 ; i<NN ; ++i ) {
         y=0.0e0;
         for ( size_t j=0 ; j<NN ; ++j ) {
            z=0.0e0;
            //for ( size_t k=0 ; k<10 ; ++k ) { tmp[count++]=wf->EvalELF(x,y,z); z+=0.1e0; }
            for ( size_t k=0 ; k<NN ; ++k ) { tmp[count++]=wf->EvalDensity(x,y,z); z+=0.1e0; }
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


