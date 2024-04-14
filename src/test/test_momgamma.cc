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
#include <complex>
#include "../common/gausswavefunction.h"
#include "../common/mytimer.h"
#include "../common/screenutils.h"

//#define SET_MY_PRECISION (cout << setprecision(10))

int main (int argc, char *argv[]) {
   bool verbose=false;
   if ( argc==2 ) {
   }
   bool passed=true;
   int N=512;
   string fname="42.wfx";
   if ( argc>1 ) {
      if ( string(argv[1])==string("-h") ) {
         cout << "Usage:\n" << "   " << argv[0] << " [-v,-h] [-n Npts]" << '\n';
         return EXIT_SUCCESS;
      }
      if ( string(argv[1])==string("-v") ) {
         verbose=true;
      }
      if ( argc==3 ) {
         cout << "Usage:\n" << "   " << argv[0] << " [-v,-h] [-n Npts]" << '\n';
         return EXIT_SUCCESS;
      }
      if ( argc==4 && string(argv[2])==string("-N") ) {
         N=std::stoi(string(argv[3]));
         cout << "Using: " << N << '\n';
      }
   }
   shared_ptr<GaussWaveFunction> wf=shared_ptr<GaussWaveFunction>(new GaussWaveFunction());
   wf->ReadFromFile(fname);
   SET_MY_PRECISION;
   double rho=wf->EvalDensity(0.0e0,0.0e0,0.0e0);
   double pidens=wf->EvalFTDensity(0.0,0.0,0.2);
   double refval1=0.028511686565958;
   double refval2=10.790258029591;
   passed=passed&&(fabs(rho-refval1)<1.0e-12);
   cout << setprecision(14);
   passed=passed&&(fabs(pidens-refval2)<1.0e-12);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "            Value         (Diff)" << '\n';
      cout << "rho: " << rho << " (" << fabs(rho-refval1) << ')' << '\n';
      cout << " pi: " << pidens << " (" << fabs(pidens-refval2) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   wf->CalcCabAAndCabB();
   if ( !wf->ihaveCABSingleSpin ) { 
      passed=false;
      ScreenUtils::SetScrRedBoldFont();
      cout << (passed? "PASSED" : "FAILED") << '\n';
      ScreenUtils::SetScrNormalFont();
      return EXIT_FAILURE;
   }
   if ( verbose ) { ScreenUtils::PrintScrCharLine('-'); }
   MyTimer timer;
   double dx=0.01e0/double(N);
   double xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalDensity(xx,0.0e0,0.0e0);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Norm Density"); }
   rho=wf->EvalDensity(xx,0.0e0,0.0e0);
   SET_MY_PRECISION;
   refval1=0.01756082669;
   rho=wf->EvalDensity(xx,0.0e0,0.0e0);
   passed=passed&&(fabs(rho-refval1)<1.0e-11);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "            Value         (Diff)" << '\n';
      cout << "rho: " << rho << " (" << fabs(rho-refval1) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); cout << std::flush; }
   }
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalFTDensity(xx,0.0e0,0.0e0);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Momentum space density"); }
   //
   double p1x=0.35e0,p1y=0.44e0,p1z=0.22e0;
   double p2x=0.42e0,p2y=0.25e0,p2z=0.36e0;
   complex<double> g1p=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p1x,p1y,p1z);
   double rhopref=wf->EvalFTDensity(p1x,p1y,p1z);
   passed=passed&&(fabs(g1p.real()-rhopref)<1.0e-14);
   passed=passed&&(fabs(g1p.imag())<1.0e-12);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      SET_MY_PRECISION;
      cout << "            Value         (Diff)" << '\n';
      cout << "rhopref1: " << rhopref << '\n';
      cout << "rhop.re: " << g1p.real() << " (" << fabs(g1p.real()-rhopref) << ')' << '\n';
      cout << "rhop.im: " << fabs(g1p.imag()) << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   //
   g1p=wf->EvalFTDensityMatrix1(p2x,p2y,p2z,p2x,p2y,p2z);
   rhopref=wf->EvalFTDensity(p2x,p2y,p2z);
   passed=passed&&(fabs(g1p.real()-rhopref)<1.0e-14);
   passed=passed&&(fabs(g1p.imag())<1.0e-12);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "rhopref2: " << rhopref << '\n';
      cout << "rhop.re: " << g1p.real() << " (" << fabs(g1p.real()-rhopref) << ')' << '\n';
      cout << "rhop.im: " << g1p.imag() << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   complex<double> g1pref(6.6274552808751,0.47050515036446);
   g1p=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z);
   passed=passed&&(fabs(g1p.real()-g1pref.real())<1.0e-12);
   passed=passed&&(fabs(g1p.imag()-g1pref.imag())<1.0e-12);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "g1p.re: " << g1p.real() << " (" << fabs(g1p.real()-g1pref.real()) << ')' << '\n';
      cout << "g1p.im: " << g1p.imag() << " (" << fabs(g1p.imag()-g1pref.imag()) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   //
   p1x=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      g1p=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z);
      p1x+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Momentum space density matrix 1"); }
   cout << (passed? "PASSED" : "FAILED") << '\n';
   if ( verbose ) { ScreenUtils::SetScrNormalFont(); }
   return EXIT_SUCCESS;
}


