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

int main (int argc, char *argv[]) {
   bool verbose=false;
   if ( argc==2 ) {
   }
   bool passed=true;
   int N=2048;
   string fname="XeF2.wfx";
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
   double rho=wf->EvalDensity(1.3e0,0.8e0,0.6e0);
   double gamm=wf->EvalDensityMatrix1(1.3e0,0.8e0,0.6e0,1.3e0,0.8e0,0.6e0);
   double pidens=wf->EvalFTDensity(0.4,0.8,0.2);
   double refval1=0.17233627471329;
   double refval2=1.354323488727;
   passed=passed&&(fabs(rho-refval1)<1.0e-12);
   cout << setprecision(14);
   passed=passed&&(fabs(pidens-refval2)<1.0e-12);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "            Value         (Diff)" << '\n';
      cout << " rho: " << rho << " (" << fabs(rho-refval1) << ')' << '\n';
      cout << "g1xx: " << rho << " (" << fabs(gamm-refval1) << ')' << '\n';
      cout << "  pi: " << pidens << " (" << fabs(pidens-refval2) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   if ( verbose ) { ScreenUtils::PrintScrCharLine('-'); }
   MyTimer timer;
   double dx=0.01e0/double(N);
   double xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalDensity(xx,0.4e0,0.2e0);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Norm Density"); }
   SET_MY_PRECISION;
   cout << setprecision(14);
   refval1=1.0695955871204;
   rho=wf->EvalDensity(xx,0.3e0,0.1e0);
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
   // ***************************************************
   double x1=0.61e0, y1=0.77e0, z1=0.85e0;
   double x2=0.43e0, y2=0.63e0, z2=0.92e0;
   SET_MY_PRECISION;
   double rho1=wf->EvalDensity(x1,y1,z1);
   double rho2=wf->EvalDensity(x2,y2,z2);
   gamm=wf->EvalDensityMatrix1(x1,y1,z1,x2,y2,z2);
   double rho2cr=0.5e0*rho1*rho2-0.25e0*gamm*gamm;
   double rho2cv=wf->EvalRho2ClosedShell(x1,y1,z1,x2,y2,z2);
   passed=passed&&(fabs(rho2cv-rho2cr)<1.0e-14);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "                   Value         (Diff)" << '\n';
      cout << "rho2clRef(x,xp): " << rho2cr << '\n';
      cout << "rho2clVal(x,xp): " << rho2cv << " (" << (rho2cv-rho2cr) << ')' << '\n';
   }
   // ***************************************************
   x1=-2.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho2cv=wf->EvalRho2ClosedShell(x1,y1,z1,x2,y2,z2);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Rho2ClosedShell"); }
   // ***************************************************
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   double p1x=0.37,p1y=0.28,p1z=0.42;
   double p2x=0.17,p2y=0.37,p2z=0.13;
   double pi1=wf->EvalFTDensity(p1x,p1y,p1z);
   double pi2=wf->EvalFTDensity(p2x,p2y,p2z);
   complex<double> g1p=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z);
   double pi2csr=0.5e0*pi1*pi2-0.25e0*norm(g1p);
   double pi2csv=wf->EvalFTRho2ClosedShell(p1x,p1y,p1z,p2x,p2y,p2z);
   passed=passed&&(fabs(pi2csv-pi2csr)<1.0e-14);
   if ( verbose ) {
      SET_MY_PRECISION;
      cout << "                   Value         (Diff)" << '\n';
      cout << "pi2clRef(p1,p2): " << pi2csr << '\n';
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "pi2clVal(p1,p2): " << pi2csv << " (" << fabs(pi2csv-pi2csr) << ')' << '\n';
   }
   // ***************************************************
   timer.Start();
   p1x=0.793;
   for ( int i=0 ; i<N ; ++i ) {
      rho2cv=wf->EvalFTRho2ClosedShell(p1x,p1y,p1z,p2x,p2y,p2z);
      p1x+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Pi2ClosedShell"); }
   // ***************************************************
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   cout << (passed? "PASSED" : "FAILED") << '\n';
   if ( verbose ) { ScreenUtils::SetScrNormalFont(); }
   return EXIT_SUCCESS;
}


