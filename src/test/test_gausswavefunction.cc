#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::scientific;
using std::setprecision;
#include <memory>
using std::shared_ptr;
#include <string>
using std::string;
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
   double refval1=0.02851168657;
   double refval2=10.79025803;
   passed=passed&&(fabs(rho-refval1)<1.0e-11);
   passed=passed&&(fabs(pidens-refval2)<5.0e-10);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "            Value         (Diff)" << '\n';
      cout << "rho: " << rho << " (" << fabs(rho-refval1) << ')' << endl;
      cout << " pi: " << pidens << " (" << fabs(pidens-refval2) << ')' << endl;
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   wf->CalcCabAAndCabB();
   if ( !wf->ihaveCABSingleSpin ) { 
      passed=false;
      ScreenUtils::SetScrRedBoldFont();
      cout << (passed? "PASSED" : "FAILED") << endl;
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
      cout << "rho: " << rho << " (" << fabs(rho-refval1) << ')' << endl;
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalFTDensity(xx,0.0e0,0.0e0);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Momentum space density"); }
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalSpinDensity(xx,0.0e0,0.0e0);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Spin Density"); }
   refval1=0.0002271651711;
   rho=wf->EvalSpinDensity(0.0e0,0.0e0,0.0e0);
   passed=passed&&(fabs(rho-refval1)<1.0e-13);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      SET_MY_PRECISION;
      cout << "            Value         (Diff)" << '\n';
      cout << "spinrho: " << rho << " (" << (rho-refval1) << ')' << endl;
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   double yy,zz;
   xx=yy=zz=0.0e0;
   double xp=0.1e0,yp=0.1e0,zp=0.2e0;
   SET_MY_PRECISION;
   double reggamma=wf->EvalDensityMatrix1(xx,yy,zz,xp,yp,zp);
   double grlgamma=wf->EvalGeneralDensityMatrix1(xx,yy,zz,xp,yp,zp,false,wf->cab);
   double g1a=wf->EvalGeneralDensityMatrix1(xx,yy,zz,xp,yp,zp,true,wf->cabA);
   double refg1a=0.013130273669108;
   double g1b=wf->EvalGeneralDensityMatrix1(xx,yy,zz,xp,yp,zp,true,wf->cabB);
   double refg1b=0.012829822024859;
   passed=passed&&(fabs(g1a-refg1a)<1.0e-14);
   passed=passed&&(fabs(g1b-refg1b)<1.0e-14);
   passed=passed&&(fabs(g1a+g1b-reggamma)<1.0e-14);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << setprecision(14);
      cout << "            Value         (Diff)" << '\n';
      cout << "    Gamma1: " << reggamma << '\n';
      cout << "GralGamma1: " << grlgamma << " (" << (grlgamma-reggamma) << ')' << '\n';
      cout << "Sum-Gamma1: " << (g1a+g1b) << " (" << (g1a+g1b-reggamma) << ')' << '\n';
      cout << "AlphGamma1: " << g1a << " (" << (refg1a-g1a) << ')' << '\n';
      cout << "BetaGamma1: " << g1b << " (" << (refg1b-g1b) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   // ***************************************************
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalDensityMatrix1(xx,yy,zz,xp,yp,zp);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Regular Gamma1"); }
   // ***************************************************
   xx=0.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalGeneralDensityMatrix1(xx,yy,zz,xp,yp,zp,false,wf->cab);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("General Gamma1"); }
   // ***************************************************
   xx=0.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalDensityMatrix1Alpha(xx,yy,zz,xp,yp,zp);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Gamma1 Alpha"); }
   // ***************************************************
   xx=0.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalDensityMatrix1Beta(xx,yy,zz,xp,yp,zp);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Gamma1 Beta"); }
   // ***************************************************
   xx=yy=zz=0.0e0;
   xp=0.1e0,yp=0.1e0,zp=0.2e0;
   double rho1=wf->EvalDensity(xx,yy,zz);
   double rho2=wf->EvalDensity(xp,yp,zp);
   double gamm=wf->EvalDensityMatrix1(xx,yy,zz,xp,yp,zp);
   double gama=wf->EvalDensityMatrix1Alpha(xx,yy,zz,xp,yp,zp);
   double gamb=wf->EvalDensityMatrix1Beta(xx,yy,zz,xp,yp,zp);
   SET_MY_PRECISION;
   double rho2or=0.5e0*(rho1*rho2-gama*gama-gamb*gamb);
   double rho2ov=wf->EvalRho2OpenShell(xx,yy,zz,xp,yp,zp);
   double rho2cr=0.5e0*rho1*rho2-0.25e0*gamm*gamm;
   double rho2cv=wf->EvalRho2ClosedShell(xx,yy,zz,xp,yp,zp);
   passed=passed&&(fabs(rho2ov-rho2or)<1.0e-14);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      cout << "                   Value         (Diff)" << '\n';
      cout << "rho2opRef(x,xp): " << rho2or << '\n';
      cout << "rho2opVal(x,xp): " << rho2ov << " (" << (rho2ov-rho2or) << ')' << '\n';
      cout << "rho2clRef(x,xp): " << rho2cr << '\n';
      cout << "rho2clVal(x,xp): " << rho2cv << " (" << (rho2cv-rho2cr) << ')' << '\n';
   }
   // ***************************************************
   xx=0.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalRho2ClosedShell(xx,yy,zz,xp,yp,zp);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Rho2ClosedShell"); }
   // ***************************************************
   xx=0.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalRho2OpenShell(xx,yy,zz,xp,yp,zp);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Rho2OpenShell"); }
   // ***************************************************
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   cout << (passed? "PASSED" : "FAILED") << '\n';
   if ( verbose ) { ScreenUtils::SetScrNormalFont(); }
   return EXIT_SUCCESS;
}


