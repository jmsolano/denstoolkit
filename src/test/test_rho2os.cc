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
   wf->CalcCabAAndCabB();
   if ( !wf->ihaveCABSingleSpin ) { 
      passed=false;
      ScreenUtils::SetScrRedBoldFont();
      cout << (passed? "PASSED" : "FAILED") << '\n';
      ScreenUtils::SetScrNormalFont();
      return EXIT_FAILURE;
   }
   SET_MY_PRECISION;
   double rho=wf->EvalDensity(1.17e0,0.93e0,0.59e0);
   double gamm=wf->EvalDensityMatrix1(1.17e0,0.93e0,0.59e0,1.17e0,0.93e0,0.59e0);
   double pidens=wf->EvalFTDensity(1.41,0.77,1.23);
   double refval1=0.012594741367475;
   double refval2=0.59069810073782;
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
   refval1=2.1742838791774;
   rho=wf->EvalFTDensity(xx,0.37e0,1.1e0);
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
      rho=wf->EvalFTDensity(xx,0.7e0,1.4e0);
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
   double gamm1=wf->EvalDensityMatrix1(x1,y1,z1,x1,y1,z1);
   double gamm2=wf->EvalDensityMatrix1(x2,y2,z2,x2,y2,z2);
   passed=passed&&(fabs(gamm1-rho1)<1.0e-14);
   passed=passed&&(fabs(gamm2-rho2)<1.0e-14);
   double rhoa1=wf->EvalDensityMatrix1Alpha(x1,y1,z1,x1,y1,z1);
   double rhoa2=wf->EvalDensityMatrix1Alpha(x2,y2,z2,x2,y2,z2);
   double rhob1=wf->EvalDensityMatrix1Beta(x1,y1,z1,x1,y1,z1);
   double rhob2=wf->EvalDensityMatrix1Beta(x2,y2,z2,x2,y2,z2);
   double gammaA=wf->EvalDensityMatrix1Alpha(x1,y1,z1,x2,y2,z2);
   double gammaB=wf->EvalDensityMatrix1Beta(x1,y1,z1,x2,y2,z2);
   double rhoaa=0.5e0*(rhoa1*rhoa2-gammaA*gammaA);
   double rhobb=0.5e0*(rhob1*rhob2-gammaB*gammaB);
   double rhoab=0.5e0*rhoa1*rhob2;
   double rhoba=0.5e0*rhoa2*rhob1;
   double rho2or=rhoaa+rhobb+rhoab+rhoba;
   double rho2ov=wf->EvalRho2OpenShell(x1,y1,z1,x2,y2,z2);
   passed=passed&&(fabs(rho2ov-rho2or)<1.0e-14);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "                   Value         (Diff)" << '\n';
      cout << "  gammxx1vsrho1: " << gamm1 << " (" << fabs(gamm1-rho1) << ')' << '\n';
      cout << "  gammxx2vsrho2: " << gamm2 << " (" << fabs(gamm2-rho2) << ')' << '\n';
      cout << "rho2opRef(x,xp): " << rho2or << '\n';
      cout << "rho2opVal(x,xp): " << rho2ov << " (" << (rho2ov-rho2or) << ')' << '\n';
   }
   // ***************************************************
   x1=-2.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho2ov=wf->EvalRho2OpenShell(x1,y1,z1,x2,y2,z2);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Rho2OpenShell"); }
   x1=-2.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rhoa1=wf->EvalDensityMatrix1Alpha(x1,y1,z1,x1,y1,z1);
      rhoa2=wf->EvalDensityMatrix1Alpha(x2,y2,z2,x2,y2,z2);
      rhob1=wf->EvalDensityMatrix1Beta(x1,y1,z1,x1,y1,z1);
      rhob2=wf->EvalDensityMatrix1Beta(x2,y2,z2,x2,y2,z2);
      gammaA=wf->EvalDensityMatrix1Alpha(x1,y1,z1,x2,y2,z2);
      gammaB=wf->EvalDensityMatrix1Beta(x1,y1,z1,x2,y2,z2);
      rhoaa=0.5e0*(rhoa1*rhoa2-gammaA*gammaA);
      rhobb=0.5e0*(rhob1*rhob2-gammaB*gammaB);
      rhoab=0.5e0*rhoa1*rhob2;
      rhoba=0.5e0*rhoa2*rhob1;
      rho2or=rhoaa+rhobb+rhoab+rhoba;
      x1+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Rho2OpenBrute"); }
   // ***************************************************
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   double p1x=1.77,p1y=2.68,p1z=1.72;
   double p2x=1.15,p2y=1.17,p2z=1.13;
   double pi1=wf->EvalFTDensity(p1x,p1y,p1z);
   double pi2=wf->EvalFTDensity(p2x,p2y,p2z);
   complex<double> g1pp1=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p1x,p1y,p1z);
   complex<double> g1pp2=wf->EvalFTDensityMatrix1(p2x,p2y,p2z,p2x,p2y,p2z);
   passed=passed&&(fabs(abs(g1pp1)-pi1)<1.0e-14);
   passed=passed&&(fabs(abs(g1pp2)-pi2)<1.0e-14);
   double pia1=abs(wf->EvalFTDensityMatrix1Alpha(p1x,p1y,p1z,p1x,p1y,p1z));
   double pia2=abs(wf->EvalFTDensityMatrix1Alpha(p2x,p2y,p2z,p2x,p2y,p2z));
   double pib1=abs(wf->EvalFTDensityMatrix1Beta(p1x,p1y,p1z,p1x,p1y,p1z));
   double pib2=abs(wf->EvalFTDensityMatrix1Beta(p2x,p2y,p2z,p2x,p2y,p2z));
   complex<double> tgamma=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z);
   complex<double> tgammaA=wf->EvalFTDensityMatrix1Alpha(p1x,p1y,p1z,p2x,p2y,p2z);
   complex<double> tgammaB=wf->EvalFTDensityMatrix1Beta(p1x,p1y,p1z,p2x,p2y,p2z);
   double piaa=0.5e0*(pia1*pia2-norm(tgammaA));
   double pibb=0.5e0*(pib1*pib2-norm(tgammaB));
   double piab=0.5e0*pia1*pib2;
   double piba=0.5e0*pia2*pib1;
   double pi2osr=piaa+pibb+piab+piba;
   double pi2osv=wf->EvalFTRho2OpenShell(p1x,p1y,p1z,p2x,p2y,p2z);
   passed=passed&&(fabs(pi2osv-pi2osr)<1.0e-14);
   if ( verbose ) {
      SET_MY_PRECISION;
      cout << "                   Value         (Diff)" << '\n';
      cout << "         piA(p1): " << pia1 << '\n';
      cout << "         piA(p2): " << pia2 << '\n';
      cout << "         piB(p1): " << pib1 << '\n';
      cout << "         piB(p2): " << pib2 << '\n';
      cout << "  tgammaA(p1,p2): " << tgammaA << '\n';
      cout << "  tgammaB(p1,p2): " << tgammaB << '\n';
      cout << "sumtgamma(p1,p2): " << (tgammaA+tgammaB) << '\n';
      cout << "   tgamma(p1,p2): " << (tgamma) << " (" << (tgammaA+tgammaB-tgamma) << ')' << '\n';
      cout << "    gammpp1vspi1: " << abs(g1pp1) << " (" << fabs(abs(g1pp1)-pi1) << ')' << '\n';
      cout << "    gammpp2vspi2: " << abs(g1pp2) << " (" << fabs(abs(g1pp2)-pi2) << ')' << '\n';
      cout << " pi2osRef(p1,p2): " << pi2osr << '\n';
      cout << " pi2osVal(p1,p2): " << pi2osv << " (" << fabs(pi2osv-pi2osr) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   // ***************************************************
   timer.Start();
   p1x=-1.731;
   for ( int i=0 ; i<N ; ++i ) {
      rho2ov=wf->EvalFTRho2OpenShell(p1x,p1y,p1z,p2x,p2y,p2z);
      p1x+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Pi2OpenShell"); }
   // ***************************************************
   timer.Start();
   p1x=-1.731;
   for ( int i=0 ; i<N ; ++i ) {
      pia1=abs(wf->EvalFTDensityMatrix1Alpha(p1x,p1y,p1z,p1x,p1y,p1z));
      pia2=abs(wf->EvalFTDensityMatrix1Alpha(p2x,p2y,p2z,p2x,p2y,p2z));
      pib1=abs(wf->EvalFTDensityMatrix1Beta(p1x,p1y,p1z,p1x,p1y,p1z));
      pib2=abs(wf->EvalFTDensityMatrix1Beta(p2x,p2y,p2z,p2x,p2y,p2z));
      tgamma=wf->EvalFTDensityMatrix1(p1x,p1y,p1z,p2x,p2y,p2z);
      tgammaA=wf->EvalFTDensityMatrix1Alpha(p1x,p1y,p1z,p2x,p2y,p2z);
      tgammaB=wf->EvalFTDensityMatrix1Beta(p1x,p1y,p1z,p2x,p2y,p2z);
      piaa=0.5e0*(pia1*pia2-norm(tgammaA));
      pibb=0.5e0*(pib1*pib2-norm(tgammaB));
      piab=0.5e0*pia1*pib2;
      piba=0.5e0*pia2*pib1;
      pi2osr=piaa+pibb+piab+piba;
      p1x+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("Pi2OpenBrute"); }
   // ***************************************************
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   cout << (passed? "PASSED" : "FAILED") << '\n';
   if ( verbose ) { ScreenUtils::SetScrNormalFont(); }
   return EXIT_SUCCESS;
}


