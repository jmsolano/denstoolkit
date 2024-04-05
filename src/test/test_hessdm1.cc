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
#include "../common/gausswavefunction.h"
#include "../common/mytimer.h"
#include "../common/screenutils.h"

//#define SET_MY_PRECISION (cout << setprecision(10))

void printHessMat(double (&h)[3][3],double (&rh)[3][3],char fd /* first derivative  */,\
      char sd /* second derivative  */);
bool testHessMat(double(&h)[3][3],double (&rh)[3][3]);
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
   refval1=0.0002271651711;
   rho=wf->EvalSpinDensity(0.0e0,0.0e0,0.0e0);
   passed=passed&&(fabs(rho-refval1)<1.0e-13);
   if ( verbose ) {
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      SET_MY_PRECISION;
      cout << "            Value         (Diff)" << '\n';
      cout << "spinrho: " << rho << " (" << (rho-refval1) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   double yy,zz;
   xx=yy=zz=0.0e0;
   double xp=0.1e0,yp=0.1e0,zp=0.2e0;
   SET_MY_PRECISION;
   // ***************************************************
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
   xx=yy=zz=0.0e0;
   xp=0.1e0,yp=0.1e0,zp=0.2e0;
   double refggm1=wf->EvalDensityMatrix1(xx,yy,zz,xp,yp,zp);
   double refggxx[3]={0.0026916196341634,-0.0047278138486487,-0.0074855108384113};
   double refggxp[3]={8.0358325728843e-05,-0.0064731533172912,-0.0081621533763259};
   double ggm1,gxx[3],gxp[3];
   wf->EvalGradDensityMatrix1(xx,yy,zz,xp,yp,zp,ggm1,gxx,gxp);
   passed=passed&&(fabs(ggm1-refggm1)<1.0e-14);
   for ( int i=0 ; i<3 ; ++i ) {
      passed=passed&&(fabs(gxx[i]-refggxx[i])<1.0e-14);
      passed=passed&&(fabs(gxp[i]-refggxp[i])<1.0e-14);
   }
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      //cout << setprecision(14);
      SET_MY_PRECISION;
      cout << "                   Value         (Diff)" << '\n';
      cout << "Gamma1Ref(xx,xp): " << refggm1 << '\n';
      cout << "   Gamma1(xx,xp): " << ggm1 << " (" << (ggm1-refggm1) << ')' << '\n';
      cout << "dxxGamma1(xx,xp): " << gxx[0] << " (" << fabs(gxx[0]-refggxx[0]) << ')' << '\n';
      cout << "dyyGamma1(xx,xp): " << gxx[1] << " (" << fabs(gxx[1]-refggxx[1]) << ')' << '\n';
      cout << "dzzGamma1(xx,xp): " << gxx[2] << " (" << fabs(gxx[2]-refggxx[2]) << ')' << '\n';
      cout << "dxpGamma1(xx,xp): " << gxp[0] << " (" << fabs(gxp[0]-refggxp[0]) << ')' << '\n';
      cout << "dypGamma1(xx,xp): " << gxp[1] << " (" << fabs(gxp[1]-refggxp[1]) << ')' << '\n';
      cout << "dzpGamma1(xx,xp): " << gxp[2] << " (" << fabs(gxp[2]-refggxp[2]) << ')' << '\n';
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalGradDensityMatrix1(xx,yy,zz,xp,yp,zp,ggm1,gxx,gxp);
      xx+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("GradDensityMatrix1"); }
   // ***************************************************
   double vxx[3]={0.0e0,0.1e0,0.3e0};
   double vxp[3]={0.1e0,0.3e0,0.2e0};
   double refhgm1,refhggxx[3],refhggxp[3];
   wf->EvalGradDensityMatrix1(vxx[0],vxx[1],vxx[2],vxp[0],vxp[1],vxp[2],refhgm1,refhggxx,refhggxp);
   double refhx1x1[3][3]={
      {-0.013758614012534,  -0.00027545483005345, -0.0063295105276852 },
      {-0.00027545483005345,-0.009875575647393,   -0.00050849817726172},
      {-0.0063295105276852, -0.00050849817726172,  0.021633349781183  }
   };
   double refhx1x2[3][3]={
      { 0.0056394992752296,-3.1556891983718e-05,-0.0038304309085288 },
      { 0.00016719544307187,0.0094585152212566, -0.00031144472434839},
      {-0.004787524771857, -0.0018910508195085 , 0.028218019818115  }
   };
   double refhx2x2[3][3]={
      {-0.01314103820661,    0.00064500491796486,-0.0074696832085443},
      { 0.00064500491796486,-0.0064380924699305, -0.0017218098583244},
      {-0.0074696832085443, -0.0017218098583244,  0.023233868275934 }
   };
   double hxx[3][3],hxp[3][3],hpp[3][3];
   wf->EvalHessDensityMatrix1(vxx,vxp,ggm1,gxx,gxp,hxx,hxp,hpp);
   passed=passed&&testHessMat(hxx,refhx1x1);
   passed=passed&&testHessMat(hxp,refhx1x2);
   passed=passed&&testHessMat(hpp,refhx2x2);
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
      cout << "                   Value         (Diff)" << '\n';
      cout << "Gamma1Ref(x1,x2): " << refhgm1 << '\n';
      cout << "   Gamma1(x1,x2): " << ggm1 << " (" << (ggm1-refhgm1) << ')' << '\n';
      cout << "dxxGamma1(x1,x2): " << gxx[0] << " (" << fabs(gxx[0]-refhggxx[0]) << ')' << '\n';
      cout << "HessDxiDxjGamma(x1,x2):\n";
      printHessMat(hxx,refhx1x1,'1','1');
      printHessMat(hxp,refhx1x2,'1','2');
      printHessMat(hpp,refhx2x2,'2','2');
      //SET_MY_PRECISION;
      if ( !passed ) { ScreenUtils::SetScrNormalFont(); }
   }
   // ***************************************************
   if ( verbose ) {
      ScreenUtils::PrintScrCharLine('-');
      if ( !passed ) { ScreenUtils::SetScrRedBoldFont(); }
   }
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalHessDensityMatrix1(vxx,vxp,ggm1,gxx,gxp,hxx,hxp,hpp);
      vxx[0]+=dx;
   }
   timer.End();
   if ( verbose ) { timer.PrintElapsedTimeMilliSec("HessDensityMatrix1"); }
   cout << (passed? "PASSED" : "FAILED") << '\n';
   if ( verbose ) { ScreenUtils::SetScrNormalFont(); }
   return EXIT_SUCCESS;
}
void printHessMat(double (&h)[3][3],double (&rh)[3][3],char fd,char sd) {
      vector<string> rowstr,colstr;
      rowstr.push_back("x"); rowstr.push_back("y"); rowstr.push_back("z");
      colstr.push_back("x"); colstr.push_back("y"); colstr.push_back("z");
      for ( size_t i=0 ; i<3 ; ++i ) { rowstr[i]+=fd; }
      for ( size_t i=0 ; i<3 ; ++i ) { colstr[i]+=sd; }
      ScreenUtils::PrintScrCharLine('=');
      for ( size_t i=0 ; i<3 ; ++i ) { cout << std::setw(20) << colstr[i]; }
      cout << setprecision(14) << '\n';
      for ( size_t i=0 ; i<3 ; ++i ) {
         cout << rowstr[i] << "  ";
         for ( size_t j=0 ; j<3 ; ++j ) {
            cout << std::setw(20) << h[i][j] << (j==2? '\n' : ' ');
         }
      }
      cout << "----------------------------------" << '\n';
      SET_MY_PRECISION;
      for ( size_t i=0 ; i<3 ; ++i ) {
         cout << rowstr[i] << "  ";
         for ( size_t j=0 ; j<3 ; ++j ) {
            cout << '(' << std::setw(16) << (h[i][j]-rh[i][j]) << ')' << (j==2? '\n' : ' ');
         }
      }
      ScreenUtils::PrintScrCharLine('=');
}
bool testHessMat(double(&h)[3][3],double (&rh)[3][3]) {
   bool res=true;
   for ( int i=0 ; i<3 ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) {
         res=res&&(fabs(h[i][j]-rh[i][j])<1.0e-14);
      }
   }
   return res;
}


