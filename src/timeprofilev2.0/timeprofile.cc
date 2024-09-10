#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <memory>
using std::shared_ptr;
#include <iomanip>
using std::scientific;
using std::setprecision;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;
#include <complex>
#include "optflags.h"
#include "screenutils.h"
#include "mytimer.h"
#include "gausswavefunction.h"

int main (int argc, char *argv[]) {
   /* ************************************************************************** */
   MyTimer globtimer;
   globtimer.Start();
   /* ************************************************************************** */
   /* Configures the program (options, variables, etc.)  */
   shared_ptr<OptionFlags> options = shared_ptr<OptionFlags>(new OptionFlags(argc,argv));
   switch ( options->GetExitCode() ) {
      case OptionFlagsBase::ExitCode::OFEC_EXITNOERR :
         std::exit(EXIT_SUCCESS);
         break;
      case OptionFlagsBase::ExitCode::OFEC_EXITERR :
         std::exit(EXIT_FAILURE);
         break;
      default :
         break;
   }
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   /* Main corpus  */
   bool verbose=true;
   int N=1024;
   string fname=argv[1];
   shared_ptr<GaussWaveFunction> wf=shared_ptr<GaussWaveFunction>(new GaussWaveFunction());
   wf->ReadFromFile(fname);
   if ( options->stpspindens && wf->ihaveSingleSpinOrbs ) {
      cout << "Setting up single-spin density matrices." << '\n';
      wf->CalcCabAAndCabB();
   }
   
   vector<string> property(0);
   vector<string> epslabels(0);
   vector<double> computeTime(0);

   // ---------------------------------------------------
   MyTimer timer;
   double rho;
   double dx=0.01e0/double(N);
   double xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalDensity(xx,0.4e0,0.2e0);
      xx+=dx;
   }
   timer.End();
   string fld2eval="Rho";
   string epsstr="'{/Symbol r}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   double reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalSpinDensity(xx,0.4e0,0.2e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="Rhos";
   epsstr="'{/Symbol r}^s'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   /*
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalLOL(xx,0.4e0,0.2e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="LOL";
   epsstr="'LOL'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalMolElecPot(xx,0.4e0,0.2e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="MEP";
   epsstr="'MEP'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // */
   // ---------------------------------------------------
   double gamma;
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      gamma=wf->EvalDensityMatrix1(xx,0.4e0,0.2e0,xx,0.91e0,0.27e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="Gamma";
   epsstr="'{/Symbol G}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   double gx1[3],gx2[3];
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalGradDensityMatrix1(xx,0.4e0,0.2e0,xx,0.91e0,0.27e0,gamma,gx1,gx2);
      xx+=dx;
   }
   timer.End();
   fld2eval="Grad6Gamma";
   epsstr="'{/Symbol \\321}_6{/Symbol G}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   gx1[0]=xx; gx1[1]=0.37e0; gx1[2]=0.23e0;
   gx2[0]=xx; gx2[1]=0.27e0; gx2[2]=0.33e0;
   for ( int i=0 ; i<N ; ++i ) {
      wf->EvalLapDensityMatrix1(gx1,gx2);
      xx+=dx;
   }
   timer.End();
   fld2eval="Lap6Gamma";
   epsstr="'{/Symbol \\321}_6^2{/Symbol G}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      gamma=wf->EvalDensityMatrix1Alpha(xx,0.4e0,0.2e0,xx,0.91e0,0.27e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="GammaAlpha";
   epsstr="'{/Symbol G}^{/Symbol a}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   double rho2;
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      gamma=wf->EvalDensityMatrix1(xx,0.4e0,0.2e0,0.91e0,xx,0.27e0);
      rho=wf->EvalDensity(xx,0.4e0,0.2e0);
      rho*=(wf->EvalDensity(0.91e0,xx,0.27e0));
      rho2=0.5e0*rho-0.25e0*gamma*gamma;
      xx+=dx;
   }
   timer.End();
   fld2eval="Rho2ClosedShellRaw";
   epsstr="'{/Symbol r}_{2,CSraw}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalRho2ClosedShell(xx,0.4e0,0.2e0,xx,0.91e0,0.27e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="Rho2ClosedShell";
   epsstr="'{/Symbol r}_{2,CS}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   double x1[3],x2[3];
   double rhoa1,rhoa2,rhob1,rhob2;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rhoa1=wf->EvalDensityMatrix1Alpha(x1[0],x1[1],x1[2],x1[0],x1[1],x1[2]);
      rhoa2=wf->EvalDensityMatrix1Alpha(x2[0],x2[1],x2[2],x2[0],x2[1],x2[2]);
      rhob1=wf->EvalDensityMatrix1Beta(x1[0],x1[1],x1[2],x1[0],x1[1],x1[2]);
      rhob2=wf->EvalDensityMatrix1Beta(x2[0],x2[1],x2[2],x2[0],x2[1],x2[2]);
      gamma=wf->EvalDensityMatrix1Alpha(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      rho2=-(wf->EvalDensityMatrix1Beta(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]));
      rho2-=gamma;
      rho2+=(rhoa1*rhoa2);
      rho2+=(rhob1*rhob2);
      rho2+=(rhoa1*rhob2);
      rho2+=(rhob1*rhoa2);
      rho*=0.5e0;
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="Rho2OpenShellRaw";
   epsstr="'{/Symbol r}_{2,OSraw}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalRho2OpenShell(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="Rho2OpenShell";
   epsstr="'{/Symbol r}_{2,OS}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      rho=wf->EvalFTDensity(xx,0.4e0,0.2e0);
      xx+=dx;
   }
   timer.End();
   fld2eval="Pi";
   epsstr="'~{/Symbol r}{0.1\\~}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   complex<double> ppii;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      ppii=wf->EvalFTDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="GammaMomSp";
   epsstr="'~{/Symbol G}{0.3\\~}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      ppii=wf->EvalFTDensityMatrix1Alpha(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="GammaAlphaMomSp";
   epsstr="'~{/Symbol G}{0.3\\~}^{/Symbol a}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   epslabels.push_back(epsstr);
   property.push_back(fld2eval);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   double pi2;
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      pi2=wf->EvalFTRho2ClosedShell(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="Pi2ClsdShMomSp";
   epsstr="'~{/Symbol r}{0.1\\~}_{2,CS}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      pi2=0.5e0*(wf->EvalFTDensity(x1[0],x1[1],x1[2]));
      pi2*=(wf->EvalFTDensity(x2[0],x2[1],x2[2]));
      ppii=wf->EvalFTDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      pi2+=(0.5e0*std::norm(ppii));
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="Pi2ClsdShMomSpRaw";
   epsstr="'~{/Symbol r}{0.1\\~}_{2,CSraw}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   for ( int i=0 ; i<N ; ++i ) {
      pi2=wf->EvalFTRho2OpenShell(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="Pi2OpenShMomSp";
   epsstr="'~{/Symbol r}{0.1\\~}_{2,OS}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------
   dx=0.01e0/double(N);
   xx=-1.0e0;
   x1[0]=xx; x1[1]=0.4e0; x1[2]=0.2e0;
   x2[0]=0.91e0; x2[1]=xx; x2[2]=0.27e0;
   timer.Start();
   double pia1,pia2,pib1,pib2;
   for ( int i=0 ; i<N ; ++i ) {
      pia1=std::abs(wf->EvalFTDensityMatrix1Alpha(x1[0],x1[1],x1[2],x1[0],x1[1],x1[2]));
      pia2=std::abs(wf->EvalFTDensityMatrix1Alpha(x2[0],x2[1],x2[2],x2[0],x2[1],x2[2]));
      pib1=std::abs(wf->EvalFTDensityMatrix1Beta(x1[0],x1[1],x1[2],x1[0],x1[1],x1[2]));
      pib2=std::abs(wf->EvalFTDensityMatrix1Beta(x2[0],x2[1],x2[2],x2[0],x2[1],x2[2]));
      pi2=-std::norm(wf->EvalFTDensityMatrix1Alpha(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]));
      pi2-=std::norm(wf->EvalFTDensityMatrix1Beta(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]));
      pi2+=(pia1*pia2);
      pi2+=(pib1*pib2);
      pi2+=(pia1*pib2);
      pi2+=(pib1*pia2);
      pi2*=0.5e0;
      x1[0]+=dx;
      x2[1]+=(0.87*dx);
   }
   timer.End();
   fld2eval="Pi2OpenShMomSpRaw";
   epsstr="'~{/Symbol r}{0.1\\~}_{2,OSraw}'";
   if ( verbose ) { timer.PrintElapsedTimeMilliSec(fld2eval); }
   reftime=timer.GetElapsedTimeSec();
   property.push_back(fld2eval);
   epslabels.push_back(epsstr);
   computeTime.push_back(reftime);
   // ---------------------------------------------------

   ofstream ofil("timeprofile.dat");
   ofil << scientific << setprecision(8);
   ofil << "#set xtics (";
   for ( size_t i=0 ; i<property.size() ; ++i ) {
      ofil << epslabels[i] << ' ' << i << (i<(property.size()-1) ? ", " : ")\n");
   }
   ofil << "#Property    ComputeTime(s)     ComputeTimeRatio(vs rho)" << '\n';
   for ( size_t i=0 ; i<property.size() ; ++i ) {
      ofil << i << ' ' << std::setw(18) << property[i] << ' ' << computeTime[i] << ' '
           << (computeTime[i]/computeTime[0]) << '\n';
   }
   ofil.close();

   /* All OK  */
   ScreenUtils::PrintHappyEnding();
   globtimer.End();
   globtimer.PrintElapsedTimeSec(string("global timer"));
   /* ************************************************************************** */
   return EXIT_SUCCESS;
}

