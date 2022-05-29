#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "helpersintegrate.h"
#include "integrator_vegas.h"
#include <climits>

shared_ptr<Integrator> FactoryIntegrator::CreateIntegrator(OptionFlags &options,\
      int argc,char* argv[],GaussWaveFunction &ugwf,BondNetWork &ubnw) {
   shared_ptr<IntegratorVegas> vegas=shared_ptr<IntegratorVegas>(new IntegratorVegas(ugwf,ubnw));
   shared_ptr<Integrator> integrator=shared_ptr<Integrator>(vegas);

   //Setting configuration parameters.
   int intervals=10;
   if ( options.setintervals ) { intervals=std::stod(string(argv[options.setintervals])); }
   size_t points=10000;
   if ( options.setpoints ) { points=std::stod(string(argv[options.setpoints])); }
   size_t iterations=20;
   if ( options.setiterations ) { iterations=std::stod(string(argv[options.setiterations])); }
   double convRate=1.0e0;
   if ( options.setconvergenceRate ) { convRate=std::stod(string(argv[options.setconvergenceRate])); }
   int terma=0;
   if ( options.settermalization ) { terma=std::stod(string(argv[options.settermalization])); }
   double tol=0.0e0;
   if ( options.settolerance ) { tol=std::stod(string(argv[options.settolerance])); }
   int stopRef=INT_MAX;
   if ( options.setstopRefinement ) { stopRef=std::stod(string(argv[options.setstopRefinement])); }
   char func='d';
   if ( options.setfunction ) { func=*argv[options.setfunction]; }
   double nelectrons=0.0e0;
   if ( func == 'd' || func == 'm' ) {
      nelectrons=ugwf.IntegralRho();
      vegas->AnalyticIntegral(nelectrons);
   }
   size_t nPntsForMax=100000;
   if ( options.setNPntsForMax ) { nPntsForMax=std::stod(string(argv[options.setNPntsForMax])); }
   double boxLimits[2]={0.0e0,0.0e0};
   if ( options.setSupBoxLimits ) {
      boxLimits[0]={std::stod(string(argv[options.setInfBoxLimits]))};
      boxLimits[1]={std::stod(string(argv[options.setSupBoxLimits]))};
   }

   //Setting integration boundaries.
   double intRmin[3]={0.0e0,0.0e0,0.0e0},intRmax[3]={0.0e0,0.0e0,0.0e0};
   double uncertainty = 3.0e0;
   if ( boxLimits[0] == 0.0e0 && boxLimits[1] == 0.0e0 ){
      if ( (func == 'm') || (func == 'T') || (func == 'k') ){
         while ( ugwf.EvalFTDensity(intRmax[0],intRmax[1],intRmax[2]) >= 1.0e-12 ) {
            for (int i=0; i<3; ++i) { intRmax[i] += 0.5e0; }
         }
         for (int i=0; i<3; ++i) { intRmin[i] = -intRmax[i]; }
      } else {
         for (int i=0;i<3;++i) {
            intRmin[i] = ubnw.bbmin[i]-uncertainty;
            intRmax[i] = ubnw.bbmax[i]+uncertainty;
         }
      }
   } else {
      for (int i=0;i<3;++i) {
         intRmin[i] = boxLimits[0];
         intRmax[i] = boxLimits[1];
      }
   }
   vegas->SetDimensions(intRmin[0],intRmin[1],intRmin[2],intRmax[0],intRmax[1],intRmax[2]);

   //Setting integration properties.
   vegas->SetIntegrand(func);
   vegas->SetIntervals(intervals);
   vegas->SetNumOfPoints(points);
   vegas->SetIterations(iterations);
   vegas->SetConvergenceRate(convRate);
   vegas->SetTermalization(terma);
   vegas->SetTolerance(tol);
   vegas->SetStopRefinement(stopRef);
   vegas->SetNSamplesToFindMaximum(nPntsForMax);
   // vegas->NormalizedEDF();
   // vegas->Relative2MaxDensity('a'); //Average of maxima.
   return integrator;
}

