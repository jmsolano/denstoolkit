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
   if ( options.vegassetinterv ) { intervals=std::stod(string(argv[options.vegassetinterv])); }
   size_t points=10000;
   if ( options.vegassetpoints ) { points=std::stod(string(argv[options.vegassetpoints])); }
   size_t iterations=20;
   if ( options.vegassetiter ) { iterations=std::stod(string(argv[options.vegassetiter])); }
   double convRate=1.0e0;
   if ( options.vegassetconvrat ) { convRate=std::stod(string(argv[options.vegassetconvrat])); }
   int therm=0;
   if ( options.vegassettherm ) { therm=std::stod(string(argv[options.vegassettherm])); }
   double tol=0.0e0;
   if ( options.vegassettol ) { tol=std::stod(string(argv[options.vegassettol])); }
   int stopRef=INT_MAX;
   if ( options.vegassetstopref ) { stopRef=std::stod(string(argv[options.vegassetstopref])); }
   char func='d';
   if ( options.integrand ) { func=*argv[options.integrand]; }
   double nelectrons=0.0e0;
   if ( func == 'd' || func == 'm' ) {
      nelectrons=ugwf.IntegralRho();
      vegas->AnalyticIntegral(nelectrons);
   }
   size_t nPntsForMax=100000;
   if ( options.vegassetnpts4max ) { nPntsForMax=std::stod(string(argv[options.vegassetnpts4max])); }
   double boxLimits[2]={0.0e0,0.0e0};
   if ( options.setupperdombox ) {
      boxLimits[0]={std::stod(string(argv[options.setlowerdombox]))};
      boxLimits[1]={std::stod(string(argv[options.setupperdombox]))};
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
   vegas->SetThermalization(therm);
   vegas->SetTolerance(tol);
   vegas->SetStopRefinement(stopRef);
   vegas->SetNSamplesToFindMaximum(nPntsForMax);
   // vegas->NormalizedEDF();
   // vegas->Relative2MaxDensity('a'); //Average of maxima.
   return integrator;
}

