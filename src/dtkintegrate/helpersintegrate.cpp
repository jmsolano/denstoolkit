/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <climits>
#include "helpersintegrate.h"
#include "../common/dtkscalarfunction3d.h"
#include "../common/integrator3d_vegas.h"
#include "../common/integrator3d_miser.h"

shared_ptr<Integrator3D> FactoryIntegrator::CreateIntegrator(OptionFlags &options,\
      int argc,char* argv[],GaussWaveFunction &ugwf,BondNetWork &ubnw) {
   char integtype='v';
   if ( options.integrator ) { integtype=argv[options.integrator][0]; }
   switch ( integtype ) {
      case 'm' :
         return CreateIntegratorMiser(options,argc,argv,ugwf,ubnw);
         break;
      case 'v' :
         return CreateIntegratorVegas(options,argc,argv,ugwf,ubnw);
         break;
      default :
         ScreenUtils::DisplayWarningMessage("Unknown integrator! Setting up Vegas integrator.");
         cout << __FILE__ << ", line: " << __LINE__ << '\n';
         break;
   }
   return CreateIntegratorVegas(options,argc,argv,ugwf,ubnw);
}
shared_ptr<Integrator3D> FactoryIntegrator::CreateIntegratorVegas(OptionFlags &options,\
      int argc,char* argv[],GaussWaveFunction &ugwf,BondNetWork &ubnw) {
   char ft='d';
   if ( options.integrand ) { ft=*argv[options.integrand]; }
   shared_ptr<DTKScalarFunction> dtkfield=shared_ptr<DTKScalarFunction>(new DTKScalarFunction(ugwf));
   dtkfield->SetScalarFunction(ft);
   shared_ptr<Function3D> the_function=shared_ptr<Function3D>(dtkfield);
   shared_ptr<Integrator3DVegas> vegas=shared_ptr<Integrator3DVegas>(new Integrator3DVegas(the_function));
   shared_ptr<Integrator3D> integrator=shared_ptr<Integrator3D>(vegas);

   //Setting configuration parameters.
   //vegas->param.integrand=ft;
   if ( ft == 'd' || ft == 'm' ) {
      vegas->AnalyticIntegral(ugwf.IntegralRho());
   }
   int intervals=10;
   if ( options.vegassetinterv ) { intervals=std::stod(string(argv[options.vegassetinterv])); }
   vegas->SetIntervals(intervals);
   size_t points=10000;
   if ( options.vegassetpoints ) { points=std::stod(string(argv[options.vegassetpoints])); }
   vegas->SetNumOfPoints(points);
   size_t iterations=20;
   if ( options.vegassetiter ) { iterations=std::stod(string(argv[options.vegassetiter])); }
   vegas->SetIterations(iterations);
   double convRate=1.0e0;
   if ( options.vegassetconvrat ) { convRate=std::stod(string(argv[options.vegassetconvrat])); }
   vegas->SetConvergenceRate(convRate);
   int therm=0;
   if ( options.vegassettherm ) { therm=std::stod(string(argv[options.vegassettherm])); }
   vegas->SetThermalization(therm);
   double tol=0.0e0;
   if ( options.vegassettol ) { tol=std::stod(string(argv[options.vegassettol])); }
   vegas->SetTolerance(tol);
   int stopRef=INT_MAX;
   if ( options.vegassetstopref ) { stopRef=std::stod(string(argv[options.vegassetstopref])); }
   vegas->SetStopRefinement(stopRef);
   size_t nPntsForMax=100000;
   if ( options.vegassetnpts4max ) { nPntsForMax=std::stod(string(argv[options.vegassetnpts4max])); }
   vegas->SetNSamplesToFindMaximum(nPntsForMax);
   vector<double> intRmin(3),intRmax(3);
   for ( size_t i=0 ; i<3 ; ++i ) { intRmin[i]=0.0e0; intRmax[i]=0.0e0; }
   FindIntegralLimits(options,argv,ugwf,ubnw,ft,intRmin,intRmax);
   vegas->SetDimensions(intRmin[0],intRmin[1],intRmin[2],intRmax[0],intRmax[1],intRmax[2]);

   // vegas->NormalizedEDF();
   // vegas->Relative2MaxDensity('a'); //Average of maxima.
   return integrator;
}
void FactoryIntegrator::FindIntegralLimits(OptionFlags &options,char*argv[],\
         GaussWaveFunction &wf,BondNetWork &bn,char ft,\
         vector<double> &rmin,vector<double> &rmax) {
   double boxLimits[2]={0.0e0,0.0e0};
   if ( options.setupperdombox ) {
      boxLimits[0]={std::stod(string(argv[options.setlowerdombox]))};
      boxLimits[1]={std::stod(string(argv[options.setupperdombox]))};
   }
   double uncertainty = 3.0e0;
   if ( boxLimits[0] == 0.0e0 && boxLimits[1] == 0.0e0 ){
      if ( (ft == 'm') || (ft == 'T') || (ft == 'k') ){
         while ( wf.EvalFTDensity(rmax[0],rmax[1],rmax[2]) >= 1.0e-12 ) {
            for (int i=0; i<3; ++i) { rmax[i] += 0.5e0; }
         }
         for (int i=0; i<3; ++i) { rmin[i]=-rmax[i]; }
      } else {
         for (int i=0;i<3;++i) {
            rmin[i] = bn.bbmin[i]-uncertainty;
            rmax[i] = bn.bbmax[i]+uncertainty;
         }
      }
   } else {
      for (int i=0;i<3;++i) {
         rmin[i] = boxLimits[0];
         rmax[i] = boxLimits[1];
      }
   }
}
shared_ptr<Integrator3D> FactoryIntegrator::CreateIntegratorMiser(OptionFlags &options,\
      int argc,char* argv[],GaussWaveFunction &ugwf,BondNetWork &ubnw) {
   char ft='d';
   if ( options.integrand ) { ft=argv[options.integrand][0]; }
   shared_ptr<DTKScalarFunction> dtkfield=shared_ptr<DTKScalarFunction>(new DTKScalarFunction(ugwf));
   dtkfield->SetScalarFunction(ft);
   shared_ptr<Function3D> the_function=shared_ptr<Function3D>(dtkfield);
   shared_ptr<Integrator3DMiser> miser=shared_ptr<Integrator3DMiser>(new Integrator3DMiser(the_function));
   shared_ptr<Integrator3D> integrator=shared_ptr<Integrator3D>(miser);
   cout << "Using " << GetFieldTypeKeyShort(ft) << '\n';


   //Setting configuration parameters.
   size_t points=500000;
   if ( options.misersetpoints ) { points=std::stoull(string(argv[options.misersetpoints])); }
   miser->SetNumPts(points);
   double tmp=0.05;
   if ( options.misersetdith ) { tmp=std::stod(string(argv[options.misersetdith])); }
   miser->SetDith(tmp);
   vector<double> intRmin(3),intRmax(3);
   for ( size_t i=0 ; i<3 ; ++i ) { intRmin[i]=0.0e0; intRmax[i]=0.0e0; }
   FindIntegralLimits(options,argv,ugwf,ubnw,ft,intRmin,intRmax);
   miser->SetXMin(intRmin);
   miser->SetXMax(intRmax);

   return integrator;
}

