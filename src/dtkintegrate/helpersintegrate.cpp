/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
#include "fldtypesdef.h"
#include "helpersintegrate.h"
#include "mymath.h"
#include "matrixvectoroperations3d.h"
#include "atomradiicust.h"
#include "dtkscalarfunction3d.h"
#include "integrator3d_vegas.h"
#include "integrator3d_miser.h"
#include "integrator3d_legsphtd.h"
#include "../common/integrator3d_diatomics.h"
#include "commonhelpers.h"

shared_ptr<Integrator3D> FactoryIntegrator::CreateIntegrator(OptionFlags &options,\
      int argc,char* argv[],GaussWaveFunction &ugwf,BondNetWork &ubnw) {
   char integtype='v';
   if ( options.integrator ) { integtype=argv[options.integrator][0]; }
   switch ( integtype ) {
      case 'd' :
         return CreateIntegratorDiatomics(options,argc,argv,ugwf,ubnw);
         break;
      case 'm' :
         return CreateIntegratorMiser(options,argc,argv,ugwf,ubnw);
         break;
      case 's' :
         return CreateIntegratorCubLegSphtDes(options,argc,argv,ugwf,ubnw);
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
   vegas->SetInternalIntegrand(ft);

   //Setting configuration parameters.
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
   //double sqrt3=sqrt(3.0e0);
   //for ( size_t i=0 ; i<3 ; ++i ) { intRmin[i]*=sqrt3; }
   //for ( size_t i=0 ; i<3 ; ++i ) { intRmax[i]*=sqrt3; }
   vegas->SetDimensions(intRmin[0],intRmin[1],intRmin[2],intRmax[0],intRmax[1],intRmax[2]);
   cout << "Integrand: " << GetFieldTypeKeyLong(ft) << '\n';

   // vegas->NormalizedEDF();
   // vegas->Relative2MaxDensity('a'); //Average of maxima.
   return integrator;
}
void FactoryIntegrator::FindIntegralLimits(OptionFlags &options,char*argv[],\
         GaussWaveFunction &wf,BondNetWork &bn,char ft,\
         vector<double> &rmin,vector<double> &rmax) {
   if ( options.setupperdombox ) {
      double boxLimits[2]={0.0e0,0.0e0};
      boxLimits[0]={std::stod(string(argv[options.setlowerdombox]))};
      boxLimits[1]={std::stod(string(argv[options.setupperdombox]))};
      for ( size_t i=0 ; i<3 ; ++i ) { rmin[i]=boxLimits[0]; }
      for ( size_t i=0 ; i<3 ; ++i ) { rmax[i]=boxLimits[1]; }
      return;
   }
   if ( Is3DMomSpaceField(ft) ) {
      while ( wf.EvalFTDensity(rmax[0],rmax[1],rmax[2]) >= 1.0e-14 ) {
         for (int i=0; i<3; ++i) { rmax[i] += 0.5e0; }
      }
      for (int i=0; i<3; ++i) { rmin[i]=-rmax[i]; }
   } else {
      while ( wf.EvalDensity(rmax[0],rmax[1],rmax[2]) >= 1.0e-14 ) {
         for (int i=0; i<3; ++i) { rmax[i] += 0.5e0; }
      }
      for (int i=0; i<3; ++i) { rmin[i]=-rmax[i]; }
   }
}
void FactoryIntegrator::DetermineDiatomicIntegralLimits(GaussWaveFunction &wf,\
      char ft,vector<double> &r0mx,vector<double> &r1mx,vector<double> &rmid) {
   vector<double> r0(3);
   for ( int i=0 ; i<3 ; ++i ) { r0[i]=wf.GetR(0,i); }
   vector<double> r1(3);
   for ( int i=0 ; i<3 ; ++i ) { r1[i]=wf.GetR(1,i); }
   double a0=GetAtomicVDWRadiusAtomicUnits(wf.atCharge[0]-1);
   double a1=GetAtomicVDWRadiusAtomicUnits(wf.atCharge[1]-1);
   CommonHelpers::DetermineDiatomicIntegralLimits(wf,'d',r0,r1,a0,a1,r0mx,r1mx,rmid);
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
   cout << "Integrand: " << GetFieldTypeKeyLong(ft) << '\n';


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
shared_ptr<Integrator3D> FactoryIntegrator::CreateIntegratorCubLegSphtDes(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw) {
   char ft='d';
   if ( options.integrand ) { ft=argv[options.integrand][0]; }
   if ( ugwf.nNuc !=1 && Is3DPosSpaceField(ft) ) {
      ScreenUtils::DisplayWarningMessage("Legendre-SphericaltDesign cubature only works with\n"
            "single atoms! The result may be mistaken.");
      cout << "Integrand: " << ft << '\n';
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   shared_ptr<DTKScalarFunction> dtkfield=shared_ptr<DTKScalarFunction>(new DTKScalarFunction(ugwf));
   dtkfield->SetScalarFunction(ft);
   shared_ptr<Function3D> the_function=shared_ptr<Function3D>(dtkfield);
   shared_ptr<Integrator3DLegSphtDes> legsphtd=shared_ptr<Integrator3DLegSphtDes>(new Integrator3DLegSphtDes(the_function));
   shared_ptr<Integrator3D> integrator=shared_ptr<Integrator3D>(legsphtd);
   cout << "Integrand: " << GetFieldTypeKeyLong(ft) << '\n';

   vector<double> intRmin(3),intRmax(3);
   for ( size_t i=0 ; i<3 ; ++i ) { intRmin[i]=0.0e0; intRmax[i]=0.0e0; }
   FindIntegralLimits(options,argv,ugwf,ubnw,ft,intRmin,intRmax);
   double a=-1.0e+50;
   for ( size_t i=0 ; i<3 ; ++i ) { if ( fabs(intRmin[i])>a ) { a=fabs(intRmin[i]); } }
   for ( size_t i=0 ; i<3 ; ++i ) { if ( fabs(intRmax[i])>a ) { a=fabs(intRmax[i]); } }
   a*=sqrt(3.0e0);
   cout << "a: " << a << '\n';
   cout << "f(" << a << ",0.0e0,0.0e0): " << dtkfield->f(a,0.0e0,0.0e0) << '\n';

   int glord=64;
   if ( options.lsptdsetol ) { glord=std::stoi(string(argv[options.lsptdsetol])); }
   int sptdord=21;
   if ( options.lsptdsetos ) { sptdord=std::stoi(string(argv[options.lsptdsetos])); }

   if ( !(legsphtd->SetupCubature(0.0e0,a,glord,sptdord)) ) {
      ScreenUtils::DisplayWarningMessage("Cubature could not be setup!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return nullptr;
   }
   //legsphtd->SetVerbosityLevel(1);
   return integrator;
}
shared_ptr<Integrator3D> FactoryIntegrator::CreateIntegratorDiatomics(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw) {
   if ( ugwf.nNuc !=2 ) {
      ScreenUtils::DisplayErrorMessage("Diatomics cubature only works with\n"
            "diatomic molelcules! The result would be wrong, thus the\n"
            "integrator will not be created.");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return nullptr;
   }
   char ft='d';
   if ( options.integrand ) { ft=argv[options.integrand][0]; }
   if ( !CommonHelpers::AtomsAreZAligned(ugwf) ) {
      ScreenUtils::DisplayErrorMessage("The atoms are not aligned along the z-axis.\n"
            "The integrator cannot be constructed...");
      return nullptr;
   }
   vector<double> intRmx0(3),intRmx1(3),intRmid(3);
   DetermineDiatomicIntegralLimits(ugwf,ft,intRmx0,intRmx1,intRmid);
   shared_ptr<DTKScalarFunction> dtkfield=shared_ptr<DTKScalarFunction>(new DTKScalarFunction(ugwf));
   dtkfield->SetScalarFunction(ft);
   shared_ptr<Function3D> the_function=shared_ptr<Function3D>(dtkfield);
   shared_ptr<Integrator3DDiatomics> diat=shared_ptr<Integrator3DDiatomics>(new Integrator3DDiatomics(the_function));
   shared_ptr<Integrator3D> integrator=shared_ptr<Integrator3D>(diat);
   cout << "Integrand: " << GetFieldTypeKeyLong(ft) << '\n';
   if ( options.verboseLevel>0 ) {
      ScreenUtils::SetScrYellowBoldFont();
      cout << "intRmx0: " << intRmx0[0] << ' ' << intRmx0[1] << ' ' << intRmx0[2] << '\n';
      cout << "intRmx1: " << intRmx1[0] << ' ' << intRmx1[1] << ' ' << intRmx1[2] << '\n';
      cout << "prop(intRmx0): " << the_function->f(intRmx0) << '\n';
      cout << "prop(intRmx1): " << the_function->f(intRmx1) << '\n';
      ScreenUtils::SetScrNormalFont();
   }

   int glradord=32;
   if ( options.diatsetrad ) { glradord=std::stoi(string(argv[options.diatsetrad])); }
   int glangord=32;
   if ( options.diatsetang ) { glangord=std::stoi(string(argv[options.diatsetang])); }

   double rt[3],xy[3];
   double vdvf=0.75;
   double r0=vdvf*GetAtomicVDWRadiusAtomicUnits(ubnw.atNum[0]);
   for ( int i=0 ; i<3 ; ++i ) { rt[i]=ubnw.R[0][i]-intRmid[i]; }
   if ( options.verboseLevel>0 ) {
      cout << "Dist(x0-c): " << magV3(rt) << ", r0: " << r0 << ", r0/D: " << (r0/magV3(rt)) << '\n';
   }
   if ( r0>magV3(rt) ) {
      r0=vdvf*magV3(rt);
      ScreenUtils::DisplayWarningMessage("The middle plane is closer than 0.75VDWRad!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   double r1=vdvf*GetAtomicVDWRadiusAtomicUnits(ubnw.atNum[1]);
   for ( int i=0 ; i<3 ; ++i ) { rt[i]=ubnw.R[1][i]-intRmid[i]; }
   if ( options.verboseLevel>0 ) {
      cout << "Dist(x0-c): " << magV3(rt) << ", r1: " << r1 << ", r1/D: " << (r1/magV3(rt)) << '\n';
   }
   if ( r1>magV3(rt) ) {
      r1=vdvf*magV3(rt);
      ScreenUtils::DisplayWarningMessage("The middle plane is closer than 0.75VDWRad!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
   }
   for ( int i=0 ; i<3 ; ++i ) { xy[i]=intRmx0[i]-intRmid[i]; }
   double RR=magV3(xy);
   for ( int i=0 ; i<3 ; ++i ) { xy[i]=intRmx1[i]-intRmid[i]; }
   if ( RR<magV3(xy) ) { RR=magV3(xy); }
   for ( int i=0 ; i<3 ; ++i ) { intRmx0[i]=ubnw.R[0][i]; } //hereon, intRmx0 is x0.
   for ( int i=0 ; i<3 ; ++i ) { intRmx1[i]=ubnw.R[1][i]; } //hereon, intRmx1 is x1.
   if ( !(diat->SetupCubature(intRmx0,intRmx1,intRmid,r0,r1,RR,glradord,glangord)) ) {
      ScreenUtils::DisplayWarningMessage("Cubature could not be setup!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return nullptr;
   }
   //diat>SetVerbosityLevel(1);
   return integrator;
}

