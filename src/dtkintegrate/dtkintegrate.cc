/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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
/*
   Adscription at the moment this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
*/
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
using std::scientific;
#include <ctime>
#include <climits>
#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "../common/integrator_vegas.h"
#include "optflags.h"
#include "../common/bondnetwork.h"
#include "mytimer.h"
#include "crtflnms.h"

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,outfilnam;
   string progname;
   OptionFlags options;
   ifstream ifile;
   ofstream ofile;
  
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   //Just to let the user know that the initial configuration is OK
  
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
  
   GaussWaveFunction gwf;

   if (!(gwf.ReadFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   BondNetWork bnw;
   bnw.ReadFromFile(infilnam);
   bnw.SetUpBNW();
   cout << "Done." << endl;
  
   IntegratorVegas integrator(gwf,bnw);

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
      nelectrons=gwf.IntegralRho();
      integrator.AnalyticIntegral(nelectrons);
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
         while ( gwf.EvalFTDensity(intRmax[0],intRmax[1],intRmax[2]) >= 1.0e-12 ) {
            for (int i=0; i<3; ++i) { intRmax[i] += 0.5e0; }
         }
         for (int i=0; i<3; ++i) { intRmin[i] = -intRmax[i]; }
      } else {
         for (int i=0;i<3;++i) {
            intRmin[i] = bnw.bbmin[i]-uncertainty;
            intRmax[i] = bnw.bbmax[i]+uncertainty;
         }
      }
   } else {
      for (int i=0;i<3;++i) {
         intRmin[i] = boxLimits[0];
         intRmax[i] = boxLimits[1];
      }
   }
   integrator.SetDimensions(intRmin[0],intRmin[1],intRmin[2],intRmax[0],intRmax[1],intRmax[2]);

   //Setting integration properties.
   integrator.SetIntegrand(func);
   integrator.SetIntervals(intervals);
   integrator.SetNumOfPoints(points);
   integrator.SetIterations(iterations);
   integrator.SetConvergenceRate(convRate);
   integrator.SetTermalization(terma);
   integrator.SetTolerance(tol);
   integrator.SetStopRefinement(stopRef);
   integrator.SetNSamplesToFindMaximum(nPntsForMax);
   // integrator.NormalizedEDF();
   // integrator.Relative2MaxDensity('a'); //Average of maxima.
   cout << scientific << setprecision(10);
   integrator.DisplayProperties();

   //Numeric integral.
   MyTimer aTim;
   aTim.Start();
   integrator.Integrate();
   aTim.End();
   aTim.PrintElapsedTimeSec(string("integration time"));

   //Print results on screen.
   integrator.DisplayResults();

   //Print results in a log file.
   cout << "\nPrinting integrand data into file " << outfilnam << " ...\n";
   ofstream ofil(outfilnam);
   if ( !ofil.good() ) {
      ScreenUtils::DisplayErrorFileNotOpen(outfilnam);
      ofil.close();
      return EXIT_FAILURE;
   }
   FileUtils::WriteHappyStart(argv,ofil,CURRENTVERSION,"JMSA/JMHP/SAFR");
   FileUtils::WriteScrStarLine(ofil);
   ofil << "Wavefunction file: " << infilnam << '\n';
   FileUtils::WriteScrStarLine(ofil);
   integrator.WriteResults(ofil);
   ofil.close();

   cout << scientific << setprecision(8);

   /* At this point the computation has ended. Usually this means no errors ocurred. */

   ScreenUtils::PrintHappyEnding();
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintScrStarLine();
   cout << setprecision(3) << "CPU Time: "
      << double( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   double end_walltime=time(NULL);
   cout << "Wall-clock time: " << double (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debuggin mode (under construction...)" << endl;
#endif
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrNormalFont();
   return EXIT_SUCCESS;
}

