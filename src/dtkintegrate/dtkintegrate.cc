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
#include <memory>
using std::shared_ptr;
#include <ctime>
#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "../common/integrator_vegas.h"
#include "optflags.h"
#include "../common/bondnetwork.h"
#include "mytimer.h"
#include "crtflnms.h"
#include "helpersintegrate.h"

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

   shared_ptr<Integrator> integrator=FactoryIntegrator::CreateIntegrator(options,argc,argv,gwf,bnw);
   cout << scientific << setprecision(10);
   integrator->DisplayProperties();

   //Numeric integral.
   MyTimer aTim;
   aTim.Start();
   integrator->Integrate();
   aTim.End();
   aTim.PrintElapsedTimeSec(string("Integration time"));

   //Display results on screen.
   integrator->DisplayResults();

   //Write results in a log file.
   cout << "\nSaving results in file '" << outfilnam << "'...\n";
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
   integrator->WriteResults(ofil);
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

