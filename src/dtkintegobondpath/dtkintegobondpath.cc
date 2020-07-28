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
   * Adscription at the moment DTK suite is started:
   * University of Guelph,
   * Guelph, Ontario, Canada.
   * May 2013
   *
   * Adscription at the moment this program is started:
   * Meritorious Autonomous University of Puebla
   * Puebla, Puebla, Mexico.
   * September 2016
*/
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <iomanip>
using std::setprecision;
#include <ctime>
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "../common/fldtypesdef.h"
#include "../common/integrateoverbondpath.h"
#include "../common/fileutils.h"
#include "optflags.h"
#include "crtflnms.h"

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,logfilnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   
   /* Processes the options from the command line.  */
   getOptions(argc,argv,options);

   /* This creates the names used.  */
   mkFileNames(argv,options,infilnam,logfilnam);

   /* Just to let the user know that the initial configuration is OK.  */
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   
   /* Loading the wave function  */
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   GaussWaveFunction gwf;
   if (!(gwf.ReadFromFile(infilnam))) {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   /* Loads the bond-network (if the wave function
    * was read, there souldn't be problems here to
    * setup the bond network).  */
   BondNetWork bnw;
   bnw.ReadFromFile(infilnam);
   bnw.SetUpBNW();

   /* Setups the critical point network object  */
   CritPtNetWork cpn(gwf,bnw);
   double stepSize=INTEGBP_BONDPATHSTEPSIZE; //defined in soldefines.h
   cpn.SetStepSizeBGP(stepSize);
   int numPtsBGPArray=INTEGBP_NUMBEROFPOINTSBGPARRAY*4;
   cpn.SetMaxGradPathNPts(numPtsBGPArray);
   cpn.SetCriticalPoints(DENS);
   cpn.SetBondPaths();

   /* Setups the IntegrateOverBondPath object  */
   vector<IntegrateOverBondPath*> integ;
   integ.push_back(new IntegrateOverBondPath(gwf,cpn,DENS));
   integ.push_back(new IntegrateOverBondPath(gwf,cpn,KEDK));
   integ.push_back(new IntegrateOverBondPath(gwf,cpn,KEDG));
   //integ.push_back(new IntegrateOverBondPath(gwf,cpn,MEPD));
   integ.push_back(new IntegrateOverBondPath(gwf,cpn,REDG));
   integ.push_back(new IntegrateOverBondPath(gwf,cpn,EDFTA));
   for ( size_t i=0 ; i<integ.size() ; ++i ) {
      cout << "Computing " << integ[i]->GetFieldTypeLabelShort() << " integrals, (" << i << " out of " << 
        integ.size() << ")..." << endl;
      integ[i]->ComputeAllIntegralsOverBondPaths();
   }
   for ( size_t i=0 ; i<integ.size() ; ++i ) {
      cout << "Total integral (" << integ[i]->GetFieldTypeLabelShort() 
           << "): " << integ[i]->GetBondPathIntegral() << endl;
   }
   /* Writes the integrals to a file  */
   double globalEnergy=0.0e0;
   if ( options.globalenergy ) {
      globalEnergy=std::stod(string(argv[options.globalenergy]));
   } 
   cout << "globalEnergy: " << globalEnergy << endl;
   ofstream ofil(logfilnam.c_str());
   for ( size_t i=0 ; i<integ.size() ; ++i ) {
      FileUtils::WriteScrStarLine(ofil);
      integ[i]->WriteIntegralValuesToFile(ofil,globalEnergy);
   }
   ofil.close();

   /* Cleans the integral vector  */
   for ( size_t i=0 ; i<integ.size() ; ++i ) {
      delete integ[i];
      integ[i]=NULL;
   }
   /* At this point the computation has ended. Usually this means no
    * errors ocurred. */
   
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintHappyEnding();
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

