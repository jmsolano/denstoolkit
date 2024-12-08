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

/*
 Adscription at the moment this program is started:
 University of Guelph,
 Guelph, Ontario, Canada.
 May 2013
 */
#include <iostream>
using std::cout;
using std::ios;
#include <cmath>
#include <string>
#include <iomanip>
using std::setprecision;
#include <ctime>
#include "screenutils.h"
#include "bondnetwork.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helpersdrawer.h"

int main (int argc, char *argv[]) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,basenam;
   string progname;
   OptionFlags options;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,basenam); //This creates the names used.
   if (options.verboseLevel && std::stoi(string(argv[options.verboseLevel]))>0)
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << "\nLoading data from file: " << infilnam << "...\n";
   
   BondNetWork bnw;
   if (!bnw.ReadFromFile(infilnam)) { return EXIT_FAILURE; } //Loading the bond-network (if the wavefunction
                               // or xyz file
                               //was read, there souldn't be problems here).
   bnw.SetUpBNW();             //To setup the bond network.

   HelpersDrawer::MakePovFile(basenam,options,argv,bnw);

   /* At this point the computation has ended. Usually this means no errors ocurred. */
   
   if (options.verboseLevel && std::stoi(string(argv[options.verboseLevel]))>0) {
      ScreenUtils::SetScrGreenBoldFont();
      ScreenUtils::PrintHappyEnding();
      ScreenUtils::PrintScrStarLine();
      cout << setprecision(3) << "CPU Time: "
         << double( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << '\n';
      double end_walltime=time(NULL);
      cout << "Wall-clock time: " << double (end_walltime-begin_walltime) << "s" << '\n';
#if DEBUG
      cout << "Debuggin mode (under construction...)" << '\n';
#endif
      ScreenUtils::PrintScrStarLine();
      ScreenUtils::SetScrNormalFont();
   }
   return EXIT_SUCCESS;
}

