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
/*
   Adscription at the moment this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
 */
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <iomanip>
using std::setprecision;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>
using std::string;
#include <cmath>

#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "wfgrid2d.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helpersplot.h"

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,outfilnam,gnpnam;
   string progname;
   OptionFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam,gnpnam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   GaussWaveFunction gwf;
   if (!(gwf.ReadFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.SetUpBNW();             //To setup the bond network.
   
   WaveFunctionGrid2D grid;    //Defining a grid object
   
   /* Looking for user grid dimensions */
   
   int nn=DEFAULTPOINTSPERDIRECTION;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&nn);
      grid.SetNPts(nn);
   } else {
      grid.SetNPts(nn);
   }
   
   /* Defining the line direction and other properties, and setting up the grid */
   
   int at1=0,at2=1,at3=2;
   if (bnw.nNuc==1) {
      at1=at2=0,at3=0;
      grid.SetUpSimplePlane(bnw,0);
      cout << "The file " << infilnam << " has only one atom" << endl;
      cout << "Using this atom to set the plane...\n";
   } else if (bnw.nNuc==2) {
      at1=0;
      at2=at3=1;
      grid.SetUpSimplePlane(bnw,0,1);
      cout << "The file " << infilnam << " has only two atoms" << endl;
      cout << "Using these atoms to set the plane...\n";
   } else {
      if (options.setats) {
         sscanf(argv[options.setats],"%d",&at1);
         sscanf(argv[options.setats+1],"%d",&at2);
         sscanf(argv[options.setats+2],"%d",&at3);
         if (at1<1||at2<1||at3<1||at1>bnw.nNuc||at2>bnw.nNuc||at3>bnw.nNuc) {
            ScreenUtils::SetScrRedBoldFont();
            cout << "Error: one of the given atoms do not exist!\n";
            ScreenUtils::SetScrNormalFont();
            exit(1);
         }
         cout << "Using atoms " << at1 << "(" << bnw.atLbl[at1-1]
         << "), " << at2 << "(" << bnw.atLbl[at2-1] << "), and "
         << at3 << "(" << bnw.atLbl[at3-1] << ") to set the plane."
         << endl;
         at1--;
         at2--;
         at3--;
         grid.SetUpSimplePlane(bnw,at1,at2,at3);
      } else {
         at1=at2=at3=0;
         grid.SetUpSimplePlane(bnw,0);
         cout << "Using the first atom to set the line...\n";
      }
   }
   //cout << "checkpoint" << endl;
   cout << "The size of the grid will be: " << grid.GetNPts(0) << "x" << grid.GetNPts(1) << endl;
   cout << "Total number of points that will be computed: " << grid.GetNPts(0)*grid.GetNPts(1) << endl;
   
   /* Setting the property to be computed */
   
   char prop;
   if (options.prop2plot) {
      prop=argv[options.prop2plot][0];
   } else {
      prop='d';
   }
   
   /* Main calculation loop, chooses between different available fields. */
   
   cout << "Evaluating and writing property..." << endl;
   cout << "(Scalar Field to plot: " << GetFieldTypeKeyLong(prop) << ")." << endl << endl;
   ScalarFieldType fldType=Char2ScalarFieldType(prop);
   grid.MakeTsv(outfilnam,gwf,fldType);
   cout << endl << "Output written in file: " << outfilnam << endl;
   
#if _HAVE_GNUPLOT_
   cout << "          Gnuplot file: " << gnpnam << endl;
   double dd=0.0e0,r[3];
   for (int i=0; i<3; i++) {
      r[i]=grid.Ca[i]-grid.Cb[i];
      dd+=r[i]*r[i];
   }
   dd=0.45e0*sqrt(dd);
   if (options.mkplt) {
      HelpersPlot::makeGnuplotFile(options,gnpnam,outfilnam,prop,dd,bnw,at1,at2,at3,grid);
   }
#endif /* _HAVE_GNUPLOT_ */

#if (defined(__APPLE__)||defined(__linux__))
   if (options.zipdat) {
      string cmdl;
      cmdl=string("gzip -9f ")+outfilnam;
      cout << "Calling gzip...";
      system(cmdl.c_str());
      cout << " Done!" << endl;
   }
#endif/* defined(__APPLE__)||defined(__linux__) */
   
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


