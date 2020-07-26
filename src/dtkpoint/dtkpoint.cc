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

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

//#ifdef DEBUG
//#undef DEBUG
//#endif
//#define DEBUG 1
//#define USEPROGRESSBAR 1
//#define EPSFORELFVALUE (2.871e-05)

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <math.h>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;
#include <ctime>

#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "optflags.h"
#include "crtflnms.h"

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   GaussWaveFunction gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   solreal Nofelec=gwf.integralRho();
   solreal totNucCharg=gwf.totalNuclearCharge();
   solreal errinteg=100.0e0*fabs(1.0e0-Nofelec/totNucCharg);
   cout << setprecision(12);
   cout << "     Total Nuclear Charge: " << totNucCharg << endl;
   cout << "N. Electrons (Integrated): " << Nofelec << endl;
   cout << setprecision(4);
   cout << "    Rel. Err. Integration: " << errinteg << "%" << endl;
   
   /* The following is left for future reference and possible implementation of 
      additional capabilities. 
   */
   
   //bondNetWork bnw;
   //bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   //bnw.setUpBNW();             //To setup the bond network.
   
   /* Setting the property to be computed */

   /*
   char prop;
   if (options.prop2plot) {
      prop=argv[options.prop2plot][0];
   } else {
      prop='a';
   }
   // */

   /* Checking whether custom fields should be computed  */

   if ( options.setscustfld ) {gwf.useScalarCustomField(true);}
   if ( options.setvcustfld ) {gwf.useVectorCustomField(true);}
   /* Openning the output log-file.  */
   
   ofile.open(outfilnam.c_str(),ios::out);
   FileUtils::WriteHappyStart(argv,ofile,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   ofile << endl << "       Wave function file: " << infilnam << endl;
   ofile << "                    Title: " << gwf.title[0] << endl;
   ofile << scientific << setprecision(12);
   ofile << "     Total Nuclear Charge: " << totNucCharg << endl;
   ofile << "N. Electrons (Integrated): " << Nofelec << endl;
   ofile << setprecision(4);
   ofile << "    Rel. Err. Integration: " << errinteg << "%" << endl;
   
   solreal x,y,z;
   x=y=z=0.0e0;
   if ((!options.rcrds)&&(!options.setat)&&(!options.crdfil)) {
      ScreenUtils::SetScrYellowBoldFont();
      cout << endl;
      cout << "No coordinates have been given, calculating all properties at the origin." << endl;
      cout << endl;
      ScreenUtils::SetScrNormalFont();
   }
   
   if (options.rcrds) {
      cout << "Reading coordinates from the command line..." << endl;
      int pos=options.rcrds;
      sscanf(argv[pos++],"%lf",&x);
      sscanf(argv[pos++],"%lf",&y);
      sscanf(argv[pos],"%lf",&z);
   }
   
   if (options.setat) {
      if (options.rcrds) {
         ScreenUtils::DisplayWarningMessage("Overwriting command line coordinates...");
      }
      int nat;
      sscanf(argv[options.setat],"%d",&nat);
      if (nat>gwf.nNuc) {
         ScreenUtils::DisplayErrorMessage("You requested a non existent atom!");
         ScreenUtils::DisplayWarningMessage("Setting R to be the zero vector.");
      } else {
         nat--;
         cout << "Setting R=R(" << gwf.atLbl[nat] << ")" << endl;
         ofile << "                        R= R(" << gwf.atLbl[nat] << ")" << endl;
         nat*=3;
         x=gwf.R[nat++];
         y=gwf.R[nat++];
         z=gwf.R[nat];
         sscanf(argv[options.setat],"%d",&nat);
      }
   }
   
   if (options.crdfil) {
      ifstream rfile;
      rfile.open(argv[options.crdfil],ios::in);
      cout << "Coordinates will be taken from " << argv[options.crdfil] << endl;
      ofile << "    File of Coords: " << argv[options.crdfil] << endl;
      if (options.setat||options.rcrds) {
         cout << "Ignoring previous coordinates..." << endl;
      }
      cout << endl;
      ofile << endl;
      while (!rfile.eof()) {
         FileUtils::DiscardComments(rfile);
         rfile >> x >> y >> z;
         gwf.displayAllFieldProperties(x,y,z);
         cout << endl;
         gwf.writeAllFieldProperties(x,y,z,ofile);
         ofile << endl;
      }
      rfile.close();
   } else {
      cout << endl;
      ofile << endl;
      gwf.displayAllFieldProperties(x,y,z);
      gwf.writeAllFieldProperties(x,y,z,ofile);
   }
   ofile.close();
   
   cout << endl << "Output written in file: " << outfilnam << endl;
   
   
   /* At this point the computation has ended. Usually this means no errors ocurred. */
   
   ScreenUtils::PrintHappyEnding();
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintScrStarLine();
   cout << setprecision(3) << "CPU Time: "
        << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   solreal end_walltime=time(NULL);
   cout << "Wall-clock time: " << solreal (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debuggin mode (under construction...)" << endl;
#endif
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrNormalFont();
   return 0;
}

