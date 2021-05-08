/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
using std::exit;
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
using std::scientific;
#include <ctime>

#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymemory.h"
#include "../common/iofuncts-wfx.h"
#include "../common/iofuncts-wfn.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "../common/solcubetools.h"
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
   int dim=0,axis=0,plane=0,npts=DEFAULTPOINTSPERDIRECTION;
   char field='d';
   double px,py,pz;
   px=py=pz=0.0e0;
   string extralbl="";
   if (options.evdim) {
      switch (argv[options.evdim][1]) {
         case '0':
            dim=0;
            sscanf (argv[options.evdim+1],"%lf",&px);
            sscanf (argv[options.evdim+2],"%lf",&py);
            sscanf (argv[options.evdim+3],"%lf",&pz);
            //cout << scientific << "px: " << px << ", py: " << py << ", pz: " << pz << endl;
            break;
         case '1':
            dim=1;
            switch (argv[options.evdim+1][0]) {
               case 'x':
                  axis=1;
                  break;
               case 'y':
                  axis=2;
                  break;
               case 'z':
                  axis=3;
                  break;
               default:
                  cout << endl;
                  ScreenUtils::DisplayErrorMessage("Invalid value for option -1...");
                  cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
                  cout << "\nto view the help menu.\n\n";
                  exit(1);
                  break;
            }
            extralbl=argv[options.evdim+1];
            extralbl.insert(0,"-P");
            break;
         case '2':
            dim=2;
            extralbl=string(argv[options.evdim+1]);
            if (!(extralbl==string("xy")||extralbl==string("xz")||extralbl==string("yz"))) {
               cout << endl;
               ScreenUtils::DisplayErrorMessage("Invalid value for option -2...");
               cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
               cout << "\nto view the help menu.\n\n";
               exit(1);
            }
            if (extralbl==string("xy")) {plane=1;}
            if (extralbl==string("xz")) {plane=2;}
            if (extralbl==string("yz")) {plane=3;}
            extralbl.insert(1,"P");
            extralbl.insert(0,"-P");
            break;
         case '3':
            dim=3;
            npts=DEFAULTNPOINTSFORCUBE;
            break;
         default:
            break;
      }
   }
   if ( options.setfld ) {field=argv[options.setfld][0];}
   if (options.setn1) {sscanf (argv[options.setn1],"%d",&npts);}
   mkFileNames(argv,options,infilnam,outfilnam,gnpnam,dim,extralbl); //This creates the names used.
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
   
   string strfield;
   switch ( field ) {
      case 'd' :
         strfield="Momentum Density";
         break;
      case 'K' :
         strfield="Kinetic Energy Density (in Momentum Space)";
         break;
      default :
         break;
   }
   
   /* Evaluates the momentum density at a single point */
   if (dim==0) {
      cout << scientific << setprecision(12);
      cout << "\nMomentum density at\nP=(" << px << "," << py << "," << pz << "): " << endl;
      cout << endl << gwf.EvalFTDensity(px,py,pz) << endl;
      cout << "\nK. E. in Momentum space at \nP=(" << px << "," << py << "," << pz << "): " << endl;
      cout << endl << gwf.EvalFTKineticEnergy(px,py,pz) << endl;
   }
   
   /* Evaluates the momentum density on a line */
   if (dim==1) {
      cout << "The size of the grid will be " << npts << endl;
      cout << "The total number of points that will be computed is " << npts << endl;
      cout << "Evaluating and writing " << strfield << " on a line..." << endl;
      HelpersPlot::MakeLineDatFile(options,outfilnam,gwf,axis,npts,field);
      if (options.mkplt) {HelpersPlot::MakeLineGnuplotFile(options,gnpnam,outfilnam,field);}
   }
   
   /* Evaluates the momentum density on a plane */
   if (dim==2) {
      cout << "The size of the grid will be " << npts << " x " << npts << endl;
      cout << "The total number of points that will be computed is " << npts*npts << endl;
      cout << "Evaluating and writing " << strfield << " on a plane..." << endl;
      HelpersPlot::MakePlaneTsvFile(options,outfilnam,gwf,plane,npts,field);
      double maxdim=DEFAULTMAXVALUEOFP;
      HelpersPlot::MakePlaneGnuplotFile(options,gnpnam,outfilnam,maxdim,field);
   }
   
   /* Evaluates the momentum density on a cube */
   if (dim==3) {
      HelpersPlot::MakeCubeFile(options,outfilnam,gwf,npts,field,strfield);
      if (options.mkplt) {
         ScreenUtils::DisplayWarningMessage("Plotting can only be performed with options -1 or -2.");
      }
   }
   
   cout << "Data saved in file: " << outfilnam << endl;
   
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


