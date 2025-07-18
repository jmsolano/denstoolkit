/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
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
#include <iomanip>
using std::setprecision;
using std::scientific;
#include <ctime>

#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymemory.h"
#include "../common/iofuncts-wfx.h"
#include "../common/iofuncts-wfn.h"
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/wfgrid3d.h"
#include "../common/vmdtools.h"
#include "optflags.h"
#include "crtflnms.h"

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,outfilnam,logfilnam;
   string progname;
   OptionFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam,logfilnam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK

   cout << endl << "Loading wave function from file: " << infilnam << endl;
   
   GaussWaveFunction gwf;
   if (!(gwf.ReadFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done" << '\n';
   if ( options.stpspindens && gwf.ihaveSingleSpinOrbs ) {
      gwf.CalcCabAAndCabB();
   }
   
   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.SetUpBNW(); //To setup the bond network.
   
   WaveFunctionGrid3D grid;
   
   /* Looking for user grid dimensions */
   
   int nn=DEFAULTPOINTSPERDIRECTION;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&nn); //same number of points per direction
      grid.SetNPts(nn);
   } else if (options.setn3) {
      int ny,nz;
      sscanf(argv[options.setn3],"%d",&nn);  //different number of points per direction
      sscanf(argv[options.setn3+1],"%d",&ny); //manually given by the user
      sscanf(argv[options.setn3+2],"%d",&nz);
      grid.SetNPts(nn,ny,nz);
   } else if (options.setsmcub1) {
      sscanf(argv[options.setsmcub1],"%d",&nn); //uses nn for determining the number of points
                     //of the largest axis. The rest will have proportional number to their length
   } else if ( options.setcentredcub ) {
      nn=DEFAULTNUMPTSCENTREDCUBE;
      grid.SetNPts(nn);
   } else {
      grid.SetNPts(nn);
   }
   if ( options.setextraspace ) { grid.SetExtraSpace(std::stod(string(argv[options.setextraspace]))); }
   
   if ((options.setsmcub)||(options.setsmcub1)) {
      grid.SetUpSmartCuboidGrid(gwf,bnw,nn);
   } else if ( options.setcentredcub ) {
     int at1,at2;
     sscanf(argv[options.setcentredcub],"%d",&at1);
     sscanf(argv[options.setcentredcub+1],"%d",&at2);
     grid.SetUpCenteredGrid(gwf,bnw,at1,at2,DEFAULTLENGTHCENTREDCUBE,nn);
   } else {
      grid.SetUpSimpleGrid(gwf,bnw);
   }
   
   cout << "The size of the grid will be: " << grid.GetNPts(0) << " x " 
        << grid.GetNPts(1) << " x " << grid.GetNPts(2) << endl;
   cout << "Total number of points that will be computed: " 
        << (grid.GetNPts(0)*grid.GetNPts(1)*grid.GetNPts(2)) << endl;
   
   /* Special configurations  */
   if ( options.configspecialnci ) {
      double ttt=std::stod(string(argv[options.configspecialnci]));
      gwf.SetNCIRhoMin(ttt);
      ttt=std::stod(string(argv[options.configspecialnci+1]));
      gwf.SetNCIRhoMax(ttt);
      ttt=std::stod(string(argv[options.configspecialnci+2]));
      gwf.SetNCISMax(ttt);
   }
#if DEBUG
      cout << "nciRhoMin: " << gwf.nciRhoMin << endl;
      cout << "nciRhoMax: " << gwf.nciRhoMax << endl;
      cout << "nciSMax: " << gwf.nciSMax << endl;
#endif
   /* Setting the property to be computed */
   
   char prop;
   if (options.prop2plot) {
      prop=argv[options.prop2plot][0];
   } else {
      prop='d';
   }
   if ( prop=='b' && (!gwf.ihaveCABSingleSpin) ) {
      ScreenUtils::DisplayErrorMessage("The alpha- and beta-spin density matrices could not\n"
            "be setup! Exiting...");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return EXIT_FAILURE;
   }
   
   /* Main calculation loop, chooses between different available fields. */
   
   cout << "Evaluating and writing property..." << endl;
   cout << "(Scalar Field to plot: " << GetFieldTypeKeyLong(prop) << ")." << endl << endl;

   grid.MakeCube(outfilnam,gwf,prop);
   /*
   switch (prop) {
      case 'd':
         grid.MakeCube(outfilnam,gwf,DENS);
         cout << endl;
         break;
      case 'g':
         grid.MakeCube(outfilnam,gwf,MGRD);
         break;
      case 'l':
         grid.MakeCube(outfilnam,gwf,LAPD);
         cout << endl;
         break;
      case 'e':
         grid.MakeCube(outfilnam,gwf,ELLPY);
         break;
      case 'E':
         grid.MakeCube(outfilnam,gwf,ELFD);
         break;
      case 'P' :
         grid.MakeCube(outfilnam,gwf,MLED);
         break;
      case 'r' :
         grid.MakeCube(outfilnam,gwf,ROSE);
         break;
      case 's' :
         grid.MakeCube(outfilnam,gwf,REDG);
         break;
      case 'S':
         grid.MakeCube(outfilnam,gwf,SENT);
         break;
      case 'L':
         grid.MakeCube(outfilnam,gwf,LOLD);
         break;
      case 'M':
         grid.MakeCube(outfilnam,gwf,MGLD);
         break;
      case 'G':
         grid.MakeCube(outfilnam,gwf,KEDG);
         break;
      case 'K':
         grid.MakeCube(outfilnam,gwf,KEDK);
         break;
      case 'u' :
         grid.MakeCube(outfilnam,gwf,SCFD);
         break;
      case 'V':
         grid.MakeCube(outfilnam,gwf,MEPD);
         break;
      case 'v':
         grid.MakeCube(outfilnam,gwf,VPED);
         break;
      case 'Z':
         grid.MakeCube(outfilnam,gwf,NCIS);
         break;
      case 'z':
         grid.MakeCube(outfilnam,gwf,NCIL);
         break;
      default:
         ScreenUtils::SetScrRedBoldFont();
         cout << "Error: The property \"" << prop << "\" does not exist!" << endl;
         ScreenUtils::SetScrNormalFont();
         exit(1);
         break;
   }
   // */
   
   cout << endl << "Output written in file: " << outfilnam << endl;
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   if (options.zipcube) {
      string cmdl;
      cmdl=string("gzip -9f ")+outfilnam;
      cout << "Calling gzip...";
      system(cmdl.c_str());
      cout << " Done!" << endl;
   }
#endif
   
   /* At this point the computation has ended. Usually this means no errors ocurred. */
   
   if (options.wrtlog) {
      ofstream lfil;
      lfil.open(logfilnam.c_str(),ios::out);
      FileUtils::WriteHappyStart(argv,lfil,CURRENTVERSION,PROGRAMCONTRIBUTORS);
      lfil << "#Wave function file name: " << endl << infilnam << endl;
      lfil << "#Number of primitives: "  << endl << gwf.nPri << endl;
      lfil << "#Grid dimensions:" << endl
           << grid.GetNPts(0) << " " << grid.GetNPts(1) << " " << grid.GetNPts(2) << endl;
      lfil << "#Total number of points in the cube:" << endl
           << (grid.GetNPts(0)*grid.GetNPts(1)*grid.GetNPts(2)) << endl;
      if ( prop=='z' || prop=='Z' ) {
         lfil << "#NCI Configuration:" << endl;
         lfil << "#nciRhoMin:\n" << gwf.nciRhoMin << endl;
         lfil << "#nciRhoMax:\n" << gwf.nciRhoMax << endl;
         lfil << "#nciSMax:\n" << gwf.nciSMax << endl;
      }
      lfil << "#CPU Time (sec):" << endl;
      lfil << scientific << setprecision(4)
           <<  double( clock () - begin_time ) / CLOCKS_PER_SEC << endl;
      lfil << "#Wall-clock Time (sec):" << endl;
      double tmp_walltime=time(NULL);
      lfil << double (tmp_walltime-begin_walltime) << endl;
      lfil.close();
   }
   if ( options.genvmdscript ) {
      VMDTools::WriteSimpleVMDScript(outfilnam,prop,options.quietrender);
   }
   
   ScreenUtils::PrintHappyEnding();
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintScrStarLine();
   cout << setprecision(3) << "CPU Time: " << double( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   double end_walltime=time(NULL);
   cout << "Wall-clock time: " << double (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debugging mode (under construction)..." << endl;
#endif
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrNormalFont();
   return EXIT_SUCCESS;
}

