/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.1.1
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
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
#include <math.h>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;
#include <ctime>

#include "../common/solscrutils.h"
#include "../common/solfileutils.h"
#include "../common/solmemhand.h"
#include "../common/iofuncts-wfx.h"
#include "../common/iofuncts-wfn.h"
#include "../common/wavefunctionclass.h"
#include "../common/bondnetwork.h"
#include "../common/wfgrid3d.h"
#include "optflags.h"
#include "crtflnms.h"



int main (int argc, char ** argv)
{
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam,logfilnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam,logfilnam); //This creates the names used.
   printHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << endl;
   
   gaussWaveFunc gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      setScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      setScrNormalFont();
      exit(1);
   }
   
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW(); //To setup the bond network.
   
   waveFunctionGrid3D grid;
   
   /* Looking for user grid dimensions */
   
   int nn=DEFAULTPOINTSPERDIRECTION;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&nn); //same number of points per direction
      grid.setNPts(nn);
   } else if (options.setn3) {
      int ny,nz;
      sscanf(argv[options.setn3],"%d",&nn);  //different number of points per direction
      sscanf(argv[options.setn3+1],"%d",&ny); //manually given by the user
      sscanf(argv[options.setn3+2],"%d",&nz);
      grid.setNPts(nn,ny,nz);
   } else if (options.setsmcub1) {
      sscanf(argv[options.setsmcub1],"%d",&nn); //uses nn for determining the number of points
                     //of the largest axis. The rest will have proportional number to their length
   } else {
      grid.setNPts(nn);
   }
   
   if ((options.setsmcub)||(options.setsmcub1)) {
      grid.setUpSmartCuboidGrid(gwf,bnw,nn);
   } else {
      grid.setUpSimpleGrid(gwf,bnw);
   }
   
   cout << "The size of the grid will be: " << grid.getNPts(0) << " x " 
        << grid.getNPts(1) << " x " << grid.getNPts(2) << endl;
   cout << "Total number of points that will be computed: " 
        << (grid.getNPts(0)*grid.getNPts(1)*grid.getNPts(2)) << endl;
   
   /* Setting the property to be computed */
   
   char prop;
   if (options.prop2plot) {
      prop=argv[options.prop2plot][0];
   } else {
      prop='d';
   }
   
   /* Main calculation loop, chooses between different available fields. */
   
   cout << "Evaluating and writing property..." << endl;
   cout << "(Scalar Field to plot: " << getFieldTypeKeyLong(prop) << ")." << endl << endl;
   switch (prop) {
      case 'd':
         grid.makeCube(outfilnam,gwf,DENS);
         cout << endl;
         break;
      case 'g':
         grid.makeCube(outfilnam,gwf,MGRD);
         break;
      case 'l':
         grid.makeCube(outfilnam,gwf,LAPD);
         cout << endl;
         break;
      case 'E':
         grid.makeCube(outfilnam,gwf,ELFD);
         break;
      case 'P' :
         grid.makeCube(outfilnam,gwf,MLED);
         break;
      case 'r' :
         grid.makeCube(outfilnam,gwf,ROSE);
         break;
      case 's' :
         grid.makeCube(outfilnam,gwf,REDG);
         break;
      case 'S':
         grid.makeCube(outfilnam,gwf,SENT);
         break;
      case 'L':
         grid.makeCube(outfilnam,gwf,LOLD);
         break;
      case 'M':
         grid.makeCube(outfilnam,gwf,MGLD);
         break;
      case 'G':
         grid.makeCube(outfilnam,gwf,KEDG);
         break;
      case 'K':
         grid.makeCube(outfilnam,gwf,KEDK);
         break;
      case 'u' :
         grid.makeCube(outfilnam,gwf,SCFD);
         break;
      case 'V':
         grid.makeCube(outfilnam,gwf,MEPD);
         break;
      default:
         setScrRedBoldFont();
         cout << "Error: The property \"" << prop << "\" does not exist!" << endl;
         setScrNormalFont();
         exit(1);
         break;
   }
   
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
      writeCommentedHappyStart(argv,lfil,CURRENTVERSION,PROGRAMCONTRIBUTORS);
      lfil << "#Wave function file name: " << endl << infilnam << endl;
      lfil << "#Number of primitives: "  << endl << gwf.nPri << endl;
      lfil << "#Grid dimensions:" << endl
           << grid.getNPts(0) << " " << grid.getNPts(1) << " " << grid.getNPts(2) << endl;
      lfil << "#Total number of points in the cube:" << endl
           << (grid.getNPts(0)*grid.getNPts(1)*grid.getNPts(2)) << endl;
      lfil << "#CPU Time (sec):" << endl;
      lfil << scientific << setprecision(4)
           <<  solreal( clock () - begin_time ) / CLOCKS_PER_SEC << endl;
      lfil << "#Wall-clock Time (sec):" << endl;
      solreal tmp_walltime=time(NULL);
      lfil << solreal (tmp_walltime-begin_walltime) << endl;
      lfil.close();
   }
   
   setScrGreenBoldFont();
   printHappyEnding();
   printScrStarLine();
   cout << setprecision(3) << "CPU Time: " << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   solreal end_walltime=time(NULL);
   cout << "Wall-clock time: " << solreal (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debugging mode (under construction)..." << endl;
#endif
   printScrStarLine();
   setScrNormalFont();
   return 0;
}

