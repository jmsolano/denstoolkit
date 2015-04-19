/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.1
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
 ----------------------
 */

#ifndef _HAVE_GNUPLOT_
#define _HAVE_GNUPLOT_ 0
#endif

#ifndef _HAVE_GRAPHICSMAGICK_
#define _HAVE_GRAPHICSMAGICK_ 0
#endif

#ifndef _HAVE_IMAGEMAGICK_
#define _HAVE_IMAGEMAGICK_ 0
#endif

#ifndef _HAVE_POVRAY_
#define _HAVE_POVRAY_ 0
#endif

#ifndef CHOOSEPOVVERSION36
#define CHOOSEPOVVERSION36 1
#endif

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
#include "../common/iofuncts-cpx.h"
#include "../common/solmath.h"
#include "../common/wavefunctionclass.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "optflags.h"
#include "crtflnms.h"
#include "custfmtmathfuncts.h"

int main (int argc, char ** argv)
{
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam,povfilnam,pngfilnam,cpxfilnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   ScalarFieldType critpttype=DENS;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   
   if (options.cptype) {
      char ccpt;
      ccpt=argv[options.cptype][0];
      switch (ccpt) {
         case 'd':
            critpttype=DENS;
            break;
         case 'L':
            critpttype=LOLD;
            options.calcbgps=0;
            break;
         default:
            displayErrorMessage("This type of field is not implemented/known...");
            exit(1);
            break;
      }
   }
   
   mkFileNames(argv,options,infilnam,outfilnam,povfilnam,
               pngfilnam,cpxfilnam,critpttype); //This creates the names used.
   printHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   gaussWaveFunc gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      setScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      setScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   cout << "nNuc: " << gwf.nNuc << " nPri: " << gwf.nPri
        << " nMOr: "  << gwf.nMOr << endl;
   
   //return 0;
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW();             //To setup the bond network.
   
   critPtNetWork cpn(gwf,bnw);
   
   switch (critpttype) {
      case DENS:
         cpn.setCriticalPoints(DENS);
         if (options.calcbgps) {cpn.setBondPaths();}
         break;
      case LOLD:
         cpn.setCriticalPoints(LOLD);
         break;
      default:
         break;
   }
   if ( options.mkextsearch ) {
      cpn.extendedSearchCPs();
   }
   if (options.calcrgps) {
      cpn.setRingPaths();
      cpn.setCagePaths();
   }

   //cpn.displayIHVCoords();
   //cpn.displayACPCoords();
   //cpn.displayBCPCoords();
   //cpn.printCPProps(gwf);
   
   cpn.writeCPProps(outfilnam,infilnam);
   ofstream lfil;
   lfil.open(outfilnam.c_str(),std::ofstream::app);
   lfil << setprecision(3) << "CPU Time: " << endl
        << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   lfil.close();
   writeCPXFile(cpxfilnam,infilnam,cpn);
   
   cout << endl << "Output written in files: " << outfilnam << ", and " << cpxfilnam << endl;
   
#if _HAVE_POVRAY_
   int cameravdir=1;
   if (options.camvdir) {sscanf(argv[options.camvdir],"%d",&cameravdir);}
   if (options.drawnuc) {
      cpn.drawNuclei(true);
   }
   string cmdl;
   if (options.mkpov||options.mkpng) {
      povRayConfProp povconf;
      if (options.drawbgps&&(!options.calcbgps)) {
         displayWarningMessage("If you want to see gradient paths, you must not use option -G");
         displayWarningMessage("Nothing to include in the pov file.");
      }
      if (options.drawbgps&&options.calcbgps) {
         cpn.drawBondGradPaths(true);
         cpn.drawBonds(false);
         if (options.bgptubes) {cpn.tubeStyleBGP(true);}
      }
      cpn.makePOVFile(povfilnam,povconf,cameravdir);
   }
   if (options.kppov) {
      cout << "           PovRay file: " << povfilnam << endl;
   }
   if (options.mkpng) {
      cout << "Calling povray..." << endl;
      cmdl=string(CMD_POVRAY);
#if (defined(__CYGWIN__))
      cmdl+=string(" /EXIT /RENDER");
#endif
      cmdl+=(string(" ")+povfilnam+string(" +ua -D +FN"));
      //cout << cmdl << endl;
      if (options.quiet) {cmdl+=string(" > /dev/null 2>&1");}
      system(cmdl.c_str());
#if (_HAVE_IMAGEMAGICK_)
      cmdl="convert ";
#elif(_HAVE_GRAPHICSMAGICK_)
      cmdl="gm convert ";
#endif
#if (((_HAVE_IMAGEMAGICK_))||(_HAVE_GRAPHICSMAGICK_))
      cmdl+=(string("-trim ")+pngfilnam+string(" ")+pngfilnam);
      system(cmdl.c_str());
#endif
      cout << "Rendering done." << endl;
   }
   if (!(options.kppov)) {
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
      cmdl="rm ";
      cmdl+=povfilnam;
      system(cmdl.c_str());
#endif
   }
#endif /* _HAVE_POVRAY_ */
   
   if (options.mkdatmat) {
      string acfname,cpfname,bpfname;
      mkDatMatFileNames(outfilnam,acfname,cpfname,bpfname);
      writeDatMatAtCrds(acfname,bnw);
      writeDatMatCritPtsCrds(cpfname,cpn);
      writeDatMatBondPathCrds(bpfname,cpn);
   }
   
   
//#if (defined(__APPLE__)||defined(__linux__))
//   if (options.zipdat) {
//      string cmdl;
//      cmdl=string("gzip -9f ")+outfilnam;
//      cout << "Calling gzip...";
//      system(cmdl.c_str());
//      cout << " Done!" << endl;
//   }
//#endif/* defined(__APPLE__)||defined(__linux__) */
   
   
   /* At this point the computation has ended. Usually this means no errors ocurred. */
   
   setScrGreenBoldFont();
   printHappyEnding();
   printScrStarLine();
   cout << setprecision(3) << "CPU Time: "
        << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   solreal end_walltime=time(NULL);
   cout << "Wall-clock time: " << solreal (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debuggin mode (under construction...)" << endl;
#endif
   printScrStarLine();
   setScrNormalFont();
   return 0;
}

