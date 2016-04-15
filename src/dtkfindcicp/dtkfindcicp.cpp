/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
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
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "optflags.h"
#include "crtflnms.h"
#include "custfmtmathfuncts.h"

void writeCICPLogFile(string lnam,GaussWaveFunction &wf,critPtNetWork &cp);


int main (int argc, char ** argv)
{
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,logfilnam,cpxfilnam,cicpfilnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   ScalarFieldType critpttype=DENS;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   
   mkFileNames(argv,options,infilnam,logfilnam,cicpfilnam,cpxfilnam); //This creates the names used.
   printHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   GaussWaveFunction gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      setScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      setScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW();             //To setup the bond network.
   
   critPtNetWork cpn(gwf,bnw);
   
   cpn.setCriticalPoints(DENS);
   cpn.setBondPaths();
   cpn.writeCPProps(logfilnam,infilnam);
   writeCPXFile(cpxfilnam,infilnam,cpn);

   /* Searching of CPs finished.  */

   ofstream lfil;
   lfil.open(logfilnam.c_str(),std::ofstream::app);
   lfil << setprecision(3) << "CPU Time: " << endl
        << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   lfil.close();
   
   writeCICPLogFile(cicpfilnam,gwf,cpn);

   cout << endl << "Output written in files: " << logfilnam
        << ", " << cpxfilnam
        << ", and " << cicpfilnam << endl;
   
   
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


void writeCICPLogFile(string lnam,GaussWaveFunction &wf,critPtNetWork &cp) {
   ofstream lfil;
   lfil.open(lnam.c_str(),ios::out);
   int nbcp=cp.nBCP;
   int acp1Idx,acp2Idx;
   solreal RACP1[3],RACP2[3];
   writeCommentedScrStarLine(lfil);
   if ( nbcp==0 ) {
      lfil << "No bond critical points found!" << endl;
      lfil.close();
      return;
   }
   solreal gamma,gg[3],gp[3],hh[3][3],hph[3][3],hp[3][3];
   writeCommentedScrCharLine(lfil,'-');
   for ( int i=0 ; i<nbcp ; ++i ) {
      acp1Idx=cp.conBCP[i][0];
      acp2Idx=cp.conBCP[i][1];
      cout << cp.RBGP[i][0][0] << " " << cp.RBGP[i][0][1] << " " << cp.RBGP[i][0][2] << endl;
      for ( int j=0 ; j<3 ; ++j ) {
         RACP1[j]=cp.RACP[acp1Idx][j];
         RACP2[j]=cp.RACP[acp2Idx][j];
      }
      lfil << "#Bond Critical Point Index: " << endl << i << endl;
      wf.evalHessDensityMatrix1(RACP1,RACP2,gamma,gg,gp,hh,hph,hp);
      lfil << "#Hessian of Density Matrix of Order 1 at CI:" << endl
           << gamma << endl;
      writeCommentedScrCharLine(lfil,'-');
   }
   writeCommentedScrStarLine(lfil);
   lfil.close();
}


