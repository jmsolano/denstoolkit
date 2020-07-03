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
#include "../common/demat1critptnetworkbp.h"
#include "optflags.h"
#include "crtflnms.h"
#include "custfmtmathfuncts.h"

void WriteLogFile(string fname,DeMat1CriticalPointNetworkBP &dcpn,bondNetWork &bn,\
      string &wfnam);

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,logfilnam,cpxfilnam,cicpfilnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   //ScalarFieldType critpttype=DENS;
   
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

   int atom1=0,atom2=1;
   if ( options.atpos ) {
      atom1=std::stoi(string(argv[options.atpos]));
      atom2=std::stoi(string(argv[options.atpos+1]));
      --atom1; --atom2;
   }
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW();             //To setup the bond network.
   
   //critPtNetWork cpn(gwf,bnw);

   DeMat1CriticalPointNetworkBP dmcpbp(gwf,bnw);
   cout << "Computing CICP eigenvalues and signatures..." << endl;
   dmcpbp.ComputeCoreInteractionCPs2D();
   dmcpbp.ComputeCoreInteractionCPs6D();
   cout << "Output saved into file " << cicpfilnam << endl;
   WriteLogFile(cicpfilnam,dmcpbp,bnw,infilnam);

   cout << "Saving cps information to " << cpxfilnam << endl;
   writeCPXFile(cpxfilnam,infilnam,*(dmcpbp.cpn));

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

void WriteLogFile(string fname,DeMat1CriticalPointNetworkBP &dcpn,bondNetWork &bn,\
      string &wfnam) {
   ofstream ofil(fname.c_str());
   if ( !ofil.good() ) {
      displayErrorMessage(string("File ")+fname+string("could not be opened!"));
      ofil.close();
      return;
   }
   /* ************************************************************************** */
   writeCommentedScrStarLine(ofil);
   centerCommentedString(string("Data obtained from file: ")+wfnam,ofil);
   if ( dcpn.differentSignaturesCICPvsNN() ) {
      centerCommentedString("Found DIFFERENT SIGNATURES [CICP] vs [Nuc-Nuc]!",ofil);
   } else {
      centerCommentedString("Found SAME SIGNATURES [CICP] vs [Nuc-Nuc]!",ofil);
   }
   writeCommentedScrStarLine(ofil);
   centerCommentedString("CICP (ACP-ACP) eigenvalues 6D",ofil);
   writeCommentedScrStarLine(ofil);
   int n=dcpn.nCICP;
   for ( int i=0 ; i<n ; ++i ) {
      ofil << dcpn.cpn->lblBCP[i] << ": ";
      for ( int j=0 ; j<6 ; ++j ) {
         ofil << dcpn.eivalCICP6D[i][j] << " ";
      }
      ofil << endl;
   }
   /* ************************************************************************** */
   writeCommentedScrStarLine(ofil);
   centerCommentedString("CICP signatures 6D",ofil);
   writeCommentedScrStarLine(ofil);
   for ( int i=0 ; i<n ; ++i ) {
      ofil << dcpn.cpn->lblBCP[i] << ": ";
      ofil << dcpn.sigCICP6D[i] << " ";
      ofil << endl;
   }
   /* ************************************************************************** */
   writeCommentedScrStarLine(ofil);
   centerCommentedString("Nuc-Nuc eigenvalues 6D",ofil);
   writeCommentedScrStarLine(ofil);
   int a1,a2;
   for ( int i=0 ; i<n ; ++i ) {
      dcpn.cpn->findTwoClosestAtomsToBCP(i,a1,a2);
      ofil << bn.atLbl[a1] << '-' << bn.atLbl[a2] << ": ";
      for ( int j=0 ; j<6 ; ++j ) {
         ofil << dcpn.eivalNN6D[i][j] << " ";
      }
      ofil << endl;
   }
   /* ************************************************************************** */
   writeCommentedScrStarLine(ofil);
   centerCommentedString("Nuc-Nuc signatures 6D",ofil);
   writeCommentedScrStarLine(ofil);
   for ( int i=0 ; i<n ; ++i ) {
      dcpn.cpn->findTwoClosestAtomsToBCP(i,a1,a2);
      ofil << bn.atLbl[a1] << '-' << bn.atLbl[a2] << ": ";
      ofil << dcpn.sigCICP6D[i] << endl;
   }
   /* ************************************************************************** */
   writeCommentedScrStarLine(ofil);
   centerCommentedString("CICP (ACP-ACP) eigenvalues 2D",ofil);
   writeCommentedScrStarLine(ofil);
   n=dcpn.nCICP;
   for ( int i=0 ; i<n ; ++i ) {
      ofil << dcpn.cpn->lblBCP[i] << ": ";
      for ( int j=0 ; j<2 ; ++j ) {
         ofil << dcpn.eivalCICP2D[i][j] << " ";
      }
      ofil << endl;
   }
   /* ************************************************************************** */
   writeCommentedScrStarLine(ofil);
   centerCommentedString("CICP signatures 2D",ofil);
   writeCommentedScrStarLine(ofil);
   for ( int i=0 ; i<n ; ++i ) {
      ofil << dcpn.cpn->lblBCP[i] << ": ";
      ofil << dcpn.sigCICP2D[i] << " ";
      ofil << endl;
   }
   writeCommentedScrStarLine(ofil);
   ofil.close();
}


