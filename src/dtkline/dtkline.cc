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
#include <ctime>

#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymemory.h"
#include "../common/iofuncts-wfx.h"
#include "../common/iofuncts-wfn.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/wfgrid1d.h"
#include "optflags.h"
#include "crtflnms.h"

#ifndef _HAVE_GNUPLOT_
#define _HAVE_GNUPLOT_ 0
#endif

void MakeGnuplotFile(string &gnpn,string &outn,char p2p,double dist,string &l1, string &l2);

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
   if ( options.stpspindens && gwf.ihaveSingleSpinOrbs ) {
      gwf.CalcCabAAndCabB();
   }
   
   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.SetUpBNW();             //To setup the bond network.
   
   WaveFunctionGrid1D grid;    //Defining a grid object
   
   /* Looking for user grid dimensions */
   
   int nn=DEFAULTPOINTSPERDIRECTION;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&nn);
      grid.SetNPts(nn);
   } else {
      grid.SetNPts(nn);
   }
   
   /* Defining the line direction and other properties, and setting up the grid */
   
   int at1=0,at2=1;
   if (bnw.nNuc==1) {
      at1=at2=0;
      cout << "The file " << infilnam << " has only one atom" << endl;
      cout << "Using this atom to set the line...\n";
      grid.SetUpSimpleLine(bnw,0);
   } else if (bnw.nNuc==2) {
      at1=0;
      at2=1;
      cout << "The file " << infilnam << " has only two atoms" << endl;
      cout << "Using these atoms to set the line...\n";
      grid.SetUpSimpleLine(bnw,0,1);
   } else {
      if (options.setats) {
         sscanf(argv[options.setats],"%d",&at1);
         sscanf(argv[options.setats+1],"%d",&at2);
         if (at1<1||at2<1) {
            ScreenUtils::SetScrRedBoldFont();
            cout << "Error: one of the given atoms do not exist!\n";
            ScreenUtils::SetScrNormalFont();
            exit(1);
         }
         cout << "Using atoms " << at1 << "(" << bnw.atLbl[at1-1]
              << ") and " << at2 << "(" << bnw.atLbl[at2-1] << ")" << " to set the line."
              << endl;
         at1--;
         at2--;
         grid.SetUpSimpleLine(bnw,at1,at2);
      } else {
         at1=at2=0;
         grid.SetUpSimpleLine(bnw,0);
         cout << "Using the first atom to set the line...\n";
      }
   }
   
   cout << "The size of the grid will be: " << grid.GetNPts() << endl;
   cout << "Total number of points that will be computed: " << grid.GetNPts() << endl;
   
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
   grid.MakeDat(outfilnam,gwf,prop);
   cout << endl << "Output written in file: " << outfilnam << endl;
   
#if _HAVE_GNUPLOT_
   cout << "          Gnuplot file: " << gnpnam << endl;
#endif /* _HAVE_GNUPLOT_ */
   
#if _HAVE_GNUPLOT_
   double dd=0.0e0,r[3];
   for (int i=0; i<3; i++) {
      r[i]=bnw.R[at2][i]-bnw.R[at1][i];
      dd+=r[i]*r[i];
   }
   dd=sqrt(dd);
   if (options.mkplt) {
      MakeGnuplotFile(gnpnam,outfilnam,prop,0.5*dd,bnw.atLbl[at1],bnw.atLbl[at2]);
   }
   if (!(options.kpgnp)) {
      cout << "Cleaning temporary files..." << endl;
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
      string cmdl="rm ";
      cmdl+=gnpnam;
      system(cmdl.c_str());
#endif//defined(__APPLE__)||defined(__linux__)
   }
#endif /* _HAVE_GNUPLOT_ */
   
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   if (options.zipdat) {
      string cmdl;
      cmdl=string("gzip -9f ")+outfilnam;
      cout << "Calling gzip...";
      system(cmdl.c_str());
      cout << " Done!" << endl;
   }
#endif/* defined(__APPLE__)||defined(__linux__) */
   
   
   /* At this point the computation has ended. Usually this means no errors ocurred. */
   
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
void MakeGnuplotFile(string &gnpn,string &outn,char p2p,double dist,string &l1,string &l2) {
   ofstream gfil;
   gfil.open(gnpn.c_str());
   
   /* Choosing the label (legend) for the plot */
   
   string plbl="";
   plbl=GnuplotFieldTitle(p2p);
   
   /* In this part the name is scanned for possible occurrings of the character '_'.
      For a proper display in the eps file, it has to be changed to "\_" */
   string line="";
   for (size_t i=0; i<outn.length(); i++) {
      if (outn[i]=='_') {line+='\\';}
      line+=outn[i];
   }
   gfil << "set title '" << line << "'" << endl;
   
   /* Adding an underscore and braces to the atome labels in order to make the numbers
      subindices */
   size_t pos;
   pos=l1.find_last_not_of("0123456789");
   line=l1;
   if (pos!=string::npos) {
      line.insert(pos+1,"_{");
      line.append("}");
   }
   gfil << "set xtics add ('" << line << "' " << (-dist) << ")" << endl;
   pos=l2.find_last_not_of("0123456789");
   line=l2;
   if (pos!=string::npos) {
      line.insert(pos+1,"_{");
      line.append("}");
   }
   gfil << "set xtics add ('" << line << "' " << (dist) << ")" << endl;
   
   /* Here is enabled the logarithmic scale in the case of G, d or g */
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   if (p2p=='G'||p2p=='d'||p2p=='g'||p2p=='r') {
      gfil << "set logscale y" << endl << "set format y \"10^{%T}\"" << endl;
      gfil << "set yrange [0.001:1000]" << endl;
   } else {
      gfil << "#set logscale y" << endl << "#set format y \"10^{%T}\"" << endl;
      gfil << "#set yrange [0.001:1000]" << endl;
   }
   
   string imgnam;
   imgnam=outn.substr(0,(outn.length()-3));
   imgnam.append("eps"); //here is set the name for the image
   string epsnam=imgnam;
   gfil << "#set output '";
   gfil << imgnam << "'" << endl; //so far, the image's name is *.eps
   imgnam=outn.substr(0,(outn.length()-3)); //Here the image name changes from eps to pdf
   imgnam.append("pdf");
   string pdfnam=imgnam;
#if (defined(__APPLE__)||defined(__linux__))
   gfil << "set output '| epstopdf --filter --outfile=";
   gfil << imgnam << "'" << endl;
#endif
#if (defined(__CYGWIN__))
   gfil << "set output '" << epsnam << "'" << endl;
#endif
   //gfil << "replot" << endl; 
   gfil << "plot '" << outn << "' with lines ls 1 lw 2 title '" << plbl << "'" << endl;
   
   gfil.close();
#if (defined(__APPLE__)||defined(__linux__))
   line=string("gnuplot ");
#endif
#if (defined(__CYGWIN__))
   line=string("wgnuplot ");
#endif
   line+=gnpn;
   cout << "Calling gnuplot... ";
   system(line.c_str());
   cout << "Done." << endl;
#if (defined(__CYGWIN__))
#if _HAVE_EPSTOPDF_
   line="epstopdf --outfile="+pdfnam+" "+epsnam;
   system(line.c_str());
   cout << "Cleaning temporary files" << endl;
   line="rm -f "+epsnam;
   system(line.c_str());
#endif
#endif
   return;
}

