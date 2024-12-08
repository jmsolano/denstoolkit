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
using std::endl;
using std::ios;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <sstream>
#include <cstdlib>
using std::exit;
#include <cmath>
#include <string>
#include <iomanip>
using std::setprecision;
#include <ctime>

#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/stringtools.h"
#include "../common/mymemory.h"
#include "../common/iofuncts-wfx.h"
#include "../common/iofuncts-wfn.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/atomcolschjmol.h"
#include "optflags.h"
#include "crtflnms.h"

void makeMolGnuplotFile(string &gname,BondNetWork &bn,bool putHs=true,const string term = "x11");

int main (int argc, char ** argv)
{
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,gnpnam;
   string progname;
   OptionFlags options;
   ifstream ifile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,gnpnam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   /*
   GaussWaveFunction gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   // */
   
   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.SetUpBNW();             //To setup the bond network.
   
   cout << " Done." << endl;
   
   string userterm="x11";
   if (options.settermgnp) {
      userterm=string(argv[options.settermgnp]);
   }
   bool drawHydrogens=true;
   if ( !(options.drawhydrogens) ) {
      drawHydrogens=false;
   }
   
   makeMolGnuplotFile(gnpnam,bnw,drawHydrogens,userterm);
   
   cout << endl << "Tip: to see the molecule type on your command line:\n" << endl;
#if (defined(__APPLE__)||defined(__linux__))
   string line="gnuplot ";
#endif
#if (defined(__CYGWIN__))
   string line="wgnuplot ";
#endif
   line+=gnpnam;
   ScreenUtils::SetScrBlueBoldFont();
   cout << line;
   ScreenUtils::SetScrNormalFont();
   cout << endl << endl;
   
   if (options.rungnp) {
      cout << "running the gnuplot file..." << endl;
      system(line.c_str());
   }
   
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
//*******************************************************************************************

void makeMolGnuplotFile(string &gname,BondNetWork &bn,bool putHs,const string term)
{
   ofstream gfil;
   gfil.open(gname.c_str(),ios::out);
   string tit=bn.title[0];
   StringTools::RemoveSpacesLeftAndRight(tit);
   gfil << "set title \"" << StringTools::GetEnhancedEpsTitle(tit) << "\\nPress any key on this window to quit...\"" << endl;
   
   gfil << "set hidden3d" << endl;
   gfil << "set border 4095" << endl;
   gfil << "unset xtics" << endl;
   gfil << "unset ytics" << endl;
   gfil << "unset ztics" << endl;
   
   gfil << "xshift=0" << endl;
   gfil << "yshift=0" << endl;
   gfil << "zshift=0" << endl;
   gfil << "set view equal xyz" << endl;
   gfil << "set xyplane relative 0" << endl;
#if (defined(__APPLE__)||defined(__linux__))
   gfil << "set term " << term << " enhanced font \"Sans,14\"" <<endl;
#endif
#if (defined(__CYGWIN__))
   gfil << "set term windows enhanced font \"Sans,14\"" <<endl;
#endif
   
   FileUtils::WriteScrCharLine(gfil,'#');
   gfil << "# Here are the labes of the atoms" << endl;
   FileUtils::WriteScrCharLine(gfil,'#');
   
   //gfil << "set label 1 at 0,0,0 'Press any key on this window to quit...' front" << endl;
   
   int effIdx=0;
   for (int i=0; i<bn.nNuc; ++i) {
      if ( (!putHs)&&(bn.atNum[effIdx]==0) ) { continue; }
      tit="set label ";
      tit+=(StringTools::GetStringFromInt(effIdx+1)+string(" \""));
      tit+=(StringTools::GetEnhancedEpsAtLbl(bn.atLbl[effIdx])+string("\" at "));
      tit+=(StringTools::GetStringFromReal(bn.R[effIdx][0])+string(","));
      tit+=(StringTools::GetStringFromReal(bn.R[effIdx][1])+string(","));
      tit+=(StringTools::GetStringFromReal(bn.R[effIdx][2])+string(" "));
      tit+=string("front offset character 0.75,0");
      gfil << tit << endl;
      ++effIdx;
   }
   FileUtils::WriteScrCharLine(gfil,'#');
   
   gfil << "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" << endl;
   gfil << "splot \"-\" using 1:2:3:(rgb($4,$5,$6)) with points "
        << "lc rgb variable lt 3 pt 7 ps 2 notitle, \\" << endl;
   gfil << "      \"-\" using 1:2:3:($4-$1):($5-$2):($6-$3) with vectors nohead"
        << " lc rgb 'black' lw 3 lt 1 notitle" << endl;
   int atnum;
   for (int i=0; i<bn.nNuc; i++) {
      atnum=bn.atNum[i];
      if ( (!putHs)&&(atnum==0) ) { continue; }
      gfil << bn.R[i][0] << " " << bn.R[i][1] << " " << bn.R[i][2] << " ";
      if (atnum==0) {
         gfil << "200 200 200" << endl;
      } else {
         gfil << GetAtomicRColorInt(atnum) << " " << GetAtomicGColorInt(atnum) << " "
              << GetAtomicBColorInt(atnum) << endl;
      }
   }
   gfil << "e" << endl;
   int k;
   for (int i=0; i<bn.nNuc; i++) {
      if ( (!putHs)&&(bn.atNum[i]==0) ) { continue; }
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bn.bNet[i][j];
         if ( (!putHs)&&(bn.atNum[k]==0) ) { continue; }
         if (k>0) {
            gfil << bn.R[i][0] << " " << bn.R[i][1] << " " << bn.R[i][2] << " ";
            gfil << bn.R[k][0] << " " << bn.R[k][1] << " " << bn.R[k][2] << endl;
         }
         
      }
   }
   gfil << "e" << endl;
#if (defined(__APPLE__)||defined(__linux__))
   gfil << "pause mouse keypress,button3,button2 ''" << endl;
#endif
#if (defined(__CYGWIN__))
   gfil << "pause mouse keypress ''" << endl;
#endif
   gfil.close();
   return;
}
