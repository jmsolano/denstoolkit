/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.1.2
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
#include <sstream>
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
#include "../common/solmath.h"
#include "../common/wavefunctionclass.h"
#include "../common/bondnetwork.h"
#include "../common/wfgrid2d.h"
#include "optflags.h"
#include "crtflnms.h"

void makeGnuplotFile(optFlags &opts, string &gnpn,string &outn,char p2p,solreal dimparam,
                     bondNetWork &bn,int a1,int a2,int a3,waveFunctionGrid2D &grd);

int main (int argc, char ** argv)
{
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam,gnpnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam,gnpnam); //This creates the names used.
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
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW();             //To setup the bond network.
   
   waveFunctionGrid2D grid;    //Defining a grid object
   
   /* Looking for user grid dimensions */
   
   int nn=DEFAULTPOINTSPERDIRECTION;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&nn);
      grid.setNPts(nn);
   } else {
      grid.setNPts(nn);
   }
   
   /* Defining the line direction and other properties, and setting up the grid */
   
   int at1=0,at2=1,at3=2;
   if (bnw.nNuc==1) {
      at1=at2=0,at3=0;
      grid.setUpSimplePlane(bnw,0);
      cout << "The file " << infilnam << " has only one atom" << endl;
      cout << "Using this atom to set the plane...\n";
   } else if (bnw.nNuc==2) {
      at1=0;
      at2=at3=1;
      grid.setUpSimplePlane(bnw,0,1);
      cout << "The file " << infilnam << " has only two atoms" << endl;
      cout << "Using these atoms to set the plane...\n";
   } else {
      if (options.setats) {
         sscanf(argv[options.setats],"%d",&at1);
         sscanf(argv[options.setats+1],"%d",&at2);
         sscanf(argv[options.setats+2],"%d",&at3);
         if (at1<1||at2<1||at3<1||at1>bnw.nNuc||at2>bnw.nNuc||at3>bnw.nNuc) {
            setScrRedBoldFont();
            cout << "Error: one of the given atoms do not exist!\n";
            setScrNormalFont();
            exit(1);
         }
         cout << "Using atoms " << at1 << "(" << bnw.atLbl[at1-1]
         << "), " << at2 << "(" << bnw.atLbl[at2-1] << "), and "
         << at3 << "(" << bnw.atLbl[at3-1] << ") to set the plane."
         << endl;
         at1--;
         at2--;
         at3--;
         grid.setUpSimplePlane(bnw,at1,at2,at3);
      } else {
         at1=at2=at3=0;
         grid.setUpSimplePlane(bnw,0);
         cout << "Using the first atom to set the line...\n";
      }
   }
   //cout << "checkpoint" << endl;
   cout << "The size of the grid will be: " << grid.getNPts(0) << "x" << grid.getNPts(1) << endl;
   cout << "Total number of points that will be computed: " << grid.getNPts(0)*grid.getNPts(1) << endl;
   
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
         grid.makeTsv(outfilnam,gwf,DENS);
         cout << endl;
         break;
      case 'g':
         grid.makeTsv(outfilnam,gwf,MGRD);
         break;
      case 'l':
         grid.makeTsv(outfilnam,gwf,LAPD);
         cout << endl;
         break;
      case 'G':
         grid.makeTsv(outfilnam,gwf,KEDG);
         break;
      case 'E':
         grid.makeTsv(outfilnam,gwf,ELFD);
         break;
      case 'K':
         grid.makeTsv(outfilnam,gwf,KEDK);
         break;
      case 'S':
         grid.makeTsv(outfilnam,gwf,SENT);
         break;
      case 'L':
         grid.makeTsv(outfilnam,gwf,LOLD);
         break;
      case 'M':
         grid.makeTsv(outfilnam,gwf,MGLD);
         break;
      case 'N':
         grid.makeTsv(outfilnam,gwf,GLOL);
         break;
      case 'p' :
         grid.makeTsv(outfilnam,gwf,LEDV);
         break;
      case 'P' :
         grid.makeTsv(outfilnam,gwf,MLED);
         break;
      case 'r' :
         grid.makeTsv(outfilnam,gwf,ROSE);
         break;
      case 's' :
         grid.makeTsv(outfilnam,gwf,REDG);
         break;
      case 'u' :
         grid.makeTsv(outfilnam,gwf,SCFD);
         break;
      case 'U' :
         grid.makeTsv(outfilnam,gwf,VCFD);
         break;
      case 'V':
         grid.makeTsv(outfilnam,gwf,MEPD);
         break;
      default:
         setScrRedBoldFont();
         cout << "Error: The property \"" << prop << "\" does not exist!" << endl;
         setScrNormalFont();
         exit(1);
         break;
   }
   
   cout << endl << "Output written in file: " << outfilnam << endl;
   
#if _HAVE_GNUPLOT_
   cout << "          Gnuplot file: " << gnpnam << endl;
   solreal dd=0.0e0,r[3];
   for (int i=0; i<3; i++) {
      r[i]=grid.Ca[i]-grid.Cb[i];
      dd+=r[i]*r[i];
   }
   dd=0.45e0*sqrt(dd);
   if (options.mkplt) {
      makeGnuplotFile(options,gnpnam,outfilnam,prop,dd,bnw,at1,at2,at3,grid);
   }
   if (!(options.kpgnp)) {
      cout << "Cleaning temporary files..." << endl;
#if (defined(__APPLE__)||defined(__linux__))
      //string cmdl="rm ";
      //cmdl+=gnpnam;
      //system(cmdl.c_str());
      system("rm contourtemp.dat");
#endif//defined(__APPLE__)||defined(__linux__)
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

//*******************************************************************************************
void makeGnuplotFile(optFlags &opts, string &gnpn,string &outn,char p2p,
                     solreal dimparam,bondNetWork &bn,int a1,int a2,int a3,waveFunctionGrid2D &grd)
{
   ofstream gfil;
   gfil.open(gnpn.c_str());
   
   /* Choosing the label (legend) for the plot  and the zrange for the plot */
   solreal minzrange,maxzrange;
   string plbl=getFieldTypeKeyShort(p2p);;
   switch (p2p) {
      case 'd':
         minzrange=0.0e0;
         maxzrange=0.6e0;
         break;
      case 'g':
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'l':
         minzrange=-2.0e0;
         maxzrange=2.0e0;
         break;
      case 'E':
         minzrange=0.0e0;
         maxzrange=1.0e0;
         break;
      case 'S':
         minzrange=-10.0e0;
         maxzrange=10.0e0;
         break;
      case 'L':
         minzrange=0.0e0;
         maxzrange=1.0e0;
         break;
      case 'M':
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'N':
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'G':
         minzrange=0.0e0;
         maxzrange=10.0e0;
         break;
      case 'K':
         minzrange=0.0e0;
         maxzrange=10.0e0;
         break;
      case 'p' :
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'P' :
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'r' :
         minzrange=-1.0e0;
         maxzrange=1.0e0;
         break;
      case 's' :
         minzrange=0.0e0;
         maxzrange=1.0e0;
         break;
      case 'u' :
         minzrange=-1.0e0;
         maxzrange=1.0e0;
         break;
      case 'U' :
         minzrange=-1.0e0;
         maxzrange=1.0e0;
         break;
      case 'V':
         minzrange=-0.6e0;
         maxzrange=0.6e0;
         break;
      default:
         setScrRedBoldFont();
         cout << "Error: The property \"" << p2p << "\" does not exist!" << endl;
         setScrNormalFont();
         exit(1);
         break;
   }
   
   gfil << "reset" << endl;
   
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
   
   gfil << "set isosample 300, 300" << endl;
   
   gfil << "set zrange [" << minzrange << ":" << maxzrange << "]" << endl;
   
   gfil << "set contour base" << endl;

   gfil << "set cntrparam bspline" << endl;
   
   //gfil << "set cntrparam levels discrete 0.001, 0.005, 0.02, 0.05, 0.1, 0.2, 0.3" <<endl;

   gfil << "set cntrparam levels incremental " << minzrange << ", " 
        << (0.10*(maxzrange-minzrange)) << ", " << maxzrange << endl;
   
   gfil << "unset surface" << endl;
   
   gfil << "set table 'contourtemp.dat'" << endl;
   
   gfil << "splot '" << outn << "'";
   if (p2p=='N' || p2p=='p' || p2p=='U') {gfil << " using 1:2:(sqrt($3*$3+$4*$4))";}
   gfil << endl;
   
   gfil << "unset table" << endl;
   
   gfil << "reset" << endl;
   
   gfil << "unset key" << endl;
   
   gfil << "set size ratio 1" << endl;
   
   gfil << "set tmargin at screen 0.99\nset bmargin at screen 0.01\nset lmargin at screen 0.01" << endl;
   
   gfil << "dimparam=" << dimparam
        << " #Decrease this number to zoom in the plot" << endl;
   
   if (p2p=='N'||p2p=='p'||p2p=='U') {
      gfil << "VMXS=dimparam/40.0 #Maximum lenght of the vectors" << endl;
   }
   
   gfil << "set notics" << endl;
   
   gfil << "set xrange[-dimparam:dimparam]\nset yrange[-dimparam:dimparam]" << endl;
   
   gfil << "minvalfield=" << minzrange << endl << "maxvalfield=" << maxzrange << endl;
   
   gfil << "set cbrange [minvalfield:maxvalfield]" << endl;
   
   gfil << "set style data lines" << endl;
   
   gfil << "set palette rgbformulae 33,13,10" << endl;
   
   if (p2p=='N'||p2p=='p'||p2p=='U') {
      //gfil << "plot '" << outn << "' using 1:2:3:4:(sqrt($3*$3+$4*$4)) "
      //<< "with vectors head size 0.1,20,60 filled lc palette";
      gfil << "plot '" << outn << "' every 1:1 using 1:2:(sqrt($3*$3+$4*$4)>VMXS? VMXS*$3/sqrt($3*$3+$4*$4) : $3):(sqrt($3*$3+$4*$4)>VMXS ? VMXS*$4/sqrt($3*$3+$4*$4) : $4):(sqrt($3*$3+$4*$4)) "
           << "with vectors head size 0.1,20,60 filled lc palette" << endl;
   } else {
      gfil << "plot '" << outn << "' with image" << endl;
   }
   
   if (opts.showcont) {gfil << "replot 'contourtemp.dat' w l lt -1 lw 1 ";}
   else {gfil << " #replot 'contourtemp.dat' w l lt -1 lw 1 #Activate this line to show contours";}
   gfil << endl;
   
   /* Here the atoms' labels are set and written to the gnuplot input file. */
   
   writeScrCharLine(gfil,'#');
   gfil << "# Here are the labes of the atoms" << endl;
   writeScrCharLine(gfil,'#');
   
   string lbl;
   int at[3];
   at[0]=a1; at[1]=a2; at[2]=a3;
   solreal xproj,yproj;
   size_t pos;
   std::ostringstream numst;
   bool IDefPlane;
   for (int i=0; i<bn.nNuc; i++) {
      xproj=0.0e0;
      yproj=0.0e0;
      IDefPlane=false;
      for (int j=0; j<3; j++) {
         xproj+=(bn.R[i][j]-grd.orig[j])*grd.dircos1[j];
         yproj+=(bn.R[i][j]-grd.orig[j])*grd.dircos2[j];
         if ((i==at[j])&&opts.showatlbl) {IDefPlane=true;}
      }
      if (opts.showallatlbl||IDefPlane) {lbl="";} else {lbl="#";}
      lbl+="set label ";
      numst.str("");
      numst << (i+1);
      lbl+=numst.str();
      lbl+=" at ";
      numst.str("");
      numst << xproj << "," << yproj;
      lbl+=numst.str();
      lbl+=" '";
      lbl+=bn.atLbl[i];
      /* Adding an underscore and braces to the atome labels in order to make the numbers
         subindices */
      pos=lbl.find_last_not_of("0123456789");
      if (pos!=string::npos) {
         lbl.insert(pos+1,"_{");
         lbl.append("}");
      }
      lbl+="' front offset character -0.5,-0.15";
      for (int k=0; k<3; k++) {
         if (at[k]==i) {lbl+=" #* This atom was used for setting the plane!";}
      }
      gfil << lbl << endl;
   }
   writeScrCharLine(gfil,'#');
   
   /* Here is enabled the logarithmic scale in the case of G, d or g */
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   
   string imgnam,epsname,extepsname,pdfname;
   imgnam=outn.substr(0,(outn.length()-4));
   imgnam.append("ext.eps"); //here is set the name for the image
   gfil << "set output '";
   gfil << imgnam << "'" << endl; //so far, the image's name is *.eps
   gfil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   gfil << "replot" << endl;
   extepsname=imgnam;
   epsname=outn.substr(0,(outn.length()-3));
   epsname.append("eps");
   pdfname=outn.substr(0,(outn.length()-3));
   pdfname.append("eps");
   gfil << "#" << endl;
   writeScrCharLine(gfil,'#');
   gfil << "#                 END OF GNUPLOT COMMANDS" << endl;
   writeScrCharLine(gfil,'#');
   gfil << "#If you want to reconstruct the plot using this file, type:" << endl
        << "#gnuplot " << gnpn << endl
        << "#epstool --copy -b " << extepsname << " " << epsname << endl
        << "#epstopdf --outfile=" << pdfname << " " << extepsname << endl;
   
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
#if _HAVE_EPSTOOL_
   line=string("epstool --copy -b ");
   line+=(imgnam+string(" "));
   imgnam=outn.substr(0,(outn.length()-3)); //Here the image name changes from eps to pdf
   imgnam.append("eps");
   line+=imgnam;
   if (opts.quiet) {line+=string(" > /dev/null 2>&1");}
   system(line.c_str());
#endif
#if _HAVE_EPSTOPDF_
   string pdfnam=outn.substr(0,(outn.length()-3));
   pdfnam.append("pdf");
   line=string("epstopdf --outfile=")+pdfnam+string(" ")+imgnam;
   if (opts.quiet) {line+=string(" > /dev/null 2>&1");}
   system(line.c_str());
#endif
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   line="rm ";
#if _HAVE_EPSTOOL_
   line+=imgnam;
#endif
   imgnam=outn.substr(0,(outn.length()-4));
   imgnam.append("ext.eps");
   line+=(string(" ")+imgnam);
   system(line.c_str());
#endif//defined(__APPLE__)||defined(__linux__)
   cout << "Done." << endl;
   return;
}


