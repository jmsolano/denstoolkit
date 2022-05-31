/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
/* optsflags.cpp
   This file contains the implementation of the definitions of the class class optflags
   
   This file contains the implementations of the definitions of the class optflags, which is
   simply a class that will contain all the flags used during the execution of the
   simulation program dess3ds.
   
   It is intended to simplify the addition of new options and new variations of calculations
   the program dess3ds makes. 
   
   ------------------------

   Juan Manuel Solano Altamirano
   Adscription at the moment this project is initiated:
   Centro de Investigaciones y Estudios Avanzados del 
   Instituto Politecnico Nacional, 
   Unidad Monterrey, Mexico.
   2011
   e-mail: jmsolanoalt@gmail.com
 
   Adscription at the moment the particular implementation for this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
*/

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <cstdlib>
using std::exit;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;
#include "optflags.h"
#include "figname.h"
#include "../common/screenutils.h"
#include "../common/gausswavefunction.h"

OptionFlags::OptionFlags() {
   infname=0;
   outfname=0;
   setn1=setn3=0;
   setsmcub=setsmcub1=0;
   configspecialnci=0;
   setsisovalue=0;
   setcolorscalesingle=setcolorscaleboth=0;
   selectpalette=setgnpangles=setviewangles=orientcam3ats=0;
   rotX=rotY=rotZ=0;
   skipcube=kppov=mkpov=mkpng=rotcam=false;
}
void getOptions(int &argc, char** &argv, OptionFlags &flags) {
   string progname;
   progname=":-)  ";
   progname+=argv[0];
   progname+="  (-:";
   size_t pos;
   pos=progname.find("./");
   if (!(pos==string::npos)) {
      progname.erase(pos,2);
   }
   if (argc<2) {
      ScreenUtils::SetScrRedBoldFont();
      cout << "\nError: Not enough arguments." << endl;
      ScreenUtils::SetScrNormalFont();
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
      exit(1);
   }
   if (string(argv[1])==string("-h")) {
      printHelpMenu(argc,argv);
      exit(1);
   }
   for (int i=1; i<argc; i++){
      if (argv[i][0] == '-'){
         switch (argv[i][1]){
            case 'a' :
               if ((i+2)>=argc) {printErrorMsg(argv,'a');}
               flags.setgnpangles=(++i);
               i++;
               break;
            case 'A' :
               if ((i+3)>=argc) {printErrorMsg(argv,'A');}
               flags.setviewangles=(++i);
               i+=2;
               break;
            case 'b' :
               flags.setcolorscalesingle=(++i);
               if (i>=argc) {printErrorMsg(argv,'b');}
               break;
            case 'B' :
               flags.setcolorscaleboth=(++i);
               i++;
               if (i>=argc) {printErrorMsg(argv,'B');}
               break;
            case 'c' :
               flags.skipcube=true;
               break;
            case 'I' :
               flags.setsisovalue=(++i);
               if (i>=argc) {printErrorMsg(argv,'I');}
               break;
            case 'k' :
               flags.kppov=true;
               break;
            case 'l' :
               flags.selectpalette=(++i);
               if (i>=argc) {printErrorMsg(argv,'l');}
               break;
            case 'n' :
               if (i>=argc) {printErrorMsg(argv,'n');}
               flags.setn1=(++i);
               break;
            case 'N' :
               if ((i+3)>=argc) {printErrorMsg(argv,'N');}
               flags.setn3=(++i);
               i+=2;
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'O' :
               if ((i+3)>=argc) {printErrorMsg(argv,'O');}
               flags.orientcam3ats=(++i);
               flags.rotcam=true;
               break;
            case 'p' :
               flags.mkpov=true;
               break;
            case 'P' :
               flags.mkpov=true;
               flags.mkpng=true;
               break;
            case 's':
               flags.setsmcub=(i);
               break;
            case 'S':
               flags.setsmcub1=(++i);
               if (i>=argc) {printErrorMsg(argv,'S');}
               break;
            case 'V':
               progname=argv[0];
               pos=progname.find("./");
               if (!(pos==string::npos)) {progname.erase(pos,2);}
               cout << progname << " " << CURRENTVERSION << endl;
               exit(0);
               break;
            case 'x' :
               flags.rotX=(++i);
               flags.rotcam=true;
               break;
            case 'y' :
               flags.rotY=(++i);
               flags.rotcam=true;
               break;
            case 'z' :
               flags.rotZ=(++i);
               flags.rotcam=true;
               break;
            case 'h':
               printHelpMenu(argc,argv);
               exit(1);
               break;
            case '-':
               processDoubleDashOptions(argc,argv,flags,i);
               break;
            default:
               cout << "\nCommand line error. Unknown switch: " << argv[i] << endl;
               cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
               exit(1);
         }
      }
   }
   return;
}
void printHelpMenu(int &argc, char** &argv) {
   string progname=argv[0];
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}
   ScreenUtils::PrintScrStarLine();
#if _SOL_USE_FIGLET_NAME_
   FigletName::PrintFigletName();
#endif
   cout << endl;
   ScreenUtils::CenterString((string(":-) ")+progname+string(" (-:")));
   cout << endl;
   ScreenUtils::CenterString("This program calculate the necessary fields to plot the NCI plot.");
   ScreenUtils::CenterString("The main data is saved into a cub file, which contains the");
   ScreenUtils::CenterString("reduced density gradient (s). This data is used to locate an");
   ScreenUtils::CenterString("isosurface, upon which the field Lambda2=l*rho is mapped");
   ScreenUtils::CenterString("with a color scale. l is the sign of second hessian");
   ScreenUtils::CenterString("eigenvalue and rho the electron density.");
   ScreenUtils::CenterString("The isosurface is composed of triangles, thus the fields");
   ScreenUtils::CenterString("l and rho are computed at the centroid of each of");
   ScreenUtils::CenterString("the triangles that conform the isosurface.");
   cout << endl;
   ScreenUtils::CenterString((string("Compilation date: ")+string(__DATE__)));
   cout << endl;
   ScreenUtils::CenterString(string("Version: ")+string(CURRENTVERSION));
   cout << endl;
   ScreenUtils::CenterString((string(":-) Created by: ")+string(PROGRAMCONTRIBUTORS)+string(" (-:")));
   cout << endl;
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrBoldFont();
   cout << "\nUsage:\n\n\t" << progname << " wf?name [option [value(s)]] ... [option [value(s)]]\n\n";
   ScreenUtils::SetScrNormalFont();
   cout << "Where wf?name is the input wfx(wfn) name, and options can be:\n\n";
   cout << "  -a aG1 aG2\tSets the gnuplot angles to be aG1 and aG2.\n"
        << "            \t  Use dtkqdmol to check and set these angles." << '\n';
   cout << "  -A aX aY aZ\tSets the view angles to be aX, aY, and aZ." << '\n';
   cout << "  -b val    \tSets the color scale to be [-val,val]." << '\n';
   cout << "  -B vMin vMax\tSets the color scale to be [vMin,vMax]." << '\n';
   cout << "  -c        \tSkips the cube calculation. Notice this will assume that the\n"
        << "            \t  cube was previously computed." << '\n';
   cout << "  -I isoval \tSets the isovalue of s to be isoval" << '\n';
   cout << "  -k        \tKeeps the pov-ray script (see also option P, below)." << '\n';
   cout << "  -l palette\tSelects the color scheme 'palette', which can be any of:\n"
        << "            \t  bentcoolwarm blues bugn gnbu greens greys inferno\n"
        << "            \t  magma moreland oranges orrd plasma pubu purples rdbu\n"
        << "            \t  rdylbu rdylgn reds spectral viridis ylgn ylgnbu\n"
        << "            \t  ylorbr ylorrd" << '\n';
   cout << "  -n  dim   \tSets the number of points per direction for the s-cube" << endl
        << "            \t  to be dim x dim x dim." << endl;
   cout << "  -N nx ny nz\tSets the individual points per direction for the s-cube" << endl
        << "            \t  to be nx x ny x nz." << endl;
   cout << "  -o outname\tSets the output file names to use 'outname' as a basename,\n"
        << "            \t  i.e., the pov and png files will be named:\n"
        << "            \t  outnameNCI.pov and outnameNCI.png, respectively." << endl;
   cout << "  -O a b c  \tOrients the POV camera, so that the atoms a, b, and c (numbering\n"
        << "            \t  according to the wf? file) are placed over the screen. The\n"
        << "            \t  final view should look like Scheme 1, below." << '\n';
   cout << "  -s        \tUses a smart cuboid for the s-cube. The number of points for the" <<endl
        << "            \t  largest direction will be " << DEFAULTPOINTSPERDIRECTION << "." << endl;
   cout << "  -S ln     \tUses a smart cuboid for the s-cube. ln is the number of points" << endl
        << "            \t  the largest axis will have. The remaining axes will have" << endl
        << "            \t  a number of points proportional to the molecule dimensions." << endl;
   cout << "  -P        \tGenerates a pov-ray script and renders it. Notice: this requires" << endl
        << "            \t   povray to be installed in your system." << endl;
   cout << "  -x alpha  \tRotates the final view by alpha degrees around the x-axis." << '\n';
   cout << "  -y beta   \tRotates the final view by beta  degrees around the y-axis." << '\n';
   cout << "  -z gamma  \tRotates the final view by gamma degrees around the z-axis." << '\n';
   cout << "  -V        \tDisplays the version of this program." << endl;
   cout << "  -h\t\tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --configure-nci rMin rMax sMax \tSet the parameters rhoMin, rhoMax," << endl
        << "             \t\t  and redDensGradMax to be rMin, rMax, and sMax," << endl
        << "             \t\t  respectively. Default values: rhoMin=" << NCIRHOMIN << ",\n"
        << "             \t\t  rhoMax=" << NCIRHOMAX << ", and redGradMax=" <<  NCISMAX << ".\n";
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
   ScreenUtils::PrintScrCharLine('-');
   cout << "            \t           a\n"
        << "            \t           |\n"
        << "            \t          y|\n"
        << "            \t           |--------c\n"
        << "            \t           |________|______\n"
        << "            \t           / b      x       \n"
        << "            \t          /\n"
        << "            \t       z / \n"
        << "  Scheme 1: View of the aligned atoms (see option -O)." << '\n';
   ScreenUtils::PrintScrCharLine('-');
   //-------------------------------------------------------------------------------------
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'A' :
         cout << "should be followed by three real numbers." << '\n';
         break;
      case 'a' :
      case 'B' :
         cout << "should be followed by two real numbers." << endl;
         break;
      case 'b' :
      case 'I' :
         cout << "should be followed by a real number." << endl;
         break;
      case 'n':
      case 'S':
         cout << "should be followed by an integer." << endl;
         break;
      case 'N':
         cout << "should be followed by three integers." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      default:
         cout << "is triggering an unknown error." << endl;
         break;
   }
   ScreenUtils::SetScrNormalFont();
   cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
   cout << "\nto view the help menu.\n\n";
   exit(1);
   return;
}
void processDoubleDashOptions(int &argc,char** &argv,OptionFlags &flags,int &pos) {
   string progname=argv[0];
   size_t progpos=progname.find("./");
   if (progpos!=string::npos) {progname.erase(progpos,2);}
   string str=argv[pos];
   str.erase(0,2);
   if (str==string("version")) {
      cout << progname << " " << CURRENTVERSION << endl;
      exit(0);
   } else if (str==string("help")) {
      printHelpMenu(argc,argv);
      exit(0);
   } else if ( str==string("configure-nci") ) {
      flags.configspecialnci=(pos+1);
      pos+=3;
      if ((pos+1)>=(argc+3)) {
         ScreenUtils::DisplayErrorMessage(string("configure-nci must be followed by 3 real numbers!"));
         exit(1);
      }
   } else {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << endl;
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   return;
}

