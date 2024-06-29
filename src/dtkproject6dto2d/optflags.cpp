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
#include "figname.h"
#include "optflags.h"
#include "fldtypesdef.h"
#include "../common/screenutils.h"

OptionFlags::OptionFlags() {
   infname=0;
   outfname=0;
   setn1=0;
   setats=0;
   setstep=0;
   uponbp=1;
   uponsl=0;
   prop2plot=0;
   zipdat=0;
   mkplt=0;
   kpgnp=1;
   quiet=1;
   showcont=0;
   showatlbls=0;
   setinccont=findcps=0;
   centredats=false;
   stpspindens=false;
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
      cout << "\nError: Not enough arguments." << '\n';
      ScreenUtils::SetScrNormalFont();
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << '\n' << "to view the help menu.\n\n";
      exit(1);
   }
   if (string(argv[1])==string("-h")) {
      printHelpMenu(argc,argv);
      exit(1);
   }
   for (int i=1; i<argc; i++){
      if (argv[i][0] == '-'){
         switch (argv[i][1]){
            case 'a':
               if ((i+2)>=argc) {printErrorMsg(argv,'a');}
               flags.setats=(++i);
               i++;
               break;
            case 'c':
               flags.showcont=i;
               break;
            case 'C' :
               if ((i+3)>=argc) {printErrorMsg(argv,'C');}
               flags.setinccont=(++i);
               i+=2;
               break;
            case 'J' :
               flags.stpspindens=true;
               break;
            case 'k':
               flags.kpgnp=i;
               break;
            case 'l':
               flags.showatlbls=i;
               break;
            case 'L':
               flags.uponbp=0;
               flags.uponsl=i;
               break;
            case 'n':
               flags.setn1=(++i);
               if (i>=argc) {printErrorMsg(argv,'n');}
               break;
            case 'p' :
               flags.prop2plot=(++i);
               if ( i>=argc ) {printErrorMsg(argv,'p');}
               break;
            case 's':
               flags.setstep=(++i);
               if (i>=argc) {printErrorMsg(argv,'s');}
               break;
            case 't' :
               flags.centredats=true;
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'P':
               flags.mkplt=i;
               break;
            case 'T' :
               flags.findcps=i;
               break;
            case 'z':
               flags.zipdat=i;
               break;
            case 'v':
               flags.quiet=0;
               break;
            case 'V':
               progname=argv[0];
               pos=progname.find("./");
               if (!(pos==string::npos)) {progname.erase(pos,2);}
               cout << progname << " " << CURRENTVERSION << endl;
               exit(0);
               break;
            case 'h':
               printHelpMenu(argc,argv);
               exit(1);
               break;
            case '-':
               processDoubleDashOptions(argc,argv,flags,i);
               break;
            default:
               cout << "\nCommand line error. Unknown switch: " << argv[i] << '\n';
               cout << "\nTry: \n\t" << argv[0] << " -h\n" << '\n' << "to view the help menu.\n\n";
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
   cout << '\n';
   ScreenUtils::CenterString((string(":-) ")+progname+string(" (-:")));
   cout << '\n';
   ScreenUtils::CenterString("This program computes a 6D field (a two-3D-coordinate dependent function)");
   ScreenUtils::CenterString("along a parametrized bond path between two atoms (optionally upon");
   ScreenUtils::CenterString("the straight line that joins those atoms).");
   ScreenUtils::CenterString("The output consists of a tsv file.");
   ScreenUtils::CenterString("The information for the calculation is obtained from a wfx(wfn) file,");
   ScreenUtils::CenterString("which is given as the input for the program.");
   ScreenUtils::CenterString("In this version, the program can project the field DM1,");
   ScreenUtils::CenterString("i.e., the Density Matrix of order 1, its gradient");
   ScreenUtils::CenterString("and its Laplacian. In future releases,");
   ScreenUtils::CenterString("the program will be able to plot other 6D fields.");
   ScreenUtils::CenterString("In addition, in the current version, the program can determine");
   ScreenUtils::CenterString("the position (over the 2D projected domain) of the");
   ScreenUtils::CenterString("critical points of DM1 (maxima, saddles, and minima).");
   ScreenUtils::CenterString("The CPs are determined upon request.");
   ScreenUtils::CenterString("(See below for the sintax.)");
   cout << '\n';
   ScreenUtils::CenterString((string("Compilation date: ")+string(__DATE__)));
   cout << '\n';
   ScreenUtils::CenterString(string("Version: ")+string(CURRENTVERSION));
   cout << '\n';
   ScreenUtils::CenterString((string(":-) Created by: ")+string(PROGRAMCONTRIBUTORS)+string(" (-:")));
   cout << '\n';
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrBoldFont();
   cout << "\nUsage:\n\n\t" << progname << " wf?name [option [value(s)]] ... [option [value(s)]]\n\n";
   ScreenUtils::SetScrNormalFont();
   cout << "Where wf?name is the input wfx(wfn) name, and options can be:\n\n"
        << "  -a a1 a2  \tDefine the atoms  (a1,a2) used to define bond path/line." << '\n'
        << "            \t  If this option is not activated, the program will " << '\n'
        << "            \t  set a1=1, a2=2." << '\n'
        << "            \t  Note: if the *.wfn (*.wfx) file has only one atom" << '\n'
        << "            \t  the program will exit and no output will be generated." << '\n';
#if _HAVE_GNUPLOT_
   cout << '\n';
   cout << "  -P        \tCreate a plot using gnuplot." << '\n'
        << "  -c        \tShow contour lines in the plot." << '\n';
   cout << "  -C s i e  \tSet the contour values (incremental style)." << '\n'
        << "            \t  s, i, and e are real numbers. s is the first" << '\n'
        << "            \t  contour value, i is the increment, and " << '\n'
        << "            \t  e is the last contour value." << '\n';
   cout << "  -J        \tSetup alpha- and beta-spin density matrices." << '\n';
   cout << "  -k     \tKeeps the *.gnp files to be used later by gnuplot." << '\n';
   cout << "  -l     \tShow labels of atoms (those set in option -a) in the plot." << '\n'
        << '\n';
#endif
   cout << "  -L        \tCalculate MD1 upon the straight line that joins the atoms" << '\n'
        << "            \t  instead of upon the bond path." << '\n';
   cout << "  -n  dim   \tSet the number of points for the tsv file per direction." << '\n'
        << "            \t  Note: for the bond path you may want to look for a good " << '\n'
        << "            \t  combination of n and the number \"step\" given in option -s," << '\n'
        << "            \t  since the number of points in the bond path will be mainly " << '\n'
        << "            \t  governed by step." << '\n';
   cout << "  -p prop   \tChoose the property to be computed. prop is a character," << '\n'
        << "            \t  which can be (D is the default value):" << '\n';
   string charFields="gnl";
   for ( size_t i=0 ; i<charFields.size() ; ++i ) {
      cout << "         \t\t" << charFields[i] << " ("
         << GetField6DTypeKeyLong(charFields[i]) << ')' << '\n';
   }
   cout << "            \t      D (Density Matrix of order 1)" << '\n'
        << "            \t      G (Gradient of the density Matrix of order 1)" << '\n'
        << "            \t      L (Laplacian of the density Matrix of order1)" << '\n';
   cout << "  -o outname\tSet the output file name." << '\n'
        << "            \t  (If not given the program will create one out of" << '\n'
        << "            \t  the input name; if given, the tsv, gnp and pdf files will" << '\n'
        << "            \t  use this name as well --but different extension--)." << '\n';
   cout << "  -s step   \tSet the stepsize for the bond path to be 'step'." << '\n'
        << "            \t  Default value: " << DEFAULTBONDPATHSTEPMD1 << '\n';
   cout << "  -t        \tTranslate the plot to the geometrical center of the" << '\n'
        << "            \t  requested atoms." << '\n';
   cout << "  -T        \tPerform the Topological analysis (find critical points)." << '\n';
   cout << '\n';
   cout << "  -v     \tVerbose (display extra information, usually output from third-" << '\n'
        << "         \t  party sofware such as gnuplot, etc.)" << '\n';
   //cout << "  -V     \tShows the current version of this program." << '\n';
#if (defined(__APPLE__)||defined(__linux__))
   cout << "  -z     \tCompress the tsv file using gzip (which must be intalled" << '\n'
        << "         \t   in your system)." << '\n';
#endif
   cout << '\n';
   cout << "  -h     \tDisplay the help menu.\n";
   cout << "  -V     \tDisplays the version of this program." << '\n';
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << '\n';
   cout << "  --version \t\tSame as -V" << '\n';
   cout << '\n';
   //-------------------------------------------------------------------------------------
#if _HAVE_GNUPLOT_
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::CenterString(string("Note that the following programs must be properly installed in your system:"));
   ScreenUtils::CenterString(string("gnuplot"));
   ScreenUtils::CenterString(string("epstool"));
   ScreenUtils::CenterString(string("epstopdf"));
#if (defined(__APPLE__)||defined(__linux__))
   ScreenUtils::CenterString(string("gzip"));
#endif
   ScreenUtils::PrintScrStarLine();
#endif
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'a':
         cout << "should be followed by two integers." << '\n';
         break;
      case 'C':
         cout << "should be followed by three real numbers." << '\n';
         break;
      case 'n':
         cout << "should be followed by an integer." << '\n';
         break;
      case 'p' :
         cout << "should be followed by a character." << '\n';
         break;
      case 's':
         cout << "should be followed by a real number." << '\n';
         break;
      case 'o':
         cout << "should be followed by a name." << '\n';
         break;
      default:
         cout << "is triggering an unknown error." << '\n';
         break;
   }
   ScreenUtils::SetScrNormalFont();
   cout << "\nTry:\n\t" << argv[0] << " -h " << '\n';
   cout << "\nto view the help menu.\n\n";
   exit(1);
   return;
}
void processDoubleDashOptions(int &argc,char** &argv,OptionFlags &flags,int pos) {
   string progname=argv[0];
   size_t progpos=progname.find("./");
   if (progpos!=string::npos) {progname.erase(progpos,2);}
   string str=argv[pos];
   str.erase(0,2);
   if (str==string("version")) {
      cout << progname << " " << CURRENTVERSION << '\n';
      exit(0);
   } else if (str==string("help")) {
      printHelpMenu(argc,argv);
      exit(0);
   } else {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << '\n';
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   return;
}

