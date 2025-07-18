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

   optsflags.cpp
   
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
#include "screenutils.h"
#include "fldtypesdef.h"


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
   showatlbls=0;
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
            case 'a':
               if ((i+2)>=argc) {printErrorMsg(argv,'a');}
               flags.setats=(++i);
               i++;
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
            case 's':
               flags.setstep=(++i);
               if (i>=argc) {printErrorMsg(argv,'s');}
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'p':
               flags.prop2plot=(++i);
               if (i>=argc) {printErrorMsg(argv,'p');}
               break;
            case 'P':
               flags.mkplt=i;
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
   FigletName::PrintFigletName();
   cout << endl;
   ScreenUtils::CenterString((string(":-) ")+progname+string(" (-:")));
   cout << endl;
   ScreenUtils::CenterString("This program computes a density field along");
   ScreenUtils::CenterString("the bond path between two atoms (optionally upon the regular straight");
   ScreenUtils::CenterString("line that joins those atoms). The output will be in a dat file.");
   ScreenUtils::CenterString("The information for the calculation is obtained from a wfx(wfn) file,");
   ScreenUtils::CenterString("which is given as the input for the program.");
   ScreenUtils::CenterString("(See below for the sintax.)");
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
   cout << "Where wf?name is the input wfx(wfn) name, and options can be:\n\n"
        << "  -a a1 a2  \tDefine the atoms  (a1,a2) used to define bond path/line." << endl
        << "            \t  If this option is not activated, the program will " << endl
        << "            \t  set a1=1, a2=2." << endl
        << "            \t  Note: if the *.wfn (*.wfx) file has only one atom" << endl
        << "            \t  the program will exit and no output will be generated." << endl;
   cout << "  -J        \tSetup alpha- and beta-spin density matrices." << '\n';
   cout << "  -L        \tCalculate the field upon the straight line that joins the atoms" << endl
        << "            \t  instead of upon the bond path." << endl;
   cout << "  -n  dim   \tSet the number of points for the dat file." << endl
        << "            \t  Note: for the bond path you may want to look for a good " << endl
        << "            \t  combination of n and the number \"step\" given in option -s," << endl
        << "            \t  since the number of points in the bond path will be mainly " << endl
        << "            \t  governed by step." << endl;
   cout << "  -o outname\tSet the output file name." << endl
        << "            \t  (If not given the program will create one out of" << endl
        << "            \t  the input name; if given, the dat, gnp and (eps)pdf files will" << endl
        << "            \t  use this name as well --but different extension--)." << endl;
   cout << "  -p prop   \tChoose the property to be computed. prop is a character," << endl
        << "            \t  which can be (d is the default value): " << endl;
   string charFields="dglKGeELMPrsSVvDbqu";
   for ( size_t i=0 ; i<charFields.size() ; ++i ) {
      cout << "         \t\t" << charFields[i] << " ("
         << GetFieldTypeKeyLong(charFields[i]) << ')' << '\n';
   }
   cout << "  -s step   \tSet the stepsize for the bond path to be 'step'." << endl
        << "            \t  Default value: " << DEFAULTBONDPATHSTEPMD1 << endl;
#if _HAVE_GNUPLOT_
   cout << endl;
   cout << "  -P     \tCreate a plot using gnuplot." << endl;
   //cout << "  -k     \tKeeps the *.gnp file to be used later by gnuplot." << endl;
   cout << "  -l     \tShow labels of atoms (those set in option -a) in the plot." << endl;
#endif
   cout << endl;
   cout << "  -v     \tVerbose (display extra information, usually output from third-" << endl
        << "         \t  party sofware such as gnuplot, etc.)" << endl;
   //cout << "  -V     \tShows the current version of this program." << endl;
#if (defined(__APPLE__)||defined(__linux__))
   cout << "  -z     \tCompress the dat file using gzip (which must be intalled" << endl
        << "         \t   in your system)." << endl;
#endif
   cout << endl;
   cout << "  -h     \tDisplay the help menu.\n";
   cout << "  -V     \tDisplays the version of this program." << endl;
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   cout << endl;
   ScreenUtils::PrintScrCharLine('-');
   cout << "            \tThe format of the dat file is:" << endl;
   cout << "            \t      L  X  Y  Z  V" << endl;
   cout << "            \t  where L is the value of the parameter that maps the bond" << endl;
   cout << "            \t  path to a line; X, Y, and Z are the actual spatial coordinates" << endl;
   cout << "            \t  of each point in the bond path; and V is the value of the" << endl;
   cout << "            \t  chosen field at the point (X,Y,Z) ---see option -p." << endl;
   ScreenUtils::PrintScrCharLine('-');
   //-------------------------------------------------------------------------------------
#if _HAVE_GNUPLOT_
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::CenterString(string("Note that the following programs must be properly installed in your system:"));
   ScreenUtils::CenterString(string("gnuplot"));
   ScreenUtils::CenterString(string("epstopdf"));
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
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
         cout << "should be followed by two integers." << endl;
         break;
      case 'n':
         cout << "should be followed by an integer." << endl;
         break;
      case 's':
         cout << "should be followed by a real number." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      case 'p':
         cout << "should be followed by a character." << endl;
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
void processDoubleDashOptions(int &argc,char** &argv,OptionFlags &flags,int pos) {
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
   } else {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << endl;
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   return;
}

