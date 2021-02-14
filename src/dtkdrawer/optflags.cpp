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
/* optflags.cpp
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
#include "../common/screenutils.h"

OptionFlags::OptionFlags() {
   infname=1;
   outfname=0;
   verboseLevel=0;
   rotatemol=align3at=rotX=rotY=rotZ=0;
   cpkview=false;
   setzoom=selectcps2draw=0;
   drawnuc=drawcps=mkpng=true;

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
               flags.drawnuc=true;
               break;
            case 'c' :
               flags.selectcps2draw=(++i);
               flags.drawcps=true;
               if (i>=argc) {printErrorMsg(argv,'c');}
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'O' :
               flags.align3at=(++i);
               if ((i+2)>=argc) {printErrorMsg(argv,'O');}
               break;
            case 's' :
               flags.cpkview=true;
               break;
            case 'h':
               printHelpMenu(argc,argv);
               exit(1);
               break;
            case 'v' :
               flags.verboseLevel=(++i);
               if (i>=argc) {printErrorMsg(argv,'v');}
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
               flags.rotatemol=true;
               break;
            case 'y' :
               flags.rotY=(++i);
               flags.rotatemol=true;
               break;
            case 'z' :
               flags.rotZ=(++i);
               flags.rotatemol=true;
               break;
            case 'Z' :
               flags.setzoom=(++i);
               if (i>=argc) {printErrorMsg(argv,'Z');}
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
   if ( flags.align3at ) { flags.rotatemol=true; }
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
   ScreenUtils::CenterString("This program draws properties of a molecule.");
   ScreenUtils::CenterString("The molecule can be parsed to the program through a");
   ScreenUtils::CenterString("file in any of the following formats:");
   ScreenUtils::CenterString("wfx, wfn, xyz");
   cout << endl;
   ScreenUtils::CenterString((string("Compilation date: ")+string(__DATE__)));
   cout << endl;
   ScreenUtils::CenterString(string("Version: ")+string(CURRENTVERSION));
   cout << endl;
   ScreenUtils::CenterString((string(":-) Created by: ")+string(PROGRAMCONTRIBUTORS)+string(" (-:")));
   cout << endl;
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrBoldFont();
   cout << "\nUsage:\n\n\t" << progname << " inputmolecule [option [value(s)]] ... [option [value(s)]]\n\n";
   ScreenUtils::SetScrNormalFont();
   cout << "Where inputmolecule is the input file name, and options can be:\n\n";
   cout << "  -a        \tDraw transparent spheres around each nucleus." << '\n';
   cout << "  -c  cpt   \tSelects the Critical points to include (by default,\n"
        << "            \t  all CPs are included. cpt is a series\n"
        << "            \t  of letters that can be:\n"
        << "            \t\ta: include ACPs\n"
        << "            \t\tb: include BCPs\n"
        << "            \t\tr: include RCPs\n"
        << "            \t\tc: include CCPs\n" << '\n';
   cout << "  -o outname\tSets the output file name.\n"
        << "            \t  (If not given the program will create one out of\n"
        << "            \t  the input name; if given, the gnp file and the pdf will\n"
        << "            \t  use this name as well --but different extension--)." << '\n';
   cout << "  -O a b c  \tOrients the molecule using three atoms (indices: a, b, c\n"
        << "            \t according to the molecule's input geometry). The order is\n"
        << "            \t  important, as the align vectors are defined as follows\n"
        << "            \t  (x, y, and z are, here, the vectors left-right,\n"
        << "            \t  bottom-up, and screen-to-user, relative to the screen):\n"
        << "            \t  x=(c-b), y=a-dotProduct(c,(b-a)), and z=crossproduct(x,y).\n"
        << "            \t  The atoms will be placed in the screen plane, and they will\n"
        << "            \t  look as in scheme 1, below.\n";
   cout << "  -s        \tSets spacefilling mode. This mode activates the spacefilling\n"
        << "            \t  view of the atoms (uses VdW atomic radius to draw the atoms)." << '\n';
   cout << "  -v verbLev\tSets the verbose level to be verboseLevel. The greater\n" 
        << "            \t  verbLev is, the greater the information displayed on the\n"
        << "            \t  screen. verboseLevel=0 minimizes the information."  << '\n';
   cout << "  -V        \tDisplays the version of this program." << endl;
   cout << "  -x alpha  \tSets the angle GNUPlotAngle1=alpha (in povray file)." << '\n';
   cout << "  -y beta   \tSets the angle YAngle=beta (in povray file)." << '\n';
   cout << "  -z gamma  \tSets the angle GNUPlotAngle2=gamma (in povray file)." << '\n';
   cout << "  -Z zmfact \tSets the zoom factor to be zmfact. This will alter the camera\n"
        << "            \t  position by a factor zmfact." << '\n';
   cout << "  -h     \tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   cout << "  --no-png  \t\tPrecludes the png rendering; povray will not be called." << endl;
   //-------------------------------------------------------------------------------------
#if _HAVE_GNUPLOT_
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::CenterString(string("Note that the following programs must be properly installed in your system:"));
   ScreenUtils::CenterString(string("gnuplot"));
   ScreenUtils::CenterString(string("epstopdf"));
#if (defined(__APPLE__)||defined(__linux__))
   ScreenUtils::CenterString(string("gzip"));
#endif
   ScreenUtils::PrintScrStarLine();
#endif
   ScreenUtils::PrintScrCharLine('-');
   cout << "            \t           a\n"
        << "            \t           |\n"
        << "            \t          y|\n"
        << "            \t           |--------c\n"
        << "            \t           |________|______\n"
        << "            \t           / b       x       \n"
        << "            \t          /\n"
        << "            \t       z / \n"
        << "  Scheme 1: View of the aligned atoms (see option -O)." << '\n';
   ScreenUtils::PrintScrCharLine('-');
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'c' :
         cout << "should be followed by a string." << '\n';
         break;
      case 'v':
         cout << "should be followed by an integer (>=0)." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      case 'O' :
         cout << "should be followed by three integers (>0 and <nAtoms)." << '\n';
         break;
      case 'x' :
      case 'y' :
      case 'z' :
      case 'Z' :
         cout << "should be followed by a real number." << '\n';
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
   } else if ( str==string("no-png") ) {
     flags.mkpng=false;
   } else {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << endl;
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   return;
}

