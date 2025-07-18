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
   infname=0;
   outfname=0;
   zipdat=0;
   mkpov=0;
   kppov=0;
   mkpng=0;
   quiet=1;
   cptype=0;
   camvdir=0;
   drawnuc=0;
   calcbgps=1;
   calcrgps=0;
   drawbgps=bgptubes=0;
   mkdatmat=mkextsearch=0;
   forcebcpconn=customseedtwoacps=forceseveralbcpconn=0;
   setbgpstep=0;
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
               flags.drawnuc=i;
               break;
            case 'c':
               flags.camvdir=(++i);
               if (i>=argc) {printErrorMsg(argv,'c');}
               break;
            case 'e' :
               flags.mkextsearch=i;
               break;
            case 'g':
               flags.drawbgps=i;
               break;
            case 'G':
               flags.calcbgps=0;
               break;
            //case 'i':
            //   flags.infname=(++i);
            //   if (i>=argc) {printErrorMsg(argv,'i');}
            //   break;
            case 'J' :
               flags.stpspindens=true;
               break;
            case 'k':
               flags.kppov=i;
               break;
            case 'm':
               flags.mkdatmat=i;
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'p':
               flags.mkpov=i;
               flags.kppov=i;
               break;
            case 'P':
               flags.mkpng=i;
               break;
            //case 'z':
            //   flags.zipdat=i;
            //   break;
            case 'r' :
               flags.calcrgps=i;
               break;
            case 't':
               flags.cptype=(++i);
               if (i>=argc) {printErrorMsg(argv,'t');}
               break;
            case 'T':
               flags.bgptubes=i;
               break;
            case 'v':
               flags.quiet=0;
               break;
            case 'h':
               printHelpMenu(argc,argv);
               exit(1);
               break;
            case 'V':
               progname=argv[0];
               pos=progname.find("./");
               if (!(pos==string::npos)) {progname.erase(pos,2);}
               cout << progname << " " << CURRENTVERSION << endl;
               exit(0);
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
   ScreenUtils::CenterString("This program looks for the critical points of LOL or Electron Density. ");
   ScreenUtils::CenterString("It is spected as input a wfx (wfn) file.");
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
        << "  -a        \tDraw transparent spheres around each nuclei in the wf? file." <<endl
        << "  -o outname\tSet the output file name (*.log)." << endl
        << "            \t  (If not given the program will create one out of" << endl
        << "            \t  the input name; if given, the pov file will" << endl
        << "            \t  use this name as well --but different extension--)." << endl;
   cout << "  -m        \tSave the coordinates of atoms, critical points and gradient" << endl
        << "            \t  paths (if calculated, see option -G) into *dat files." << endl;
#if _HAVE_POVRAY_
   cout << "  -p     \tCreate a pov file. No image will be generated." << endl;
   cout << "  -P     \tCreate a png image using povray (calls povray)," << endl
        << "         \t  deletes pov file." << endl;
   cout << "  -c p   \tSet the camera vector direction for the povray scene." << endl
        << "         \t  p is an integer which can take the values" << endl
        << "         \t  001,010,100,101,110,111; the camera will be located" << endl
        << "         \t  somewhere in the direction (x,y,z), where xyz is one of "<< endl
        << "         \t  the above values." << endl;
   cout << "  -g     \tDraw the gradient paths in the POV file." << endl
        << "  -T     \tSet the style for the gradient paths to be tubes." << endl;
   cout << "  -G     \tSkip the calculation of the bond gradient paths." << endl;
   cout << "  -r     \tPerform the search of ring and cage gradient paths." << endl;
   cout << "  -k     \tKeeps the pov file if option -P is used." << endl;
#endif
   cout << "  -t cpt \tSet the type of critical point to be searched. cpt can be "<< endl
        << "         \t  one of the following fields (Density is the default):" << endl
        << "         \t\td (Density)" << endl
        << "         \t\tL (Localized Orbital Locator -LOL-)" << endl;
   cout << "  -e     \tPerform an extended search of critical points. This" << endl
        << "         \t  will take some more time, but it could find more CPs" << endl
        << "         \t  than the simple search." << endl;
   cout << "  -J     \tSetup alpha- and beta-spin density matrices." << '\n';
   cout << "  -v     \tVerbose mode (displays additional information (for example the " << endl
        << "         \t  output of povray, etc." <<endl;
   //#if (defined(__APPLE__)||defined(__linux__))
   //   cout << "  -z     \tCompress the cube file using gzip (which must be intalled" << endl
   //        << "         \t   in your system)." << endl;
   //#endif
   cout << "  -V        \tDisplay the version of this program." << endl;
   cout << "  -h\t\tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --force-bcp-connectivity bcp acp1 acp2 \tTry to force the connectivity" << endl
        << "            \t\t  of the bcp to be acpA and acpB. Needed to correct some" << endl
        << "            \t\t  cases wherein bond paths are not properly drawn." << endl
        << "            \t\t  The option will try to force the bond path to be" << endl
        << "            \t\t  connected between acpA-bcp-acpB." << endl;
   cout << "  --force-bcp-connectivities n bcp1 acp1A acp1B bcp2 acp2A acp2B" << endl
        << "            \t\t  Try to connect more than one bcp (see --force-bcp-connectivity)" << endl
        << "            \t\t  n is the number of bcp that will be tested for connection." << endl
        << "            \t\t  Connectibities will be tried to be acpiA-bcpi-acpiB," << endl
        << "            \t\t  here i is the number 1, 2, ..., n." << endl;
   cout << "  --add-seed-twoacps acp1 acp2\n"
        << "            \t\tPerform a critical point search around" << endl
        << "            \t\t  the middle point between acp1 and acp2." << endl;
   cout << "  --set-bgp-step dh\n"
        << "            \t\tSet the bond-gradient path step to be dh. Default: "
        << CPNW_DEFAULTGRADIENTPATHS << '\n';
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
#if _HAVE_POVRAY_
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::CenterString(string("Note that the following programs must be properly installed in your system:"));
   ScreenUtils::CenterString(string("povray"));
#if _HAVE_IMAGEMAGICK_
   ScreenUtils::CenterString(string("imagemagick"));
#endif
#if _HAVE_GRAPHICSMAGICK_
   ScreenUtils::CenterString(string("graphicsmagick"));
#endif
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
      case 'c':
         cout << "should be followed by an integer." << endl;
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
   } else if ( str==string("add-seed-twoacps") ) {
     if ( (pos+1)>=(argc+2) ) {
        ScreenUtils::DisplayErrorMessage(string("add-seed-twoacps must be followed by 2 integers!"));
        exit(1);
     }
     flags.customseedtwoacps=(pos+1);
   } else if ( str==string("force-bcp-connectivity") ) {
      if ((pos+1)>=(argc+3)) {
         ScreenUtils::DisplayErrorMessage(string("force-bcp-connectivity must be followed by 3 integers!"));
         exit(1);
      }
      flags.forcebcpconn=(pos+1);
   } else if ( str==string("force-bcp-connectivities") ) {
      ScreenUtils::DisplayWarningMessage("In this version, the user is responsible of providing\n"
            "the correct number of arguments.");
      flags.forceseveralbcpconn=(pos+1);
   } else if ( str==string("set-bgp-step") ) {
      if ( (pos+1)>argc ) {
         ScreenUtils::DisplayErrorMessage("This option must be followed by a number.");
      }
      flags.setbgpstep=(++pos);
   } else {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << endl;
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   return;
}

