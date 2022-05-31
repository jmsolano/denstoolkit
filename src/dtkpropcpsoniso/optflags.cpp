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
/* optflags.cpp
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
   ------------------------
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

OptionFlags::OptionFlags() {
   infname=0;
   outfname=0;
   prop2eval=isoprop=verboselevel=0;
   setcentat=setdirat1=setdirat2=setdirat3=0;
   setviewangles=setgnpangles=0;
   setisovalue=0;
   refinemesh=0;
   isofromcube=0;
   mkpov=kppov=mkpng=false;
   transparentiso=false;
   drawiso=true;
   cpkview=true;
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
            case 'c':
               flags.setcentat=(++i);
               if (i>=argc) {printErrorMsg(argv,'c');}
               break;
            case 'C' :
               flags.setdirat3=(++i);
               if ( (i+2)>=argc ) { printErrorMsg(argv,'C'); }
               i+=2;
               break;
            case 'H' :
               flags.drawiso=false;
               break;
            case 'I' :
               flags.isoprop=(++i);
               if ( (i+1)>=argc ) { printErrorMsg(argv,'I'); }
               ++i;
               break;
            case 'k' :
               flags.kppov=true;
               break;
            case 'l':
               flags.setdirat1=(++i);
               if (i>=argc) {printErrorMsg(argv,'l');}
               break;
            case 'm' :
               flags.setdirat2=(++i);
               if ( (i+1)>=argc ) { printErrorMsg(argv,'m'); }
               ++i;
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'p' :
               flags.prop2eval=(++i);
               if (i>=argc) {printErrorMsg(argv,'p');}
               break;
            case 'P' :
               flags.mkpov=true;
               flags.mkpng=true;
               break;
            case 'r' :
               flags.refinemesh=(++i);
               if (i>=argc) {printErrorMsg(argv,'r');}
               break;
            case 's' :
               flags.cpkview=false;
               break;
            case 't' :
               flags.transparentiso=true;
               break;
            case 'h':
               printHelpMenu(argc,argv);
               exit(EXIT_SUCCESS);
               break;
            case 'v' :
               flags.verboselevel=(++i);
               if (i>=argc) {printErrorMsg(argv,'v');}
               break;
            case 'V':
               progname=argv[0];
               pos=progname.find("./");
               if (!(pos==string::npos)) {progname.erase(pos,2);}
               cout << progname << " " << CURRENTVERSION << endl;
               exit(EXIT_SUCCESS);
               break;
            case '-':
               processDoubleDashOptions(argc,argv,flags,i);
               break;
            default:
               cout << "\nCommand line error. Unknown switch: " << argv[i] << endl;
               cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
               exit(EXIT_FAILURE);
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
   ScreenUtils::CenterString("This program searches for the critical points of a field");
   ScreenUtils::CenterString("evaluated at an isosurface.");
   ScreenUtils::CenterString("Both property evaluated and the isosurface can be");
   ScreenUtils::CenterString("any of the fields implemented in DTK.");
   ScreenUtils::CenterString("The isosurface can be constructed in two forms.");
   ScreenUtils::CenterString("a) The isosurface is evaluated around a single atom");
   ScreenUtils::CenterString("chosen through options (see below), and only on");
   ScreenUtils::CenterString("a cap pointing backwards another point; this point");
   ScreenUtils::CenterString("is selected through several options, below.");
   ScreenUtils::CenterString("b) The isosurface is extracted from a guassian cube,");
   ScreenUtils::CenterString("which can be obtained with dtkcube. In this variant,");
   ScreenUtils::CenterString("the isosurface is not partitioned, but computed through");
   ScreenUtils::CenterString("the complete volume contained in the cube file.");
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
   cout << "  -c cAtNum \tSets the atom cAtNum to be the center around which the\n"
        << "            \t  cap isosurface is computed." << '\n';
   cout << "  -H        \tHides the isosurface in png image (pov file)." << '\n';
   cout << "  -k        \tKeeps the pov-ray script (see also option P, below)." << '\n';
   cout << "  -l dAtNum \tSets the direction point to be the coordinates of atom\n"
        << "            \t  dAtNum. See the scheme at the end of this help menu." << '\n';
   cout << "  -m a1 a2  \tSets the direction point to be the middle point located\n"
        << "            \t  between the coordinates of atoms a1 and a2." << '\n'
        << "            \t  See scheme at the end of this menu.\n";
   cout << "  -C a1 a2 a3\tSets the direction point to be the centroid of the\n"
        << "            \t  triangle formed by the coordinates of atoms\n"
        << "            \t  a1, a2, and a3.\n"
        << "            \t  See scheme at the end of this menu." << '\n';
   cout << "  -p prop   \tSets the property evaluated at the isosurface to be\n"
        << "            \t  prop, which is a char that can be any of the listed\n"
        << "            \t  fields enumerated below for option -I." << '\n'
        << "            \t  (This option is included for future implementations,\n"
        << "            \t  and in this version, only rho[d] is implemented).\n";
   cout << "  -I prop v \tSets the field to compute the isosurface to be prop,\n"
        << "            \t  and the isosurface value to be v. prop is a char,\n"
        << "            \t  which can be (this is valid for options -I and -p)\n"
        << "            \t  (This option is included for future implementations,\n"
        << "            \t  and in this version, only MEP[V] is implemented):\n"
        << "         \t\td (Density)\n"
        << "         \t\tg (Magnitude of the Gradient of the Density)\n"
        << "         \t\tl (Laplacian of density)\n"
        << "         \t\tK (Kinetic Energy Density K)\n"
        << "         \t\tG (Kinetic Energy Density G)\n"
        << "         \t\tE (Electron Localization Function -ELF-)\n"
        << "         \t\tL (Localized Orbital Locator -LOL-)\n"
        << "         \t\tM (Magnitude of the Gradient of LOL)\n"
        << "         \t\tP (Magnitude of the Localized Electrons Detector -LED-)\n"
        << "         \t\tr (Region of Slower Electrons -RoSE-)\n"
        << "         \t\ts (Reduced Density Gradient -s-)\n"
        << "         \t\tS (Shannon Entropy Density)\n"
        << "         \t\tV (Molecular Electrostatic Potential)\n";
   cout << "  -o basename\tSet the base name for the output files to be basename." <<'\n'
        << "            \t  (If not given the program will use the wf?name as the" <<'\n'
        << "            \t  the base name.)" << '\n';
   cout << "  -P        \tGenerates a pov-ray script and renders it. Notice: this requires" << endl
        << "            \t   povray to be installed in your system." << endl;
   cout << "  -r rlev   \tRefines the isosurface cap mesh. rlev is a number used to\n"
        << "            \t  perform a series of subdivisions of an initial icosahedron \n"
        << "            \t  (rlev=4 by default). Each iteration increases\n"
        << "            \t  the number of triangles by a factor of 4, and it must be >0.\n"
        << "            \t  rlev>8 is not recommended as it may cause numerical issues." << '\n';
   cout << "  -s        \tSets spacefilling mode. This mode activates the spacefilling\n"
        << "            \t  view of the atoms (uses VdW atomic radius to draw the atoms)." << '\n';
   cout << "  -t        \tDraw transparent isosurface." << '\n';
   cout << "  -v vrbslev\tSets the verbose level to be vrbslev. This prints information according to\n"
        << "            \t  the level chosen. Default: vrbslev=0." << '\n';
   cout << "  -V     \tDisplays the version of this program." << '\n';
   cout << "  -h     \tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --use-cube cub \t\tUses the cube file 'cub' to generate the isosurface,\n"
        << "                 \t\t  as opposed to generate a cap around a selected atom\n"
        << "                 \t\t  and direction. This option overrides all options\n"
        << "                 \t\t  related to the cap. E.g. the mesh cannot be refined,\n"
        << "                 \t\t  because the mesh is determined from the cube sampling;\n"
        << "                 \t\t  " << '\n';
   cout << "  --help    \t\tSame as -h" << '\n';
   cout << "  --version \t\tSame as -V" << '\n';
   //-------------------------------------------------------------------------------------
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::PrintScrCharLine('-');
   cout << "Scheme of the cap drawn when the isosurface is option a).\n"
        <<"(See description of the program following DTK logo, above.)\n"
        << "Below, r(c) is the radius vector" << '\n';
   cout << "of the center atom (given through option -c), and r(d)" << '\n'
        << "is the radius vector of the direction point (given through options\n"
        << "-l, -m, or -C).\n"
        << "The vector r(c)-r(d) is used to determine the orientation of the\n"
        << "isosurface cap, which is located in the opposite direction of r(d),\n"
        << "relative to r(c), i.e. the cap points in the same direction as\n"
        << "r(c)-r(d):\n";
   ScreenUtils::PrintScrCharLine('-');
   cout << "                                     _" << '\n';
   cout << "                                       -" << '\n';
   cout << "                                         \\" << '\n';
   cout << "                                           \\" << '\n';
   cout << "                                            |" << '\n';
   cout << "       o  -------[r(c)-r(d)]------> o        |<---Isosurface CAP" << '\n';
   cout << "       |                            |       |" << '\n';
   cout << "       |                   r(c)_____|      /" << '\n';
   cout << "       |                                 /" << '\n';
   cout << "     r(d)                            _ -\n" << '\n';
   ScreenUtils::PrintScrCharLine('-');
   ScreenUtils::PrintScrStarLine();
   //-------------------------------------------------------------------------------------
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case '2' :
         cout << "should be followed by a real number." << '\n';
         break;
      case 'A' :
         cout << "should be followed by three real numbers." << '\n';
         break;
      case 'a' :
         cout << "should be followed by two real numbers." << endl;
         break;
      case 'c':
      case 'l':
      case 'r':
      case 'v':
         cout << "should be followed by an integer." << endl;
         break;
      case 'm':
         cout << "should be followed by two integers." << endl;
         break;
      case 'C' :
         cout << "should be followed by three integers." << endl;
         break;
      case 'i':
      case 'o':
      case '1': // see long option "use-cube"
         cout << "should be followed by a string." << endl;
         break;
      case 'p':
         cout << "should be followed by a character." << endl;
         break;
      case 'I' :
         cout << "should be followed by a character and a real number." << '\n';
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
   } else if (str==string("use-cube")) {
      flags.isofromcube=(++pos);
      if (pos>=argc) {printErrorMsg(argv,'1');}
   } else if (str==string("set-isovalue")) {
      flags.setisovalue=(++pos);
      if (pos>=argc) {printErrorMsg(argv,'2');}
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

