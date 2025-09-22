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
#include "fldtypesdef.h"

OptionFlags::OptionFlags() {
   infname=0;
   outfname=0;
   prop2plot=0;
   setn1=0;
   setn3=0;
   setsmcub=0;
   setsmcub1=0;
   setdelta1=0;
   setdelta3=0;
   setcentredcub=0;
   zipcube=0;
   wrtlog=0;
   configspecialnci=0;
   genvmdscript=quietrender=false;
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
            case 'c' :
               if ((i+2)>=argc) {printErrorMsg(argv,'c');}
               flags.setcentredcub=(++i);
               ++i;
               break;
            case 'd' :
               if (i>=argc) {printErrorMsg(argv,'d');}
               flags.setdelta1=(++i);
               break;
            case 'D' :
               if ((i+3)>=argc) {printErrorMsg(argv,'D');}
               flags.setdelta3=(++i);
               i+=2;
               break;
            case 'J' :
               flags.stpspindens=true;
               break;
            case 'l':
               flags.wrtlog=i;
               break;
            case 'n':
               if (i>=argc) {printErrorMsg(argv,'n');}
               flags.setn1=(++i);
               break;
            case 'N':
               if ((i+3)>=argc) {printErrorMsg(argv,'N');}
               flags.setn3=(++i);
               i+=2;
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
               flags.genvmdscript=true;
               break;
            case 'q':
               flags.quietrender=true;
               break;
            case 's':
               flags.setsmcub=(i);
               break;
            case 'S':
               flags.setsmcub1=(++i);
               if (i>=argc) {printErrorMsg(argv,'S');}
               break;
            case 'z':
               flags.zipcube=i;
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
#if _SOL_USE_FIGLET_NAME_
   FigletName::PrintFigletName();
#endif
   cout << endl;
   ScreenUtils::CenterString((string(":-) ")+progname+string(" (-:")));
   cout << endl;
   ScreenUtils::CenterString("This program creates a cub file. The information for");
   ScreenUtils::CenterString("the calculation is obtained from a wfx(wfn) file, which is");
   ScreenUtils::CenterString("given as the input for the program.");
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
   cout << "Where wf?name is the input wfx(wfn) name, and options can be:\n\n";
   cout << "  -d  dd    \tSet the delta between coordinates per direction for the cube" << endl
        << "            \t  to be dd: dx=dy=dz=dd. I.e., the grid is constructed so that\n"
        << "            \t  the distance between points is dd on the x-axis, dd on the y-axis\n"
        << "            \t  and dd on the z-axis." << endl;
   cout << "  -D dx dy dz\tSet the deltas between coordinates per direction for the cube" << endl
        << "            \t  to be dx, dy, and dz. I.e., the grid is constructed so that\n"
        << "            \t  the distance between points is dx on the x-axis, dy on the y-axis,\n"
        << "            \t  and dz on the z-axis." << endl;
   cout << "  -J        \tSetup alpha- and beta-spin density matrices." << '\n';
   cout << "  -l        \tWrites cpu time, input/output information etc. on a log file" << endl;
   cout << "  -n  dim   \tSets the number of points per direction for the cube" << endl
        << "            \t  to be dim x dim x dim." << endl;
   cout << "  -N nx ny nz\tSets the individual points per direction for the cube" << endl
        << "            \t  to be nx x ny x nz." << endl;
   cout << "  -o outname\tSets the output file name." << endl;
   cout << "  -s        \tUses a smart cuboid for the grid. The number of points for the" <<endl
        << "            \t  largest direction will be " << DEFAULTPOINTSPERDIRECTION << "." << endl;
   cout << "  -S ln     \tUses a smart cuboid for the grid. ln is the number of points" << endl
        << "            \t  the largest axis will have. The remaining axes will have" << endl
        << "            \t  a number of points proportional to its length." << endl;
   cout << "  -c a1 a2  \tUses a cube centred at the midpoint of the atoms a1 and a2." << endl;
   cout << "  -p prop\tChooses the property to be computed. prop is a character," << endl 
        << "         \t  which can be (d is the default value): " << endl;
   string charFields="dglKGeELMPrsSVDbquUvzZ";
   for ( size_t i=0 ; i<charFields.size() ; ++i ) {
      cout << "         \t\t" << charFields[i] << " ("
         << GetFieldTypeKeyLong(charFields[i]) << ')' << '\n';
   }
   cout << "  -P     \tGenerates a VMD script to render the field. Notice: this requires VMD and the internal" << endl
        << "         \t  Tachyon render (usually it is included automatically in VMD) to be installed" << endl
        << "         \t  in your system. The VMD script should run with 'vmd -e filename.vmd'" << endl;
   cout << "  -q     \tMakes the vmd script to close VMD automatically, after rendering the field." << endl;
#if (defined(__APPLE__)||defined(__linux__))
   cout << "  -z     \tCompresses the cube file using gzip (which must be installed" << endl
        << "         \t   in your system)." << endl;
#endif
   cout << "  -V        \tDisplays the version of this program." << endl;
   cout << "  -h\t\tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --configure-nci rMin rMax sMax \tSet the parameters rhoMin, rhoMax," << endl
        << "             \t\t  and redGradMax to be rMin, rMax, and sMax, respectively." << endl
        << "             \t\t  This option only affects NCI cubes (see properties z and" << endl
        << "             \t\t  Z in \"-p\" option. Default values: rhoMin=" << NCIRHOMIN << "," << endl
        << "             \t\t  rhoMax= " << NCIRHOMAX << ", and redGradMax= " <<  NCISMAX << endl;
   cout << "  --set-extra-space X \tSet the space around the most-external atoms to be X a.u.\n"
        << "             \t\t  By default, the cube is computed within a grid that includes\n"
        << "             \t\t  all atoms, plus an extra space of X=" << EXTRASPACECUBEFACTOR << " a.u.\n";
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'c' :
         cout << "should be followed by two integers." << endl;
         break;
      case 'p':
         cout << "should be followed by a character." << endl;
         break;
      case 'n':
         cout << "should be followed by an integer." << endl;
         break;
      case 'N':
         cout << "should be followed by three integers." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      case 'S':
         cout << "should be followed by an integer" << endl;
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
   } else if ( str==string("set-extra-space") ) {
     flags.setextraspace=(++pos);
      if ((pos)>=(argc)) {
         ScreenUtils::DisplayErrorMessage(string("set-extra-space must be followed by 1 real number!"));
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

