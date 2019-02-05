/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.3.1
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
#ifndef _OPTSFLAGS_CPP
#define _OPTSFLAGS_CPP

#include "optflags.h"
#include "../common/solscrutils.h"

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <cstdlib>
using std::exit;
using namespace std;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;
#include "../common/gausswavefunction.h"

//**************************************************************************************************

optFlags::optFlags()
{
   infname=0;
   outfname=0;
   prop2plot=0;
   setn1=0;
   setn3=0;
   setsmcub=0;
   setsmcub1=0;
   zipcube=0;
   wrtlog=0;
   configspecialnci=0;
}


//**************************************************************************************************

void getOptions(int &argc, char** &argv, optFlags &flags)
{
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
      setScrRedBoldFont();
      cout << "\nError: Not enough arguments." << endl;
      setScrNormalFont();
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
}//end updateFlags


//**************************************************************************************************

void printHelpMenu(int &argc, char** &argv)
{
   string progname=argv[0];
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}
   printScrStarLine();
#if _SOL_USE_FIGLET_NAME_
   printFigletName();
#endif
   cout << endl;
   centerString((string(":-) ")+progname+string(" (-:")));
   cout << endl;
   centerString("This program creates a cub file. The information for");
   centerString("the calculation is obtained from a wfx(wfn) file, which is");
   centerString("given as the input for the program.");
   centerString("(See below for the sintax.)");
   cout << endl;
   centerString((string("Compilation date: ")+string(__DATE__)));
   cout << endl;
   centerString(string("Version: ")+string(CURRENTVERSION));
   cout << endl;
   centerString((string(":-) Created by: ")+string(PROGRAMCONTRIBUTORS)+string(" (-:")));
   cout << endl;
   printScrStarLine();
   setScrBoldFont();
   cout << "\nUsage:\n\n\t" << progname << " wf?name [option [value(s)]] ... [option [value(s)]]\n\n";
   setScrNormalFont();
   cout << "Where wf?name is the input wfx(wfn) name, and options can be:\n\n"
        << "  -l        \tWrite cpu time, input/output information etc. on a log file" << endl
        << "  -n  dim   \tSet the number of points per direction for the cube" << endl
        << "            \t  to be dim x dim x dim." << endl
        << "  -N nx ny nz\tSet the individual points per direction for the cube" << endl
        << "            \t  to be nx x ny x nz." << endl
        << "  -o outname\tSet the output file name." << endl
        << "  -s        \tUse a smart cuboid for the grid. The number of points for the" <<endl
        << "            \t  largest direction will be " << DEFAULTPOINTSPERDIRECTION << "." << endl;
   cout << "  -S ln     \tUse a smart cuboid for the grid. ln is the number of points" << endl
        << "            \t  the largest axis will have. The remaining axes will have" << endl
        << "            \t  a number of points proportional to its length." << endl
        << "  -p prop\tChoose the property to be computed. prop is a character," << endl 
        << "         \t  which can be (d is the default value): " << endl
        << "         \t\td (Density)" << endl
        << "         \t\tg (Magnitude of the Gradient of the Density)" << endl
        << "         \t\tl (Laplacian of density)" << endl
        << "         \t\tK (Kinetic Energy Density K)" << endl
        << "         \t\tG (Kinetic Energy Density G)" << endl
        << "         \t\tE (Electron Localization Function -ELF-)" << endl
        << "         \t\tL (Localized Orbital Locator -LOL-)" << endl
        << "         \t\tM (Magnitude of the gradient of LOL)" << endl
        << "         \t\tP (Magnitude of Localized Electrons Detector -LED-)" << endl
        << "         \t\tr (Region of Slow electrons -RoSE-)" << endl
        << "         \t\ts (Reduced Density Gradient -s-)" << endl
        << "         \t\tS (Shannon Entropy Density)" << endl;
   cout << "         \t\tV (Molecular Electrostatic Potential)" << endl;
   cout << "         \t\tu (Scalar Custom Field)" << endl;
   cout << "         \t\tv (Virial Potential Energy Density)" << endl;
   cout << "         \t\tZ (Non Covalent Interactions(NCI) -- s)" << endl;
   cout << "         \t\tz (Non Covalent Interactions(NCI) -- Rho)" << endl;
#if (defined(__APPLE__)||defined(__linux__))
   cout << "  -z     \tCompress the cube file using gzip (which must be installed" << endl
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
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
}//end printHelpMenu

//**************************************************************************************************
void printErrorMsg(char** &argv,char lab)
{
   setScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
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
   setScrNormalFont();
   cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
   cout << "\nto view the help menu.\n\n";
   exit(1);
   return;
}
//**************************************************************************************************
void processDoubleDashOptions(int &argc,char** &argv,optFlags &flags,int &pos)
{
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
         displayErrorMessage(string("configure-nci must be followed by 3 real numbers!"));
         exit(1);
      }
   } else {
      setScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << endl;
      setScrNormalFont();
      exit(1);
   }
   return;
}
//**************************************************************************************************
//**************************************************************************************************
#endif //_OPTSFLAGS_CPP
