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
   infname=0;
   outfname=0;
   setconvergencerate=0;
   setpoints=0;
   setiterations=0;
   setincrements=0;
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
   //if (argc<2) { //orig line (in DTK)
   if (argc<1) {
      ScreenUtils::SetScrRedBoldFont();
      cout << "\nError: Not enough arguments." << endl;
      ScreenUtils::SetScrNormalFont();
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
      exit(1);
   }
   //if ( string(argv[1])==string("-h")) { //Orig line (in DTK)
   if (argc>1 && string(argv[1])==string("-h")) {
      printHelpMenu(argc,argv);
      exit(1);
   }
   for (int i=1; i<argc; i++){
      if (argv[i][0] == '-'){
         switch (argv[i][1]){
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'n':
               flags.setincrements=(++i);
               if (i>=argc) {printErrorMsg(argv,'n');}
	       break;
            case 'p':
               flags.setpoints=(++i);
               if (i>=argc) {printErrorMsg(argv,'p');}
	       break;
            case 'r' :
               flags.setconvergencerate=(++i);
               if (i>=argc) {printErrorMsg(argv,'r');}
               break;
            case 'i' :
               flags.setiterations=(++i);
               if (i>=argc) {printErrorMsg(argv,'i');}
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
   ScreenUtils::CenterString("This program computes the integral of a property,");
   ScreenUtils::CenterString("using the method Vegas-MonteCarlo");
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
   cout << "  -t tol    \tSets the tolerance to be tol." << '\n';
   cout << "  -o outname\tSet the output file name." << endl
        << "            \t  (If not given the program will create one out of" << endl
        << "            \t  the input name; if given, the gnp file and the pdf will" << endl
        << "            \t  use this name as well --but different extension--)." << endl;
   cout << "  -V        \tDisplays the version of this program." << endl;
   cout << "  -h     \tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      case 't' :
         cout << "should be followed by a small real number." << '\n';
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

