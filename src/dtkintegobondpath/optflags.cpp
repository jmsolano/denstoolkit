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

/*

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

//**************************************************************************************************

optFlags::optFlags() {
   infname=0;
   outfname=0;
   globalenergy=0;
}


//**************************************************************************************************

void getOptions(int &argc, char** &argv, optFlags &flags) {
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
            case 'e' :
               flags.globalenergy=(++i);
               if (i>=argc) {printErrorMsg(argv,'e');}
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
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
}//end updateFlags


//**************************************************************************************************

void printHelpMenu(int &argc, char** &argv) {
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
   centerString("This program computes a line integral of several scalar fields.");
   centerString("The fields are some of the regular scalar fields implemented");
   centerString("in this version of DensToolKit.");
   centerString("The program will save the information, both total integral");
   centerString("upon all bond paths, i.e. the sum of all individual integrals");
   centerString("over the bond paths.");
   centerString("It also saves the individual information of all integrals.");
   centerString("This information is saved into a *.log file.");
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
        << "  -o outname\tSet the output file name." << endl
        << "            \t  (If not given the program will create one out of" << endl
        << "            \t  the input name; if given, the gnp file and the pdf will" << endl
        << "            \t  use this name as well --but different extension--)." << endl;
   cout << "  -e globEn \tSet the global energy to be globEn. If this is given," << endl
        << "            \t  then all contributions (path integrals) will contain" << endl
        << "            \t  extra columns with entries proportional to globEn." << endl;
   cout << "  -V        \tDisplays the version of this program." << endl;
   cout << "  -h     \tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
/*
#if _HAVE_GNUPLOT_
   printScrStarLine();
   centerString(string("Note that the following programs must be properly installed in your system:"));
   centerString(string("gnuplot"));
   centerString(string("epstopdf"));
#if (defined(__APPLE__)||defined(__linux__))
   centerString(string("gzip"));
#endif
   printScrStarLine();
#endif
*/
}//end printHelpMenu

//**************************************************************************************************
void printErrorMsg(char** &argv,char lab) {
   setScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'e' :
         cout << "should be followed by a real number." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
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
void processDoubleDashOptions(int &argc,char** &argv,optFlags &flags,int pos) {
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
