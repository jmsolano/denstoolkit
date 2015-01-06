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
   Instituto Politécnico Nacional, 
   Unidad Monterrey, Mexico.
   2011
   e-mail: jmsolanoalt@gmail.com
 
   Adscription at the moment the particular implementation for this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
   ------------------------

	 This code is free code; you can redistribute it and/or
	 modify it under the terms of the GNU General Public License
	 as published by the Free Software Foundation; either version 2
	 of the License, or (at your option) any later version.

	 This program is distributed in the hope that it will be useful,
	 but WITHOUT ANY WARRANTY; without even the implied warranty of
	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
	 along with this program; if not, write to the Free Software 
	 Foundation, Inc., 59 Temple Place - Suite 330, 
	 Boston, MA  02111-1307, USA.

   WWW:  http://www.gnu.org/copyleft/gpl.html
	
	----------------------
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

optFlags::optFlags()
{
   infname=0;
   outfname=0;
   setn1=0;
   evdim=0;
   zipdat=0;
   mkplt=0;
   kpgnp=0;
   quiet=1;
   setfld=0;
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
            case '0':
               flags.evdim=(i);
               if ((i+3)>=argc) {printErrorMsg(argv,'0');}
               i+=3;
               break;
            case '1':
               flags.evdim=(i);
               if ((++i)>=argc) {printErrorMsg(argv,'1');}
               break;
            case '2':
               flags.evdim=(i);
               if ((++i)>=argc) {printErrorMsg(argv,'2');}
               break;
            case '3':
               flags.evdim=i;
               break;
            case 'k':
               flags.kpgnp=i;
               break;
            case 'n':
               if (i>=argc) {printErrorMsg(argv,'n');}
               flags.setn1=(++i);
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'p' :
               flags.setfld=(++i);
               if ( i>=argc ) {printErrorMsg(argv,'p');}
               break;
            case 'P':
               flags.mkplt=i;
               break;
            case 'v':
               flags.quiet=0;
               break;
            case 'z':
               flags.zipdat=i;
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
   centerString("This program creates a cub/tsv/dat file and a gnp file. The information for");
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
        << "  -0 x y z  \tEvaluate the momentum density at the point (x,y,z)." << endl
        << "            \t  Here x,y, and z are numbers. No file is created." << endl
        << "  -1  x     \tEvaluate the momentum density on a line. \"x\" sets the coordinate" << endl
        << "            \t  direction, and can take the values x, y or z." << endl
        << "            \t  Here x is a character. It creates a *.dat file." << endl
        << "  -2  xy    \tEvaluate the momentum density on a plane. \"xy\" sets the plane of" << endl
        << "            \t  interest, and can take the values xy,xz,yz." << endl
        << "            \t  Here xy are two characters. It creates a *.tsv file." << endl
        << "  -3        \tEvaluate the momentum density on a cube." << endl
        << "            \t  It creates a *.cub file." << endl
        << "  -n  dim   \tSet the number of points for the cub/tsv/dat file per direction" << endl
        << "  -o outname\tSet the output file name." << endl
        << "            \t  (If not given the program will create one out of" << endl
        << "            \t  the input name; if given, the dat/tsv/cub/gnp/pdf files will" << endl
        << "            \t  use this name as well --but different extension--)." << endl;
   cout << "  -p prop   \tChoose the property to be computed. prop is a character," << endl
        << "            \t  which can be (d is the default value):" << endl
        << "            \t     d Density (momentum density)" << endl
        << "            \t     K Kinetic Energy Density (in momentum space)" << endl;;
#if _HAVE_GNUPLOT_
   cout << "  -P     \tCreate a plot using gnuplot. (Only works with options -1 or -2)" << endl
        << "  -k     \tKeeps the *.gnp file to be used later by gnuplot." << endl;
#endif
   cout << "  -v     \tVerbose (display extra information, usually output from third-" << endl
        << "         \t  party sofware such as gnuplot, etc.)" << endl;
#if (defined(__APPLE__)||defined(__linux__))
   cout << "  -z     \tCompress the tsv file using gzip (which must be intalled" << endl
        << "         \t   in your system)." << endl;
#endif
   cout << "  -V        \tDisplays the version of this program." << endl;
   cout << "  -h     \tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   //-------------------------------------------------------------------------------------
#if _HAVE_GNUPLOT_
   printScrStarLine();
   centerString(string("Note that the following programs must be properly installed in your system:"));
   centerString(string("gnuplot"));
   centerString(string("epstool"));
   centerString(string("epstopdf"));
#if (defined(__APPLE__)||defined(__linux__))
   centerString(string("gzip"));
#endif
   printScrStarLine();
#endif
}//end printHelpMenu

//**************************************************************************************************
void printErrorMsg(char** &argv,char lab)
{
   setScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case '0':
         cout << "should be followed by three real numbers." << endl;
         break;
      case '1':
         cout << "should be followed by a character." << endl;
         break;
      case '2':
         cout << "should be followed by a string of two characters" << endl;
         break;
      case 'n':
         cout << "should be followed by an integer." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      case 'p' :
         cout << "should be followed by a character." << endl;
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
void processDoubleDashOptions(int &argc,char** &argv,optFlags &flags,int pos)
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
