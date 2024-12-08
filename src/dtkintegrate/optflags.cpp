/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
#include "../common/basegausslegendre.h"
#include "../common/basesphtdesign.h"
#include "fldtypesdef.h"

OptionFlags::OptionFlags() {
   infname = outfname = 0;
   integrator = integrand = 0;
   vegassetconvrat = vegassetpoints = vegassetiter    = vegassetinterv = 0;
   vegassettherm   = vegassettol    = vegassetstopref = vegassetnpts4max = 0;
   misersetpoints = misersetdith = 0;
   lsptdsetol = lsptdsetos = 0;
   setlowerdombox = setupperdombox = 0;
   verboseLevel=0;
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
   if ( string(argv[1])==string("-h")) {
      printHelpMenu(argc,argv);
      exit(1);
   }
   for (int i=1; i<argc; i++){
      if (argv[i][0] == '-'){
         switch (argv[i][1]){
            case 'i' :
               flags.integrator=(++i);
               if (i>=argc) {printErrorMsg(argv,'i');}
               break;
            case 'J' :
               flags.stpspindens=true;
               break;
            case 'o':
               flags.outfname=(++i);
               if (i>=argc) {printErrorMsg(argv,'o');}
               break;
            case 'p' :
               flags.integrand=(++i);
               if (i>=argc) {printErrorMsg(argv,'p');}
               break;
            case 'b' :
               flags.setlowerdombox=(++i);
               flags.setupperdombox=(++i);
               if (i>=argc) {printErrorMsg(argv,'b');}
               break;
            case 'h':
               printHelpMenu(argc,argv);
               exit(1);
               break;
            case 'v' :
               if ((i+1)>=argc) {printErrorMsg(argv,'v');}
               flags.verboseLevel=std::stoi(string(argv[++i]));
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
   ScreenUtils::CenterString("using the method Las Vegas+");
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
   cout << "  -p prop   \tChoose the property to be integrated. prop is a character,\n";
   cout << "            \t  which can be (d is the default value):" << '\n';
   string charFields="dmlESTgLGKekMVPsrvDbqQu";
   for ( size_t i=0 ; i<charFields.size() ; ++i ) {
      cout << "         \t\t" << charFields[i] << " ("
         << GetFieldTypeKeyLong(charFields[i]) << ')' << '\n';
   }
   cout << "  -b L U    \tSet the integration domain to be a cube, whose" << '\n';
   cout << "            \t  left,lower,back corner is (L,L,L) and right,upper,front corner" << '\n';
   cout << "            \t  is (U,U,U). Setting L=U=0 instruct the program to" << '\n';
   cout << "            \t  set the boundaries automatically." << '\n';
   cout << "  -i t      \tSet the integrator type to be t. Here, t is a char of the following list:\n"
        << "            \t\tv Vegas-Monte Carlo\n"
        << "            \t\tm Miser-Monte Carlo\n"
        << "            \t\ts Gauss-Legendre-Spherical-t design cubature.\n"
        << "            \t\t    Warning: this only works for single atoms when the field's\n"
        << "            \t\t    domain is real space, or for any molecule if the field's\n"
        << "            \t\t    domain is momentum space." << '\n';
   cout << "  -J        \tSetup alpha- and beta-spin density matrices." << '\n';
   cout << "  -o outname \tSet the output file name." << endl
        << "             \t  (If not given, the program will create one out of" << endl
        << "             \t  the input name; if given, the gnp file and the pdf will" << endl
        << "             \t  use this name as well --but different extension--)." << endl;
   cout << "  -V         \tDisplay the version of this program." << endl;
   cout << "  -h         \tDisplay the help menu.\n\n";
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   ScreenUtils::PrintScrCharLine('-');
   ScreenUtils::CenterString("Specific options for Vegas-Monte Carlo integrator.");
   ScreenUtils::PrintScrCharLine('-');
   //-------------------------------------------------------------------------------------
   cout << "  --vegas-interv  N  \tSet the number of intervals to be N. N is the number of\n"
        << "                     \t  subregions (intervals) to which the original domain is\n"
        << "                     \t  divided. Default: N=10." << '\n';
   cout << "  --vegas-npts NMCP  \tSet the number of MC-points to sample per iteration\n"
        << "                     \t  to be NMCP. Default: NMCP=10,000." << '\n';
   cout << "  --vegas-iter I     \tSet the maximum number of iterations to be I." << '\n';
   cout << "                     \t  Default: I=20." << '\n';
   cout << "  --vegas-conv-rat A \tSet the damping parameter to be A." << '\n';
   cout << "                     \t  Default: A=1.0." << '\n';
   cout << "  --vegas-therm T    \tSet the first T iterations to be ignored for computing\n"
        << "                     \t  the expected integral. Default value: T=0." << '\n';
   cout << "  --vegas-tol tol    \tSet the tolerance for considering an optimal grid to be\n"
        << "                     \t  tol. Default value: tol=0.0." << '\n';
   cout << "  --vegas-stop-ref X \tSet the num. of iterations where the grid will be\n"
        << "                     \t  refined. Usually, X is not defined." << '\n';
   cout << "  --vegas-searchmax N\tUse N sample points to find the global maximum of prop." << '\n';
   cout << "                     \t  This option is applied only momentum space functions." << '\n';
   cout << "                     \t  Default: N=1,000,000." << '\n';
   //-------------------------------------------------------------------------------------
   ScreenUtils::PrintScrCharLine('-');
   ScreenUtils::CenterString("Specific options for Miser-Monte Carlo integrator.");
   ScreenUtils::PrintScrCharLine('-');
   cout << "  --miser-npts NMCP  \tSet the number of MC-points to sample per iteration\n"
        << "                     \t  to be NMCP. Default: NMCP=500,000." << '\n';
   cout << "  --miser-dith d     \tSet dith=d. Default: d=0.05." << '\n';
   ScreenUtils::PrintScrCharLine('-');
   //-------------------------------------------------------------------------------------
   ScreenUtils::CenterString("Specific options for Gauss-Legendre-Spherical-t Cubature integrator.");
   ScreenUtils::PrintScrCharLine('-');
   cout << "  --lsptd-ol OGL     \tUse OLG order for the Gauss-Legendre terms (radial,\n"
        << "                     \t   default: 64).\n";
   cout << "                     \t  Options: " << BaseGaussLegendre::GetAvailableOrders() << '\n';
   cout << "  --lsptd-os OSTD    \tUse OSTD order for the Sperical-t design terms (angular,\n"
        << "                     \t   default: 21).\n";
   cout << "                     \t  Options: " << BaseSphericalTDesign::GetAvailablePrecisions() << '\n';
   ScreenUtils::PrintScrCharLine('-');
   ScreenUtils::CenterString("Specific options for diatomic systems.");
   ScreenUtils::PrintScrCharLine('-');
   ScreenUtils::CenterString("Cubatures for diatomic systems are built using Nr (Na) points");
   ScreenUtils::CenterString("for the radial (angular) part (see dtkmanual for specific details)");
   ScreenUtils::CenterString("For 3D integrals, the number of abscissas is (3*Na*Nr).");
   ScreenUtils::CenterString("For 6D integrals, the number of abscissas is (9*Na*Na*Nr*Nr/16).");
   ScreenUtils::CenterString("Since the most expensive computation is related to 6D integrals,");
   ScreenUtils::CenterString("In this version, only Nr and Na can be altered and must be");
   ScreenUtils::CenterString("even. If an odd number is provided, DTK will guess a close");
   ScreenUtils::CenterString("even number.");
   cout << "  --diatomic-rord Nr \tUse Nr points for the radial part of the cubatures.\n"
        << "                     \t  This option applies to diatomic systems, Default: Nr=32.\n"
        << "                     \t  Nr is the number of points used for the radial part (see above)." << '\n';
   cout << "                     \t  Options: " << BaseGaussLegendre::GetAvailableOrders() << '\n';
   cout << "  --diatomic-aord Na \tUse Na points for the angular part of the cubatures.\n"
        << "                     \t  This option applies to diatomic systems, Default: Na=32.\n"
        << "                     \t  Na is the number of points used for the angular part (see above)." << '\n';
   cout << "                     \t  Options: " << BaseGaussLegendre::GetAvailableOrders() << '\n';
   ScreenUtils::PrintScrCharLine('-');
   //-------------------------------------------------------------------------------------
}
void printErrorMsg(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" ";
   switch (lab) {
      case 'b' :
         cout << " should be followed by two real numbers." << '\n';
         break;
      case 'o':
         cout << " should be followed by a name." << endl;
         break;
      case 'p' :
         cout << " should be followed by a char." << '\n';
         break;
      case 't' :
         cout << "should be followed by a small real number." << '\n';
         break;
      case 'v' :
         cout << "should be followed by a positive integer." << '\n';
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
   static const string mstbf=" must be followed by ";
   string str=argv[pos];
   str.erase(0,2);
   if (str==string("version")) {
      cout << progname << " " << CURRENTVERSION << endl;
      exit(0);
   } else if (str==string("help")) {
      printHelpMenu(argc,argv);
      exit(0);
   } else if ( str==string("vegas-interv") ) {
      flags.vegassetinterv=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("a real number!"));
         exit(1);
      }
   } else if ( str==string("vegas-npts") ) {
      flags.vegassetpoints=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("vegas-iter") ) {
      flags.vegassetiter=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("vegas-conv-rat") ) {
      flags.vegassetconvrat=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("a real number!"));
         exit(1);
      }
   } else if ( str==string("vegas-therm") ) {
      flags.vegassettherm=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("vegas-tol") ) {
      flags.vegassettol=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("a real number!"));
         exit(1);
      }
   } else if ( str==string("vegas-stop-ref") ) {
      flags.vegassetstopref=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("vegas-searchmax") ) {
      flags.vegassetnpts4max=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("miser-npts") ) {
      flags.misersetpoints=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("miser-dith") ) {
      flags.misersetdith=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("a real number!"));
         exit(1);
      }
   } else if ( str==string("lsptd-ol") ) {
      flags.lsptdsetol=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("lsptd-os") ) {
      flags.lsptdsetos=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("diatomic-rord") ) {
      flags.diatsetrad=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
         exit(1);
      }
   } else if ( str==string("diatomic-aord") ) {
      flags.diatsetang=(++pos);
      if (pos>=argc) {
         ScreenUtils::DisplayErrorMessage(str+mstbf+string("an integer!"));
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
void checkOptionsConsistency(int &argc,char **&argv,OptionFlags &flags) {
   if ( argv[flags.integrator][0]=='d' && (flags.setlowerdombox || flags.setupperdombox) ) {
      ScreenUtils::DisplayWarningMessage("Diatomics integrator and option -b is "
            "not compatible!\nIgnoring option -b...");
   }
}

