#ifndef _OPTIONFLAGSBASE_CPP
#define _OPTIONFLAGSBASE_CPP
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;
#include "optionflagsbase.h"
#include "screenutils.h"

/* It is important to initialize to 0 all flags, otherwise segmentation faults
 * may occur.  */
OptionFlagsBase::OptionFlagsBase() {
   inFileName=1;
   haveFigletName=false;
   figletName="";
   minArgs=2;
   exitcode=ExitCode::OFEC_CONTINUE;
   verboseLevel=0;
}
string OptionFlagsBase::GetRawProgramName(char** &argv) {
  string res(argv[0]);
  size_t pos=res.find_last_of('/');
  if ( pos!=string::npos ) {
     res.erase(0,pos+1);
  }
  return res;
}
string OptionFlagsBase::GetSmileyProgramName(char** &argv) {
   string res=":-)  "+GetRawProgramName(argv)+"  (-:";
   return res;
}
bool OptionFlagsBase::BaseProcessOptions(int &argc,char** &argv) {
   string progname=GetSmileyProgramName(argv);
   if (argc<minArgs) {
      ScreenUtils::DisplayErrorMessage("Not enough arguments.");
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
      exitcode=OptionFlagsBase::ExitCode::OFEC_EXITERR;
      return false;
   }
   //*
   for (int i=1; i<argc; i++){
      if (argv[i][0] == '-'){
         switch (argv[i][1]){
            case 'i' :
               inFileName=(++i);
               if (i>=argc) {BasePrintErrorMessage(argv,'i');}
               break;
            case 'o' :
               outFileName=(++i);
               if (i>=argc) {BasePrintErrorMessage(argv,'o');}
               break;
            case 'h':
               //PrintHelpMenu(argc,argv);
               exitcode=OptionFlagsBase::ExitCode::OFEC_EXITNOERR;
               break;
            case 'v' :
               verboseLevel=(++i);
               if (i>=argc) {BasePrintErrorMessage(argv,'v');}
               break;
            case 'V':
               progname=GetRawProgramName(argv);
               cout << progname << " " << CURRENTVERSION << endl;
               exitcode=OptionFlagsBase::ExitCode::OFEC_EXITNOERR;
               break;
            case '-':
               BaseProcessDoubleDashOption(argc,argv,i);
               break;
            default:
               exitcode=OptionFlagsBase::ExitCode::OFEC_CONTINUE;
         }
      }
   }
   // */
   return true;
}
void OptionFlagsBase::PrintFigletName(void) {
   cout << figletName << endl << endl;
}
void OptionFlagsBase::PrintHelpMenu(int &argc,char** &argv) {
   // This is the base help menu. If you want to overload it, copy/paste
   // the code of this function, and modify it accordingly.
   string progname=GetRawProgramName(argv);
   ScreenUtils::PrintScrStarLine();
   if ( haveFigletName ) {
      PrintFigletName();
   }
   ScreenUtils::CenterString((string(":-) ")+progname+string(" (-:")));
   cout << endl;
   ScreenUtils::CenterString("This program ...");
   ScreenUtils::CenterString("See below for the sintax.");
   cout << endl;
   ScreenUtils::CenterString((string("Compilation date: ")+string(__DATE__)));
   cout << endl;
   ScreenUtils::CenterString(string("Version: ")+string(CURRENTVERSION));
   cout << endl;
   ScreenUtils::CenterString((string(":-) Created by: ")+string(PROGRAMCONTRIBUTORS)+string(" (-:")));
   cout << endl;
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrBoldFont();
   cout << "\nUsage:\n\n\t" << progname << " name.dat [option [value(s)]] ... [option [value(s)]]\n\n";
   ScreenUtils::SetScrNormalFont();
   cout << "Where name.dat is the input dat-file name, and options can be:\n\n";
   cout << "  -v VL  \tSet the verbose level to be VL." << endl;
   cout << "  -V     \tDisplays the version of this program." << endl;
   cout << "  -h     \tDisplay the help menu.\n\n";
   //--------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   cout << endl;
   //--------------------------------------------------------------------------
}
void OptionFlagsBase::BasePrintErrorMessage(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"" << lab << "\" should be followed by ";
   switch (lab) {
      case 'i':
      case 'o':
         cout << "a name." << endl;
         break;
      case 'v' :
         cout << "an integer." << endl;
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
void OptionFlagsBase::BaseProcessDoubleDashOption(int &argc,char** &argv,int &pos) {
   string str=argv[pos];
   str.erase(0,2);
   if (str==string("version")) {
      cout << CURRENTVERSION << endl;
      exitcode=OptionFlagsBase::ExitCode::OFEC_EXITNOERR;
   } else if (str==string("help")) {
      PrintHelpMenu(argc,argv);
      exitcode=OptionFlagsBase::ExitCode::OFEC_EXITNOERR;
   } else {
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: Unrecognized option '" << argv[pos] << "'" << endl;
      ScreenUtils::SetScrNormalFont();
      exitcode=OptionFlagsBase::ExitCode::OFEC_EXITERR;
   }
   return;
}
#endif //_OPTIONFLAGSBASE_CPP
