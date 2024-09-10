#include "optflags.h"
#include "screenutils.h"
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <cstdlib>
using namespace std;
#include <fstream>
using std::ifstream;
#include <string>
using std::string;

OptionFlags::OptionFlags() : OptionFlagsBase() {
   minArgs=2; //Change if needed
   stpspindens=false;
}
OptionFlags::OptionFlags(int &argc,char** &argv) : OptionFlags() {
   /* Remember to initialize local short ints before calling Init()!  */
   Init(argc,argv);
}
OptionFlags::~OptionFlags() {
   
}
void OptionFlags::Init(int &argc, char** &argv) {
   if (!BaseProcessOptions(argc,argv)) {return;}
   GetProgramNames(argv[0]);
   if (argc>1 && string(argv[1])==string("-h")) {
      PrintHelpMenu(argc,argv);
      exitcode=OptionFlagsBase::ExitCode::OFEC_EXITNOERR;
      return;
   }
   if (argc<minArgs) {
      ScreenUtils::SetScrRedBoldFont();
      cout << "\nError: Not enough arguments." << endl;
      ScreenUtils::SetScrNormalFont();
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
      exitcode=OptionFlagsBase::ExitCode::OFEC_EXITERR;
      return;
   }
   //outFileName and inFileName are processed in ProcessBaseOptions
   //cases h an V are also processed there.
   //See OptionFlagsBase class ProcessBaseOptions
   for (int i=1; i<argc; i++){
      if (argv[i][0] == '-'){
         switch (argv[i][1]){
            //add cases here
            //help, and version are handled in BaseProcessOptions
            case 'J' :
               stpspindens=true;
               break;
            case '-':
               BaseProcessDoubleDashOption(argc,argv,i);
               ProcessDoubleDashOption(argc,argv,i);
               break;
            default:
               cout << "\nCommand line error. Unknown switch: " << argv[i] << endl;
               cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
               exitcode=OptionFlagsBase::ExitCode::OFEC_EXITERR;
         }
      }
   }
   return;
}
void OptionFlags::PrintHelpMenu(int argc, char** &argv) {
   ScreenUtils::PrintScrStarLine();
   cout << endl;
   ScreenUtils::CenterString(smileyprogramname);
   cout << endl;
   ScreenUtils::CenterString("This program...");
   cout << endl;
   ScreenUtils::CenterString((string("Compilation date: ")+string(__DATE__)));
   cout << endl;
   ScreenUtils::CenterString(string("Version: ")+string(CURRENTVERSION));
   cout << endl;
   ScreenUtils::CenterString((string(":-) Created by: ")+string(PROGRAMCONTRIBUTORS)+string(" (-:")));
   cout << endl;
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrBoldFont();
   cout << "\nUsage:\n\n\t" << rawprogramname << " [option [value(s)]] ... [option [value(s)]]\n\n";
   ScreenUtils::SetScrNormalFont();
   cout << "Here options can be:\n\n";
   cout << "  -o outfname\tSets the output file name to be outfname." << endl;
   cout << "  -V         \tDisplays the version of this program." << endl;
   cout << "  -h         \tDisplay the help menu.\n\n";
   //-------------------------------------------------------------------------------------
   cout << "  --help    \t\tSame as -h" << endl;
   cout << "  --version \t\tSame as -V" << endl;
   cout << endl;
   ScreenUtils::PrintScrStarLine();
   //-------------------------------------------------------------------------------------
}
void OptionFlags::PrintErrorMessage(char** &argv,char lab) {
   ScreenUtils::SetScrRedBoldFont();
   cout << "\nError: the option \"-" << lab << "\" ";
   switch (lab) {
      case 'a' :
      case 'b' :
      case 'm' :
      case 'r' :
         cout << "should be followed by a real number." << endl;
         break;
      case 'o':
         cout << "should be followed by a name." << endl;
         break;
      case 'S' :
         cout << "should be followed by a string." << endl;
         break;
      default:
         cout << "is triggering an unknown error." << endl;
         break;
   }
   ScreenUtils::SetScrNormalFont();
   cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
   cout << "\nto view the help menu.\n\n";
   exitcode=OptionFlagsBase::ExitCode::OFEC_EXITERR;
   return;
}
void OptionFlags::ProcessDoubleDashOption(int &argc,char** &argv,int &pos) {
   string str=argv[pos];
   str.erase(0,2);
   if (str==string("version")) {
      cout << rawprogramname << " " << CURRENTVERSION << endl;
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
void OptionFlags::GetProgramNames(char* argv0) {
   string progname=argv0;
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}
   rawprogramname=progname;
   smileyprogramname=":-)  ";
   smileyprogramname+=rawprogramname;
   smileyprogramname+="  (-:";

}


