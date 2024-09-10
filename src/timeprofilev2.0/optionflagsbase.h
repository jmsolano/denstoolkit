#ifndef _OPTIONFLAGSBASE_H
#define _OPTIONFLAGSBASE_H
#include <string>
using std::string;

class OptionFlagsBase {
   /* ************************************************************************** */
public: 
   enum class ExitCode {
      OFEC_CONTINUE, OFEC_EXITNOERR, OFEC_EXITERR
   };
   OptionFlagsBase();
   virtual bool BaseProcessOptions(int &argc,char** &argv);
   virtual bool ProcessOptions(int &argc,char** &argv) {return true;}
   virtual void PrintHelpMenu(int &argc,char** &argv);
   void BaseProcessDoubleDashOption(int &argc,char** &argv,int &pos);
   virtual void ProcessDoubleDashOption(int &argc,char** &argv,int &pos) {}
   static string GetRawProgramName(char** &argv);
   static string GetSmileyProgramName(char** &argv);
   void SetMinNumberOfArguments(int nn) {minArgs=nn;}
   void PrintFigletName(void);
   OptionFlagsBase::ExitCode GetExitCode() {return exitcode;}
   /* ************************************************************************** */
   unsigned short int inFileName,outFileName;
   unsigned short int verboseLevel;
   /* ************************************************************************** */
protected:
   void BasePrintErrorMessage(char** &argv,char lab);
   bool haveFigletName;
   string figletName;
   int minArgs;
   OptionFlagsBase::ExitCode exitcode;
   string rawprogramname,smileyprogramname;
};
   /* ************************************************************************** */

#endif //_OPTIONFLAGSBASE_H

