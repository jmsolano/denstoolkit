#ifndef _OPTSFLAGS_H
#define _OPTSFLAGS_H

#include <string>
using std::string;

class OptionFlags {
public: 
   OptionFlags();//default constructor, initialize all the flags to convenient (default) values.
   unsigned short int infname,outfname,integrand,setlowerdombox,setupperdombox;
   unsigned short int vegassetpoints,vegassetiter,vegassetconvrat,vegassetinterv,vegassetnpts4max;
   unsigned short int vegassettherm,vegassettol,vegassetstopref;
};//end class optsFlags
void printErrorMsg(char** &argv,char lab);
void printHelpMenu(int &argc, char** &argv);//self-described
void getOptions(int &argc, char** &argv, OptionFlags &flags);//this function will assign the values to
                                                           //all the flags. Implementation is in optsflags.cpp
void processDoubleDashOptions(int &argc,char** &argv,OptionFlags &flags,int pos);
#endif //_OPTSFLAGS_H


