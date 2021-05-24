#ifndef _OPTSFLAGS_H
#define _OPTSFLAGS_H

#include <string>
using std::string;

class OptionFlags {
public: 
   OptionFlags();//default constructor, initialize all the flags to convenient (default) values.
   unsigned short int infname,outfname;
   unsigned short int setpoints,setiterations,setconvergenceRate,setintervals;
   unsigned short int settermalization,settolerance,setstopRefinement,setfunction;
};//end class optsFlags
void printErrorMsg(char** &argv,char lab);
void printHelpMenu(int &argc, char** &argv);//self-described
void getOptions(int &argc, char** &argv, OptionFlags &flags);//this function will assign the values to
                                                           //all the flags. Implementation is in optsflags.cpp
void processDoubleDashOptions(int &argc,char** &argv,OptionFlags &flags,int pos);
#endif //_OPTSFLAGS_H


