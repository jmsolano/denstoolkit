#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "figname.h"
#include "solscrutils.h"

void printFigletName(void){
   static const int FN_NROWS=5;
   static const string figname[FN_NROWS]={
      "|____/ \\___|_| |_|___/|_|\\___/ \\___/|_|_|\\_\\_|\\__|",\
      "| |_| |  __/ | | \\__ \\| | (_) | (_) | | . \\| | |_ ",\
      "| | | |/ _ \\ '_ \\/ __|| |/ _ \\ / _ \\| | ' /| | __|",\
      "|  _ \\  ___ _ __  __|_   _|__   ___ | | |/ (_) |_ ",\
      " ____                _____           _ _  ___ _   "
   };
   for (int i=0; i<FN_NROWS; i++) {centerString(figname[FN_NROWS-i-1]);}
}

