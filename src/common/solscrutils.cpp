/*
 *  solmiscutil.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 12/11/10.
 *  Copyright 2010 Cinvestav, Unidad Monterrey. All rights reserved.
 ---------------------------------------------------------------------
   Further development:
      University of Guelph,
      Guelph, Ontario, Canada.
      2013
 *
 */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef _SOL_SCR_UTILS_CPP_
#define _SOL_SCR_UTILS_CPP_
#include "solscrutils.h"

//**************************************************************************************************
//**************************************************************************************************
void centerString(const std::string &s)
{
   int pos=(int)((80-s.length())/2);
   for(int i=0;i<pos;i++){std::cout<<" ";}
   std::cout<<s<<std::endl;
}
//**************************************************************************************************
void centerString(const char* word)
{
   centerString(std::string(word));
}
//**************************************************************************************************
void printProgressBar(int perc){
   std::string b;
   for(int i=0; i<50; i++){
      if(i<(perc/2)) {
         b.replace(i,1,"=");
      } else if(i==(perc/2)){
         b.replace(i,1,"|");
      } else {
         b.replace(i,1," ");
      }
   }
   std::cout << "\r" << "[" << b << "]";
   std::cout.width(3);
   std::cout << perc << "%" << std::flush;
   return;
}
//**************************************************************************************************
void printHappyEnding(void)
{
   std::cout << "\n\n";
   centerString("(-: Normal termination! :-)");
   std::cout << "\n\n";
   //system("date \"+\%Y\%m\%d-\%H\%M\"");
}
//**************************************************************************************************
bool isDigit(char c)
{
   if ((c=='0')||(c=='1')||(c=='2')||(c=='3')||(c=='4')||
       (c=='5')||(c=='6')||(c=='7')||(c=='8')||(c=='9')) {
      return true;
   } else {
      return false;
   }

}
//**************************************************************************************************
void printScrStarLine(void)
{
   for (int i=0; i<80; i++) {std::cout << '*';}
   std::cout << std::endl;
}
//**************************************************************************************************
void printScrCharLine(char h)
{
   for (int i=0; i<80; i++) {std::cout << h;}
   std::cout << std::endl;
}
//**************************************************************************************************
void setScrGreenBoldFont(void)
{
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[32m";
#endif
   return;
}
//**************************************************************************************************
void setScrRedBoldFont(void)
{
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[31m";
#endif
   return;
}
//**************************************************************************************************
void setScrYellowBoldFont(void)
{
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[33m";
#endif
   return;
}
//**************************************************************************************************
void setScrBlueBoldFont(void)
{
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[1;34m";
#endif
   return;
}
//**************************************************************************************************
void setScrBoldFont(void)
{
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[1m";
#endif
   return;
}
//**************************************************************************************************
void setScrNormalFont(void)
{
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[m";
#endif
   return;
}
//**************************************************************************************************
void printHappyStart(char ** (&argv),const char *vers,const char *contrib)
{
   std::string progname=argv[0];
   progname.append("  (-:");
   progname.insert(0,":-) ");
   size_t pos=progname.find("./");
   if (pos!=std::string::npos) {progname.erase(pos,2);}
   setScrGreenBoldFont();
   printScrStarLine();
   std::cout << std::endl;
   centerString(progname);
   std::cout << std::endl;
   progname=__DATE__;
   progname.insert(0,"Compilation date: ");
   centerString(progname);
   std::cout << std::endl;
   progname="Version: ";
   //progname.append(CURRENTVERSION);
   progname.append(vers);
   centerString(progname);
   std::cout << std::endl;
   progname=":-)    Created by ";
   //progname.append(PROGRAMCONTRIBUTORS);
   progname.append(contrib);
   progname.append("    (-:");
   centerString(progname);
   std::cout << std::endl;
   printScrStarLine();
   setScrNormalFont();
   return;
}
//**************************************************************************************************
void displayErrorMessage(const std::string &s)
{
   setScrRedBoldFont();
   std::cout << "Error: " << s << std::endl;
   setScrNormalFont();
   return;
}
//**************************************************************************************************
void displayErrorMessage(const char* word)
{
   displayErrorMessage(std::string(word));
   return;
}
//**************************************************************************************************
void displayWarningMessage(const std::string &s)
{
   setScrYellowBoldFont();
   std::cout << "Warning: " << s << std::endl;
   setScrNormalFont();
   return;
}
//**************************************************************************************************
void displayWarningMessage(const char* word)
{
   displayWarningMessage(std::string(word));
   return;
}
//**************************************************************************************************
void printBetweenStarLines(const std::string &s)
{
   printScrStarLine();
   centerString(s);
   printScrStarLine();
   return;
}
//**************************************************************************************************
void printBetweenStarLines(const char* word)
{
   printBetweenStarLines(std::string(word));
   return;
}
//**************************************************************************************************
void printV3Comp(const solreal (&v)[3])
{
   for (int i=0; i<3; i++) {std::cout << v[i] << " ";}
   std::cout << std::endl;
   return;
}
//**************************************************************************************************
void printV3Comp(const std::string &s,const solreal (&v)[3])
{
   std::cout << s;
   printV3Comp(v);
   return;
}
//**************************************************************************************************
void printV3Comp(const char* word,const solreal (&v)[3])
{
   std::cout << word;
   printV3Comp(v);
   return;
}
//**************************************************************************************************
void printV3Comp(const int (&v)[3])
{
   for (int i=0; i<3; i++) {std::cout << v[i] << " ";}
   std::cout << std::endl;
   return;
}
//**************************************************************************************************
void printV3Comp(const std::string &s,const int (&v)[3])
{
   std::cout << s;
   printV3Comp(v);
   return;
}
//**************************************************************************************************
void printV3Comp(const char* word,const int (&v)[3])
{
   std::cout << word;
   printV3Comp(v);
   return;
}
//**************************************************************************************************
void printM3x3Comp(const solreal (&m)[3][3])
{
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {std::cout << m[i][j] << " ";}
      std::cout << std::endl;
   }
   //std::cout << std::endl;
   return;
}
//**************************************************************************************************
void printM3x3Comp(const std::string &s,const solreal (&m)[3][3])
{
   std::cout << s;
   printM3x3Comp(m);
   return;
}
//**************************************************************************************************
void printM3x3Comp(const char* word,const solreal (&m)[3][3])
{
   std::cout << word;
   printM3x3Comp(m);
   return;
}
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#endif //_SOL_SCR_UTILS_CPP_