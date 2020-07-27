/* screenutils.cpp
   ---------------------------------------------------------------------
   Created by Juan Manuel Solano Altamirano on 12/11/10.
   Copyright 2010 Cinvestav, Unidad Monterrey. All rights reserved.
   ---------------------------------------------------------------------
   Further development:
      University of Guelph,
      Guelph, Ontario, Canada.
      2013

      Meritorious Autonomous University of Puebla
      Puebla, Puebla, Mexico.
      2016
*/
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "screenutils.h"

#ifdef _SOL_USE_FIGLET_NAME_
#include "figname.h"
#endif

void ScreenUtils::CenterString(const std::string &s) {
   int pos=(int)((80-s.length())/2);
   for(int i=0;i<pos;i++){std::cout<<" ";}
   std::cout<<s<<std::endl;
}
void ScreenUtils::CenterString(const char* word) {
   CenterString(std::string(word));
}
void ScreenUtils::PrintProgressBar(int perc){
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
void ScreenUtils::PrintHappyEnding(void) {
   std::cout << "\n\n";
   SetScrGreenBoldFont();
   CenterString("(-: Normal termination! :-)");
   std::cout << "\n\n";
   SetScrNormalFont();
   //system("date \"+\%Y\%m\%d-\%H\%M\"");
}
bool ScreenUtils::IsDigit(char c) {
   if ((c=='0')||(c=='1')||(c=='2')||(c=='3')||(c=='4')||
       (c=='5')||(c=='6')||(c=='7')||(c=='8')||(c=='9')) {
      return true;
   } else {
      return false;
   }
}
void ScreenUtils::PrintScrStarLine(void) {
   for (int i=0; i<80; i++) {std::cout << '*';}
   std::cout << std::endl;
}
void ScreenUtils::PrintScrCharLine(char h) {
   for (int i=0; i<80; i++) {std::cout << h;}
   std::cout << std::endl;
}
void ScreenUtils::SetScrGreenBoldFont(void) {
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[32m";
#endif
   return;
}
void ScreenUtils::SetScrRedBoldFont(void) {
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[31m";
#endif
   return;
}
void ScreenUtils::SetScrYellowBoldFont(void) {
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[33m";
#endif
   return;
}
void ScreenUtils::SetScrBlueBoldFont(void) {
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[1;34m";
#endif
   return;
}
void ScreenUtils::SetScrBoldFont(void) {
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[1m";
#endif
   return;
}
void ScreenUtils::SetScrNormalFont(void) {
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[m";
#endif
   return;
}
void ScreenUtils::PrintHappyStart(char ** (&argv),const char *vers,const char *contrib) {
   std::string progname=argv[0];
   progname.append("  (-:");
   progname.insert(0,":-) ");
   size_t pos=progname.find("./");
   if (pos!=std::string::npos) {progname.erase(pos,2);}
   SetScrGreenBoldFont();
   PrintScrStarLine();
#ifdef _SOL_USE_FIGLET_NAME_
#if _SOL_USE_FIGLET_NAME_
   FigletName::PrintFigletName();
#endif
#endif
   std::cout << std::endl;
   CenterString(progname);
   std::cout << std::endl;
   progname=__DATE__;
   progname.insert(0,"Compilation date: ");
   CenterString(progname);
   std::cout << std::endl;
   progname="Version: ";
   //progname.append(CURRENTVERSION);
   progname.append(vers);
   CenterString(progname);
   std::cout << std::endl;
   progname=":-)    Created by ";
   //progname.append(PROGRAMCONTRIBUTORS);
   progname.append(contrib);
   progname.append("    (-:");
   CenterString(progname);
   std::cout << std::endl;
   PrintScrStarLine();
   SetScrNormalFont();
   return;
}
void ScreenUtils::DisplayErrorMessage(const std::string &s) {
   SetScrRedBoldFont();
   std::cout << "Error: " << s;
   SetScrNormalFont();
   std::cout << std::endl;
   return;
}
void ScreenUtils::DisplayErrorMessage(const char* word) {
   DisplayErrorMessage(std::string(word));
   return;
}
void ScreenUtils::DisplayWarningMessage(const std::string &s) {
   SetScrYellowBoldFont();
   std::cout << "Warning: " << s;
   SetScrNormalFont();
   std::cout << std::endl;
   return;
}
void ScreenUtils::DisplayWarningMessage(const char* word) {
   DisplayWarningMessage(std::string(word));
   return;
}
void ScreenUtils::DisplayGreenMessage(const std::string &s) {
   SetScrGreenBoldFont();
   std::cout << s;
   SetScrNormalFont();
   std::cout << std::endl;
   return;
}
void ScreenUtils::DisplayGreenMessage(const char *word) {
   DisplayGreenMessage(std::string(word));
   return;
}
void ScreenUtils::PrintBetweenStarLines(const std::string &s) {
   PrintScrStarLine();
   CenterString(s);
   PrintScrStarLine();
   return;
}
void ScreenUtils::PrintBetweenStarLines(const char* word) {
   PrintBetweenStarLines(std::string(word));
   return;
}
void ScreenUtils::PrintV3Comp(const double (&v)[3]) {
   for (int i=0; i<3; i++) {std::cout << v[i] << " ";}
   std::cout << std::endl;
   return;
}
void ScreenUtils::PrintV3Comp(const std::string &s,const double (&v)[3]) {
   std::cout << s;
   PrintV3Comp(v);
   return;
}
void ScreenUtils::PrintV3Comp(const char* word,const double (&v)[3]) {
   std::cout << word;
   PrintV3Comp(v);
   return;
}
void ScreenUtils::PrintV3Comp(const int (&v)[3]) {
   for (int i=0; i<3; i++) {std::cout << v[i] << " ";}
   std::cout << std::endl;
   return;
}
void ScreenUtils::PrintV3Comp(const std::string &s,const int (&v)[3]) {
   std::cout << s;
   PrintV3Comp(v);
   return;
}
void ScreenUtils::PrintV3Comp(const char* word,const int (&v)[3]) {
   std::cout << word;
   PrintV3Comp(v);
   return;
}
void ScreenUtils::PrintM3x3Comp(const double (&m)[3][3]) {
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {std::cout << m[i][j] << " ";}
      std::cout << std::endl;
   }
   //std::cout << std::endl;
   return;
}
void ScreenUtils::PrintM3x3Comp(const std::string &s,const double (&m)[3][3]) {
   std::cout << s;
   PrintM3x3Comp(m);
   return;
}
void ScreenUtils::PrintM3x3Comp(const char* word,const double (&m)[3][3]) {
   std::cout << word;
   PrintM3x3Comp(m);
   return;
}
void ScreenUtils::PrintFancyMemoryUsage(size_t memus_,std::string msg) {
   std::string memunit=" B";
   if ( memus_<=1024 ) {std::cout << msg << memus_  << memunit << std::endl; return;}
   double memus=double(memus_)/1024.0e0;
   if ( memus<=(1024.0e0) ) {
      memunit=" KB";
      std::cout << msg << memus << memunit << std::endl;
      return;
   }
   memus/=1024.0e0;
   if ( memus<=(1024.0e0) ) {
      memunit=" MB";
      std::cout << msg << memus << memunit << std::endl;
      return;
   }
   memus/=1024.0e0;
   if ( memus<=(1024.0e0) ) {
      memunit=" GB";
      std::cout << msg << memus << memunit << std::endl;
      return;
   }
   memus/=1024.0e0;
   if ( memus<=(1024.0e0) ) {
      memunit=" TB";
      std::cout << msg << memus << memunit << std::endl;
      return;
   }
   std::cout << "More than 1 PetaByte!" << std::endl;
#if DEBUG
   DISPLAYDEBUGINFOFILELINE;
#endif /* ( DEBUG ) */
}
void ScreenUtils::PrintFancyMemoryUsage(int memus_,std::string msg) {
   if ( memus_<0 ) {DisplayErrorMessage("Memory usage cannot be negative!"); return;}
   else {PrintFancyMemoryUsage(size_t(memus_),msg);}
}
void ScreenUtils::PrintM2x2Comp(const double (&m)[2][2]) {
   for (int i=0; i<2; i++) {
      for (int j=0; j<2; j++) {std::cout << m[i][j] << " ";}
      std::cout << std::endl;
   }
   return;
}
void ScreenUtils::PrintM2x2Comp(const std::string &s,const double (&m)[2][2]) {
   std::cout << s;
   PrintM2x2Comp(m);
   return;
}
void ScreenUtils::PrintM2x2Comp(const char* word,const double (&m)[2][2]) {
   std::cout << word;
   PrintM2x2Comp(m);
   return;
}

