/*  screenutils.h
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
#ifndef _SCREEN_UTILS_H_
#define _SCREEN_UTILS_H_
#include <iomanip>
#include <iostream>
#ifndef CURRENTVERSION
#define CURRENTVERSION "0.0.0"
#endif

#ifndef PROGRAMCONTRIBUTORS
#define PROGRAMCONTRIBUTORS "JMSA"
#endif

/* ************************************************************************** */
class ScreenUtils {
/* ************************************************************************** */
public:
   ScreenUtils() {}
/* ************************************************************************************** */
   /* ************************************************************************************** */
   static void CenterString(const std::string &s);
   /* ************************************************************************************** */
   static void CenterString(const char* word);
   /* ************************************************************************************** */
   static void PrintProgressBar(int perc);
   /* ************************************************************************************** */
   static void PrintHappyEnding(void);
   /* ************************************************************************************** */
   static bool IsDigit(char c);
   /* ************************************************************************************** */
   static void PrintScrStarLine(void);
   /* ************************************************************************************** */
   static void PrintScrCharLine(char h);
   /* ************************************************************************************** */
   static void SetScrGreenBoldFont(void);
   /* ************************************************************************************** */
   static void SetScrRedBoldFont(void);
   /* ************************************************************************************** */
   static void SetScrYellowBoldFont(void);
   /* ************************************************************************************** */
   static void SetScrBlueBoldFont(void);
   /* ************************************************************************************** */
   static void SetScrBoldFont(void);
   /* ************************************************************************************** */
   static void SetScrNormalFont(void);
   /* ************************************************************************************** */
   static void PrintHappyStart(char ** (&argv),const char *vers,const char *contrib);
   /* ************************************************************************************** */
   static void DisplayErrorMessage(const std::string &s);
   /* ************************************************************************************** */
   static void DisplayErrorMessage(const char* word);
   /* ************************************************************************************** */
   static void DisplayWarningMessage(const std::string &s);
   /* ************************************************************************************** */
   static void DisplayWarningMessage(const char* word);
   /* ************************************************************************************** */
   static void DisplayGreenMessage(const std::string &s);
   /* ************************************************************************************** */
   static void DisplayGreenMessage(const char* word);
   /* ************************************************************************************** */
   static void PrintBetweenStarLines(const std::string &s);
   /* ************************************************************************************** */
   static void PrintBetweenStarLines(const char* word);
   /* ************************************************************************************** */
   static void PrintV3Comp(const double (&v)[3]);
   /* ************************************************************************************** */
   static void PrintV3Comp(const std::string &s,const double (&v)[3]);
   /* ************************************************************************************** */
   static void PrintV3Comp(const char* word,const double (&v)[3]);
   /* ************************************************************************************** */
   static void PrintM3x3Comp(const double (&m)[3][3]);
   /* ************************************************************************************** */
   static void PrintM3x3Comp(const std::string &s,const double (&m)[3][3]);
   /* ************************************************************************************** */
   static void PrintM3x3Comp(const char* word,const double (&m)[3][3]);
   /* ************************************************************************************** */
   static void PrintV3Comp(const int (&v)[3]);
   /* ************************************************************************************** */
   static void PrintV3Comp(const std::string &s,const int (&v)[3]);
   /* ************************************************************************************** */
   static void PrintV3Comp(const char* word,const int (&v)[3]);
   /* ************************************************************************************** */
   static void PrintFancyMemoryUsage(size_t memus_,std::string msg="Memory usage: ");
   /* ************************************************************************************** */
   static void PrintFancyMemoryUsage(int memus_,std::string msg="Memory usage: ");
   /* ************************************************************************************** */
   static void PrintM2x2Comp(const double (&m)[2][2]);
   /* ************************************************************************************** */
   static void PrintM2x2Comp(const std::string &s,const double (&m)[2][2]);
   /* ************************************************************************************** */
   static void PrintM2x2Comp(const char* word,const double (&m)[2][2]);
   /* ************************************************************************************** */
protected:
   /* ************************************************************************** */
};
#endif //_SCREEN_UTILS_H_
