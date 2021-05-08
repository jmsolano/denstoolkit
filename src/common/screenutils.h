/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
#include <string>

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
   static void CenterString(const std::string &s);
   static void CenterString(const char* word);
   /* ************************************************************************************** */
   static void PrintProgressBar(int perc);
   static void PrintHappyEnding(void);
   static bool IsDigit(char c);
   /* ************************************************************************************** */
   static void PrintScrStarLine(void);
   static void PrintScrCharLine(char h);
   /* ************************************************************************************** */
   static void SetScrGreenBoldFont(void);
   static void SetScrRedBoldFont(void);
   static void SetScrYellowBoldFont(void);
   static void SetScrBlueBoldFont(void);
   static void SetScrBoldFont(void);
   static void SetScrNormalFont(void);
   /* ************************************************************************************** */
   static void PrintHappyStart(char ** (&argv),const char *vers,const char *contrib);
   /* ************************************************************************************** */
   static void DisplayErrorMessage(const std::string &s);
   static void DisplayErrorMessage(const char* word);
   static void DisplayWarningMessage(const std::string &s);
   static void DisplayErrorFileNotOpen(const std::string &s);
   static void DisplayErrorFileNotOpen(const char* word);
   static void DisplayWarningMessage(const char* word);
   static void DisplayGreenMessage(const std::string &s);
   static void DisplayGreenMessage(const char* word);
   /* ************************************************************************************** */
   static void PrintBetweenStarLines(const std::string &s);
   static void PrintBetweenStarLines(const char* word);
   /* ************************************************************************************** */
   static void PrintV3Comp(const double (&v)[3]);
   static void PrintV3Comp(const std::string &s,const double (&v)[3]);
   static void PrintV3Comp(const char* word,const double (&v)[3]);
   static void PrintM3x3Comp(const double (&m)[3][3]);
   static void PrintM3x3Comp(const std::string &s,const double (&m)[3][3]);
   static void PrintM3x3Comp(const char* word,const double (&m)[3][3]);
   static void PrintV3Comp(const int (&v)[3]);
   static void PrintV3Comp(const std::string &s,const int (&v)[3]);
   static void PrintV3Comp(const char* word,const int (&v)[3]);
   static void PrintM2x2Comp(const double (&m)[2][2]);
   static void PrintM2x2Comp(const std::string &s,const double (&m)[2][2]);
   static void PrintM2x2Comp(const char* word,const double (&m)[2][2]);
   /* ************************************************************************************** */
   static void PrintFancyMemoryUsage(size_t memus_,std::string msg="Memory usage: ");
   static void PrintFancyMemoryUsage(int memus_,std::string msg="Memory usage: ");
   /* ************************************************************************************** */
protected:
   /* ************************************************************************** */
};
#endif //_SCREEN_UTILS_H_
