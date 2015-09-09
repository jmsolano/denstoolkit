/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */

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
#ifndef _SOL_SCR_UTIL_H_
#define _SOL_SCR_UTIL_H_
#include <iomanip>
#include <iostream>
#ifndef CURRENTVERSION
#define CURRENTVERSION "0.0.0"
#endif

#ifndef PROGRAMCONTRIBUTORS
#define PROGRAMCONTRIBUTORS "JMSA"
#endif

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

/* ************************************************************************************** */
#ifdef _SOL_USE_FIGLET_NAME_
void printFigletName(void);
#endif
/* ************************************************************************************** */
void centerString(const std::string &s);
/* ************************************************************************************** */
void centerString(const char* word);
/* ************************************************************************************** */
void printProgressBar(int perc);
/* ************************************************************************************** */
void printHappyEnding(void);
/* ************************************************************************************** */
bool isDigit(char c);
/* ************************************************************************************** */
void printScrStarLine(void);
/* ************************************************************************************** */
void printScrCharLine(char h);
/* ************************************************************************************** */
void setScrGreenBoldFont(void);
/* ************************************************************************************** */
void setScrRedBoldFont(void);
/* ************************************************************************************** */
void setScrYellowBoldFont(void);
/* ************************************************************************************** */
void setScrBlueBoldFont(void);
/* ************************************************************************************** */
void setScrBoldFont(void);
/* ************************************************************************************** */
void setScrNormalFont(void);
/* ************************************************************************************** */
void printHappyStart(char ** (&argv),const char *vers,const char *contrib);
/* ************************************************************************************** */
void displayErrorMessage(const std::string &s);
/* ************************************************************************************** */
void displayErrorMessage(const char* word);
/* ************************************************************************************** */
void displayWarningMessage(const std::string &s);
/* ************************************************************************************** */
void displayWarningMessage(const char* word);
/* ************************************************************************************** */
void displayGreenMessage(const std::string &s);
/* ************************************************************************************** */
void displayGreenMessage(const char* word);
/* ************************************************************************************** */
void printBetweenStarLines(const std::string &s);
/* ************************************************************************************** */
void printBetweenStarLines(const char* word);
/* ************************************************************************************** */
void printV3Comp(const solreal (&v)[3]);
/* ************************************************************************************** */
void printV3Comp(const std::string &s,const solreal (&v)[3]);
/* ************************************************************************************** */
void printV3Comp(const char* word,const solreal (&v)[3]);
/* ************************************************************************************** */
void printM3x3Comp(const solreal (&m)[3][3]);
/* ************************************************************************************** */
void printM3x3Comp(const std::string &s,const solreal (&m)[3][3]);
/* ************************************************************************************** */
void printM3x3Comp(const char* word,const solreal (&m)[3][3]);
/* ************************************************************************************** */
void printV3Comp(const int (&v)[3]);
/* ************************************************************************************** */
void printV3Comp(const std::string &s,const int (&v)[3]);
/* ************************************************************************************** */
void printV3Comp(const char* word,const int (&v)[3]);
/* ************************************************************************************** */
void printFancyMemoryUsage(size_t memus_,std::string msg="Memory usage: ");
/* ************************************************************************************** */
void printFancyMemoryUsage(int memus_,std::string msg="Memory usage: ");
/* ************************************************************************************** */
void printM2x2Comp(const solreal (&m)[2][2]);
/* ************************************************************************************** */
void printM2x2Comp(const std::string &s,const solreal (&m)[2][2]);
/* ************************************************************************************** */
void printM2x2Comp(const char* word,const solreal (&m)[2][2]);
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#endif //_SOL_SCR_CUTIL_H_
