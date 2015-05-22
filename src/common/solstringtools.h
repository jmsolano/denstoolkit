/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.1.1
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



#ifndef _SOL_STRING_TOOLS_H_
#define _SOL_STRING_TOOLS_H_

#include <iostream>
#include <string>
using std::string;

/* ******************************************************************************************* */
string getEnhancedEpsAtLbl(const string &instr);
/*
   This function add an underscore and braces around numbers in the string instr.
   The function should receive strings such as "Ca12", "H1", "C32"... etc.
   Different styles may not work...
 */
/* ******************************************************************************************* */
string getEnhancedEpsAtLbl(const char *inword); //Overloaded function, see above.
/* ******************************************************************************************* */
string getEnhancedEpsTitle(const string &instr);
/* This function adds a backslash wherever an underscore appears in instr */
/* ******************************************************************************************* */
void removeSpacesLeftAndRight(string &str); //self descriptive (tab character is considered a space)
/* ******************************************************************************************* */
void removeSpacesLeft(string &str); //self descriptive (tab character is considered a space)
/* ******************************************************************************************* */
void removeSpacesRight(string &str); //self descriptive (tab character is considered a space)
/* ******************************************************************************************* */
string getStringFromReal(const solreal number); //alternative implementation of to_string (C++11)
/* ******************************************************************************************* */
string getStringFromInt(const int number); //alternative implementation of to_string (C++11)
/* ******************************************************************************************* */
string getFilledStringFromInt(const int number,const int width,char filler='0');
/* ******************************************************************************************* */
string genStrRandSeq(const int len);
/* ******************************************************************************************* */
void replaceTabsForSpaces(string &s);
/* ******************************************************************************************* */
void removeRedundantSpaces(string &s);
/* ******************************************************************************************* */
string getStrFromRealForFileNaming(solreal number,int prev=2,int post=2);
/* ******************************************************************************************* */
string getFirstChunk(const string &line,char delim=' ');
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */


#endif /* defined(_SOL_STRING_TOOLS_H_) */
