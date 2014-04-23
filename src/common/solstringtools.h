//
//  solstringtools.h
//  
//
//  Created by Juan Manuel Solano on 2013-10-06.
//
//

#ifndef _SOL_STRING_TOOLS_H_
#define _SOL_STRING_TOOLS_H_

#include <iostream>
#include <string>
using std::string;

//**************************************************************************************************
string getEnhancedEpsAtLbl(const string &instr);
/*
   This function add an underscore and braces around numbers in the string instr.
   The function should receive strings such as "Ca12", "H1", "C32"... etc.
   Different styles may not work...
 */
//**************************************************************************************************
string getEnhancedEpsAtLbl(const char *inword); //Overloaded function, see above.
//**************************************************************************************************
string getEnhancedEpsTitle(const string &instr);
/* This function adds a backslash wherever an underscore appears in instr */
//**************************************************************************************************
void removeSpacesLeftAndRight(string &str); //self descriptive (tab character is considered a space)
//**************************************************************************************************
void removeSpacesLeft(string &str); //self descriptive (tab character is considered a space)
//**************************************************************************************************
void removeSpacesRight(string &str); //self descriptive (tab character is considered a space)
//**************************************************************************************************
string getStringFromReal(const solreal number); //alternative implementation of to_string (C++11)
//**************************************************************************************************
string getStringFromInt(const int number); //alternative implementation of to_string (C++11)
//**************************************************************************************************
string genStrRandSeq(const int len);
//**************************************************************************************************
void replaceTabsForSpaces(string &s);
//**************************************************************************************************
void removeRedundantSpaces(string &s);
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************


#endif /* defined(_SOL_STRING_TOOLS_H_) */
