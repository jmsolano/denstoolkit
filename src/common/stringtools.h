#ifndef _STRING_TOOLS_H_
#define _STRING_TOOLS_H_

#include <iostream>
#include <string>
using std::string;

/* ************************************************************************** */
class StringTools {
/* ************************************************************************** */
public:
/* ******************************************************************************************* */
   /** This function add an underscore and braces around numbers in the string instr.
      The function should receive strings such as "Ca12", "H1", "C32"... etc.
      Different styles may not work... */
   static string GetEnhancedEpsAtLbl(const string &instr);
   /* ******************************************************************************************* */
   static string GetEnhancedEpsAtLbl(const char *inword); //Overloaded function, see above.
   /* ******************************************************************************************* */
   static string GetEnhancedEpsTitle(const string &instr);
   /* This function adds a backslash wherever an underscore appears in instr */
   /* ******************************************************************************************* */
   static void RemoveSpacesLeftAndRight(string &str); //self descriptive (tab character is considered a space)
   /* ******************************************************************************************* */
   static void RemoveSpacesLeft(string &str); //self descriptive (tab character is considered a space)
   /* ******************************************************************************************* */
   static void RemoveSpacesRight(string &str); //self descriptive (tab character is considered a space)
   /* ******************************************************************************************* */
   static string GetStringFromReal(const double number); //alternative implementation of to_string (C++11)
   /* ******************************************************************************************* */
   static string GetStringFromInt(const int number); //alternative implementation of to_string (C++11)
   /* ******************************************************************************************* */
   static string GetFilledStringFromInt(const int number,const int width,char filler='0');
   /* ******************************************************************************************* */
   static string GenerateStrRandSeq(const int len);
   /* ******************************************************************************************* */
   static void ReplaceTabsForSpaces(string &s);
   /* ******************************************************************************************* */
   static void RemoveRedundantSpaces(string &s);
   /* ******************************************************************************************* */
   static string GetStrFromRealForFileNaming(double number,int prev=2,int post=2);
   /* ******************************************************************************************* */
   static string GetFirstChunk(const string &line,char delim=' ');
   /* ******************************************************************************************* */
   static string GetFirstChunkAndDeleteFromLine(string &line,char delim=' ');
   /* ************************************************************************** */
   static string RemoveAllDigits(const string &line);
   /* ************************************************************************** */
protected:
   /* ************************************************************************** */
};
/* ************************************************************************** */

/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */


#endif /* defined(_STRING_TOOLS_H_) */
