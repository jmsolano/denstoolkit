

#ifndef _SOL_STRING_TOOLS_CPP_
#define _SOL_STRING_TOOLS_CPP_

#include "solstringtools.h"
#include <sstream>
#include <iomanip>
#include <cstdlib>


/* ******************************************************************************************* */
string getEnhancedEpsAtLbl(const string &instr)
{
   string outstr=instr;
   size_t pos=outstr.find_last_not_of("0123456789");
   if (pos!=string::npos) {
      outstr.insert(pos+1,"_{");
      outstr.append("}");
   }
   return outstr;
}
/* ******************************************************************************************* */
string getEnhancedEpsAtLbl(const char *inword)
{
   return getEnhancedEpsAtLbl(string(inword));
}
/* ******************************************************************************************* */
string getEnhancedEpsTitle(const string &instr)
{
   string outstr="";
   for (size_t i=0; i<instr.length(); i++) {
      if (instr[i]=='_') {outstr+='\\';}
      outstr+=instr[i];
   }
   return outstr;
}
/* ******************************************************************************************* */
void removeSpacesLeftAndRight(string &str)
{
   removeSpacesLeft(str);
   removeSpacesRight(str);
   return;
}
/* ******************************************************************************************* */
void removeSpacesLeft(string &str)
{
   while ((str.length()>0)&&(str[0]==' '||str[0]=='\t')) {
      str.erase(0,1);
   }
   return;
}
/* ******************************************************************************************* */
void removeSpacesRight(string &str)
{
   int len=str.length()-1;
   while (len>=0&&(str[len]==' '||str[len]=='\t')) {
      str.erase(len,1);
      len--;
   }
   return;
}
/* ******************************************************************************************* */
string getStringFromReal(const solreal number)
{
   std::ostringstream numstr;
   numstr.str("");
   numstr << number;
   return numstr.str();
}
/* ******************************************************************************************* */
string getStringFromInt(const int number)
{
   std::ostringstream numstr;
   numstr.str("");
   numstr << number;
   return numstr.str();
}/* ******************************************************************************************* */
string getFilledStringFromInt(const int number,const int width,char filler)
{
   std::ostringstream numstr;
   numstr.str("");
   numstr << std::setfill(filler) << std::setw(width) << number;
   return numstr.str();
}
/* ******************************************************************************************* */
string genStrRandSeq(const int len)
{
   static const int dm1=61;
   static const char alphanum[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   srand (clock());
   string s="";
   for (int i=0; i<len; ++i) {s+=alphanum[rand()%(dm1)];}
   return s;
}
/* ******************************************************************************************* */
void replaceTabsForSpaces(string &s)
{
   for (int i=0; i<s.length(); i++) {if (s[i]=='\t') {s[i]=' ';}}
   return;
}
/* ******************************************************************************************* */
void removeRedundantSpaces(string &s)
{
   std::size_t pos=s.find("  ");
   while (pos!=string::npos) {
      s.erase(pos,1);
      pos=s.find("  ");
   }
   return;
}
/* ******************************************************************************************* */
string getStrFromRealForFileNaming(solreal number,int prev,int post)
{
   std::ostringstream numstr;
   numstr.str("");
   numstr.setf(std::ios::fixed);
   numstr << std::setw(prev+post+1) << 
      std::setprecision(post) << std::setfill('0') << number;
   std::string finstr=numstr.str();
   size_t pos=finstr.find_first_of('.');
   if ( pos!=std::string::npos ) {finstr[pos]='p';}
   return finstr;
}
/* ******************************************************************************************* */
string getFirstChunk(const string &line,char delim)
{
   size_t pos=line.find_first_of(delim);
   string res=line.substr(0,pos);
   return res;
}
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */

#endif /* defined(_SOL_STRING_TOOLS_CPP_) */

