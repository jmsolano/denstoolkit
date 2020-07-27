#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "stringtools.h"

string StringTools::GetEnhancedEpsAtLbl(const string &instr) {
   string outstr=instr;
   size_t pos=outstr.find_last_not_of("0123456789");
   if (pos!=string::npos) {
      outstr.insert(pos+1,"_{");
      outstr.append("}");
   }
   return outstr;
}
string StringTools::GetEnhancedEpsAtLbl(const char *inword) {
   return GetEnhancedEpsAtLbl(string(inword));
}
string StringTools::GetEnhancedEpsTitle(const string &instr) {
   string outstr="";
   for (size_t i=0; i<instr.length(); i++) {
      if (instr[i]=='_') {outstr+='\\';}
      outstr+=instr[i];
   }
   return outstr;
}
void StringTools::RemoveSpacesLeftAndRight(string &str) {
   RemoveSpacesLeft(str);
   RemoveSpacesRight(str);
   return;
}
void StringTools::RemoveSpacesLeft(string &str) {
   while ((str.length()>0)&&(str[0]==' '||str[0]=='\t')) {
      str.erase(0,1);
   }
   return;
}
void StringTools::RemoveSpacesRight(string &str) {
   int len=str.length()-1;
   while (len>=0&&(str[len]==' '||str[len]=='\t')) {
      str.erase(len,1);
      len--;
   }
   return;
}
string StringTools::GetStringFromReal(const double number) {
   std::ostringstream numstr;
   numstr.str("");
   numstr << number;
   return numstr.str();
}
string StringTools::GetStringFromInt(const int number) {
   std::ostringstream numstr;
   numstr.str("");
   numstr << number;
   return numstr.str();
}
string StringTools::GetFilledStringFromInt(const int number,const int width,char filler) {
   std::ostringstream numstr;
   numstr.str("");
   numstr << std::setfill(filler) << std::setw(width) << number;
   return numstr.str();
}
string StringTools::GenerateStrRandSeq(const int len) {
   static const int dm1=61;
   static const char alphanum[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   srand (clock());
   string s="";
   for (int i=0; i<len; ++i) {s+=alphanum[rand()%(dm1)];}
   return s;
}
void StringTools::ReplaceTabsForSpaces(string &s) {
   for (size_t i=0; i<s.length(); i++) {if (s[i]=='\t') {s[i]=' ';}}
   return;
}
void StringTools::RemoveRedundantSpaces(string &s) {
   std::size_t pos=s.find("  ");
   while (pos!=string::npos) {
      s.erase(pos,1);
      pos=s.find("  ");
   }
   return;
}
string StringTools::GetStrFromRealForFileNaming(double number,int prev,int post) {
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
string StringTools::GetFirstChunk(const string &line,char delim) {
   size_t pos=line.find_first_of(delim);
   string res=line.substr(0,pos);
   return res;
}
string StringTools::GetFirstChunkAndDeleteFromLine(string &line,char delim) {
   size_t pos=line.find_first_of(delim);
   if ( pos==string::npos ) {
      return string("");
   }
   string res=line.substr(0,pos);
   line.erase(0,res.length()+1);
   return res;
}
string StringTools::RemoveAllDigits(const string &line) {
   string res=line;
   res.erase(std::remove_if(res.begin(), res.end(), &isdigit),res.end());
   return res;
}


