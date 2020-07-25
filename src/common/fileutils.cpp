#ifndef _FILEUTILS_CPP_
#define _FILEUTILS_CPP_
#include <cstdlib>
#include <cstdio>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ifstream;
#include "fileutils.h"
#include "screenutils.h"
#include "stringtools.h"

#ifndef DEBUG
#define DEBUG 0
#endif

FileUtils::FileUtils() {
}
int FileUtils::CountNumberOfNonCommentLines(const string &fname) {
   ifstream ifil(fname.c_str());
   if ( !ifil.good() ) {
       cout << "Error while openning the file " << fname << endl;
       ifil.close();
       return -1;
   }
   DiscardComments(ifil);
   string line;
   int nl=0;
   while (!ifil.eof()) {
      std::getline(ifil,line);
      if ( line.length()>0 ) { ++nl; }
   }
   ifil.close();
   return nl;
}
int FileUtils::CountNumberOfLines(const string &fname) {
   FILE* ifil = std::fopen(fname.c_str(), "r");
   if(!ifil) {
#if DEBUG
      cout << "Error while openning file " << fname << endl;
#endif /* ( DEBUG ) */
      std::fclose(ifil);
      return -1;
   }

   int c; // note: int, not char, required to handle EOF
   int nl=0;
   while ((c = std::fgetc(ifil)) != EOF) { // standard C I/O file reading loop
      if ( c=='\n' ) { ++nl; }
   }
   std::fclose(ifil);
   return nl;
}
int FileUtils::CountNumberOfLinesThatContainsString(const string &fname,const string &str) {
   ifstream ifil(fname.c_str());
   if(!ifil.good()) {
#if DEBUG
      cout << "Error while openning file " << fname << endl;
#endif /* ( DEBUG ) */
      ifil.close();
      return -1;
   }
   string line;
   int nl=0;
   while (!ifil.eof()) {
      std::getline(ifil,line);
      if ( line.length()>0 ) {
         if ( line.find(str)!=std::string::npos ) {
            ++nl;
         }
      }
   }
   ifil.close();
   return nl;
}
int FileUtils::CountNumberOfLinesInTSV(const string &fname) {
   ifstream ifil(fname.c_str());
   if(!ifil.good()) {
#if DEBUG
      cout << "Error while openning file " << fname << endl;
#endif /* ( DEBUG ) */
      ifil.close();
      return -1;
   }
   string line;
   int nl=0;
   while (!ifil.eof()) {
      std::getline(ifil,line);
      if ( line.length()>0 ) { ++nl; }
   }
   ifil.close();
   return nl;
}
int FileUtils::GetBlockLengthInTSV(const string &fname) {
   ifstream ifil(fname.c_str());
   if(!ifil.good()) {
#if DEBUG
      cout << "Error while openning file " << fname << endl;
#endif /* ( DEBUG ) */
      ifil.close();
      return -1;
   }
   string line;
   int nb=0;
   DiscardComments(ifil);
   while (!ifil.eof()) {
      std::getline(ifil,line);
      if ( line.length()>0 ) { ++nb; } else {break;}
   }
   ifil.close();
   return nb;
}
void FileUtils::GenerateRandomTmpFileName(string &s, const int len) {
   static const int dm1=61;
   static const char alphanum[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   srand (clock());
   s="";
   for (int i=0; i<len; ++i) {s+=alphanum[rand()%(dm1)];}
}
int FileUtils::CountColumnsInFile(const string &fnam) {
   string line;
   StringTools st;
   ifstream ifil;
   ifil.open(fnam.c_str(),std::ios::in);
   DiscardComments(ifil);
   getline(ifil,line);
   ifil.close();
   st.ReplaceTabsForSpaces(line);
   st.RemoveSpacesLeft(line);
   line.append(" ");
   st.RemoveRedundantSpaces(line);
   int numtabs=0;
   for (size_t i=0; i<line.length(); i++) {if (line[i]==' ') {numtabs++;}}
   return numtabs;
}
void FileUtils::DiscardComments(ifstream &ifil) {
   std::string tl;
   int prevpos;
   bool procnext=true;
   while (procnext) {
      prevpos=ifil.tellg();
      getline(ifil,tl);
      while (tl.length()==0) {
         getline(ifil,tl);
         if (ifil.eof()) {return;}
      }
      while((tl[0]==' ')||(tl[0]=='\t')){tl.erase(0,1);}
      if (tl[0]!='#') {
         procnext=false;
         ifil.seekg(prevpos);
      }
   }
}
void FileUtils::ReadXYColumns(string fname,vector<double> &x,vector<double> &y,\
      int nx,int ny) {
   int nr=CountNumberOfNonCommentLines(fname);
   int nc=CountColumnsInFile(fname);
   vector<vector<double> > m(nr);
   for ( int i=0 ; i<nr ; ++i ) {
      m[i].resize(nc);
   }
   ifstream ifil(fname.c_str());
   if ( !ifil.good() ) {
       cout << "Error while openning the file " << fname << endl;
       ifil.close();
       return;
   }
   DiscardComments(ifil);
   for ( int i=0 ; i<nr ; ++i ) {
      for ( int j=0 ; j<nc ; ++j ) {
         ifil >> m[i][j];
      }
   }
   ifil.close();
   if ( int(x.size())!=nr ) { x.resize(nr); }
   if ( int(y.size())!=nr ) { y.resize(nr); }
   for ( int i=0 ; i<nr ; ++i ) {
      x[i]=m[i][nx-1];
      y[i]=m[i][ny-1];
   }
}
vector<double> FileUtils::ReadSingleColumn(const string &fname,\
      int nocol) {
   int nr=CountNumberOfNonCommentLines(fname);
   int nc=CountColumnsInFile(fname);
   vector<vector<double> > m(nr);
   if ( nocol<1 || nocol>nc ) {
      ScreenUtils::DisplayErrorMessage("Invalid column number");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
   }
   for ( int i=0 ; i<nr ; ++i ) { m[i].resize(nc); }
   ifstream ifil(fname.c_str());
   vector<double> col(0);
   if ( !ifil.good() ) {
       cout << "Error while openning the file " << fname << endl;
       ifil.close();
       return col;
   }
   DiscardComments(ifil);
   for ( int i=0 ; i<nr ; ++i ) {
      for ( int j=0 ; j<nc ; ++j ) {
         ifil >> m[i][j];
      }
   }
   ifil.close();
   col.resize(nr);
   for ( int i=0 ; i<nr ; ++i ) { col[i]=m[i][nocol-1]; }
   return col;
}
vector<vector<double> > FileUtils::ReadDataMultiColumsFromFile(string fname,\
         size_t colx,size_t coly,double minx,double maxx,double miny,double maxy) {
   vector<vector<double> > data(0);
   size_t cols=size_t(CountColumnsInFile(fname));
   size_t rows=size_t(CountNumberOfNonCommentLines(fname));
   vector<double> tmp(cols);
   ifstream ifil(fname.c_str());
   if ( !ifil.good() ) {
      ifil.close();
      ScreenUtils::DisplayErrorMessage(string("Could not open the file '")+fname+string("'!"));
      return data;
   }
   DiscardComments(ifil);
   //Reads the data
   for ( size_t i=0 ; i<rows ; ++i ) {
      for ( size_t j=0 ; j<cols ; ++j ) { ifil >> tmp[j]; }
      if ( (tmp[colx-1]<= maxx) && (tmp[colx-1]>=minx) && (tmp[coly]<=maxy) && (tmp[coly]>=miny) ) {
         data.push_back(tmp);
      }
   }
   ifil.close();
   return data;
}
void FileUtils::WriteCenteredString(ofstream &ofil,const string &s,bool comm) {
   int pos=int((80-s.length())/2);
   if (comm) { ofil << "#"; }
   for( int i=0; i<pos; ++i ){ ofil << " "; }
   ofil << s << endl;
}
void FileUtils::WriteScrCharLine(ofstream &ofil,const char h,bool comm) {
   if ( comm ) { ofil << "#"; }
   for ( int i=0 ; i<80 ; ++i ) { ofil << h; }
   ofil << endl;
}
void FileUtils::WriteHappyStart(char ** (&argv),ofstream &ofil,string vers,string progcont,bool comm) {
   string progname=argv[0];
   progname.append("  (-:");
   progname.insert(0,":-) ");
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}
   WriteScrStarLine(ofil,comm);
   if (comm) { ofil << "#"; }
   ofil << endl;
   WriteCenteredString(ofil,progname,comm);
   if (comm) { ofil << "#"; }
   ofil << endl;
   progname=__DATE__;
   progname.insert(0,"Compilation date: ");
   WriteCenteredString(ofil,progname,comm);
   if (comm) { ofil << "#"; }
   ofil << endl;
   progname="Version: ";
   progname.append(vers);
   WriteCenteredString(ofil,progname,comm);
   if (comm) { ofil << "#"; }
   ofil << endl;
   progname=":-)    Created by ";
   progname.append(progcont);
   progname.append("    (-:");
   WriteCenteredString(ofil,progname,comm);
   if (comm) { ofil << "#"; }
   ofil << endl;
   WriteScrStarLine(ofil,comm);
}
void FileUtils::WriteV3Components(ofstream &ofil,const double (&v)[3]) {
   for (int i=0; i<3; i++) {ofil << v[i] << " ";}
   ofil << endl;
}
void FileUtils::WriteV3Components(ofstream &ofil,const std::string &s,const double (&v)[3]){
   ofil << s;
   FileUtils::WriteV3Components(ofil,v);
}
void FileUtils::ReplaceExtensionOfFileName(string &orig,\
      const string thenewext) {
   size_t pos=orig.find_last_of('.')+1;
   orig=orig.substr(0,pos)+thenewext;
}
void FileUtils::InsertAtEndOfFileName(string &orig,const string str2insrt) {
   size_t pos=orig.find_last_of('.');
   orig.insert(pos,str2insrt);
}
bool FileUtils::ExtensionMatches(const string &fname,const string ext) {
   size_t nExt=ext.size();
   size_t pos=fname.find_last_of('.')+1;
   return (fname.substr(pos,nExt)==ext);
}


#endif  /* _FILEUTILS_CPP_ */

