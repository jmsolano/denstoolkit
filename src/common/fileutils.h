#ifndef _FILEUTILS_H_
#define _FILEUTILS_H_
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <vector>
using std::vector;

/* ************************************************************************** */
class FileUtils {
/* ************************************************************************** */
public:
   FileUtils();
/* ************************************************************************** */
   static int CountNumberOfLines(const string &fname);
   static int CountNumberOfNonCommentLines(const string &fname);
   static int CountNumberOfLinesThatContainsString(const string &fname,const string &str);
   static int CountNumberOfLinesInTSV(const string &fname);
   static int GetBlockLengthInTSV(const string &fname);
/* ************************************************************************** */
   static int CountColumnsInFile(const string &fnam);
/* ************************************************************************** */
   static void GenerateRandomTmpFileName(string &s, const int len);
/* ************************************************************************** */
   static void DiscardComments(ifstream &ifil);
/* ************************************************************************** */
   /** Reads a single column of a file.  */
   static vector<double> ReadSingleColumn(const string &fname,int nocol=1);
   /** As the name suggests, reads two columns (from a dat file).
    * The column numbers are given through nx and ny (warning: the
    * first column is 1!)  */
   static void ReadXYColumns(string fname,vector<double> &x,vector<double> &y,int nx=1,int ny=2);
   /** Reads data from a dat file. The data can be selected using range
    * criteria. The ranges can be selected up to two variables (or columns).
    * The columns that will be used for selecting the data are given through
    * colx and coly (warning: the first column is 1!).
    * The ranges are given as [minx,maxx] and [miny,maxy]
    * (square brackets here indicate closed sets).*/
   static vector<vector<double> > ReadDataMultiColumsFromFile(string fname,\
         size_t colx=0,size_t coly=1,\
         double minx=-1.0e+50,double maxx=1.0e+50,\
         double miny=-1.0e+50,double maxy=1.0e+50);
/* ************************************************************************** */
   static void WriteCenteredString(ofstream &ofil,const string &s,bool comm=true);
   static void WriteCenteredString(ofstream &ofil,const char* word,bool comm=true) {
      return WriteCenteredString(ofil,string(word),comm);
   }
   static void WriteScrCharLine(ofstream &ofil,const char h,bool comm=true);
   static void WriteScrStarLine(ofstream &ofil,bool comm=true) {return WriteScrCharLine(ofil,'*',comm);}
   static void WriteHappyStart(char ** (&argv),ofstream &ofil,string vers,string progcont,bool comm=true);
/* ************************************************************************** */
   static void WriteV3Components(ofstream &ofil,const double (&v)[3]);
   static void WriteV3Components(ofstream &ofil,const std::string &s,const double (&v)[3]);
   /** Writes the matrix into the file named fname.
    * sep is the character used for separating columns.
    * scient is a boolean to write the data using scientific notation.
    * for separating the matrix values.
    * hdr is a string that will be printed
    * out at the beginning of the file if it is not empty. Warning: The function
    * prints a '#' character for the first line of hdr, however, if the hdr contains
    * newlines, the user must also include the respective '#' after each new line.  */
   static void SaveMatrix(const string &fname,const vector<vector<double> > &mat,\
         const string &hdr="",const bool scient=true,const char sep=' ');
/* ************************************************************************** */
   static void ReplaceExtensionOfFileName(string &orig,const string thenewext);
   static void RemoveExtensionFromFileName(string &str);
/* ************************************************************************** */
   static void InsertAtEndOfFileName(string &orig,const string str2insrt);
/* ************************************************************************** */
   static bool ExtensionMatches(const string &fname,const string ext);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _FILEUTILS_H_ */

