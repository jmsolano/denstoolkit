/*
   This file contains a lot of subroutines, for example to change the extension of a file-name,
   or to add some characters to the name of a file, etc.
   
   Created by: Juan Manuel Solano Altamirano
               Present institution at the moment of the creation of this file:
               Centro de Investigaciones y Estudios Avanzados del IPN,
               Unidad Monterrey.
               México.
               email: jmsolanoalt@gmail.com (prefered)
                      jmsolano@venus.ifuap.buap.mx (alternative)
               
   Created:    04 April 2010
*/

/*
  
  int countchofast(char * str) returns the number of characters in a string.
  int fnamesize(char *filename) returns the size of a file-name, without extension
         nor the period.
  int fextsize (char * filename) returns the size of the extension of a file-name
         without the period.
  void addext (char * origname, char *ext) add the extension to a name (also adds the '.')
  void chgext (char * origname, char * newext, char * newname) changes the extension
         of a file-name into another. Further details below.
  void addch2endname(char * origname, char * sttoadd, char * newname) adds the
       string "sttoadd" to the end of a file-name (and before the extension).
  void addch2begname(char * origname, char * sttoadd, char * newname) adds the
       string "sttoadd" to the begining of a file-name.
  int countcolumns (char *filename) counts the number of columns in the data-file
       named "filename". It assumes that the file ONLY uses the tab character
       in between the columns and not at the end of the row (line).
  int countrows (char * filename) counts the number of rows in the data-file
       named "filename". It assumes that there is a new-line character
       right after the last number (before the EOF character).
  
  
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  addext() adds the extension of a file-name. It receives a char* (the file-name
  with no extension) and a char* (the extension). This function also adds the
  character '.' to the name, for example(if name="filename" and ext="txt"): 
            addext(name,ext)
  then the name string is modified such that after this function takes the value
  name="filename.txt"
  
  chgext () changes the extension of a file-name. It receives a char*
  (the file-name, including the old extension), a char * (the new extension)
  and a char * (where the new name with the new extension will be put). It copies the name
  but changes its extension. NOTE: this routine only supports names up to 128 characters
  (included the extension)
  
  addch2endname () adds a series of characters to the name of a file (it keeps its 
  root name and extension). It adds the characters to the end of the name.
  
  countchofast (char *str) count the number of characters that are in the string whose 
  pointer is char *str
*/

#ifndef _SOLFILEUTILS_H
#define _SOLFILEUTILS_H

//# include <iostream.h>
# include <stdio.h>
# define MAX_ROWS (3000)
# define MAX_COLS (3000)
#include <cstdlib>
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;

//**************************************************************************************************
int countchofast(char * str);
//**************************************************************************************************
int fnamesize(char *filename);
//**************************************************************************************************
int fextsize (char * filename);
//**************************************************************************************************
void addext (char* origname, char* ext);
//**************************************************************************************************
void chgext(char * origname, char * newext, char * newname);
//**************************************************************************************************
void addch2endname(char * origname, char * sttoadd, char * newname);
//**************************************************************************************************
void addch2begname(char * origname, char * sttoadd, char * newname);
//**************************************************************************************************
int countcolumns (char *filename);
//**************************************************************************************************
int countrows (char * filename);
//**************************************************************************************************
void discardComments(ifstream &ifil);
//**************************************************************************************************
void writeHappyStart(char ** (&argv),ofstream &ofil,string vers,string progcont);
//**************************************************************************************************
void writeCommentedHappyStart(char ** (&argv),ofstream &ofil,string vers,string progcont);
//**************************************************************************************************
void writeScrStarLine(ofstream &ofil);
//**************************************************************************************************
void writeCommentedScrStarLine(ofstream &ofil);
//**************************************************************************************************
void writeScrCharLine(ofstream &ofil,char h);
//**************************************************************************************************
void writeCommentedScrCharLine(ofstream &ofil,char h);
//**************************************************************************************************
void centerString(const string &s,ofstream &ofil);
//**************************************************************************************************
void centerCommentedString(const string &s,ofstream &ofil);
//**************************************************************************************************
void centerString(const char* word,ofstream &ofil);
//**************************************************************************************************
void centerCommentedString(const char* word,ofstream &ofil);
//**************************************************************************************************
void writeV3Comp(ofstream &ofil,const solreal (&v)[3]);
//**************************************************************************************************
void writeV3Comp(ofstream &ofil,const std::string &s,const solreal (&v)[3]);
//**************************************************************************************************
void writeV3Comp(ofstream &ofil,const char* word,const solreal (&v)[3]);
//**************************************************************************************************
void writeM3x3Comp(ofstream &ofil,const solreal (&m)[3][3]);
//**************************************************************************************************
void writeM3x3Comp(ofstream &ofil,const std::string &s,const solreal (&m)[3][3]);
//**************************************************************************************************
void writeM3x3Comp(ofstream &ofil,const char* word,const solreal (&m)[3][3]);
//**************************************************************************************************
void genRandomTmpFileName(string &s, const int len);
//**************************************************************************************************
int countColumnsInFile(string fnam);
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
# endif //_SOLFILEUTILS_H
