/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.0
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
   This file contains a lot of subroutines, for example to change the extension of a file-name,
   or to add some characters to the name of a file, etc.
   
   Created by: Juan Manuel Solano Altamirano
               Present institution at the moment of the creation of this file:
               Centro de Investigaciones y Estudios Avanzados del IPN,
               Unidad Monterrey.
               Mexico.
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

#ifndef _SOLFILEUTILS_CPP_
#define _SOLFILEUTILS_CPP_

#include "solfileutils.h"
#include "solscrutils.h"
#include "solstringtools.h"
# include <stdio.h>
# define MAX_ROWS (3000)
# define MAX_COLS (3000)
#include <cstdlib>
#include <string>


/* ******************************************************************************************* */
int countchofast(char * str){
   
   int size=0;
   char chtemp='a';
   
   while (chtemp!='\0') {
      //printf("\n%d",size);
      chtemp=str[size];
      //printf("%c",str[size]);
      size++;
      if (size==2000) {
         chtemp='\0';
         printf("\nWarning: size is larger than 2000. Possible error.");
      }
   }
   return (size-1); //size counts the '\0' character. Must be discounted.
}
/* ******************************************************************************************* */
int fnamesize(char *filename) {

   int size=0;
   int i=0;
   char chtemp='a';
   
   while (chtemp != '\0') {
      chtemp=filename[i];
      //printf("%c",chtemp);
      if (chtemp == '.') {
         size=i;
         //printf("\nsize= %d\n", size);
      }
      i++;
      //printf("%d",i);
      if (i==128) {
         chtemp='\0';
         printf("\nOnly names of 128 characters are permited.\n");
         //system("pause");
         std::exit(1);
      }
   }
   return size;

}
/* ******************************************************************************************* */
int fextsize (char * filename) {

   int extsize=0;
   int namesize=0;
   int i=0;
   char chtemp='a';
   
   namesize=countchofast(filename);
   //printf("namesize: %d",namesize);
   
   while (chtemp!='.') {
      chtemp=filename[namesize-i];
      if (chtemp=='.') {
         extsize=i;
      }
      i++;
   }
   //printf("extsize= %d",(extsize-1));
   return (extsize-1); //extsize includes the '.' character. Must be removed.
}
/* ******************************************************************************************* */
void addext (char* origname, char* ext)
{
   int namesize=0,extsize=0;
   int i;
   
   namesize=countchofast(origname);
   extsize=countchofast(ext);
   
   if ((namesize+extsize)>127){
      printf("\nThe new name is larger than 128 characters long...\n");
      std::exit(1);
   }
   
   origname[namesize]='.';
   for(i=0; i<extsize; i++){
      origname[namesize+i+1]=ext[i];
   }
   
   origname[namesize+extsize+1]='\0';
   
}
/* ******************************************************************************************* */
void chgext(char * origname, char * newext, char * newname){

   int NameSize,ExtSize,i;
   char chtemp;
   
   chtemp='a';
   NameSize=0;
   i=0;
   
   //printf("\nOrigname (fh): %s",origname);
   //printf("\nNewext (fh): %s", newext);
   
   while (chtemp != '\0') {
      chtemp=origname[i];
      //printf("%c",chtemp);
      if (chtemp == '.') {
         NameSize=i;
         //printf("\nNameSize= %d\n", NameSize);
      }
      i++;
      //printf("%d",i);
      if (i==124) {
         chtemp='\0';
         printf("\nOnly names of 128 characters are permited.\n");
         //system("pause");
         std::exit(1);
      }
   }
   ExtSize=i-NameSize-2;
   //printf("ExtSize: %d",ExtSize);
   
   for (i=0; i<NameSize; i++) {
      newname[i]=origname[i];
   }
   newname[NameSize]='.';
   for(i=0; (i<(NameSize+2+ExtSize)); i++) {
      newname[NameSize+1+i]=newext[i];
   }
   newname[NameSize+1+ExtSize]='\0';
   //printf("\nNew name (from .h): %s",newname);
}

/* ******************************************************************************************* */
void addch2endname(char * origname, char * sttoadd, char * newname){

   int i,CompNameSize,NameSize,ExtSize,sttoaddSize;
   
   CompNameSize=countchofast(origname);
   NameSize=fnamesize(origname);
   ExtSize=fextsize(origname);
   sttoaddSize=countchofast(sttoadd);
   
   for (i=0; i<NameSize; i++) {
      newname[i]=origname[i];
   }
   
   for (i=0; i<sttoaddSize; i++) {
      newname[NameSize+i]=sttoadd[i];
   }
   newname[NameSize+sttoaddSize]='.';
   for (i=0; i<ExtSize; i++) {
      newname[NameSize+sttoaddSize+1+i]=origname[NameSize+1+i];
   }
   newname[CompNameSize+sttoaddSize]='\0';
}
/* ******************************************************************************************* */
void addch2begname(char * origname, char * sttoadd, char * newname){

   int i,CompNameSize;
   int sttoaddSize;
   
   CompNameSize=countchofast(origname);
   sttoaddSize=countchofast(sttoadd);
   for (i=0; i<sttoaddSize; i++) {
      newname[i]=sttoadd[i];
   }
   for (i=0; i<CompNameSize; i++) {
      newname[sttoaddSize+i]=origname[i];
   }
   newname[CompNameSize+sttoaddSize]='\0';
}
/* ******************************************************************************************* */
int countcolumns (char *filename){
   
   //int j=0;
   int numtabs=0;
   FILE *datafile;
   char chtemp='a';
   
   if ((datafile=fopen(filename,"r"))==NULL) {
      printf("Error opening the file: %s\n",filename);
      //system("pause");
      std::exit (1);
   }
   
   while ((chtemp!='\n')&&(chtemp!='\r')&&(numtabs<MAX_COLS)&&(chtemp!=EOF)) {
         //if the OS is windows '\n' is defined as the TWO
         //characters '0D 0A' (in hex format), so one cannot count
         // '\n'-characters. Instead, one can count the '\r' characters.
         // If is unix-linux, '\r' chars are not used and can
         // be counted the '\n' (defined as 0A) chars.
      chtemp=fgetc(datafile);
      //j++;
      if ((chtemp == '\t') || (chtemp == ' ') || (chtemp == ',')) {
         numtabs++;
      }
      if(numtabs==MAX_COLS) {
         printf("\nThe number of columns is greater than %d...\n",MAX_COLS);
         printf("Terminating program...");
         fclose(datafile);
         std::exit(1);
      }
   }
   fclose(datafile);
   return (numtabs+1); //it's assumed that there is no '\t' character after
                       //the last column.

}
/* ******************************************************************************************* */
int countrows (char * filename) {
   
   int numnewlines=0;
   FILE *datafile;
   char chtemp='a';
   
   if ((datafile=fopen(filename,"r"))==NULL) {
      printf("Error opening the file: %s\n",filename);
      //system("pause");
      std::exit (1);
   }
   
   while ((numnewlines<MAX_ROWS)&&(chtemp != EOF)) {
      chtemp=fgetc(datafile);
      if ((chtemp == '\n') || (chtemp == '\r') ) { //if is running on windows '\n' is defined as the TWO
                                                   //characters '0D 0A' (in hex format), so one cannot count
                                                   // '\n'-characters. Instead, one can count the '\r' characters.
                                                   // If is running on unix-linux, '\r' chars are not used and can
                                                   // be counted the '\n' (defined as 0A) chars.
         numnewlines++;
      }
      if(numnewlines==MAX_ROWS) {
         printf("The number of rows is greater than %d...\n",MAX_ROWS);
         printf("Terminating program...\n");
         fclose(datafile);
         std::exit (1);
      }
   }
   fclose(datafile);
   return numnewlines;

}
/* ******************************************************************************************* */
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
/* ******************************************************************************************* */
void discardComments(ifstream &ifil)
{
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
      //cout << tl << endl;
   }
}
/* ******************************************************************************************* */
void writeHappyStart(char ** (&argv),ofstream &ofil,string vers,string progcont)
{
   string progname=argv[0];
   progname.append("  (-:");
   progname.insert(0,":-) ");
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}
   writeScrStarLine(ofil);
   ofil << endl;
   centerString(progname,ofil);
   ofil << endl;
   progname=__DATE__;
   progname.insert(0,"Compilation date: ");
   centerString(progname,ofil);
   ofil << endl;
   progname="Version: ";
   progname.append(vers);
   centerString(progname,ofil);
   ofil << endl;
   progname=":-)    Created by ";
   progname.append(progcont);
   progname.append("    (-:");
   centerString(progname,ofil);
   ofil << endl;
   writeScrStarLine(ofil);
   return;
}
/* ******************************************************************************************* */
void writeCommentedHappyStart(char ** (&argv),ofstream &ofil,string vers,string progcont)
{
   string progname=argv[0];
   progname.append("  (-:");
   progname.insert(0,":-) ");
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}
   writeCommentedScrStarLine(ofil);
   ofil << "#" << endl;
   centerCommentedString(progname,ofil);
   ofil << "#" << endl;
   progname=__DATE__;
   progname.insert(0,"Compilation date: ");
   centerCommentedString(progname,ofil);
   ofil << "#" << endl;
   progname="Version: ";
   progname.append(vers);
   centerCommentedString(progname,ofil);
   ofil << "#" << endl;
   progname=":-)    Created by ";
   progname.append(progcont);
   progname.append("    (-:");
   centerCommentedString(progname,ofil);
   ofil << "#" << endl;
   writeCommentedScrStarLine(ofil);
   return;
}
/* ******************************************************************************************* */
void writeScrStarLine(ofstream &ofil)
{
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   return;
}
/* ******************************************************************************************* */
void writeCommentedScrStarLine(ofstream &ofil)
{
   ofil << "#";
   for (int i=0; i<80; i++) {ofil << '*';}
   ofil << endl;
   return;
}
/* ******************************************************************************************* */
void writeScrCharLine(ofstream &ofil,char h)
{
   for (int i=0; i<80; i++) {ofil << h;}
   ofil << std::endl;
   return;
}
/* ******************************************************************************************* */
void writeCommentedScrCharLine(ofstream &ofil,char h)
{
   ofil << "#";
   for (int i=0; i<80; i++) {ofil << h;}
   ofil << std::endl;
   return;
}
/* ******************************************************************************************* */
void centerCommentedString(const string &s,ofstream &ofil)
{
   int pos=(int)((80-s.length())/2);
   ofil<<"#";
   for(int i=0;i<pos;i++){ofil<<" ";}
   ofil<<s<<endl;
}
/* ******************************************************************************************* */
void centerString(const string &s,ofstream &ofil)
{
   int pos=(int)((80-s.length())/2);
   for(int i=0;i<pos;i++){ofil<<" ";}
   ofil<<s<<endl;
}
/* ******************************************************************************************* */
void centerString(const char* word,ofstream &ofil)
{
   centerString(string(word),ofil);
   return;
}
/* ******************************************************************************************* */
void centerCommentedString(const char* word,ofstream &ofil)
{
   centerCommentedString(string(word),ofil);
   return;
}
/* ******************************************************************************************* */
void writeV3Comp(ofstream &ofil,const solreal (&v)[3])
{
   for (int i=0; i<3; i++) {ofil << v[i] << " ";}
   ofil << endl;
   return;
}
/* ******************************************************************************************* */
void writeV3Comp(ofstream &ofil,const std::string &s,const solreal (&v)[3])
{
   ofil << s;
   writeV3Comp(ofil,v);
   return;
}
/* ******************************************************************************************* */
void writeV3Comp(ofstream &ofil,const char* word,const solreal (&v)[3])
{
   ofil << word;
   writeV3Comp(ofil,v);
   return;
}
/* ******************************************************************************************* */
void writeM3x3Comp(ofstream &ofil,const solreal (&m)[3][3])
{
   for (int i=0; i<3; i++) {
      //for (int j=0; j<3; j++) {ofil << m[i][j] << " ";}
      ofil << m[i][0] << " " << m[i][1] << " " << m[i][2];
      ofil << endl;
   }
   return;
}
/* ******************************************************************************************* */
void writeM3x3Comp(ofstream &ofil,const std::string &s,const solreal (&m)[3][3])
{
   ofil << s;
   writeM3x3Comp(ofil,m);
   return;
}
/* ******************************************************************************************* */
void writeM3x3Comp(ofstream &ofil,const char* word,const solreal (&m)[3][3])
{
   ofil << word;
   writeM3x3Comp(ofil,m);
   return;
}
/* ******************************************************************************************* */
void genRandomTmpFileName(string &s, const int len) {
   static const int dm1=61;
   static const char alphanum[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   srand (clock());
   s="";
   for (int i=0; i<len; ++i) {s+=alphanum[rand()%(dm1)];}
   return;
}
/* ******************************************************************************************* */
int countColumnsInFile(string fnam)
{
   string line;
   ifstream ifil;
   ifil.open(fnam.c_str(),std::ios::in);
   discardComments(ifil);
   getline(ifil,line);
   ifil.close();
   removeSpacesLeft(line);
   line.append(" ");
   removeRedundantSpaces(line);
   int numtabs=0;
   for (size_t i=0; i<line.length(); i++) {if (line[i]==' ') {numtabs++;}}
   return numtabs;
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
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
/* ******************************************************************************************* */
# endif //_SOLFILEUTILS_CPP_
