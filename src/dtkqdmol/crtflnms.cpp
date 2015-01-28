/*

   crtflnms.cpp
   
   This file contains the implementation of the definitions to modify or to create names for the 
   different files used to record the information extracted from the program.

   ------------------------

   Juan Manuel Solano Altamirano
   Adscription at the moment this project is initiated:
   Centro de Investigaciones y Estudios Avanzados del 
   Instituto Politécnico Nacional, 
   Unidad Monterrey, Mexico.
   2011
   e-mail: jmsolanoalt@gmail.com
   
   Adscription at the moment the particular implementation for this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
*/
#ifndef _CRTFLNMS_CPP
#define _CRTFLNMS_CPP
#include "crtflnms.h"
#include "../common/solfileutils.h"
#include "../common/solscrutils.h"

void mkFileNames(char ** (&argv), optFlags &opts, string &i_fn,string &g_fn)
{
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   /*
      If you need more names to be created by this function, you need to add the new
      string in the arguments list here and in the corresponding header file.
    */
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   i_fn=string(argv[1]);
   size_t pos;
   //string sl="wfn";
   //pos=i_fn.find(sl);
   if (!((i_fn.find("wfn")!=string::npos)||
       (i_fn.find("WFN")!=string::npos)||
       (i_fn.find("wfx")!=string::npos)||
       (i_fn.find("WFX")!=string::npos))) {
      setScrRedBoldFont();
      cout << "\nError: the file " << i_fn << " is not a valid wave function file." << endl << endl;
      setScrNormalFont();
      exit(1);
   }
   g_fn=i_fn.substr(0,(i_fn.length()-3));
   g_fn.append("gnp");
   pos=g_fn.find_last_of('.');
   if (pos!=string::npos) {
      string plbl="QDMol";
      g_fn.insert(pos,plbl);
   }
   if (opts.outfname) {
      g_fn=argv[opts.outfname];
      g_fn.append(".gnp");
   }
   return;
}

#endif //_CRTFLNMS_CPP
