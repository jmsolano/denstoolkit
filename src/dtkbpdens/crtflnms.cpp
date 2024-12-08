/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/
/* crtflnms.cpp
   
   This file contains the implementation of the definitions to modify or to create names for the 
   different files used to record the information extracted from the program.

   ------------------------

   Juan Manuel Solano Altamirano
   Adscription at the moment this project is initiated:
   Centro de Investigaciones y Estudios Avanzados del 
   Instituto Politecnico Nacional, 
   Unidad Monterrey, Mexico.
   2011
   e-mail: jmsolanoalt@gmail.com
   
   Adscription at the moment the particular implementation for this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
*/
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include "crtflnms.h"
#include "../common/fileutils.h"
#include "../common/screenutils.h"
#include "../common/fldtypesdef.h"

void mkFileNames(char ** (&argv), OptionFlags &opts, string &i_fn, string &o_fn,
                 string &g_fn,string &l_fn) {
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
      ScreenUtils::SetScrRedBoldFont();
      cout << "\nError: the file " << i_fn << " is not a valid wave function file." << endl << endl;
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   o_fn=i_fn.substr(0,(i_fn.length()-3));
   g_fn=o_fn;
   l_fn=o_fn;
   g_fn.append("gnp");
   l_fn.append("log");
   o_fn.append("dat");
   
   pos=o_fn.find_last_of('.');
   if (pos!=string::npos) {
      string plbl;
      if (opts.uponbp) {plbl="BP";}
      if (opts.uponsl) {plbl="SL";}
      o_fn.insert(pos,plbl);
      g_fn.insert(pos,plbl);
      l_fn.insert(pos,plbl);
   }
   
   char prop;
   if (opts.prop2plot) {
      prop=argv[opts.prop2plot][0];
   } else {
      prop='d';
   }
   pos=o_fn.find_last_of('.');
   if (pos!=string::npos) {
      string plbl=GetFieldTypeKeyShort(prop);
      /*
      switch (prop) {
         case 'd':
            plbl="Rho";
            break;
         case 'g':
            plbl="MagGradRho";
            break;
         case 'l':
            plbl="LapRho";
            break;
         case 'E':
            plbl="ELF";
            break;
         case 'L':
            plbl="LOL";
            break;
         case 'M':
            plbl="MagGradLOL";
            break;
         case 'S':
            plbl="ShannEnt";
            break;
         case 'G':
            plbl="KinEnerDensG";
            break;
         case 'K':
            plbl="KinEnerDensK";
            break;
         default:
            plbl="";
            break;
      }
      // */
      o_fn.insert(pos,plbl);
      g_fn.insert(pos,plbl);
      l_fn.insert(pos,plbl);
   }
   
   if (opts.outfname) {
      o_fn=argv[opts.outfname];
      g_fn=o_fn;
      l_fn=g_fn;
      o_fn.append(".dat");
      g_fn.append(".gnp");
      l_fn.append(".log");
   }
   return;
}

