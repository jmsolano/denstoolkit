/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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

/*

   optsflags.h
   The actual implementation of the code is in the file optsflags.cpp
   
   This file contains the interface definitions of the class optflags, which is
   simply a class that will contain all the flags used during the execution of the
   main program.
   
   It is intended to simplify the addition of new options and new variations of calculations
   the program makes. 

   ------------------------

   Juan Manuel Solano Altamirano
   Adscription at the moment the template for this class is initiated:
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

#ifndef _OPTSFLAGS_H
#define _OPTSFLAGS_H
#include "../common/fldtypesdef.h"
#include <string>
using std::string;

class OptionFlags
{
public: 
   OptionFlags();//default constructor, initialize all the flags to convenient (default) values.
   unsigned short int infname,outfname;
   unsigned short int zipdat,quiet;
   unsigned short int atpos;
   bool stpspindens;
};//end class optsFlags
void printErrorMsg(char** &argv,char lab);
void printHelpMenu(int &argc, char** &argv);//self-described
void getOptions(int &argc, char** &argv, OptionFlags &flags);//this function will assign the values to
                                                           //all the flags. Implementation is in optsflags.cpp
void processDoubleDashOptions(int &argc,char** &argv,OptionFlags &flags,int pos);
#endif //_OPTSFLAGS_H


