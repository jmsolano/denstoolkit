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
   Instituto Politécnico Nacional, 
   Unidad Monterrey, Mexico.
   2011
   e-mail: jmsolanoalt@gmail.com
   
   Adscription at the moment the particular implementation for this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
   ------------------------

	 This code is free code; you can redistribute it and/or
	 modify it under the terms of the GNU General Public License
	 as published by the Free Software Foundation; either version 2
	 of the License, or (at your option) any later version.

	 This program is distributed in the hope that it will be useful,
	 but WITHOUT ANY WARRANTY; without even the implied warranty of
	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
	 along with this program; if not, write to the Free Software 
	 Foundation, Inc., 59 Temple Place - Suite 330, 
	 Boston, MA  02111-1307, USA.

   WWW:  http://www.gnu.org/copyleft/gpl.html
	
	----------------------
*/

#ifndef _OPTSFLAGS_H
#define _OPTSFLAGS_H

#include <string>
using std::string;

class optFlags
{
public: 
   optFlags();//default constructor, initialize all the flags to convenient (default) values.
   unsigned short int infname,outfname,setn1,setats,setstep;
   unsigned short int uponbp,uponsl,prop2plot;
   unsigned short int zipdat,mkplt,kpgnp,quiet,showcont,showatlbls;
   unsigned short int setinccont,findcps;
};//end class optsFlags
void printErrorMsg(char** &argv,char lab);
void printHelpMenu(int &argc, char** &argv);//self-described
void getOptions(int &argc, char** &argv, optFlags &flags);//this function will assign the values to
                                                           //all the flags. Implementation is in optsflags.cpp
void processDoubleDashOptions(int &argc,char** &argv,optFlags &flags,int pos);
#endif //_OPTSFLAGS_H


