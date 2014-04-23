/*

   crtflnms.h
   The actual implementation of the code is in the file crtflnms.cpp
   
   This file contains the interface definitions to modify or to create names for the 
   different files used to record the information extrated from the program. All the details of
   the functions contained in this file and the *.cpp version is below. (Here it is just provided
   a brief description of the function.)

   ------------------------

   Juan Manuel Solano Altamirano
   Adscription at the moment this project is initiated:
   Centro de Investigaciones y Estudios Avanzados del 
   Instituto Politécnico Nacional, 
   Unidad Monterrey, México.
   e-mail: jmsolanoalt@gmail.com
 
   Adscription at the moment the particular implementation for this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   
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

#ifndef _CRTFLNMS_H
#define _CRTFLNMS_H

#include "optflags.h"
#include <string>
using std::string;

/** This function makes the names out of a base name. 
   argv is the usual argv from the command line, and pos is the 
   position (or arg number) of the base name in the argv list.
   If no base name is given there, then the DEFAULTBASENAME will be
   used. The latter has been defined at the begining of the main program.
 */
void mkFileNames(char ** (&argv), optFlags &opts, string &i_fn, string &o_fn,
                 string &g_fn,string &l_fn);

#endif //_CRTFLNMS_H


