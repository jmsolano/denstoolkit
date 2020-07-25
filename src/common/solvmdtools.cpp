/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
#include "solvmdtools.h"
#include "fileutils.h"
#include "screenutils.h"

void VMDTools::writeSimpleVMDScript(string cubename,const char prop2rend,bool addquiet) {
   string vmdname=cubename;
   string tganame=cubename;
   FileUtils::ReplaceExtensionOfFileName(vmdname,"vmd");
   FileUtils::ReplaceExtensionOfFileName(tganame,"tga");
   double isovalue=getDefaultIsolvalueForCube(prop2rend);
   ofstream ofil(vmdname.c_str());
   if ( !ofil.good() ) {
      ScreenUtils::DisplayErrorMessage(string("Could not open the file '")+vmdname+string("'."));
   }
   ofil << "#VMD file created with DensToolKit." << endl;
   ofil << "#If you want to keep VMD's display on, comment or remove the last line 'quit'." << endl;
   ofil << "display resetview" << endl;
   ofil << "display projection Orthographic" << endl;
   ofil << "display nearclip set 0.000000" << endl;
   ofil << "display depthcue Off" << endl;
   ofil << "color Display Background white" << endl;
   ofil << "axes location Off" << endl;
   ofil << "mol addrep 0" << endl;
   ofil << "mol new {" << cubename
        << "} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }" << endl;
   ofil << "mol modstyle 0 0 Licorice 0.300000 12.000000 12.000000" << endl;
   ofil << "mol addrep 0" << endl;
   ofil << "mol modrep 1 0" << endl;
   ofil << "mol modstyle 1 0 Isosurface " << isovalue << " 0 0 0 1 1" << endl;
   ofil << "mol color Name" << endl;
   ofil << "mol material Opaque" << endl;
   ofil << "mol representation Isosurface " << isovalue << " 0 0 0 1 1" << endl;
   ofil << "render TachyonInternal " << tganame << endl;
   ofil << "#" << endl;
   if ( addquiet ) { ofil << "quit" << endl; }
   ofil << "#" << endl;
   ofil.close();
}

