/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
using std::ifstream;
#include "fileutils.h"
#include "screenutils.h"
#include "gaussiancube.h"

GaussianCube::GaussianCube() {
   x0.resize(3);
   dx.resize(3);
   for ( size_t i=0 ; i<3 ; ++i ) { dx[i].resize(3); }
   nx=ny=nz=0;
   useBohrUnits=true;
   cubeLoaded=false;
   data.clear();
}
GaussianCube::GaussianCube(const string &cubname) : GaussianCube() {
   if ( !FileUtils::ExtensionMatches(cubname,"cub") ) {
      ScreenUtils::DisplayWarningMessage(string("The file '")+cubname+
            string("' has no *cub or *cube extension."));
   }
   ifstream ifil(cubname.c_str());
   if ( !ifil.good() ) {
      ScreenUtils::DisplayErrorFileNotOpen(cubname);
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      ifil.close();
   } else {
      mol.ReadFromFile(ifil); //Recall this function sets buffer pos=0 (rewind file) after reading molecule's info.
      string tmp;
      for ( size_t i=0 ; i<2 ; ++i ) { std::getline(ifil,tmp); }
      int itmp;
      ifil >> itmp;
      for ( size_t i=0 ; i<3 ; ++i ) { ifil >> x0[i]; }
      ifil >> itmp;
      nx=(itmp >= 0 ? itmp : -itmp);
      if ( itmp<0 ) { useBohrUnits=false; }
      for ( size_t i=0 ; i<3 ; ++i ) { ifil >> dx[0][i]; }
      ifil >> itmp;
      ny=(itmp >= 0 ? itmp : -itmp);
      if ( itmp<0 ) { useBohrUnits=false; }
      for ( size_t i=0 ; i<3 ; ++i ) { ifil >> dx[1][i]; }
      ifil >> itmp;;
      nz=(itmp >= 0 ? itmp : -itmp);
      if ( itmp<0 ) { useBohrUnits=false; }
      for ( size_t i=0 ; i<3 ; ++i ) { ifil >> dx[2][i]; }
      nyz=ny*nz;
      if ( !useBohrUnits ) {
         ScreenUtils::DisplayWarningMessage("negative cube, see Paul Bourke webpage about cube format.");
      }
      for ( size_t i=0 ; i<=mol.Size() ; ++i ) { //required to read an empty line...
         std::getline(ifil,tmp);
      }
      data.resize(nx*ny*nz);
      for ( size_t i=0 ; i<nx ; ++i ) {
         for ( size_t j=0 ; j<ny ; ++j ) {
            for ( size_t k=0 ; k<nz ; ++k ) { ifil >> data[i*nyz+j*nz+k]; }
         }
      }
      ifil.close();
      cubeLoaded=true;
   }
}
GaussianCube::~GaussianCube() {
   data.clear();
   x0.clear();
}

