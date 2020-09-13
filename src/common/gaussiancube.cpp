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
   data.clear();
}
GaussianCube::GaussianCube(const string &cubname) : GaussianCube() {
   if ( !FileUtils::ExtensionMatches(cubname,"cub") ) {
      ScreenUtils::DisplayWarningMessage(string("The file '")+cubname+
            string("' has no *cub or *cube extension."));
   }
   ifstream ifil(cubname.c_str());
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
}
GaussianCube::~GaussianCube() {
   data.clear();
   x0.clear();
}

