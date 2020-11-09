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
/*
   Adscription at the moment this program is started:
   University of Guelph,
   Guelph, Ontario, Canada.
   May 2013
*/

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
using std::scientific;
#include <ctime>
#include "screenutils.h"
#include "fileutils.h"
#include "mymath.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helperspropcpsoniso.h"
#include "symmetricsurfacegrid.h"
#include "atomradiicust.h"

#ifndef PROPCPSONISOCPSIZE
#define PROPCPSONISOCPSIZE 0.10
#endif

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,logfilnam,povfilnam;
   string progname;
   OptionFlags options;
   ifstream ifile;
   ofstream ofile;

   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,logfilnam,povfilnam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   //Just to let the user know that the initial configuration is OK

   cout << endl << "Loading wave function from file: " << infilnam << "... ";

   GaussWaveFunction gwf;
   if (!(gwf.ReadFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(EXIT_FAILURE);
   }
   cout << "Done." << endl;

   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
   //was read, there souldn't be problems here.
   bnw.SetUpBNW();             //To setup the bond network.

   /* Setting the properties to be computed */

   char evprop='V';
   if (options.prop2eval) { evprop=argv[options.prop2eval][0]; }
   char isoprop='d';
   if ( options.isoprop ) { isoprop=argv[options.isoprop][0]; }
   cout << "Logfile: " << logfilnam << '\n';

   if ( HelpersPropCPsOnIso::HaveIncompatibleOptions(argc,argv,options) ) { return EXIT_FAILURE; }

   int cAt;
   vector<double> xc(3),xd(3);
   HelpersPropCPsOnIso::GetCenterIndexAndVectors(argv,options,bnw,cAt,xc,xd);
   cout << "xc: " << xc[0] << ' ' << xc[1] << ' ' << xc[2] << '\n';
   cout << "xd: " << xd[0] << ' ' << xd[1] << ' ' << xd[2] << '\n';
   SymmetricSurfaceGrid grid;
   int refCapMeshLevel=4;
   if ( options.refinemesh ) {
      refCapMeshLevel=std::stoi(string(argv[options.refinemesh]));
      if ( refCapMeshLevel<0 ) { refCapMeshLevel=0; }
      if ( refCapMeshLevel>8 ) {
         ScreenUtils::DisplayWarningMessage("This level of refinement may cause numerical issues.");
      }
   }
   grid. SetupSphereIcosahedron(refCapMeshLevel);
   grid.Translate(xc);
   grid.Scale(GetAtomicVDWRadius(bnw.atNum[cAt]));
   /*
      cout << "Meshgrid contains " << grid.vertex.size() << " vertices and "
      << grid.tface.size() << " before trimming." << '\n'; // */
   grid.TrimFacesCentroidDotProdGreaterThanZero(xd);
   /*
      cout << "Meshgrid contains " << grid.vertex.size() << " vertices and "
      << grid.tface.size() << " after trimming." << '\n'; // */

   /*
      for ( size_t i=0 ; i<grid.vertex.size() ; ++i ) {
      cout << grid.vertex[i][0] << ' ' << grid.vertex[i][1] << ' ' << grid.vertex[i][2] << ' '
      << gwf.EvalMolElecPot(grid.vertex[i][0],grid.vertex[i][1],grid.vertex[i][2]) << '\n';
      }
   // */

   /* Looking for partial isosurface  */
   //grid.DisplayFaces();

   /*
      cout << "Meshgrid contains " << grid.vertex.size() << " vertices and "
      << grid.tface.size() <<  " faces before projecting." << '\n'; // */
   cout << "Computing cap isosurface...\n";
   HelpersPropCPsOnIso::ProjectGridOntoIsosurface(gwf,grid,isoprop,0.001);
   cout << "Done\n";
   /*cout << "Meshgrid contains " << grid.vertex.size() << " vertices and "
     << grid.tface.size() <<  " faces after projecting." << '\n'; // */

   /*
      for ( size_t i=0 ; i<grid.vertex.size() ; ++i ) {
      cout << grid.vertex[i][0] << ' ' << grid.vertex[i][1] << ' ' << grid.vertex[i][2] << ' '
      << gwf.EvalMolElecPot(grid.vertex[i][0],grid.vertex[i][1],grid.vertex[i][2]) << '\n';
      }
   // */
   cout << "Seeking critical points on cap isosurface...\n";
   vector<vector<double> > rcp;
   vector<double> vcp;
   vector<size_t> poscp;
   vector<int> sigcp;
   bool foundcps=HelpersPropCPsOnIso::SearchCPs(grid,gwf,\
         rcp,poscp,sigcp,vcp,'V');
   cout << "Done.\n";
   if ( foundcps ) {
      string cptype;
      ScreenUtils::PrintScrCharLine('+');
      for ( size_t i=0 ; i<rcp.size() ; ++i ) {
         if ( sigcp[i]==-3 ) { cptype="Maximum"; }
         else if ( sigcp[i]==3 ) { cptype="Minimum"; }
         else {cptype="Saddle critical point";}
         ScreenUtils::PrintScrCharLine('-');
         cout << cptype << " found at: "
              << rcp[i][0] << ' ' << rcp[i][1] << ' ' <<  rcp[i][2] << ' '
              << " : " << vcp[i] << '\n';
         ScreenUtils::PrintScrCharLine('-');
         gwf.DisplayAllFieldProperties(rcp[i][0],rcp[i][1],rcp[i][2]);
      }
      ScreenUtils::PrintScrCharLine('+');
   }
   /* Rendering  */

   POVRayConfiguration pvp;
   if ( HelpersPropCPsOnIso::ComputeNormalsAtVertices(grid,gwf,isoprop) ) {
      grid.UseNormals(true);
   }
   double cv[3]={0.0e0,0.0e0,0.0e0};
   if ( options.setgnpangles ) {
      cv[0]=std::stod(string(argv[options.setgnpangles]));
      cv[2]=std::stod(string(argv[options.setgnpangles+1]));
   }
   if ( options.setviewangles ) {
      cv[0]=std::stod(string(argv[options.setviewangles  ]));
      cv[1]=std::stod(string(argv[options.setviewangles+1]));
      cv[2]=std::stod(string(argv[options.setviewangles+2]));
   }
   pvp.SetAngView(cv[0],cv[1],cv[2]);
   cout << "Generating povray file (" << povfilnam << ")" << '\n';
   vector<vector<double> > sp(0);
   vector<double> tcp(8);
   for ( size_t i=0 ; i<rcp.size() ; ++i ) {
      tcp[0]=rcp[i][0]; tcp[1]=rcp[i][1]; tcp[2]=rcp[i][2]; tcp[3]=PROPCPSONISOCPSIZE;
      tcp[4]=(sigcp[i]==3 ? 0.5529412 : 0.00000e0);
      tcp[5]=0.00000e0;
      tcp[6]=(sigcp[i]==-3 ? 0.5529412 : 0.00000e0);
      tcp[7]=vcp[i];
      sp.push_back(tcp);
   }
   HelpersPropCPsOnIso::MakePovFile(povfilnam,options,pvp,bnw,grid,sp);

   /* At this point the computation has ended. Usually this means no errors ocurred. */

   ScreenUtils::PrintHappyEnding();
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintScrStarLine();
   cout << setprecision(3) << "CPU Time: "
      << double( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   double end_walltime=time(NULL);
   cout << "Wall-clock time: " << double (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debuggin mode (under construction...)" << endl;
#endif
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrNormalFont();
   return EXIT_SUCCESS;
}

