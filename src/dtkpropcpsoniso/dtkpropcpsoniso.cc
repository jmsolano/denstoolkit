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
#include <memory>
using std::shared_ptr;
#include <ctime>
#include "screenutils.h"
#include "fileutils.h"
#include "stringtools.h"
#include "mymath.h"
#include "gausswavefunction.h"
#include "fldtypesdef.h"
#include "bondnetwork.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helperspropcpsoniso.h"
#include "symmetricsurfacegrid.h"
#include "atomradiicust.h"
#include "matrixvectoroperations3d.h"

#ifndef PROPCPSONISOCPSIZE
#define PROPCPSONISOCPSIZE 0.10
#endif

/** Computes the angle between two 3D vectors. Notice that r is normalized!  */
double GetAngleBetweenVectors(vector<double> &xc,vector<double> &xx,vector<double> &xd);

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,logfilnam,povfilnam;
   string progname;
   OptionFlags options;
   ifstream ifile;
   ofstream ofile;

   getOptions(argc,argv,options); //This processes the options from the command line.
   int verboseLevel=0;
   if ( options.verboselevel ) { verboseLevel=std::stoi(string(argv[options.verboselevel])); }
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
   if ( options.stpspindens && gwf.ihaveSingleSpinOrbs ) {
      gwf.CalcCabAAndCabB();
   }

   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
   //was read, there souldn't be problems here.
   bnw.SetUpBNW();             //To setup the bond network.

   /* Setting the properties to be computed */

   char evprop='V';
   if (options.prop2eval) { evprop=argv[options.prop2eval][0]; }
   if ( evprop=='b' && (!gwf.ihaveCABSingleSpin) ) {
      ScreenUtils::DisplayErrorMessage("The alpha- and beta-spin density matrices could not\n"
            "be setup! Exiting...");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return EXIT_FAILURE;
   }
   char isoprop='d';
   if ( options.isoprop ) { isoprop=argv[options.isoprop][0]; }
   cout << "Logfile: " << logfilnam << '\n';

   if ( HelpersPropCPsOnIso::HaveIncompatibleOptions(argc,argv,options) ) { return EXIT_FAILURE; }

   /* Building the mesh  */
   double isovalue=0.001;
   if ( options.setisovalue ) {
      isovalue=std::stod(string(argv[options.setisovalue]));
   }
   if ( options.estimpkbaminesprim ) {
      isoprop='d';
      isovalue=0.0062;
      evprop='V';
      if ( options.setisovalue || options.isoprop || options.prop2eval ) {
         ScreenUtils::DisplayWarningMessage("Estimating pKb of a primary amine, setting\n"
               "isovalue rho=0.0062 a.u. and property evaluated at iso surf: MEP.");
      }
   }
   if ( options.estimpkbaminessec ) {
      isoprop='d';
      isovalue=0.0108e0;
      evprop='V';
      if ( options.setisovalue || options.isoprop || options.prop2eval ) {
         ScreenUtils::DisplayWarningMessage("Estimating pKb of a secondary amine, setting\n"
               "isovalue rho=0.0108 a.u. and property evaluated at iso surf: MEP.");
      }
   }
   if ( options.estimpkbaminester ) {
      isoprop='d';
      isovalue=0.0113e0;
      evprop='V';
      if ( options.setisovalue || options.isoprop || options.prop2eval ) {
         ScreenUtils::DisplayWarningMessage("Estimating pKb of a tertiary amine, setting\n"
               "isovalue rho=0.0113 a.u. and property evaluated at iso surf: MEP.");
      }
   }
   shared_ptr<MeshGrid> grid=HelpersPropCPsOnIso::BuildMesh(argc,\
         argv,options,gwf,bnw,isovalue);

   /* Computing critical points on the isosurface.  */

   cout << "Seeking critical points on" << (options.isofromcube? " " : " cap " ) << "isosurface...\n";
   if ( verboseLevel>0 ) {
      cout << "The isosurface has " << grid->vertex.size() << " vertices and "
           << grid->face.size() << " faces." << '\n';
   }
   vector<vector<double> > rcp;
   vector<double> vcp;
   vector<size_t> poscp;
   vector<int> sigcp;
   bool foundcps=HelpersPropCPsOnIso::SearchCPsIso(grid,gwf,rcp,poscp,sigcp,vcp,evprop);
   cout << "Done.\n";
   if ( foundcps && verboseLevel>0 ) {
      string cptype;
      ScreenUtils::PrintScrCharLine('+');
      for ( size_t i=0 ; i<rcp.size() ; ++i ) {
         if ( sigcp[i]==-3 ) { cptype="Maximum"; }
         else if ( sigcp[i]==3 ) { cptype="Minimum"; }
         else {cptype="Saddle critical point";}
         ScreenUtils::PrintScrCharLine('-');
         cout << scientific << setprecision(12);
         cout << cptype << " found at: "
              << rcp[i][0] << ' ' << rcp[i][1] << ' ' <<  rcp[i][2] << ' '
              << " : " << vcp[i] << " a.u., " << vcp[i]*627.5 << " kcal/mol, "
              << vcp[i]*2625.5 << " kJ/mol";
         if ( verboseLevel>1 ) {
            cout << " (" << grid->vneigh2v[poscp[i]].size() << " neighbs)" << '\n';
         }
         cout << '\n';
         ScreenUtils::PrintScrCharLine('-');
         if ( verboseLevel>2 ) {
            gwf.DisplayAllFieldProperties(rcp[i][0],rcp[i][1],rcp[i][2]);
         }
      }
      ScreenUtils::PrintScrCharLine('+');
   }
   string palettename="moreland";
   /* Rendering  */

   POVRayConfiguration pvp;
   if ( HelpersPropCPsOnIso::ComputeNormalsAtVertices(grid,gwf,isoprop) ) {
      grid->UseNormals(true);
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
      tcp[0]=rcp[i][0]; tcp[1]=rcp[i][1]; tcp[2]=rcp[i][2]; tcp[3]=PROPCPSONISOCPSIZE; // sphere coords and radius
      tcp[4]=(sigcp[i]==3 ? 0.5529412 : 0.00000e0); // red component (RGB)
      tcp[5]=0.00000e0; // green component (RGB)
      tcp[6]=(sigcp[i]==-3 ? 0.5529412 : 0.00000e0); // blue component (RGB)
      tcp[7]=vcp[i]; //field value at the sphere centre.
      sp.push_back(tcp); 
   }
   HelpersPropCPsOnIso::MakePovFile(povfilnam,options,pvp,bnw,grid,sp,palettename);

   /* At this point the computation has ended. Usually this means no errors ocurred. */
   if ( rcp.size()!=poscp.size() || rcp.size()!=sigcp.size() || rcp.size() !=vcp.size() ) {
      ScreenUtils::DisplayErrorMessage("Incorrect sizes!");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
   } else {
      ofstream ofil(logfilnam);
      int atnum;
      FileUtils::WriteHappyStart(argv,ofil,CURRENTVERSION,PROGRAMCONTRIBUTORS,true);
      ofil << scientific << setprecision(12);
      ofil << "#The isosurface was computed using ";
      if ( options.isofromcube ) {
         ofil << "the cube file '" << argv[options.isofromcube]
              << "',\n#  with a " << GetFieldTypeKeyLong(isoprop) << " isovalue of " << isovalue << '\n';
      } else {
         ofil << "a cap surrounding atom " << argv[options.setcentat] << '\n';
      }
      ofil << "#Below, (x,y,z) is the position of each critical point (CP) found upon the surface.\n"
           << "#sigma is the signature of the CP; +2 indicates a local minimum, and -2 indicates a local maximum.\n"
           << "#atomIdxInWF? is the atom index in the wfn (or wfx) file, i.e. the number of\n"
           << "#  the atom in the list of atoms, of the closest atom to the CP.\n"
           << "#atomSymbol is the atomic symbol of the closest atom to the CP\n"
           << "#dist2cent is the distance to the associated atom center (atomIdxInWF)\n"
           << "#angFromDir is the angle relative to the cap's direction." << '\n';
      ofil << "#        x                  y                  z        sigma       "
           << GetFieldTypeKeyShort(evprop) << "   atomIdxInWF? atomSymbol   dist2cent   angFromDir" << '\n';
      int cat;
      double dist2c,angrcpdir;
      vector<double> Rc(3),xd(3),xc(3);
      HelpersPropCPsOnIso::GetCenterIndexAndVectors(argv,options,bnw,cat,xc,xd);
      for ( size_t i=0 ; i<sp.size() ; ++i ) {
         atnum=HelpersPropCPsOnIso::FindClosestAtom(bnw,rcp[i]);
         for ( int j=0 ; j<3 ; ++j ) { Rc[j]=bnw.R[atnum][j]; }
         dist2c=0.0e0;
         for ( int j=0 ; j<3 ; ++j ) { dist2c+=((rcp[i][j]-Rc[j])*(rcp[i][j]-Rc[j])); }
         dist2c=sqrt(dist2c);
         angrcpdir=GetAngleBetweenVectors(xc,rcp[i],xd);
         ofil << rcp[i][0] << ' ' << rcp[i][1] << ' ' << rcp[i][2] << (sigcp[i]>0 ? " +" : " ")
              << sigcp[i] << ' ' << vcp[i] << ' '
              << (atnum+1) << ' ' << StringTools::RemoveAllDigits(gwf.atLbl[atnum])
              << ' ' << dist2c << ' ' << angrcpdir << '\n';
      }
      ofil.close();
   }

   if ( options.estimpkbaminesprim ) {
      HelpersPropCPsOnIso::EstimatepKbPrimaryAmine(vcp,sigcp);
      HelpersPropCPsOnIso::RequestCitation(argc,argv,options);
   }
   if ( options.estimpkbaminessec ) {
      HelpersPropCPsOnIso::EstimatepKbSecondaryAmine(vcp,sigcp);
      HelpersPropCPsOnIso::RequestCitation(argc,argv,options);
   }
   if ( options.estimpkbaminester ) {
      HelpersPropCPsOnIso::EstimatepKbTertiaryAmine(vcp,sigcp);
      HelpersPropCPsOnIso::RequestCitation(argc,argv,options);
   }

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
double GetAngleBetweenVectors(vector<double> &xc,vector<double> &xx,vector<double> &xd) {
   vector<double> r(3);
   for ( size_t i=0 ; i<3 ; ++i ) { r[i]=xx[i]-xc[i]; }
   MatrixVectorOperations3D::Normalize(r);
   double res=MatrixVectorOperations3D::InnerProduct(r,xd);
   return (acos(res)*180.0e0/M_PI);
}

