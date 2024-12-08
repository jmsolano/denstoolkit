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
/*
   Adscription at the moment this program is started:
   Meritorious Autonomous Unversity of Puebla
   Puebla, Pue., Mexico
   2020
 */
#include <iostream>
using std::cout;
#include <cstdlib>
using std::exit;
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
#include <memory>
using std::shared_ptr;
#include <ctime>
#include "screenutils.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "wfgrid3d.h"
#include "gaussiancube.h"
#include "isosurface.h"
#include "palette.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helpersmapfieldonisosurf.h"
int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const double begin_walltime = time(NULL);
   string infilnam,cubfnam,povfnam;
   string progname;

   OptionFlags options;
   getOptions(argc,argv,options); //This processes the options from the command line.
   getWFName(argv,options,infilnam);
   //Just to let the user know that the initial configuration is OK:
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS);

   cout << '\n' << "Loading wave function from file: " << infilnam << '\n';
   
   GaussWaveFunction gwf;
   if (!(gwf.ReadFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << '\n';
   if ( options.stpspindens && gwf.ihaveSingleSpinOrbs ) {
      gwf.CalcCabAAndCabB();
   }

   BondNetWork bnw;
   bnw.ReadFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.SetUpBNW(); //To setup the bond network.
   
   WaveFunctionGrid3D grid;
   
   /* Looking for user grid dimensions */
   
   //*
   int nn=DEFAULTPOINTSPERDIRECTION;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&nn); //same number of points per direction
      grid.SetNPts(nn);
   } else if (options.setn3) {
      int ny,nz;
      sscanf(argv[options.setn3],"%d",&nn);  //different number of points per direction
      sscanf(argv[options.setn3+1],"%d",&ny); //manually given by the user
      sscanf(argv[options.setn3+2],"%d",&nz);
      grid.SetNPts(nn,ny,nz);
   } else if (options.setsmcub1) {
      sscanf(argv[options.setsmcub1],"%d",&nn); //uses nn for determining the number of points
                     //of the largest axis. The rest will have proportional number to their length
   } else {
      grid.SetNPts(nn);
   }
   // */
   
   if ((options.setsmcub)||(options.setsmcub1)) {
      grid.SetUpSmartCuboidGrid(gwf,bnw,nn);
   } else {
      grid.SetUpSimpleGrid(gwf,bnw);
   }

   /* Special configurations  */
   if ( options.plotnci && options.plotdori ) {
      ScreenUtils::DisplayErrorMessage("--dori and --nci options are mutually exclusive,\n"
            "please, choose only one.");
      return EXIT_FAILURE;
   }

   /* Set the isosurface property to be computed */
   char isoprop='z';
   /* Set the property to be mapped */
   char mapprop='Z';

   if (options.plotnci) {
      if ( options.configspecialnci ) {
         double ttt=std::stod(string(argv[options.configspecialnci]));
         gwf.SetNCIRhoMin(ttt);
         ttt=std::stod(string(argv[options.configspecialnci+1]));
         gwf.SetNCIRhoMax(ttt);
         ttt=std::stod(string(argv[options.configspecialnci+2]));
         gwf.SetNCISMax(ttt);
      }
      cout << "nciRhoMin: " << gwf.nciRhoMin << '\n';
      cout << "nciRhoMax: " << gwf.nciRhoMax << '\n';
      cout << "nciSMax: " << gwf.nciSMax << '\n';
   } else if ( options.plotdori ) {
      isoprop='D';
      gwf.SetNCIRhoMin(-0.02);
      gwf.SetNCIRhoMax(0.02);
      cout << "Lambda2Min: " << gwf.nciRhoMin << '\n';
      cout << "Lambda2Max: " << gwf.nciRhoMax << '\n';
   } else {
      if ( options.setisoprop ) {
         isoprop=argv[options.setisoprop][0];
      }
      if ( options.setmapprop ) {
         mapprop=argv[options.setmapprop][0];
      }
   }
   mkFileNames(argv,options,infilnam,cubfnam,povfnam,isoprop,mapprop); //This creates the names used.
                                                       //
   /* Computes cube (if requested). */
   if ( !options.skipcube ) {
      if ( isoprop=='b' && (!gwf.ihaveCABSingleSpin) ) {
         ScreenUtils::DisplayErrorMessage("The alpha- and beta-spin density matrices could not\n"
               "be setup! Exiting...");
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         return EXIT_FAILURE;
      }
      cout << "The size of the grid will be: " << grid.GetNPts(0) << " x "
           << grid.GetNPts(1) << " x " << grid.GetNPts(2) << '\n';
      cout << "Total number of points that will be computed: "
           << (grid.GetNPts(0)*grid.GetNPts(1)*grid.GetNPts(2)) << '\n';
      cout << "Evaluating and writing " << GetFieldTypeKeyLong(isoprop) << ")." << '\n' << '\n';
      ScalarFieldType ft=Char2ScalarFieldType(isoprop);
      grid.MakeCube(cubfnam,gwf,ft);
   }
   cout << '\n' << GetFieldTypeKeyLong(isoprop) << " cube file: " << cubfnam << '\n';
   
   cout << "Extracting isosurface..." << std::endl;
   GaussianCube gc(cubfnam);
   if ( !gc.CubeLoaded() ) {
      ScreenUtils::DisplayErrorMessage("Please deactivate the option -c");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return EXIT_FAILURE;
   }
   Isosurface iso;
   //iso.UseTetrahedrons(true);
   double isovalue=GetDefaultIsolvalueForCube(mapprop);
   if ( options.setsisovalue ) {
      isovalue=std::stod(string(argv[options.setsisovalue]));
   }
   iso.ExtractMarchingCubes(gc,isovalue);
   if ( mapprop=='0' ) {
      iso.UseColorMap(false);
      ScreenUtils::SetScrGreenBoldFont();
      cout << "Drawing the isosurface without color map." << '\n';
      ScreenUtils::SetScrNormalFont();
      double rr=0.7,gg=0.7,bb=0.7;
      if ( options.setgrayrgb ) {
         rr=gg=bb=std::stod(string(argv[options.setgrayrgb]));
      }
      if ( options.setfullrgb ) {
         rr=std::stod(string(argv[options.setfullrgb+0]));
         gg=std::stod(string(argv[options.setfullrgb+1]));
         bb=std::stod(string(argv[options.setfullrgb+2]));
      }
      iso.SetRGB(rr,gg,bb);
   } else {
      if ( mapprop=='b' && (!gwf.ihaveCABSingleSpin) ) {
         ScreenUtils::DisplayErrorMessage("The alpha- and beta-spin density matrices could not\n"
               "be setup! Exiting...");
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         return EXIT_FAILURE;
      }
      cout << "Computing " << GetFieldTypeKeyLong(mapprop) << " at mesh vertices..." << '\n';
      HelpersMapFieldOnIsoSurf::ComputeFieldAtVertices(gwf,iso,mapprop);
   }
   cout << "Computing Normals at mesh vertices..." << '\n';
   HelpersMapFieldOnIsoSurf::ComputeNormalsAtVertices(gwf,iso,isoprop);
   iso.UseNormals(true);

   /* At this point the computation has ended. Usually this means no errors ocurred. */

   /* Rendering  */
   double Lmin=-0.02e0;
   double Lmax=0.02e0;
   if ( mapprop!='0' ) {
      cout << GetFieldTypeKeyLong(mapprop) << " limits over iso of '" << cubfnam << "':\n";
      cout << "    " << GetFieldTypeKeyLong(mapprop) << " min: " << iso.MinP2Map() << '\n';
      cout << "    " << GetFieldTypeKeyLong(mapprop) << " max: " << iso.MaxP2Map() << '\n';
      if ( options.setcolorscalesingle ) {
         Lmax=std::stod(string(argv[options.setcolorscalesingle]));
         Lmin=-Lmax;
      }
      if ( options.setcolorscaleboth ) {
         Lmin=std::stod(string(argv[options.setcolorscaleboth]));
         Lmax=std::stod(string(argv[options.setcolorscaleboth+1]));
      }
   }
   POVRayConfiguration pvp;
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
   cout << "Generating povray file (" << povfnam << ")" << '\n';
   cout << "Color scale: [" << Lmin << ',' << Lmax << ']' << '\n';
   iso.SetMinP2Map(Lmin);
   iso.SetMaxP2Map(Lmax);
   string palname="greens";
   if ( options.selectpalette ) {
      palname=argv[options.selectpalette];
   }
   HelpersMapFieldOnIsoSurf::MakePovFile(povfnam,pvp,bnw,iso,options,argv);

   /* In principle, everything went OK. Exiting.  */
   ScreenUtils::PrintHappyEnding();
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintScrStarLine();
   cout << setprecision(3) << "CPU Time: " << double( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << '\n';
   double end_walltime=time(NULL);
   cout << "Wall-clock time: " << double (end_walltime-begin_walltime) << "s" << '\n';
#if DEBUG
   cout << "Debugging mode (under construction)..." << '\n';
#endif
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrNormalFont();
   return EXIT_SUCCESS;
}

