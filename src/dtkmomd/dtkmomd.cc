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
#include <sstream>
#include <cstdlib>
using std::exit;
#include <math.h>
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;
#include <ctime>

#include "../common/screenutils.h"
#include "../common/fileutils.h"
#include "../common/mymemory.h"
#include "../common/iofuncts-wfx.h"
#include "../common/iofuncts-wfn.h"
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
//#include "../common/bondnetwork.h"
#include "../common/solcubetools.h"
#include "optflags.h"
#include "crtflnms.h"

void makeLineGnuplotFile(optFlags &opts, string &gnpn,string &outn,char thefield);
void makePlaneGnuplotFile(optFlags &opts, string &gnpn,string &outn,solreal dimparam,\
      char thefield);
void makeLineDatFile(optFlags &opts,string &datnam,GaussWaveFunction &wf,int theaxis,\
      int npts,char thefield);
void makePlaneTsvFile(optFlags &opts,string &tsvnan,GaussWaveFunction &wf,int theplane,\
      int npts,char thefield);
void makeCubeFile(optFlags &opts,string &cubnam,GaussWaveFunction &wf,int npts,\
      char thefield,string &strfield);

int main (int argc, char ** argv)
{
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam,gnpnam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   int dim=0,axis=0,plane=0,npts=DEFAULTPOINTSPERDIRECTION;
   char field='d';
   solreal px,py,pz;
   px=py=pz=0.0e0;
   string extralbl="";
   if (options.evdim) {
      switch (argv[options.evdim][1]) {
         case '0':
            dim=0;
            sscanf (argv[options.evdim+1],"%lf",&px);
            sscanf (argv[options.evdim+2],"%lf",&py);
            sscanf (argv[options.evdim+3],"%lf",&pz);
            //cout << scientific << "px: " << px << ", py: " << py << ", pz: " << pz << endl;
            break;
         case '1':
            dim=1;
            switch (argv[options.evdim+1][0]) {
               case 'x':
                  axis=1;
                  break;
               case 'y':
                  axis=2;
                  break;
               case 'z':
                  axis=3;
                  break;
               default:
                  cout << endl;
                  ScreenUtils::DisplayErrorMessage("Invalid value for option -1...");
                  cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
                  cout << "\nto view the help menu.\n\n";
                  exit(1);
                  break;
            }
            extralbl=argv[options.evdim+1];
            extralbl.insert(0,"-P");
            break;
         case '2':
            dim=2;
            extralbl=string(argv[options.evdim+1]);
            if (!(extralbl==string("xy")||extralbl==string("xz")||extralbl==string("yz"))) {
               cout << endl;
               ScreenUtils::DisplayErrorMessage("Invalid value for option -2...");
               cout << "\nTry:\n\t" << argv[0] << " -h " << endl;
               cout << "\nto view the help menu.\n\n";
               exit(1);
            }
            if (extralbl==string("xy")) {plane=1;}
            if (extralbl==string("xz")) {plane=2;}
            if (extralbl==string("yz")) {plane=3;}
            extralbl.insert(1,"P");
            extralbl.insert(0,"-P");
            break;
         case '3':
            dim=3;
            npts=DEFAULTNPOINTSFORCUBE;
            break;
         default:
            break;
      }
   }
   if ( options.setfld ) {field=argv[options.setfld][0];}
   if (options.setn1) {sscanf (argv[options.setn1],"%d",&npts);}
   mkFileNames(argv,options,infilnam,outfilnam,gnpnam,dim,extralbl); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   GaussWaveFunction gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   string strfield;
   switch ( field ) {
      case 'd' :
         strfield="Momentum Density";
         break;
      case 'K' :
         strfield="Kinetic Energy Density (in Momentum Space)";
         break;
      default :
         break;
   }

   /* This evaluates the momentum density at a single point */
   
   if (dim==0) {
      cout << scientific << setprecision(12);
      cout << "\nMomentum density at\nP=(" << px << "," << py << "," << pz << "): " << endl;
      cout << endl << gwf.evalFTDensity(px,py,pz) << endl;
      cout << "\nK. E. in Momentum space at \nP=(" << px << "," << py << "," << pz << "): " << endl;
      cout << endl << gwf.evalFTKineticEnergy(px,py,pz) << endl;
   }
   
   /* This evaluates the momentum density on a line */
   
   if (dim==1) {
      cout << "The size of the grid will be " << npts << endl;
      cout << "The total number of points that will be computed is " << npts << endl;
      cout << "Evaluating and writing " << strfield << " on a line..." << endl;
      makeLineDatFile(options,outfilnam,gwf,axis,npts,field);
      if (options.mkplt) {makeLineGnuplotFile(options,gnpnam,outfilnam,field);}
   }
   
   /* This evaluates the momentum density on a plane */
   
   if (dim==2) {
      cout << "The size of the grid will be " << npts << " x " << npts << endl;
      cout << "The total number of points that will be computed is " << npts*npts << endl;
      cout << "Evaluating and writing " << strfield << " on a plane..." << endl;
      makePlaneTsvFile(options,outfilnam,gwf,plane,npts,field);
      solreal maxdim=DEFAULTMAXVALUEOFP;
      makePlaneGnuplotFile(options,gnpnam,outfilnam,maxdim,field);
   }
   
   /* This evaluates the momentum density on a cube */
   
   if (dim==3) {
      makeCubeFile(options,outfilnam,gwf,npts,field,strfield);
      if (options.mkplt) {
         ScreenUtils::DisplayWarningMessage("Plotting can only be performed with options -1 or -2.");
      }
   }
   
   cout << "Data saved in file: " << outfilnam << endl;
   
   /* At this point the computation has ended. Usually this means no errors ocurred. */
   
   ScreenUtils::PrintHappyEnding();
   ScreenUtils::SetScrGreenBoldFont();
   ScreenUtils::PrintScrStarLine();
   cout << setprecision(3) << "CPU Time: "
        << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   solreal end_walltime=time(NULL);
   cout << "Wall-clock time: " << solreal (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debuggin mode (under construction...)" << endl;
#endif
   ScreenUtils::PrintScrStarLine();
   ScreenUtils::SetScrNormalFont();
   return 0;
}

//**************************************************************************************************

void makeLineGnuplotFile(optFlags &opts, string &gnpn,string &outn,char thefield)
{
   ofstream gfil;
   gfil.open(gnpn.c_str());
   
   /* Choosing the label (legend) for the plot */
   
   string plbl;
   switch ( thefield ) {
      case 'd' :
         plbl="{/Symbol r}(p)";
         break;
      case 'K' :
         plbl="{/Bold K}(p)";
         break;
      default :
         break;
   }
   
   /* In this part the name is scanned for possible occurrings of the character '_'.
    For a proper display in the eps file, it has to be changed to "\_" */
   string line="";
   for (size_t i=0; i<outn.length(); i++) {
      if (outn[i]=='_') {line+='\\';}
      line+=outn[i];
   }
   gfil << "set title '" << line << "'" << endl;
   
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   
   string imgnam;
   imgnam=outn.substr(0,(outn.length()-3));
   imgnam.append("eps"); //here is set the name for the image
   string epsnam=imgnam;
   gfil << "#set output '";
   gfil << imgnam << "'" << endl; //so far, the image's name is *.eps
   imgnam=outn.substr(0,(outn.length()-3)); //Here the image name changes from eps to pdf
   imgnam.append("pdf");
   string pdfnam=imgnam;
#if (defined(__APPLE__)||defined(__linux__))
   gfil << "set output '| epstopdf --filter --outfile=";
   gfil << imgnam << "'" << endl;
#endif
#if (defined(__CYGWIN__))
   gfil << "set output '" << epsnam << "'" << endl;
#endif
   //gfil << "replot" << endl;
   gfil << "plot '" << outn << "' with lines ls 1 lw 2 title '" << plbl << "'" << endl;
   
   gfil.close();
#if (defined(__APPLE__)||defined(__linux__))
   line=string("gnuplot ");
#endif
#if (defined(__CYGWIN__))
   line=string("wgnuplot ");
#endif
   line+=gnpn;
   cout << "Calling gnuplot... ";
   system(line.c_str());
   cout << "Done." << endl;
#if(defined(__CYGWIN__))
#if _HAVE_EPSTOPDF_
   line="epstopdf --outfile="+pdfnam+" "+epsnam;
   cout << "Cleaning temporary files" << endl;
   line="rm -f "+epsnam;
   system(line.c_str());
#endif
#endif
   return;
}

//**************************************************************************************************

void makePlaneGnuplotFile(optFlags &opts, string &gnpn,string &outn,solreal dimparam,char thefield)
{
   ofstream gfil;
   gfil.open(gnpn.c_str());
   
   /* Choosing the label (legend) for the plot */
   
   string plbl="";
   switch ( thefield ) {
      case 'd':
         plbl="{/Symbol r}(p)";
         break;
      case 'K':
         plbl="{/Bold K}(p)";
         break;
      default :
         break;
   }
   
   gfil << "reset" << endl;
   
   /* In this part the name is scanned for possible occurrings of the character '_'.
    For a proper display in the eps file, it has to be changed to "\_" */
   string line="";
   for (size_t i=0; i<outn.length(); i++) {
      if (outn[i]=='_') {line+='\\';}
      line+=outn[i];
   }
   gfil << "set title '" << line << "'" << endl;
   
   /* Adding an underscore and braces to the atome labels in order to make the numbers
    subindices */
   
   gfil << "set isosample 300, 300" << endl;
   
   //gfil << "set zrange [" << minzrange << ":" << maxzrange << "]" << endl;
   
   gfil << "set contour base" << endl;
   
   gfil << "set cntrparam levels discrete 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0" <<endl;
   
   gfil << "unset surface" << endl;
   
   gfil << "set table 'contourtemp.dat'" << endl;
   
   gfil << "splot '" << outn << "'" << endl;
   
   gfil << "unset table" << endl;
   
   gfil << "reset" << endl;
   
   gfil << "unset key" << endl;
   
   gfil << "set size ratio 1" << endl;
   
   gfil << "set tmargin at screen 0.99\nset bmargin at screen 0.01\nset lmargin at screen 0.01" << endl;
   
   gfil << "dimparam=" << dimparam
   << " #Decrease this number to zoom in the plot" << endl;
   
   gfil << "set notics" << endl;
   
   gfil << "set xrange[-dimparam:dimparam]\nset yrange[-dimparam:dimparam]" << endl;
   
   //gfil << "minvalfield=" << minzrange << endl << "maxvalfield=" << maxzrange << endl;
   
   //gfil << "set cbrange [minvalfield:maxvalfield]" << endl;
   
   gfil << "set style data lines" << endl;
   
   gfil << "set palette rgbformulae 33,13,10" << endl;
   
   gfil << "plot '" << outn << "' with image" << endl;
   gfil << "replot 'contourtemp.dat' w l lt -1 lw 1 #Activate this line to show contours" << endl;
   
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   
   string imgnam,epsname,extepsname,pdfname;
   imgnam=outn.substr(0,(outn.length()-4));
   imgnam.append("ext.eps"); //here is set the name for the image
   gfil << "set output '";
   gfil << imgnam << "'" << endl; //so far, the image's name is *.eps
   gfil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   gfil << "replot" << endl;
   extepsname=imgnam;
   epsname=outn.substr(0,(outn.length()-3));
   epsname.append("eps");
   pdfname=outn.substr(0,(outn.length()-3));
   pdfname.append("eps");
   gfil << "#" << endl;
   FileUtils::WriteScrCharLine(gfil,'#');
   gfil << "#                 END OF GNUPLOT COMMANDS" << endl;
   FileUtils::WriteScrCharLine(gfil,'#');
   gfil << "#If you want to reconstruct the plot using this file, type:" << endl
   << "#gnuplot " << gnpn << endl
   << "#epstool --copy -b " << extepsname << " " << epsname << endl
   << "#epstopdf --outfile=" << pdfname << " " << epsname << endl;
   
   gfil.close();
#if (defined(__APPLE__)||defined(__linux__))
   line=string("gnuplot ");
#endif
#if (defined(__CYGWIN__))
   line=string("wgnuplot ");
#endif
   line+=gnpn;
   cout << "Calling gnuplot... ";
   system(line.c_str());
#if _HAVE_EPSTOOL_
   line=string("epstool --copy -b ");
   line+=(imgnam+string(" "));
   imgnam=outn.substr(0,(outn.length()-3)); //Here the image name changes from eps to pdf
   imgnam.append("eps");
   line+=imgnam;
   if (opts.quiet) {line+=string(" > /dev/null 2>&1");}
   system(line.c_str());
#endif
#if _HAVE_EPSTOPDF_
   string pdfnam=outn.substr(0,(outn.length()-3));
   pdfnam.append("pdf");
   line=string("epstopdf --outfile=")+pdfnam+string(" ")+imgnam;
   if (opts.quiet) {line+=string(" > /dev/null 2>&1");}
   system(line.c_str());
#endif
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   line="rm ";
#if _HAVE_EPSTOOL_
   line+=imgnam;
#endif
   imgnam=outn.substr(0,(outn.length()-4));
   imgnam.append("ext.eps");
   line+=(string(" ")+imgnam);
   system(line.c_str());
#endif//defined(__APPLE__)||defined(__linux__)
   cout << "Done." << endl;
   return;
}

/* ************************************************************************************ */
void makeLineDatFile(optFlags &opts,string &datnam,GaussWaveFunction &wf,int theaxis,\
      int npts,char thefield)
{
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   ofstream ofile;
   ofile.open(datnam.c_str(),ios::out);
   solreal dx,dy,dz,px,py,pz;
   dx=dy=dz=0.0e0;
   px=py=pz=0.0e0;
   switch (theaxis) {
      case 1:
         dx=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         px=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  ofile << px << " " << wf.evalFTDensity(px,py,pz) << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  ofile << px << " " << wf.evalFTKineticEnergy(px,py,pz) << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         
         break;
      case 2:
         dy=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         py=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  ofile << py << " " << wf.evalFTDensity(px,py,pz) << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            case 'K' :
              for (int i=0; i<npts; i++) {
                  ofile << py << " " << wf.evalFTKineticEnergy(px,py,pz) << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
      case 3:
         dz=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         pz=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  ofile << pz << " " << wf.evalFTDensity(px,py,pz) << endl;
                  pz+=dz;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            case 'K' :
              for (int i=0; i<npts; i++) {
                  ofile << pz << " " << wf.evalFTKineticEnergy(px,py,pz) << endl;
                  pz+=dz;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
              break;
            default :
               break;
         }
         
         break;
      default:
         break;
   }
   cout << endl;
   ofile.close();

}
/* ************************************************************************************ */
void makePlaneTsvFile(optFlags &opts,string &tsvnam,GaussWaveFunction &wf,int theplane,\
      int npts,char thefield)
{
   ofstream ofile;
   ofile.open(tsvnam.c_str(),ios::out);
   solreal dx,dy,dz,px,py,pz;
   dx=dy=dz=0.0e0;
   px=py=pz=0.0e0;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch (theplane) {
      case 1:
         dx=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         dy=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         px=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  py=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << py << "\t" << wf.evalFTDensity(px,py,pz) << endl;
                     py+=dy;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  py=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << py << "\t" << wf.evalFTKineticEnergy(px,py,pz) << endl;
                     py+=dy;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
      case 2:
         dx=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         dz=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         px=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << pz << "\t" << wf.evalFTDensity(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << pz << "\t" << wf.evalFTKineticEnergy(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
      case 3:
         dy=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         dz=2.0e0*DEFAULTMAXVALUEOFP/solreal(npts-1);
         py=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << py << "\t" << pz << "\t" << wf.evalFTDensity(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << py << "\t" << pz << "\t" << wf.evalFTKineticEnergy(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
         default:
            break;
      }
      cout << endl;
      ofile.close();
}
/* ************************************************************************************ */
void makeCubeFile(optFlags &opts,string &cubnam,GaussWaveFunction &wf,int npts,\
      char thefield,string &strfield)
{
   string comments="#Property: ";
   switch ( thefield ) {
      case 'd' :
         comments+="Momentum Density.";
         break;
      case 'K' :
         comments+="Kinetic Energy Density (in Momentum Space).";
         break;
      default :
         break;
   }
   solreal px,py,pz;
   px=py=pz=0.0e0;
   int boxnpts[3];
   for (int i=0; i<3; i++) {boxnpts[i]=npts;}
   solreal xin[3],delta[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         delta[i][j]=0.0e0;
      }
      delta[i][i]=2.0e0*DEFAULTMAXVALUEOFP/solreal(boxnpts[i]-1);
      xin[i]=-DEFAULTMAXVALUEOFP;
   }
   ofstream ofile;
   ofile.open(cubnam.c_str(),ios::out);
   writeCubeHeader(ofile,wf.title[0],comments,boxnpts,xin,delta,wf.nNuc,wf.atCharge,wf.R);
   solreal *prop1d;
   MyMemory::Alloc1DRealArray("prop1d",boxnpts[2],prop1d);
   cout << "The size of the grid will be " << boxnpts[0] << " x "
      << boxnpts[1] << " x " << boxnpts[2] << endl;
   cout << "The total number of points that will be computed is "
      << boxnpts[0]*boxnpts[1]*boxnpts[2] << endl;
   cout << "Evaluating and writing " << strfield << " on a cube..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch ( thefield ) {
      case 'd' :
         px=xin[0];
         for (int i=0; i<boxnpts[0]; i++) {
            py=xin[1];
            for (int j=0; j<boxnpts[1]; j++) {
               pz=xin[2];
               for (int k=0; k<boxnpts[2]; k++) {
                  prop1d[k]=wf.evalFTDensity(px,py,pz);
                  //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
                  pz+=delta[2][2];
               }
               writeCubeProp(ofile,boxnpts[2],prop1d);
               py+=delta[1][1];
            }
            px+=delta[0][0];
#if USEPROGRESSBAR
            ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((boxnpts[0]-1))));
#endif
         }
         cout << endl;
         break;
      case 'K' :
         px=xin[0];
         for (int i=0; i<boxnpts[0]; i++) {
            py=xin[1];
            for (int j=0; j<boxnpts[1]; j++) {
               pz=xin[2];
               for (int k=0; k<boxnpts[2]; k++) {
                  prop1d[k]=wf.evalFTKineticEnergy(px,py,pz);
                  //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
                  pz+=delta[2][2];
               }
               writeCubeProp(ofile,boxnpts[2],prop1d);
               py+=delta[1][1];
            }
            px+=delta[0][0];
#if USEPROGRESSBAR
            ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((boxnpts[0]-1))));
#endif
         }
         cout << endl;
         break;
      default :
         break;
   }
   //writeCubeProp(ofstream &ofil,int dim,solreal* (&prop));
   ofile.close();
   MyMemory::Dealloc1DRealArray(prop1d);
}
/* ************************************************************************************ */

