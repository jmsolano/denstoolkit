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
 October 2013
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

#include "../common/solscrutils.h"
#include "../common/solfileutils.h"
#include "../common/solmemhand.h"
#include "../common/solmath.h"
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "../common/solstringtools.h"
#include "../common/fldtypesdef.h"
#include "solgnuplottools.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helpersplot.h"

int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam,gnpnam,lognam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam,gnpnam,lognam); //This creates the names used.
   printHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   /* Stablishes the field property to be evaluated. */
   
   char prop='d';
   if (options.prop2plot) { prop=argv[options.prop2plot][0]; }

   /* Checks for non valid property fields */
   
   if (!(prop=='d'||prop=='g'||prop=='l'||prop=='E'||prop=='L'||prop=='M'\
         ||prop=='S'||prop=='G'||prop=='K'||prop=='V'\
         ||prop=='P'||prop=='r'||prop=='s'||prop=='u'\
         ||prop=='e')) {
      displayErrorMessage("Non valid field type");
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
      exit(1);
   }
   
   /* Loads the wave function */
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   GaussWaveFunction gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      setScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      setScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   
   if (gwf.nNuc==1) {
      displayWarningMessage("This file contains only one atom... There are no bond paths...");
      cout << "Nothing to do!" << endl;
      gwf.~GaussWaveFunction();
      exit(0);
   }
   
   /* Sets the atoms for defining the line/bond path */
   
   int at1=0,at2=1;
   if (options.setats) {
      sscanf(argv[options.setats],"%d",&at1);
      at1--;
      sscanf(argv[options.setats+1],"%d",&at2);
      at2--;
   }
   if (at1>=gwf.nNuc||at2>=gwf.nNuc||at1<0||at2<0) {
      displayErrorMessage("Requesting a non existent atom!");
      gwf.~GaussWaveFunction();
      exit(1);
   }
   
   /* Sets the bond network of the molecule */
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW();             //To setup the bond network.
   
   /* Defines the main critical point network object. */
   
   critPtNetWork cpn(gwf,bnw);
   
   int dimarr=300;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&dimarr);
      if (dimarr<=0) {
         displayErrorMessage("Please provide a positive number for the number of points!");
         gwf.~GaussWaveFunction();
         bnw.~bondNetWork();
         cpn.~critPtNetWork();
         exit(1);
      }
   }
   
   /* Sets the final file names. */
   
   string lbls="-"+gwf.atLbl[at1]+"-"+gwf.atLbl[at2];
   size_t pos=gnpnam.find_last_of('.');
   if (pos!=string::npos) {
      gnpnam.insert(pos,lbls);
      outfilnam.insert(pos,lbls);
      lognam.insert(pos,lbls);
   }
   
   /* Declares an array to store the coordinates of the bond path points. */
   solreal **rbgp=NULL;
   alloc2DRealArray(string("rbgp"),dimarr,3,rbgp);
   
   int nbgppts=1;
   solreal dl=DEFAULTBONDPATHSTEPMD1;
   solreal robcp[3];
   solreal x1[3],x2[3],dx[3];
   solreal lenline=0.0e0;
   solreal dist=0.0e0;
   
   if (options.setstep) {sscanf(argv[options.setstep],"%lg",&dl);}
   
   /* ==================================================================== */
   
   if (options.uponbp) {
      /* Computes the gradient path */
      nbgppts=cpn.findSingleRhoBondGradientPathRK5(at1,at2,dl,dimarr,rbgp,robcp);
      //cout << "npts: " << nbgppts << endl;
      for (int i=0; i<3; i++) {x1[i]=rbgp[0][i];}
      for (int i=1; i<nbgppts; i++) {
         dist=0.0e0;
         for (int j=0; j<3; j++) {
            x2[j]=rbgp[i][j];
            dist+=((x2[j]-x1[j])*(x2[j]-x1[j]));
            x1[j]=x2[j];
         }
         lenline+=sqrt(dist);
      }
      cout << "The bond path consists of " << nbgppts << " points" << endl;
      cout << "The value of the step: " << dl << endl;
   }
   
   if (options.uponsl) {
      nbgppts=dimarr;
      lenline=0.0e0;
      for (int i=0; i<3; i++) {
         lenline+=((bnw.R[at1][i]-bnw.R[at2][i])*(bnw.R[at1][i]-bnw.R[at2][i]));
      }
      lenline=sqrt(lenline);
      dl=lenline/solreal(dimarr-1);
      for (int i=0; i<3; i++) {
         dx[i]=(bnw.R[at2][i]-bnw.R[at1][i])/solreal(dimarr-1);
         rbgp[0][i]=bnw.R[at1][i];
         robcp[i]=0.0e0;
      }
      for (int i=1; i<dimarr; i++) {for (int j=0; j<3; j++) {rbgp[i][j]=rbgp[i-1][j]+dx[j];}}
      cout << "The line that joins the atoms consists of " << nbgppts << " points." << endl;
   }
   
   /* Opens the dat file */
   
   ofile.open(outfilnam.c_str(),ios::out);
   ofile << scientific << setprecision(12);
   cout << scientific << setprecision(12);
   
   cout << "Evaluating " << getFieldTypeKeyLong(prop) << " upon ";
   if (options.uponbp) {cout << "the bond path:" << endl;}
   if (options.uponsl) {cout << "the bond line:" << endl;}

#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   
   solreal xx[3],xt[3],tmpval,pp,pbcp,pmax,pmin,maxval=-1.0e+50,minval=1.0e+50,bcpval;
   solreal xmax[3],xmin[3];
   pp=0.0e0;
   if (options.uponbp) {
      if ((rbgp[0][0]==robcp[0])&&(rbgp[0][1]==robcp[1])&&(rbgp[0][2]==robcp[2])) {pbcp=pp;}
   }
   for (int k=0; k<3; k++) {xx[k]=rbgp[0][k];}
   tmpval=HelpersPlot::EvalFieldProperty(prop,xx,gwf);
   if (tmpval<minval) {minval=tmpval; pmin=pp; for(int k=0; k<3; k++) {xmin[k]=xx[k];}}
   if (tmpval>maxval) {maxval=tmpval; pmax=pp; for(int k=0; k<3; k++) {xmax[k]=xx[k];}}
   ofile << pp << " ";
   for (int k=0; k<3; k++) {ofile << xx[k] << " ";}
   ofile << tmpval << endl;
   for (int i=1; i<nbgppts; i++) {
      dist=0.0e0;
      for (int k=0; k<3; k++) {
         xt[k]=xx[k]; xx[k]=rbgp[i][k];
         xt[k]-=xx[k]; dist+=(xt[k]*xt[k]);
      }
      pp+=sqrt(dist);
      tmpval=HelpersPlot::EvalFieldProperty(prop,xx,gwf);
      ofile << pp << " ";
      for (int k=0; k<3; k++) {ofile << xx[k] << " ";}
      ofile << tmpval << endl;
      if (options.uponbp) {
         if ((xx[0]==robcp[0])&&(xx[1]==robcp[1])&&(xx[2]==robcp[2])) {pbcp=pp; bcpval=tmpval;}
      }
      if (tmpval<minval) {minval=tmpval; pmin=pp; for(int k=0; k<3; k++) {xmin[k]=xx[k];}}
      if (tmpval>maxval) {maxval=tmpval; pmax=pp; for(int k=0; k<3; k++) {xmax[k]=xx[k];}}
#if USEPROGRESSBAR
      printProgressBar(int(100.0e0*solreal(i)/solreal((nbgppts-1))));
#endif
   }
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   
   /* Closes the dat file */
   ofile.close();
   
   /* Displays the information of min/max of prop */
   
   cout << scientific << setprecision(12) << endl;
   printScrCharLine('-');
   cout << "The maximum value (within the evaluated line) is located at:" << endl;
   printV3Comp("R(maxval): ",xmax);
   cout << "p(maxval): " << pmax << endl;
   cout << getFieldTypeKeyShort(prop) << "(max): " << maxval << endl;
   printScrCharLine('-');
   cout << "The minimum value (within the evaluated line) is located at:" << endl;
   printV3Comp("R(minval): ",xmin);
   cout << "p(minval): " << pmin << endl;
   cout << getFieldTypeKeyShort(prop) << "(min): " << minval << endl;
   printScrCharLine('-');
   if (options.uponbp) {
      cout << "The value at BCP is:" << endl;
      printV3Comp("R(bcp): ",robcp);
      cout << "p(pcp): " << pbcp << endl;
      cout << getFieldTypeKeyShort(prop) << "(bcp): " << bcpval << endl;
   }
   printScrCharLine('-');
   
   /* Writes the information of the field property to the log file */

   ofstream logfil;

   logfil.open(lognam.c_str(),ios::out);
   writeCommentedHappyStart(argv,logfil,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   logfil << scientific << setprecision(12);
   writeCommentedScrCharLine(logfil,'-');
   logfil << "#The maximum value (within the evaluated line) is located at:" << endl;
   writeV3Comp(logfil,"#R(maxval):\n",xmax);
   logfil << "#p(maxval):\n" << pmax << endl;
   logfil << "#" << getFieldTypeKeyShort(prop) << "(max):\n" << maxval << endl;
   writeCommentedScrCharLine(logfil,'-');
   logfil << "#The minimum value (within the evaluated line) is located at:" << endl;
   writeV3Comp(logfil,"#R(minval):\n",xmin);
   logfil << "#p(minval):\n" << pmin << endl;
   logfil << "#" << getFieldTypeKeyShort(prop) << "(min):\n" << minval << endl;
   writeCommentedScrCharLine(logfil,'-');
   if (options.uponbp) {
      logfil << "#The value at BCP is:" << endl;
      writeV3Comp(logfil,"R(bcp):\n",robcp);
      logfil << "#p(pcp):\n" << pbcp << endl;
      logfil << getFieldTypeKeyShort(prop) << "(bcp):\n" << bcpval << endl;
      writeCommentedScrCharLine(logfil,'-');
   }
   
   logfil.close();

   /* Writes the gnuplot script */
   
   HelpersPlot::MakeGnuplotFile(options,gnpnam,outfilnam,prop,bnw,lenline,minval,maxval,\
         at1,at2,pbcp,rbgp);
   
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   if (options.zipdat) {
      string cmdl;
      cmdl=string("gzip -9f ")+outfilnam;
      cout << "Calling gzip...";
      system(cmdl.c_str());
      cout << " Done!" << endl;
   }
#endif/* defined(__APPLE__)||defined(__linux__) */
   
   dealloc2DRealArray(rbgp,dimarr);
   
   /* At this point the computation has ended. Usually this means that no errors ocurred. */
   
   setScrGreenBoldFont();
   printHappyEnding();
   printScrStarLine();
   cout << setprecision(3) << "CPU Time: "
        << solreal( clock () - begin_time ) / CLOCKS_PER_SEC << "s" << endl;
   solreal end_walltime=time(NULL);
   cout << "Wall-clock time: " << solreal (end_walltime-begin_walltime) << "s" << endl;
#if DEBUG
   cout << "Debuggin mode (under construction...)" << endl;
#endif
   printScrStarLine();
   setScrNormalFont();
   return 0;
}


