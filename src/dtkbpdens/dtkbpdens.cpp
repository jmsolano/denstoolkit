/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.0
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
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
#include "../common/wavefunctionclass.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "../common/solstringtools.h"
#include "../common/fldtypesdef.h"
#include "optflags.h"
#include "crtflnms.h"

solreal evalFieldProperty(char prop,solreal (&x)[3],gaussWaveFunc &wf);


int main (int argc, char ** argv)
{
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
   
   /* Stablish the field property to be evaluated. */
   
   char prop;
   if (options.prop2plot) {
      prop=argv[options.prop2plot][0];
   } else {
      prop='d';
   }
   
   /* Checking for non valid property fields */
   
   if (!(prop=='d'||prop=='g'||prop=='l'||prop=='E'||prop=='L'||prop=='M'\
         ||prop=='S'||prop=='G'||prop=='K'||prop=='V'\
         ||prop=='P'||prop=='r'||prop=='s'||prop=='u')) {
      displayErrorMessage("Non valid field type");
      cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
      exit(1);
   }
   
   /* Loading the wave function */
   
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   
   gaussWaveFunc gwf;
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
      gwf.~gaussWaveFunc();
      exit(0);
   }
   
   /* Setting the atoms for defining the line/bond path */
   
   int at1,at2;
   at1=0;
   at2=1;
   if (options.setats) {
      sscanf(argv[options.setats],"%d",&at1);
      at1--;
      sscanf(argv[options.setats+1],"%d",&at2);
      at2--;
   }
   if (at1>=gwf.nNuc||at2>=gwf.nNuc||at1<0||at2<0) {
      displayErrorMessage("Requesting a non existent atom!");
      gwf.~gaussWaveFunc();
      exit(1);
   }
   
   
   /* Setting the bond network of the molecule */
   
   bondNetWork bnw;
   bnw.readFromFile(infilnam); //Loading the bond-network (if the wave function
                               //was read, there souldn't be problems here.
   bnw.setUpBNW();             //To setup the bond network.
   
   /* Defining the main critical point network object. */
   
   critPtNetWork cpn(gwf,bnw);
   
   
   
   int dimarr=200;
   if (options.setn1) {
      sscanf(argv[options.setn1],"%d",&dimarr);
      if (dimarr<=0) {
         displayErrorMessage("Please provide a positive number for the number of points!");
         gwf.~gaussWaveFunc();
         bnw.~bondNetWork();
         cpn.~critPtNetWork();
         exit(1);
      }
   }
   
   /* Setting the final file names. */
   
   string lbls="-"+gwf.atLbl[at1]+"-"+gwf.atLbl[at2];
   size_t pos=gnpnam.find_last_of('.');
   if (pos!=string::npos) {
      gnpnam.insert(pos,lbls);
      outfilnam.insert(pos,lbls);
      lognam.insert(pos,lbls);
   }
   
   /* Declaring an array to store the coordinates of the bond path points. */
   solreal **rbgp=NULL;
   alloc2DRealArray(string("rbgp"),dimarr,3,rbgp);
   
   int nbgppts;
   solreal dl=DEFAULTBONDPATHSTEPMD1;
   solreal robcp[3];
   solreal x1[3],x2[3],dx[3];
   solreal lenline=0.0e0;
   solreal dist=0.0e0;
   
   if (options.setstep) {sscanf(argv[options.setstep],"%lg",&dl);}
   
   /* ==================================================================== */
   
   if (options.uponbp) {
      /* Compute the gradient path */
      nbgppts=cpn.findSingleRhoGradientPathRK5(at1,at2,dl,dimarr,rbgp,robcp);
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
   
   /* Open the dat file */
   
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
   tmpval=evalFieldProperty(prop,xx,gwf);
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
      tmpval=evalFieldProperty(prop,xx,gwf);
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
   
   /* Close the dat file */
   ofile.close();
   
   /* Display the information of min/max of prop */
   
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
   
   /* Writing the information of the field property to the log file */
   
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
   
   /* Writing the gnuplot script */
   
   ofstream gfil;
   string line;
   gfil.open(gnpnam.c_str(),ios::out);
   
   if (maxval>1000.e0) {maxval=1000.0e0;}
   
   gfil << "set datafile missing \"inf\"" << endl;
   gfil << "namedatfile='" << outfilnam << "'" << endl;
   gfil << "xmin2plot=0.0" << endl;
   gfil << "xmax2plot=" << lenline << endl;
   gfil << "ymin2plot=" << (0.0e0<minval ? 0.0e0 : minval) << endl;
   gfil << "ymax2plot=" << maxval*1.05e0 << endl;
   gfil << "set xrange [xmin2plot:xmax2plot]" << endl;
   gfil << "set yrange [ymin2plot:ymax2plot]" << endl;
   if (!options.showatlbls) {gfil << "#";}
   solreal tmpdist=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {tmpdist+=((bnw.R[at1][i]-rbgp[0][i])*(bnw.R[at1][i]-rbgp[0][i]));}
   tmpdist=sqrt(tmpdist);
   int tmpa1=at1,tmpa2=at2;
   if ( tmpdist>0.1 ) {tmpa1=at2; tmpa2=at1;}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bnw.atLbl[tmpa1]) << "' at "
   << 0.05e0*(lenline) << "," << 0.1e0*(maxval) << " front offset character 0,0.75" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bnw.atLbl[tmpa2]) << "' at "
   << 0.95e0*(lenline) << "," << 0.1e0*(maxval) << " front offset character -2,0.75" << endl;
   gfil << "set arrow 1 from " << 0.05e0*(lenline) << "," << 0.1e0*(maxval) << " to "
        << "0.0,ymin2plot" << endl;
   gfil << "set arrow 2 from " << 0.95e0*(lenline) << "," << 0.1e0*(maxval) << " to "
        << lenline << ",ymin2plot" << endl;
   if (options.uponbp) {
      if (!options.showatlbls) {gfil << "#";}
      gfil << "set label 3 'BCP' at "
      << pbcp << "," << 0.05e0*(maxval) << " rotate by 90 front" << endl;
      if (!options.showatlbls) {gfil << "#";}
      gfil << "set xtics add ('' " << pbcp << ")" << endl;
   }
   gfil << "set title '" << getEnhancedEpsTitle(outfilnam) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   string epsnam=gnpnam.substr(0,(gnpnam.length()-3));
   string pdfnam=epsnam;
   epsnam.append("eps");
   pdfnam.append("pdf");
   gfil << "set output '" << epsnam << "'" << endl;
   gfil << "plot namedatfile using 1:5 with lines lw 2 title '"
        << gnuplotFieldTitle(prop) << "'" << endl;
   gfil << "#" << endl;
   writeScrCharLine(gfil,'#');
   gfil << "#                 END OF GNUPLOT COMMANDS" << endl;
   writeScrCharLine(gfil,'#');
   gfil << "#If you want to reconstruct the plots using this file, type:" << endl
   << "#gnuplot " << gnpnam << endl
   << "#epstopdf --outfile=" << pdfnam << " " << epsnam << endl;
   gfil.close();
   
#if _HAVE_GNUPLOT_
   if (options.mkplt) {
      cout << "Preparing the plot..." << endl;
#if (defined(__APPLE__)||defined(__linux__))
      line="gnuplot "+gnpnam;
#endif
#if (defined(__CYGWIN__))
      line="wgnuplot "+gnpnam;
#endif
      if (options.quiet) {line+=" > /dev/null 2>&1";}
      system(line.c_str());
#if _HAVE_EPSTOPDF_
      line= "epstopdf "+epsnam+" --outfile="+pdfnam;
      if (options.quiet) {line+=" > /dev/null 2>&1";}
      system(line.c_str());
      line= "rm -f "+epsnam;
      if (options.quiet) {line+=" > /dev/null 2>&1";}
      system(line.c_str());
#endif
      if (!(options.kpgnp)) {
         line="rm -f "+gnpnam;
         if (options.quiet) {line+=" > /dev/null 2>&1";}
         system(line.c_str());
      }
   }
#endif
   
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

solreal evalFieldProperty(char prop,solreal (&x)[3],gaussWaveFunc &wf)
{
   solreal res;
   switch (prop) {
      case 'd':
         res=wf.evalDensity(x[0],x[1],x[2]);
         break;
      case 'g':
         res=wf.evalMagGradRho(x[0],x[1],x[2]);
         break;
      case 'l':
         res=wf.evalLapRho(x[0],x[1],x[2]);
         break;
      case 'E':
         res=wf.evalELF(x[0],x[1],x[2]);
         break;
      case 'L':
         res=wf.evalLOL(x[0],x[1],x[2]);
         break;
      case 'M':
         res=wf.evalMagGradLOL(x[0],x[1],x[2]);
         break;
      case 'P' :
         res=wf.evalMagLED(x[0],x[1],x[2]);
         break;
      case 'r' :
         res=wf.evalRoSE(x[0],x[1],x[2]);
         break;
      case 's' :
         res=wf.evalReducedDensityGradient(x[0],x[1],x[2]);
         break;
      case 'S':
         res=wf.evalShannonEntropy(x[0],x[1],x[2]);
         break;
      case 'G':
         res=wf.evalKineticEnergyG(x[0],x[1],x[2]);
         break;
      case 'K':
         res=wf.evalKineticEnergyK(x[0],x[1],x[2]);
         break;
      case 'u' :
         res=wf.evalCustomScalarField(x[0],x[1],x[2]);
         break;
      case 'V':
         res=wf.evalMolElecPot(x[0],x[1],x[2]);
         break;
      default:
#if DEBUG
         cout << "Unknown field type! Returning zero..." << endl;
         DISPLAYDEBUGINFOFILELINE;
#endif
         res=0.0e0;
         break;
   }
   return res;
}



