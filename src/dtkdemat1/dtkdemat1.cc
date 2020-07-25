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
 ------------------------
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
#include "../common/mymath.h"
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "../common/demat1critptnetworksl.h"
#include "../common/stringtools.h"
#include "../common/gnuplottools.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helperplots.h"


/* computeUVProjection evaluates the projection of $\nabla\gamma(x,x')$, $\nabla'\gamma(x,x')$
 * upon the coordinates (u,v)  */
void computeUVProjection(solreal (&x1)[3],solreal (&x2)[3],\
      solreal (&g)[3],solreal (&gp)[3],solreal (&uv)[2]);



int main (int argc, char ** argv) {
   const clock_t begin_time = clock();
   const solreal begin_walltime = time(NULL);
   string infilnam,outfilnam,o1dfilnam,o1sfilnam,basegnpnam,lognam;
   string progname;
   optFlags options;
   ifstream ifile;
   ofstream ofile,o1dfile,o1sfile;
   
   getOptions(argc,argv,options); //This processes the options from the command line.
   mkFileNames(argv,options,infilnam,outfilnam,o1dfilnam,o1sfilnam,basegnpnam,lognam); //This creates the names used.
   ScreenUtils::PrintHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   /* Loading the wave function */
   cout << endl << "Loading wave function from file: " << infilnam << "... ";
   GaussWaveFunction gwf;
   if (!(gwf.readFromFile(infilnam))) { //Loading the wave function
      ScreenUtils::SetScrRedBoldFont();
      cout << "Error: the wave function could not be loaded!\n";
      ScreenUtils::SetScrNormalFont();
      exit(1);
   }
   cout << "Done." << endl;
   if (gwf.nNuc==1) {
      ScreenUtils::DisplayWarningMessage("This file contains only one atom... There are no bond paths...");
      cout << "Nothing to do!" << endl;
      gwf.~GaussWaveFunction();
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
      ScreenUtils::DisplayErrorMessage("Requesting a non existent atom!");
      gwf.~GaussWaveFunction();
      exit(1);
   }

   char prop='D';
   if ( options.prop2plot ) {prop=argv[options.prop2plot][0];}
   
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
         ScreenUtils::DisplayErrorMessage("Please provide a positive number for the number of points!");
         gwf.~GaussWaveFunction();
         bnw.~bondNetWork();
         cpn.~critPtNetWork();
         exit(1);
      }
   }
   
   /* Setting the final file names. */
   string lbls="-"+gwf.atLbl[at1]+"-"+gwf.atLbl[at2];
   size_t pos=basegnpnam.find_last_of('.');
   if (pos!=string::npos) {
      basegnpnam.insert(pos,lbls);
      outfilnam.insert(pos,lbls);
      lognam.insert(pos,lbls);
      o1dfilnam.insert(pos,lbls);
      o1sfilnam.insert(pos,lbls);
   }
   pos=o1dfilnam.find_last_of('.');
   if (pos!=string::npos) {
      o1dfilnam.insert(pos,"-1D1");
      o1sfilnam.insert(pos,"-1D2");
   }
   
   /* Declaring an array to store the coordinates of the bond path points. */
   solreal **rbgp=NULL;
   MyMemory::Alloc2DRealArray(string("rbgp"),dimarr,3,rbgp);
   
   int nbgppts;
   solreal dl=DEFAULTBONDPATHSTEPMD1;
   solreal robcp[3];
   solreal x1[3],x2[3],dx[3];
   solreal xbeg[3],xend[3];
   solreal lenline=0.0e0;
   solreal dist=0.0e0;
   
   if (options.setstep) {sscanf(argv[options.setstep],"%lg",&dl);}
   
   /* ==================================================================== */
   
   if (options.uponbp) {
      /* Compute the gradient path */
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
      lenline*=2.0e0;
      dl=lenline/solreal(dimarr-1);
      for (int i=0; i<3; i++) {
         dx[i]=2.0e0*(bnw.R[at2][i]-bnw.R[at1][i])/solreal(dimarr-1);
         rbgp[0][i]=0.5e0*(3.0e0*bnw.R[at1][i]-bnw.R[at2][i]);
         robcp[i]=0.0e0;
      }
      for (int i=1; i<dimarr; i++) {for (int j=0; j<3; j++) {rbgp[i][j]=rbgp[i-1][j]+dx[j];}}
      cout << "The line that joins the atoms consists of " << nbgppts << " points." << endl;
   }

   for ( int i=0 ; i<3 ; i++ ) {
      xbeg[i]=rbgp[0][i];
      xend[i]=rbgp[nbgppts-1][i];
   }
   
   
   /* Open the tsv and dat files */
   
   ofile.open(outfilnam.c_str(),ios::out);
   o1dfile.open(o1dfilnam.c_str(),ios::out);
   o1sfile.open(o1sfilnam.c_str(),ios::out);
   ofile << scientific << setprecision(12);
   o1dfile << scientific << setprecision(12);
   o1sfile << scientific << setprecision(12);
   cout << scientific << setprecision(12);
   
   solreal p1,p2,pbcp,xt[3];
   solreal md1tmp,md1max=-1.0e+50,md1min=1.0e+50,rhomin=1.0e+50,diagmax=-1.0e+50;
   solreal x1max[3],x1min[3],x2max[3],x2min[3],p1max,p1min,p2max,p2min,p1dmax,p2dmax;
   solreal xrmin[3],xd1max[3],xd2max[3];
   solreal gmd1max=-1.0e+50,gmd1min=1.0e+50;
   solreal ggradmagmax=-1.0e+50,ggradmagmin=1.0e+50;
   
   /* Evaluating and writing the density matrix of order 1 into the tsv file
      On the fly, determining min/max of MD1 and the respective
      coordinates of such min/max. The same for the global min/max.
   */
   
   cout << (string("Evaluating ")+string(prop=='G' ? "Grad" : "the ")+string("MD1 "));
   if (options.uponbp) {cout << "upon the bond path... " << endl;}
   if (options.uponsl) {cout << "upon the straight line that joins the selected atoms..." << endl;}
   cout << "Progress: " << endl;
   
   if ( options.centredats ) {
      p1=-0.5e0*lenline;
      p2=-0.5e0*lenline;
   } else {
      p1=0.0e0;
      p2=0.0e0;
   }
   solreal gg[3],gp[3],proj[2],magproj;
   
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   /* Upon the bond path... */
   if (options.uponbp) {
      for (int i=0; i<nbgppts; i++) {
         for (int k=0; k<3; k++) {x1[k]=rbgp[i][k];}
         if ((rbgp[i][0]==robcp[0])&&(rbgp[i][1]==robcp[1])&&(rbgp[i][2]==robcp[2])) {pbcp=p1;}
         for (int k=0; k<3; k++) {x2[k]=rbgp[0][k];}
         if (options.centredats) { p2=-0.5e0*lenline; } else { p2=0.0e0; }
         if ( prop=='G' ) {
            gwf.evalGradDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],md1tmp,gg,gp);
         } else {
            md1tmp=gwf.evalDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
         }
         if (md1tmp>md1max) {
            md1max=md1tmp;
            for (int k=0; k<3; k++) {
               x1max[k]=x1[k];
               x2max[k]=x2[k];
            }
            p1max=p1;
            p2max=p2;
         }
         if (md1tmp<md1min) {
            md1min=md1tmp;
            for (int k=0; k<3; k++) {
               x1min[k]=x1[k];
               x2min[k]=x2[k];
            }
            p1min=p1;
            p2min=p2;
         }
         if ( md1tmp>gmd1max ) {gmd1max=md1tmp;}
         if ( md1tmp<gmd1min ) {gmd1min=md1tmp;}
         ofile << p1 << "\t" << p2 << "\t" << md1tmp;
         if ( prop=='G' ) {
            computeUVProjection(xbeg,xend,gg,gp,proj);
            magproj=(proj[0]*proj[0]+proj[1]*proj[1]);
            if ( magproj>ggradmagmax ) {ggradmagmax=magproj;}
            if ( magproj<ggradmagmin ) {ggradmagmin=magproj;}
            for ( int ss=0 ; ss<2 ; ss++ ) {ofile << "\t" << proj[ss];}
            ofile << endl;
         } else {
            ofile << endl;
         }
         if (i==0) {
            o1dfile << p1 << " " << md1tmp;
            if ( prop=='G' ) {
               o1dfile << endl;
            } else {
               o1dfile << endl;
            }
         }
         for (int j=1; j<nbgppts; j++) {
            dist=0.0e0;
            for (int k=0; k<3; k++) {
               xt[k]=x2[k];
               x2[k]=rbgp[j][k];
               dist+=((x2[k]-xt[k])*(x2[k]-xt[k]));
            }
            p2+=sqrt(dist);
            if ( prop=='G' ) {
               gwf.evalGradDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],md1tmp,gg,gp);
            } else {
               md1tmp=gwf.evalDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
            }
            if (i==j) {
               if (rhomin>md1tmp) {
                  rhomin=md1tmp;
                  for (int k=0; k<3; k++) {xrmin[k]=x1[k];}
               }
               o1dfile << p1 << " " << md1tmp;
               if ( prop == 'G' ) {
                  o1dfile << endl;
               } else {
                  o1dfile << endl;
               }
            }
            if (i==(nbgppts-j)) {
               if (diagmax<md1tmp) {
                  diagmax=md1tmp;
                  for (int k=0; k<3; k++) {
                     xd1max[k]=x1[k];
                     xd2max[k]=x2[k];
                  }
                  p1dmax=p1;
                  p2dmax=p2;
               }
            }
            if ( md1tmp>gmd1max ) {gmd1max=md1tmp;}
            if ( md1tmp<gmd1min ) {gmd1min=md1tmp;}
            ofile << p1 << "\t" << p2 << "\t" << md1tmp;
            if ( prop=='G' ) {
               computeUVProjection(xbeg,xend,gg,gp,proj);
               magproj=(proj[0]*proj[0]+proj[1]*proj[1]);
               if ( magproj>ggradmagmax ) {ggradmagmax=magproj;}
               if ( magproj<ggradmagmin ) {ggradmagmin=magproj;}
               for ( int ss=0 ; ss<2 ; ss++ ) {ofile << "\t" << proj[ss];}
               ofile << endl;
            } else {
               ofile << endl;
            }
            if (j==(nbgppts-1-i)) {
               o1sfile << p1 << " " << md1tmp;
               if ( prop=='G' ) {
                  o1sfile << endl;
               } else {
                  o1sfile << endl;
               }
            }
            if (md1tmp>md1max) {
               md1max=md1tmp;
               for (int k=0; k<3; k++) {
                  x1max[k]=x1[k];
                  x2max[k]=x2[k];
               }
               p1max=p1;
               p2max=p2;
            }
            if (md1tmp<md1min) {
               md1min=md1tmp;
               for (int k=0; k<3; k++) {
                  x1min[k]=x1[k];
                  x2min[k]=x2[k];
               }
               p1min=p1;
               p2min=p2;
            }
         }
         ofile << endl;
         dist=0.0e0;
         if (i<nbgppts) {
            for (int k=0; k<3; k++) {
               xt[k]=rbgp[i+1][k];
               dist+=((x1[k]-xt[k])*(x1[k]-xt[k]));
            }
         }
         p1+=sqrt(dist);
#if USEPROGRESSBAR
         ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((nbgppts-1))));
#endif
      }
   }
   /* Upon the line that joins the atoms */
   if (options.uponsl) {
      for (int i=0; i<nbgppts; i++) {
         for (int k=0; k<3; k++) {x1[k]=rbgp[i][k];}
         for (int k=0; k<3; k++) {x2[k]=rbgp[0][k];}
         if (options.centredats) { p2=-0.5e0*lenline; } else { p2=0.0e0; }
         if ( prop=='G' ) {
            gwf.evalGradDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],md1tmp,gg,gp);
         } else {
            md1tmp=gwf.evalDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
         }
         if (md1tmp>md1max) {
            md1max=md1tmp;
            for (int k=0; k<3; k++) {
               x1max[k]=x1[k];
               x2max[k]=x2[k];
            }
            p1max=p1;
            p2max=p2;
         }
         if (md1tmp<md1min) {
            md1min=md1tmp;
            for (int k=0; k<3; k++) {
               x1min[k]=x1[k];
               x2min[k]=x2[k];
            }
            p1min=p1;
            p2min=p2;
         }
         if ( md1tmp>gmd1max ) {gmd1max=md1tmp;}
         if ( md1tmp<gmd1min ) {gmd1min=md1tmp;}
         ofile << p1 << "\t" << p2 << "\t" << md1tmp;
         if ( prop=='G' ) {
            computeUVProjection(xbeg,xend,gg,gp,proj);
            magproj=(proj[0]*proj[0]+proj[1]*proj[1]);
            if ( magproj>ggradmagmax ) {ggradmagmax=magproj;}
            if ( magproj<ggradmagmin ) {ggradmagmin=magproj;}
            for ( int ss=0 ; ss<2 ; ss++ ) {ofile << "\t" << proj[ss];}
            ofile << endl;
         } else {
            ofile << endl;
         }
         if (i==0) {
            o1dfile << p1 << " " << md1tmp;
            if ( prop=='G' ) {
               o1dfile << endl;
            } else {
               o1dfile << endl;
            }
         }
         for (int j=1; j<nbgppts; j++) {
            for (int k=0; k<3; k++) {x2[k]=rbgp[j][k];}
            p2+=dl;
            if ( prop=='G' ) {
               gwf.evalGradDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],md1tmp,gg,gp);
            } else {
               md1tmp=gwf.evalDensityMatrix1(x1[0],x1[1],x1[2],x2[0],x2[1],x2[2]);
            }
            if (i==j) {
               if (rhomin>md1tmp) {
                  rhomin=md1tmp;
                  for (int k=0; k<3; k++) {xrmin[k]=x1[k];}
               }
               o1dfile << p1 << " " << md1tmp;
               if ( prop == 'G' ) {
                  o1dfile << endl;
               } else {
                  o1dfile << endl;
               }
            }
            if (i==(nbgppts-j)) {
               if (diagmax<md1tmp) {
                  diagmax=md1tmp;
                  for (int k=0; k<3; k++) {
                     xd1max[k]=x1[k];
                     xd2max[k]=x2[k];
                  }
                  p1dmax=p1;
                  p2dmax=p2;
               }
            }
            if (j==(nbgppts-1-i)) {
               o1sfile << p1 << " " << md1tmp;
               if ( prop=='G' ) {
                  o1sfile << endl;
               } else {
                  o1sfile << endl;
               }
            }
            if ( md1tmp>gmd1max ) {gmd1max=md1tmp;}
            if ( md1tmp<gmd1min ) {gmd1min=md1tmp;}
            ofile << p1 << "\t" << p2 << "\t" << md1tmp;
            if ( prop == 'G' ) {
               computeUVProjection(xbeg,xend,gg,gp,proj);
               magproj=(proj[0]*proj[0]+proj[1]*proj[1]);
               if ( magproj>ggradmagmax ) {ggradmagmax=magproj;}
               if ( magproj<ggradmagmin ) {ggradmagmin=magproj;}
               for ( int ss=0 ; ss<2 ; ss++ ) {ofile << "\t" << proj[ss];}
               ofile << endl;
            } else {
               ofile << endl;
            }
            if (md1tmp>md1max) {
               md1max=md1tmp;
               for (int k=0; k<3; k++) {
                  x1max[k]=x1[k];
                  x2max[k]=x2[k];
               }
               p1max=p1;
               p2max=p2;
            }
            if (md1tmp<md1min) {
               md1min=md1tmp;
               for (int k=0; k<3; k++) {
                  x1min[k]=x1[k];
                  x2min[k]=x2[k];
               }
               p1min=p1;
               p2min=p2;
            }
         }
         ofile << endl;
         p1+=dl;
#if USEPROGRESSBAR
         ScreenUtils::PrintProgressBar(int(100.0e0*solreal(i)/solreal((nbgppts-1))));
#endif
      }
   }
   
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(100);
   cout << endl;
#endif
   
   /* Close the tsv and dat files */
   ofile.close();
   o1dfile.close();
   o1sfile.close();
   
   /* Finding critical points (if requested)  */

   DeMat1CriticalPointNetworkSL cp(&gwf,at1,at2);
   if ( options.findcps ) {
      ScreenUtils::PrintBetweenStarLines("Searching Critical Points");
      cp.setGammaCriticalPoints();
      ScreenUtils::PrintScrStarLine();
   }

   /* Display the information of min/max of MD1 */
   
   cout << scientific << setprecision(12) << endl;
   ScreenUtils::PrintScrCharLine('-');
   //cout << "md1max: " << md1max << ", md1min: " << md1min << endl << endl;
   cout << "The maximum value (within the evaluated line) is located at:" << endl;
   ScreenUtils::PrintV3Comp("x1(max): ",x1max);
   ScreenUtils::PrintV3Comp("x2(max): ",x2max);
   cout << "(p1: " << p1max << ")\n(p2: " << p2max << ")" << endl;
   cout << "MD1(max): " << gwf.evalDensityMatrix1(x1max[0],x1max[1],x1max[2],
                                             x2max[0],x2max[1],x2max[2]);
   cout << endl;
   ScreenUtils::PrintScrCharLine('-');
   cout << "The minimum value (within the evaluated line) is located at:" << endl;
   ScreenUtils::PrintV3Comp("x1(min): ",x1min);
   ScreenUtils::PrintV3Comp("x2(min): ",x2min);
   cout << "(p1: " << p1min << ")\n(p2: " << p2min << ")" << endl;
   cout << "MD1(min): " << gwf.evalDensityMatrix1(x1min[0],x1min[1],x1min[2],
                                             x2min[0],x2min[1],x2min[2]);
   cout << endl;
   ScreenUtils::PrintScrCharLine('-');
   solreal md1lmin;
   if (options.uponbp) {
      md1lmin=gwf.evalDensityMatrix1(robcp[0],robcp[1],robcp[2],robcp[0],robcp[1],robcp[2]);
      cout << "The bond critical point is located at:" << endl;
      ScreenUtils::PrintV3Comp("R(BCP): ",robcp);
      cout << "p1(p2): " << pbcp << endl;
      cout << "MD1(BCP): " << md1lmin << endl;
   }
   if (options.uponsl) {
      cout << "The minimum value of the density is located at" << endl;
      cout << "(considering only the evaluated points)" << endl;
      md1lmin= gwf.evalDensityMatrix1(xrmin[0],xrmin[1],xrmin[2],xrmin[0],xrmin[1],xrmin[2]);
      ScreenUtils::PrintV3Comp("R(rho,min): ",xrmin);
      cout << "MD1(rho,min): " << md1lmin << endl;
   }
   ScreenUtils::PrintScrCharLine('-');
   solreal md1dmax;
   md1dmax=gwf.evalDensityMatrix1(xd1max[0],xd1max[1],xd1max[2],xd2max[0],xd2max[1],xd2max[2]);
   cout << "The maximum value of MD1 at the diagonal (90 degrees from the rho line)" << endl;
   cout << "is located at" << endl;
   ScreenUtils::PrintV3Comp("x1(diag,max): ",xd1max);
   ScreenUtils::PrintV3Comp("x2(diag,max): ",xd2max);
   cout << "(p1: " << p1dmax << ")\n(p2: " << p2dmax << ")" << endl;
   cout << "MD1(diag,max): " << md1dmax << endl;
   if ( options.findcps ) {
      ScreenUtils::PrintScrCharLine('-');
      cp.displayCPsInfo();
   }
   ScreenUtils::PrintScrCharLine('-');
   cout << "The value of MD1 at the point e --aka cicp-- is:" << endl;
   cout << "MD1(cicp): " << gwf.evalDensityMatrix1(\
         bnw.R[at1][0],bnw.R[at1][1],bnw.R[at1][2],\
         bnw.R[at2][0],bnw.R[at2][1],bnw.R[at2][2]) << endl;
   ScreenUtils::PrintScrCharLine('-');
   cout << "Logfile: " << lognam << endl;
   ScreenUtils::PrintScrCharLine('-');
   cout << endl;
   

   /* Writing the information of MD1 (min/max) to the log file */
   
   ofstream logfil;
   
   logfil.open(lognam.c_str(),ios::out);
   FileUtils::WriteHappyStart(argv,logfil,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   logfil << scientific << setprecision(12);
   FileUtils::WriteScrCharLine(logfil,'-');
   logfil << "#The maximum value (within the evaluated line) is located at:" << endl;
   FileUtils::WriteV3Components(logfil,"#x1(max):\n",x1max);
   FileUtils::WriteV3Components(logfil,"#x2(max):\n",x2max);
   logfil << "#(p1):\n" << p1max << "\n#(p2):\n" << p2max << endl;
   logfil << "#MD1(max):\n" << gwf.evalDensityMatrix1(x1max[0],x1max[1],x1max[2],
                                                    x2max[0],x2max[1],x2max[2]);
   logfil << endl;
   FileUtils::WriteScrCharLine(logfil,'-');
   logfil << "#The minimum value (within the evaluated line) is located at:" << endl;
   FileUtils::WriteV3Components(logfil,"#x1(min):\n",x1min);
   FileUtils::WriteV3Components(logfil,"#x2(min):\n",x2min);
   logfil << "#(p1):\n" << p1min << "\n#(p2):\n" << p2min << endl;
   logfil << "#MD1(min):\n" << gwf.evalDensityMatrix1(x1min[0],x1min[1],x1min[2],
                                                    x2min[0],x2min[1],x2min[2]);
   logfil << endl;
   FileUtils::WriteScrCharLine(logfil,'-');
   if (options.uponbp) {
      logfil << "#The bond critical point is located at:" << endl;
      FileUtils::WriteV3Components(logfil,"#R(BCP):\n",robcp);
      logfil << "#p1(p2):\n" << pbcp << endl;
      logfil << "#MD1(BCP):\n" << md1lmin << endl;
   }
   if (options.uponsl) {
      logfil << "#The minimum value of the density is located at" << endl;
      logfil << "#(considering only the evaluated points)" << endl;
      FileUtils::WriteV3Components(logfil,"#R(rho,min):\n",xrmin);
      logfil << "#MD1(rho,min):\n" << md1lmin << endl;
   }
   FileUtils::WriteScrCharLine(logfil,'-');
   logfil << "#The maximum value of MD1 at the diagonal (90 degrees from the rho line)" << endl;
   logfil << "#is located at" << endl;
   FileUtils::WriteV3Components(logfil,"#x1(diag,max):\n",xd1max);
   FileUtils::WriteV3Components(logfil,"#x2(diag,max):\n",xd2max);
   logfil << "#(p1):\n" << p1dmax << "\n#(p2):\n" << p2dmax << endl;
   logfil << "#MD1(diag,max):\n" << md1dmax << endl;
   /* Writing the information of MD1 (min/max) to the log file */
   if ( options.findcps ) {
      FileUtils::WriteScrCharLine(logfil,'-');
      cp.writeCPsInfo(logfil);
   }
   /* Writes the value of gamma at cicp  */
   FileUtils::WriteScrCharLine(logfil,'-');
   logfil << "#The value of MD1 at the point e --aka cicp-- is:" << endl;
   logfil << "#MD1(CICP): " << endl << gwf.evalDensityMatrix1(\
         bnw.R[at1][0],bnw.R[at1][1],bnw.R[at1][2],\
         bnw.R[at2][0],bnw.R[at2][1],bnw.R[at2][2]) << endl;
   FileUtils::WriteScrCharLine(logfil,'-');

   logfil.close();
   
   /* Makes the plots */
   
   string line,tmpnam;
   solreal minval,maxval;
   solreal range=fabs(md1dmax-md1lmin);
#if DEBUG
   if ( md1lmin<0.0e0 ) {
      ScreenUtils::DisplayWarningMessage(string("md1lmin: "+getStringFromReal(md1lmin)));
   }
   cout << "gmin: " << gmd1min << ", gmax: " << gmd1max << endl;
   cout << "md1dmax: " << md1dmax << endl;
   cout << "md1lmin: " << md1lmin << endl;
   cout << "range: " << range << endl;
#endif
   minval=round(110*gmd1min)/100.0e0;
   maxval=round(110*(md1lmin+range))/100.0e0;
#if DEBUG
   cout << "minval: " << minval << endl;
   cout << "maxval: " << maxval << endl;
#endif /* ( DEBUG ) */
   if ((fabs(md1lmin-md1min)>range)&&(range<1.0e-02)) {
      range=2.0e0*fabs(md1lmin-md1min);
      minval=round(110*gmd1min)/100.0e0;
      maxval=round(110*(md1min+range))/100.0e0;
   }
#if DEBUG
   cout << "minval: " << minval << endl;
   cout << "maxval: " << maxval << endl;
#endif /* ( DEBUG ) */
   HelperPlot::generate3DPlot(options,outfilnam,minval,maxval,lenline,nbgppts);
   HelperPlot::generateMainDiagPlot(options,o1dfilnam,bnw,at1,at2,0.0e0,(maxval-minval),lenline,range);
   HelperPlot::generateSecDiagPlot(options,o1sfilnam,bnw,at1,at2,minval,maxval,lenline,range);
   HelperPlot::generateHeatMap(options,argv,outfilnam,bnw,cp,rbgp,nbgppts,minval,maxval,\
         lenline,md1lmin,md1dmax,at1,at2);
   if ( prop=='G' ) {
      HelperPlot::generateVectorField(options,argv,outfilnam,bnw,cp,rbgp,nbgppts,minval,maxval,ggradmagmin,ggradmagmax,\
            lenline,md1lmin,md1dmax,at1,at2);
   }
   
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   if (options.zipdat) {
      cout << "gzipping tsv..." << endl;
      line="gzip -9f "+outfilnam;
      system(line.c_str());
   }
#endif
   
   /* Freeing memory space */
   
   MyMemory::Dealloc2DRealArray(rbgp,dimarr);
   
   /* At this point the computation has ended. Usually this means that no errors ocurred. */
   
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


void computeUVProjection(solreal (&x1)[3],solreal (&x2)[3],\
      solreal (&g)[3],solreal (&gp)[3],solreal (&uv)[2])
{
   solreal tmp1=0.0e0,tmp2=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {
      tmp1+=((x2[i]-x1[i])*g[i]);
      tmp2+=((x2[i]-x1[i])*gp[i]);
   }
   uv[0]=tmp1;
   uv[1]=tmp2;
}

