/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
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

#include "../common/solscrutils.h"
#include "../common/solfileutils.h"
#include "../common/solmemhand.h"
#include "../common/solmath.h"
#include "../common/gausswavefunction.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"
#include "../common/demat1critptnetworksl.h"
#include "../common/solstringtools.h"
#include "../common/solgnuplottools.h"
#include "optflags.h"
#include "crtflnms.h"
#include "helperplots.h"


/* computeUVProjection evaluates the projection of $\nabla\gamma(x,x')$, $\nabla'\gamma(x,x')$
 * upon the coordinates (u,v)  */
void computeUVProjection(solreal (&x1)[3],solreal (&x2)[3],\
      solreal (&g)[3],solreal (&gp)[3],solreal (&uv)[2]);
/*
void generateMainDiagPlot(optFlags &options,const string &datname,\
      bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
      solreal linelength,solreal frange);
void generateSecDiagPlot(optFlags &options,const string &datname,\
      bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
      solreal linelength,solreal frange);
void generate3DPlot(optFlags &options,const string &tsvname,solreal minval2plot,\
      solreal maxval2plot,solreal linelength,int nptsinline);
void generateHeatMap(optFlags &options,char *argv[],const string &tsvname,\
      bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),int nptx,solreal minval2plot,\
      solreal maxval2plot,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
void generateVectorField(optFlags &options,char *argv[],const string &tsvname,\
      bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),int nptx,solreal minval2plot,\
      solreal maxval2plot,solreal maggradmin,solreal maggradmax,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
// */
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
   printHappyStart(argv,CURRENTVERSION,PROGRAMCONTRIBUTORS); //Just to let the user know that the initial configuration is OK
   
   /* Loading the wave function */
   
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
         displayErrorMessage("Please provide a positive number for the number of points!");
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
   alloc2DRealArray(string("rbgp"),dimarr,3,rbgp);
   
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
   
   p1=0.0e0;
   p2=0.0e0;
   solreal gg[3],gp[3],proj[2],magproj;
   
#if USEPROGRESSBAR
   printProgressBar(0);
#endif
   /* Upon the bond path... */
   if (options.uponbp) {
      for (int i=0; i<nbgppts; i++) {
         for (int k=0; k<3; k++) {x1[k]=rbgp[i][k];}
         if ((rbgp[i][0]==robcp[0])&&(rbgp[i][1]==robcp[1])&&(rbgp[i][2]==robcp[2])) {pbcp=p1;}
         for (int k=0; k<3; k++) {x2[k]=rbgp[0][k];}
         p2=0.0e0;
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
         printProgressBar(int(100.0e0*solreal(i)/solreal((nbgppts-1))));
#endif
      }
   }
   /* Upon the line that joins the atoms */
   if (options.uponsl) {
      for (int i=0; i<nbgppts; i++) {
         for (int k=0; k<3; k++) {x1[k]=rbgp[i][k];}
         for (int k=0; k<3; k++) {x2[k]=rbgp[0][k];}
         p2=0.0e0;
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
         printProgressBar(int(100.0e0*solreal(i)/solreal((nbgppts-1))));
#endif
      }
   }
   
#if USEPROGRESSBAR
   printProgressBar(100);
   cout << endl;
#endif
   
   /* Close the tsv and dat files */
   ofile.close();
   o1dfile.close();
   o1sfile.close();
   
   /* Finding critical points (if requested)  */

   DeMat1CriticalPointNetworkSL cp(&gwf,at1,at2);
   if ( options.findcps ) {
      printBetweenStarLines("Searching Critical Points");
      cp.setGammaCriticalPoints();
      printScrStarLine();
   }

   /* Display the information of min/max of MD1 */
   
   cout << scientific << setprecision(12) << endl;
   printScrCharLine('-');
   //cout << "md1max: " << md1max << ", md1min: " << md1min << endl << endl;
   cout << "The maximum value (within the evaluated line) is located at:" << endl;
   printV3Comp("x1(max): ",x1max);
   printV3Comp("x2(max): ",x2max);
   cout << "(p1: " << p1max << ")\n(p2: " << p2max << ")" << endl;
   cout << "MD1(max): " << gwf.evalDensityMatrix1(x1max[0],x1max[1],x1max[2],
                                             x2max[0],x2max[1],x2max[2]);
   cout << endl;
   printScrCharLine('-');
   cout << "The minimum value (within the evaluated line) is located at:" << endl;
   printV3Comp("x1(min): ",x1min);
   printV3Comp("x2(min): ",x2min);
   cout << "(p1: " << p1min << ")\n(p2: " << p2min << ")" << endl;
   cout << "MD1(min): " << gwf.evalDensityMatrix1(x1min[0],x1min[1],x1min[2],
                                             x2min[0],x2min[1],x2min[2]);
   cout << endl;
   printScrCharLine('-');
   solreal md1lmin;
   if (options.uponbp) {
      md1lmin=gwf.evalDensityMatrix1(robcp[0],robcp[1],robcp[2],robcp[0],robcp[1],robcp[2]);
      cout << "The bond critical point is located at:" << endl;
      printV3Comp("R(BCP): ",robcp);
      cout << "p1(p2): " << pbcp << endl;
      cout << "MD1(BCP): " << md1lmin << endl;
   }
   if (options.uponsl) {
      cout << "The minimum value of the density is located at" << endl;
      cout << "(considering only the evaluated points)" << endl;
      md1lmin= gwf.evalDensityMatrix1(xrmin[0],xrmin[1],xrmin[2],xrmin[0],xrmin[1],xrmin[2]);
      printV3Comp("R(rho,min): ",xrmin);
      cout << "MD1(rho,min): " << md1lmin << endl;
   }
   printScrCharLine('-');
   solreal md1dmax;
   md1dmax=gwf.evalDensityMatrix1(xd1max[0],xd1max[1],xd1max[2],xd2max[0],xd2max[1],xd2max[2]);
   cout << "The maximum value of MD1 at the diagonal (90 degrees from the rho line)" << endl;
   cout << "is located at" << endl;
   printV3Comp("x1(diag,max): ",xd1max);
   printV3Comp("x2(diag,max): ",xd2max);
   cout << "(p1: " << p1dmax << ")\n(p2: " << p2dmax << ")" << endl;
   cout << "MD1(diag,max): " << md1dmax << endl;
   if ( options.findcps ) {
      printScrCharLine('-');
      cp.displayCPsInfo();
   }
   printScrCharLine('-');
   cout << "The value of MD1 at the point e --aka cicp-- is:" << endl;
   cout << "MD1(cicp): " << gwf.evalDensityMatrix1(\
         bnw.R[at1][0],bnw.R[at1][1],bnw.R[at1][2],\
         bnw.R[at2][0],bnw.R[at2][1],bnw.R[at2][2]) << endl;
   printScrCharLine('-');
   cout << "Logfile: " << lognam << endl;
   printScrCharLine('-');
   cout << endl;
   

   /* Writing the information of MD1 (min/max) to the log file */
   
   ofstream logfil;
   
   logfil.open(lognam.c_str(),ios::out);
   writeCommentedHappyStart(argv,logfil,CURRENTVERSION,PROGRAMCONTRIBUTORS);
   logfil << scientific << setprecision(12);
   writeCommentedScrCharLine(logfil,'-');
   logfil << "#The maximum value (within the evaluated line) is located at:" << endl;
   writeV3Comp(logfil,"#x1(max):\n",x1max);
   writeV3Comp(logfil,"#x2(max):\n",x2max);
   logfil << "#(p1):\n" << p1max << "\n#(p2):\n" << p2max << endl;
   logfil << "#MD1(max):\n" << gwf.evalDensityMatrix1(x1max[0],x1max[1],x1max[2],
                                                    x2max[0],x2max[1],x2max[2]);
   logfil << endl;
   writeCommentedScrCharLine(logfil,'-');
   logfil << "#The minimum value (within the evaluated line) is located at:" << endl;
   writeV3Comp(logfil,"#x1(min):\n",x1min);
   writeV3Comp(logfil,"#x2(min):\n",x2min);
   logfil << "#(p1):\n" << p1min << "\n#(p2):\n" << p2min << endl;
   logfil << "#MD1(min):\n" << gwf.evalDensityMatrix1(x1min[0],x1min[1],x1min[2],
                                                    x2min[0],x2min[1],x2min[2]);
   logfil << endl;
   writeCommentedScrCharLine(logfil,'-');
   if (options.uponbp) {
      logfil << "#The bond critical point is located at:" << endl;
      writeV3Comp(logfil,"#R(BCP):\n",robcp);
      logfil << "#p1(p2):\n" << pbcp << endl;
      logfil << "#MD1(BCP):\n" << md1lmin << endl;
   }
   if (options.uponsl) {
      logfil << "#The minimum value of the density is located at" << endl;
      logfil << "#(considering only the evaluated points)" << endl;
      writeV3Comp(logfil,"#R(rho,min):\n",xrmin);
      logfil << "#MD1(rho,min):\n" << md1lmin << endl;
   }
   writeCommentedScrCharLine(logfil,'-');
   logfil << "#The maximum value of MD1 at the diagonal (90 degrees from the rho line)" << endl;
   logfil << "#is located at" << endl;
   writeV3Comp(logfil,"#x1(diag,max):\n",xd1max);
   writeV3Comp(logfil,"#x2(diag,max):\n",xd2max);
   logfil << "#(p1):\n" << p1dmax << "\n#(p2):\n" << p2dmax << endl;
   logfil << "#MD1(diag,max):\n" << md1dmax << endl;
   /* Writing the information of MD1 (min/max) to the log file */
   if ( options.findcps ) {
      writeCommentedScrCharLine(logfil,'-');
      cp.writeCPsInfo(logfil);
   }
   /* Writes the value of gamma at cicp  */
   writeCommentedScrCharLine(logfil,'-');
   logfil << "#The value of MD1 at the point e --aka cicp-- is:" << endl;
   logfil << "#MD1(CICP): " << endl << gwf.evalDensityMatrix1(\
         bnw.R[at1][0],bnw.R[at1][1],bnw.R[at1][2],\
         bnw.R[at2][0],bnw.R[at2][1],bnw.R[at2][2]) << endl;
   writeCommentedScrCharLine(logfil,'-');

   logfil.close();
   
   /* Makes the plots */
   
   string line,tmpnam;
   solreal minval,maxval;
   solreal range=fabs(md1dmax-md1lmin);
#if DEBUG
   if ( md1lmin<0.0e0 ) {
      displayWarningMessage(string("md1lmin: "+getStringFromReal(md1lmin)));
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
   generate3DPlot(options,outfilnam,minval,maxval,lenline,nbgppts);
   generateMainDiagPlot(options,o1dfilnam,bnw,at1,at2,minval,maxval,lenline,range);
   generateSecDiagPlot(options,o1sfilnam,bnw,at1,at2,minval,maxval,lenline,range);
   generateHeatMap(options,argv,outfilnam,bnw,cp,rbgp,minval,maxval,\
         lenline,md1lmin,md1dmax,at1,at2);
   generateVectorField(options,argv,outfilnam,bnw,cp,rbgp,minval,maxval,ggradmagmin,ggradmagmax,\
         lenline,md1lmin,md1dmax,at1,at2);
   
#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   if (options.zipdat) {
      cout << "gzipping tsv..." << endl;
      line="gzip -9f "+outfilnam;
      system(line.c_str());
   }
#endif
   
   /* Freeing memory space */
   
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
/*
void generateMainDiagPlot(optFlags &options,const string &datname,\
      bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
      solreal linelength,solreal frange)
{
   string gnpname=datname;
   replaceExtensionOfFileName(gnpname,string("gnp"));
   string pdfname=gnpname;
   replaceExtensionOfFileName(pdfname,string("pdf"));
   string epsname=gnpname;
   replaceExtensionOfFileName(epsname,string("eps"));
   //----------------------------------------------------
   ofstream gfil(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << datname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   gfil << "set xrange[0:dimparam]" << endl;
   gfil << "set yrange [minval2plot:maxval2plot]" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 from " << 0.15e0*(linelength) << "," << (minval2plot+0.15*frange)
   << " to " << 0.25e0*linelength << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 from " << 0.85e0*(linelength) << "," << (minval2plot+0.15*frange)
        << " to " << 0.75e0*linelength << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << (0.25*linelength) << ",minval2plot "
        << "to " << (0.25*linelength) << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << (0.75*linelength) << ",minval2plot "
        << "to " << (0.75*linelength) << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
   << 0.15e0*(linelength) << "," << (minval2plot+frange*0.15e0) << " front offset character -2,0.00" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << 0.85e0*(linelength) << "," << (minval2plot+frange*0.15e0) << " front offset character 0.3,0.00" << endl;
   gfil << "set title '" << getEnhancedEpsTitle(datname) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot '" << datname << "' w lines lw 2 notitle"  << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void generateSecDiagPlot(optFlags &options,const string &datname,\
      bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
      solreal linelength,solreal frange)
{
   string gnpname=datname;
   replaceExtensionOfFileName(gnpname,string("gnp"));
   string pdfname=gnpname;
   replaceExtensionOfFileName(pdfname,string("pdf"));
   string epsname=gnpname;
   replaceExtensionOfFileName(epsname,string("eps"));
   //----------------------------------------------------
   ofstream gfil(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << datname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   gfil << "set xrange[0:dimparam]" << endl;
   gfil << "set yrange [minval2plot:maxval2plot]" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 from " << 0.15e0*(linelength) << "," << (minval2plot+0.15*frange)
   << " to " << 0.25e0*linelength << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 from " << 0.85e0*(linelength) << "," << (minval2plot+0.15*frange)
        << " to " << 0.75e0*linelength << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << (0.25*linelength) << ",minval2plot "
        << "to " << (0.25*linelength) << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << (0.75*linelength) << ",minval2plot "
        << "to " << (0.75*linelength) << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
   << 0.15e0*(linelength) << "," << (minval2plot+frange*0.15e0) << " front offset character -2,0.00" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << 0.85e0*(linelength) << "," << (minval2plot+frange*0.15e0) << " front offset character 0.3,0.00" << endl;
   gfil << "set title '" << getEnhancedEpsTitle(datname) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot '" << datname << "' w lines lw 2 notitle"  << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}

void generate3DPlot(optFlags &options,const string &tsvname,\
      solreal minval2plot,solreal maxval2plot,solreal linelength,int nptsinline)
{
   string gnpname=tsvname;
   insertAtEndOfFileName(gnpname,string("-3D"));
   replaceExtensionOfFileName(gnpname,string("gnp"));
   string pdfname=gnpname;
   replaceExtensionOfFileName(pdfname,string("pdf"));
   string epsname=gnpname;
   replaceExtensionOfFileName(epsname,string("eps"));
   //----------------------------------------------------
   ofstream gfil(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << tsvname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   gfil << "set isosample 300, 300" << endl;
   gfil << "set cntrparam cubicspline" << endl;
   gfil << "set zrange [minval2plot:maxval2plot]" << endl;
   gfil << "set cbrange [minval2plot:maxval2plot]" << endl;
   gfil << "unset contour" << endl;
   gfil << "set pm3d depthorder hidden3d 1" << endl;
   gfil << "set palette rgbformulae 33,13,10" << endl;
   gfil << "set style line 1 linecolor rgb \"#444444\"" << endl;
   gfil << "unset title" << endl;
   gfil << "set xyplane 0" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "#set output '|epstopdf --filter --outfile=" << pdfname << "'" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "splot namedatfile ";
   if ((nptsinline>80&&nptsinline<=120)) {gfil << "every 2:2 ";}
   if ((nptsinline>120&&nptsinline<=160)) {gfil << "every 3:3 ";}
   if ((nptsinline>160&&nptsinline<240)) {gfil << "every 4:4 ";}
   if ((nptsinline>240&&nptsinline<300)) {gfil << "every 5:5 ";}
   gfil << "using 1:2:($3>maxval2plot? maxval2plot:($3<minval2plot ? minval2plot : $3)) "
   << "with pm3d title '" << getEnhancedEpsTitle(tsvname) << "'" << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void generateHeatMap(optFlags &options,char *argv[],const string &tsvname,bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,
      solreal **xx,int nptx,solreal minval2plot,\
      solreal maxval2plot,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2)
{
   string gnpname=tsvname;
   insertAtEndOfFileName(gnpname,string("-2D"));
   replaceExtensionOfFileName(gnpname,string("gnp"));
   string pdfname=gnpname;
   replaceExtensionOfFileName(pdfname,string("pdf"));
   string epsname=gnpname;
   replaceExtensionOfFileName(epsname,string("eps"));
   //----------------------------------------------------
   ofstream gfil;
   gfil.open(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << tsvname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   gfil << "set isosample 300, 300" << endl;
   gfil << "set cntrparam cubicspline" << endl;
   gfil << "set zrange [minval2plot:maxval2plot]" << endl;
   gfil << "set cbrange [minval2plot:maxval2plot]" << endl;
   gfil << "unset contour" << endl;
   gfil << "set pm3d depthorder hidden3d 1" << endl;
   gfil << "set palette rgbformulae 33,13,10" << endl;
   gfil << "set style line 1 linecolor rgb \"#444444\"" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "unset key" << endl;
   gfil << "set contour base" << endl;
   int noconts=DEFAULTNUMBEROFCONTOURLINES;
   solreal dcont=(maxval2plot-minval2plot)/solreal(noconts);
   solreal contval=md1lmin+0.25e0*dcont;
   gfil << "set cntrparam cubicspline" << endl;
   gfil << "set cntrparam levels ";
   if ( !(options.setinccont) ) {
      gfil << "discrete " << contval;
      gfil << ", " << (md1dmax-0.5*dcont);
      contval=minval2plot+dcont;
      for (int i=1; i<noconts; i++) {gfil << ", " << contval; contval+=dcont;}
      gfil << endl;
   } else {
      gfil << "incremental " << argv[options.setinccont] << ","
           << argv[options.setinccont+1] << "," << argv[options.setinccont+2] << endl;
   }
   gfil << "unset surface" << endl;
   gfil << "set table 'contourtemp.dat'" << endl;
   gfil << "splot namedatfile" << endl;
   gfil << "unset table" << endl;
   gfil << "set xrange[0:dimparam]" << endl;
   gfil << "set yrange [minval2plot:maxval2plot]" << endl;
   gfil << "unset title" << endl;
   gfil << "set yrange[0:dimparam]" << endl;
   gfil << "set tmargin at screen 0.98\nset bmargin at screen 0.01\nset lmargin at screen 0.02" << endl;
   gfil << "set size ratio 1" << endl;
   gfil << "set notics" << endl;
   gfil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   if (!((xx[0][0]==bn.R[idx1][0])&&(xx[0][1]==bn.R[idx1][1])&&(xx[0][2]==bn.R[idx1][2]))) {
      int tmpaaa=idx1;
      idx1=idx2;
      idx2=tmpaaa;
   }
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set style fill solid 1.0 border lt -1" << endl;
   gfil << "set object 1 circle at 0.25*dimparam,0.25*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   gfil << "set object 2 circle at 0.75*dimparam,0.75*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   if ( options.findcps ) {
      gfil << "set style fill solid 1.0 border lt -1" << endl;
      int nobjs=3;
      for ( int i=0 ; i<cp.nACP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (0.25e0+0.5e0*cp.RACP[i][0]) << "*dimparam),("
            << (0.25e0+0.5e0*cp.RACP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'white'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblACP[i] << "' at (" << (0.25e0+0.5e0*cp.RACP[i][0]) << "*dimparam),("
                << (0.25e0+0.5e0*cp.RACP[i][1]) << "*dimparam) front center";
         gfil << " offset character " << (i==0? "-2" : "-1");
         gfil << ",0 font \",8\""<< endl;
      }
      nobjs+=cp.nACP;
      for ( int i=0 ; i<cp.nSCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (0.25e0+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
            << (0.25e0+0.5e0*cp.RSCP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'blue'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblSCP[i] << "' at (" << (0.25e0+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
                << (0.25e0+0.5e0*cp.RSCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
      nobjs+=cp.nSCP;
      for ( int i=0 ; i<cp.nRCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (0.25e0+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
            << (0.25e0+0.5e0*cp.RRCP[i][1]) << "*dimparam) front size (dimparam*0.008) "
            << "fc rgb 'green'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblRCP[i] << "' at (" << (0.25e0+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
                << (0.25e0+0.5e0*cp.RRCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
   }

   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 nohead lt 2 lc rgb 'black' lw 1 from " << (0.25*linelength) << 
           ",0 to " << (0.25*linelength) << ",dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 nohead lt 2 lc rgb 'black' lw 1 from " << (0.75*linelength) << 
           ",0 to " << (0.75*linelength) << ",dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from 0," << (0.25*linelength) << 
           " to dimparam," << (0.25*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from 0," << (0.75*linelength) << 
           " to dimparam," << (0.75*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
      << 0.25e0*(linelength) << "," << 0.25*(linelength) << " front offset character -2.2,-1" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << 0.75e0*(linelength) << "," << 0.75*(linelength) << " front offset character 0.5,0.8" << endl;
gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot namedatfile with image notitle";
   if (!options.showcont) {gfil << "#";}
   gfil << ",'contourtemp.dat' w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void generateVectorField(optFlags &options,char *argv[],const string &tsvname,bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,
      solreal **xx,int nptx,solreal minval2plot,\
      solreal maxval2plot,solreal maggradmin,solreal maggradmax,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2)
{
   char prop='D';
   if ( options.prop2plot ) {prop=argv[options.prop2plot][0];}
   string gnpname=tsvname;
   insertAtEndOfFileName(gnpname,string("-2DV"));
   replaceExtensionOfFileName(gnpname,string("gnp"));
   string pdfname=gnpname;
   replaceExtensionOfFileName(pdfname,string("pdf"));
   string epsname=gnpname;
   replaceExtensionOfFileName(epsname,string("eps"));
   //----------------------------------------------------
   ofstream gfil;
   gfil.open(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << tsvname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   gfil << "set isosample 300, 300" << endl;
   gfil << "set cntrparam cubicspline" << endl;
   gfil << "set zrange [minval2plot:maxval2plot]" << endl;
   gfil << "set cbrange [minval2plot:maxval2plot]" << endl;
   gfil << "unset contour" << endl;
   gfil << "set pm3d depthorder hidden3d 1" << endl;
   gfil << "set palette rgbformulae 33,13,10" << endl;
   gfil << "set style line 1 linecolor rgb \"#444444\"" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "unset key" << endl;
   gfil << "set contour base" << endl;
   int noconts=DEFAULTNUMBEROFCONTOURLINES;
   solreal dcont=(maxval2plot-minval2plot)/solreal(noconts);
   solreal contval=md1lmin+0.25e0*dcont;
   gfil << "set cntrparam cubicspline" << endl;
   gfil << "set cntrparam levels ";
   if ( !(options.setinccont) ) {
      gfil << "discrete " << contval;
      gfil << ", " << (md1dmax-0.5*dcont);
      contval=minval2plot+dcont;
      for (int i=1; i<noconts; i++) {gfil << ", " << contval; contval+=dcont;}
      gfil << endl;
   } else {
      gfil << "incremental " << argv[options.setinccont] << ","
           << argv[options.setinccont+1] << "," << argv[options.setinccont+2] << endl;
   }
   gfil << "unset surface" << endl;
   gfil << "set table 'contourtemp.dat'" << endl;
   gfil << "splot namedatfile" << endl;
   gfil << "unset table" << endl;
   gfil << "set xrange[0:dimparam]" << endl;
   gfil << "set yrange [minval2plot:maxval2plot]" << endl;
   gfil << "unset title" << endl;
   gfil << "set yrange[0:dimparam]" << endl;
   gfil << "set tmargin at screen 0.98\nset bmargin at screen 0.01\nset lmargin at screen 0.02" << endl;
   gfil << "set size ratio 1" << endl;
   gfil << "set notics" << endl;
   gfil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   if (!((xx[0][0]==bn.R[idx1][0])&&(xx[0][1]==bn.R[idx1][1])&&(xx[0][2]==bn.R[idx1][2]))) {
      int tmpaaa=idx1;
      idx1=idx2;
      idx2=tmpaaa;
   }
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set style fill solid 1.0 border lt -1" << endl;
   gfil << "set object 1 circle at 0.25*dimparam,0.25*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   gfil << "set object 2 circle at 0.75*dimparam,0.75*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   if ( options.findcps ) {
      gfil << "set style fill solid 1.0 border lt -1" << endl;
      int nobjs=3;
      for ( int i=0 ; i<cp.nACP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (0.25e0+0.5e0*cp.RACP[i][0]) << "*dimparam),("
            << (0.25e0+0.5e0*cp.RACP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'white'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblACP[i] << "' at (" << (0.25e0+0.5e0*cp.RACP[i][0]) << "*dimparam),("
                << (0.25e0+0.5e0*cp.RACP[i][1]) << "*dimparam) front center";
         gfil << " offset character " << (i==0? "-2" : "-1");
         gfil << ",0 font \",8\""<< endl;
      }
      nobjs+=cp.nACP;
      for ( int i=0 ; i<cp.nSCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (0.25e0+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
            << (0.25e0+0.5e0*cp.RSCP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'blue'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblSCP[i] << "' at (" << (0.25e0+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
                << (0.25e0+0.5e0*cp.RSCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
      nobjs+=cp.nSCP;
      for ( int i=0 ; i<cp.nRCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (0.25e0+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
            << (0.25e0+0.5e0*cp.RRCP[i][1]) << "*dimparam) front size (dimparam*0.008) "
            << "fc rgb 'green'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblRCP[i] << "' at (" << (0.25e0+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
                << (0.25e0+0.5e0*cp.RRCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
   }

   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 nohead lt 2 lc rgb 'black' lw 1 from " << (0.25*linelength) << 
           ",0 to " << (0.25*linelength) << ",dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 nohead lt 2 lc rgb 'black' lw 1 from " << (0.75*linelength) << 
           ",0 to " << (0.75*linelength) << ",dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from 0," << (0.25*linelength) << 
           " to dimparam," << (0.25*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from 0," << (0.75*linelength) << 
           " to dimparam," << (0.75*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
      << 0.25e0*(linelength) << "," << 0.25*(linelength) << " front offset character -2.2,-1" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << 0.75e0*(linelength) << "," << 0.75*(linelength) << " front offset character 0.5,0.8" << endl;
gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot namedatfile with image notitle";
   if (!options.showcont) {gfil << "#";}
   gfil << ",'contourtemp.dat' w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   //
   if ( prop=='G' ) {
#if DEBUG
      cout << "maggradmin: " << maggradmin << ", maggradmax: " << maggradmax << endl;
#endif
      minval2plot=round(90*sqrt(maggradmin))/100.0e0;
      if ( (10.0e0*maxval2plot)<(sqrt(maggradmax)) ) {
         maxval2plot*=10.0e0;
      } else {
         maxval2plot=round(100*sqrt(maggradmax))/100.0e0;
      }
      gfil << "minval2plot=" << minval2plot << endl;
      gfil << "maxval2plot=" << maxval2plot << endl;
      gfil << "set zrange [minval2plot:maxval2plot]" << endl;
      gfil << "set cbrange [minval2plot:maxval2plot]" << endl;
      gfil << "VMXS=dimparam/40.0 #Maximum lenght of the vectors" << endl;
      if ( options.findcps ) {
         int nobjs=3;
         for ( int i=0 ; i<cp.nACP ; i++ ) {
            gfil << "set object " << getStringFromInt(i+nobjs) 
                 << " circle at (" << (0.25e0+0.5e0*cp.RACP[i][0]) << ")*dimparam,("
                 << (0.25e0+0.5e0*cp.RACP[i][1]) << ")*dimparam front size dimparam*0.008 "
                 << "fc rgb 'white'" << endl;
         }
         nobjs+=cp.nACP;
         for ( int i=0 ; i<cp.nSCP ; i++ ) {
            gfil << "set object " << getStringFromInt(i+nobjs) 
               << " circle at (" << (0.25e0+0.5e0*cp.RSCP[i][0]) << ")*dimparam,("
               << (0.25e0+0.5e0*cp.RSCP[i][1]) << ")*dimparam front size dimparam*0.008 "
               << "fc rgb 'blue'" << endl;
         }
         nobjs+=cp.nSCP;
         for ( int i=0 ; i<cp.nRCP ; i++ ) {
            gfil << "set object " << getStringFromInt(i+nobjs) 
               << " circle at (" << (0.25e0+0.5e0*cp.RRCP[i][0]) << ")*dimparam,("
               << (0.25e0+0.5e0*cp.RRCP[i][1]) << ")*dimparam front size dimparam*0.008 "
               << "fc rgb 'green'" << endl;
            ++nobjs;
         }
      }
      gfil << "set output '" << epsname << "'" << endl;
      gfil << "plot namedatfile u 1:2:(sqrt($4*$4+$5*$5)>VMXS? VMXS*$4/sqrt($4*$4+$5*$5) : $4):(sqrt($4*$4+$5*$5)>VMXS ? VMXS*$5/sqrt($4*$4+$5*$5) : $5):(sqrt($4*$4+$5*$5)) "
         << "with vectors head size 0.1,20,60 filled lc palette";
      if (!options.showcont) {gfil << "#";}
      gfil << ",'contourtemp.dat' w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   }
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}

// */

