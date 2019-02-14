#ifndef _HELPERPLOTS_CPP_
#define _HELPERPLOTS_CPP_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "helperplots.h"
#include "../common/solstringtools.h"
#include "../common/solgnuplottools.h"
#include "../common/solfileutils.h"


void HelperPlot::generateMainDiagPlot(optFlags &options,const string &datname,\
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
   solreal at1pos,at2pos;
   if ( options.centredats ) {
      at1pos=-0.25*linelength;
      at2pos=0.25*linelength;
   } else {
      at1pos=0.25*linelength;
      at2pos=0.75*linelength;
   }
   solreal arroffset=0.1*linelength;
   solreal lbloffset=0.15*linelength;
   ofstream gfil(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << datname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   if ( options.centredats ) {
      gfil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      gfil << "set xrange[0:dimparam]" << endl;
   }
   gfil << "set yrange [minval2plot:maxval2plot]" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 from " << (at1pos-arroffset) << "," << (minval2plot+0.15*frange)
   << " to " << at1pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 from " << at2pos+arroffset << "," << (minval2plot+0.15*frange)
        << " to " << at2pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << at1pos << ",minval2plot "
        << "to " << at1pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << at2pos << ",minval2plot "
        << "to " << at2pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
   << (at1pos-lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character -2,0.00" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << (at2pos+lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character 0.3,0.00" << endl;
   gfil << "set title '" << getEnhancedEpsTitle(datname) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot '" << datname << "' w lines lw 2 lc 'black' notitle"  << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void HelperPlot::generateSecDiagPlot(optFlags &options,const string &datname,\
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
   solreal at1pos,at2pos;
   if ( options.centredats ) {
      at1pos=-0.25*linelength;
      at2pos=0.25*linelength;
   } else {
      at1pos=0.25*linelength;
      at2pos=0.75*linelength;
   }
   solreal arroffset=0.1*linelength;
   solreal lbloffset=0.15*linelength;
   ofstream gfil(gnpname.c_str(),ios::out);
   gfil << "namedatfile='" << datname << "'" << endl;
   gfil << "minval2plot=" << minval2plot << endl;
   gfil << "maxval2plot=" << maxval2plot << endl;
   gfil << "dimparam=" << linelength << endl;
   if ( options.centredats ) {
      gfil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      gfil << "set xrange[0:dimparam]" << endl;
   }
   gfil << "set yrange [minval2plot:maxval2plot]" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 from " << (at1pos-arroffset) << "," << (minval2plot+0.15*frange)
   << " to " << at1pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 from " << at2pos+arroffset << "," << (minval2plot+0.15*frange)
        << " to " << at2pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << at1pos << ",minval2plot "
        << "to " << at1pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << at2pos << ",minval2plot "
        << "to " << at2pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
   << (at1pos-lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character -2,0.00" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << (at2pos+lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character 0.3,0.00" << endl;
   gfil << "set title '" << getEnhancedEpsTitle(datname) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot '" << datname << "' w lines lw 2 lc 'black' notitle"  << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}

void HelperPlot::generate3DPlot(optFlags &options,const string &tsvname,\
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
   int theevery=2;
   if ((nptsinline>80&&nptsinline<=120)) {theevery=2;}
   if ((nptsinline>120&&nptsinline<=160)) {theevery=3;}
   if ((nptsinline>160&&nptsinline<=240)) {theevery=4;}
   if ((nptsinline>240&&nptsinline<=300)) {theevery=5;}
   if ( nptsinline>300 ) { theevery=int(nptsinline/60); }
   gfil << "every " << theevery << ":" << theevery << " ";
   gfil << "using 1:2:($3>maxval2plot? maxval2plot:($3<minval2plot ? minval2plot : $3)) "
   << "with pm3d ls 1 title '" << getEnhancedEpsTitle(tsvname) << "'" << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void HelperPlot::generateHeatMap(optFlags &options,char *argv[],const string &tsvname,bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,
      solreal **xx,int nptsinline,solreal minval2plot,\
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
   gfil << "unset title" << endl;
   if ( options.centredats ) {
      gfil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
      gfil << "set yrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      gfil << "set xrange[0:dimparam]" << endl;
      gfil << "set yrange[0:dimparam]" << endl;
   }
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
   solreal at1relpos,at2relpos,cpoffset;
   if ( options.centredats ) {
      at1relpos=-0.25;
      at2relpos=0.25;
      cpoffset=-0.25;
   } else {
      at1relpos=0.25;
      at2relpos=0.75;
      cpoffset=0.25;
   }
   gfil << "set object 1 circle at " << at1relpos << "*dimparam,"
        << at1relpos << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   gfil << "set object 2 circle at " << at2relpos << "*dimparam,"
        << at2relpos << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   if ( options.findcps ) {
      gfil << "set style fill solid 1.0 border lt -1" << endl;
      int nobjs=3;
      for ( int i=0 ; i<cp.nACP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'white'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblACP[i] << "' at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front center";
         gfil << " offset character " << (i==0? "-2" : "-1");
         gfil << ",0 font \",8\""<< endl;
      }
      nobjs+=cp.nACP;
      for ( int i=0 ; i<cp.nSCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'blue'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblSCP[i] << "' at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
      nobjs+=cp.nSCP;
      for ( int i=0 ; i<cp.nRCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front size (dimparam*0.008) "
            << "fc rgb 'green'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblRCP[i] << "' at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
   }
   solreal lowoffset=0.0,uppoffset=1.0;
   if ( options.centredats ) { lowoffset=-0.5; uppoffset=0.5; }
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 nohead lt 2 lc rgb 'black' lw 1 from " << (at1relpos*linelength) << 
           "," << lowoffset << "*dimparam to " << (at1relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 nohead lt 2 lc rgb 'black' lw 1 from " << (at2relpos*linelength) << 
           "," << lowoffset << "*dimparam to " << (at2relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at1relpos*linelength) << 
           " to " << uppoffset << "*dimparam," << (at1relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at2relpos*linelength) << 
           " to dimparam*" << uppoffset << "," << (at2relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
      << at1relpos*(linelength) << "," << at1relpos*(linelength) << " front offset character -2.2,-1" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << at2relpos*(linelength) << "," << at2relpos*(linelength) << " front offset character 0.5,0.8" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot namedatfile ";
   int theevery=2;
   if ((nptsinline>80&&nptsinline<=120)) {theevery=2;}
   if ((nptsinline>120&&nptsinline<=160)) {theevery=3;}
   if ((nptsinline>160&&nptsinline<=240)) {theevery=4;}
   if ((nptsinline>240&&nptsinline<=300)) {theevery=5;}
   if ( nptsinline>300 ) { theevery=int(nptsinline/60); }
   gfil << "every " << theevery << ":" << theevery << " ";
   gfil << " with image notitle";
   if (!options.showcont) {gfil << "#";}
   gfil << ",'contourtemp.dat' w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void HelperPlot::generateVectorField(optFlags &options,char *argv[],const string &tsvname,bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,
      solreal **xx,int nptsinline,solreal minval2plot,\
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
   gfil << "unset title" << endl;
   if ( options.centredats ) {
      gfil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
      gfil << "set yrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      gfil << "set xrange[0:dimparam]" << endl;
      gfil << "set yrange[0:dimparam]" << endl;
   }
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
   solreal at1relpos,at2relpos,cpoffset;
   if ( options.centredats ) {
      at1relpos=-0.25;
      at2relpos=0.25;
      cpoffset=-0.25;
   } else {
      at1relpos=0.25;
      at2relpos=0.75;
      cpoffset=0.25;
   }
   gfil << "set object 1 circle at " << at1relpos << "*dimparam," << at1relpos 
        << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   gfil << "set object 2 circle at " << at2relpos << "*dimparam," << at2relpos
        << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   if ( options.findcps ) {
      gfil << "set style fill solid 1.0 border lt -1" << endl;
      int nobjs=3;
      for ( int i=0 ; i<cp.nACP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'white'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblACP[i] << "' at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front center";
         gfil << " offset character " << (i==0? "-2" : "-1");
         gfil << ",0 font \",8\""<< endl;
      }
      nobjs+=cp.nACP;
      for ( int i=0 ; i<cp.nSCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'blue'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblSCP[i] << "' at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
      nobjs+=cp.nSCP;
      for ( int i=0 ; i<cp.nRCP ; i++ ) {
         gfil << "set object " << getStringFromInt(i+nobjs) 
            << " circle at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front size (dimparam*0.008) "
            << "fc rgb 'green'" << endl;
         gfil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblRCP[i] << "' at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front center";
         gfil << " offset character -1,0 font \",8\""<< endl;
      }
   }
   solreal lowoffset=0.0,uppoffset=1.0;
   if ( options.centredats ) { lowoffset=-0.5; uppoffset=0.5; }
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 1 nohead lt 2 lc rgb 'black' lw 1 from " << (at1relpos*linelength) << 
           "," << lowoffset << "*dimparam to " << (at1relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 2 nohead lt 2 lc rgb 'black' lw 1 from " << (at2relpos*linelength) << 
           "," << lowoffset << "*dimparam to " << (at2relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at1relpos*linelength) << 
           " to dimparam*" << uppoffset << "," << (at1relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at2relpos*linelength) << 
           " to dimparam*" << uppoffset << "," << (at2relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
      << (at1relpos*linelength) << "," << (at1relpos*linelength) << " front offset character -2.2,-1" << endl;
   if (!options.showatlbls) {gfil << "#";}
   gfil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << (at2relpos*linelength) << "," << (at2relpos*linelength) << " front offset character 0.5,0.8" << endl;
   gfil << "set output '" << epsname << "'" << endl;
   gfil << "plot namedatfile ";
   int theevery=2;
   if ((nptsinline>80&&nptsinline<=120)) {theevery=2;}
   if ((nptsinline>120&&nptsinline<=160)) {theevery=3;}
   if ( nptsinline>160 ) { theevery=int(nptsinline/40); }
   gfil << "every " << theevery << ":" << theevery << " ";
   gfil << "with image notitle";
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
                 << " circle at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << ")*dimparam,("
                 << (cpoffset+0.5e0*cp.RACP[i][1]) << ")*dimparam front size dimparam*0.008 "
                 << "fc rgb 'white'" << endl;
         }
         nobjs+=cp.nACP;
         for ( int i=0 ; i<cp.nSCP ; i++ ) {
            gfil << "set object " << getStringFromInt(i+nobjs) 
               << " circle at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << ")*dimparam,("
               << (cpoffset+0.5e0*cp.RSCP[i][1]) << ")*dimparam front size dimparam*0.008 "
               << "fc rgb 'blue'" << endl;
         }
         nobjs+=cp.nSCP;
         for ( int i=0 ; i<cp.nRCP ; i++ ) {
            gfil << "set object " << getStringFromInt(i+nobjs) 
               << " circle at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << ")*dimparam,("
               << (cpoffset+0.5e0*cp.RRCP[i][1]) << ")*dimparam front size dimparam*0.008 "
               << "fc rgb 'green'" << endl;
            ++nobjs;
         }
      }
      gfil << "set output '" << epsname << "'" << endl;
      gfil << "plot namedatfile ";
      gfil << "every " << theevery << ":" << theevery << " ";
      gfil << "u 1:2:(sqrt($4*$4+$5*$5)>VMXS? VMXS*$4/sqrt($4*$4+$5*$5) : $4):(sqrt($4*$4+$5*$5)>VMXS ? VMXS*$5/sqrt($4*$4+$5*$5) : $5):(sqrt($4*$4+$5*$5)) "
         << "with vectors head size 0.1,20,60 filled lc palette";
      if (!options.showcont) {gfil << "#";}
      gfil << ",'contourtemp.dat' w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   }
   gfil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}

#endif  /* _HELPERPLOTS_CPP_ */

