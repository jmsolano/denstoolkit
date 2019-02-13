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
      solreal **xx,solreal minval2plot,\
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
      solreal **xx,solreal minval2plot,\
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

#endif  /* _HELPERPLOTS_CPP_ */

