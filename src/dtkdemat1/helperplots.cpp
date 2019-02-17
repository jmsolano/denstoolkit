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
   string gnpname,epsname,pdfname;
   generateGNPEPSAndPDFNamesFromDATORTSV(datname,gnpname,epsname,pdfname);
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
   ofstream ofil(gnpname.c_str(),ios::out);
   addHeaderInfo2GNP(ofil,minval2plot,maxval2plot,linelength,datname,string("y"));
   if ( options.centredats ) {
      ofil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      ofil << "set xrange[0:dimparam]" << endl;
   }
   ofil << "set yrange [minval2plot:maxval2plot]" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 1 from " << (at1pos-arroffset) << "," << (minval2plot+0.15*frange)
   << " to " << at1pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 2 from " << at2pos+arroffset << "," << (minval2plot+0.15*frange)
        << " to " << at2pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << at1pos << ",minval2plot "
        << "to " << at1pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << at2pos << ",minval2plot "
        << "to " << at2pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
   << (at1pos-lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character -2,0.00" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << (at2pos+lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character 0.3,0.00" << endl;
   ofil << "set title '" << getEnhancedEpsTitle(datname) << "'" << endl;
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   ofil << "set output '" << epsname << "'" << endl;
   ofil << "plot '" << datname << "' w lines lw 2 lc 'black' notitle"  << endl;
   ofil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void HelperPlot::generateSecDiagPlot(optFlags &options,const string &datname,\
      bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
      solreal linelength,solreal frange)
{
   string gnpname,epsname,pdfname;
   generateGNPEPSAndPDFNamesFromDATORTSV(datname,gnpname,epsname,pdfname);
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
   ofstream ofil(gnpname.c_str(),ios::out);
   addHeaderInfo2GNP(ofil,minval2plot,maxval2plot,linelength,datname,string("y"));
   if ( options.centredats ) {
      ofil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      ofil << "set xrange[0:dimparam]" << endl;
   }
   ofil << "set yrange [minval2plot:maxval2plot]" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 1 from " << (at1pos-arroffset) << "," << (minval2plot+0.15*frange)
   << " to " << at1pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 2 from " << at2pos+arroffset << "," << (minval2plot+0.15*frange)
        << " to " << at2pos << "," << (minval2plot+0.10*frange) << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << at1pos << ",minval2plot "
        << "to " << at1pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << at2pos << ",minval2plot "
        << "to " << at2pos << ",maxval2plot" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
   << (at1pos-lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character -2,0.00" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << (at2pos+lbloffset) << "," << (minval2plot+frange*0.15e0) << " front offset character 0.3,0.00" << endl;
   ofil << "set title '" << getEnhancedEpsTitle(datname) << "'" << endl;
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   ofil << "set output '" << epsname << "'" << endl;
   ofil << "plot '" << datname << "' w lines lw 2 lc 'black' notitle"  << endl;
   ofil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}

void HelperPlot::generate3DPlot(optFlags &options,const string &tsvname,\
      solreal minval2plot,solreal maxval2plot,solreal linelength,int nptsinline)
{
   string gnpname=tsvname,epsname,pdfname;
   insertAtEndOfFileName(gnpname,string("-3D"));
   replaceExtensionOfFileName(gnpname,string("gnp"));
   generateEPSAndPDFNamesFromGNP(gnpname,epsname,pdfname);
   //----------------------------------------------------
   ofstream ofil(gnpname.c_str(),ios::out);
   addHeaderInfo2GNP(ofil,minval2plot,maxval2plot,linelength,tsvname,string("z"));
   ofil << "set isosample 300, 300" << endl;
   ofil << "set cntrparam cubicspline" << endl;
   ofil << "set zrange [minval2plot:maxval2plot]" << endl;
   ofil << "set cbrange [minval2plot:maxval2plot]" << endl;
   ofil << "unset contour" << endl;
   ofil << "set pm3d depthorder hidden3d 1" << endl;
   ofil << "set palette rgbformulae 33,13,10" << endl;
   ofil << "set style line 1 linecolor rgb \"#444444\"" << endl;
   ofil << "unset title" << endl;
   ofil << "set xyplane 0" << endl;
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   ofil << "#set output '|epstopdf --filter --outfile=" << pdfname << "'" << endl;
   ofil << "set output '" << epsname << "'" << endl;
   ofil << "splot namedatfile ";
   int theevery=2;
   if ((nptsinline>80&&nptsinline<=120)) {theevery=2;}
   if ((nptsinline>120&&nptsinline<=160)) {theevery=3;}
   if ((nptsinline>160&&nptsinline<=240)) {theevery=4;}
   if ((nptsinline>240&&nptsinline<=300)) {theevery=5;}
   if ( nptsinline>300 ) { theevery=int(nptsinline/60); }
   ofil << "every " << theevery << ":" << theevery << " ";
   ofil << "using 1:2:($3>maxval2plot? maxval2plot:($3<minval2plot ? minval2plot : $3)) "
   << "with pm3d ls 1 title '" << getEnhancedEpsTitle(tsvname) << "'" << endl;
   ofil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void HelperPlot::generateHeatMap(optFlags &options,char *argv[],const string &tsvname,bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,
      solreal **xx,int nptsinline,solreal minval2plot,\
      solreal maxval2plot,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2)
{
   string gnpname=tsvname,epsname,pdfname;
   insertAtEndOfFileName(gnpname,string("-2D"));
   replaceExtensionOfFileName(gnpname,string("gnp"));
   generateEPSAndPDFNamesFromGNP(gnpname,epsname,pdfname);
   //----------------------------------------------------
   ofstream ofil;
   ofil.open(gnpname.c_str(),ios::out);
   addHeaderInfo2GNP(ofil,minval2plot,maxval2plot,linelength,tsvname,string("z"));
   ofil << "set isosample 300, 300" << endl;
   ofil << "set cntrparam cubicspline" << endl;
   ofil << "set zrange [minval2plot:maxval2plot]" << endl;
   ofil << "set cbrange [minval2plot:maxval2plot]" << endl;
   ofil << "unset contour" << endl;
   ofil << "set pm3d depthorder hidden3d 1" << endl;
   ofil << "set palette rgbformulae 33,13,10" << endl;
   ofil << "set style line 1 linecolor rgb \"#444444\"" << endl;
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   ofil << "unset key" << endl;
   ofil << "set contour base" << endl;
   int noconts=DEFAULTNUMBEROFCONTOURLINES;
   solreal dcont=(maxval2plot-minval2plot)/solreal(noconts);
   solreal contval=md1lmin+0.25e0*dcont;
   ofil << "set cntrparam cubicspline" << endl;
   ofil << "set cntrparam levels ";
   if ( !(options.setinccont) ) {
      ofil << "discrete " << contval;
      ofil << ", " << (md1dmax-0.5*dcont);
      contval=minval2plot+dcont;
      for (int i=1; i<noconts; i++) {ofil << ", " << contval; contval+=dcont;}
      ofil << endl;
   } else {
      ofil << "incremental " << argv[options.setinccont] << ","
           << argv[options.setinccont+1] << "," << argv[options.setinccont+2] << endl;
   }
   ofil << "unset surface" << endl;
   string contourtempname=generateStrRandSeq(32);
   ofil << "contourtempname='" << contourtempname << "'" << endl;
   ofil << "set table contourtempname" << endl;
   ofil << "splot namedatfile" << endl;
   ofil << "unset table" << endl;
   ofil << "unset title" << endl;
   if ( options.centredats ) {
      ofil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
      ofil << "set yrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      ofil << "set xrange[0:dimparam]" << endl;
      ofil << "set yrange[0:dimparam]" << endl;
   }
   ofil << "set tmargin at screen 0.98\nset bmargin at screen 0.01\nset lmargin at screen 0.02" << endl;
   ofil << "set size ratio 1" << endl;
   ofil << "set notics" << endl;
   ofil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   if (!((xx[0][0]==bn.R[idx1][0])&&(xx[0][1]==bn.R[idx1][1])&&(xx[0][2]==bn.R[idx1][2]))) {
      int tmpaaa=idx1;
      idx1=idx2;
      idx2=tmpaaa;
   }
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set style fill solid 1.0 border lt -1" << endl;
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
   ofil << "set object 1 circle at " << at1relpos << "*dimparam,"
        << at1relpos << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   ofil << "set object 2 circle at " << at2relpos << "*dimparam,"
        << at2relpos << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   if ( options.findcps ) {
      ofil << "set style fill solid 1.0 border lt -1" << endl;
      int nobjs=3;
      for ( int i=0 ; i<cp.nACP ; i++ ) {
         ofil << "set object " << getStringFromInt(i+nobjs)
            << " circle at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'white'" << endl;
         ofil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblACP[i] << "' at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front center";
         ofil << " offset character " << (i==0? "-2" : "-1");
         ofil << ",0 font \",8\""<< endl;
      }
      nobjs+=cp.nACP;
      for ( int i=0 ; i<cp.nSCP ; i++ ) {
         ofil << "set object " << getStringFromInt(i+nobjs)
            << " circle at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'blue'" << endl;
         ofil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblSCP[i] << "' at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front center";
         ofil << " offset character -1,0 font \",8\""<< endl;
      }
      nobjs+=cp.nSCP;
      for ( int i=0 ; i<cp.nRCP ; i++ ) {
         ofil << "set object " << getStringFromInt(i+nobjs)
            << " circle at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front size (dimparam*0.008) "
            << "fc rgb 'green'" << endl;
         ofil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblRCP[i] << "' at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front center";
         ofil << " offset character -1,0 font \",8\""<< endl;
      }
   }
   solreal lowoffset=0.0,uppoffset=1.0;
   if ( options.centredats ) { lowoffset=-0.5; uppoffset=0.5; }
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 1 nohead lt 2 lc rgb 'black' lw 1 from " << (at1relpos*linelength) <<
           "," << lowoffset << "*dimparam to " << (at1relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 2 nohead lt 2 lc rgb 'black' lw 1 from " << (at2relpos*linelength) <<
           "," << lowoffset << "*dimparam to " << (at2relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at1relpos*linelength) <<
           " to " << uppoffset << "*dimparam," << (at1relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at2relpos*linelength) <<
           " to dimparam*" << uppoffset << "," << (at2relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
      << at1relpos*(linelength) << "," << at1relpos*(linelength) << " front offset character -2.2,-1" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << at2relpos*(linelength) << "," << at2relpos*(linelength) << " front offset character 0.5,0.8" << endl;
   ofil << "#set output '|epstopdf --filter --outfile=" << pdfname << "'" << endl;
   ofil << "set output '" << epsname << "'" << endl;
   ofil << "plot namedatfile ";
   int theevery=2;
   if ((nptsinline>80&&nptsinline<=120)) {theevery=2;}
   if ((nptsinline>120&&nptsinline<=160)) {theevery=3;}
   if ((nptsinline>160&&nptsinline<=240)) {theevery=4;}
   if ((nptsinline>240&&nptsinline<=300)) {theevery=5;}
   if ( nptsinline>300 ) { theevery=int(nptsinline/60); }
   ofil << "every " << theevery << ":" << theevery << " ";
   ofil << " with image notitle";
   if (!options.showcont) {ofil << "#";}
   ofil << ",\\\ncontourtempname w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   ofil << "\n\nsystem('rm -f '.contourtempname)" << endl;
   ofil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
   //contourtempname=string("rm -f ")+contourtempname;
   //system(contourtempname.c_str());
}
void HelperPlot::generateVectorField(optFlags &options,char *argv[],const string &tsvname,bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,
      solreal **xx,int nptsinline,solreal minval2plot,\
      solreal maxval2plot,solreal maggradmin,solreal maggradmax,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2)
{
   char prop='D';
   if ( options.prop2plot ) {prop=argv[options.prop2plot][0];}
   string gnpname=tsvname,epsname,pdfname;
   insertAtEndOfFileName(gnpname,string("-2DV"));
   replaceExtensionOfFileName(gnpname,string("gnp"));
   generateEPSAndPDFNamesFromGNP(gnpname,epsname,pdfname);
   //----------------------------------------------------
   ofstream ofil;
   ofil.open(gnpname.c_str(),ios::out);
   addHeaderInfo2GNP(ofil,minval2plot,maxval2plot,linelength,tsvname,string("colorbar"));
   ofil << "set isosample 300, 300" << endl;
   ofil << "set cntrparam cubicspline" << endl;
   ofil << "set zrange [minval2plot:maxval2plot]" << endl;
   ofil << "set cbrange [minval2plot:maxval2plot]" << endl;
   ofil << "unset contour" << endl;
   ofil << "set pm3d depthorder hidden3d 1" << endl;
   ofil << "set palette rgbformulae 33,13,10" << endl;
   ofil << "set style line 1 linecolor rgb \"#444444\"" << endl;
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   ofil << "unset key" << endl;
   ofil << "set contour base" << endl;
   int noconts=DEFAULTNUMBEROFCONTOURLINES;
   solreal dcont=(maxval2plot-minval2plot)/solreal(noconts);
   solreal contval=md1lmin+0.25e0*dcont;
   ofil << "set cntrparam cubicspline" << endl;
   ofil << "set cntrparam levels ";
   if ( !(options.setinccont) ) {
      ofil << "discrete " << contval;
      ofil << ", " << (md1dmax-0.5*dcont);
      contval=minval2plot+dcont;
      for (int i=1; i<noconts; i++) {ofil << ", " << contval; contval+=dcont;}
      ofil << endl;
   } else {
      ofil << "incremental " << argv[options.setinccont] << ","
           << argv[options.setinccont+1] << "," << argv[options.setinccont+2] << endl;
   }
   ofil << "unset surface" << endl;
   string contourtempname=generateStrRandSeq(32);
   ofil << "contourtempname='" << contourtempname << "'" << endl;
   ofil << "set table contourtempname" << endl;
   ofil << "splot namedatfile" << endl;
   ofil << "unset table" << endl;
   ofil << "unset title" << endl;
   if ( options.centredats ) {
      ofil << "set xrange[-0.5*dimparam:0.5*dimparam]" << endl;
      ofil << "set yrange[-0.5*dimparam:0.5*dimparam]" << endl;
   } else {
      ofil << "set xrange[0:dimparam]" << endl;
      ofil << "set yrange[0:dimparam]" << endl;
   }
   ofil << "set tmargin at screen 0.98\nset bmargin at screen 0.01\nset lmargin at screen 0.02" << endl;
   ofil << "set size ratio 1" << endl;
   ofil << "set notics" << endl;
   ofil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   if (!((xx[0][0]==bn.R[idx1][0])&&(xx[0][1]==bn.R[idx1][1])&&(xx[0][2]==bn.R[idx1][2]))) {
      int tmpaaa=idx1;
      idx1=idx2;
      idx2=tmpaaa;
   }
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set style fill solid 1.0 border lt -1" << endl;
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
   ofil << "set object 1 circle at " << at1relpos << "*dimparam," << at1relpos
        << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   ofil << "set object 2 circle at " << at2relpos << "*dimparam," << at2relpos
        << "*dimparam front size dimparam*0.012 fc rgb 'black'" << endl;
   if ( options.findcps ) {
      ofil << "set style fill solid 1.0 border lt -1" << endl;
      int nobjs=3;
      for ( int i=0 ; i<cp.nACP ; i++ ) {
         ofil << "set object " << getStringFromInt(i+nobjs)
            << " circle at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'white'" << endl;
         ofil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblACP[i] << "' at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RACP[i][1]) << "*dimparam) front center";
         ofil << " offset character " << (i==0? "-2" : "-1");
         ofil << ",0 font \",8\""<< endl;
      }
      nobjs+=cp.nACP;
      for ( int i=0 ; i<cp.nSCP ; i++ ) {
         ofil << "set object " << getStringFromInt(i+nobjs)
            << " circle at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front size dimparam*0.008 "
            << "fc rgb 'blue'" << endl;
         ofil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblSCP[i] << "' at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RSCP[i][1]) << "*dimparam) front center";
         ofil << " offset character -1,0 font \",8\""<< endl;
      }
      nobjs+=cp.nSCP;
      for ( int i=0 ; i<cp.nRCP ; i++ ) {
         ofil << "set object " << getStringFromInt(i+nobjs)
            << " circle at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
            << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front size (dimparam*0.008) "
            << "fc rgb 'green'" << endl;
         ofil << "set label " << getStringFromInt(i+nobjs) << " '"
                 << cp.lblRCP[i] << "' at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << "*dimparam),("
                << (cpoffset+0.5e0*cp.RRCP[i][1]) << "*dimparam) front center";
         ofil << " offset character -1,0 font \",8\""<< endl;
      }
   }
   solreal lowoffset=0.0,uppoffset=1.0;
   if ( options.centredats ) { lowoffset=-0.5; uppoffset=0.5; }
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 1 nohead lt 2 lc rgb 'black' lw 1 from " << (at1relpos*linelength) <<
           "," << lowoffset << "*dimparam to " << (at1relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 2 nohead lt 2 lc rgb 'black' lw 1 from " << (at2relpos*linelength) <<
           "," << lowoffset << "*dimparam to " << (at2relpos*linelength) << "," << uppoffset << "*dimparam front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 3 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at1relpos*linelength) <<
           " to dimparam*" << uppoffset << "," << (at1relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set arrow 4 nohead lt 2 lc rgb 'black' lw 1 from " << lowoffset << "*dimparam," << (at2relpos*linelength) <<
           " to dimparam*" << uppoffset << "," << (at2relpos*linelength) << " front" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 1 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx2 : idx1)]) << "' at "
      << (at1relpos*linelength) << "," << (at1relpos*linelength) << " front offset character -2.2,-1" << endl;
   if (!options.showatlbls) {ofil << "#";}
   ofil << "set label 2 '" << getEnhancedEpsAtLbl(bn.atLbl[(options.uponsl? idx1 : idx2)]) << "' at "
   << (at2relpos*linelength) << "," << (at2relpos*linelength) << " front offset character 0.5,0.8" << endl;
   ofil << "set output '" << epsname << "'" << endl;
   ofil << "plot namedatfile ";
   int theevery=2;
   if ((nptsinline>80&&nptsinline<=120)) {theevery=2;}
   if ((nptsinline>120&&nptsinline<=160)) {theevery=3;}
   if ( nptsinline>160 ) { theevery=int(nptsinline/40); }
   ofil << "every " << theevery << ":" << theevery << " ";
   ofil << "with image notitle";
   if (!options.showcont) {ofil << "#";}
   ofil << ",contourtempname w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
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
      ofil << "minval2plot=" << minval2plot << endl;
      ofil << "maxval2plot=" << maxval2plot << endl;
      ofil << "set zrange [minval2plot:maxval2plot]" << endl;
      ofil << "set cbrange [minval2plot:maxval2plot]" << endl;
      ofil << "VMXS=dimparam/40.0 #Maximum lenght of the vectors" << endl;
      if ( options.findcps ) {
         int nobjs=3;
         for ( int i=0 ; i<cp.nACP ; i++ ) {
            ofil << "set object " << getStringFromInt(i+nobjs)
                 << " circle at (" << (cpoffset+0.5e0*cp.RACP[i][0]) << ")*dimparam,("
                 << (cpoffset+0.5e0*cp.RACP[i][1]) << ")*dimparam front size dimparam*0.008 "
                 << "fc rgb 'white'" << endl;
         }
         nobjs+=cp.nACP;
         for ( int i=0 ; i<cp.nSCP ; i++ ) {
            ofil << "set object " << getStringFromInt(i+nobjs)
               << " circle at (" << (cpoffset+0.5e0*cp.RSCP[i][0]) << ")*dimparam,("
               << (cpoffset+0.5e0*cp.RSCP[i][1]) << ")*dimparam front size dimparam*0.008 "
               << "fc rgb 'blue'" << endl;
         }
         nobjs+=cp.nSCP;
         for ( int i=0 ; i<cp.nRCP ; i++ ) {
            ofil << "set object " << getStringFromInt(i+nobjs)
               << " circle at (" << (cpoffset+0.5e0*cp.RRCP[i][0]) << ")*dimparam,("
               << (cpoffset+0.5e0*cp.RRCP[i][1]) << ")*dimparam front size dimparam*0.008 "
               << "fc rgb 'green'" << endl;
            ++nobjs;
         }
      }
      ofil << "#set output '|epstopdf --filter --outfile=" << pdfname << "'" << endl;
      ofil << "set output '" << epsname << "'" << endl;
      ofil << "plot \\\nnamedatfile ";
      ofil << "every " << theevery << ":" << theevery << " ";
      ofil << "u 1:2:(sqrt($4*$4+$5*$5)>VMXS? VMXS*$4/sqrt($4*$4+$5*$5) : $4):(sqrt($4*$4+$5*$5)>VMXS ? VMXS*$5/sqrt($4*$4+$5*$5) : $5):(sqrt($4*$4+$5*$5)) "
         << "with vectors head size 0.1,20,60 filled lc palette";
      if (!options.showcont) {ofil << "#";}
      ofil << ",\\\ncontourtempname w l lt -1 lw 1 notitle #Activate this line to show contours" << endl;
   }
   ofil << "\n\nsystem('rm -f '.contourtempname)" << endl;
   ofil.close();
   bool rmgnp=!(options.kpgnp);
   renderGnpFile(gnpname,rmgnp);
   gnuplottools_eps2pdf(epsname);
}
void HelperPlot::addHeaderInfo2GNP(ofstream &ofil,solreal minval,solreal maxval,\
         solreal dimpar,string datortsv,string axis)
{
   ofil << "#\n#File generated by DensToolKit (dtkdemat1)\n#" << endl;
   ofil << "namedatfile='" << datortsv << "'" << endl;
   ofil << "#\n#Sets the minimum and maximum value to show in\n#the " << axis << " axis of the plot." << endl;
   ofil << "minval2plot=" << minval<< endl;
   ofil << "maxval2plot=" << maxval<< endl;
   ofil << "dimparam=" << dimpar << endl;
}
void HelperPlot::generateGNPEPSAndPDFNamesFromDATORTSV(const string &dnam,\
      string &gnam,string &enam,string &pnam)
{
   gnam=enam=pnam=dnam;
   replaceExtensionOfFileName(gnam,string("gnp"));
   generateEPSAndPDFNamesFromGNP(gnam,enam,pnam);
}
void HelperPlot::generateEPSAndPDFNamesFromGNP(const string &gnam,string &enam,string &pnam)
{
   enam=pnam=gnam;
   replaceExtensionOfFileName(pnam,string("pdf"));
   replaceExtensionOfFileName(enam,string("eps"));
}
#endif  /* _HELPERPLOTS_CPP_ */

