#ifndef _GNUPLOT_TOOLS_CPP_
#define _GNUPLOT_TOOLS_CPP_
#include "gnuplottools.h"
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <sstream>
using std::stringstream;
#include "fileutils.h"
#include "stringtools.h"
#include "screenutils.h"
#ifndef DEBUG
#define DEBUG 0
#endif
#ifndef _HAVE_EPSTOOL_
#define _HAVE_EPSTOOL_ 1
#endif
#ifndef _HAVE_EPSTOPDF_
#define _HAVE_EPSTOPDF_ 1
#endif

/* ************************************************************************** */
GnuplotTools::GnuplotTools() {
   init();
}
/* ************************************************************************** */
void GnuplotTools::init(void) {
   plotType=string("l");
   term="postscript eps enhanced color linewidth 2";
   termtype="pdf";
}
void GnuplotTools::GenerateGraphName(string dtnm) {
   graphnam=dtnm;
   if ( termtype==string("pdf") ) {
      FileUtils::ReplaceExtensionOfFileName(graphnam,string("pdf"));
      graphnam=string("|epstopdf --filter --outfile=")+graphnam;
   } else if ( termtype==string("png") ) {
      FileUtils::ReplaceExtensionOfFileName(graphnam,string("png"));
   }
}
void GnuplotTools::SetPNGTerminal(string ss) {
   term=string("png enhanced linewidth 2 ")+ss;
   termtype="png";
}
/* ************************************************************************** */
void GnuplotTools::MakeSimplePlotFromDat(string datnam,int cx,int cy) {
   MakeLogLogPlotFromDat(datnam,false,false,cx,cy);
}
/* ************************************************************************** */
void GnuplotTools::MakeSimplePlotFromMultiYDat(string datnam) {
   MakeLogLogPlotFromMultiYDat(datnam,false,false);
}
/* ************************************************************************** */
string GnuplotTools::GetGnuplotRGBColorString(const int r, const int g,\
      const int b) {
#if DEBUG
   ScreenUtils scrut;
   if (r>255||g>255||b>255||r<0||g<0||b<0) {
      scrut.DisplayErrorMessage("Bad rgb component...");
      DISPLAYDEBUGINFOFILELINE;
      return string("'#000000'");
   }
#endif
   stringstream oo;
   oo.fill('0'); oo.width(2); oo << std::hex << r;
   oo.fill('0'); oo.width(2); oo << std::hex << g;
   oo.fill('0'); oo.width(2); oo << std::hex << b;
   string res="'#";
   res+=(oo.str()+"'");
   return res;
}
/* ************************************************************************** */
/* **************************************************************************************** */
void GnuplotTools::MakeLogLogPlotFromMultiYDat(string datnam,bool setxlog,\
      bool setylog) {
   string gnpnam;
   FileUtils fut;
   StringTools st;
   fut.GenerateRandomTmpFileName(gnpnam,24);
   gnpnam.append(".gnp");
   GenerateGraphName(datnam);
   int nocols=fut.CountColumnsInFile(datnam);
   ofstream gfil;
   gfil.open(gnpnam.c_str(),ios::out);
   gfil << "set title '" << st.GetEnhancedEpsTitle(datnam) << "'" << endl;
   gfil << "set key horizontal bmargin center samplen 2 reverse" << endl;
   if (setxlog) {gfil << "set logscale x; set format x \"10^{%T}\"" << endl;}
   if (setylog) {gfil << "set logscale y; set format y \"10^{%T}\"" << endl;}
   gfil << "set terminal postscript eps enhanced color linewidth 2 "
        << "dashlength 4 fontscale 1.75\n";
   gfil << "set output '" << graphnam << "' " 
        << endl;
   gfil << "plot '" << datnam << "' using 1:2 w " << plotType 
        << " lw 2 title '2'";
   for (int i=3; i<=nocols; i++) {
      gfil << "\\\n  , '" << datnam << "' using 1:" << i << " w " << plotType 
           << " lw 2 title '" << i << "'";
   }
   gfil << endl;
   gfil.close();
   string l;
   l="gnuplot "+gnpnam;
   system(l.c_str());
   l="rm -f "+gnpnam;
   system(l.c_str());
   return;
}
/* ************************************************************************** */
void GnuplotTools::MakeLinLogPlotFromMultiYDat(string datnam) {
   MakeLogLogPlotFromMultiYDat(datnam,true,false);
}
/* ************************************************************************** */
void GnuplotTools::MakeLinLogPlotFromDat(string datnam,int cx,int cy) {
   MakeLogLogPlotFromDat(datnam,true,false,cx,cy);
}
/* ************************************************************************** */
void GnuplotTools::MakeLogLinPlotFromDat(string datnam,int cx,int cy) {
   MakeLogLogPlotFromDat(datnam,false,true,cx,cy);
}
/* ************************************************************************** */
void GnuplotTools::MakeLogLogPlotFromDat(string datnam,bool setxlog,\
      bool setylog,int cx,int cy) {
   string gnpnam;
   FileUtils fut;
   StringTools st;
   fut.GenerateRandomTmpFileName(gnpnam,24);
   gnpnam.append(".gnp");
   GenerateGraphName(datnam);
   ofstream gfil;
   gfil.open(gnpnam.c_str(),ios::out);
   if (setxlog) {gfil << "set logscale x; set format x \"10^{%T}\"" << endl;}
   if (setylog) {gfil << "set logscale y; set format y \"10^{%T}\"" << endl;}
   gfil << "set title '" << st.GetEnhancedEpsTitle(datnam) << "'" << endl;
   gfil << "set terminal " << term
        << " dashlength 4 fontscale 1.75\n";
   gfil << "set output '" << graphnam << "' "
        << endl;
   gfil << "plot '" << datnam << "' using " << cx << ":" << cy << " w "
        << plotType << " lw 2 notitle" << endl;
   gfil.close();
   string l;
   l="gnuplot "+gnpnam;
   system(l.c_str());
   l="rm -f "+gnpnam;
   system(l.c_str());
   return;
}
/* ************************************************************************** */
void GnuplotTools::MakeLogLogPlotFromDat(string datnam,int cx,int cy) {
   MakeLogLogPlotFromDat(datnam,true,true,cx,cy);
}
/* ************************************************************************** */
void GnuplotTools::MakeLogLogPlotFromMultiYDat(string datnam) {
   MakeLogLogPlotFromMultiYDat(datnam,true,true);
}
/* ************************************************************************** */
void GnuplotTools::MakeLogLinPlotFromMultiYDat(string datnam) {
   MakeLogLogPlotFromMultiYDat(datnam,false,true);
}
/* ************************************************************************** */
void GnuplotTools::MakeSimple3DPlotFromTSV(string tsvnam) {
   string gnpnam;
   FileUtils fut;
   StringTools st;
   fut.GenerateRandomTmpFileName(gnpnam,24);
   gnpnam.append(".gnp");
   GenerateGraphName(tsvnam);
   FileUtils::InsertAtEndOfFileName(graphnam,string("3D"));
   ofstream gfil;
   gfil.open(gnpnam.c_str(),ios::out);
   gfil << "set title '" << st.GetEnhancedEpsTitle(tsvnam) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color linewidth 2 "
        << "dashlength 4 fontscale 1.75\n";
   gfil << "set output '" << graphnam << "' "
        << endl;
   gfil << "set palette rgbformulae 33,13,10" << endl;
   gfil << "set linetype 1 lc rgb '#999999'" << endl;
   gfil << "set pm3d depthorder hidden3d 1" << endl;
   gfil << "set xyplane 0" << endl;
   gfil << "splot '" << tsvnam << "' w pm3d notitle" << endl;
   gfil.close();
   string l;
   l="gnuplot "+gnpnam;
   system(l.c_str());
   l="rm -f "+gnpnam;
   system(l.c_str());
   return;
}
/* ************************************************************************** */
void GnuplotTools::MakeSimpleHeatMapFromTSV(string tsvnam) {
   string gnpnam;
   FileUtils fut;
   StringTools st;
   fut.GenerateRandomTmpFileName(gnpnam,24);
   gnpnam.append(".gnp");
   GenerateGraphName(tsvnam);
   FileUtils::InsertAtEndOfFileName(graphnam,string("HM"));
   ofstream gfil;
   gfil.open(gnpnam.c_str(),ios::out);
   gfil << "set title '" << st.GetEnhancedEpsTitle(tsvnam) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color linewidth 2 "
        << "dashlength 4 fontscale 1.75\n";
   gfil << "set output '" << graphnam << "' "
        << endl;
   gfil << "set palette rgbformulae 33,13,10" << endl;
   gfil << "plot '" << tsvnam << "' w image notitle" << endl;
   gfil.close();
   string l;
   l="gnuplot "+gnpnam;
   system(l.c_str());
   l="rm -f "+gnpnam;
   system(l.c_str());
   return;
}
void GnuplotTools::eps2pdf(const string &epsname,bool quiet) {
   cout << "Trimming eps (epstool)..." << endl;
   string tmpepsname=StringTools::GenerateStrRandSeq(16);
   string cmdline;
#if ( _HAVE_EPSTOOL_)
   tmpepsname+=".eps";
   cmdline="epstool --copy -b ";
   cmdline+=epsname;
   cmdline+=" ";
   cmdline+=tmpepsname;
   if ( quiet ) { cmdline+=" > /dev/null 2>&1"; }
   system(cmdline.c_str());
#endif
   //_HAVE_EPSTOOL_
#if (_HAVE_EPSTOPDF_)
   cout << "Converting eps to pdf (epstopdf)..." << endl;
   string pdfname=epsname;
   FileUtils::ReplaceExtensionOfFileName(pdfname,"pdf");
   cmdline="epstopdf --outfile=";
   cmdline+=pdfname;
   cmdline+=" ";
   cmdline+=tmpepsname;
   if ( quiet ) { cmdline+=" > /dev/null 2>&1"; }
   system(cmdline.c_str());
   cmdline="rm "+tmpepsname;
   system(cmdline.c_str());
#endif
   // _HAVE_EPSTOPDF_
}
void GnuplotTools::RenderGnpFile(string gnpname,bool rmGnpFile) {
#if (defined(__CYGWIN__))
   string cmdline="wgnuplot ";
#else
   string cmdline="gnuplot ";
#endif
   cmdline+=gnpname;
   cout << "Calling gnuplot; cmd: '" << cmdline << "'" << endl;
   system(cmdline.c_str());
   if ( rmGnpFile ) {
      cmdline="rm "+gnpname;
      system(cmdline.c_str());
   }
}
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */


#endif /* defined(_GNUPLOT_TOOLS_CPP_) */
