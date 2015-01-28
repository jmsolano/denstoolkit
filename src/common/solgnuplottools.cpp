

#ifndef _SOL_GNUPLOT_TOOLS_CPP_
#define _SOL_GNUPLOT_TOOLS_CPP_

#include "solgnuplottools.h"
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <sstream>
using std::stringstream;
#include "solstringtools.h"
#include "solfileutils.h"
#include "solscrutils.h"

/* **************************************************************************************** */
void mkSimplePlotFromDat(string datnam,int cx,int cy)
{
   mkLogLogPlotFromDat(datnam,false,false,cx,cy);
}
/* **************************************************************************************** */
void mkSimplePlotFromMultiYDat(string datnam)
{
   mkLogLogPlotFromMultiYDat(datnam,false,false);
}
/* **************************************************************************************** */
string getGnuplotRGBColorString(const int r, const int g,const int b)
{
#if DEBUG
   if (r>255||g>255||b>255||r<0||g<0||b<0) {
      displayErrorMessage("Bad rgb component...");
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
/* **************************************************************************************** */
void mkLogLogPlotFromMultiYDat(string datnam,bool setxlog,bool setylog)
{
   string gnpnam,pdfnam;
   genRandomTmpFileName(gnpnam,24);
   gnpnam.append(".gnp");
   pdfnam=datnam.substr(0,datnam.length()-3);
   pdfnam.append("pdf");
   int nocols=countColumnsInFile(datnam);
   ofstream gfil;
   gfil.open(gnpnam.c_str(),ios::out);
   gfil << "set title '" << getEnhancedEpsTitle(datnam) << "'" << endl;
   gfil << "set key horizontal bmargin center samplen 2 reverse" << endl;
   if (setxlog) {gfil << "set logscale x" << endl;}
   if (setylog) {gfil << "set logscale y" << endl;}
   gfil << "set terminal postscript eps enhanced color linewidth 2 dashlength 4 fontscale 1.75\n";
   gfil << "set output '|epstopdf --filter --outfile=" << pdfnam << "' " << endl;
   gfil << "plot '" << datnam << "' using 1:2 w l lw 2 title '2'";
   for (int i=3; i<=nocols; i++) {
      gfil << "\\\n  , '" << datnam << "' using 1:" << i << " w l lw 2 title '" << i << "'";
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
/* **************************************************************************************** */
void mkLinLogPlotFromMultiYDat(string datnam)
{
   mkLogLogPlotFromMultiYDat(datnam,true,false);
}
/* **************************************************************************************** */
void mkLinLogPlotFromDat(string datnam,int cx,int cy)
{
   mkLogLogPlotFromDat(datnam,true,false,cx,cy);
}
/* **************************************************************************************** */
void mkLogLinPlotFromDat(string datnam,int cx,int cy)
{
   mkLogLogPlotFromDat(datnam,false,true,cx,cy);
}
/* **************************************************************************************** */
void mkLogLogPlotFromDat(string datnam,bool setxlog,bool setylog,int cx,int cy)
{
   string gnpnam,pdfnam;
   genRandomTmpFileName(gnpnam,24);
   gnpnam.append(".gnp");
   pdfnam=datnam.substr(0,datnam.length()-3);
   pdfnam.append("pdf");
   ofstream gfil;
   gfil.open(gnpnam.c_str(),ios::out);
   if (setxlog) {gfil << "set logscale x" << endl;}
   if (setylog) {gfil << "set logscale y" << endl;}
   gfil << "set title '" << getEnhancedEpsTitle(datnam) << "'" << endl;
   gfil << "set terminal postscript eps enhanced color linewidth 2 dashlength 4 fontscale 1.75\n";
   gfil << "set output '|epstopdf --filter --outfile=" << pdfnam << "' " << endl;
   gfil << "plot '" << datnam << "' using " << cx << ":" << cy << " w l lw 2 notitle" << endl;
   gfil.close();
   string l;
   l="gnuplot "+gnpnam;
   system(l.c_str());
   l="rm -f "+gnpnam;
   system(l.c_str());
   return;
}
/* **************************************************************************************** */
void mkLogLogPlotFromDat(string datnam,int cx,int cy)
{
   mkLogLogPlotFromDat(datnam,true,true,cx,cy);
}
/* **************************************************************************************** */
void mkLogLogPlotFromMultiYDat(string datnam)
{
   mkLogLogPlotFromMultiYDat(datnam,true,true);
}
/* **************************************************************************************** */
void mkLogLinPlotFromMultiYDat(string datnam)
{
   mkLogLogPlotFromMultiYDat(datnam,false,true);
}
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */
/* **************************************************************************************** */

#endif /* defined(_SOL_GNUPLOT_TOOLS_CPP_) */
