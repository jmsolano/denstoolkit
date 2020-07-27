#ifndef _GNUPLOT_TOOLS_H_
#define _GNUPLOT_TOOLS_H_

#include <iostream>
#include <string>
using std::string;
#include <fstream>
using std::ofstream;

/* ************************************************************************** */
class GnuplotTools {
/* ************************************************************************** */
public:
   GnuplotTools();
   void init();
   void SetPNGTerminal(string ss="");
   /* ************************************************************************** */
   void MakeSimplePlotFromDat(string datnam,int cx=1,int cy=2);
   void MakeSimplePlotFromMultiYDat(string datnam);
   void MakeSimple3DPlotFromTSV(string tsvnam);
   void MakeSimpleHeatMapFromTSV(string tsvnam);
   string GetGnuplotRGBColorString(const int r, const int g,const int b);
   /* ************************************************************************** */
   void MakeLogLinPlotFromMultiYDat(string datnam);
   void MakeLinLogPlotFromMultiYDat(string datnam);
   void MakeLinLogPlotFromDat(string datnam,int cx=1,int cy=2);
   void MakeLogLinPlotFromDat(string datnam,int cx=1,int cy=2);
   void MakeLogLogPlotFromDat(string datnam,bool setxlog,bool setylog,int cx=1,int cy=2);
   void MakeLogLogPlotFromMultiYDat(string datnam,bool setxlog,bool setylog);
   void MakeLogLogPlotFromDat(string datnam,int cx=1,int cy=2);
   void MakeLogLogPlotFromMultiYDat(string datnam);
   void SetPlotType(string pt) {plotType=pt;}
   void SetPlotType(const char *pt) {SetPlotType(string(pt));}
   static void eps2pdf(const string &epsname,bool quiet=true);
   static void RenderGnpFile(string gnpname,bool rmGnpFile=true);
   static void AddCommandsToRemoveTemporaryFileFromGnuplotScript(ofstream &ofil,const string &f2rm);
   /* ************************************************************************** */
protected:
   string plotType;
   string term,graphnam;
   string termtype;
   void GenerateGraphName(string dtnm);
   /* ************************************************************************** */
};
/* ************************************************************************** */

#endif /* defined(_GNUPLOT_TOOLS_H_) */
