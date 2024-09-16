/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.1
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
