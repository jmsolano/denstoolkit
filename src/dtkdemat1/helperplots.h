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

#ifndef _HELPERPLOTS_H_
#define _HELPERPLOTS_H_
#include "optflags.h"
#include <fstream>
using std::ofstream;
#include "../common/bondnetwork.h"
#include "../common/demat1critptnetworksl.h"


/* ************************************************************************** */
class HelperPlot {
/* ************************************************************************** */
public:
   static void generateMainDiagPlot(optFlags &options,const string &datname,\
         bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
         solreal linelength,solreal frange);
   static void generateSecDiagPlot(optFlags &options,const string &datname,\
         bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
         solreal linelength,solreal frange);
   static void generate3DPlot(optFlags &options,const string &tsvname,solreal minval2plot,\
         solreal maxval2plot,solreal linelength,int nptsinline);
   static void generateHeatMap(optFlags &options,char *argv[],const string &tsvname,\
         bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),int nptsinline,solreal minval2plot,\
         solreal maxval2plot,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
   static void generateVectorField(optFlags &options,char *argv[],const string &tsvname,\
         bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),int nptsinline,solreal minval2plot,\
         solreal maxval2plot,solreal maggradmin,solreal maggradmax,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
   static void addHeaderInfo2GNP(ofstream &ofil,solreal minval,solreal maxval,\
         solreal dimpar,string datortsv,string axis="y");
   static void generateGNPEPSAndPDFNamesFromDATORTSV(const string &dnam,\
         string &gnam,string &enam,string &pnam);
   static void generateEPSAndPDFNamesFromGNP(const string &gnam,\
         string &enam,string &pnam);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERPLOTS_H_ */

