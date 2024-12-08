/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
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
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "helpersplot.h"
#include "gnuplottools.h"
#include "stringtools.h"
#include "fldtypesdef.h"
#include "fileutils.h"


void HelpersPlot::MakeGnuplotFile(OptionFlags &opts,string &gnpnam,string &datnam,char p2p,\
      BondNetWork &bn,double lenline,double minval,double maxval,\
      int at1,int at2,double pbcp,double** (&rbgp)) {
   ofstream ofil;
   string line;
   ofil.open(gnpnam.c_str(),std::ios::out);
   
   if (maxval>1000.e0) {maxval=1000.0e0;}
   
   ofil << "set datafile missing \"inf\"" << endl;
   ofil << "namedatfile='" << datnam << "'" << endl;
   ofil << "xmin2plot=0.0" << endl;
   ofil << "xmax2plot=" << lenline << endl;
   ofil << "ymin2plot=" << (0.0e0<minval ? 0.0e0 : minval) << endl;
   ofil << "ymax2plot=" << maxval*1.05e0 << endl;
   ofil << "set xrange [xmin2plot:xmax2plot]" << endl;
   ofil << "set yrange [ymin2plot:ymax2plot]" << endl;
   double tmpdist=0.0e0;
   for ( int i=0 ; i<3 ; i++ ) {tmpdist+=((bn.R[at1][i]-rbgp[0][i])*(bn.R[at1][i]-rbgp[0][i]));}
   tmpdist=sqrt(tmpdist);
   int tmpa1=at1,tmpa2=at2;
   if ( tmpdist>0.1 ) {tmpa1=at2; tmpa2=at1;}
   if (!opts.showatlbls) {ofil << "#";}
   ofil << "set label 1 '" << StringTools::GetEnhancedEpsAtLbl(bn.atLbl[tmpa1]) << "' at "
   << 0.05e0*(lenline) << "," << 0.1e0*(maxval) << " front offset character 0,0.75" << endl;
   if (!opts.showatlbls) {ofil << "#";}
   ofil << "set label 2 '" << StringTools::GetEnhancedEpsAtLbl(bn.atLbl[tmpa2]) << "' at "
   << 0.95e0*(lenline) << "," << 0.1e0*(maxval) << " front offset character -2,0.75" << endl;
   if (!opts.showatlbls) {ofil << "#";}
   ofil << "set arrow 1 from " << 0.05e0*(lenline) << "," << 0.1e0*(maxval) << " to "
        << "0.0,ymin2plot" << endl;
   if (!opts.showatlbls) {ofil << "#";}
   ofil << "set arrow 2 from " << 0.95e0*(lenline) << "," << 0.1e0*(maxval) << " to "
        << lenline << ",ymin2plot" << endl;
   if (opts.uponbp) {
      if (!opts.showatlbls) {ofil << "#";}
      ofil << "set label 3 'BCP' at "
      << pbcp << "," << 0.05e0*(maxval) << " rotate by 90 front" << endl;
      if (!opts.showatlbls) {ofil << "#";}
      ofil << "set xtics add ('' " << pbcp << ")" << endl;
   }
   ofil << "set title '" << StringTools::GetEnhancedEpsTitle(datnam) << "'" << endl;
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   string epsnam=gnpnam.substr(0,(gnpnam.length()-3));
   string pdfnam=epsnam;
   epsnam.append("eps");
   pdfnam.append("pdf");
#if _HAVE_EPSTOPDF_
   ofil << "set output '|epstopdf --filter --outfile=" << pdfnam << "'" << endl;
#else
   ofil << "set output '" << epsnam << "'" << endl;
#endif
   ofil << "plot namedatfile using 1:5 with lines lw 2 title '"
        << GnuplotFieldTitle(p2p) << "'" << endl;
   ofil << "#" << endl;
   FileUtils::WriteScrCharLine(ofil,'#');
   ofil << "#                 END OF GNUPLOT COMMANDS" << endl;
   FileUtils::WriteScrCharLine(ofil,'#');
   ofil << "#If you want to reconstruct the plots using this file, type:" << endl
   << "#gnuplot " << gnpnam << endl;
#if !(_HAVE_EPSTOPDF_)
   ofil << "#epstopdf --outfile=" << pdfnam << " " << epsnam << endl;
#endif
   ofil.close();
   if ( opts.mkplt ) {
      GnuplotTools::RenderGnpFile(gnpnam,(!opts.kpgnp));
   }
}
double HelpersPlot::EvalFieldProperty(char prop,double (&x)[3],GaussWaveFunction &wf) {
   double res;
   switch (prop) {
      case 'd':
         res=wf.EvalDensity(x[0],x[1],x[2]);
         break;
      case 'g':
         res=wf.EvalMagGradRho(x[0],x[1],x[2]);
         break;
      case 'l':
         res=wf.EvalLapRho(x[0],x[1],x[2]);
         break;
      case 'e':
         res=wf.EvalEllipticity(x[0],x[1],x[2]);
         break;
      case 'E':
         res=wf.EvalELF(x[0],x[1],x[2]);
         break;
      case 'L':
         res=wf.EvalLOL(x[0],x[1],x[2]);
         break;
      case 'M':
         res=wf.EvalMagGradLOL(x[0],x[1],x[2]);
         break;
      case 'P' :
         res=wf.EvalMagLED(x[0],x[1],x[2]);
         break;
      case 'r' :
         res=wf.EvalRoSE(x[0],x[1],x[2]);
         break;
      case 's' :
         res=wf.EvalReducedDensityGradient(x[0],x[1],x[2]);
         break;
      case 'S':
         res=wf.EvalShannonEntropy(x[0],x[1],x[2]);
         break;
      case 'G':
         res=wf.EvalKineticEnergyG(x[0],x[1],x[2]);
         break;
      case 'K':
         res=wf.EvalKineticEnergyK(x[0],x[1],x[2]);
         break;
      case 'u' :
         res=wf.EvalCustomScalarField(x[0],x[1],x[2]);
         break;
      case 'V':
         res=wf.EvalMolElecPot(x[0],x[1],x[2]);
         break;
      case 'v' :
         res=wf.EvalVirialPotentialEnergyDensity(x[0],x[1],x[2]);
         break;
      case 'D' :
         res=wf.EvalDORI(x[0],x[1],x[2]);
         break;
      case 'b' :
         res=wf.EvalSpinDensity(x[0],x[1],x[2]);
         break;
      case 'q' :
         res=wf.EvalOneElecDisequilibrium(x[0],x[1],x[2]);
         break;
      default:
#if DEBUG
         cout << "Unknown field type! Returning zero..." << endl;
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
#endif
         res=0.0e0;
         break;
   }
   return res;
}

