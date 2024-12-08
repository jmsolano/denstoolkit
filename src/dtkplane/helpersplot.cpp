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
#include "fileutils.h"
#include "gnuplottools.h"

void HelpersPlot::makeGnuplotFile(OptionFlags &opts, string &gnpnam,string &tsvnam,char p2p,
      double dimparam,BondNetWork &bn,int a1,int a2,int a3,WaveFunctionGrid2D &grd) {
   ofstream gfil;
   gfil.open(gnpnam.c_str());
   
   /* Choosing the label (legend) for the plot  and the zrange for the plot */
   double minzrange,maxzrange;
   string plbl=GetFieldTypeKeyShort(p2p);;
   switch (p2p) {
      case 'b' :
         minzrange=0.0e+00;
         maxzrange=1.0e-01;
         break;
      case 'd':
      case 'q':
         minzrange=0.0e0;
         maxzrange=0.6e0;
         break;
      case 'g':
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'l':
         minzrange=-2.0e0;
         maxzrange=2.0e0;
         break;
      case 'e':
      case 'E':
      case 'D' :
         minzrange=0.0e0;
         maxzrange=1.0e0;
         break;
      case 'S':
         minzrange=-10.0e0;
         maxzrange=10.0e0;
         break;
      case 'L':
         minzrange=0.0e0;
         maxzrange=1.0e0;
         break;
      case 'M':
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'N':
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'G':
      case 'K':
         minzrange=0.0e0;
         maxzrange=10.0e0;
         break;
      case 'p' :
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'P' :
         minzrange=0.0e0;
         maxzrange=2.0e0;
         break;
      case 'r' :
         minzrange=-1.0e0;
         maxzrange=1.0e0;
         break;
      case 's' :
         minzrange=0.0e0;
         maxzrange=1.0e0;
         break;
      case 'u' :
         minzrange=-1.0e0;
         maxzrange=1.0e0;
         break;
      case 'U' :
         minzrange=-1.0e0;
         maxzrange=1.0e0;
         break;
      case 'V':
         minzrange=-0.6e0;
         maxzrange=0.6e0;
         break;
      case 'v' :
         minzrange=-1.0e0;
         maxzrange=0.0e0;
         break;
      default:
         ScreenUtils::SetScrRedBoldFont();
         cout << "Error: The property \"" << p2p << "\" does not exist!" << endl;
         cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
         ScreenUtils::SetScrNormalFont();
         exit(1);
         break;
   }
   
   gfil << "reset" << endl;
   
   /* In this part the name is scanned for possible occurrings of the character '_'.
    For a proper display in the eps file, it has to be changed to "\_" */
   string line="";
   for (size_t i=0; i<tsvnam.length(); i++) {
      if (tsvnam[i]=='_') {line+='\\';}
      line+=tsvnam[i];
   }
   gfil << "set title '" << line << "'" << endl;
   
   /* Adding an underscore and braces to the atome labels in order to make the numbers
    subindices */
   
   gfil << "set isosample 300, 300" << endl;
   gfil << "set zrange [" << minzrange << ":" << maxzrange << "]" << endl;
   gfil << "set contour base" << endl;
   gfil << "set cntrparam bspline" << endl;
   //gfil << "set cntrparam levels discrete 0.001, 0.005, 0.02, 0.05, 0.1, 0.2, 0.3" <<endl;
   gfil << "set cntrparam levels incremental " << minzrange << ", " 
        << (0.10*(maxzrange-minzrange)) << ", " << maxzrange << endl;
   gfil << "unset surface" << endl;

   string cntname="contourtemp.dat";
#if ( defined(__APPLE__)||defined(__linux__) )
   FileUtils::GenerateRandomTmpFileName(cntname, 20);
#endif
   gfil << "set table '" << cntname << "'" << endl;
   gfil << "splot '" << tsvnam << "'";
   if (p2p=='N' || p2p=='p' || p2p=='U') {gfil << " using 1:2:(sqrt($3*$3+$4*$4))";}
   gfil << endl;
   gfil << "unset table" << endl;
   gfil << "reset" << endl;
   gfil << "unset key" << endl;
   gfil << "set size ratio 1" << endl;
   gfil << "set tmargin at screen 0.99\nset bmargin at screen 0.01\nset lmargin at screen 0.01" << endl;
   gfil << "dimparam=" << dimparam
        << " #Decrease this number to zoom in the plot" << endl;
   if (p2p=='N'||p2p=='p'||p2p=='U') {
      gfil << "VMXS=dimparam/40.0 #Maximum lenght of the vectors" << endl;
   }
   gfil << "set notics" << endl;
   gfil << "set xrange[-dimparam:dimparam]\nset yrange[-dimparam:dimparam]" << endl;
   gfil << "minvalfield=" << minzrange << endl << "maxvalfield=" << maxzrange << endl;
   gfil << "set cbrange [minvalfield:maxvalfield]" << endl;
   gfil << "set style data lines" << endl;
   gfil << "set palette rgbformulae 33,13,10" << endl;
   
   /* Here the atoms' labels are set and written to the gnuplot input file. */
   
   FileUtils::WriteScrCharLine(gfil,'#');
   gfil << "# Here are the labes of the atoms" << endl;
   FileUtils::WriteScrCharLine(gfil,'#');
   
   string lbl;
   int at[3];
   at[0]=a1; at[1]=a2; at[2]=a3;
   double xproj,yproj;
   size_t pos;
   std::ostringstream numst;
   bool IDefPlane;
   for (int i=0; i<bn.nNuc; i++) {
      xproj=0.0e0;
      yproj=0.0e0;
      IDefPlane=false;
      for (int j=0; j<3; j++) {
         xproj+=(bn.R[i][j]-grd.orig[j])*grd.dircos1[j];
         yproj+=(bn.R[i][j]-grd.orig[j])*grd.dircos2[j];
         if ((i==at[j])&&opts.showatlbl) {IDefPlane=true;}
      }
      if (opts.showallatlbl||IDefPlane) {lbl="";} else {lbl="#";}
      lbl+="set label ";
      numst.str("");
      numst << (i+1);
      lbl+=numst.str();
      lbl+=" at ";
      numst.str("");
      numst << xproj << "," << yproj;
      lbl+=numst.str();
      lbl+=" '";
      lbl+=bn.atLbl[i];
      /* Adding an underscore and braces to the atome labels in order to make the numbers
         subindices */
      pos=lbl.find_last_not_of("0123456789");
      if (pos!=string::npos) {
         lbl.insert(pos+1,"_{");
         lbl.append("}");
      }
      lbl+="' front offset character -0.5,-0.15";
      for (int k=0; k<3; k++) {
         if (at[k]==i) {lbl+=" #* This atom was used for setting the plane!";}
      }
      gfil << lbl << endl;
   }
   FileUtils::WriteScrCharLine(gfil,'#');

   /* Here is enabled the logarithmic scale in the case of G, d or g */
   gfil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   
   
   string basenam,epsnam;
   basenam=tsvnam.substr(0,(tsvnam.length()-4));
   gfil << "set output '";
   gfil << basenam << ".eps'" << endl; //so far, the image's name is *.eps
   gfil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;
   if (p2p=='N'||p2p=='p'||p2p=='U') {
      //gfil << "plot '" << tsvnam << "' using 1:2:3:4:(sqrt($3*$3+$4*$4)) "
      //<< "with vectors head size 0.1,20,60 filled lc palette";
      gfil << "plot '" << tsvnam << "' every 1:1 using 1:2:(sqrt($3*$3+$4*$4)>VMXS? VMXS*$3/sqrt($3*$3+$4*$4) : $3):(sqrt($3*$3+$4*$4)>VMXS ? VMXS*$4/sqrt($3*$3+$4*$4) : $4):(sqrt($3*$3+$4*$4)) "
           << "with vectors head size 0.1,20,60 filled lc palette";
   } else {
      gfil << "plot '" << tsvnam << "' with image,\\\n";
   }
   
   if (!opts.showcont) {gfil << "#";}
   gfil << "'" << cntname << "' w l lt -1 lw 1  #Activate this line to show contours" << endl;
   
    // Prints some guidelines to produce a trimmed pdf from the eps image.
   epsnam=tsvnam.substr(0,(tsvnam.length()-3));
   epsnam.append("eps");
#if ( defined(__APPLE__)||defined(__linux__) )
   gfil << "\nsystem '" << "rm " << cntname << "'" << endl;
#endif
   gfil << "#" << endl << endl << "q" << endl << endl;
   HelpersPlot::addGnuplotRenderingHelpEPS2PDF(gfil,gnpnam,epsnam);
   gfil.close();

   if ( opts.mkplt ) {
      GnuplotTools::RenderGnpFile(gnpnam,(!opts.kpgnp));
      GnuplotTools::eps2pdf(epsnam);
   }
   cout << "Done." << endl;
   return;
}
void HelpersPlot::addGnuplotRenderingHelpEPS2PDF(ofstream &ofil,const string &gnpnam,\
      const string &epsnam) {
   FileUtils::WriteScrCharLine(ofil,'#');
   ofil << "#                 END OF GNUPLOT COMMANDS" << endl;
   FileUtils::WriteScrCharLine(ofil,'#');
   ofil << "#If you want to reconstruct the plot using this file, type:" << endl
        << "#gnuplot " << gnpnam << endl
        << "#dtkepstopdf " << epsnam << endl;
}

