#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "helpersplot.h"
#include "solcubetools.h"
#include "mymemory.h"
#include "../common/stringtools.h"
#include "../common/gnuplottools.h"

void HelpersPlot::MakeLineGnuplotFile(optFlags &opts, string &gnpn,string &outn,char thefield) {
   ofstream ofil;
   ofil.open(gnpn.c_str());
   
   /* Choosing the label (legend) for the plot */
   string plbl;
   switch ( thefield ) {
      case 'd' :
         plbl="{/Symbol r}(p)";
         break;
      case 'K' :
         plbl="{/Bold K}(p)";
         break;
      default :
         break;
   }
   
   string line=StringTools::GetEnhancedEpsTitle(outn);
   ofil << "set title '" << line << "'" << endl;
   
   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;
   
   string pdfname=gnpn;
   FileUtils::ReplaceExtensionOfFileName(pdfname,string("pdf"));

   ofil << "set output '|epstopdf --filter --outfile=" << pdfname << "'" << endl;
   ofil << "plot '" << outn << "' with lines ls 1 lw 2 title '" << plbl << "'" << endl;
   
   ofil.close();
   GnuplotTools::RenderGnpFile(gnpn,!(opts.kpgnp));
   return;
}
void HelpersPlot::MakePlaneGnuplotFile(optFlags &opts, string &gnpn,string &outn,double dimparam,char thefield) {
   ofstream ofil;
   ofil.open(gnpn.c_str());
   
   /* Choosing the label (legend) for the plot */
   
   string plbl="";
   switch ( thefield ) {
      case 'd':
         plbl="{/Symbol r}(p)";
         break;
      case 'K':
         plbl="{/Bold K}(p)";
         break;
      default :
         break;
   }
   
   ofil << "reset" << endl;
   
   /* In this part the name is scanned for possible occurrings of the character '_'.
    For a proper display in the eps file, it has to be changed to "\_" */
   string line="";
   for (size_t i=0; i<outn.length(); i++) {
      if (outn[i]=='_') {line+='\\';}
      line+=outn[i];
   }
   ofil << "set title '" << line << "'" << endl;
   
   /* Adding an underscore and braces to the atome labels in order to make the numbers
    subindices */

   string contourtemp=StringTools::GenerateStrRandSeq(32);
   
   ofil << "set isosample 300, 300" << endl;
   ofil << "set contour base" << endl;
   ofil << "set cntrparam levels discrete 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0" <<endl;
   ofil << "unset surface" << endl;
   ofil << "set table '" << contourtemp << "'" << endl;
   ofil << "splot '" << outn << "'" << endl;
   ofil << "unset table" << endl;
   ofil << "reset" << endl;
   ofil << "unset key" << endl;
   ofil << "set size ratio 1" << endl;
   ofil << "set tmargin at screen 0.99\nset bmargin at screen 0.01\nset lmargin at screen 0.01" << endl;
   ofil << "dimparam=" << dimparam
   << " #Decrease this number to zoom in the plot" << endl;
   ofil << "set notics" << endl;
   ofil << "set xrange[-dimparam:dimparam]\nset yrange[-dimparam:dimparam]" << endl;
   ofil << "set style data lines" << endl;
   ofil << "set palette rgbformulae 33,13,10" << endl;
   ofil << "set cbtics #deactivate if you do not want tics on the colorbox scale" << endl;

   ofil << "set terminal postscript eps enhanced color fontscale 1.75 lw 2 dashlength 4" << endl;

   string epsname=gnpn;
   FileUtils::ReplaceExtensionOfFileName(epsname,string("eps"));
   string pdfname=gnpn;
   FileUtils::ReplaceExtensionOfFileName(pdfname,string("pdf"));

   ofil << "set output '" << epsname << "'" << endl;

   ofil << "plot '" << outn << "' with image,\\" << endl;
   if ( !opts.mkplt ) { ofil << '#'; }
   ofil << "'" << contourtemp << "' w l lt -1 lw 1 #Activate this line to show contours" << endl;
   
   ofil << endl << endl;
   GnuplotTools::AddCommandsToRemoveTemporaryFileFromGnuplotScript(ofil,contourtemp);
   ofil << "q" << endl << endl;
   FileUtils::WriteScrCharLine(ofil,'#');
   ofil << "#                 END OF GNUPLOT COMMANDS" << endl;
   FileUtils::WriteScrCharLine(ofil,'#');
   ofil << "#If you want to reconstruct the plot using this file, type:" << endl
   << "#gnuplot " << gnpn << endl
   << "#dtkeps2pdf " << epsname << endl;
   ofil.close();

   GnuplotTools::RenderGnpFile(gnpn,(!(opts.kpgnp)));
   GnuplotTools::eps2pdf(epsname);

#if (defined(__APPLE__)||defined(__linux__)||defined(__CYGWIN__))
   line=string("rm ")+epsname;
   system(line.c_str());
#endif
   cout << "Done." << endl;
   return;
}
void HelpersPlot::MakeLineDatFile(optFlags &opts,string &datnam,GaussWaveFunction &wf,int theaxis,\
      int npts,char thefield) {
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   ofstream ofile;
   ofile.open(datnam.c_str(),std::ios::out);
   double dx,dy,dz,px,py,pz;
   dx=dy=dz=0.0e0;
   px=py=pz=0.0e0;
   switch (theaxis) {
      case 1:
         dx=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         px=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  ofile << px << " " << wf.evalFTDensity(px,py,pz) << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  ofile << px << " " << wf.evalFTKineticEnergy(px,py,pz) << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         
         break;
      case 2:
         dy=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         py=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  ofile << py << " " << wf.evalFTDensity(px,py,pz) << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            case 'K' :
              for (int i=0; i<npts; i++) {
                  ofile << py << " " << wf.evalFTKineticEnergy(px,py,pz) << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
      case 3:
         dz=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         pz=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  ofile << pz << " " << wf.evalFTDensity(px,py,pz) << endl;
                  pz+=dz;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            case 'K' :
              for (int i=0; i<npts; i++) {
                  ofile << pz << " " << wf.evalFTKineticEnergy(px,py,pz) << endl;
                  pz+=dz;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
              break;
            default :
               break;
         }
         break;
      default:
         break;
   }
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(100);
#endif
   cout << endl;
   ofile.close();
}
void HelpersPlot::MakePlaneTsvFile(optFlags &opts,string &tsvnam,GaussWaveFunction &wf,int theplane,\
      int npts,char thefield) {
   ofstream ofile;
   ofile.open(tsvnam.c_str(),std::ios::out);
   double dx,dy,dz,px,py,pz;
   dx=dy=dz=0.0e0;
   px=py=pz=0.0e0;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch (theplane) {
      case 1:
         dx=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         dy=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         px=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  py=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << py << "\t" << wf.evalFTDensity(px,py,pz) << endl;
                     py+=dy;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  py=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << py << "\t" << wf.evalFTKineticEnergy(px,py,pz) << endl;
                     py+=dy;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
      case 2:
         dx=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         dz=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         px=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << pz << "\t" << wf.evalFTDensity(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << px << "\t" << pz << "\t" << wf.evalFTKineticEnergy(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  px+=dx;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
      case 3:
         dy=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         dz=2.0e0*DEFAULTMAXVALUEOFP/double(npts-1);
         py=-1.0e0*DEFAULTMAXVALUEOFP;
         switch ( thefield ) {
            case 'd' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << py << "\t" << pz << "\t" << wf.evalFTDensity(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            case 'K' :
               for (int i=0; i<npts; i++) {
                  pz=-1.0e0*DEFAULTMAXVALUEOFP;
                  for (int j=0; j<npts; j++) {
                     ofile << py << "\t" << pz << "\t" << wf.evalFTKineticEnergy(px,py,pz) << endl;
                     pz+=dz;
                  }
                  ofile << endl;
                  py+=dy;
#if USEPROGRESSBAR
                  ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((npts-1))));
#endif
               }
               break;
            default :
               break;
         }
         break;
         default:
            break;
      }
      cout << endl;
      ofile.close();
}
void HelpersPlot::MakeCubeFile(optFlags &opts,string &cubnam,GaussWaveFunction &wf,int npts,\
      char thefield,string &strfield) {
   string comments="#Property: ";
   switch ( thefield ) {
      case 'd' :
         comments+="Momentum Density.";
         break;
      case 'K' :
         comments+="Kinetic Energy Density (in Momentum Space).";
         break;
      default :
         break;
   }
   double px,py,pz;
   px=py=pz=0.0e0;
   int boxnpts[3];
   for (int i=0; i<3; i++) {boxnpts[i]=npts;}
   double xin[3],delta[3][3];
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         delta[i][j]=0.0e0;
      }
      delta[i][i]=2.0e0*DEFAULTMAXVALUEOFP/double(boxnpts[i]-1);
      xin[i]=-DEFAULTMAXVALUEOFP;
   }
   ofstream ofile;
   ofile.open(cubnam.c_str(),std::ios::out);
   writeCubeHeader(ofile,wf.title[0],comments,boxnpts,xin,delta,wf.nNuc,wf.atCharge,wf.R);
   double *prop1d;
   MyMemory::Alloc1DRealArray("prop1d",boxnpts[2],prop1d);
   cout << "The size of the grid will be " << boxnpts[0] << " x "
      << boxnpts[1] << " x " << boxnpts[2] << endl;
   cout << "The total number of points that will be computed is "
      << boxnpts[0]*boxnpts[1]*boxnpts[2] << endl;
   cout << "Evaluating and writing " << strfield << " on a cube..." << endl;
#if USEPROGRESSBAR
   ScreenUtils::PrintProgressBar(0);
#endif
   switch ( thefield ) {
      case 'd' :
         px=xin[0];
         for (int i=0; i<boxnpts[0]; i++) {
            py=xin[1];
            for (int j=0; j<boxnpts[1]; j++) {
               pz=xin[2];
               for (int k=0; k<boxnpts[2]; k++) {
                  prop1d[k]=wf.evalFTDensity(px,py,pz);
                  //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
                  pz+=delta[2][2];
               }
               writeCubeProp(ofile,boxnpts[2],prop1d);
               py+=delta[1][1];
            }
            px+=delta[0][0];
#if USEPROGRESSBAR
            ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((boxnpts[0]-1))));
#endif
         }
         cout << endl;
         break;
      case 'K' :
         px=xin[0];
         for (int i=0; i<boxnpts[0]; i++) {
            py=xin[1];
            for (int j=0; j<boxnpts[1]; j++) {
               pz=xin[2];
               for (int k=0; k<boxnpts[2]; k++) {
                  prop1d[k]=wf.evalFTKineticEnergy(px,py,pz);
                  //if (prop1d[k]<1.0e-20) {prop1d[k]=0.0e0;}
                  pz+=delta[2][2];
               }
               writeCubeProp(ofile,boxnpts[2],prop1d);
               py+=delta[1][1];
            }
            px+=delta[0][0];
#if USEPROGRESSBAR
            ScreenUtils::PrintProgressBar(int(100.0e0*double(i)/double((boxnpts[0]-1))));
#endif
         }
         cout << endl;
         break;
      default :
         break;
   }
   //writeCubeProp(ofstream &ofil,int dim,double* (&prop));
   ofile.close();
   MyMemory::Dealloc1DRealArray(prop1d);
}

