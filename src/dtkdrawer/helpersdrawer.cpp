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
#include <iomanip>
using std::scientific;
using std::setprecision;
#include <vector>
using std::vector;
#include "gausswavefunction.h"
#include "critptnetwork.h"
#include "commonhelpers.h"
#include "helpersdrawer.h"
#include "povraytools.h"
#include "fileutils.h"
#include "../common/matrixvectoroperations3d.h"

bool HelpersDrawer::MakePovFile(const string &basename,const OptionFlags &options,\
      char *argv[],BondNetWork &bn) {
   string inname=argv[options.infname];
   if ( FileUtils::ExtensionMatches(inname,"cpx") ) {
      return MakeCriticalPointsPovFile(basename,options,argv,bn,inname);
   }
   string povname=basename+string(".pov");
   ofstream ofil(povname);
   if ( !ofil.good() ) {
      string msg="The file '";
      msg+=povname;
      msg+="' could not be opened!";
      ScreenUtils::DisplayErrorMessage(msg);
      ofil.close();
      return false;
   }
   WritePovHeader(ofil);
   if (!(bn.ImStp())) {bn.SetUpBNW();}
   if (!(bn.ballAndStickMode)) {bn.drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
   bn.CalcViewRadius();

   POVRayConfiguration pvp;
   SetupPovConf(pvp,bn,options,argv);
   CommonHelpers::WriteAngleDeclarations(ofil,pvp);
   bool tmpbool=true;
   ofil << "#declare DrawStandardBonds=" << tmpbool << ";" << '\n';
   //tmpbool=false;
   //ofil << "#declare DrawIsosurface=" << tmpbool << ";" << '\n';
   ofil << "#declare TransmitAtomSphere=0.0;" << '\n';
   ofil << "#declare TransmitStdBondCylinder=0.0;" << '\n';
   /*
   ofil << "#declare TransmitIsosurface=" << ( options.transparentiso ? "0.4" : "0.0")
        << "; // set this between 0 (opaque) and 1 (completely transparent)." << '\n';
   // */
   WriteCameraSettings(ofil,pvp,options.align3at);
   WriteLightSources(ofil,pvp);
   WriteBondNetworkElements(ofil,bn,pvp,options,argv);
   ofil.close();
   tmpbool=options.verboseLevel && std::stoi(string(argv[options.verboseLevel]))>0;
   int ww=1200;
   if ( options.setpngwidth ) { ww=std::stoi(string(argv[options.setpngwidth])); }
   CommonHelpers::RenderPovfile(povname,tmpbool,ww);
   return true;
}
bool HelpersDrawer::MakeCriticalPointsPovFile(const string &basename,\
      const OptionFlags &options,char *argv[],BondNetWork &bn,\
      const string cpname) {
   GaussWaveFunction wf;
   CritPtNetWork cp(wf,bn);
   cp.ReadFromFile(cpname);
   cp.DrawBondGradPaths(true);
   string povname=basename+string(".pov");
   POVRayConfiguration pvc;
   ofstream ofil;
   ofil.open(povname.c_str(),std::ios::out);
   if (!(ofil.good())) {
      cout << "Error: File " << povname << "could not be opened...\n";
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   WritePovHeader(ofil);
   if (!(bn.ImStp())) {bn.SetUpBNW();}
   if (!(bn.ballAndStickMode)) {bn.drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
   bn.CalcViewRadius();

   POVRayConfiguration pvp;
   SetupPovConf(pvp,bn,options,argv);
   CommonHelpers::WriteAngleDeclarations(ofil,pvp);
   bool tmpbool=true;
   if (options.drawnuc) {tmpbool=true;} else {tmpbool=false;}
   ofil << "#declare DrawAtomTranspSpheres=" << (tmpbool?"true":"false") << ";" << '\n';
   tmpbool=options.cylbond;
   ofil << "#declare DrawStandardBonds=" << (tmpbool?"true":"false") << ';' << '\n';
   ofil << "#declare TransmitStdBondCylinder=0.7;" << '\n';
   bool drawacps=DrawCPofType('a',options,argv);
   bool drawbcps=DrawCPofType('b',options,argv);
   bool drawrcps=DrawCPofType('r',options,argv);
   bool drawccps=DrawCPofType('c',options,argv);
   bool drawbgps=DrawGPofType('b',options,argv);
   bool drawrgps=DrawGPofType('r',options,argv);
   bool drawcgps=DrawGPofType('c',options,argv);

   ofil << "#declare DrawAttractorCriticalPoints=" << (drawacps?"true":"false") << ';' << '\n';
   ofil << "#declare DrawBondCriticalPoints=" << (drawbcps?"true":"false") << ';' << '\n';
   ofil << "#declare DrawRingCriticalPoints=" << (drawrcps?"true":"false") << ';' << '\n';
   ofil << "#declare DrawCageCriticalPoints=" << (drawccps?"true":"false") << ';' << '\n';
   tmpbool=((options.drawgps)&&(options.tubestyle));
   ofil << "#declare DrawGradientPathSpheres=" << (tmpbool? "false":"true") << ";" << '\n';
   ofil << "#declare DrawGradientPathTubes=" << (tmpbool? "true":"false") << ";" << '\n';
   ofil << "//  Activation of \"DrawGradientPathSpheres\" requires deactivation of "
      << "\n//  \"DrawGradientPathTubes\", and vice versa." << '\n';
   double allcprad=bn.drawAtSize*CPNW_ATOMCRITICALPOINTSIZEFACTOR;
   ofil << "#declare TransmitAtomSphere=0.7;" << '\n';
   ofil << "#declare RadiusAllCriticalPoints=" << allcprad << ";" << '\n';
   ofil << "#declare ColorACP=rgb <0.0,0.0,0.0>;" << '\n';
   ofil << "#declare RadiusACP=RadiusAllCriticalPoints;" << '\n';
   ofil << "#declare ColorBCP=rgb <0.0,0.6,1.0>;" << '\n';
   ofil << "#declare RadiusBCP=RadiusAllCriticalPoints;" << '\n';
   ofil << "#declare ColorRCP=rgb <1.0,1.0,0.0>;" << '\n';
   ofil << "#declare RadiusRCP=RadiusAllCriticalPoints;" << '\n';
   ofil << "#declare ColorCCP=rgb <1.0,0.0,0.0>;" << '\n';
   ofil << "#declare RadiusCCP=RadiusAllCriticalPoints;" << '\n';
   ofil << "#declare ColorBGradPath=rgb <0.0,0.2,1.0>;" << '\n';
   ofil << "#declare ColorRGradPath=rgb <0.0,0.8,0.0>;" << '\n';
   ofil << "#declare ColorCGradPath=rgb <1.0,0.5,0.0>;" << '\n';
   ofil << "#default { finish { specular 0.3 roughness 0.03 phong .1 } }" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "//For the colors, instead of rgb < ... >, you may want to try Red, Yellow, ..." << '\n';
   ofil << "//  or any of the colors defined in \"colors.inc\"" << '\n';
   ofil << "//" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "// END OF CUSTOM OPTIONS" << '\n';
   FileUtils::WriteScrCharLine(ofil ,'/',false);
#if DEBUG
   ScreenUtils::DisplayWarningMessage("In this version, calling HelpersDrawer::MakeCriticalPointsPovFile(...)\n\
         will overwrite the original coordinates of the critical points\n\
         and the coordinates on the bondnetwork object as well.");
#endif

   WriteCameraSettings(ofil,pvp,options.align3at);
   WriteLightSources(ofil,pvp);
   WriteBondNetworkElements(ofil,bn,pvp,options,argv);

   FileUtils::WriteScrCharLine(ofil,'/',false);
   if (drawacps && cp.IKnowACPs()) {
      ofil << "#if(DrawAttractorCriticalPoints)" << '\n';
      int nACP=cp.nACP;
      for (int i=0; i<nACP; i++) {
         HelpersPOVRay::WriteSphere(ofil,0,cp.RACP[i][0],cp.RACP[i][1],cp.RACP[i][2],"RadiusACP",
               "ColorACP");
      }
      ofil << "#end\n//end if DrawAttractorCriticalPoints" << '\n';
      FileUtils::WriteScrCharLine(ofil,'/',false);
   }
   if (drawbcps && cp.IKnowBCPs()) {
      ofil << "#if(DrawBondCriticalPoints)" << '\n';
      int nBCP=cp.nBCP;
      for (int i=0; i<nBCP; i++) {
         HelpersPOVRay::WriteSphere(ofil,0,cp.RBCP[i][0],cp.RBCP[i][1],cp.RBCP[i][2],"RadiusBCP",
               "ColorBCP");
      }
      ofil << "#end\n//end if DrawBondCriticalPoints" << '\n';
      FileUtils::WriteScrCharLine(ofil,'/',false);
   }
   if (drawrcps && cp.IKnowRCPs()) {
      ofil << "#if(DrawRingCriticalPoints)" << '\n';
      int nRCP=cp.nRCP;
      for (int i=0; i<nRCP; i++) {
         HelpersPOVRay::WriteSphere(ofil,0,cp.RRCP[i][0],cp.RRCP[i][1],cp.RRCP[i][2],"RadiusRCP",
               "ColorRCP");
      }
      ofil << "#end\n//end if DrawBondCriticalPoints" << '\n';
      FileUtils::WriteScrCharLine(ofil,'/',false);
   }
   if (drawccps && cp.IKnowCCPs()) {
      ofil << "#if(DrawCageCriticalPoints)" << '\n';
      int nCCP=cp.nCCP;
      for (int i=0; i<nCCP; i++) {
         HelpersPOVRay::WriteSphere(ofil,0,cp.RCCP[i][0],cp.RCCP[i][1],cp.RCCP[i][2],"RadiusCCP",
               "ColorCCP");
      }
      ofil << "#end\n//end if DrawBondCriticalPoints" << '\n';
      FileUtils::WriteScrCharLine(ofil,'/',false);
   }

   FileUtils::WriteScrCharLine(ofil,'/',false);
   if (drawbgps&&cp.IKnowBGPs()) {
      double gprad=0.06;
      int npts,nBCP=cp.nBCP;
      ofil << "#if(DrawGradientPathSpheres)" << '\n';
      ofil << "union {" << '\n';
      for (int i=0; i<nBCP; ++i) {
         npts=cp.conBCP[i][2];
         HelpersPOVRay::WriteSphere(ofil,1,cp.RBGP[i][0][0],cp.RBGP[i][0][1],cp.RBGP[i][0][2], \
               gprad,"ColorBGradPath");
         for (int j=1; j<npts; j++) {
            HelpersPOVRay::WriteSphere(ofil,1,cp.RBGP[i][j][0],cp.RBGP[i][j][1],cp.RBGP[i][j][2], \
                  gprad,"ColorBGradPath");
         }
      }
      ofil << "}" << '\n';
      ofil << "#end\n//end if DrawGradientPathSpheres" << '\n';
      ofil << "#if(DrawGradientPathTubes)" << '\n';
      ofil << "union {" << '\n';
      for (int i=0; i<nBCP; i++) {
         npts=cp.conBCP[i][2];
         HelpersPOVRay::WriteSphere(ofil,1,cp.RBGP[i][0][0],cp.RBGP[i][0][1],cp.RBGP[i][0][2], \
               gprad,"ColorBGradPath");
         for (int j=1; j<npts; j++) {
            HelpersPOVRay::WriteSphere(ofil,1,cp.RBGP[i][j][0],cp.RBGP[i][j][1],cp.RBGP[i][j][2], \
                  gprad,"ColorBGradPath");
            HelpersPOVRay::WriteCylinder(ofil,1,\
                  cp.RBGP[i][j][0],cp.RBGP[i][j][1],cp.RBGP[i][j][2], \
                  cp.RBGP[i][j-1][0],cp.RBGP[i][j-1][1],cp.RBGP[i][j-1][2], \
                  gprad,"ColorBGradPath");
         }
      }
      ofil << "}" << '\n';
      ofil << "#end\n//end if DrawGradientPathTubes" << '\n';
   }
   if ( drawrgps&&cp.IKnowRGPs() ) {
      double gprad=0.045;
      int npts,currBcpPos,nRCP=cp.nRCP;
      ofil << "#if(DrawGradientPathTubes)" << '\n';
      ofil << "union {" << '\n';
      for (int rcpIdx=0; rcpIdx<nRCP; ++rcpIdx) {
         currBcpPos=0;
         while ( cp.conRCP[rcpIdx][0][currBcpPos]>=0 ) {
            npts=cp.conRCP[rcpIdx][1][currBcpPos];
            HelpersPOVRay::WriteSphere(ofil,1,cp.RRGP[rcpIdx][currBcpPos][0][0],\
                  cp.RRGP[rcpIdx][currBcpPos][0][1],\
                  cp.RRGP[rcpIdx][currBcpPos][0][2],gprad,"ColorRGradPath");
            for ( int j=1 ; j<npts ; ++j ) {
               HelpersPOVRay::WriteSphere(ofil,1,cp.RRGP[rcpIdx][currBcpPos][j][0],\
                     cp.RRGP[rcpIdx][currBcpPos][j][1],\
                     cp.RRGP[rcpIdx][currBcpPos][j][2],gprad,"ColorRGradPath");
               if ( j%2 ==0 ) {
               HelpersPOVRay::WriteCylinder(ofil,1,\
                     cp.RRGP[rcpIdx][currBcpPos][j][0],\
                     cp.RRGP[rcpIdx][currBcpPos][j][1],\
                     cp.RRGP[rcpIdx][currBcpPos][j][2],\
                     cp.RRGP[rcpIdx][currBcpPos][j-1][0],\
                     cp.RRGP[rcpIdx][currBcpPos][j-1][1],\
                     cp.RRGP[rcpIdx][currBcpPos][j-1][2],\
                     gprad,"ColorRGradPath");
               }
            }
            ++currBcpPos;
         }
      }
      ofil << "}" << '\n';
      ofil << "#end\n//end if DrawGradientPathTubes" << '\n';
      //-------------------
      ofil << "#if(DrawGradientPathSpheres)" << '\n';
      ofil << "union {" << '\n';
      for (int rcpIdx=0; rcpIdx<nRCP; ++rcpIdx) {
         currBcpPos=0;
         while ( cp.conRCP[rcpIdx][0][currBcpPos]>=0 ) {
            npts=cp.conRCP[rcpIdx][1][currBcpPos];
            for ( int j=0 ; j<npts ; ++j ) {
               HelpersPOVRay::WriteSphere(ofil,1,cp.RRGP[rcpIdx][currBcpPos][j][0],\
                     cp.RRGP[rcpIdx][currBcpPos][j][1],\
                     cp.RRGP[rcpIdx][currBcpPos][j][2],gprad,"ColorRGradPath");
            }
            ++currBcpPos;
         }
      }
      ofil << "}" << '\n';
      ofil << "#end\n//end if DrawGradientPathSpheres" << '\n';
   }
   if ( drawcgps&&cp.IKnowCGPs() ) {
      double gprad=0.045;
      int npts,currRcpPos,nCCP=cp.nCCP;
      ofil << "#if(DrawGradientPathTubes)" << '\n';
      ofil << "union {" << '\n';
      for (int ccpIdx=0; ccpIdx<nCCP; ++ccpIdx) {
         currRcpPos=0;
         while ( cp.conCCP[ccpIdx][0][currRcpPos]>=0 ) {
            npts=cp.conCCP[ccpIdx][1][currRcpPos];
            HelpersPOVRay::WriteSphere(ofil,1,cp.RCGP[ccpIdx][currRcpPos][0][0],\
                  cp.RCGP[ccpIdx][currRcpPos][0][1],\
                  cp.RCGP[ccpIdx][currRcpPos][0][2],gprad,"ColorCGradPath");
            for ( int j=1 ; j<npts ; ++j ) {
               if (j%3 != 1) {
               HelpersPOVRay::WriteSphere(ofil,1,cp.RCGP[ccpIdx][currRcpPos][j][0],\
                     cp.RCGP[ccpIdx][currRcpPos][j][1],\
                     cp.RCGP[ccpIdx][currRcpPos][j][2],gprad,"ColorCGradPath");
               }
               if (j%3 ==0) {
               HelpersPOVRay::WriteCylinder(ofil,1,\
                     cp.RCGP[ccpIdx][currRcpPos][j][0],\
                     cp.RCGP[ccpIdx][currRcpPos][j][1],\
                     cp.RCGP[ccpIdx][currRcpPos][j][2],\
                     cp.RCGP[ccpIdx][currRcpPos][j-1][0],\
                     cp.RCGP[ccpIdx][currRcpPos][j-1][1],\
                     cp.RCGP[ccpIdx][currRcpPos][j-1][2],\
                     gprad,"ColorCGradPath");
               }
            }
            ++currRcpPos;
         }
      }
      ofil << "}" << '\n';
      ofil << "#end\n//end if DrawGradientPathTubes" << '\n';
      //-------------------
      ofil << "#if(DrawGradientPathSpheres)" << '\n';
      ofil << "union {" << '\n';
      for (int ccpIdx=0; ccpIdx<nCCP; ++ccpIdx) {
         currRcpPos=0;
         while ( cp.conCCP[ccpIdx][0][currRcpPos]>=0 ) {
            npts=cp.conCCP[ccpIdx][1][currRcpPos];
            for ( int j=0 ; j<npts ; ++j ) {
               HelpersPOVRay::WriteSphere(ofil,1,cp.RCGP[ccpIdx][currRcpPos][j][0],\
                     cp.RCGP[ccpIdx][currRcpPos][j][1],\
                     cp.RCGP[ccpIdx][currRcpPos][j][2],gprad,"ColorCGradPath");
            }
            ++currRcpPos;
         }
      }
      ofil << "}" << '\n';
      ofil << "#end\n//end if DrawGradientPathSpheres" << '\n';
   }

   ofil.close();

   //Rendering
   tmpbool=options.verboseLevel && std::stoi(string(argv[options.verboseLevel]))>0;
   int ww=1200;
   if ( options.setpngwidth ) { ww=std::stoi(string(argv[options.setpngwidth])); }
   CommonHelpers::RenderPovfile(povname,tmpbool,ww);
   return true;
}
bool HelpersDrawer::AlignMolecule(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]) {
   if ( options.align3at ) { return AlignMolecule3Atoms(pvp,bn,options,argv); }
   return false;
}
bool HelpersDrawer::AlignMolecule3Atoms(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]) {
   int aIdx=std::stoi(string(argv[options.align3at  ]))-1;
   int bIdx=std::stoi(string(argv[options.align3at+1]))-1;
   int cIdx=std::stoi(string(argv[options.align3at+2]))-1;
   vector<double> a(3),b(3),c(3),tmp(3);
   for ( size_t i=0 ; i<3 ; ++i ) { a[i]=bn.R[aIdx][i]; }
   for ( size_t i=0 ; i<3 ; ++i ) { b[i]=bn.R[bIdx][i]; }
   for ( size_t i=0 ; i<3 ; ++i ) { c[i]=bn.R[cIdx][i]; }
#if 1
   vector<vector<double> > M=MatrixVectorOperations3D::GetRotationMatrix2AlignPassive(a,b,c);
   pvp.SelectStandardCameraVectors();
   pvp.ApplyRotationMatrixToCameraAndLightSources(M);
#else
   vector<vector<double> > M=MatrixVectorOperations3D::GetRotationMatrix2AlignActive(a,b,c);
   tmp[0]=0.0e0; tmp[1]=0.0e0; tmp[2]=1.0e0;
   vector<double> vpos=MatrixVectorOperations3D::MatrixVectorProduct(M,tmp);
   for ( int i=0 ; i<bn.nNuc ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) { tmp[j]=bn.R[i][j]; }
      vpos=MatrixVectorOperations3D::MatrixVectorProduct(M,tmp);
      for ( int j=0 ; j<3 ; ++j ) { bn.R[i][j]=vpos[j]; }
   }
#endif
   return true;
}
void HelpersDrawer::SetupPovConf(POVRayConfiguration &pvp,BondNetWork &bn,\
      const OptionFlags &options,char* argv[]) {
   double camdist=2.5e0;
   if ( options.setzoom ) { camdist*=(std::stod(string(argv[options.setzoom]))); }
   for (int i=0; i<3; i++) {pvp.locCam[i]=0.0e0;}
   pvp.locCam[2]=1.0e0;
   double zsep=1.5e0;
   pvp.lightSource[1][0]=zsep;
   pvp.lightSource[1][1]=zsep;
   pvp.lightSource[1][2]=1.0e0;
   pvp.AddLightSource(zsep,-zsep,1.0e0);
   pvp.AddLightSource(-zsep,zsep,1.0e0);
   pvp.AddLightSource(-zsep,-zsep,1.0e0);
   pvp.AddLightSource(0.0e0,1.0e0,0.0e0);
   pvp.AddLightSource(-1.0e0,0.0e0,0.0e0);
   for (int i=1; i<pvp.nLightSources; i++) {
      for (int j=0; j<3; j++) {pvp.lightSource[i][j]*=(bn.rView*4.0e0);}
   }
   for ( size_t i=0 ; i<3 ; ++i ) { pvp.vecAngView[i]=0.0e0; }
   if ( options.rotatemol ) {
      AlignMolecule(pvp,bn,options,argv);
      double angle;
      if ( options.rotX ) {
         angle=std::stod(string(argv[options.rotX]));
         CommonHelpers::RotateCameraAroundRight(pvp,angle);
      }
      if ( options.rotY ) {
         angle=std::stod(string(argv[options.rotY]));
         CommonHelpers::RotateCameraAroundUp(pvp,angle);
      }
      if ( options.rotZ ) {
         angle=std::stod(string(argv[options.rotZ]));
         CommonHelpers::RotateCameraAroundLocCam(pvp,angle);
      }
   }
   for (int i=0; i<3; i++) {
      pvp.locCam[i]*=bn.rView*camdist;
      for (int j=0; j<2; j++) {
         pvp.lightSource[j][i]*=(bn.rView*2.0e0);
      }
   }
   
}
void HelpersDrawer::WriteCameraSettings(ofstream &ofil,POVRayConfiguration &pvc,\
      bool lookat) {
   ofil << "camera {";
   ofil << "\n  location ";
   HelpersPOVRay::WriteVector(ofil,pvc.locCam[0],pvc.locCam[1],pvc.locCam[2]);
   ofil << "\n  up ";
   HelpersPOVRay::WriteVector(ofil,pvc.vecUp[0],pvc.vecUp[1],pvc.vecUp[2]);
   ofil << '\n' << "  right ";
   HelpersPOVRay::WriteVector(ofil,pvc.vecRight[0],pvc.vecRight[1],pvc.vecRight[2]);
   if ( lookat ) {
      ofil << '\n' << "  direction ";
      HelpersPOVRay::WriteVector(ofil,pvc.vecDir[0],pvc.vecDir[1],pvc.vecDir[2]);
   } else {
      ofil << '\n' << "  look_at ";
      HelpersPOVRay::WriteVector(ofil,pvc.lookAtCam[0],pvc.lookAtCam[1],pvc.lookAtCam[2]);
   }
   ofil << '\n' << "  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >";
   ofil << '\n' << "}" << '\n';
}
void HelpersDrawer::WriteLightSources(ofstream &ofil,POVRayConfiguration &pvc) {
   ofil << "#default { finish { specular 0.2 roughness 0.03 phong .1 } }" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "// END OF CUSTOM OPTIONS" << '\n';
   FileUtils::WriteScrCharLine(ofil ,'/',false);
   ofil << "global_settings { ambient_light White }" << '\n';
   ofil << "\nbackground { color < 0, 0.5, 0.7 > }\n" << '\n';
   for (int i=0; i<pvc.nLightSources; i++) {
      pvc.WriteLightSource(ofil,i,0.5,"  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   }
}
void HelpersDrawer::WritePovHeader(ofstream &ofil) {
#if CHOOSEPOVVERSION36
   ofil << "#version 3.6; //Unless you know what you are doing, do not modify this line..." << '\n';
#endif
   ofil << "#include \"colors.inc\"" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "//" << '\n';
   ofil << "//Code generated by HelpersDrawerNCI::MakePovFile (part of DTK)." << '\n';
   ofil << "//Below you can find some options to be parsed to povray" << '\n';
   ofil << "//set your custom values." << '\n';
   ofil << "//You can reconstruct the image using the script dtkpov2png" << '\n';
   ofil << "//" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
}
void HelpersDrawer::WriteBondNetworkElements(ofstream &ofil,BondNetWork &bn,\
      POVRayConfiguration &pvc,const OptionFlags &options,char *argv[]) {
   CommonHelpers::PutNuclei(ofil,bn,pvc.currIndLev,"TransmitAtomSphere",options.cpkview);
   ofil << "#if(DrawStandardBonds)" << '\n';
   CommonHelpers::PutBonds(ofil,bn,pvc.currIndLev,"TransmitStdBondCylinder");
   ofil << "#end\n//end if DrawStandardBonds" << '\n';
}
bool HelpersDrawer::DrawCPofType(char t,const OptionFlags &options,char *argv[]) {
   if ( options.drawcps ) {
      if ( !options.selectcps2draw ) {
         return true;
      } else {
         string list=argv[options.selectcps2draw];
         char c=std::tolower(t);
         for ( size_t i=0 ; i<list.size() ; ++i ) { list[i]=std::tolower(list[i]); }
         if ( list.find(c)==string::npos ) {
            return false;
         } else {
            return true;
         }
      }
   }
   return false;
}
bool HelpersDrawer::DrawGPofType(char t,const OptionFlags &options,char *argv[]) {
   if ( options.drawgps ) {
      if ( !options.selectgps2draw ) {
         return true;
      } else {
         string list=argv[options.selectgps2draw];
         char c=std::tolower(t);
         for ( size_t i=0 ; i<list.size() ; ++i ) { list[i]=std::tolower(list[i]); }
         if ( list.find(c)==string::npos ) {
            return false;
         } else {
            return true;
         }
      }
   }
   return false;
}


