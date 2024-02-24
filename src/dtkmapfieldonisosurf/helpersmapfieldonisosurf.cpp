/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
#include <iomanip>
using std::scientific;
using std::setprecision;
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;
#include <memory>
using std::shared_ptr;
#include "screenutils.h"
#include "fileutils.h"
#include "helpersmapfieldonisosurf.h"
#include "commonhelpers.h"
#include "palette.h"
#include "matrixvectoroperations3d.h"
#include "dtkscalarfunction3d.h"

bool HelpersMapFieldOnIsoSurf::ComputeFieldAtCentroids(GaussWaveFunction &wf,Isosurface &iso,\
      const char prop) {
   size_t n=iso.face.size();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   iso.ComputeCentroids();
   DTKScalarFunction f(wf);
   f.SetScalarFunction(prop);
   vector<double> v(n);
   double x,y,z;
   for ( size_t i=0 ; i<n ; ++i ) {
      x=iso.centroid[i][0];
      y=iso.centroid[i][1];
      z=iso.centroid[i][2];
      v[i]=f.f(x,y,z);
   }
   iso.SetProperty2Map(v);
   iso.UseColorMap(true);
   return true;
}
bool HelpersMapFieldOnIsoSurf::ComputeFieldAtVertices(GaussWaveFunction &wf,\
      Isosurface &iso,const char prop) {
   size_t n=iso.vertex.size();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   vector<double> v(n);
   DTKScalarFunction f(wf);
   f.SetScalarFunction(prop);
   double x,y,z;
   for ( size_t i=0 ; i<n ; ++i ) {
      x=iso.vertex[i][0];
      y=iso.vertex[i][1];
      z=iso.vertex[i][2];
      v[i]=f.f(x,y,z);
   }
   iso.SetProperty2Map(v);
   iso.UseColorMap(true);
   return true;
}
bool HelpersMapFieldOnIsoSurf::ComputeNormalsAtVertices(GaussWaveFunction &wf,\
      Isosurface &iso,const char prop) {
   size_t n=iso.vertex.size();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   if ( iso.vertex.size()!=iso.normal.size() ) {
      for ( size_t i=0 ; i<iso.normal.size() ; ++i ) { iso.normal[i].clear(); }
      iso.normal.resize(iso.vertex.size());
      for ( size_t i=0 ; i<iso.normal.size() ; ++i ) { iso.normal[i].resize(3); }
   }
   vector<double> x(3);
   DTKScalarFunction f(wf);
   f.SetScalarFunction(prop);
   double gs[3];
   for ( size_t i=0 ; i<n ; ++i ) {
      x[0]=iso.vertex[i][0]; x[1]=iso.vertex[i][1]; x[2]=iso.vertex[i][2];
      f.gradf(x[0],x[1],x[2],gs);
      iso.normal[i][0]=gs[0]; iso.normal[i][1]=gs[1]; iso.normal[i][2]=gs[2];
   }
   iso.NormaliseNormals();
   return true;
}
bool HelpersMapFieldOnIsoSurf::ComputeTriangleNormals(GaussWaveFunction &wf,Isosurface &iso) {
   size_t nv=iso.vertex.size();
   if ( nv==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   double x[3],gs[3];
   for ( size_t i=0 ; i<nv ; ++i ) {
      x[0]=iso.vertex[i][0]; x[1]=iso.vertex[i][1]; x[2]=iso.vertex[i][2];
      wf.EvalGradReducedDensityGradient(x[0],x[1],x[2],gs);
      iso.normal[i][0]=gs[0]; iso.normal[i][1]=gs[1]; iso.normal[i][2]=gs[2];
      //cout << sqrt(gs[0]*gs[0]+gs[1]*gs[1]+gs[2]*gs[2]) << '\n';
   }
   iso.NormaliseNormals();
   return true;
}
bool HelpersMapFieldOnIsoSurf::MakePovFile(const string &povname,POVRayConfiguration &pvp,BondNetWork &bn,Isosurface &iso,\
      OptionFlags &options,char *argv[]) {
   ofstream ofil(povname.c_str());
   if ( !ofil.good() ) {
      string msg="The file '";
      msg+=povname;
      msg+="' could not be opened!";
      ScreenUtils::DisplayErrorMessage(msg);
      ofil.close();
      return false;
   }
#if CHOOSEPOVVERSION36
   ofil << "#version 3.6; //Unless you know what you are doing, do not modify this line..." << '\n';
#endif
   ofil << "#include \"colors.inc\"" << endl;
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "//" << endl;
#if DEBUG
   ofil << "//Code generated by HelpersMapFieldOnIsoSurf::MakePovFile (part of DTK)." << endl;
#endif
   ofil << "//Below you can find some options to be parsed to povray" << endl;
   ofil << "//set your custom values." << endl;
   ofil << "//You can reconstruct the image using the script dtkpov2png" << endl;
   ofil << "//" << endl;
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "#declare GNUPlotAngle1=" << pvp.vecAngView[0] << ";" << endl;
   ofil << "#declare GNUPlotAngle2=" << pvp.vecAngView[2] << ";" << endl;
   ofil << "#declare YAngle=" << pvp.vecAngView[1] << ";" << endl;
   bool tmpbool=true;
   ofil << "#declare DrawStandardBonds=" << tmpbool << ";" << endl;
   ofil << "#declare TransmitAtomSphere=0.0;" << endl;
   ofil << "#declare TransmitStdBondCylinder=0.0;" << endl;
   ofil << "#default { finish { specular 0.3 roughness 0.03 phong 0.1 } }" << endl;
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "// END OF CUSTOM OPTIONS" << endl;
   FileUtils::WriteScrCharLine(ofil ,'/',false);
   if (!(bn.ImStp())) {bn.SetUpBNW();}
   if (!(bn.ballAndStickMode)) {bn.drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
#if DEBUG
   ScreenUtils::DisplayWarningMessage("In this version, calling CritPtNetWork::makePovFILE(...)\n\
         will overwrite the original coordinates of the critical points\n\
         and the coordinates on the bondnetwork object as well.");
#endif
   CenterMolecule(bn,iso);
   bn.CalcViewRadius();
   double camdist=2.5e0;
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
   if ( options.rotcam ) {
      HelpersMapFieldOnIsoSurf::AlignMolecule(pvp,bn,options,argv);
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
   pvp.inccolors=false;
   ofil << "global_settings { ambient_light White }" << endl;
   ofil << "\nbackground { color < 0, 0.5, 0.7 > }\n" << endl;
   //pvp.WriteLightSource(ofil,0,0.75," rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   for (int i=0; i<pvp.nLightSources; i++) {
      pvp.WriteLightSource(ofil,i,0.5," rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   }
   ofil << scientific << setprecision(10);
   ofil << "camera {" << endl;
   //ofil << "  orthographic" << '\n';
   ofil << "\n  location ";
   HelpersPOVRay::WriteVector(ofil,pvp.locCam[0],pvp.locCam[1],pvp.locCam[2]);
   ofil << "\n  up ";
   HelpersPOVRay::WriteVector(ofil,pvp.vecUp[0],pvp.vecUp[1],pvp.vecUp[2]);
   ofil << '\n' << "  right ";
   HelpersPOVRay::WriteVector(ofil,pvp.vecRight[0],pvp.vecRight[1],pvp.vecRight[2]);
   if ( options.orientcam3ats ) {
      ofil << '\n' << "  direction ";
      HelpersPOVRay::WriteVector(ofil,pvp.vecDir[0],pvp.vecDir[1],pvp.vecDir[2]);
   } else {
      ofil << endl << "  look_at ";
      HelpersPOVRay::WriteVector(ofil,pvp.lookAtCam[0],pvp.lookAtCam[1],pvp.lookAtCam[2]);
   }
   ofil << endl << "  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >";
   ofil << endl << "}" << endl;
   CommonHelpers::PutNuclei(ofil,bn,pvp.currIndLev,"TransmitAtomSphere");
   ofil << "#if(DrawStandardBonds)" << endl;
   CommonHelpers::PutBonds(ofil,bn,pvp.currIndLev,"TransmitStdBondCylinder");
   ofil << "#end\n//end if DrawStandardBonds" << endl;
   shared_ptr<Palette> pal=std::make_shared<Palette>();
   string palname="greens";
   if ( options.selectpalette ) { palname=argv[options.selectpalette]; }
   pal->SelectPalette(palname);
   if ( iso.UseNormals() ) {
      HelpersIsosurface::AddIsosurfacePOVMeshWithNormals(ofil,iso,pal,pvp.currIndLev);
   } else {
      HelpersIsosurface::AddIsosurfacePOVMeshNoNormals(ofil,iso,pal,pvp.currIndLev);
   }
   ofil.close();
   if ( options.mkpng ) {
      string cmd="dtkpov2png "+povname+" 2>/dev/null";
      system(cmd.c_str());
   }
   return true;
}
void HelpersMapFieldOnIsoSurf::CenterMolecule(BondNetWork &bn,Isosurface &iso) {
   vector<double> trn(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      trn[i]=-0.5e0*(bn.rmax[i]+bn.rmin[i]);
   }
   bn.CenterMolecule();
   iso.Translate(trn);
   cout << "isocenter: " << iso.center[0] << " " << iso.center[1] << " " << iso.center[2] << '\n';
}
bool HelpersMapFieldOnIsoSurf::AlignMolecule(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]) {
   if ( options.orientcam3ats ) { return AlignMolecule3Atoms(pvp,bn,options,argv); }
   return false;
}
bool HelpersMapFieldOnIsoSurf::AlignMolecule3Atoms(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]) {
   int aIdx=std::stoi(string(argv[options.orientcam3ats  ]))-1;
   int bIdx=std::stoi(string(argv[options.orientcam3ats+1]))-1;
   int cIdx=std::stoi(string(argv[options.orientcam3ats+2]))-1;
   vector<double> a(3),b(3),c(3),tmp(3);
   for ( size_t i=0 ; i<3 ; ++i ) { a[i]=bn.R[aIdx][i]; }
   for ( size_t i=0 ; i<3 ; ++i ) { b[i]=bn.R[bIdx][i]; }
   for ( size_t i=0 ; i<3 ; ++i ) { c[i]=bn.R[cIdx][i]; }
   vector<vector<double> > M=MatrixVectorOperations3D::GetRotationMatrix2AlignPassive(a,b,c);
   pvp.SelectStandardCameraVectors();
   pvp.ApplyRotationMatrixToCameraAndLightSources(M);
   return true;
}

