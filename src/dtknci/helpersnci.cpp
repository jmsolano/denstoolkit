#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <vector>
using std::vector;
#include <fstream>
using std::ofstream;
#include <memory>
using std::shared_ptr;
#include "screenutils.h"
#include "fileutils.h"
#include "helpersnci.h"
#include "commonhelpers.h"
#include "palette.h"

bool HelpersNCI::ComputeLambdaOnCentroids(GaussWaveFunction &wf,Isosurface &iso) {
   size_t n=iso.TrianglesSize();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   iso.ComputeCentroids();
   vector<double> v(n);
   double x,y,z;
   for ( size_t i=0 ; i<n ; ++i ) {
      x=iso.centroid[i][0];
      y=iso.centroid[i][1];
      z=iso.centroid[i][2];
      v[i]=wf.EvalNCILambda(x,y,z);
   }
   iso.SetProperty2Map(v);
   iso.UseColorMap(true);
   return true;
}
bool HelpersNCI::ComputeLambdaOnVertices(GaussWaveFunction &wf,Isosurface &iso) {
   size_t n=iso.vertices.size();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   vector<double> v(n);
   double x,y,z;
   for ( size_t i=0 ; i<n ; ++i ) {
      x=iso.vertices[i][0];
      y=iso.vertices[i][1];
      z=iso.vertices[i][2];
      v[i]=wf.EvalNCILambda(x,y,z);
   }
   iso.SetProperty2Map(v);
   iso.UseColorMap(true);
   return true;
}
bool HelpersNCI::ComputeNormalsAtVertices(GaussWaveFunction &wf,Isosurface &iso) {
   size_t n=iso.vertices.size();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   if ( iso.vertices.size()!=iso.normals.size() ) {
      for ( size_t i=0 ; i<iso.normals.size() ; ++i ) { iso.normals[i].clear(); }
      iso.normals.resize(iso.vertices.size());
      for ( size_t i=0 ; i<iso.normals.size() ; ++i ) { iso.normals[i].resize(3); }
   }
   vector<double> x(3);
   double gs[3];
   for ( size_t i=0 ; i<n ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) {
         x=iso.vertices[i];
         wf.EvalGradReducedDensityGradient(x[0],x[1],x[2],gs);
         iso.normals[i][0]=gs[0];
         iso.normals[i][1]=gs[1];
         iso.normals[i][2]=gs[2];
         //cout << sqrt(gs[0]*gs[0]+gs[1]*gs[1]+gs[2]*gs[2]) << '\n';
      }
   }
   return true;
}
bool HelpersNCI::ComputeTriangleNormals(GaussWaveFunction &wf,Isosurface &iso) {
   size_t n=iso.TrianglesSize();
   if ( n==0 ) {
      ScreenUtils::DisplayErrorMessage("There are zero triangles in the isosurface!"
            "\nNothing to be done.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return false;
   }
   vector<double> x(3),nv(3);
   double gs[3];
   for ( size_t i=0 ; i<n ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) {
         x=iso.Vertex(i,j);
         wf.EvalGradReducedDensityGradient(x[0],x[1],x[2],gs);
         nv[0]=gs[0]; nv[1]=gs[1]; nv[2]=gs[2];
         iso.SetNormal(i,j,nv);
         //cout << sqrt(gs[0]*gs[0]+gs[1]*gs[1]+gs[2]*gs[2]) << '\n';
      }
   }
   return true;
}
bool HelpersNCI::MakePovFile(const string &povname,POVRayConfiguration &pvp,BondNetWork &bn,Isosurface &iso,\
      const string &palname,bool render) {
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
   ofil << "//Code generated by HelpersNCI::MakePovFile (part of DTK)." << endl;
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
   ofil << "#default { finish { specular 0.2 roughness 0.03 phong .1 } }" << endl;
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
   //cout << "rView: " << bn.rView << endl;
   double camdist=2.5e0;
   for (int i=0; i<3; i++) {pvp.locCam[i]=0.0e0;}
   pvp.locCam[2]=camdist;
   for (int i=0; i<3; i++) {
      pvp.locCam[i]*=bn.rView;
      for (int j=0; j<2; j++) {
         pvp.lightSource[j][i]*=(bn.rView*2.0e0);
      }
   }
   pvp.inccolors=false;
   //pvp.writeHeader(ofil,false);
   ofil << "global_settings { ambient_light White }" << endl;
   ofil << "\nbackground { color < 0, 0.5, 0.7 > }\n" << endl;
   double zsep=0.5e0;
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
   for (int i=0; i<pvp.nLightSources; i++) {
      pvp.WriteLightSource(ofil,i,0.5,"  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   }
   ofil << "camera {" << endl;
   //ofil << "  orthographic" << '\n';
   ofil << "  up < 0, 1, 0 >" << endl;
   ofil << "  right < -4/3, 0, 0 >" << endl;
   ofil << "  location ";
   HelpersPOVRay::WriteVector(ofil,pvp.locCam[0],pvp.locCam[1],pvp.locCam[2]);
   ofil << endl << "  look_at ";
   HelpersPOVRay::WriteVector(ofil,pvp.lookAtCam[0],pvp.lookAtCam[1],pvp.lookAtCam[2]);
   ofil << endl << "  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >";
   ofil << endl << "}" << endl;
   CommonHelpers::PutNuclei(ofil,bn,pvp.currIndLev,"TransmitAtomSphere");
   ofil << "#if(DrawStandardBonds)" << endl;
   CommonHelpers::PutBonds(ofil,bn,pvp.currIndLev,"TransmitStdBondCylinder");
   ofil << "#end\n//end if DrawStandardBonds" << endl;
   shared_ptr<Palette> pal=std::make_shared<Palette>();
   pal->SelectPalette(palname);
   if ( iso.UseNormals() ) {
      HelpersIsosurface::AddIsosurfacePOVMeshWithNormals(ofil,iso,pal,pvp.currIndLev);
   } else {
      HelpersIsosurface::AddIsosurfacePOVMeshNoNormals(ofil,iso,pal,pvp.currIndLev);
   }
   ofil.close();
   if ( render ) {
      string cmd="dtkpov2png "+povname+" 2>/dev/null";
      system(cmd.c_str());
   }
   return true;
}
void HelpersNCI::CenterMolecule(BondNetWork &bn,Isosurface &iso) {
   vector<double> trn(3);
   for ( size_t i=0 ; i<3 ; ++i ) {
      trn[i]=0.5e0*(bn.rmax[i]+bn.rmin[i]);
   }
   bn.CenterMolecule();
   iso.Translate(trn);
}
