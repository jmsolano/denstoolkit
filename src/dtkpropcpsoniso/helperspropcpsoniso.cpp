#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <fstream>
using std::ofstream;
#include <vector>
using std::vector;
#include "helperspropcpsoniso.h"
#include "screenutils.h"
#include "fileutils.h"
#include "commonhelpers.h"
#include "matrixvectoroperations3d.h"
#include "palette.h"
#include "atomradiicust.h"
#include "gaussiancube.h"

#ifndef PROPCPSONISOCOORDSEPS2
#define PROPCPSONISOCOORDSEPS2 1.0e-12
#endif
#ifndef PROPCPSONISOISOPRECISIONFACTOR
#define PROPCPSONISOISOPRECISIONFACTOR 1.0e-4
#endif
#ifndef PROPCPSONISOINISTEP
#define PROPCPSONISOINISTEP 0.2
#endif
#ifndef DEBUG
#define DEBUG 0
#endif

void HelpersPropCPsOnIso::GetCenterIndexAndVectors( char *argv[],\
      const OptionFlags &options,const BondNetWork &bn,int &cat,\
      vector<double> &xc,vector<double> &xd) {
   cat=0;
   if ( options.setcentat ) { cat=std::stoi(string(argv[options.setcentat]))-1; }
   for ( int i=0 ; i<3 ; ++i ) { xc[i]=bn.R[cat][i]; }
   if ( options.setdirat1 ) {
      int a1=std::stoi(string(argv[options.setdirat1]))-1;
      for ( int i=0 ; i<3 ; ++i ) { xd[i]=(bn.R[cat][i]-bn.R[a1][i]); }
   } else if ( options.setdirat2 ) {
      int a1=std::stoi(string(argv[options.setdirat2]))-1;
      int a2=std::stoi(string(argv[options.setdirat2+1]))-1;
      for ( int i=0 ; i<3 ; ++i ) {
         xd[i]=(bn.R[cat][i]-0.5e0*(bn.R[a1][i]+bn.R[a2][i]));
      }
   } else if ( options.setdirat3 ) {
      int a1=std::stoi(string(argv[options.setdirat3]))-1;
      int a2=std::stoi(string(argv[options.setdirat3+1]))-1;
      int a3=std::stoi(string(argv[options.setdirat3+2]))-1;
      for ( int i=0 ; i<3 ; ++i ) {
         xd[i]=(bn.R[cat][i]-(bn.R[a1][i]+bn.R[a2][i])+bn.R[a3][i])/3.0e0;
      }
   }
   MatrixVectorOperations3D::Normalize(xd);
}
bool HelpersPropCPsOnIso::HaveIncompatibleOptions(int argc,char *argv[],OptionFlags &opt) {
   int count=0;
   if ( opt.setdirat1 ) { ++count; }
   if ( opt.setdirat2 ) { ++count; }
   if ( opt.setdirat3 ) { ++count; }
   if ( count>1 ) {
      string msg="Direction point is defined more than once.\n"
         "Please use only one option. Command line:\n";
      for ( int i=0 ; i<argc ; ++i ) {
         msg+=string(argv[i])+string(" ");
      }
      ScreenUtils::DisplayErrorMessage(msg);
      return true;
   }
   return false;
}
void HelpersPropCPsOnIso::ProjectGridOntoIsosurface(GaussWaveFunction &wf,\
      shared_ptr<SymmetricSurfaceGrid> g,const char prop,const double iso) {
   double f;
   vector<double> r(3);
   vector<size_t> v2erase;
   v2erase.clear();
   for ( size_t i=0 ; i<g->vertex.size() ; ++i ) {
      //cout << "(" << wf.EvalDensity(g->vertex[i][0],g->vertex[i][1],g->vertex[i][2]) << ")";
      f=SearchValueAlongLineDescending(g->center,g->vertex[i],wf,r,prop,iso);
      if ( f>0.0e0 ) {
         g->vertex[i][0]=r[0]; g->vertex[i][1]=r[1]; g->vertex[i][2]=r[2];
      } else {
         //ScreenUtils::DisplayWarningMessage(string("Erasing vertex ")+std::to_string(i));
         v2erase.push_back(i);
      }
   }
   for ( size_t i=0 ; i<v2erase.size() ; ++i ) {
      cout <<  v2erase[i] << '\n';
   }
   g->RemoveFacesUsingVertices(v2erase);
   g->ComputeCentroids();
   for ( size_t i=0 ; i<g->centroid.size() ; ++i ) {
      f=SearchValueAlongLineDescending(g->center,g->centroid[i],wf,r,prop,iso);
      if ( f>0.0e0 ) {
         g->centroid[i][0]=r[0]; g->centroid[i][1]=r[1]; g->centroid[i][2]=r[2];
      }
   }
}
double HelpersPropCPsOnIso::SearchValueAlongLineDescending(const vector<double> &c,\
      const vector<double> &r0,GaussWaveFunction &wf,vector<double> &ri,const char prop,const double iso) {
   double (GaussWaveFunction::* f)(double,double,double);
   switch ( prop ) {
      case 'd' :
         f=&GaussWaveFunction::EvalDensity;
         break;
      case 'V' :
         f=&GaussWaveFunction::EvalMolElecPot;
         break;
      default :
         f=&GaussWaveFunction::EvalDensity;
         break;
   }
   double h[3];
   h[0]=PROPCPSONISOINISTEP*(r0[0]-c[0]);
   h[1]=PROPCPSONISOINISTEP*(r0[1]-c[1]);
   h[2]=PROPCPSONISOINISTEP*(r0[2]-c[2]);
   vector<double> rl(3),rr(3);
   rl[0]=r0[0]-h[0]; rl[1]=r0[1]-h[1]; rl[2]=r0[2]-h[2];
   rr[0]=r0[0];      rr[1]=r0[1];      rr[2]=r0[2];
   double vl=(wf.*f)(rl[0],rl[1],rl[2]);
   double vr=(wf.*f)(rr[0],rr[1],rr[2]);
   //cout << "Before while: vl: " << vl << ", vr: " << vr << '\n';
   double dv=vr-vl;
   size_t count=0;
   while ( vr>iso && dv<0.0e0 && vr>0.000001 ) {
      rl[0]+=h[0]; rl[1]+=h[1]; rl[2]+=h[2];
      rr[0]+=h[0]; rr[1]+=h[1]; rr[2]+=h[2];
      vl=(wf.*f)(rl[0],rl[1],rl[2]);
      vr=(wf.*f)(rr[0],rr[1],rr[2]);
      dv=vr-vl;
      ++count;
   }
   if ( dv>0.0e0 ) {
      ScreenUtils::DisplayWarningMessage("dv>0!");
      cout << "After while(count=" << count << "): vl=" << vl << ", vr="
           << vr << '\n';
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return (-1.0e0);
   }
   vector<double> rm(3);
   rm[0]=0.5e0*(rl[0]+rr[0]); rm[1]=0.5e0*(rl[1]+rr[1]); rm[2]=0.5e0*(rl[2]+rr[2]);
   double vm=(wf.*f)(rm[0],rm[1],rm[2]);
   double dist2=MatrixVectorOperations3D::Distance2(rl,rr);
   //cout << "vl: " << vl << ", vr: " << vr << ",vm: " << vm << ", dist2: " << dist2 << '\n';
   double epsiso=iso*PROPCPSONISOISOPRECISIONFACTOR;
   //cout << "epsiso: " << epsiso << '\n';
   while ( dist2>PROPCPSONISOCOORDSEPS2 ) {
      //Check if iso is close enough to vm
      if ( fabs(vm-iso)<=epsiso ) { break; }
      if ( iso<vm ) { // If vm is smaller, ignore left half
         rl[0]=rm[0]; rl[1]=rm[1]; rl[2]=rm[2];
         vl=vm;
      } else { // If vm is greater, ignore right half
         rr[0]=rm[0]; rr[1]=rm[1]; rr[2]=rm[2];
         vr=vm;
      }
      rm[0]=0.5e0*(rl[0]+rr[0]); rm[1]=0.5e0*(rl[1]+rr[1]); rm[2]=0.5e0*(rl[2]+rr[2]);
      vm=(wf.*f)(rm[0],rm[1],rm[2]);
      dist2=MatrixVectorOperations3D::Distance2(rl,rr);
   }
   ri[0]=rm[0]; ri[1]=rm[1]; ri[2]=rm[2];
   //cout << "iso: " << vm << ", dist2: " << dist2 << '\n';
   vector<double> r(3);
   for ( size_t i=0 ; i<3 ; ++i ) { r[i]=rm[i]-c[i]; }
   return MatrixVectorOperations3D::Norm(r);
}
bool HelpersPropCPsOnIso::SearchCPsCap(shared_ptr<MeshGrid> g,\
      GaussWaveFunction &wf,vector<vector<double> > &rcp,vector<size_t>  &poscp,\
      vector<int> &sigcp,vector<double> &valcp,const char prop) {
   double (GaussWaveFunction::* f)(double,double,double);
   switch ( prop ) {
      case 'd' :
         f=&GaussWaveFunction::EvalDensity;
         break;
      case 'V' :
         f=&GaussWaveFunction::EvalMolElecPot;
         break;
      default :
         f=&GaussWaveFunction::EvalDensity;
         break;
   }
   size_t nf=g->face.size(),nv=g->vertex.size();
   vector<double> propAtCentroid(nf);
   vector<double> propAtVertex(nv);
   for ( size_t i=0 ; i<nf ; ++i ) {
      propAtCentroid[i]=(wf.*f)(g->centroid[i][0],g->centroid[i][1],g->centroid[i][2]);
   }
   for ( size_t i=0 ; i<nv ; ++i ) {
      propAtVertex[i]=(wf.*f)(g->vertex[i][0],g->vertex[i][1],g->vertex[i][2]);
   }
   for ( size_t i=0 ; i<nv ; ++i ) { g->value[i]=propAtVertex[i]; }
   size_t p0,p1,p2;
   double v0,v1,v2,vc;
   for ( size_t i=0 ; i<rcp.size() ; ++i ) { rcp[i].clear(); }
   rcp.clear();
   poscp.clear();
   valcp.clear();
   sigcp.clear();
   for ( size_t i=0 ; i<nf ; ++i ) {
      p0=g->face[i][0]; p1=g->face[i][1]; p2=g->face[i][2];
      vc=propAtCentroid[i]; v0=propAtVertex[p0]; v1=propAtVertex[p1]; v2=propAtVertex[p2];
      if ( vc>v0 && vc>v1 && vc>v2 ) {
         rcp.push_back(g->centroid[i]);
         poscp.push_back(i);
         valcp.push_back(propAtCentroid[i]);
         sigcp.push_back(-3);
      }
      if ( vc<v0 && vc<v1 && vc<v2 ) {
         rcp.push_back(g->centroid[i]);
         poscp.push_back(i);
         valcp.push_back(propAtCentroid[i]);
         sigcp.push_back(3);
      }
   }
   if ( valcp.size()==0 ) { return false; }
   return true;
}
bool HelpersPropCPsOnIso::SearchCPsIso(shared_ptr<MeshGrid> g,GaussWaveFunction &wf,\
      vector<vector<double> > &rcp,vector<size_t> &poscp,vector<int> &sigcp,\
      vector<double> &valcp,const char prop) {
   g->FindAllVertexNeighbours();
   double (GaussWaveFunction::* f)(double,double,double);
   switch ( prop ) {
      case 'd' :
         f=&GaussWaveFunction::EvalDensity;
         break;
      case 'V' :
         f=&GaussWaveFunction::EvalMolElecPot;
         break;
      default :
         f=&GaussWaveFunction::EvalDensity;
         break;
   }
   size_t nv=g->vertex.size();
   for ( size_t i=0 ; i<nv ; ++i ) {
      g->value[i]=(wf.*f)(g->vertex[i][0],g->vertex[i][1],g->vertex[i][2]);
   }
   for ( size_t i=0 ; i<rcp.size() ; ++i ) { rcp[i].clear(); }
   rcp.clear();
   poscp.clear();
   valcp.clear();
   sigcp.clear();
   bool immax,immin;
   size_t nnei,pos;
   double valatvtx,valatnei;
   for ( size_t i=0 ; i<nv ; ++i ) {
      immax=immin=true;
      nnei=g->vneigh2v[i].size();
      //cout << "nnei(" << i << "): " << nnei;
      valatvtx=g->value[i];
      if ( nnei<4 ) { continue; }
      for ( size_t j=0 ; j<nnei ; ++j ) {
         pos=g->vneigh2v[i][j];
         //cout << ' ' << pos;
         valatnei=g->value[pos];
         immax=immax&&(valatnei<valatvtx);
         immin=immin&&(valatnei>valatvtx);
      }
      //cout << '\n';
      if ( immax ) {
         rcp.push_back(g->vertex[i]);
         poscp.push_back(i);
         valcp.push_back(valatvtx);
         sigcp.push_back(-3);
      }
      if ( immin ) {
         rcp.push_back(g->vertex[i]);
         poscp.push_back(i);
         valcp.push_back(valatvtx);
         sigcp.push_back(3);
      }
   }
   if ( valcp.size()==0 ) { return false; }
   return true;
}
bool HelpersPropCPsOnIso::MakePovFile(const string &povname,OptionFlags &options,POVRayConfiguration &pvp,\
      BondNetWork &bn,shared_ptr<MeshGrid> grid,vector<vector<double> > &sp,const string &palname) {
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
   ofil << "#include \"colors.inc\"" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "//" << '\n';
#if DEBUG
   ofil << "//Code generated by HelpersNCI::MakePovFile (part of DTK)." << '\n';
#endif
   ofil << "//Below you can find some options to be parsed to povray" << '\n';
   ofil << "//set your custom values." << '\n';
   ofil << "//You can reconstruct the image using the script dtkpov2png" << '\n';
   ofil << "//" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "#declare GNUPlotAngle1=" << pvp.vecAngView[0] << ";" << '\n';
   ofil << "#declare GNUPlotAngle2=" << pvp.vecAngView[2] << ";" << '\n';
   ofil << "#declare YAngle=" << pvp.vecAngView[1] << ";" << '\n';
   bool tmpbool=true;
   ofil << "#declare DrawStandardBonds=" << tmpbool << ";" << '\n';
   tmpbool=options.drawiso;
   ofil << "#declare DrawIsosurface=" << tmpbool << ";" << '\n';
   ofil << "#declare TransmitAtomSphere=0.0;" << '\n';
   ofil << "#declare TransmitStdBondCylinder=0.0;" << '\n';
   ofil << "#declare TransmitCritPts=0.0;" << '\n';
   ofil << "#declare TransmitIsosurface=0.2;" << '\n';
   ofil << "#default { finish { specular 0.2 roughness 0.03 phong .1 } }" << '\n';
   FileUtils::WriteScrCharLine(ofil,'/',false);
   ofil << "// END OF CUSTOM OPTIONS" << '\n';
   FileUtils::WriteScrCharLine(ofil ,'/',false);
   if (!(bn.ImStp())) {bn.SetUpBNW();}
   if (!(bn.ballAndStickMode)) {bn.drawAtSize*=AUTOMATICBALLANDSTICKRATIO;}
#if DEBUG
   ScreenUtils::DisplayWarningMessage("In this version, calling CritPtNetWork::makePovFILE(...)\n\
         will overwrite the original coordinates of the critical points\n\
         and the coordinates on the bondnetwork object as well.");
#endif
   //CenterMolecule(bn,iso);
   bn.CalcViewRadius();
   //cout << "rView: " << bn.rView << '\n';
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
   ofil << "global_settings { ambient_light White }" << '\n';
   ofil << "\nbackground { color < 0, 0.5, 0.7 > }\n" << '\n';
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
   for (int i=0; i<pvp.nLightSources; i++) {
      pvp.WriteLightSource(ofil,i,0.5,"  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >\n");
   }
   ofil << "camera {" << '\n';
   //ofil << "  orthographic" << '\n';
   ofil << "  up < 0, 1, 0 >" << '\n';
   ofil << "  right < -4/3, 0, 0 >" << '\n';
   ofil << "  location ";
   HelpersPOVRay::WriteVector(ofil,pvp.locCam[0],pvp.locCam[1],pvp.locCam[2]);
   ofil << '\n' << "  look_at ";
   HelpersPOVRay::WriteVector(ofil,pvp.lookAtCam[0],pvp.lookAtCam[1],pvp.lookAtCam[2]);
   ofil << '\n' << "  rotate < GNUPlotAngle1, YAngle, GNUPlotAngle2 >";
   ofil << '\n' << "}" << '\n';
   CommonHelpers::PutNuclei(ofil,bn,pvp.currIndLev,"TransmitAtomSphere");
   ofil << "#if(DrawStandardBonds)" << '\n';
   CommonHelpers::PutBonds(ofil,bn,pvp.currIndLev,"TransmitStdBondCylinder");
   ofil << "#end\n//end if DrawStandardBonds" << '\n';
   vector<double> rgb{0.1,0.1,1.0};
   vector<vector<double> > textures=ComputeTextures(grid,-0.025,0.01,palname);
   if ( grid->UseNormals() ) {
      //HelpersPOVRay::WriteMesh2SingleRGB(ofil,grid->vertex,grid->normal,grid->face,0,rgb,"transmit TransmitIsosurface");
      HelpersPOVRay::WriteMesh2WithTextures(ofil,grid->vertex,grid->normal,textures,grid->face,0,"transmit TransmitIsosurface");
   } else {
      HelpersPOVRay::WriteMesh2SingleRGB(ofil,grid->vertex,grid->face,0,rgb,"transmit TransmitIsosurface");
   }
   CommonHelpers::PutSpecialSpheres(ofil,pvp.currIndLev,sp,"TransmitCritPts");
   ofil.close();
   if ( options.mkpng ) {
      string cmd="dtkpov2png "+povname+" 2>/dev/null";
      system(cmd.c_str());
   }
   return true;
}
bool HelpersPropCPsOnIso::ComputeNormalsAtVertices(shared_ptr<MeshGrid> g,GaussWaveFunction &wf,\
      const char prop) {
   void (GaussWaveFunction::* f)(double,double,double,double&,double&,double&,double&);
   switch ( prop ) {
      case 'd' :
         f=&GaussWaveFunction::EvalRhoGradRho;
         break;
      default :
         f=&GaussWaveFunction::EvalRhoGradRho;
         break;
   }
   //double vl=(wf.*f)(rl[0],rl[1],rl[2]);
   double rho;
   size_t nv=g->vertex.size();
   if ( nv!=g->normal.size() ) {
      ScreenUtils::DisplayWarningMessage("Needs resizing normals!");
      return false;
   }
   for ( size_t i=0 ; i<nv ; ++i ) {
      (wf.*f)(g->vertex[i][0],g->vertex[i][1],g->vertex[i][2],rho,\
            g->normal[i][0],g->normal[i][1],g->normal[i][2]);
   }
   g->NormaliseNormals();
   return true;
}
vector<vector<double> > HelpersPropCPsOnIso::ComputeTextures(shared_ptr<MeshGrid> g,const double valmin,\
      const double valmax,const string &palname) {
   size_t nv=g->vertex.size();
   double tmp;
   double vmin=1.0e+50,vmax=-1.0e+50;
   double mean=0.0e0l;
   for ( size_t i=0 ; i<g->value.size() ; ++i ) {
      tmp=g->value[i];
      mean+=tmp;
      if ( tmp<vmin ) { vmin=tmp; }
      if ( tmp>vmax ) { vmax=tmp; }
   }
   cout << "vecvmin: " << vmin << ", vecvmax: " << vmax << ", mean: " << (mean/double(nv)) << '\n';
   if ( valmin>vmin ) { vmin=valmin; }
   if ( valmax<vmax ) { vmax=valmax; }
   cout << "valmin: " << valmin << ", valmax: " << valmax << '\n';
   cout << "vmin: " << vmin << ", vmax: " << vmax << '\n';
   double oodeltatimes255=255.0e0/(vmax-vmin);
   vector<vector<double> > t(nv);
   for ( size_t i=0 ; i<nv ; ++i ) { t[i].resize(3); }
   Palette pal;
   pal.SelectPalette(palname);
   size_t pos;
   for ( size_t i=0 ; i<nv ; ++i ) {
      tmp=g->value[i];
      tmp-=vmin;
      tmp*=oodeltatimes255;
      if ( tmp<0.0e0 ) { tmp=0.0e0; }
      if ( tmp>255.0e0 ) { tmp=255.0e0; }
      pos=size_t(tmp);
      pal.GetRGB(pos,t[i][0],t[i][1],t[i][2]);
      //cout << tmp << ", rgb: " << t[i][0] << ' ' << t[i][1] << ' ' << t[i][2] << '\n';
   }
   return t;
}
shared_ptr<MeshGrid> HelpersPropCPsOnIso::BuildCapMesh(int argc,\
      char *argv[],const OptionFlags &opt,GaussWaveFunction &wf,\
      BondNetWork &bn) {
   char isoprop='d';
   if ( opt.isoprop ) { isoprop=argv[opt.isoprop][0]; }
   int cAt;
   vector<double> xc(3),xd(3);
   HelpersPropCPsOnIso::GetCenterIndexAndVectors(argv,opt,bn,cAt,xc,xd);
   cout << "xc: " << xc[0] << ' ' << xc[1] << ' ' << xc[2] << '\n';
   cout << "xd: " << xd[0] << ' ' << xd[1] << ' ' << xd[2] << '\n';
   shared_ptr<SymmetricSurfaceGrid> grid=shared_ptr<SymmetricSurfaceGrid>(new SymmetricSurfaceGrid);
   int refCapMeshLevel=4;
   if ( opt.refinemesh ) {
      refCapMeshLevel=std::stoi(string(argv[opt.refinemesh]));
      if ( refCapMeshLevel<0 ) {
         ScreenUtils::DisplayWarningMessage("Requesting a negative number of refining iterations.");
         refCapMeshLevel=0;
      }
      if ( refCapMeshLevel>8 ) {
         ScreenUtils::DisplayWarningMessage("This level of refinement may cause numerical issues.");
      }
   }
   grid-> SetupSphereIcosahedron(refCapMeshLevel);
   grid->Translate(xc);
   grid->ScaleVertices(GetAtomicVDWRadius(bn.atNum[cAt]));
   grid->TrimFacesCentroidDotProdGreaterThanZero(xd);

   /* Looking for partial isosurface  */
   cout << "Computing cap isosurface...\n";
   HelpersPropCPsOnIso::ProjectGridOntoIsosurface(wf,grid,isoprop,0.001);
   cout << "Done\n";
   return grid;
}
shared_ptr<MeshGrid> HelpersPropCPsOnIso::BuildMeshFromCube(int argc,\
      char *argv[],const OptionFlags &opt,GaussWaveFunction &wf,\
      BondNetWork &bn) {
   if ( !opt.isofromcube ) {
      ScreenUtils::DisplayErrorMessage("Requesting isofromcube when isofromcube==false.");
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return nullptr;
   }
   string cubfnam=argv[opt.isofromcube];
   cout << "Extracting isosurface..." << std::endl;
   GaussianCube gc(cubfnam);
   if ( !gc.CubeLoaded() ) {
      ScreenUtils::DisplayErrorMessage(string("The file '")+cubfnam\
            +string("' could not be opened."));
      cout << __FILE__ << ", line: " << __LINE__ << '\n';
      return nullptr;
   }
   shared_ptr<Isosurface> iso=shared_ptr<Isosurface>(new Isosurface);
   //iso->UseTetrahedrons(true);
   double isovalue=0.001e0;
   if ( opt.setisovalue ) {
      isovalue=std::stod(string(argv[opt.setisovalue]));
   }
   iso->ExtractMarchingCubes(gc,isovalue);
   return iso;
}
shared_ptr<MeshGrid> HelpersPropCPsOnIso::BuildMesh(int argc,\
      char *argv[],const OptionFlags &opt,GaussWaveFunction &wf,\
      BondNetWork &bn) {
   if ( opt.isofromcube ) {
      return BuildMeshFromCube(argc,argv,opt,wf,bn);
   } else {
      return BuildCapMesh(argc,argv,opt,wf,bn);
   }
   ScreenUtils::DisplayErrorMessage("For some reason, HelpersPropCPsOnIso::BuildMesh failed!");
   cout << __FILE__ << ", line: " << __LINE__ << '\n';
   return nullptr;
}

