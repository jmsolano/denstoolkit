#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <iomanip>
using std::scientific;
using std::setprecision;
#include "povraytools.h"
#include "mymemory.h"

POVRayConfiguration::POVRayConfiguration() {
   angView=67.5e0;
   vecAngView[0]=angView; vecAngView[1]=vecAngView[2]=0.0e0;
   skyCam[0]=0.0e0; skyCam[1]=0.0e0; skyCam[2]=1.0e0;
   locCam[0]=0.0e0; locCam[1]=0.0e0; locCam[2]=2.3e0;
   lookAtCam[0]=0.0e0; lookAtCam[1]=0.0e0; lookAtCam[2]=0.0e0;
   vecUp[0]=0.0e0; vecUp[1]=1.0e0; vecUp[2]=0.0e0;
   vecRight[0]=-4.0e0/3.0e0; vecRight[1]=0.0e0; vecRight[2]=0.0e0;
   vecDir[0]=0.0e0; vecDir[1]=0.0e0; vecDir[2]=1.0e0;
   nLightSources=2;
   MyMemory::Alloc2DRealArray("lightSource",nLightSources,3,lightSource);
   defLightSource[0]=1.0; defLightSource[1]=1.0; defLightSource[2]=1.0;
   for (int i=0; i<3; i++) {
      lightSource[0][i]=defLightSource[i];
      lightSource[1][i]=defLightSource[i];
   }
   lightSource[1][2]*=(-1.0);
   backGroundColor[0]=0.0; backGroundColor[1]=0.5; backGroundColor[2]=0.7;
   currIndLev=0;
   shine=false;
   inccolors=setversion36=true;
}
POVRayConfiguration::~POVRayConfiguration() {
   MyMemory::Dealloc2DRealArray(lightSource,nLightSources);
}
bool POVRayConfiguration::WriteHeader(ofstream &ofil,bool placecam) {
   string thetabs;
   if ( setversion36 ) {
      ofil << "#version 3.6; "
           << "//Unless you know what you are doing, do not modify "
           << "this line..." << '\n';
   }
   WriteIncColors(ofil);
   if (shine) {
      ofil << '\n' << "#default { finish { specular 0.3 roughness 0.03 phong .1 } }" << '\n';
   }
   ofil << '\n' << "global_settings { ambient_light White }" << '\n';
   ofil << '\n' << "light_source {" << '\n';
   currIndLev++; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << thetabs;
   HelpersPOVRay::WriteVector(ofil,lightSource[0][0],lightSource[0][1],lightSource[0][2]);
   ofil << '\n' << thetabs << "color White*0.5" << '\n';
   currIndLev--; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << '}' << '\n';
   for (int i=1; i<nLightSources; i++) {
      ofil << '\n' << "light_source {" << '\n';
      currIndLev++; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
      ofil << thetabs;
      HelpersPOVRay::WriteVector(ofil,lightSource[i][0],lightSource[i][1],lightSource[i][2]);
      ofil << '\n' << thetabs << "color White*0.5" << '\n';
      currIndLev--; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
      ofil << '}' << '\n';
   }
   ofil << '\n' << "background { color ";
   HelpersPOVRay::WriteVector(ofil,backGroundColor[0],backGroundColor[1],backGroundColor[2]);
   ofil << " }" << '\n';
   if (placecam) {WritePlaceCamera(ofil);}
   return false;
}
void POVRayConfiguration::SetBGColor(double rgbr,double rgbg,double rgbb) {
   backGroundColor[0]=rgbr; backGroundColor[1]=rgbg; backGroundColor[2]=rgbb;
   return;
}
void POVRayConfiguration::SetAngView(double av) {
   angView=av;
   return;
}
void POVRayConfiguration::SetAngView(double avx,double avy,double avz) {
   vecAngView[0]=avx; vecAngView[1]=avy; vecAngView[2]=avz;
}
void POVRayConfiguration::SetSkyCam(double xx,double yy,double zz) {
   skyCam[0]=xx; skyCam[1]=yy; skyCam[2]=zz;
   return;
}
void POVRayConfiguration::SetLocCam(double xx,double yy,double zz) {
   locCam[0]=xx; locCam[1]=yy; locCam[2]=zz;
   return;
}
void POVRayConfiguration::SetLookAtCam(double xx,double yy,double zz) {
   lookAtCam[0]=xx; lookAtCam[1]=yy; lookAtCam[2]=zz;
   return;
}
void POVRayConfiguration::AddLightSource(double xx,double yy,double zz) {
   double **cpls;
   int cpnls;
   MyMemory::Alloc2DRealArray("cpls",nLightSources,3,cpls);
   for (int i=0; i<nLightSources; i++) {
      for (int j=0; j<3; j++) {
         cpls[i][j]=lightSource[i][j];
      }
   }
   cpnls=nLightSources;
   MyMemory::Dealloc2DRealArray(lightSource,nLightSources);
   nLightSources=cpnls+1;
   MyMemory::Alloc2DRealArray(string("lightSource"),nLightSources,3,lightSource);
   for (int i=0; i<cpnls; i++) {
      for (int j=0; j<3; j++) {
         lightSource[i][j]=cpls[i][j];
      }
   }
   lightSource[cpnls][0]=xx; lightSource[cpnls][1]=yy; lightSource[cpnls][2]=zz;
   MyMemory::Dealloc2DRealArray(cpls,cpnls);
   cpnls=0;
   cpls=NULL;
   return;
}
void POVRayConfiguration::WriteIncColors(ofstream &ofil) {
   if (inccolors) {
      ofil << "//To include colors, textures, etc..." << '\n';
      ofil << "#include \"colors.inc\"" << '\n';
   }
   return;
}
void POVRayConfiguration::WritePlaceCamera(ofstream &ofil) {
   string thetabs;
   ofil << '\n' << "//Place the camera:" << '\n' << "camera{" << '\n';
   currIndLev++;
   thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << thetabs << "sky ";
   HelpersPOVRay::WriteVector(ofil,skyCam[0],skyCam[1],skyCam[2]);
   ofil << '\n' << thetabs << "location ";
   HelpersPOVRay::WriteVector(ofil,locCam[0],locCam[1],locCam[2]);
   ofil << '\n' << thetabs << "look_at ";
   HelpersPOVRay::WriteVector(ofil,lookAtCam[0],lookAtCam[1],lookAtCam[2]);
   ofil << '\n';
   currIndLev--;
   ofil << '}' << '\n';
   return;
}
void POVRayConfiguration::WriteLightSource(ofstream &ofil,int is,double intens,const char* opts) {
   string thetabs;
   ofil << '\n' << "light_source {" << '\n';
   currIndLev++; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << thetabs;
   HelpersPOVRay::WriteVector(ofil,lightSource[is][0],lightSource[is][1],lightSource[is][2]);
   ofil << '\n' << thetabs << "color White*" << intens << '\n' << thetabs << opts << '\n';
   currIndLev--; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << '}' << '\n';
}
void POVRayConfiguration::ScaleLightSources(double scfactor) {
   for ( int i=0 ; i<nLightSources ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {lightSource[i][k]*=scfactor;}
   }
}
void HelpersPOVRay::WriteIndTabs(ofstream &ofil, int nt) {
   for (int i=0; i<nt; i++) {
      for (int j=0; j<SPTSIZEINDENT; j++) {ofil << ' ';}
   }
}
/*
void HelpersPOVRay::WriteVector(ofstream &ofil,double xx,double yy, double zz) {
   ofil << "< " << xx << ", " << yy << ", " << zz << " >";
}
void WriteVector(ofstream &ofil,int xx,int yy, int zz) {
   ofil << "< " << xx << ", " << yy << ", " << zz << " >";
}
// */
string HelpersPOVRay::IndTabsStr(int nt) {
   string s="";
   for (int i=0; i<nt; i++) {s.append("  ");}
   return s;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr,
                    double cr,double cg, double cb) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << '\n';
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << '}' << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx, double yy, double zz, double rr,
                    const char * str) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << '\n';
   ofil << thetabs << "pigment { " << str << " }" << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx, double yy, double zz, const char* strrad,
                    const char * strcol) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << strrad << '\n';
   ofil << thetabs << "pigment { " << strcol << " }" << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteTransparentSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr,
                    double cr,double cg, double cb,double trc) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << '\n';
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << "transmit "  << trc << '}' << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteTransparentSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr,
                    double cr,double cg, double cb,string trnsmStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << '\n';
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << "transmit "  << trnsmStr << '}' << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteTriangle(ofstream &ofil,int nt) {
   if (nt==0) {
      //std::cout << "The function WriteTriangle is under construction, nothing done..." << '\n';
      ofil << "//Attempt to use function WriteTriangle..." << '\n';
   }
   return false;
}
bool HelpersPOVRay::WriteCylinder(ofstream &ofil,int nt, 
                      double xa, double ya, double za, 
                      double xb, double yb, double zb, double rr,
                      double cr, double cg, double cb,string pigmentStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "cylinder { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xa,ya,za);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,xb,yb,zb);
   ofil << ", " << rr << '\n';
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << " " << pigmentStr << " " << '\n';
   ofil << " }" << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteCylinder(ofstream &ofil,int nt,
                      double xa, double ya, double za,
                      double xb, double yb, double zb, double rr,
                      const char* str) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "cylinder { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xa,ya,za);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,xb,yb,zb);
   ofil << ", " << rr << '\n';
   ofil << thetabs << "pigment { " << str << " }" << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteCylinder(ofstream &ofil,int nt, 
                      double xa, double ya, double za, 
                      double xb, double yb, double zb, double rr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "cylinder { " << '\n';
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xa,ya,za);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,xb,yb,zb);
   ofil << ", " << rr << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteArrow(ofstream &ofil,int nt) {
   if (nt==0) {
      //std::cout << "The function WriteTriangle is under construction, nothing done..." << '\n';
      ofil << "//Attempt to use function WriteArrow..." << '\n';
   }
   return false;
}
bool HelpersPOVRay::WriteTriangle(ofstream &ofil,int nt,\
      double x1,double y1,double z1,double x2,double y2,double z2,
      double x3,double y3,double z3,double cr,double cg,double cb,
      const string &pigmentStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "triangle { " << '\n';
   ++indlev; thetabs=IndTabsStr(indlev);
   ofil << scientific << setprecision(8);
   ofil << thetabs;
   WriteVector(ofil,x1,y1,z1);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,x2,y2,z2);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,x3,y3,z3);
   ofil << '\n' << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << (pigmentStr.size()>0 ? (string(" ")+pigmentStr+string(" ")) : " ") << '}' << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSmoothTriangle(ofstream &ofil,int nt,\
         const vector<vector<double> > &v,const vector<vector<double> > &n,
         double cr,double cg,double cb,const string &pigmentStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "smooth_triangle { " << '\n';
   ++indlev; thetabs=IndTabsStr(indlev);
   ofil << scientific << setprecision(12);
   ofil << thetabs;
   WriteVector(ofil,v[0][0],v[0][1],v[0][2]);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,n[0][0],n[0][1],n[0][2]);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,v[1][0],v[1][1],v[1][2]);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,n[1][0],n[1][1],n[1][2]);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,v[2][0],v[2][1],v[2][2]);
   ofil << "," << '\n' << thetabs;
   WriteVector(ofil,n[2][0],n[2][1],n[2][2]);
   ofil << '\n' << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << (pigmentStr.size()>0 ? (string(" ")+pigmentStr+string(" ")) : " ") << '}' << '\n';
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << '}' << '\n';
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteMesh2SingleRGB(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<size_t> > &f,const int nt,vector<double> rgb,\
         const string &trnsmStr) {
   vector<vector<double> > n(0);
   return WriteMesh2SingleRGB(ofil,v,n,f,nt,rgb,trnsmStr);
}
bool HelpersPOVRay::WriteMesh2WithTextures(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<double> > &t,const vector<vector<size_t> > &f,\
         const int nt,const string &trnsmStr) {
   vector<vector<double> > n(0);
   return WriteMesh2WithTextures(ofil,v,n,t,f,nt,trnsmStr);
}
bool HelpersPOVRay::WriteMesh2SingleRGB(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<double> > &n,const vector<vector<size_t> > &f,\
         const int nt,vector<double> rgb,const string &trnsmStr) {
   int indlev=nt;
   size_t nvm1=v.size()-1;
   size_t nfm1=f.size()-1;
   bool writeTransmit=(trnsmStr.size()>0);
   string thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "mesh2 {" << '\n';
   thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "vertex_vectors {\n";
   thetabs=IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      WriteVector(ofil,v[i][0],v[i][1],v[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   WriteVector(ofil,v[nvm1][0],v[nvm1][1],v[nvm1][2]);
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << "}\n";//end of vertex_vectors
   if ( n.size()>0 ) {
      ofil << thetabs << "normal_vectors {\n";
      thetabs=IndTabsStr(indlev);
      ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
      for ( size_t i=0 ; i<nvm1 ; ++i ) {
         WriteVector(ofil,n[i][0],n[i][1],n[i][2]);
         ofil << ",";
         if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
      }
      WriteVector(ofil,n[nvm1][0],n[nvm1][1],n[nvm1][2]);
      thetabs=IndTabsStr(--indlev);
      ofil << thetabs << '}' << '\n';//end of normal_vectors
   }
   thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "face_indices {\n";
   thetabs=IndTabsStr(indlev);
   ofil << thetabs << (nfm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nfm1 ; ++i ) {
      WriteVector(ofil,f[i][0],f[i][1],f[i][2]);
      ofil << ",";
      if ( (i%5) == 4 ) { ofil << '\n' << thetabs; }
   }
   WriteVector(ofil,f[nfm1][0],f[nfm1][1],f[nfm1][2]);
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << '}' << '\n';//end of face_indices
   thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,rgb[0],rgb[1],rgb[2]);
   if ( writeTransmit ) { ofil << ' ' << trnsmStr << ' '; }
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << '}' << '\n';//end of pigment
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << '}' << '\n'; //end of mesh2
   return true;
}
bool HelpersPOVRay::WriteMesh2WithTextures(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<double> > &n,const vector<vector<double> > &t,const vector<vector<size_t> > &f,\
         const int nt,const string &trnsmStr) {
   int indlev=nt;
   size_t nvm1=v.size()-1;
   size_t nfm1=f.size()-1;
   bool writeTransmit=(trnsmStr.size()>0);
   string thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "mesh2 {" << '\n';
   thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "vertex_vectors {\n";
   thetabs=IndTabsStr(indlev);
   ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nvm1 ; ++i ) {
      WriteVector(ofil,v[i][0],v[i][1],v[i][2]);
      ofil << ",";
      if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
   }
   WriteVector(ofil,v[nvm1][0],v[nvm1][1],v[nvm1][2]);
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << "}\n";//end of vertex_vectors
   if ( n.size()>0 ) {
      ofil << thetabs << "normal_vectors {\n";
      thetabs=IndTabsStr(indlev);
      ofil << thetabs << (nvm1+1) << ",\n" << thetabs;
      for ( size_t i=0 ; i<nvm1 ; ++i ) {
         WriteVector(ofil,n[i][0],n[i][1],n[i][2]);
         ofil << ",";
         if ( (i%3) == 2 ) { ofil << '\n' << thetabs; }
      }
      WriteVector(ofil,n[nvm1][0],n[nvm1][1],n[nvm1][2]);
      thetabs=IndTabsStr(--indlev);
      ofil << thetabs << '}' << '\n';//end of normal_vectors
   }
   ofil << thetabs << "texture_list {\n";
   thetabs=IndTabsStr(indlev);
   size_t nnvv=nvm1+1;
   ofil << thetabs << (nnvv) << ",\n";
   for ( size_t i=0 ; i<nnvv ; ++i ) {
      ofil << thetabs << "texture{pigment{rgb ";
      HelpersPOVRay::WriteVector(ofil,t[i][0],t[i][1],t[i][2]);
      if ( writeTransmit ) { ofil << ' ' << trnsmStr; }
      ofil << " }}" << (i<nvm1? ',' : ' ') << "\n";
   }
   thetabs=HelpersPOVRay::IndTabsStr(--indlev);
   ofil << thetabs << "}" << endl;//end of texture_list
   thetabs=IndTabsStr(indlev++);
   ofil << thetabs << "face_indices {\n";
   thetabs=IndTabsStr(indlev);
   ofil << thetabs << (nfm1+1) << ",\n" << thetabs;
   for ( size_t i=0 ; i<nfm1 ; ++i ) {
      WriteVector(ofil,f[i][0],f[i][1],f[i][2]);
      ofil << ", " << f[i][0] << ',' << f[i][1] << ',' << f[i][2] << ", ";
      if ( (i%5) == 4 ) { ofil << '\n' << thetabs; }
   }
   WriteVector(ofil,f[nfm1][0],f[nfm1][1],f[nfm1][2]);
   ofil << ", " <<  f[nfm1][0] << ',' << f[nfm1][1] << ',' << f[nfm1][2] << '\n';
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << '}' << '\n';//end of face_indices
   thetabs=IndTabsStr(indlev++);
   thetabs=IndTabsStr(--indlev);
   ofil << thetabs << '}' << '\n'; //end of mesh2
   return true;
}








