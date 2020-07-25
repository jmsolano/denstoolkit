#ifndef _POVRAYTOOLS_CPP_
#define _POVRAYTOOLS_CPP_
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
   skyCam[0]=0.0e0; skyCam[1]=0.0e0; skyCam[2]=1.0e0;
   locCam[0]=0.0e0; locCam[1]=0.0e0; locCam[2]=2.3e0;
   lookAtCam[0]=0.0e0; lookAtCam[1]=0.0e0; lookAtCam[2]=0.0e0;
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
           << "this line..." << endl;
   }
   WriteIncColors(ofil);
   if (shine) {
      ofil << endl << "#default { finish { specular 0.3 roughness 0.03 phong .1 } }" << endl;
   }
   ofil << endl << "global_settings { ambient_light White }" << endl;
   ofil << endl << "light_source {" << endl;
   currIndLev++; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << thetabs;
   HelpersPOVRay::WriteVector(ofil,lightSource[0][0],lightSource[0][1],lightSource[0][2]);
   ofil << endl << thetabs << "color White*0.5" << endl;
   currIndLev--; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << "}" << endl;
   for (int i=1; i<nLightSources; i++) {
      ofil << endl << "light_source {" << endl;
      currIndLev++; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
      ofil << thetabs;
      HelpersPOVRay::WriteVector(ofil,lightSource[i][0],lightSource[i][1],lightSource[i][2]);
      ofil << endl << thetabs << "color White*0.5" << endl;
      currIndLev--; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
      ofil << "}" << endl;
   }
   ofil << endl << "background { color ";
   HelpersPOVRay::WriteVector(ofil,backGroundColor[0],backGroundColor[1],backGroundColor[2]);
   ofil << " }" << endl;
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
      ofil << "//To include colors, textures, etc..." << endl;
      ofil << "#include \"colors.inc\"" << endl;
   }
   return;
}
void POVRayConfiguration::WritePlaceCamera(ofstream &ofil) {
   string thetabs;
   ofil << endl << "//Place the camera:" << endl << "camera{" << endl;
   currIndLev++;
   thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << thetabs << "sky ";
   HelpersPOVRay::WriteVector(ofil,skyCam[0],skyCam[1],skyCam[2]);
   ofil << endl << thetabs << "location ";
   HelpersPOVRay::WriteVector(ofil,locCam[0],locCam[1],locCam[2]);
   ofil << endl << thetabs << "look_at ";
   HelpersPOVRay::WriteVector(ofil,lookAtCam[0],lookAtCam[1],lookAtCam[2]);
   ofil << endl;
   currIndLev--;
   ofil << "}" << endl;
   return;
}
void POVRayConfiguration::WriteLightSource(ofstream &ofil,int is,double intens,const char* opts) {
   string thetabs;
   ofil << endl << "light_source {" << endl;
   currIndLev++; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << thetabs;
   HelpersPOVRay::WriteVector(ofil,lightSource[is][0],lightSource[is][1],lightSource[is][2]);
   ofil << endl << thetabs << "color White*" << intens << endl << thetabs << opts << endl;
   currIndLev--; thetabs=HelpersPOVRay::IndTabsStr(currIndLev);
   ofil << "}" << endl;
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
void HelpersPOVRay::WriteVector(ofstream &ofil,double xx,double yy, double zz) {
   ofil << "< " << xx << ", " << yy << ", " << zz << " >";
}
void WriteVector(ofstream &ofil,int xx,int yy, int zz) {
   ofil << "< " << xx << ", " << yy << ", " << zz << " >";
}
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
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << "}" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx, double yy, double zz, double rr,
                    const char * str) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { " << str << " }" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx, double yy, double zz, const char* strrad,
                    const char * strcol) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << strrad << endl;
   ofil << thetabs << "pigment { " << strcol << " }" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteTransparentSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr,
                    double cr,double cg, double cb,double trc) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << "transmit "  << trc << "}" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteTransparentSphere(ofstream &ofil,int nt,
                    double xx,double yy,double zz, double rr,
                    double cr,double cg, double cb,string trnsmStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << "transmit "  << trnsmStr << "}" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteTriangle(ofstream &ofil,int nt) {
   if (nt==0) {
      //std::cout << "The function WriteTriangle is under construction, nothing done..." << endl;
      ofil << "//Attempt to use function WriteTriangle..." << endl;
   }
   return false;
}
bool HelpersPOVRay::WriteCylinder(ofstream &ofil,int nt, 
                      double xa, double ya, double za, 
                      double xb, double yb, double zb, double rr,
                      double cr, double cg, double cb,string pigmentStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "cylinder { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xa,ya,za);
   ofil << "," << endl << thetabs;
   WriteVector(ofil,xb,yb,zb);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << " " << pigmentStr << " " << endl;
   ofil << " }" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteCylinder(ofstream &ofil,int nt,
                      double xa, double ya, double za,
                      double xb, double yb, double zb, double rr,
                      const char* str) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "cylinder { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xa,ya,za);
   ofil << "," << endl << thetabs;
   WriteVector(ofil,xb,yb,zb);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { " << str << " }" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteCylinder(ofstream &ofil,int nt, 
                      double xa, double ya, double za, 
                      double xb, double yb, double zb, double rr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "cylinder { " << endl;
   indlev++; thetabs=IndTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   WriteVector(ofil,xa,ya,za);
   ofil << "," << endl << thetabs;
   WriteVector(ofil,xb,yb,zb);
   ofil << ", " << rr << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
bool HelpersPOVRay::WriteArrow(ofstream &ofil,int nt) {
   if (nt==0) {
      //std::cout << "The function WriteTriangle is under construction, nothing done..." << endl;
      ofil << "//Attempt to use function WriteArrow..." << endl;
   }
   return false;
}
bool HelpersPOVRay::WriteTriangle(ofstream &ofil,int nt,\
      double x1,double y1,double z1,double x2,double y2,double z2,
      double x3,double y3,double z3,double cr,double cg,double cb,
      const string &pigmentStr) {
   int indlev=nt;
   string thetabs=IndTabsStr(indlev);
   ofil << thetabs << "triangle { " << endl;
   ++indlev; thetabs=IndTabsStr(indlev);
   ofil << scientific << setprecision(8);
   ofil << thetabs;
   WriteVector(ofil,x1,y1,z1);
   ofil << "," << endl << thetabs;
   WriteVector(ofil,x2,y2,z2);
   ofil << "," << endl << thetabs;
   WriteVector(ofil,x3,y3,z3);
   ofil << endl << thetabs << "pigment { rgb ";
   WriteVector(ofil,cr,cg,cb);
   ofil << (pigmentStr.size()>0 ? (string(" ")+pigmentStr+string(" ")) : " ") << "}" << endl;
   indlev--; thetabs=IndTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
#endif//_POVRAYTOOLS_CPP_

