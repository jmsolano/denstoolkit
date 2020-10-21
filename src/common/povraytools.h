#ifndef _POVRAYTOOLS_H_
#define _POVRAYTOOLS_H_
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>
using std::string;
#include <vector>
using std::vector;

#ifndef SPTSIZEINDENT
#define SPTSIZEINDENT 2
#endif
/* ************************************************************************** */
class POVRayConfiguration {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   POVRayConfiguration();
   ~POVRayConfiguration();
   bool WriteHeader(ofstream &ofil,bool placecam = true);
   void SetBGColor(double rgbr,double rgbg,double rgbb);
   void SetAngView(double av);
   void SetSkyCam(double xx,double yy,double zz);
   void SetLocCam(double xx,double yy,double zz);
   void SetLookAtCam(double xx,double yy,double zz);
   void AddLightSource(double xx,double yy,double zz);
   void WriteIncColors(ofstream &ofil);
   void WritePlaceCamera(ofstream &ofil);
   void WriteLightSource(ofstream &ofil,int is,double intens,const char* opts);
   void ScaleLightSources(double scfactor);
   void UnsetVersion36() {setversion36=false;};
/* ************************************************************************** */
   int nLightSources;
   double **lightSource;
   double backGroundColor[3];
   double angView;
   double skyCam[3];
   double locCam[3];
   double lookAtCam[3];
   int currIndLev;
   double defLightSource[3];
   bool shine,inccolors,setversion36;
/* ************************************************************************** */
};

/* ************************************************************************** */
class HelpersPOVRay {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static void WriteIndTabs(ofstream &ofil, int nt);
   template <class Tp> static inline void WriteVector(ofstream &ofil,Tp xx,Tp yy,Tp zz)
   {ofil << "< " << xx << ", " << yy << ", " << zz << " >";}
   //static void WriteVector(ofstream &ofil,double xx,double yy, double zz);
   //static void WriteVector(ofstream &ofil,int xx,int yy,int zz);
   //static void WriteVector(ofstream &ofil,size_t xx,size_t yy,size_t zz);
   static string IndTabsStr(int nt);
   static bool WriteSphere(ofstream &ofil,int nt,
         double xx, double yy, double zz, double rr,
         double cr, double cg, double cb);
   static bool WriteSphere(ofstream &ofil,int nt,
         double xx, double yy, double zz, double rr,
         const char *str);
   static bool WriteSphere(ofstream &ofil,int nt,
         double xx, double yy, double zz, const char* strrad,
         const char *strcol);
   static bool WriteSphere(ofstream &ofil,int nt,
         double xx,double yy,double zz,double rr);
   static bool WriteTransparentSphere(ofstream &ofil,int nt,
         double xx, double yy, double zz, double rr,
         double cr, double cg, double cb,double trc);
   static bool WriteTransparentSphere(ofstream &ofil,int nt,
         double xx, double yy, double zz, double rr,
         double cr, double cg, double cb,string trnsmStr);
   static bool WriteTriangle(ofstream &ofil,int nt);
   static bool WriteCylinder(ofstream &ofil,int nt, 
         double xa, double ya, double za, 
         double xb, double yb, double zb, double rr,
         double cr, double cg, double cb,string pigmentStr=string(""));
   static bool WriteCylinder(ofstream &ofil,int nt,
         double xa, double ya, double za,
         double xb, double yb, double zb, double rr,
         const char * str);
   static bool WriteCylinder(ofstream &ofil,int nt,
         double xa, double ya, double za,
         double xb, double yb, double zb, double rr);
   static bool WriteArrow(ofstream &ofil,int nt);
   static bool WriteTriangle(ofstream &ofil,int nt,\
         double x1,double y1,double z1,
         double x2,double y2,double z2,
         double x3,double y3,double z3,
         double cr,double cg,double cb,const string &pigmentStr="");
   static bool WriteSmoothTriangle(ofstream &ofil,int nt,\
         const vector<vector<double> > &v,const vector<vector<double> > &n,
         double cr,double cg,double cb,const string &pigmentStr="");
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif//_POVRAYTOOLS_H_



