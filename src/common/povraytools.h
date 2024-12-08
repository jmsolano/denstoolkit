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
   void SetAngView(double avx,double avy,double avz);
   void SetSkyCam(double xx,double yy,double zz);
   void SetLocCam(double xx,double yy,double zz);
   void SetLookAtCam(double xx,double yy,double zz);
   void AddLightSource(double xx,double yy,double zz);
   void WriteIncColors(ofstream &ofil);
   void WritePlaceCamera(ofstream &ofil);
   void WriteLightSource(ofstream &ofil,int is,double intens,const char* opts);
   void ScaleLightSources(double scfactor);
   void UnsetVersion36() {setversion36=false;};
   void SetupUsualLightSources(const double rv);
   void SelectStandardCameraVectors(const double zsep=1.0e0);
   /** This function multiply locCam, up, right, dir, and lightSources by M.
    * Hence, this function should be called AFTER setting these vectors.  */
   void ApplyRotationMatrixToCameraAndLightSources(const vector<vector<double> > &M);
/* ************************************************************************** */
   int nLightSources;
   double **lightSource;
   double backGroundColor[3];
   double angView;
   double vecAngView[3];
   double skyCam[3];
   double locCam[3];
   double lookAtCam[3];
   double vecUp[3];
   double vecRight[3];
   double vecDir[3];
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
   static bool WriteMesh2SingleRGB(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<double> > &n,const vector<vector<size_t> > &f,\
         const int nt,vector<double> rgb,const string &trnsmStr="");
   /** Writes a mesh2 povray object with (v)ertices, (n)normals, (t)extures,
    * and (f)aces.  */
   static bool WriteMesh2WithTextures(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<double> > &n,const vector<vector<double> > &t,const vector<vector<size_t> > &f,\
         const int nt,const string &trnsmStr="");
   static bool WriteMesh2WithTextures(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<double> > &t,const vector<vector<size_t> > &f,\
         const int nt,const string &trnsmStr="");
   static bool WriteMesh2SingleRGB(ofstream &ofil,const vector<vector<double> > &v,\
         const vector<vector<size_t> > &f,const int nt,vector<double> rgb,\
         const string &trnsmStr="");
   static void ApplyRotationMatrix(const vector<vector<double> > &M,double*v);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif//_POVRAYTOOLS_H_



