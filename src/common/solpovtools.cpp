/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.1.2
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */



#ifndef _SOLPOVTOOLS_CPP_
#define _SOLPOVTOOLS_CPP_

#include "solpovtools.h"
#include "solmemhand.h"
#include <iomanip>
using std::setprecision;
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
povRayConfProp::povRayConfProp()
{
   angView=67.5e0;
   skyCam[0]=0.0e0; skyCam[1]=0.0e0; skyCam[2]=1.0e0;
   locCam[0]=0.0e0; locCam[1]=0.0e0; locCam[2]=2.3e0;
   lookAtCam[0]=0.0e0; lookAtCam[1]=0.0e0; lookAtCam[2]=0.0e0;
   nLightSources=2;
   alloc2DRealArray("lightSource",nLightSources,3,lightSource);
   defLightSource[0]=1.0; defLightSource[1]=1.0; defLightSource[2]=1.0;
   for (int i=0; i<3; i++) {
      lightSource[0][i]=defLightSource[i];
      lightSource[1][i]=defLightSource[i];
   }
   lightSource[1][2]*=(-1.0);
   backGroundColor[0]=0.0; backGroundColor[1]=0.5; backGroundColor[2]=0.7;
   currIndLev=0;
   shine=false;
   inccolors=true;
}
//**********************************************************************************************
povRayConfProp::~povRayConfProp()
{
   dealloc2DRealArray(lightSource,nLightSources);
}
//**********************************************************************************************
bool povRayConfProp::writeHeader(ofstream &ofil,bool placecam)
{
   string thetabs;
   writeIncColors(ofil);
   if (shine) {
      ofil << endl << "#default { finish { specular 0.3 roughness 0.03 phong .1 } }" << endl;
   }
   ofil << endl << "global_settings { ambient_light White }" << endl;
   ofil << endl << "light_source {" << endl;
   currIndLev++; thetabs=indTabsStr(currIndLev);
   ofil << thetabs;
   writePoVVector(ofil,lightSource[0][0],lightSource[0][1],lightSource[0][2]);
   ofil << endl << thetabs << "color White*0.5" << endl;
   currIndLev--; thetabs=indTabsStr(currIndLev);
   ofil << "}" << endl;
   for (int i=1; i<nLightSources; i++) {
      ofil << endl << "light_source {" << endl;
      currIndLev++; thetabs=indTabsStr(currIndLev);
      ofil << thetabs;
      writePoVVector(ofil,lightSource[i][0],lightSource[i][1],lightSource[i][2]);
      ofil << endl << thetabs << "color White*0.5" << endl;
      currIndLev--; thetabs=indTabsStr(currIndLev);
      ofil << "}" << endl;
   }
   ofil << endl << "background { color ";
   writePoVVector(ofil,backGroundColor[0],backGroundColor[1],backGroundColor[2]);
   ofil << " }" << endl;
   if (placecam) {writePlaceCamera(ofil);}
   return false;
}
//**********************************************************************************************
void povRayConfProp::setBGColor(solreal rgbr,solreal rgbg,solreal rgbb)
{
   backGroundColor[0]=rgbr; backGroundColor[1]=rgbg; backGroundColor[2]=rgbb;
   return;
}
//**********************************************************************************************
void povRayConfProp::setAngView(solreal av)
{
   angView=av;
   return;
}
//**********************************************************************************************
void povRayConfProp::setSkyCam(solreal xx,solreal yy,solreal zz)
{
   skyCam[0]=xx; skyCam[1]=yy; skyCam[2]=zz;
   return;
}
//**********************************************************************************************
void povRayConfProp::setLocCam(solreal xx,solreal yy,solreal zz)
{
   locCam[0]=xx; locCam[1]=yy; locCam[2]=zz;
   return;
}
//**********************************************************************************************
void povRayConfProp::setLookAtCam(solreal xx,solreal yy,solreal zz)
{
   lookAtCam[0]=xx; lookAtCam[1]=yy; lookAtCam[2]=zz;
   return;
}
//**********************************************************************************************
void povRayConfProp::addLightSource(solreal xx,solreal yy,solreal zz)
{
   solreal **cpls;
   int cpnls;
   alloc2DRealArray("cpls",nLightSources,3,cpls);
   for (int i=0; i<nLightSources; i++) {
      for (int j=0; j<3; j++) {
         cpls[i][j]=lightSource[i][j];
      }
   }
   cpnls=nLightSources;
   dealloc2DRealArray(lightSource,nLightSources);
   nLightSources=cpnls+1;
   alloc2DRealArray(string("lightSource"),nLightSources,3,lightSource);
   for (int i=0; i<cpnls; i++) {
      for (int j=0; j<3; j++) {
         lightSource[i][j]=cpls[i][j];
      }
   }
   lightSource[cpnls][0]=xx; lightSource[cpnls][1]=yy; lightSource[cpnls][2]=zz;
   dealloc2DRealArray(cpls,cpnls);
   cpnls=0;
   cpls=NULL;
   return;
}
//**************************************************************************************************
void povRayConfProp::writeIncColors(ofstream &ofil)
{
   if (inccolors) {
      ofil << "//To include colors, textures, etc..." << endl;
      ofil << "#include \"colors.inc\"" << endl;
   }
   return;
}
//**************************************************************************************************
void povRayConfProp::writePlaceCamera(ofstream &ofil)
{
   string thetabs;
   ofil << endl << "//Place the camera:" << endl << "camera{" << endl;
   currIndLev++;
   thetabs=indTabsStr(currIndLev);
   ofil << thetabs << "sky ";
   writePoVVector(ofil,skyCam[0],skyCam[1],skyCam[2]);
   ofil << endl << thetabs << "location ";
   writePoVVector(ofil,locCam[0],locCam[1],locCam[2]);
   ofil << endl << thetabs << "look_at ";
   writePoVVector(ofil,lookAtCam[0],lookAtCam[1],lookAtCam[2]);
   ofil << endl;
   currIndLev--;
   ofil << "}" << endl;
   return;
}
//**************************************************************************************************
//*
void povRayConfProp::writeLightSource(ofstream &ofil,int is,solreal intens,const char* opts)
{
   string thetabs;
   ofil << endl << "light_source {" << endl;
   currIndLev++; thetabs=indTabsStr(currIndLev);
   ofil << thetabs;
   writePoVVector(ofil,lightSource[is][0],lightSource[is][1],lightSource[is][2]);
   ofil << endl << thetabs << "color White*" << intens << endl << thetabs << opts << endl;
   currIndLev--; thetabs=indTabsStr(currIndLev);
   ofil << "}" << endl;
}
// */
//**************************************************************************************************
void povRayConfProp::scaleLightSources(solreal scfactor)
{
   for ( int i=0 ; i<nLightSources ; i++ ) {
      for ( int k=0 ; k<3 ; k++ ) {lightSource[i][k]*=scfactor;}
   }
}
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**************************************************************************************************
//**********************************************************************************************
void writeIndTabs(ofstream &ofil, int nt)
{
   for (int i=0; i<nt; i++) {
      for (int j=0; j<SPTSIZEINDENT; j++) {ofil << ' ';}
   }
}
//**********************************************************************************************
void writePoVVector(ofstream &ofil,solreal xx,solreal yy, solreal zz)
{
   ofil << "< " << xx << ", " << yy << ", " << zz << " >";
}
//**********************************************************************************************
void writePoVVector(ofstream &ofil,int xx,int yy, int zz)
{
   ofil << "< " << xx << ", " << yy << ", " << zz << " >";
}
//**********************************************************************************************
string indTabsStr(int nt)
{
   string s="";
   for (int i=0; i<nt; i++) {s.append("  ");}
   return s;
}
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx,solreal yy,solreal zz, solreal rr,
                    solreal cr,solreal cg, solreal cb)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   writePoVVector(ofil,cr,cg,cb);
   ofil << "}" << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx, solreal yy, solreal zz, solreal rr,
                    const char * str)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { " << str << " }" << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx, solreal yy, solreal zz, const char* strrad,
                    const char * strcol)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xx,yy,zz);
   ofil << ", " << strrad << endl;
   ofil << thetabs << "pigment { " << strcol << " }" << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx,solreal yy,solreal zz, solreal rr)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVTransparentSphere(ofstream &ofil,int nt,
                    solreal xx,solreal yy,solreal zz, solreal rr,
                    solreal cr,solreal cg, solreal cb,solreal trc)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "sphere { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xx,yy,zz);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   writePoVVector(ofil,cr,cg,cb);
   ofil << "transmit "  << trc << "}" << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVTriangle(ofstream &ofil,int nt)
{
   if (nt==0) {
      //std::cout << "The function writePOVTriangle is under construction, nothing done..." << endl;
      ofil << "//Attempt to use function writePOVTriangle..." << endl;
   }
   return false;
}
//**********************************************************************************************
bool writePOVCylinder(ofstream &ofil,int nt, 
                      solreal xa, solreal ya, solreal za, 
                      solreal xb, solreal yb, solreal zb, solreal rr,
                      solreal cr, solreal cg, solreal cb)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "cylinder { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xa,ya,za);
   ofil << "," << endl << thetabs;
   writePoVVector(ofil,xb,yb,zb);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { rgb ";
   writePoVVector(ofil,cr,cg,cb);
   ofil << " }" << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVCylinder(ofstream &ofil,int nt,
                      solreal xa, solreal ya, solreal za,
                      solreal xb, solreal yb, solreal zb, solreal rr,
                      const char* str)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "cylinder { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xa,ya,za);
   ofil << "," << endl << thetabs;
   writePoVVector(ofil,xb,yb,zb);
   ofil << ", " << rr << endl;
   ofil << thetabs << "pigment { " << str << " }" << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVCylinder(ofstream &ofil,int nt, 
                      solreal xa, solreal ya, solreal za, 
                      solreal xb, solreal yb, solreal zb, solreal rr)
{
   int indlev=nt;
   string thetabs=indTabsStr(indlev);
   ofil << thetabs << "cylinder { " << endl;
   indlev++; thetabs=indTabsStr(indlev);
   ofil << thetabs;
   ofil << scientific << setprecision(8);
   writePoVVector(ofil,xa,ya,za);
   ofil << "," << endl << thetabs;
   writePoVVector(ofil,xb,yb,zb);
   ofil << ", " << rr << endl;
   indlev--; thetabs=indTabsStr(indlev);
   ofil << thetabs << "}" << endl;
   ofil.unsetf(ios::scientific);
   return true;
}
//**********************************************************************************************
bool writePOVArrow(ofstream &ofil,int nt)
{
   if (nt==0) {
      //std::cout << "The function writePOVTriangle is under construction, nothing done..." << endl;
      ofil << "//Attempt to use function writePOVArrow..." << endl;
   }
   return false;
}
//**************************************************************************************************
//**************************************************************************************************
#endif//_SOLPOVTOOLS_CPP_

