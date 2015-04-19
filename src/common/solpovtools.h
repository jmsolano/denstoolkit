/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.1
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



#ifndef _SOLPOVTOOLS_H_
#define _SOLPOVTOOLS_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#include <string>
using namespace std;
#include <fstream>
using std::ifstream;
using std::ofstream;

#ifndef SPTSIZEINDENT
#define SPTSIZEINDENT 2
#endif

//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
class povRayConfProp
//**********************************************************************************************
{
public:
   //***********************************************************************************************
   int nLightSources;
   solreal **lightSource;
   solreal backGroundColor[3];
   solreal angView;
   solreal skyCam[3];
   solreal locCam[3];
   solreal lookAtCam[3];
   int currIndLev;
   solreal defLightSource[3];
   bool shine,inccolors;
   //**********************************************************************************************
   povRayConfProp();
   //**********************************************************************************************
   ~povRayConfProp();
   //**********************************************************************************************
   bool writeHeader(ofstream &ofil,bool placecam = true);
   //**********************************************************************************************
   void setBGColor(solreal rgbr,solreal rgbg,solreal rgbb);
   //**********************************************************************************************
   void setAngView(solreal av);
   //**********************************************************************************************
   void setSkyCam(solreal xx,solreal yy,solreal zz);
   //**********************************************************************************************
   void setLocCam(solreal xx,solreal yy,solreal zz);
   //**********************************************************************************************
   void setLookAtCam(solreal xx,solreal yy,solreal zz);
   //**********************************************************************************************
   void addLightSource(solreal xx,solreal yy,solreal zz);
   //**********************************************************************************************
   void writeIncColors(ofstream &ofil);
   //**********************************************************************************************
   void writePlaceCamera(ofstream &ofil);
   //**********************************************************************************************
   void writeLightSource(ofstream &ofil,int is,solreal intens,const char* opts);
   //**********************************************************************************************
   void scaleLightSources(solreal scfactor);
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
   //**********************************************************************************************
};
//**********************************************************************************************
void writeIndTabs(ofstream &ofil, int nt);
//**********************************************************************************************
void writePoVVector(ofstream &ofil,solreal xx,solreal yy, solreal zz);
//**********************************************************************************************
void writePoVVector(ofstream &ofil,int xx, int yy, int zz);
//**********************************************************************************************
string indTabsStr(int nt);
//**********************************************************************************************
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                 solreal xx, solreal yy, solreal zz, solreal rr,
                 solreal cr, solreal cg, solreal cb);
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx, solreal yy, solreal zz, solreal rr,
                    const char *str);
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx, solreal yy, solreal zz, const char* strrad,
                    const char *strcol);
//**********************************************************************************************
bool writePOVSphere(ofstream &ofil,int nt,
                    solreal xx,solreal yy,solreal zz,solreal rr);
//**********************************************************************************************
bool writePOVTransparentSphere(ofstream &ofil,int nt,
                    solreal xx, solreal yy, solreal zz, solreal rr,
                    solreal cr, solreal cg, solreal cb,solreal trc);
//**********************************************************************************************
bool writePOVTriangle(ofstream &ofil,int nt);
//**********************************************************************************************
bool writePOVCylinder(ofstream &ofil,int nt, 
                      solreal xa, solreal ya, solreal za, 
                      solreal xb, solreal yb, solreal zb, solreal rr,
                      solreal cr, solreal cg, solreal cb);
//**********************************************************************************************
bool writePOVCylinder(ofstream &ofil,int nt,
                      solreal xa, solreal ya, solreal za,
                      solreal xb, solreal yb, solreal zb, solreal rr,
                      const char * str);
//**********************************************************************************************
bool writePOVCylinder(ofstream &ofil,int nt,
                      solreal xa, solreal ya, solreal za,
                      solreal xb, solreal yb, solreal zb, solreal rr);
//**********************************************************************************************
bool writePOVArrow(ofstream &ofil,int nt);
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
#endif//_SOLPOVTOOLS_H_



