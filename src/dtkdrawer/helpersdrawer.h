/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
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
#ifndef _HELPERSDRAWER_H_
#define _HELPERSDRAWER_H_
#include <string>
using std::string;
#include "optflags.h"
#include "bondnetwork.h"

/* ************************************************************************** */
class HelpersDrawer {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static bool MakePovFile(const string &basename,const OptionFlags &options,\
         char *argv[],BondNetWork &bn);
   static bool MakeCriticalPointsPovFile(const string &basename,const OptionFlags &options,\
         char *argv[],BondNetWork &bn,const string cpname);
   static bool AlignMolecule(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]);
   static bool AlignMolecule3Atoms(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char *argv[]);
   static void SetupPovConf(POVRayConfiguration &pvp,BondNetWork &bn,\
         const OptionFlags &options,char* argv[]);
   static void WriteCameraSettings(ofstream &ofil,POVRayConfiguration &pvc,\
         bool lookat=false);
   static void WriteLightSources(ofstream &ofil,POVRayConfiguration &pvc);
   static void WritePovHeader(ofstream &ofil);
   static void WriteBondNetworkElements(ofstream &ofil,BondNetWork &bn,\
         POVRayConfiguration &pvc,const OptionFlags &options,char *argv[]);
   static bool DrawCPofType(char t,const OptionFlags &options,char *argv[]);
   static bool DrawGPofType(char t,const OptionFlags &options,char *argv[]);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSDRAWER_H_ */

