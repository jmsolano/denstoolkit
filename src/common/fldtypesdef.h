/*
 *  fldtypesdef.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 10/05/13.
 *  Copyright 2013 . All rights reserved.
 *
 */


#ifndef _FLDTYPESDEF_H_
#define _FLDTYPESDEF_H_
#include <string>
using std::string;
//**********************************************************************************************
enum ScalarFieldType {
   NONE,DENS,MGRD,LAPD,LOLD,ELFD,SENT,KEDK,KEDG,MGLD,GLOL,MEPD,LEDV,MLED,REDG,ROSE
};
//**********************************************************************************************
inline string getFieldTypeKeyShort(const char prop)
{
   string plbl="";
   switch (prop) {
      case 'd':
         plbl="Rho";
         break;
      case 'g':
         plbl="MagGradRho";
         break;
      case 'l':
         plbl="LapRho";
         break;
      case 'E':
         plbl="ELF";
         break;
      case 'L':
         plbl="LOL";
         break;
      case 'M':
         plbl="MagGradLOL";
         break;
      case 'N':
         plbl="GradLOL";
         break;
      case 'p' :
         plbl="LED";
         break;
      case 'P' :
         plbl="MagLED";
         break;
      case 'r' :
         plbl="RoSE";
         break;
      case 's' :
         plbl="RedDensGrad";
         break;
      case 'S':
         plbl="ShannEnt";
         break;
      case 'G':
         plbl="KinEnerDensG";
         break;
      case 'K':
         plbl="KinEnerDensK";
         break;
      case 'V':
         plbl="MEP";
         break;
      default:
         plbl="";
         break;
   }
   return plbl;
}
//**********************************************************************************************
inline string getFieldTypeKeyLong(const char prop)
{
   string plbl="";
   switch (prop) {
      case 'd':
         plbl="Electron Density --Rho--";
         break;
      case 'g':
         plbl="Magnitude of the Gradient of the Electron Density";
         break;
      case 'l':
         plbl="Laplacian of Electron Density";
         break;
      case 'E':
         plbl="Electron Localization Function --ELF--";
         break;
      case 'L':
         plbl="Localized Orbital Locator --LOL--";
         break;
      case 'M':
         plbl="Magnitude of the Gradrient of LOL";
         break;
      case 'N':
         plbl="Gradient of LOL";
         break;
      case 'p' :
         plbl="Localized Electrons Detector --LED--";
         break;
      case 'P' :
         plbl="Magnitude of Localized Electrons Detector";
         break;
      case 'r' :
         plbl="Region of Slow Electrons --RoSE--";
         break;
      case 's' :
         plbl="Reduced Density Gradient --s--";
         break;
      case 'S':
         plbl="Shannon-Entropy Density";
         break;
      case 'G':
         plbl="Kinetic Energy Density G";
         break;
      case 'K':
         plbl="Kinetic Energy Density K";
         break;
      case 'V':
         plbl="Molecular Electrostatic Potential";
         break;
      default:
         plbl="";
         break;
   }
   return plbl;
}
//**********************************************************************************************
inline string gnuplotFieldTitle(const char p2p)
{
   string plbl;
   switch (p2p) {
      case 'd':
         plbl=string("{/Symbol r}");
         break;
      case 'g':
         plbl=string("|{/Symbol \\321 r}|");
         break;
      case 'l':
         plbl=string("{/Symbol \\321}^2{/Symbol r}");
         break;
      case 'E':
         plbl=string("ELF");
         break;
      case 'p' :
         plbl=string("~{/Bold P}{0.6\\176} [LED]");
         break;
      case 'P' :
         plbl=string("|~{/Bold P}{0.6\\176}| [|LED|]");
         break;
      case 'r' :
         plbl=string("RoSE");
         break;
      case 's' :
         plbl=string("s");
         break;
      case 'S':
         plbl=string("S_{/Symbol r}");
         break;
      case 'L':
         plbl=string("LOL");
         break;
      case 'M':
         plbl=string("|{/Symbol \\321}LOL|");
         break;
      case 'G':
         plbl=string("{/Bold G}");
         break;
      case 'K':
         plbl=string("{/Bold K}");
         break;
      case 'V':
         plbl=string("M.E.P.");
         break;
      default:
         plbl="";
         break;
   }
   return plbl;
}
//**********************************************************************************************
#endif//_FLDTYPESDEF_H_

