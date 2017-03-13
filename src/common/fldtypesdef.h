/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.2.0
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



#ifndef _FLDTYPESDEF_H_
#define _FLDTYPESDEF_H_
#include <string>
using std::string;
//**********************************************************************************************
enum ScalarFieldType {
   NONE,\
   DENS,/* Electron density (Rho)  */\
   MGRD,/* MagGradRho Density  */\
   LAPD,/* Laplacian Density  */\
   LOLD,/* LOL Density  */\
   ELFD,/* ELF Density  */\
   SENT,/* Shannon Entropy Density  */\
   KEDK,/* Kinetic Energy Density K  */\
   KEDG,/* Kinetic Energy Density G  */\
   MGLD,/* MagGradLOL Density  */\
   GLOL,/* Grad LOL  */\
   MEPD,/* Molecular Electrostatic Potential Density  */\
   LEDV,/* LEDVector  */\
   MLED,/* MagLEDVector  */\
   REDG,/* Reduced Density Gradient  */\
   ROSE,/* Region of Slow Electrons  */\
   VPED,/* Potential Energy Density */\
   NCIS,/* Non Covalent Interactions (NCI) to Reduced Density Gradient */\
   NCIL,/* Non Covalent Interactions (NCI) to Rho */\
   EDFTA, /* DFT Exchange and Corrrelation Energy ($(-3/4)(3\rho/\pi)^{1/3}$) */\
   SCFD, /* Scalar Custom Field Density */\
   VCFD /* Vector Custom Field Density */
};
//**********************************************************************************************
inline char convertScalarFieldType2Char(ScalarFieldType fftt)
{
   char res='d';
   switch ( fftt ) {
      case DENS :
         res='d';
         break;
      case MGRD :
         res='g';
         break;
      case LAPD :
         res='l';
         break;
      case ELFD :
         res='E';
         break;
      case LOLD :
         res='L';
         break;
      case MGLD :
         res='M';
         break;
      case SENT :
         res='S';
         break;
      case KEDK :
         res='K';
         break;
      case KEDG :
         res='G';
         break;
      case GLOL :
         res='M';
         break;
      case MEPD :
         res='V';
         break;
      case LEDV :
         res='p';
         break;
      case MLED :
         res='P';
         break;
      case ROSE :
         res='r';
         break;
      case REDG :
         res='s';
         break;
      case SCFD :
         res='u';
         break;
      case VCFD :
         res='U';
         break;
      case EDFTA :
         res= 'a';
         break;
      case VPED :
         res= 'v';
         break;
      case NCIS :
         res= 'z';
         break;
      case NCIL :
         res= 'Z';
         break;
      case NONE :
      default :
         res='0';
         break;
   }
   return res;
}
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
      case 'a' :
         plbl="ExDFTa";
         break;
      case 'u' :
         plbl="ScalarCustFld";
         break;
      case 'U' :
         plbl="VectorCustFld";
         break;
      case 'V':
         plbl="MEP";
         break;
      case 'v':
         plbl="V.P.E.D.";
         break;
      case 'z':
         plbl="NCIRedDensGrad";
         break;
      case 'Z':
         plbl="NCIRho";
         break;
      default:
         plbl="Unknown";
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
      case 'a' :
         plbl="Exchange And Correlation Energy DFTa";
         break;
      case 'u':
         plbl="Scalar Custom Field";
         break;
      case 'U' :
         plbl="Vector Custom Field";
         break;
      case 'V':
         plbl="Molecular Electrostatic Potential";
         break;
      case 'v':
         plbl="Virial Potential Energy Density";
         break;
      case 'Z':
         plbl="Non Covalent Interactions - Reduced Density Gradient";
         break;
      case 'z':
         plbl="Non Covalent Interactions - Density";
         break;
      default:
         plbl="Unknown Field Type!";
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
      case 'a' :
         plbl="E_{x}(DFT a)";
         break;
      case 'u' :
         plbl=string("S.C.F.");
         break;
      case 'U' :
         plbl=string("V.C.F.");
         break;
      case 'V':
         plbl=string("M.E.P.");
         break;
      case 'v':
         plbl=string("V.P.E.D.");
         break;
      case 'Z':
         plbl=string("NCI -- s");
         break;
      case 'z':
         plbl=string("NCI -- Rho");
         break;
      default:
         plbl="Unknown";
         break;
   }
   return plbl;
}
//**********************************************************************************************
#endif//_FLDTYPESDEF_H_

