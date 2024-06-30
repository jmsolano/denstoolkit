/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
#ifndef _FLDTYPESDEF_H_
#define _FLDTYPESDEF_H_
#include <string>
using std::string;
/* 3D Fields  */
enum ScalarFieldType {
   NONE,\
   DENS,/* Electron density (Rho)  */\
   DENSM,/* Electron density (Rho) in Momentum Space  */\
   MGRD,/* MagGradRho Density  */\
   LAPD,/* Laplacian Density  */\
   LOLD,/* LOL Density  */\
   ELFD,/* ELF Density  */\
   SENT,/* Shannon Entropy Density  */\
   SENTM,/* Shannon Entropy Density in Momentum Space  */\
   KEDK,/* Kinetic Energy Density K  */\
   KEDKM,/* Kinetic Energy Density K in Momentum Space  */\
   KEDG,/* Kinetic Energy Density G  */\
   MGLD,/* MagGradLOL Density  */\
   GLOL,/* Grad LOL  */\
   MEPD,/* Molecular Electrostatic Potential Density  */\
   LEDV,/* LEDVector  */\
   MLED,/* MagLEDVector  */\
   REDG,/* Reduced Density Gradient  */\
   ROSE,/* Region of Slow Electrons  */\
   VPED,/* Virial Potential Energy Density */\
   NCIS,/* Non Covalent Interactions (NCI) -- Reduced Density Gradient */\
   NCIL,/* Non Covalent Interactions (NCI) -- Lambda2 */\
   EDFTA, /* DFT Exchange and Corrrelation Energy ($(-3/4)(3\rho/\pi)^{1/3}$) */\
   ELLPY, /*!< Ellipticity  */\
   DORI, /*!< Density Overlap Regions Indicator  */\
   SPND, /*!< [SP]i[N] [D]ensity  */\
   RHO2, /*!< \f$\rho^2\f$ (one eletron density disequilibrium)  */\
   SCFD, /* Scalar Custom Field Density */\
   VCFD /* Vector Custom Field Density */
};
inline char ConvertScalarFieldType2Char(ScalarFieldType fftt) {
   char res='d';
   switch ( fftt ) {
      case DENS :
         res='d';
         break;
      case DENSM :
         res='m';
         break;
      case MGRD :
         res='g';
         break;
      case LAPD :
         res='l';
         break;
      case LOLD :
         res='L';
         break;
      case ELFD :
         res='E';
         break;
      case SENT :
         res='S';
         break;
      case SENTM :
         res='T';
         break;
      case KEDK :
         res='K';
         break;
      case KEDKM :
         res='k';
         break;
      case KEDG :
         res='G';
         break;
      case MGLD :
         res='M';
         break;
      case GLOL :
         res='N';
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
      case REDG :
         res='s';
         break;
      case ROSE :
         res='r';
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
      case EDFTA :
         res= 'a';
         break;
      case ELLPY :
         res= 'e';
         break;
      case DORI :
         res= 'D';
         break;
      case SPND :
         res= 'b';
         break;
      case RHO2 :
         res= 'q';
         break;
      case SCFD :
         res='u';
         break;
      case VCFD :
         res='U';
         break;
      case NONE :
      default :
         res='0';
         break;
   }
   return res;
}
inline ScalarFieldType Char2ScalarFieldType(const char prop) {
   ScalarFieldType res=ScalarFieldType::NONE;
   switch (prop) {
      case 'd':
         res=DENS;
         break;
      case 'm':
         res=DENSM;
         break;
      case 'g' :
         res=MGRD;
         break;
      case 'l' :
         res=LAPD;
         break;
      case 'L' :
         res=LOLD;
         break;
      case 'E' :
         res=ELFD;
         break;
      case 'S' :
         res=SENT;
         break;
      case 'T' :
         res=SENTM;
         break;
      case 'K' :
         res=KEDK;
         break;
      case 'k' :
         res=KEDKM;
         break;
      case 'G' :
         res=KEDG;
         break;
      case 'M' :
         res=MGLD;
         break;
      case 'N' :
         res=GLOL;
         break;
      case 'V' :
         res=MEPD;
         break;
      case 'p' :
         res=LEDV;
         break;
      case 'P' :
         res=MLED;
         break;
      case 's' :
         res=REDG;
         break;
      case 'r' :
         res=ROSE;
         break;
      case 'v' :
         res=VPED;
         break;
      case 'z' :
         res=NCIS;
         break;
      case 'Z' :
         res=NCIL;
         break;
      case 'a' :
         res=EDFTA;
         break;
      case 'e' :
         res=ELLPY;
         break;
      case 'D' :
         res=DORI;
         break;
      case 'b' :
         res=SPND;
         break;
      case 'q' :
         res=RHO2;
         break;
      case 'u' :
         res=SCFD;
         break;
      case 'U' :
         res=VCFD;
         break;
      default :
         res=NONE;
         break;
   }
   return res;
}
inline string GetFieldTypeKeyShort(const char prop) {
   string plbl="";
   switch (prop) {
      case 'd':
         plbl="Rho";
         break;
      case 'm':
         plbl="RhoMom";
         break;
      case 'g':
         plbl="MagGradRho";
         break;
      case 'l':
         plbl="LapRho";
         break;
      case 'L':
         plbl="LOL";
         break;
      case 'E':
         plbl="ELF";
         break;
      case 'S':
         plbl="ShannEnt";
         break;
      case 'T':
         plbl="ShannEntMom";
         break;
      case 'K':
         plbl="KinEnerDensK";
         break;
      case 'k':
         plbl="KinEnerDensKMom";
         break;
      case 'G':
         plbl="KinEnerDensG";
         break;
      case 'M':
         plbl="MagGradLOL";
         break;
      case 'N':
         plbl="GradLOL";
         break;
      case 'V':
         plbl="MEP";
         break;
      case 'p' :
         plbl="LED";
         break;
      case 'P' :
         plbl="MagLED";
         break;
      case 's' :
         plbl="RedDensGrad";
         break;
      case 'r' :
         plbl="RoSE";
         break;
      case 'v':
         plbl="VPED";
         break;
      case 'z':
         plbl="NCIRedDensGrad";
         break;
      case 'Z':
         plbl="NCILambda2";
         break;
      case 'a' :
         plbl="ExDFTa";
         break;
      case 'e':
         plbl="Ellipticity";
         break;
      case 'D' :
         plbl="DORI";
         break;
      case 'b' :
         plbl="SpinDensity";
         break;
      case 'q' :
         plbl="OneElecDiseq";
         break;
      case 'u' :
         plbl="ScalarCustFld";
         break;
      case 'U' :
         plbl="VectorCustFld";
         break;
      default:
         plbl="Unknown";
         break;
   }
   return plbl;
}
inline string GetFieldTypeKeyLong(const char prop) {
   string plbl="";
   switch (prop) {
      case 'd':
         plbl="Electron Density --Rho--";
         break;
      case 'm':
         plbl="Electron Density --Rho-- in Momentum Space";
         break;
      case 'g':
         plbl="Magnitude of the Gradient of the Electron Density";
         break;
      case 'l':
         plbl="Laplacian of Electron Density";
         break;
      case 'L':
         plbl="Localized Orbital Locator --LOL--";
         break;
      case 'E':
         plbl="Electron Localization Function --ELF--";
         break;
      case 'S':
         plbl="Shannon-Entropy Density";
         break;
      case 'T':
         plbl="Shannon-Entropy Density in Momentum Space";
         break;
      case 'K':
         plbl="Kinetic Energy Density K";
         break;
      case 'k':
         plbl="Kinetic Energy Density K in Momentum Space";
         break;
      case 'G':
         plbl="Kinetic Energy Density G";
         break;
      case 'M':
         plbl="Magnitude of the Gradrient of LOL";
         break;
      case 'N':
         plbl="Gradient of LOL";
         break;
      case 'V':
         plbl="Molecular Electrostatic Potential";
         break;
      case 'p' :
         plbl="Localized Electrons Detector --LED--";
         break;
      case 'P' :
         plbl="Magnitude of Localized Electrons Detector";
         break;
      case 's' :
         plbl="Reduced Density Gradient --s--";
         break;
      case 'r' :
         plbl="Region of Slow Electrons --RoSE--";
         break;
      case 'v':
         plbl="Virial Potential Energy Density";
         break;
      case 'z':
         plbl="Non Covalent Interactions - Reduced Density Gradient";
         break;
      case 'Z':
         plbl="Non Covalent Interactions - Density";
         break;
      case 'a' :
         plbl="Exchange And Correlation Energy DFTa";
         break;
      case 'e':
         plbl="Ellipticity";
         break;
      case 'D' :
         plbl="Density Overlap Regions Indicator --DORI--";
         break;
      case 'b' :
         plbl="Spin density";
         break;
      case 'q' :
         plbl="One electron disequilibrium --Rho squared--";
         break;
      case 'u':
         plbl="Scalar Custom Field";
         break;
      case 'U' :
         plbl="Vector Custom Field";
         break;
      default:
         plbl="Unknown Field Type!";
         break;
   }
   return plbl;
}
inline string GnuplotFieldTitle(const char p2p) {
   string plbl;
   switch (p2p) {
      case 'd':
         plbl=string("{/Symbol r}");
         break;
      case 'm':
         plbl=string("{/Symbol p}");
         break;
      case 'g':
         plbl=string("|{/Symbol \\321 r}|");
         break;
      case 'l':
         plbl=string("{/Symbol \\321}^2{/Symbol r}");
         break;
      case 'e':
         plbl=string("{/Symbol e}");
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
      case 'T':
         plbl=string("{/Symbol P}_{/Symbol r}");
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
      case 'k':
         plbl=string("{/Bold K}_P");
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
      case 'z':
         plbl=string("NCI -- s");
         break;
      case 'Z':
         plbl=string("NCI -- {/Symbol L}_2");
         break;
      case 'D' :
         plbl=string("DORI");
         break;
      case 'b' :
         plbl=string("{/Symbol r}_{/Symbol a}-{/Symbol r}_{/Symbol b}");
         break;
      case 'q' :
         plbl=string("{/Symbol r}^2");
         break;
      default:
         plbl="Unknown";
         break;
   }
   return plbl;
}
inline double GetDefaultIsolvalueForCube(const char p2p) {
   double isoval=0.01e0;
   switch ( p2p ) {
      case 'd':
         isoval=0.01e0;
         break;
      case 'D' :
         isoval=0.90e0;
         break;
      case 's':
         isoval=0.2e0;
         break;
      case 'z' :
         isoval=0.055e0;
         break;
      case 'Z' :
         isoval=0.55e0;
         break;
      default :
         isoval=0.01e0;
         break;
   }
   return isoval;
}
/* 6D fields  */
enum ScalarField6DType {
   NONE6D,\
   DM1P, /*!< [D]ensity [M]atrix of order [1]; 6D field  */\
   GDM1P, /*!< [G]radient of [D]ensity [M]atrix of order [1]; 6D vector field  */\
   LDM1P, /*!< [L]aplacian of [D]ensity [M]atrix of order [1]; 6D field  */\
   ADM1P, /*!< [A]lpha [D]ensity [M]atrix of order [1]; 6D field  */\
   BDM1P, /*!< [B]eta [D]ensity [M]atrix of order [1]; 6D field  */\
   DM1M, /*!< [D]ensity [M]atrix of order [1] in [M]omentum spac3; 6D field  */\
   ADM1M, /*!< [A]lpha [D]ensity [M]atrix of order [1] in [M]omentum space; 6D field  */\
   BDM1M, /*!< [A]lpha [D]ensity [M]atrix of order [1] in [M]omentum space; 6D field  */\
   CPDFP, /*!< [C]losed-shell[P]air [D]ensity [F]unction; 6D field  */\
   CPDFM, /*!< [C]losed-shell [P]air [D]ensity [F]unction in [M]omentum space; 6D field  */\
   OPDFP, /*!< [O]pen-shell [P]air [D]ensity [F]unction; 6D field  */\
   OPDFM, /*!< [O]pen-shell [P]air [D]ensity [F]unction in [M]omentum space; 6D field  */\
   SCF6, /* Scalar Custom Field 6D */\
   VCF6 /* Vector Custom Field 6D */
};
inline char ScalarFieldType6D2Char(ScalarField6DType fftt) {
   char res='0';
   switch ( fftt ) {
      case DM1P  :
         res='g';
         break;
      case GDM1P :
         res='n';
         break;
      case LDM1P :
         res='l';
         break;
      case DM1M :
         res='G';
         break;
      case ADM1P :
         res='a';
         break;
      case ADM1M:
         res='A';
         break;
      case BDM1P:
         res='b';
         break;
      case BDM1M:
         res='B';
         break;
      case CPDFP:
         res='c';
         break;
      case CPDFM:
         res='C';
         break;
      case OPDFP :
         res='o';
         break;
      case OPDFM:
         res='O';
         break;
      case SCF6 :
         res='u';
         break;
      case VCF6 :
         res='U';
         break;
      case NONE6D :
      default   :
         break;
   }
   return res;
}
inline ScalarField6DType Char2ScalarField6DType(const char prop) {
   ScalarField6DType res=ScalarField6DType::NONE6D;
   switch (prop) {
      case 'a' :
         res=ADM1P;
         break;
      case 'A' :
         res=ADM1M;
         break;
      case 'b' :
         res=BDM1P;
         break;
      case 'B' :
         res=BDM1M;
         break;
      case 'c' :
         res=CPDFP;
         break;
      case 'C' :
         res=CPDFM;
         break;
      case 'g' :
         res=DM1P;
         break;
      case 'G' :
         res=DM1M ;
         break;
      case 'l' :
         res=LDM1P;
         break;
      case 'n' :
         res=GDM1P;
         break;
      case 'o' :
         res=OPDFP;
         break;
      case 'O' :
         res=OPDFM;
         break;
      case 'u' :
         res=SCF6;
         break;
      case 'U' :
         res=VCF6;
         break;
      case '0' :
      default  :
         break;
   }
   return res;
}
inline string GetField6DTypeKeyShort(const char prop) {
   string plbl="";
   switch (prop) {
      case 'a' :
         plbl="AlphaDM1";
         break;
      case 'A' :
         plbl="AlphaDM1MomSp";
         break;
      case 'b' :
         plbl="BetaDM1";
         break;
      case 'B' :
         plbl="BetaDM1MomSp";
         break;
      case 'c' :
         plbl="CSPairDens";
         break;
      case 'C' :
         plbl="CSPairDensMomSp";
         break;
      case 'g':
         plbl="DM1";
         break;
      case 'G':
         plbl="DM1MomSp";
         break;
      case 'l' :
         plbl="LapDM1";
         break;
      case 'n' :
         plbl="GradDM1";
         break;
      case 'o' :
         plbl="OSPairDens";
         break;
      case 'O' :
         plbl="OSPairDensMomSp";
         break;
      case 'u' :
         plbl="ScalarCustFld6D";
         break;
      case 'U' :
         plbl="VectorCustFld6D";
         break;
      default:
         plbl="Unknown6DFld";
         break;
   }
   return plbl;
}
inline string GetField6DTypeKeyLong(const char prop) {
   string plbl="";
   switch (prop) {
      case 'a' :
         plbl="Alpha density matrix of order 1";
         break;
      case 'A' :
         plbl="Alpha density matrix of order 1 in momentum space";
         break;
      case 'b' :
         plbl="Beta density matrix of order 1";
         break;
      case 'B' :
         plbl="Beta density matrix of order 1 in momentum space";
         break;
      case 'c' :
         plbl="Closed-shell pair density function";
         break;
      case 'C' :
         plbl="Closed-shell pair density function in momentum space";
         break;
      case 'g':
         plbl="Density Matrix of order 1";
         break;
      case 'G':
         plbl="Density matrix of order 1 in momentum space";
         break;
      case 'n' :
         plbl="Gradient of Density Matrix of order 1";
         break;
      case 'l' :
         plbl="Laplacian of Density Matrix of order 1";
         break;
      case 'o' :
         plbl="Open-shell pair density function";
         break;
      case 'O' :
         plbl="Open-shell pair density function in momentum space";
         break;
      case 'u' :
         plbl="Scalar custom field 6D";
         break;
      case 'U' :
         plbl="Vector custom field 6D";
         break;
      default:
         plbl="Unknown 6D field";
         break;
   }
   return plbl;
}
inline string GnuplotField6DTitle(const char prop) {
   string plbl="";
   switch (prop) {
      case 'a' :
         plbl="{/Symbol G}_1^{/Symbol a}";
         break;
      case 'A' :
         plbl="~{/Symbol G}{0.3\\~}_1^{/Symbol a}";
         break;
      case 'b' :
         plbl="{/Symbol G}_1^{/Symbol b}";
         break;
      case 'B' :
         plbl="~{/Symbol G}{0.3\\~}_1^{/Symbol b}";
         break;
      case 'c' :
         plbl="{/Symbol r}_{2,cs}";
         break;
      case 'C' :
         plbl="~{/Symbol r}{0.1\\~}_{2,cs}";
         break;
      case 'g':
         plbl="{/Symbol G}_1";
         break;
      case 'G':
         plbl="~{/Symbol G}{0.3\\~}_1";
         break;
      case 'n' :
         plbl="{/Symbol \\321 G}_1";
         break;
      case 'l' :
         plbl="{/Symbol \\321}^2{/Symbol G}_1";
         break;
      case 'o' :
         plbl="{/Symbol r}_{2,os}";
         break;
      case 'O' :
         plbl="~{/Symbol r}{0.1\\~}_{2,os}";
         break;
      case 'u' :
         plbl=string("S.C.F. 6D");
         break;
      case 'U' :
         plbl=string("V.C.F. 6D");
         break;
      default:
         plbl="Unknown 6D Field";
         break;
   }
   return plbl;
}
inline bool Is6DMomSpaceField(const char prop) {
   bool res=false;
   switch ( prop ) {
      case 'G' :
      case 'H' :
      case 'P' :
      case 'Q' :
      case 'U' :
         res=true;
         break;
      default :
         res=false;
         break;
   }
   return res;
}
inline bool Is6DPosSpaceField(const char prop) { return !Is6DMomSpaceField(prop); }
#endif//_FLDTYPESDEF_H_

