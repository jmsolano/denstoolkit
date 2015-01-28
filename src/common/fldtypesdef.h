

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
   REDG,/* Reduced Density Grandient  */\
   ROSE, /* Region of Slow Electrons  */
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
      case 'u' :
         plbl="ScalarCustFld";
         break;
      case 'U' :
         plbl="VectorCustFld";
         break;
      case 'V':
         plbl="MEP";
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
      case 'u':
         plbl="Scalar Custom Field";
         break;
      case 'U' :
         plbl="Vector Custom Field";
         break;
      case 'V':
         plbl="Molecular Electrostatic Potential";
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
      case 'u' :
         plbl=string("S.C.F.");
         break;
      case 'U' :
         plbl=string("V.C.F.");
         break;
      case 'V':
         plbl=string("M.E.P.");
         break;
      default:
         plbl="Unknown";
         break;
   }
   return plbl;
}
//**********************************************************************************************
#endif//_FLDTYPESDEF_H_

