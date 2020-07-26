/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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
#ifndef _CUSTFMTMATHFUNCTS_CPP_
#define _CUSTFMTMATHFUNCTS_CPP_
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ofstream;
#include <iomanip>
using std::setprecision;
using std::scientific;
using std::setw;
#include "custfmtmathfuncts.h"

void writeDatMatAtCrds(string &acfn,bondNetWork &bn) {
   ofstream ofil;
   size_t pos;
   string lbl;
   ofil.open(acfn.c_str(),std::ios::out);
   ofil << scientific << setprecision(12);
   for (int i=0; i<bn.nNuc; i++) {
      lbl=bn.atLbl[i];
      pos=lbl.find_last_not_of("0123456789");
      ofil << lbl.substr(0,(++pos));
      ofil << " " << (bn.atNum[i]+1);
      for (int j=0; j<3; j++) {ofil << " " << setw(19) << bn.R[i][j];}
      ofil << endl;
   }
   ofil.close();
   return;
}
void writeDatMatCritPtsCrds(string &cpfn,critPtNetWork &cp) {
   ofstream ofil;
   ofil.open(cpfn.c_str(),std::ios::out);
   ofil << scientific << setprecision(12);
   for (int i=0; i<cp.nACP; i++) {
      ofil << "acp";
      for (int j=0; j<3; j++) {ofil << " " << setw(19) << cp.RACP[i][j];}
      ofil << endl;
   }
   for (int i=0; i<cp.nBCP; i++) {
      ofil << "bcp";
      for (int j=0; j<3; j++) {ofil << " " << setw(19) << cp.RBCP[i][j];}
      ofil << endl;
   }
   for (int i=0; i<cp.nRCP; i++) {
      ofil << "rcp";
      for (int j=0; j<3; j++) {ofil << " " << setw(19) << cp.RRCP[i][j];}
      ofil << endl;
   }
   for (int i=0; i<cp.nCCP; i++) {
      ofil << "ccp";
      for (int j=0; j<3; j++) {ofil << " " << setw(19) << cp.RCCP[i][j];}
      ofil << endl;
   }
   ofil.close();
   return;
}
void writeDatMatBondPathCrds(string &bpfn,critPtNetWork &cp) {
   ofstream ofil;
   ofil.open(bpfn.c_str(),std::ios::out);
   ofil << scientific << setprecision(12);
   int npts,nbps;
   nbps=cp.nBCP;
   for (int i=0; i<nbps; i++) {
      npts=cp.conBCP[i][2];
      for (int j=0; j<npts; j++) {
         ofil << cp.RBGP[i][j][0] << " " << cp.RBGP[i][j][1] << " " << cp.RBGP[i][j][2] << endl;
      }
   }
   ofil.close();
   return;
}

#endif /* defined(_CUSTFMTMATHFUNCTS_CPP_) */

