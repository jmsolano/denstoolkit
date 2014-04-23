//
//  custfmtmathfuncts.cpp
//  
//
//  Created by Juan Manuel Solano on 2013-10-28.
//
//

#ifndef _CUSTFMTMATHFUNCTS_CPP_
#define _CUSTFMTMATHFUNCTS_CPP_

#include "custfmtmathfuncts.h"
#include <fstream>
using std::ofstream;
#include <iomanip>
using std::setprecision;
using std::scientific;
using std::setw;

void writeDatMatAtCrds(string &acfn,bondNetWork &bn)
{
   ofstream ofil;
   size_t pos;
   string lbl;
   ofil.open(acfn.c_str(),ios::out);
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

void writeDatMatCritPtsCrds(string &cpfn,critPtNetWork &cp)
{
   ofstream ofil;
   ofil.open(cpfn.c_str(),ios::out);
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

void writeDatMatBondPathCrds(string &bpfn,critPtNetWork &cp)
{
   ofstream ofil;
   ofil.open(bpfn.c_str(),ios::out);
   ofil << scientific << setprecision(12);
   int npts,nbps;
   nbps=cp.nBCP;
   for (int i=0; i<nbps; i++) {
      npts=cp.atBCP[i][2];
      for (int j=0; j<npts; j++) {
         ofil << cp.RBGP[i][j][0] << " " << cp.RBGP[i][j][1] << " " << cp.RBGP[i][j][2] << endl;
      }
   }
   ofil.close();
   return;
}


#endif /* defined(_CUSTFMTMATHFUNCTS_CPP_) */

