#include <cstdlib>
#include <iostream>
using std::cout;
#include "commonhelpers.h"
#include "atomcolschjmol.h"
#include "atomradiicust.h"
#include "povraytools.h"

void CommonHelpers::PutNuclei(ofstream &ofil,BondNetWork &bn,int ntbs,\
      const string trnsmat) {
   int atomn;
   double atrad;
   int nt=ntbs;
   ofil << "//          Nuclei of molecule." << '\n';
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union {" << '\n';
   for (int i=0; i<bn.nNuc; i++) {
      atomn=bn.atNum[i];
      atrad=bn.drawAtSize;
      HelpersPOVRay::WriteTransparentSphere(ofil,nt,bn.R[i][0],bn.R[i][1],bn.R[i][2],atrad,
            GetAtomicRColorReal(atomn),GetAtomicGColorReal(atomn),
            GetAtomicBColorReal(atomn),trnsmat);
   }
   ofil << HelpersPOVRay::IndTabsStr(--nt) << "}\n";
}
void CommonHelpers::PutBonds(ofstream &ofil,BondNetWork &bn,const int ntbs,\
         const string trnsmbnd) {
   string pigmstr="transmit ";
   pigmstr+=trnsmbnd;
   int nt=ntbs;
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union{" << '\n';
   int k=0,atni,atnk;
   double startpt[3],frak1;
   for (int i=0; i<bn.nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bn.bNet[i][j];
         atni=bn.atNum[i];
         atnk=bn.atNum[k];
         frak1=GetAtomicVDWRadius(atni)/(GetAtomicVDWRadius(atni)+GetAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=bn.R[i][l]*(1.0e0-frak1)+bn.R[k][l]*frak1;
         }
         if (k>0) {
            HelpersPOVRay::WriteCylinder(ofil,nt,
                  bn.R[i][0],bn.R[i][1],bn.R[i][2],
                  startpt[0],startpt[1],startpt[2],
                  bn.drawStickSize,
                  GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
                  GetAtomicBColorReal(atni),pigmstr);
            HelpersPOVRay::WriteCylinder(ofil,nt,
                  startpt[0],startpt[1],startpt[2],
                  bn.R[k][0],bn.R[k][1],bn.R[k][2],
                  bn.drawStickSize,
                  GetAtomicRColorReal(atnk),GetAtomicGColorReal(atnk),
                  GetAtomicBColorReal(atnk),pigmstr);
         }
      }
      HelpersPOVRay::WriteSphere(ofil,nt,bn.R[i][0],bn.R[i][1],bn.R[i][2],
            bn.drawStickSize,
            GetAtomicRColorReal(atni),GetAtomicGColorReal(atni),
            GetAtomicBColorReal(atni));
   }
   ofil << HelpersPOVRay::IndTabsStr(--nt) << "}" << '\n';
}

