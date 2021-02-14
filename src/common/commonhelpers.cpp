#include <cstdlib>
#include <iostream>
using std::cout;
#include "commonhelpers.h"
#include "atomcolschjmol.h"
#include "atomradiicust.h"

void CommonHelpers::PutNuclei(ofstream &ofil,BondNetWork &bn,int ntbs,\
      const string trnsmat,bool cpkview) {
   int atomn;
   double atrad;
   int nt=ntbs;
   ofil << "//          Nuclei of molecule." << '\n';
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union {" << '\n';
   for (int i=0; i<bn.nNuc; i++) {
      atomn=bn.atNum[i];
      if ( cpkview ) {
         atrad=GetAtomicVDWRadius(atomn)*AUTOMATICSPACEFILLINGRATIO;
      } else {
         atrad=bn.drawAtSize;
      }
      HelpersPOVRay::WriteTransparentSphere(ofil,nt,bn.R[i][0],bn.R[i][1],bn.R[i][2],atrad,
            GetAtomicRColorReal(atomn),GetAtomicGColorReal(atomn),
            GetAtomicBColorReal(atomn),trnsmat);
   }
   ofil << HelpersPOVRay::IndTabsStr(--nt) << "}\n";
}
void CommonHelpers::PutSpecialSpheres(ofstream &ofil,int ntbs,\
      const vector<vector<double> > &sp,const string trnsmat) {
   if ( sp[0].size()<7 ) {
      ScreenUtils::DisplayErrorMessage("CommonHelpers::PutSpecialSpheres needs 7 componets per sp[i] vector!");
      return;
   }
   double atrad;
   int nt=ntbs;
   ofil << "//          Special additional points." << '\n';
   ofil << HelpersPOVRay::IndTabsStr(nt++) << "union {" << '\n';
   for (int i=0; i<sp.size(); i++) {
      atrad=sp[i][3];
      HelpersPOVRay::WriteTransparentSphere(ofil,nt,sp[i][0],sp[i][1],sp[i][2],atrad,
            sp[i][4],sp[i][5],sp[i][6],trnsmat);
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
void CommonHelpers::WriteAngleDeclarations(ofstream &ofil,POVRayConfiguration &pvc) {
   ofil << "#declare GNUPlotAngle1=" << pvc.vecAngView[0]
        << "; // Equivalent to XAngle." << '\n';
   ofil << "#declare GNUPlotAngle2=" << pvc.vecAngView[2]
        << "; // Equivalent to ZAngle." << '\n';
   ofil << "#declare YAngle=" << pvc.vecAngView[1] << ";" << '\n';
}
void CommonHelpers::RenderPovfile(const string &povname,bool verbose) {
   string cmd="dtkpov2png "+povname;
   if ( !verbose ) { cmd+=" 2>/dev/null"; }
   system(cmd.c_str());
   if (verbose) {
      cout << "Rendering done." << '\n';
   }
}

