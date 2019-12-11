#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
#include "solvmdtools.h"
#include "solfileutils.h"
#include "solscrutils.h"

void VMDTools::writeSimpleVMDScript(string cubename,const char prop2rend,bool addquiet) {
   string vmdname=cubename;
   string tganame=cubename;
   replaceExtensionOfFileName(vmdname,"vmd");
   replaceExtensionOfFileName(tganame,"tga");
   double isovalue=getDefaultIsolvalueForCube(prop2rend);
   ofstream ofil(vmdname.c_str());
   if ( !ofil.good() ) {
      displayErrorMessage(string("Could not open the file '")+vmdname+string("'."));
   }
   ofil << "#VMD file created with DensToolKit." << endl;
   ofil << "#If you want to keep VMD's display on, comment or remove the last line 'quit'." << endl;
   ofil << "display resetview" << endl;
   ofil << "display projection Orthographic" << endl;
   ofil << "display nearclip set 0.000000" << endl;
   ofil << "display depthcue Off" << endl;
   ofil << "color Display Background white" << endl;
   ofil << "axes location Off" << endl;
   ofil << "mol addrep 0" << endl;
   ofil << "mol new {" << cubename
        << "} type {cube} first 0 last -1 step 1 waitfor 1 volsets {0 }" << endl;
   ofil << "mol modstyle 0 0 Licorice 0.300000 12.000000 12.000000" << endl;
   ofil << "mol addrep 0" << endl;
   ofil << "mol modrep 1 0" << endl;
   ofil << "mol modstyle 1 0 Isosurface " << isovalue << " 0 0 0 1 1" << endl;
   ofil << "mol color Name" << endl;
   ofil << "mol material Opaque" << endl;
   ofil << "mol representation Isosurface " << isovalue << " 0 0 0 1 1" << endl;
   ofil << "render TachyonInternal " << tganame << endl;
   ofil << "#" << endl;
   if ( addquiet ) { ofil << "quit" << endl; }
   ofil << "#" << endl;
   ofil.close();
}

