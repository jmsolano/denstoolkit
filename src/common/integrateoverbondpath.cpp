#ifndef _INTEGRATEOVERBONDPATH_CPP_
#define _INTEGRATEOVERBONDPATH_CPP_

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "critptnetwork.h"
#include "integrateoverbondpath.h"
#include "solscrutils.h"

IntegrateOverBondPath::IntegrateOverBondPath(critPtNetWork &ucpn) {
   init();
   if ( !ucpn.iKnowBCPs() ) {
      displayErrorMessage("First seek the critical points!");
      displayWarningMessage("critPtNetWork pointer is set to null!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   if ( !ucpn.iKnowBGPs() ) {
      displayErrorMessage("First compute the bond paths!");
      displayWarningMessage("critPtNetWork pointer is set to null!");
      cout << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   cp=&ucpn;
}
void IntegrateOverBondPath::init(void) {
   cp=NULL;
}

#endif  /* _INTEGRATEOVERBONDPATH_CPP_ */

