//
//  custfmtmathfuncts.h
//  
//
//  Created by Juan Manuel Solano on 2013-10-28.
//
//

#ifndef _CUSTFMTMATHFUNCTS_H_
#define _CUSTFMTMATHFUNCTS_H_

#include <iostream>
#include <string>
using std::string;
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"


void writeDatMatAtCrds(string &acfn,bondNetWork &bn);
void writeDatMatCritPtsCrds(string &cpfn,critPtNetWork &cp);
void writeDatMatBondPathCrds(string &bpfn,critPtNetWork &cp);

#endif /* defined(_CUSTFMTMATHFUNCTS_H_) */
