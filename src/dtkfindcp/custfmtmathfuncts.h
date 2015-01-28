

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
