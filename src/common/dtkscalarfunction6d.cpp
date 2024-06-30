#include <cstdlib>
#include <iostream>
using std::cout;
#include "dtkscalarfunction6d.h"
#include "screenutils.h"
#include "fldtypesdef.h"

bool DTKScalarFunction6D::prntunknfld1st=true;
DTKScalarFunction6D::DTKScalarFunction6D() {
   wf=nullptr;
   sft=ScalarField6DType::NONE6D;
}
DTKScalarFunction6D::DTKScalarFunction6D(GaussWaveFunction &ugwf) : DTKScalarFunction6D() {
   wf=&ugwf;
   sft=ScalarField6DType::DM1P;
}
void DTKScalarFunction6D::SetScalarFunction(char t) {
   sft=Char2ScalarField6DType(t);
}
