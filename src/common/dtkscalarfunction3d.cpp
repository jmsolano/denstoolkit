#include <cstdlib>
#include <iostream>
using std::cout;
#include "dtkscalarfunction3d.h"

bool DTKScalarFunction::prntunknfld1st=true;
DTKScalarFunction::DTKScalarFunction() : Function3D() {
   wf=nullptr;
   sft=ScalarFieldType::NONE;
}
DTKScalarFunction::DTKScalarFunction(GaussWaveFunction &ugwf) : DTKScalarFunction() {
   wf=&ugwf;
   sft=ScalarFieldType::DENS;
}
void DTKScalarFunction::SetScalarFunction(char t) {
   sft=Char2ScalarFieldType(t);
}
