#include <cstdlib>
#include <iostream>
using std::cout;
#include "integrator.h"

Integrator::Integrator() {
   wf=nullptr;
   bn=nullptr;
   integral=0.0e0;
}
void Integrator::DisplayProperties() {
   cout << "Base integrator, nothing to display..." << '\n';
}
void Integrator::DisplayResults() {
   DisplayProperties();
}
void Integrator::WriteResults(ofstream &ofil) {
   ofil << "#Base integrator, nothing to write!..." << '\n';
}


