#ifndef _HELPERSINTEGRATE_H_
#define _HELPERSINTEGRATE_H_
#include <memory>
using std::shared_ptr;
#include "optflags.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "integrator.h"

/* ************************************************************************** */
class FactoryIntegrator {
/* ************************************************************************** */
public:
   static shared_ptr<Integrator> CreateIntegrator(OptionFlags &options,\
         int argc, char *argv[],GaussWaveFunction &ugwf,\
         BondNetWork &ubnw);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSINTEGRATE_H_ */


