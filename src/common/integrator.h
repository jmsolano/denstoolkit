#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_
#include "gausswavefunction.h"
#include "bondnetwork.h"

/* ************************************************************************** */
class Integrator {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   Integrator();
   virtual void Integrate() = 0;
   virtual double Integral() = 0;
   virtual void DisplayProperties();
   virtual void DisplayResults();
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   double integral;
   GaussWaveFunction *wf;
   BondNetWork *bn;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATOR_H_ */

