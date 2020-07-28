#ifndef _HELPERSPLOT_H_
#define _HELPERSPLOT_H_
#include "bondnetwork.h"
#include "gausswavefunction.h"
#include "optflags.h"

/* ************************************************************************** */
class HelpersPlot {
/* ************************************************************************** */
public:
/* ************************************************************************** */
static void MakeGnuplotFile(optFlags &opts,string &gnpnam,string &datnam,char p2p,\
      BondNetWork &bn,double lenline,double minval,double maxval,\
      int at1,int at2,double pbcp,double** (&rbgp));
static double EvalFieldProperty(char prop,double (&x)[3],GaussWaveFunction &wf);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSPLOT_H_ */

