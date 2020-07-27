#ifndef _HELPERSPLOT_H_
#define _HELPERSPLOT_H_
#include "screenutils.h"
#include "fileutils.h"
#include "optflags.h"
#include "gausswavefunction.h"

/* ************************************************************************** */
class HelpersPlot {
/* ************************************************************************** */
public:
/* ************************************************************************** */
static void MakeLineGnuplotFile(optFlags &opts, string &gnpn,string &outn,char thefield);
static void MakePlaneGnuplotFile(optFlags &opts, string &gnpn,string &outn,double dimparam,\
      char thefield);
static void MakeLineDatFile(optFlags &opts,string &datnam,GaussWaveFunction &wf,int theaxis,\
      int npts,char thefield);
static void MakePlaneTsvFile(optFlags &opts,string &tsvnan,GaussWaveFunction &wf,int theplane,\
      int npts,char thefield);
static void MakeCubeFile(optFlags &opts,string &cubnam,GaussWaveFunction &wf,int npts,\
      char thefield,string &strfield);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSPLOT_H_ */

