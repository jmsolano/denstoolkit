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
static void MakeLineGnuplotFile(OptionFlags &opts, string &gnpn,string &outn,char thefield);
static void MakePlaneGnuplotFile(OptionFlags &opts, string &gnpn,string &outn,double dimparam,\
      char thefield);
static void MakeLineDatFile(OptionFlags &opts,string &datnam,GaussWaveFunction &wf,int theaxis,\
      int npts,char thefield);
static void MakePlaneTsvFile(OptionFlags &opts,string &tsvnan,GaussWaveFunction &wf,int theplane,\
      int npts,char thefield);
static void MakeCubeFile(OptionFlags &opts,string &cubnam,GaussWaveFunction &wf,int npts,\
      char thefield,string &strfield);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSPLOT_H_ */

