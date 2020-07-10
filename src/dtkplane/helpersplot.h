#ifndef _HELPERSPLOT_H_
#define _HELPERSPLOT_H_
#include <fstream>
using std::ofstream;
#include "optflags.h"
#include "bondnetwork.h"
#include "wfgrid2d.h"

/* ************************************************************************** */
class HelpersPlot {
/* ************************************************************************** */
public:
/* ************************************************************************** */
static void makeGnuplotFile(optFlags &opts, string &gnpnam,string &tsvnam,char p2p,solreal dimparam,
                     bondNetWork &bn,int a1,int a2,int a3,waveFunctionGrid2D &grd);
static void addGnuplotRenderingHelpEPS2PDF(ofstream &ofil,const string &gpnnam,const string &epsnam);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERSPLOT_H_ */

