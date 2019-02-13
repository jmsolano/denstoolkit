#ifndef _HELPERPLOTS_H_
#define _HELPERPLOTS_H_
#include "optflags.h"
#include "../common/bondnetwork.h"
#include "../common/demat1critptnetworksl.h"

/* ************************************************************************** */
//class HelperPlot {
/* ************************************************************************** */
//public:
   void generateMainDiagPlot(optFlags &options,const string &datname,\
         bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
         solreal linelength,solreal frange);
   void generateSecDiagPlot(optFlags &options,const string &datname,\
         bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
         solreal linelength,solreal frange);
   void generate3DPlot(optFlags &options,const string &tsvname,solreal minval2plot,\
         solreal maxval2plot,solreal linelength,int nptsinline);
   void generateHeatMap(optFlags &options,char *argv[],const string &tsvname,\
         bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),solreal minval2plot,\
         solreal maxval2plot,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
   void generateVectorField(optFlags &options,char *argv[],const string &tsvname,\
         bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),solreal minval2plot,\
         solreal maxval2plot,solreal maggradmin,solreal maggradmax,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
/* ************************************************************************** */
//protected:
/* ************************************************************************** */
//};
/* ************************************************************************** */


#endif  /* _HELPERPLOTS_H_ */

