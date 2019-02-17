#ifndef _HELPERPLOTS_H_
#define _HELPERPLOTS_H_
#include "optflags.h"
#include <fstream>
using std::ofstream;
#include "../common/bondnetwork.h"
#include "../common/demat1critptnetworksl.h"


/* ************************************************************************** */
class HelperPlot {
/* ************************************************************************** */
public:
   static void generateMainDiagPlot(optFlags &options,const string &datname,\
         bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
         solreal linelength,solreal frange);
   static void generateSecDiagPlot(optFlags &options,const string &datname,\
         bondNetWork &bn,int idx1,int idx2,solreal minval2plot,solreal maxval2plot,\
         solreal linelength,solreal frange);
   static void generate3DPlot(optFlags &options,const string &tsvname,solreal minval2plot,\
         solreal maxval2plot,solreal linelength,int nptsinline);
   static void generateHeatMap(optFlags &options,char *argv[],const string &tsvname,\
         bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),int nptsinline,solreal minval2plot,\
         solreal maxval2plot,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
   static void generateVectorField(optFlags &options,char *argv[],const string &tsvname,\
         bondNetWork &bn,DeMat1CriticalPointNetworkSL &cp,solreal **(xx),int nptsinline,solreal minval2plot,\
         solreal maxval2plot,solreal maggradmin,solreal maggradmax,solreal linelength,solreal md1lmin,solreal md1dmax,int idx1,int idx2);
   static void addHeaderInfo2GNP(ofstream &ofil,solreal minval,solreal maxval,\
         solreal dimpar,string datortsv,string axis="y");
   static void generateGNPEPSAndPDFNamesFromDATORTSV(const string &dnam,\
         string &gnam,string &enam,string &pnam);
   static void generateEPSAndPDFNamesFromGNP(const string &gnam,\
         string &enam,string &pnam);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _HELPERPLOTS_H_ */

