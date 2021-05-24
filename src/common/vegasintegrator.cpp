#include <cmath>
#include <vector>
using std::vector;
#include <algorithm>
using std::numeric_limits;
#include <iostream>
using std::cout;

#include "vegasintegrator.h"
#include "gausswavefunction.h"


VegasIntegrator::VegasIntegrator() {
   wf=nullptr;
   integral=0.0e0;
   countIter=countEval=0;
   stopIterating=false;
}
VegasIntegrator::VegasIntegrator(GaussWaveFunction &uwf) : VegasIntegrator() {
   wf=&uwf;
}
double VegasIntegrator::Function(double x,double y,double z){
   if (param.function == "rho") return wf->EvalDensity(x,y,z);
   else if (param.function == "lap_rho") return wf->EvalLapRho(x,y,z);
   else if (param.function == "ELF") return wf->EvalELF(x,y,z);
   else if (param.function == "shannon_rho") return wf->EvalShannonEntropy(x,y,z);
   else if (param.function == "shannon_ftrho") return wf->EvalMomentumShannonEntropy(x,y,z);
   else if (param.function == "maggrad_rho") return wf->EvalMagGradRho(x,y,z);
   else if (param.function == "LOL") return wf->EvalLOL(x,y,z);
   else if (param.function == "kineticdens_g") return wf->EvalKineticEnergyG(x,y,z);
   else if (param.function == "kineticdens_k") return wf->EvalKineticEnergyK(x,y,z);
   else if (param.function == "ellipticity") return wf->EvalEllipticity(x,y,z);
   else if (param.function == "ftrho") return wf->EvalFTDensity(x,y,z);
   else if (param.function == "ftkinetic_k") return wf->EvalFTKineticEnergy(x,y,z);
   else if (param.function == "maggrad_LOL") return wf->EvalMagGradLOL(x,y,z);
   else if (param.function == "MEP") return wf->EvalMolElecPot(x,y,z);
   else if (param.function == "mag_LED") return wf->EvalMagLED(x,y,z);
   else if (param.function == "reduced_densgrad") return wf->EvalReducedDensityGradient(x,y,z);
   else if (param.function == "RoSE") return wf->EvalRoSE(x,y,z);
   else if (param.function == "customfield") return wf->EvalCustomScalarField(x,y,z);
   else if (param.function == "potentialdens") return wf->EvalVirialPotentialEnergyDensity(x,y,z);
   else if (param.function == "reduced_densgrad_NCI") return wf->EvalNCIs(x,y,z);
   else return wf->EvalDensity(x,y,z);
}
void VegasIntegrator::Integrate(void) {
   vector<vector<double> > interval(3,vector<double>(param.numOfIntervals+1));
   vector<vector<double> > meanIntegral(3,vector<double>(param.numOfIntervals));

   while (repeatIntegral == true || (stopIterating == false && countIter < param.iterations)) {
      weightedAverage=inverseVariance=chiSquare=0;
      countIter=0;
      repeatIntegral = stopIterating = false;
      for (int j=0; j<3; j++) {
	 for (int i=0; i<param.numOfIntervals+1; i++) interval[j][i] = xMin[j]+i*1./param.numOfIntervals*width[j];
      }

      while (repeatIntegral == false && stopIterating == false && countIter < param.iterations) {
	 countIter++;
	 // cout << "Iteration: " << countIter << endl;

	 MonteCarloIntegration(interval,meanIntegral);
	 // cout << sampling.simple << " " << sampling.square << endl;

	 varPerIt = VariancePerIteration(sampling.simple,sampling.square);

	 if (countIter > param.termalization) {
	    weightedAverage += sampling.simple/varPerIt;
	    inverseVariance += 1./varPerIt;

	    integral = weightedAverage/inverseVariance;
	    variance = 1./inverseVariance;
	 }

	 // cout << "sampling= " << sampling.simple << endl;
	 // cout << "integral= " << integral << endl;

	 stopIterating = AlteratesAverageIntegral(meanIntegral);
	 AlteratesIncrements(interval,meanIntegral);

	 if (countIter > param.termalization) {
	    chiSquare += ChiSquare(sampling.simple,varPerIt,integral);
	    if (chiSquare > countIter) {
	       // cout << "Chi Square" << endl;
	       // getchar();
	       repeatIntegral = true;
	    }
	 }
      }
   }

   countEval = param.numOfPoints*countIter;
}
void VegasIntegrator::MonteCarloIntegration(vector<vector<double> > interval,vector<vector<double> > &meanIntegral) {
   int floorYNg[3];
   double decimalYNg,yNg,deltaXi,jacobian;
   double xi[3];
   vector<vector<int> > countPointsSampled(3,vector<int>(param.numOfIntervals,0));

   sampling.simple = sampling.square = 0;
   for (int j=0; j<3; j++) {
      for (int i=0; i<param.numOfIntervals; i++) meanIntegral[j][i] = 0;
   }

   for (long int k=0; k<param.numOfPoints; k++) {
      jacobian = 1;
      for (int j=0; j<3; j++) {
	 yNg = dis(engine)*param.numOfIntervals;
	 floorYNg[j] = int(yNg);
	 decimalYNg = yNg-floorYNg[j];

	 deltaXi = interval[j][floorYNg[j]+1]-interval[j][floorYNg[j]];
	 xi[j] = interval[j][floorYNg[j]]+deltaXi*decimalYNg;

	 jacobian *= deltaXi*param.numOfIntervals;

	 countPointsSampled[j][floorYNg[j]]++;
      }

      sampling.simple += Function(xi[0],xi[1],xi[2])*jacobian;
      sampling.square += Power(Function(xi[0],xi[1],xi[2])*jacobian,2);
      for (int j=0; j<3; j++) meanIntegral[j][floorYNg[j]] += Power(Function(xi[0],xi[1],xi[2])*jacobian,2);
      // sampling.simple += wf->EvalDensity(xi[0],xi[1],xi[2])*jacobian;
      // sampling.square += Power(wf->EvalDensity(xi[0],xi[1],xi[2])*jacobian,2);
      // for (int j=0; j<3; j++) meanIntegral[j][floorYNg[j]] += Power(wf->EvalDensity(xi[0],xi[1],xi[2])*jacobian,2);
   }

   sampling.simple /= param.numOfPoints;
   sampling.square /= param.numOfPoints;

   for (int j=0; j<3; j++) {
      for (int i=0; i<param.numOfIntervals; i++) {
	 if (countPointsSampled[j][i] > 0) meanIntegral[j][i] /= countPointsSampled[j][i];
	 // if (countPointsSampled[j][i] > 0) meanIntegral[j][i] /= param.numOfPoints*1./param.numOfIntervals;
	 // if (meanIntegral[j][i] == 0) meanIntegral[j][i] = numeric_limits<double>::min();
      }
   }

   // cout << "points sampled per interval: ";
   // for (int j=0; j<3; j++) {
      // for (int i=0; i<param.numOfIntervals; i++) cout << countPointsSampled[j][i] << " ";
      // cout << endl;
   // }
}
bool VegasIntegrator::AlteratesAverageIntegral(vector<vector<double> > &meanIntegral) {
   double sumInferior=0,sumSuperior=0,convRate=param.convergenceRate;
   double sumMeanIntegral=0;
   vector<double> meanIntegralCopy(param.numOfIntervals);

   for (int j=0; j<3; j++){
      sumInferior += *min_element(meanIntegral[j].begin(),meanIntegral[j].end());
      sumSuperior += *max_element(meanIntegral[j].begin(),meanIntegral[j].end());

      // cout << *min_element(meanIntegral[j].begin(),meanIntegral[j].end()) << " " << *max_element(meanIntegral[j].begin(),meanIntegral[j].end()) << endl;
   }
   sumInferior /= 3;
   sumSuperior /= 3;

   if (fabs(sumSuperior-sumInferior) <= param.tolerance && countIter > param.termalization) return true;
   else {
      for (int j=0; j<3; j++) {
	 for (int i=0; i<param.numOfIntervals; i++) sumMeanIntegral += meanIntegral[j][i];

	 for (int i=0; i<param.numOfIntervals; i++) {
	    if (i == 0) meanIntegralCopy[i] = (7*meanIntegral[j][i]+meanIntegral[j][i+1])/8;
	    else if (i == param.numOfIntervals-1) meanIntegralCopy[i] = (meanIntegral[j][i-1]+7*meanIntegral[j][i])/8;
	    else meanIntegralCopy[i] = (meanIntegral[j][i-1]+6*meanIntegral[j][i]+meanIntegral[j][i+1])/8;

	    meanIntegralCopy[i] /= sumMeanIntegral;
	 }

	 if (countIter > param.noMoreRefinement) convRate = 0;
	 for (int i=0; i<param.numOfIntervals; i++) meanIntegral[j][i] = pow((1-meanIntegralCopy[i])/log(1./meanIntegralCopy[i]),convRate);
      }

      return false;
   }
}
void VegasIntegrator::AlteratesIncrements(vector<vector<double> > &interval,vector<vector<double> > meanIntegral) {
   int k;
   double deltaMeanIntegral,accumMeanIntegral,deltaXk,sumMeanIntegral;
   vector<vector<double> > intervalCopy(interval);

   for (int j=0; j<3; j++) {
      sumMeanIntegral = 0;

      for (int i=0; i<param.numOfIntervals; i++) sumMeanIntegral += meanIntegral[j][i];

      accumMeanIntegral = k = 0;
      deltaMeanIntegral = sumMeanIntegral/param.numOfIntervals;
      for (int i=1; i<param.numOfIntervals; i++) {
	 while (accumMeanIntegral < deltaMeanIntegral) {
	    accumMeanIntegral += meanIntegral[j][k];
	    k++;
	 }
	 accumMeanIntegral -= deltaMeanIntegral;
	 deltaXk = interval[j][k]-interval[j][k-1];

	 intervalCopy[j][i] = interval[j][k]-accumMeanIntegral/meanIntegral[j][k-1]*deltaXk;
	 // cout << k << " ";
      }
      // cout << endl;
   }
   interval = intervalCopy;
}
double VegasIntegrator::ChiSquare(double integralPerIteration,double variancePerIteration,double estimatedIntegral) {
   double chiSq;

   chiSq = Power(integralPerIteration-estimatedIntegral,2)/variancePerIteration;
   // chi_sq = Power(integralPerIteration-estimatedIntegral,2)*Power(integralPerIteration,2)/(Power(estimatedIntegral,2)*variancePerIteration);

   return chiSq;
}
double VegasIntegrator::VariancePerIteration(double sample,double squareSample) {
   return (squareSample-sample*sample)/param.numOfPoints;
}
double VegasIntegrator::ComputesRelativeError(double analyticIntegral,double estimatedIntegral) {
   return fabs((analyticIntegral-estimatedIntegral)/analyticIntegral)*100;
}
double VegasIntegrator::RelativeError(void) {
   return ComputesRelativeError(param.analyticInt,integral);
}
double VegasIntegrator::Power(double x,int n) {
   double y;
   if (n == 0) return 1;
   else {
      y = x;
      for (int i=1; i<n; i++) y *= x;
      return y;
   }
}
