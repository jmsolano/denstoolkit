#include <cmath>
#include <vector>
using std::vector;
#include <algorithm>
using std::numeric_limits;
#include <iostream>
using std::cout;
using std::endl;
#include <climits>
#include <cstdlib>
#include <string>
using std::string;
#include <cmath>
#include "fldtypesdef.h"
#include "gausswavefunction.h"
#include "bondnetwork.h"
#include "integrator_vegas.h"


IntegratorVegas::IntegratorVegas() : Integrator() {
   countIter=countEval=0;
   stopIterating=false;
   normConstant=maxDensity=maxMomDensity=0;

   for (int j=0; j<3; j++){
      xMin[j] = 0; 
      xMax[j] = 1;
      width[j] = xMax[j]-xMin[j];

      criticalPoint[j] = 0;
   }

   param.integrand = 'd';
   param.iterations = 20;
   param.numOfIntervals = 10;
   param.numOfPoints = 1000;
   param.analyticInt = 0;
   param.relativeError = false;
   param.convergenceRate = 1.0;
   param.termalization = 0; 
   param.tolerance = 0;
   param.noMoreRefinement = INT_MAX;
   param.nPointsForMax = 100000;
}
IntegratorVegas::IntegratorVegas(GaussWaveFunction &uwf,BondNetWork &ubnw) : IntegratorVegas() {
   wf=&uwf;
   bnw=&ubnw;
}
void IntegratorVegas::AnalyticIntegral(double analyticResult) {
   param.analyticInt = analyticResult;
   param.relativeError = true;
}
void IntegratorVegas::SetDimensions(double xLeft,double yLeft,double zLeft,double xRight,double yRight,double zRight) {
   if (xLeft >= xRight || yLeft >= yRight || zLeft >= zRight) {
      cout << "Dimensions added in wrong order or limits are equal." << endl;
      exit (1);
   }

   xMin[0] = xLeft;
   xMin[1] = yLeft;
   xMin[2] = zLeft;

   xMax[0] = xRight;
   xMax[1] = yRight;
   xMax[2] = zRight;

   for (int j=0; j<3; j++) width[j] = xMax[j]-xMin[j];
}
void IntegratorVegas::DisplayProperties(void) {
   printf("\nLeft limit: (%lf,%lf,%lf)",xMin[0],xMin[1],xMin[2]);
   printf("\nRight limit: (%lf,%lf,%lf)",xMax[0],xMax[1],xMax[2]);
   cout << "\nIntegrand: " << GetFieldTypeKeyLong(param.integrand);
   if ( param.integrand == 'u' ) printf("\nNormalization constant: %s",normConstant ? std::to_string(normConstant).c_str() : "Function not normalized.");
   else printf("\nN electrons (integrated): %s",normConstant ? std::to_string(normConstant).c_str() : "No electrons calculated.");
   if ( param.integrand == 'm' || param.integrand == 'T' || param.integrand == 'k' ) {
      printf("\nDensity in (0,0,0): %lf",wf->EvalFTDensity(0,0,0)); 
      printf("\nMaximum of density (computed): %s",maxMomDensity ? std::to_string(maxMomDensity).c_str() : "There is no chosen value.");
      printf("\nOne critical point (computed): (%lf,%lf,%lf)",criticalPoint[0],criticalPoint[1],criticalPoint[2]);
      printf("\nNum of points to find global maximum: %ld",param.nPointsForMax);
   }else printf("\nMaximum of density (computed): %s",maxDensity ? std::to_string(maxDensity).c_str() : "There is no chosen value.");
   printf("\nConvergence rate: %lf",param.convergenceRate);
   printf("\nNumber of intervals: %d",param.numOfIntervals);
   printf("\nNumber of points to sample: %ld",param.numOfPoints);
   printf("\nMax number of iterations: %ld",param.iterations);
   printf("\nTermalization: %ld",param.termalization);
   printf("\nStop refinement after iteration: %ld",param.noMoreRefinement);
   printf("\nTolerance: %lf",param.tolerance);
   if (param.relativeError == true) printf("\nAnalytic integral: %lg",param.analyticInt);
   cout << "\n" << endl;
}
void IntegratorVegas::Relative2MaxDensity(char choice) {
   double evalDensity;

   if ( normConstant == 0 ) NormalizedEDF(); //Define normConstant

   if ( maxDensity == 0 ){
      if ( (param.integrand == 'm') || (param.integrand == 'T') | (param.integrand == 'k') ){
	 SearchForMaximum();
	 maxDensity = ( wf->EvalFTDensity(0,0,0) > maxMomDensity ) ? wf->EvalFTDensity(0,0,0) : maxMomDensity;
      }else{
	 switch ( choice ) {
	    case 'a': /* Average of rho max */
	       for (int i=0; i<bnw->nNuc; i++) maxDensity += wf->EvalDensity(bnw->R[i][0],bnw->R[i][1],bnw->R[i][2]);
	       maxDensity /= bnw->nNuc;
	    case 'g': /* Value of global maximum/maxima */
	       for (int i=0; i<bnw->nNuc; i++) {
		  evalDensity = wf->EvalDensity(bnw->R[i][0],bnw->R[i][1],bnw->R[i][2]);
		  if ( evalDensity  >= maxDensity ) maxDensity = evalDensity;
	       }
	 }
      }
   }
}
void IntegratorVegas::SearchForMaximum(void) {
   double functionImage;
   double radialFactor=0.5*width[0];
   double azimuthalFactor=2*M_PI+1e-4; //1e-4 because random numbers appear from 0 to 0.9999 and we want an inclusive interval.
   double polarFactor=M_PI+1e-4;
   double point[3]={dis(engine)*radialFactor,dis(engine)*azimuthalFactor,dis(engine)*polarFactor};

   for ( int i=0; i<param.nPointsForMax; i++ ) {
      functionImage = wf->EvalFTDensity(point[0]*cos(point[1])*sin(point[2]),
					point[0]*sin(point[1])*sin(point[2]),
					point[0]*cos(point[2]));
      if ( functionImage > maxMomDensity ) {
	 maxMomDensity = functionImage;
	 for ( int j=0; j<3; j++ ) criticalPoint[j] = point[j];
	 radialFactor = point[0]+1e-4;
      }
      point[0] = dis(engine)*radialFactor;
      point[1] = dis(engine)*azimuthalFactor;
      point[2] = dis(engine)*polarFactor;
   }
}
double IntegratorVegas::Integral(void) {
   if ( normConstant > 0 && maxDensity == 0 ){
      switch ( param.integrand ) {
	 case 'd' : /* Electron density (Rho)  */
	    integral = ( integral-0.5 >= int(integral) ) ? int(integral+1) : int(integral);
	    return integral*1./normConstant;
	 case 'm' : /* Electron density (Rho) in Momentum Space  */
	    integral = ( integral-0.5 >= int(integral) ) ? int(integral+1) : int(integral);
	    return integral*1./normConstant;
	 case 'S' : /* Shannon Entropy Density  */
	    return integral*1./normConstant+log(normConstant);
	 case 'T' : /* Shannon Entropy Density in Momentum Space  */
	    return integral*1./normConstant+log(normConstant);
	 default :
	    integral = ( integral-0.5 >= int(integral) ) ? int(integral+1) : int(integral);
	    return integral*1./normConstant; 
      }
   }else if ( maxDensity > 0 ){
      switch ( param.integrand ) {
	 case 'd' : /* Electron density (Rho)  */
	    return integral*1./maxDensity;
	 case 'm' : /* Electron density (Rho) in Momentum Space  */
	    return integral*1./maxDensity;
	 case 'S' : /* Shannon Entropy Density  */
	    return integral*1./normConstant+log(maxDensity);
	 case 'T' : /* Shannon Entropy Density in Momentum Space  */
	    return integral*1./normConstant+log(maxDensity);
	 default :
	    return integral*1./maxDensity; 
      }
   }

   return integral;
}
double IntegratorVegas::Integrand(double x,double y,double z) {
   switch ( param.integrand ) {
      case 'd' : /* Electron density (Rho)  */
	 return wf->EvalDensity(x,y,z);
      case 'm' : /* Electron density (Rho) in Momentum Space  */
	 return wf->EvalFTDensity(x,y,z);
      case 'g' : /* MagGradRho Density  */
	 return wf->EvalMagGradRho(x,y,z);
      case 'l' : /* Laplacian Density  */
	 return wf->EvalLapRho(x,y,z);
      case 'L' : /* LOL Density  */
	 return wf->EvalLOL(x,y,z);
      case 'E' : /* ELF Density  */
	 return wf->EvalELF(x,y,z);
      case 'S' : /* Shannon Entropy Density  */
	 return wf->EvalShannonEntropy(x,y,z);
      case 'T' : /* Shannon Entropy Density in Momentum Space  */
	 return wf->EvalMomentumShannonEntropy(x,y,z);
      case 'K' : /* Kinetic Energy Density K  */
	 return wf->EvalKineticEnergyK(x,y,z);
      case 'k' : /* Kinetic Energy Density K in Momentum Space  */
	 return wf->EvalFTKineticEnergy(x,y,z);
      case 'G' : /* Kinetic Energy Density G  */
	 return wf->EvalKineticEnergyG(x,y,z);
      case 'M' : /* MagGradLOL Density  */
	 return wf->EvalMagGradLOL(x,y,z);
      case 'V' : /* Molecular Electrostatic Potential Density  */
	 return wf->EvalMolElecPot(x,y,z);
      case 'P' : /* MagLEDVector  */
	 return wf->EvalMagLED(x,y,z);
      case 's' : /* Reduced Density Gradient  */
	 return wf->EvalReducedDensityGradient(x,y,z);
      case 'r' : /* Region of Slow Electrons  */
	 return wf->EvalRoSE(x,y,z);
      case 'v' : /* Potential Energy Density */
	 return wf->EvalVirialPotentialEnergyDensity(x,y,z);
      case 'z' : /* Non Covalent Interactions (NCI) -- Reduced Density Gradient */
	 return wf->EvalNCIs(x,y,z);
      case 'Z' : /* Non Covalent Interactions (NCI) -- Rho */
	 return wf->EvalNCILambda(x,y,z);
      case 'e' : /*!< Ellipticity  */
	 return wf->EvalEllipticity(x,y,z);
      case 'u' : /* Scalar Custom Field Density */
	 return wf->EvalCustomScalarField(x,y,z);
      default :
	 return wf->EvalDensity(x,y,z);
   }
}
void IntegratorVegas::Integrate(void) {
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
   repeatIntegral = true;
}
void IntegratorVegas::MonteCarloIntegration(vector<vector<double> > interval,vector<vector<double> > &meanIntegral) {
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

      sampling.simple += Integrand(xi[0],xi[1],xi[2])*jacobian;
      sampling.square += Power(Integrand(xi[0],xi[1],xi[2])*jacobian,2);
      for (int j=0; j<3; j++) meanIntegral[j][floorYNg[j]] += Power(Integrand(xi[0],xi[1],xi[2])*jacobian,2);
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
bool IntegratorVegas::AlteratesAverageIntegral(vector<vector<double> > &meanIntegral) {
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
void IntegratorVegas::AlteratesIncrements(vector<vector<double> > &interval,vector<vector<double> > meanIntegral) {
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
void IntegratorVegas::PrintInLogFile(string inFileName,string outFileName) {
   ofstream outFile(outFileName);
      outFile << "File read: " << inFileName
	      << "\n***********************************************************************************"
	      << "\n***********************************************************************************"
	      << "\nIntegral properties\n"
	      << "\nLeft limit: (" << xMin[0] << "," << xMin[1] << "," << xMin[2] << ")"
	      << "\nRight limit: (" << xMax[0] << "," << xMax[1] << "," << xMax[2] << ")"
	      << "\nIntegrand: " << GetFieldTypeKeyLong(param.integrand);
      if ( param.integrand == 'u' ) outFile << "\nNormalization constant: " << normConstant ? std::to_string(normConstant).c_str() : "Function not normalized";
      else outFile << "\nN electrons (integrated): " << normConstant ? std::to_string(normConstant).c_str() : "No electrons calculated.";
      outFile << "\nConvergence rate: " << param.convergenceRate;
      if ( param.integrand == 'm' || param.integrand == 'T' || param.integrand == 'k' ) {
         outFile << "\nDensity in (0,0,0): " << wf->EvalFTDensity(0,0,0);
         outFile << "\nMaximum of Density (computed): " << maxMomDensity ? std::to_string(maxMomDensity).c_str() : "There is no chosen value";
	 outFile << "\nOne critical point (computed): (" << criticalPoint[0] << "," << criticalPoint[1] << "," << criticalPoint[2] << ")"
		 << "\nNum of points to find global maximum: " << param.nPointsForMax;
      } else outFile << "\nMaximum of Density (computed): " << maxDensity ? std::to_string(maxDensity).c_str() : "There is no chosen value";
      outFile << "\nN. intervals: " << param.numOfIntervals
	      << "\nN. points to sample: " << param.numOfPoints
	      << "\nMaximum number of iterations: " << param.iterations
	      << "\nTermalization: " << param.termalization
	      << "\nStop refinement after iteration: " << param.noMoreRefinement
	      << "\nTolerance: " << param.tolerance;
      if (param.relativeError == true) outFile << "\nAnalytic integral: " << param.analyticInt;
      outFile << "\n***********************************************************************************"
	      << "\nResults\n"
	      << "\nN. integrand evaluations: " << CountEvaluations()
	      << "\nN. iterations: " << CountIterations()
	      << "\nIntegral: " << Integral();
      if (param.integrand == 'd' || param.integrand == 'm') {
	 outFile << "\nN. Electrons (Integrated): " 
	         << (Integral()-0.5 >= int(Integral()) ? int(Integral()+1) : int(Integral()));
	 // if (nelectrons > 0) cout << "Relerr(%) = " << integrator.RelativeError() << '\n';
      } 
      outFile << "\nVariance: " << Variance()
	      << "\n***********************************************************************************"
	      << "\n";
   outFile.close();
}
double IntegratorVegas::ChiSquare(double integralPerIteration,double variancePerIteration,double estimatedIntegral) {
   double chiSq;

   chiSq = Power(integralPerIteration-estimatedIntegral,2)/variancePerIteration;
   // chi_sq = Power(integralPerIteration-estimatedIntegral,2)*Power(integralPerIteration,2)/(Power(estimatedIntegral,2)*variancePerIteration);

   return chiSq;
}
double IntegratorVegas::VariancePerIteration(double sample,double squareSample) {
   return (squareSample-sample*sample)/param.numOfPoints;
}
double IntegratorVegas::ComputesRelativeError(double analyticIntegral,double estimatedIntegral) {
   return fabs((analyticIntegral-estimatedIntegral)/analyticIntegral)*100;
}
double IntegratorVegas::RelativeError(void) {
   return ComputesRelativeError(param.analyticInt,integral);
}
double IntegratorVegas::Power(double x,int n) {
   double y;
   if (n == 0) return 1;
   else {
      y = x;
      for (int i=1; i<n; i++) y *= x;
      return y;
   }
}
