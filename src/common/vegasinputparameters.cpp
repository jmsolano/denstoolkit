#include <iostream>
using std::cout;
using std::endl;
#include <climits>
#include <string>
using std::string;

#include "vegasinputparameters.h"


VegasInputParameters::VegasInputParameters() {
   for (int j=0; j<3; j++){
      xMin[j] = 0; 
      xMax[j] = 1;
      width[j] = xMax[j]-xMin[j];
   }

   param.function = "rho";
   param.iterations = 20;
   param.numOfIntervals = 10;
   param.numOfPoints = 1000;
   param.analyticInt = 0;
   param.relativeError = false;
   param.convergenceRate = 1.0;
   param.termalization = 0;
   param.tolerance = 100;
   param.noMoreRefinement = INT_MAX;
}
void VegasInputParameters::AnalyticIntegral(double analyticResult) {
   param.analyticInt = analyticResult;
   param.relativeError = true;
}
void VegasInputParameters::SetDimensions(double xLeft,double yLeft,double zLeft,double xRight,double yRight,double zRight) {
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
void VegasInputParameters::DisplayProperties(void) {
   printf("\nLeft limit: (%lf,%lf,%lf)",xMin[0],xMin[1],xMin[2]);
   printf("\nRight limit: (%lf,%lf,%lf)",xMax[0],xMax[1],xMax[2]);
   cout << "\nIntegrand: " << param.function;
   printf("\nConvergence rate: %lf",param.convergenceRate);
   printf("\nNumber of intervals: %d",param.numOfIntervals);
   printf("\nNumber of points to sample: %ld",param.numOfPoints);
   printf("\nMax number of iterations: %ld",param.iterations);
   printf("\nTermalization: %ld",param.termalization);
   printf("\nStop refinement after iteration: %ld",param.noMoreRefinement);
   printf("\nTolerance: %lf",param.tolerance);
   if (param.relativeError == true) printf("\nAnalytic integral: %lf",param.analyticInt);
   cout << "\n" << endl;
}
