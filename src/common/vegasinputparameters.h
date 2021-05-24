#ifndef _INPUTPARAMETERS_H
#define _INPUTPARAMETERS_H

#include <string>
using std::string;

/* *********************************************************************************** */
class VegasInputParameters {
/* *********************************************************************************** */
public:
/* *********************************************************************************** */
   VegasInputParameters();
   /* *********************************************************************************** */
   /** Sets the integration region limits for a function \f$f\f$ of type \f$f:R^3 -> R\f$. */
   void SetDimensions(double xLeft,double yLeft,double zLeft,double xRight,double yRight,double zRight);
   /**  */
   void SetFunction(string func){param.function = func;}
   /** Sets the iterations limit. */
   void SetIterations(double it){param.iterations = it;}
   /** Sets the grid size (number of intervals for the integration grid). */
   void SetIntervals(int intervals){param.numOfIntervals = intervals;}
   /**
    * Sets the convergence rate (damping parameter).
    * If one doesn't know its optimal value, one can set it to 1.
    * To set it to 0 means there is no refinement.
    * This parameter helps Las Vegas method to smooth the grid refinement each iteration. However,
    * its value, set typically between 0 and 2, depends a lot on the integrand and it is
    * established after trials and errors.
   */
   void SetConvergenceRate(double conv){param.convergenceRate = conv;}
   /** Number of Monte Carlo points to sample randomly over the integrand each iteration. */
   void SetNumOfPoints(double nPoints){param.numOfPoints = nPoints;}
   /**
    * Sets the number of iterations where the cumulative expected value (integral) is not computed.
    * During the first iterations, Las Vegas Method receives few information about the integrand, causing
    * possible poor expected values. In order not to bias the expected value, it can be computed afer the
    * grid has been optimally refined (after 10-20 iterations).
   */
   void SetTermalization(double terma){param.termalization = terma;}
   /**
    * Sets the confidence interval limits related to the optimal grid. 
    * The reliable meassurement to use for finding the optimal grid is the average value for each
    * interval: All of them have to tend to a same equal value.
    * The tolerance established the minimum distance between the maximum and minimum average values before 
    * the grid gets refined, so that if both values reach the tolerance, the grid reaches its optimal form.
   */
   void SetTolerance(double tol){param.tolerance = tol;}
   /**
    * Sets the number of iterations where the grid will be continually refined.
    * As the accumulated expected value might get biased, one can avoid that happening by making 
    * MC integrations without refinement (after the grid becomes refined enough), so as the 
    * accumulated expected value starts to beheave as a normal distribution.
   */
   void SetStopRefinement(double stopRef){param.noMoreRefinement = stopRef;}
   /**
    * Displays Las Vegas method's parameters: Integration region, convergence rate, grid size (number of
    * intervals), number of Monte Carlo points and maximum number of iterations.
   */
   void DisplayProperties(void);
   void AnalyticIntegral(double analyticResult);
/* *********************************************************************************** */
protected:
/* *********************************************************************************** */
   struct Parameters{
      int numOfIntervals;
      long int iterations,numOfPoints,termalization,noMoreRefinement;
      double analyticInt,convergenceRate,tolerance;
      bool relativeError,printVar,printNumPoints;
      string function;
   }param;
		
   double xMin[3],xMax[3],width[3];
/* *********************************************************************************** */
};
/* *********************************************************************************** */
#endif /* _H_VEGASINPUTPARAMETERS_H_ */
