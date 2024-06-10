/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/
#ifndef _INTEGRATOR3D_VEGAS_H_
#define _INTEGRATOR3D_VEGAS_H_
#include <vector>
using std::vector;
#include <random>
using std::random_device;
using std::mt19937_64;
using std::uniform_real_distribution;
#include <cmath>
#include <string>
using std::string;
#include "gausswavefunction.h"
#include "fldtypesdef.h"
#include "bondnetwork.h"
#include "integrator3d.h"

/* *********************************************************************************** */
class Integrator3DVegas : public Integrator3D {
   /* *********************************************************************************** */
public:
   /* *********************************************************************************** */
   Integrator3DVegas();
   Integrator3DVegas(shared_ptr<Function3D> i);
   Integrator3DVegas(GaussWaveFunction &uwf,BondNetWork &ubnw);
   /* *********************************************************************************** */
   /** Computes integral through Las Vegas method. */
   void ComputeIntegral(void);
   /** Shows the integral (expected value) once Las Vegas integration has ended. */
   double Integral(void);
   /** Displays Las Vegas method's parameters: Integration region, integrand, convergence rate, grid size 
    * (number of intervals), number of Monte Carlo points, maximum number of iterations, thermalization
    * (num. of iterations to ignore statistical computing), num. of iterations before stopping refinement 
    * and refinement tolerance. */
   void DisplayProperties(void);
   /** Displays the results of the integration, as well as additional information  */
   void DisplayResults();
   size_t NumberOfEvaluations() { return size_t(countEval);}
   /** Sets the integration region limits for a function \f$f\f$ of type \f$f:R^3 -> R\f$. */
   void SetDimensions(double xLeft,double yLeft,double zLeft,double xRight,double yRight,double zRight);
   /** Sets the iterations limit. */
   void SetIterations(double it){param.iterations = it;}
   /** Sets the grid size (number of intervals for the integration grid). */
   void SetIntervals(int intervals){param.numOfIntervals = intervals;}
   /** Sets the convergence rate (damping parameter).
    * If one doesn't know its optimal value, one can set it to 1.
    * To set it to 0 means there is no refinement.
    * This parameter helps Las Vegas method to smooth the grid refinement each iteration. However,
    * its value, set typically between 0 and 2, depends a lot on the integrand and it is
    * established after trials and errors. */
   void SetConvergenceRate(double conv){param.convergenceRate = conv;}
   /** Number of Monte Carlo points to sample randomly over the integrand each iteration. */
   void SetNumOfPoints(double nPoints){param.numOfPoints = nPoints;}
   /** Sets the number of iterations where the cumulative expected value (integral) is not computed.
    * During the first iterations, Las Vegas Method receives few information about the integrand, causing
    * possible poor expected values. In order not to bias the expected value, it can be computed afer the
    * grid has been optimally refined (after 10-20 iterations). */
   void SetThermalization(double therm){param.thermalization = therm;}
   /** Sets the confidence interval limits related to the optimal grid. 
    * The reliable meassurement to use for finding the optimal grid is the average value for each
    * interval: All of them have to tend to a same equal value.
    * The tolerance established the minimum distance between the
    * maximum and minimum average values before the grid gets refined,
    * so that if both values reach the tolerance, the grid reaches its optimal form. */
   void SetTolerance(double tol){param.tolerance = tol;}
   /** Sets the number of iterations where the grid will be continually refined.
    * As the accumulated expected value might get biased, one can avoid that happening by making 
    * MC integrations without refinement (after the grid becomes refined enough), so as the 
    * accumulated expected value starts to beheave as a normal distribution. */
   void SetStopRefinement(double stopRef){param.noMoreRefinement = stopRef;}
   /** Set the internal integrand char.  */
   void SetInternalIntegrand(char ft) { param.integrand = ft; }
   /** Computes the normalization constant of the Electron Density, so that functions related with the
    * normalized Electron Density can be integrated. Note that if you ask for the Electron Density,
    * DTK will give you approximately 1. */
   void NormalizedEDF(void){if ( normConstant == 0 ) normConstant = round(wf->IntegralRho());}
   /** Computes the global maximum/maxima of the Electron Density (choice = 'g') or an average of the 
    * maxima (choice = 'a') in position space. In momentum space, its critical point is found near the 
    * origin or in it, so it is set in the origin.
    * Note that if you ask for the Electron Density, DTK will give you \f$\rho/\rho_{max}\f$. */
   //void Relative2MaxDensity(char choice);
   /** Shows the integral variance (Las Vegas method). */
   double Variance(void) {return fabs(variance);}
   /** Shows the maximum value of electron density. */
   double MaxDensity(void) {return maxDensity;}
   /** Shows the normalization constant of integrand. */
   double NormConstant(void) {return normConstant;}
   /** Shows the total number of evaluations during Las Vegas integration. */
   long int CountEvaluations(void) {return countEval;}
   /** Shows the number of iterations runned during Las Vegas integration. */
   long int CountIterations(void) {return countIter;}
   /** Sets the number of points to find the global maximum of the electron density.
    * This method only works for functions in momentum space (Shannon Entropy, Electron
    * Density and Kinnetic Energy). */
   void SetNSamplesToFindMaximum(double NSamples){param.nPointsForMax = NSamples;}
   /** Print all information related to the integral: Input data and output data. */
   void WriteProperties(ofstream &ofil);
   void WriteResults(ofstream &ofil);
   double RelativeError(void);
   void AnalyticIntegral(double analyticResult);
   /* *********************************************************************************** */
protected:
   /* *********************************************************************************** */
   // Use random_device to generate a seed for Mersenne twister engine.
   random_device rd{};
   // Use Mersenne twister engine to generate pseudo-random numbers.
   mt19937_64 engine{rd()};
   // "Filter" MT engine's output to generate pseudo-random double values,
   // **uniformly distributed** on the closed interval [0, 1).
   uniform_real_distribution<double> dis{0.0,0.9999};

   GaussWaveFunction* wf;
   BondNetWork* bnw;
   double integral,variance,maxDensity,maxMomDensity,normConstant;
   double weightedAverage,inverseVariance,chiSquare,varPerIt;
   double xMin[3],xMax[3],width[3],criticalPoint[3];
   long int countEval,countIter;
   bool stopIterating,repeatIntegral;

   struct Parameters{
      int numOfIntervals;
      long int iterations,numOfPoints,thermalization,noMoreRefinement,nPointsForMax;
      double analyticInt,convergenceRate,tolerance;
      bool relativeError,printVar,printNumPoints;
      char integrand;
   } param;
   struct Sampling{
      double simple,square;
   } sampling;
   double ChiSquare(double integralPerIteration,double variancePerIteration,double estimatedIntegral);
   double ComputesRelativeError(double analyticIntegral,double estimatedIntegral);
   double VariancePerIteration(double sample,double squareSample);
   //double Integrand(double x,double y,double z);
   void SearchForMaximum(void);
   void MonteCarloIntegration(vector<vector<double> > interval,vector<vector<double> > &meanIntegral);
   void AlteratesIncrements(vector<vector<double> > &interval,vector<vector<double> > meanIntegral);
   bool AlteratesAverageIntegral(vector<vector<double> > &meanIntegral);
};

#endif /* _INTEGRATOR3D_VEGAS_H_ */

