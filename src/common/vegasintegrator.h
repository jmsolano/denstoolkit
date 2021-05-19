#ifndef _VEGASINTEGRATOR_H_
#define _VEGASINTEGRATOR_H_

#include <vector>
using std::vector;
#include <random>
using std::random_device;
using std::mt19937_64;
using std::uniform_real_distribution;
#include <cmath>

#include "vegasinputparameters.h"
#include "gausswavefunction.h"

/* **************************************************************** */
class VegasIntegrator : public VegasInputParameters {
/* **************************************************************** */
public:
/* **************************************************************** */
   VegasIntegrator();
   VegasIntegrator(GaussWaveFunction &uwf);
   /* **************************************************************** */
   double RelativeError(void);
   /** Shows the integral variance (Las Vegas method). */
   double Variance(void) {return fabs(variance);}
   /** Shows the integral (expected value) once Las Vegas integration has ended. */
   double Integral(void) {return integral;}
   /** Shows the total number of evaluations during Las Vegas integration. */
   long int CountEvaluations(void) {return countEval;}
   /** Shows the number of iterations runned during Las Vegas integration. */
   long int CountIterations(void) {return countIter;}
   /** Computes integral through Las Vegas method. */
   void Integrate(void);
/* **************************************************************** */
protected:
/* **************************************************************** */
   // Use random_device to generate a seed for Mersenne twister engine.
   random_device rd{};
   // Use Mersenne twister engine to generate pseudo-random numbers.
   mt19937_64 engine{rd()};
   // "Filter" MT engine's output to generate pseudo-random double values,
   // **uniformly distributed** on the closed interval [0, 1).
   uniform_real_distribution<double> dis{0.0,0.9999};

   GaussWaveFunction* wf;
   double integral,variance;
   double weightedAverage,inverseVariance,chiSquare,varPerIt;
   long int countEval,countIter;
   bool stopIterating,repeatIntegral;

   struct Sampling{
      double simple,square;
   }sampling;

   double Power(double x,int n);
   double ChiSquare(double integralPerIteration,double variancePerIteration,double estimatedIntegral);
   double ComputesRelativeError(double analyticIntegral,double estimatedIntegral);
   double VariancePerIteration(double sample,double squareSample);
   void MonteCarloIntegration(vector<vector<double> > interval,vector<vector<double> > &meanIntegral);
   void AlteratesIncrements(vector<vector<double> > &interval,vector<vector<double> > meanIntegral);
   bool AlteratesAverageIntegral(vector<vector<double> > &meanIntegral);
};

#endif /* _VEGASINTEGRATOR_H_ */

