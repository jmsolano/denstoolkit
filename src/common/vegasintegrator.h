#ifndef _VEGASINTEGRATOR_H_
#define _VEGASINTEGRATOR_H_

#include "vegasinputparameters.h"
#include "gausswavefunction.h"
#include <vector>
#include <random>

using namespace std;

/* **************************************************************** */
class VegasIntegrator : public VegasInputParameters{
/* **************************************************************** */
   public:
/* **************************************************************** */
      VegasIntegrator();
      VegasIntegrator(GaussWaveFunction &uwf);
      double RelativeError(void);
      double Variance(void){return abs(variance);}
      double Integral(void){return integral;}
      long int CountEvaluations(void){return count_eval;}
      long int CountIterations(void){return count_iter;}
      void Integrate(void);
/* **************************************************************** */
   protected:
      // Use random_device to generate a seed for Mersenne twister engine.
      random_device rd{};
      // Use Mersenne twister engine to generate pseudo-random numbers.
      mt19937_64 engine{rd()};
      // "Filter" MT engine's output to generate pseudo-random double values,
      // **uniformly distributed** on the closed interval [0, 1).
      uniform_real_distribution<double> dis{0.0,0.9999};

      GaussWaveFunction* wf;
      double integral,variance;
      double weighted_average,inverse_variance,chisquare,var_per_it;
      long int count_eval,count_iter;
      bool stop_iterating,repeat_integral;

      struct Sampling{
	 double simple,square;
      }sampling;

      double Power(double x,int n);
      double ChiSquare(double integral_per_iteration,double variance_per_iteration,double estimated_integral);
      double ComputesRelativeError(double analytic_integral,double estimated_integral);
      double VariancePerIteration(double sample,double square_sample);
      void MonteCarloIntegration(vector<vector<double> > interval,vector<vector<double> >& d_i);
      void AlteratesIncrements(vector<vector<double> >& interval,vector<vector<double> > d_i);
      bool AlteratesAverageValue(vector<vector<double> >& d_i);
};

#endif /* _VEGASINTEGRATOR_H_ */

