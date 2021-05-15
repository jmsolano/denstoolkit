#ifndef _INPUTPARAMETERS_H
#define _INPUTPARAMETERS_H

class VegasInputParameters{
   protected:
   struct Parameters{
      int num_of_intervals;
      long int iterations,num_of_points,termalization,no_more_refinement;
      double analytic_int,convergence_rate,tolerance;
      bool relative_error,print,print_var,print_num_points;
   }param;
		
   double xmin[3],xmax[3],width[3];
   
   void Help(void);

   public:
   VegasInputParameters(void);
   void SetDimensions(double xleft,double yleft,double zleft,double xright,double yright,double zright);
   void SetDimensions(double xleft,double yleft,double xright,double yright);
   void SetDimensions(double xleft,double xright);
   void SetIterations(double it){param.iterations = it;}
   void SetIncrements(int intervals){param.num_of_intervals = intervals;}
   void SetConvergenceRate(double conv){param.convergence_rate = conv;}
   void SetNumOfPoints(double npoints){param.num_of_points = npoints;}
   void SetTermalization(double terma){param.termalization = terma;}
   void SetTolerance(double tol){param.tolerance = tol;}
   void AnalyticIntegral(double analytic_result);
   void DisplayProperties(void);
};

#endif
