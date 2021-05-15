#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "vegasintegrator.h"
#include "gausswavefunction.h"

using namespace std;
/* ********************************************************************************************************************* */
VegasIntegrator::VegasIntegrator(){
    wf=nullptr;
    integral=0.0e0;
    count_iter=count_eval=0;
    stop_iterating = false;
}
/* ********************************************************************************************************************* */
VegasIntegrator::VegasIntegrator(GaussWaveFunction &uwf) : VegasIntegrator(){
    wf = &uwf;
}
/* ********************************************************************************************************************* */
void VegasIntegrator::Integrate(void){
    vector<vector<double> > interval(3,vector<double>(param.num_of_intervals+1));
    vector<vector<double> > d_i(3,vector<double>(param.num_of_intervals));

    while (repeat_integral == true || (stop_iterating == false && count_iter < param.iterations)){
	weighted_average=inverse_variance=chisquare=0;
	count_iter=0;
	repeat_integral = stop_iterating = false;
	for (int j=0; j<3; j++){
	    for (int i=0; i<param.num_of_intervals+1; i++) interval[j][i] = xmin[j]+i*1./param.num_of_intervals*width[j];
	}

	while (repeat_integral == false && stop_iterating == false && count_iter < param.iterations){
	    count_iter++;
	    // cout << "Iteration: " << count_iter << endl;

	    MonteCarloIntegration(interval,d_i);
	    // cout << sampling.simple << " " << sampling.square << endl;

	    var_per_it = VariancePerIteration(sampling.simple,sampling.square);

	    if (count_iter > param.termalization){
		weighted_average += sampling.simple/var_per_it;
		inverse_variance += 1./var_per_it;

		integral = weighted_average/inverse_variance;
		variance = 1./inverse_variance;
	    }

	    // cout << "sampling= " << sampling.simple << endl;
	    // cout << "integral= " << integral << endl;

	    stop_iterating = AlteratesAverageValue(d_i);
	    AlteratesIncrements(interval,d_i);

	    if (count_iter > param.termalization){
		chisquare += ChiSquare(sampling.simple,var_per_it,integral);
		if (chisquare > count_iter){
		    // cout << "Chi Square" << endl;
		    // getchar();
		    repeat_integral = true;
		}
	    }
	}
    }

    count_eval = param.num_of_points*count_iter;
}
/* ********************************************************************************************************************* */
void VegasIntegrator::MonteCarloIntegration(vector<vector<double> > interval,vector<vector<double> >& d_i){
    int floor_yng[3];
    double decimal_yng,yng,delta_x_i,jacobian_i;
    double x_i[3];
    vector<vector<int> > count_points_sampled(3,vector<int>(param.num_of_intervals,0));

    sampling.simple = sampling.square = 0;
    for (int j=0; j<3; j++){
	for (int i=0; i<param.num_of_intervals; i++) d_i[j][i] = 0;
    }

    for (long int k=0; k<param.num_of_points; k++){
	jacobian_i = 1;
	for (int j=0; j<3; j++){
	    yng = dis(engine)*param.num_of_intervals;
	    floor_yng[j] = int(yng);
	    decimal_yng = yng-floor_yng[j];
	    
	    delta_x_i = interval[j][floor_yng[j]+1]-interval[j][floor_yng[j]];
	    x_i[j] = interval[j][floor_yng[j]]+delta_x_i*decimal_yng;

	    jacobian_i *= delta_x_i*param.num_of_intervals;

	    count_points_sampled[j][floor_yng[j]]++;
	}

	sampling.simple += wf->EvalDensity(x_i[0],x_i[1],x_i[2])*jacobian_i;
	sampling.square += Power(wf->EvalDensity(x_i[0],x_i[1],x_i[2])*jacobian_i,2);
	for (int j=0; j<3; j++) d_i[j][floor_yng[j]] += Power(wf->EvalDensity(x_i[0],x_i[1],x_i[2])*jacobian_i,2);
    }

    sampling.simple /= param.num_of_points;
    sampling.square /= param.num_of_points;

    for (int j=0; j<3; j++){
	for (int i=0; i<param.num_of_intervals; i++){
	    if (count_points_sampled[j][i] > 0) d_i[j][i] /= count_points_sampled[j][i];
	    // if (count_points_sampled[j][i] > 0) d_i[j][i] /= param.num_of_points*1./param.num_of_intervals;
	    // if (d_i[j][i] == 0) d_i[j][i] = numeric_limits<double>::min();
	}
    }

    // cout << "points sampled per interval: ";
    // for (int j=0; j<3; j++){
	// for (int i=0; i<param.num_of_intervals; i++) cout << count_points_sampled[j][i] << " ";
	// cout << endl;
    // }
}
/* ********************************************************************************************************************* */
bool VegasIntegrator::AlteratesAverageValue(vector<vector<double> >& d_i){
    double sum_inferior=0,sum_superior=0,conv_rate=param.convergence_rate;
    double sum_d_i=0;
    vector<double> d_i_copy(param.num_of_intervals);

    for (int j=0; j<3; j++){
	sum_inferior += *min_element(d_i[j].begin(),d_i[j].end());
	sum_superior += *max_element(d_i[j].begin(),d_i[j].end());

	// cout << *min_element(d_i[j].begin(),d_i[j].end()) << " " << *max_element(d_i[j].begin(),d_i[j].end()) << endl;
    }
    sum_inferior /= 3;
    sum_superior /= 3;

    if (fabs(sum_superior-sum_inferior) <= param.tolerance && count_iter > param.termalization) return true;
    else{
	for (int j=0; j<3; j++){
	    for (int i=0; i<param.num_of_intervals; i++) sum_d_i += d_i[j][i];

	    for (int i=0; i<param.num_of_intervals; i++){
		if (i == 0) d_i_copy[i] = (7*d_i[j][i]+d_i[j][i+1])/8;
		else if (i == param.num_of_intervals-1) d_i_copy[i] = (d_i[j][i-1]+7*d_i[j][i])/8;
		else d_i_copy[i] = (d_i[j][i-1]+6*d_i[j][i]+d_i[j][i+1])/8;

		d_i_copy[i] /= sum_d_i;
	    }

	    if (count_iter > param.no_more_refinement) conv_rate = 0;
	    for (int i=0; i<param.num_of_intervals; i++) d_i[j][i] = pow((1-d_i_copy[i])/log(1./d_i_copy[i]),conv_rate);
	}

	return false;
    }
}
/* ********************************************************************************************************************* */
void VegasIntegrator::AlteratesIncrements(vector<vector<double> >& interval,vector<vector<double> > d_i){
    int k;
    double delta_d,s_d,delta_x_k,sum_d_i;
    vector<vector<double> > interval_copy(interval);

    for (int j=0; j<3; j++){
	sum_d_i = 0;

	for (int i=0; i<param.num_of_intervals; i++) sum_d_i += d_i[j][i];
	
	s_d = k = 0;
	delta_d = sum_d_i/param.num_of_intervals;
	for (int i=1; i<param.num_of_intervals; i++){
	    while(s_d < delta_d){
		s_d += d_i[j][k];
		k++;
	    }
	    s_d -= delta_d;
	    delta_x_k = interval[j][k]-interval[j][k-1];

	    interval_copy[j][i] = interval[j][k]-s_d/d_i[j][k-1]*delta_x_k;
	    // cout << k << " ";
	}
	// cout << endl;
    }
    interval = interval_copy;
}
/* ********************************************************************************************************************* */
double VegasIntegrator::ChiSquare(double integral_per_iteration,double variance_per_iteration,double estimated_integral){
    double chi_sq;

    chi_sq = Power(integral_per_iteration-estimated_integral,2)/variance_per_iteration;
    // chi_sq = Power(integral_per_iteration-estimated_integral,2)*Power(integral_per_iteration,2)/(Power(estimated_integral,2)*variance_per_iteration);

    return chi_sq;
}
/* ********************************************************************************************************************* */
double VegasIntegrator::VariancePerIteration(double sample,double square_sample){
    return (square_sample-sample*sample)/param.num_of_points;
}
/* ********************************************************************************************************************* */
double VegasIntegrator::ComputesRelativeError(double analytic_integral,double estimated_integral){
    return fabs((analytic_integral-estimated_integral)/analytic_integral)*100;
}
/* ********************************************************************************************************************* */
double VegasIntegrator::RelativeError(void){
    return ComputesRelativeError(param.analytic_int,integral);
}
/* ********************************************************************************************************************* */
double VegasIntegrator::Power(double x,int n){
    double y;
    if (n == 0) return 1;
    else{
	y = x;
	for (int i=1; i<n; i++) y *= x;
	return y;
    }
}
