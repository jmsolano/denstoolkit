#include "vegasinputparameters.h"
#include <iostream>
#include <string>
#include <climits>

using namespace std;

VegasInputParameters::VegasInputParameters(void){
   for (int j=0; j<3; j++){
      xmin[j] = 0; 
      xmax[j] = 1;
      width[j] = xmax[j]-xmin[j];
   }

   param.iterations = 20; //Maximum number of iterations to compute the integral.
   param.num_of_intervals = 10; //Number of intervals for the mesh.
   param.num_of_points = 1000; //Number of points to integer.
   param.analytic_int = 0; //Analytic integral. If not established, programe doesn't compare both the analytic one and the numeric one.
   param.relative_error = false; //Print relative error?
   param.convergence_rate = 1.0;
   param.termalization = 0;
   param.tolerance = 100;
   param.no_more_refinement = INT_MAX;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VegasInputParameters::AnalyticIntegral(double analytic_result){
   param.analytic_int = analytic_result;
   param.relative_error = true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VegasInputParameters::SetDimensions(double xleft,double yleft,double zleft,double xright,double yright,double zright){
   if (xleft > xright || yleft > yright || zleft > zright){
      cout << "Dimensions added in wrong order." << endl;
      exit (1);
   }

   xmin[0] = xleft;
   xmin[1] = yleft;
   xmin[2] = zleft;

   xmax[0] = xright;
   xmax[1] = yright;
   xmax[2] = zright;

   for (int j=0; j<3; j++) width[j] = xmax[j]-xmin[j];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VegasInputParameters::SetDimensions(double xleft,double yleft,double xright,double yright){
   if (xleft > xright || yleft > yright){
      cout << "Dimensions added in wrong order." << endl;
      exit (1);
   }

   xmin[0] = xleft;
   xmin[1] = yleft;

   xmax[0] = xright;
   xmax[1] = yright;

   for (int j=0; j<2; j++) width[j] = xmax[j]-xmin[j];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VegasInputParameters::SetDimensions(double xleft,double xright){
   if (xleft > xright){
      cout << "Dimensions added in wrong order." << endl;
      exit (1);
   }

   xmin[0] = xleft;
   xmax[0] = xright;
   width[0] = xmax[0]-xmin[0];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Imprime guía de ayuda.

void VegasInputParameters::Help(void){
   cout << "\nUso: ./lasvegas [<guía>]/[<parametros obligatorios>] [<parametros complementarios>]" << endl;
   cout << endl;
   cout << "Parametros obligatorios:" << endl;
   cout << "En los primeros seis argumentos se indica la region de integracion." << endl;
   cout << "El orden es el siguiente:" << endl;
   cout << "xmin ymin zmin xmax ymax zmax" << endl;
   cout << endl;
   cout << "Opciones complementarias:" << endl;
   cout << "-np\tNumero de puntos para integrar (10,000 por ejemplo). Por defecto son 10,000." << endl;
   cout << "-it\tMaximo numero de iteraciones a realizar (10,000 por ejemplo). Por defecto son el equivalente a MAX_INT." << endl;
   cout << "-inc\tNumero de intervalos para componer el mallado (75 por ejemplo; recomendado entre 50 y 100). Por defecto son 50." << endl;
   cout << "-rate\tRadio de convergencia de la integral (1.5 por ejemplo; recomendado entre 1.0 y 2.0). Por defecto es 2.0." << endl;
   cout << "-aint\tValor de integral analitica (para obtener error relativo)." << endl;
   cout << "-pvar\tImprimir varianza de la integral." << endl;
   cout << "-p\tImprimir parametros de entrada." << endl;
   cout << endl;
   cout << "La funcion a integrar se agrega en el archivo 'function.cpp'." << endl;
   cout << endl;
   cout << "Para abrir la guia, coloque -h.\n" << endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Imprime los datos de entrada aportados por la línea de comando.

void VegasInputParameters::DisplayProperties(void){
   printf("\nLeft_limit: (%lf,%lf,%lf)",xmin[0],xmin[1],xmin[2]);
   printf("\nRight_limit: (%lf,%lf,%lf)",xmax[0],xmax[1],xmax[2]);
   printf("\nConvergence_rate: %lf",param.convergence_rate);
   printf("\nNumber_of_intervals: %d",param.num_of_intervals);
   printf("\nNumber_of_points_to_sample: %ld",param.num_of_points);
   printf("\nMax_number_of_iterations: %ld",param.iterations);
   if (param.relative_error == true) printf("\nAnalytic_integral: %lf",param.analytic_int);
   cout << "\n" << endl;
}
