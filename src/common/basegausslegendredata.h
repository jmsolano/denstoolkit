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
//
//  basegausslegendredata.h
//  
//
//  Created by Juan Manuel Solano on 2013-10-18.
//  Numerical values were taken from 
//  http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/
//
//  Refactored on 2016-09-19, by Juan Manuel Solano Altamirano
//  e-mail: jmsolanoalt@gmail.com
//
//
#ifndef _BASEGAUSSLEGENDREDATA_H_
#define _BASEGAUSSLEGENDREDATA_H_
#include <vector>
using std::vector;
#include <string>
using std::string;

/* ************************************************************************** */
class BaseGaussLegendreData {
/* ************************************************************************** */
public:
   BaseGaussLegendreData();
   ~BaseGaussLegendreData();
   static void GetWeightsAndAbscissas(vector<double> &ww,vector<double> &xx,int nn);
   static void GetWeightsAndAbscissas(vector<double> &ww,vector<double> &xx,\
         const double xmin,const double xmax,int nn);
   static void GetZero2InfWeightsAndAbscissas(vector<double> &ww,vector<double> &xx,\
         int nn);
   static string GetAvailableOrders();
   /** Returns the size of the internal data arrays. This is
    * not to be confused with the number of abscissas/weights.  */
   static size_t GetSizeArray(int nn);
   static const double* GetBaseAbscissas(int nn);
   static const double* GetBaseWeights(int nn);
/* Computes weights and abscissas, and it was obtained from the open source:
 * http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/ */
   static void gauss_legendre_tbl(int n, double* x, double* w, double eps);
   /** This function selects the number of points used in the quadrature, and
    * it copies the respective values into the internal arrays.  */
   bool SelectAbscissasAndWeigths(int nn);
   /** This function assign absisas and weights onto the internal arrays
    * (x and w, respectively).  */
   bool AssignValues(const double* origx,const double* origw);
   static bool AssignAbscissas(const double* orig,const int nn,vector<double> &dest);
   static bool AssignWeights(const double* orig,const int nn,vector<double> &dest);
   inline int Dim() {return dim;}
   double *x,*w;
   double *myx,*myw;
/* ************************************************************************** */
protected:
   int dim;
   bool AllocArrays();
   static vector<double> GetAbscissas(int nn);
   static vector<double> GetAbscissas(const double xmin,const double xmax,int nn);
   static vector<double> GetWeights(int nn);
   static vector<double> GetWeights(const double xmin,const double xmax,int nn);
   const static int dimBasGauLeg02,dimBasGauLeg04,dimBasGauLeg06,dimBasGauLeg08;
   const static int dimBasGauLeg11,dimBasGauLeg19;
   const static int dimBasGauLeg10,dimBasGauLeg12,dimBasGauLeg14,dimBasGauLeg16;
   const static int dimBasGauLeg18,dimBasGauLeg20,dimBasGauLeg32;
   const static double xxxBasGauLeg02[1],wwwBasGauLeg02[1];
   const static double xxxBasGauLeg04[2],wwwBasGauLeg04[2];
   const static double xxxBasGauLeg06[3],wwwBasGauLeg06[3];
   const static double xxxBasGauLeg08[4],wwwBasGauLeg08[4];
   const static double xxxBasGauLeg10[5],wwwBasGauLeg10[5];
   const static double xxxBasGauLeg11[6],wwwBasGauLeg11[6];
   const static double xxxBasGauLeg12[6],wwwBasGauLeg12[6];
   const static double xxxBasGauLeg14[7],wwwBasGauLeg14[7];
   const static double xxxBasGauLeg16[8],wwwBasGauLeg16[8];
   const static double xxxBasGauLeg18[9],wwwBasGauLeg18[9];
   const static double xxxBasGauLeg19[10],wwwBasGauLeg19[10];
   const static double xxxBasGauLeg20[10],wwwBasGauLeg20[10];
   const static double xxxBasGauLeg32[16],wwwBasGauLeg32[16];
/* ************************************************************************** */
   // n = 2
};
/* ************************************************************************** */
#endif
//_BASEGAUSSLEGENDREDATA_H_
