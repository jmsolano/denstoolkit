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
#ifndef _BASESPHTDESIGN_H_
#define _BASESPHTDESIGN_H_
#include <vector>
using std::vector;
#include <string>
using std::string;
#define SPHERICAL_T_DESIGN_MAX_N_RULES_DEFINED 60

/* ************************************************************************** */
/** This class provides abscissas and weigths of symmetrical spherical t-designs.
   *  Tabulated values can be found in
   *  https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html
   *  Some information and details can be found in
   *  ---> Beentjes, C. H. L., Quadrature on a Spherical Surface, Mathematical
   *         Institute, University of Oxford, UK, 2015. */
class BaseSphericalTDesign {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   static vector<vector<double> > GetAbscissas(int nn);
   static vector<double> GetWeights(int nn);
   static const double* GetBaseAbscissas(int nn);
   /** Returns the size of the internal data array. This is not to be confused
    * with the number of abscissas/weigths.  */
   static size_t GetSizeArray(int nn);
   static string GetAvailableOrders();
   static string GetAvailableRules();
   static string GetAvailablePrecisions();
   static int order_table (int rule);
   static int precision_table(int rule);
   static bool available_table (int rule);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
/*          read !cat ../gentsphdess/DIMBASESPHTDESDEFS                       */
   constexpr static size_t dim_ss003_00006=3;
   constexpr static size_t dim_ss005_00012=6;
   constexpr static size_t dim_ss007_00032=16;
   constexpr static size_t dim_ss009_00048=24;
   constexpr static size_t dim_ss011_00070=35;
   constexpr static size_t dim_ss013_00094=47;
   constexpr static size_t dim_ss015_00120=60;
   constexpr static size_t dim_ss017_00156=78;
   constexpr static size_t dim_ss019_00192=96;
   constexpr static size_t dim_ss021_00234=117;
   constexpr static size_t dim_ss031_00498=249;
   constexpr static size_t dim_ss039_00782=391;
   constexpr static size_t dim_ss045_01038=519;
   constexpr static size_t dim_ss055_01542=771;
   constexpr static size_t dim_ss063_02018=1009;
   constexpr static size_t dim_ss071_02558=1279;
   constexpr static size_t dim_ss077_03006=1503;
   constexpr static size_t dim_ss089_04008=2004;
   constexpr static size_t dim_ss101_05154=2577;
   constexpr static size_t dim_ss111_06218=3109;
/* ************************************************************************** */
/*          read !cat ../gentsphdess/PTRBASESPHTDESDEFS                       */
   const static double ss003_00006[9];
   const static double ss005_00012[18];
   const static double ss007_00032[48];
   const static double ss009_00048[72];
   const static double ss011_00070[105];
   const static double ss013_00094[141];
   const static double ss015_00120[180];
   const static double ss017_00156[234];
   const static double ss019_00192[288];
   const static double ss021_00234[351];
   const static double ss031_00498[747];
   const static double ss039_00782[1173];
   const static double ss045_01038[1557];
   const static double ss055_01542[2313];
   const static double ss063_02018[3027];
   const static double ss071_02558[3837];
   const static double ss077_03006[4509];
   const static double ss089_04008[6012];
   const static double ss101_05154[7731];
   const static double ss111_06218[9327];
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _BASESPHTDESIGN_H_ */

