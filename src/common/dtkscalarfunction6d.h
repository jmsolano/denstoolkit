/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
#ifndef _DTKSCALARFUNCTION6D_H_
#define _DTKSCALARFUNCTION6D_H_
#include <vector>
using std::vector;
#include "../common/fldtypesdef.h"
#include "../common/gausswavefunction.h"

/* ************************************************************************** */
class DTKScalarFunction6D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   DTKScalarFunction6D();
   DTKScalarFunction6D(GaussWaveFunction &ugwf);
   void SetScalarFunction(char t);
   inline double f(const vector<double> &x) { return f(x[0],x[1],x[2],x[3],x[4],x[5]); }
   inline double f(const double x1,const double y1,const double z1,\
                   const double x2,const double y2,const double z2) {
      switch ( sft ) {
         case ScalarField6DType::DM1P:
            return wf->EvalDensityMatrix1(x1,y1,z1,x2,y2,z2);
            break;
         case ScalarField6DType::LDM1P:
            return wf->EvalLapDensityMatrix1(x1,y1,z1,x2,y2,z2);
            break;
         case ScalarField6DType::CPDFP:
            return wf->EvalRho2ClosedShell(x1,y1,z1,x2,y2,z2);
            break;
         case ScalarField6DType::OPDFP:
            return wf->EvalRho2OpenShell(x1,y1,z1,x2,y2,z2);
            break;
         default :
            if ( prntunknfld1st ) {
               printf("\033[31mUnknown or unimplemented 6D field!\n");
               printf("%s, line: %d\033[0m\n",__FILE__,__LINE__);
               prntunknfld1st=false;
            }
            return 0.0e0;
            break;
      }
      return 0.0e0;
   }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   ScalarField6DType sft;
   GaussWaveFunction* wf;
   static bool prntunknfld1st; /*!< Whether to print error messages.  */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DTKSCALARFUNCTION6D_H_ */

