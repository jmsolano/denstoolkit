/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.1.2
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */



#ifndef _SOLCUBETOOLS_CPP_
#define _SOLCUBETOOLS_CPP_

#include "solcubetools.h"

//**********************************************************************************************
void writeCubeHeader(ofstream &ofil,string &t1,string &t2,int (&bdim)[3],
                solreal (&x0)[3],solreal (&dx)[3][3],int nat,solreal* (&atchrg),solreal* (&x))
{
   ofil.width(60);
   ofil.setf(std::ios::left);
   ofil << t1 << endl;
   ofil << t2 << endl;
   ofil.setf(std::ios::right);
   ofil.width(5);
   ofil << nat;
   ofil << setprecision(6);
   ofil.setf( std::ios::fixed, std:: ios::floatfield );
   for (int i=0; i<3; i++) {
      ofil.width(12);
      ofil << x0[i];
   }
   ofil << endl;
   for (int i=0; i<3; i++) {
      ofil.width(5);
      ofil << bdim[i];
      for (int j=0; j<3; j++) {
         ofil.width(12);
         ofil << dx[i][j];
      }
      ofil << endl;
   }
   for (int i=0; i<nat; i++) {
      ofil.width(5);
      ofil << int(atchrg[i]);
      ofil.width(12);
      ofil << atchrg[i];
      for (int j=0; j<3; j++) {
         ofil.width(12);
         ofil << x[3*i+j];
      }
      ofil << endl;
   }
   return;
}
//**********************************************************************************************
void writeCubeProp(ofstream &ofil,int dim,solreal* (&prop))
{
   int count;
   count=0;
   for (int k=0; k<dim; k++) {
      ofil.width(13);
      ofil.fill(' ');
      ofil << scientific << setprecision(5) << prop[k];
      count++;
      if (count==6) {
         ofil << endl;
         count=0;
         ofil.setf( std::ios::fixed, std::ios::floatfield );
      }
   }
   if (count!=0) {
      ofil << endl;
      count=0;
   }
   return;
}
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************
//**********************************************************************************************

#endif//_SOLCUBETOOLS_CPP_

