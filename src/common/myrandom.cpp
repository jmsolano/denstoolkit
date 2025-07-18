/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
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
#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "myrandom.h"


#define MYRANDOM_IM1 2147483563
#define MYRANDOM_IM2 2147483399
#define MYRANDOM_AM (1.0/MYRANDOM_IM1)
#define MYRANDOM_IMM1 (MYRANDOM_IM1-1)
#define MYRANDOM_IA1 40014
#define MYRANDOM_IA2 40692
#define MYRANDOM_IQ1 53668
#define MYRANDOM_IQ2 52774
#define MYRANDOM_IR1 12211
#define MYRANDOM_IR2 3791
#define MYRANDOM_NTAB 32
#define MYRANDOM_NDIV (1+MYRANDOM_IMM1/MYRANDOM_NTAB)
#define MYRANDOM_EPS 1.2e-7
#define MYRANDOM_RNMX (1.0-MYRANDOM_EPS)

long MyRandom::idum=time(NULL);
MyRandom::MyRandom() {
}
MyRandom::MyRandom(long seed) : MyRandom() {
   idum=seed;
}
double MyRandom::rannr() {
   int j;
   long k;
   static long idum2=123456789;
   static long iy=0;
   static long iv[MYRANDOM_NTAB];
   double temp;
   
   if (idum <= 0) {
      if (-(idum) < 1) idum=1;
      else idum = -(idum);
      idum2=(idum);
      for (j=MYRANDOM_NTAB+7;j>=0;j--) {
         k=(idum)/MYRANDOM_IQ1;
         idum=MYRANDOM_IA1*(idum-k*MYRANDOM_IQ1)-k*MYRANDOM_IR1;
         if (idum < 0) {idum += MYRANDOM_IM1;}
         if (j < MYRANDOM_NTAB) {iv[j] = idum;}
      }
      iy=iv[0];
   }
   k=(idum)/MYRANDOM_IQ1;
   idum=MYRANDOM_IA1*(idum-k*MYRANDOM_IQ1)-k*MYRANDOM_IR1;
   if (idum < 0) idum += MYRANDOM_IM1; 
   k=idum2/MYRANDOM_IQ2;
   idum2=MYRANDOM_IA2*(idum2-k*MYRANDOM_IQ2)-k*MYRANDOM_IR2; 
   if (idum2 < 0) idum2 += MYRANDOM_IM2;
   j=iy/MYRANDOM_NDIV;
   iy=iv[j]-idum2;
   iv[j] = idum;
   if (iy < 1) iy += MYRANDOM_IMM1;
   if ((temp=double(MYRANDOM_AM*iy)) > MYRANDOM_RNMX) {
      return MYRANDOM_RNMX;
   } else {return temp;}
}
double MyRandom::gaussdev() {
   static int iset=0;
   static double gset;
   double fac,rsq,v1,v2;
   if (idum < 0) iset=0;
   if (iset == 0) {
      do {
         v1=2.0*rannr()-1.0;
         v2=2.0*rannr()-1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else {
      iset=0;
      return gset;
   }
}

