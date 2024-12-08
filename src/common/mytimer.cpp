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
#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <iomanip>
#include "mytimer.h"

MyTimer::MyTimer() {
}
void MyTimer::Start(void) {
   gettimeofday(&start,NULL);
}
void MyTimer::End(void) {
   gettimeofday(&end,NULL);
}
double MyTimer::GetElapsedTimeMilliSec() {
   double elaps=GetCPUSecond(end)-GetCPUSecond(start);
   return (elaps*1000.0e0);
}
double MyTimer::GetElapsedTimeSec() {
   return (GetCPUSecond(end)-GetCPUSecond(start));
}
void MyTimer::PrintElapsedTimeMilliSec(string msg) {
   double elaps=GetCPUSecond(end)-GetCPUSecond(start);
   elaps*=1000.0e0;
   int len=msg.length();
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[32m";                                                                                                                                                 
#endif
   std::cout << std::setprecision(6) << "Elapsed Time";
   std::cout << (len>0? " (": "")
      << msg << (len>0? ")":"") << ": " << elaps << "ms" << std::endl;
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[m";
#endif
}
void MyTimer::PrintElapsedTimeSec(string msg) {
   double elaps=GetCPUSecond(end)-GetCPUSecond(start);
   int len=msg.length();
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[32m";                                                                                                                                                 
#endif
   std::cout << std::setprecision(6) << "Elapsed Time";
   std::cout << (len>0? " (": "")
      << msg << (len>0? ")":"") << ": " << elaps << "s" << std::endl;
#if (defined(__APPLE__))||(defined __linux__)||(defined(__CYGWIN__))
   std::cout << "\033[m";
#endif
}


