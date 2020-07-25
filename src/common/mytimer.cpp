#ifndef _MYTIMER_CPP_
#define _MYTIMER_CPP_

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <iomanip>
#include "mytimer.h"

/* ************************************************************************** */
MyTimer::MyTimer() {
   
}
/* ************************************************************************** */
void MyTimer::Start(void) {
   gettimeofday(&start,NULL);
}
/* ************************************************************************** */
void MyTimer::End(void) {
   gettimeofday(&end,NULL);
}
/* ************************************************************************** */
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
/* ************************************************************************** */
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
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */


#endif  /* _MYTIMER_CPP_ */

