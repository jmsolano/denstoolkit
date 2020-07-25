#ifndef _MYTIMER_H_
#define _MYTIMER_H_

#include <sys/time.h>
#include <string>
using std::string;
/* ************************************************************************** */
class MyTimer {
/* ************************************************************************** */
public:
   MyTimer();
   void Start(void);
   void End(void);
   void PrintElapsedTimeMilliSec(string msg="");
   void PrintElapsedTimeSec(string ms="");
/* ************************************************************************** */
protected:
   inline double GetCPUSecond(timeval &tp) {
      return (double(tp.tv_sec)+double(tp.tv_usec)*1.0e-6);
   }
   timeval start,end;
/* ************************************************************************** */
};
/* ************************************************************************** */
#endif  /* _MYTIMER_H_ */

