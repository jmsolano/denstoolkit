#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
#include <algorithm>
#include "colorutils.h"

void ColorUtils::rgb2hsl(const double r,const double g,const double b,\
      double &h,double &s,double &l) {
   double themin,themax,delta;
   themin=std::min(r,std::min(g,b));
   themax=std::max(r,std::max(g,b));
   delta=themax-themin;
   l=0.5e0*(themin+themax);
   s=0.0e0;
   if ( (l>0.0e0)&&(l<1.0e0) ) { s=delta/(l<0.5?(2.0e0*l):(2.0e0-2.0e0*l)); }
   h=0.0e0;
   if ( delta>0 ) {
      if ( themax==r && themax!=g ) { h+=((g-b)/delta); }
      if ( themax==g && themax!=b ) { h+=(2.0e0+(b-r)/delta); }
      if ( themax==b && themax!=r ) { h+=(4.0e0+(r-g)/delta); }
      h*=60.0e0;
   }
}
void ColorUtils::hsl2rgb(const double h,const double s,const double l,\
      double &r,double &g,double &b) {
   double c1h=h,c1s=s,c1l=l;
   double satr,satg,satb;

   while ( c1h<0.0e0 ) { c1h+=360.0e0; }
   while ( c1h>360.0e0) { c1h-=360.0e0; }

   if ( c1h<120.0e0 ) {
      satr=(120.0e0-c1h)*oo60;
      satg=c1h*oo60;
      satb=0.0e0;
   } else if (c1h < 240) {
      satr=0.0e0;
      satg=(240.0e0-c1h)*oo60;
      satb=(c1h-120.0e0)*oo60;
   } else {
      satr=(c1h-240.0e0)*oo60;
      satg=0.0e0;
      satb=(360.0e0-c1h)*oo60;
   }
   satr=std::min(satr,1.0e0);
   satg=std::min(satg,1.0e0);
   satb=std::min(satb,1.0e0);

   double ctmpr=2*c1s*satr+(1-c1s);
   double ctmpg=2*c1s*satg+(1-c1s);
   double ctmpb=2*c1s*satb+(1-c1s);

   if ( c1l< 0.5e0 ) {
      r=c1l*ctmpr;
      g=c1l*ctmpg;
      b=c1l*ctmpb;
   } else {
      r=(1.0e0-c1l)*ctmpr+2.0e0*c1l-1.0e0;
      g=(1.0e0-c1l)*ctmpg+2.0e0*c1l-1.0e0;
      b=(1.0e0-c1l)*ctmpb+2.0e0*c1l-1.0e0;
   }
}
