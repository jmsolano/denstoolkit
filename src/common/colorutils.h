#ifndef _COLORUTILS_H_
#define _COLORUTILS_H_
#include <vector>
using std::vector;
#include <cmath>

/* ************************************************************************** */
class ColorUtils {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   /* Calculate hsl from rgb
      Hue is in degrees
      Lightness is between 0 and 1
      Saturation is between 0 and 1
      This code is a refactor of the code taken from:
      http://paulbourke.net/miscellaneous/colourspace/ */
   static void rgb2hsl(const double r,const double g,const double b,\
         double &h,double &s,double &l);
   static inline double rgbDouble(int v) {return (double(v)*oo255);}
   static inline int rgbInt(const double v) { return v<=1.0e0 ? round(255.0e0*v) : 1.0e0;} 
   /** Calculate rgb from hsl, reverse of rgb2hsl(...)
      Hue is in degrees
      Lightness is between 0 and 1
      Saturation is between 0 and 1
      This code is a refactor of the code taken from:
      http://paulbourke.net/miscellaneous/colourspace/ */
   static void hsl2rgb(const double h,const double s,const double l,\
         double &r,double &g,double &b);
/* ************************************************************************** */
   static constexpr double oo255=1.0e0/255.0e0;
   static constexpr double oo60=1.0e0/60.0e0;
protected:
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _COLORUTILS_H_ */

