

#ifndef _RADII_SCHEME_JMHP_CPP_
#define _RADII_SCHEME_JMHP_CPP_

#include "atomradiicust.h"

solreal getAtomicVDWRadius(int atn)
{
   static const solreal atomicRadius[MAXDEFINEDATOMICRADII]={
      0.37e0, 0.70e0, 1.23e0, 0.89e0, 0.80e0,
      0.79e0, 0.74e0, 0.74e0, 0.72e0, 0.70e0,
      1.57e0, 1.36e0, 1.25e0, 1.17e0, 1.10e0,
      1.04e0, 0.99e0, 0.70e0, 2.03e0, 1.74e0,
      1.44e0, 1.32e0, 1.22e0, 1.17e0, 1.16e0,
      1.16e0, 1.15e0, 1.17e0, 1.25e0, 1.25e0,
      1.22e0, 1.21e0, 1.17e0, 0.70e0, 1.24e0,
      1.91e0, 1.62e0, 1.45e0, 1.34e0, 1.29e0,
      1.29e0, 1.24e0, 1.25e0, 1.28e0, 1.34e0,
      1.41e0, 1.50e0, 1.40e0, 1.41e0, 1.37e0,
      1.33e0, 0.70e0, 1.33e0, 1.98e0, 1.69e0,
      1.69e0, 1.69e0, 1.69e0, 1.69e0, 1.69e0,
      1.69e0, 1.69e0, 1.69e0, 1.69e0, 1.69e0,
      1.69e0, 1.69e0, 1.69e0, 1.69e0, 1.69e0,
      1.69e0, 1.44e0, 1.34e0, 1.30e0, 1.28e0,
      1.26e0, 1.29e0, 1.34e0, 1.44e0, 1.55e0,
      1.54e0, 1.52e0, 1.52e0, 1.40e0, 0.70e0,
      2.40e0, 2.00e0, 1.90e0, 1.90e0, 1.90e0,
      1.90e0, 1.90e0, 0.70e0, 0.26e0, 0.80e0,
      0.80e0, 0.80e0, 0.80e0, 0.80e0, 0.80e0,
      0.80e0, 0.80e0, 0.80e0, 0.80e0, 0.80e0,
      0.80e0, 0.80e0, 0.80e0, 0.80e0, 0.80e0,
      0.80e0, 0.80e0
   };
   return atomicRadius[atn];
}

#endif /* defined(_RADII_SCHEME_JMHP_H_) */
