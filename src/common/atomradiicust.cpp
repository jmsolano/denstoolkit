/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.1
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
