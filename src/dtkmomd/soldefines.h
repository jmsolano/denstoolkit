/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.0
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


#ifndef _SOLDEFINES_H_
#define _SOLDEFINES_H_


#ifndef _HAVE_DEF_REAL_TYPE_
#define _HAVE_DEF_REAL_TYPE_
typedef double real;
typedef double solreal;
//typedef float real;
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef PARALLELISEDTK
#define PARALLELISEDTK 0
#endif

#define USEPROGRESSBAR 1
#define CURRENTVERSION "1.0.0"
#define PROGRAMCONTRIBUTORS "JMSA/JMHP"
#define EPSFORELFVALUE (2.871e-05)
#define DEFAULTPOINTSPERDIRECTION (200)
#define DEFAULTNPOINTSFORCUBE (80)
#define DEFAULTMAXVALUEOFP (5.0e0)

#define DISPLAYDEBUGINFOFILELINE (std::cout << __FILE__ << ", line: " << __LINE__ << std::endl)

#if DEBUG
#define _SOL_USE_SAFE_CHECKS_ 1
#else
#define _SOL_USE_SAFE_CHECKS_ 0
#endif

#define _HAVE_GNUPLOT_ 1

#define _HAVE_EPSTOOL_ 1

#define _HAVE_EPSTOPDF_ 1

#define _SOL_USE_FIGLET_NAME_ 1

#define _MAX_MEM_ALLOWANCE_ (1024*1024*1024)

#endif /* defined(_SOLDEFINES_H_) */
