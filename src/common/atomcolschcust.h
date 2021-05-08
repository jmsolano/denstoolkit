/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
#ifndef _COL_SCHEME_JMHP_H_
#define _COL_SCHEME_JMHP_H_

#define _HAVE_SELECTED_ATOM_PALETTE_ 1
#define MAXDEFINEDATOMICCOLORS 19

/** This table contains the colors for the atoms. The numbers are RGB integers (\f$0-255\f$) */
static int atomicColorInt[MAXDEFINEDATOMICCOLORS][3]={
   {255, 255, 255},
   {255, 192, 203},
   {178, 34, 34},
   {0, 255, 0},
   {211, 211, 211},
   {143, 143, 255},
   {255, 0, 0},
   {218, 165, 32},
   {0, 0, 255},
   {34, 139, 34},
   {0, 0, 255},
   {34, 139, 34},
   {128, 128, 144},
   {218, 165, 32},
   {255, 165, 0},
   {255, 200, 50},
   {0, 255, 0},
   {255, 20, 147},
   {255, 20, 147}
};

/** This table contains the colors for the atoms. The numbers are RGB floats (\f$0.0-1.0\f$) */
static double atomicColor[MAXDEFINEDATOMICCOLORS][3]={
   {1.00000, 1.00000, 1.00000},
   {1.00000, 0.75294, 0.79608},
   {0.69804, 0.13333, 0.13333},
   {0.00000, 1.00000, 0.00000},
   {0.82745, 0.82745, 0.82745},
   {0.56078, 0.56078, 1.00000},
   {1.00000, 0.00000, 0.00000},
   {0.85490, 0.64706, 0.12549},
   {0.00000, 0.00000, 1.00000},
   {0.13333, 0.54510, 0.13333},
   {0.00000, 0.00000, 1.00000},
   {0.13333, 0.54510, 0.13333},
   {0.50196, 0.50196, 0.56471},
   {0.85490, 0.64706, 0.12549},
   {1.00000, 0.64706, 0.00000},
   {1.00000, 0.78431, 0.19608},
   {0.00000, 1.00000, 0.00000},
   {1.00000, 0.07843, 0.57647},
   {1.00000, 0.07843, 0.57647}
};

#endif /* defined(_COL_SCHEME_JMHP_H_) */
