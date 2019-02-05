/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.3.1
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



#ifndef _COL_SCHEME_JMOL_H_
#define _COL_SCHEME_JMOL_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#define _HAVE_SELECTED_ATOM_PALETTE_ 1

solreal getAtomicRGBColorReal(int nat,int rgb);
solreal getAtomicRColorReal(int nat);
solreal getAtomicGColorReal(int nat);
solreal getAtomicBColorReal(int nat);

int getAtomicRGBColorInt(int nat,int rgb);
int getAtomicRColorInt(int nat);
int getAtomicGColorInt(int nat);
int getAtomicBColorInt(int nat);

void getAtomicRGBColorsReal(int nat,solreal &rr,solreal &gg,solreal &bb);
void getAtomicRGBColorsInt(int nat,int &rr,int &gg,int &bb);


#endif /* defined(_COL_SCHEME_JMOL_H_) */
