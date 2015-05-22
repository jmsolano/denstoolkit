/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.1.1
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



#ifndef _SOLMEMHANDLE_H_
#define _SOLMEMHANDLE_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
#endif

#include <iostream>

#include <cstring>
using std::string;

/* ********************************************************************************** */
/* ********************************************************************************** */
bool alloc1DRealArray(string ptrname,const int n,solreal* &thptr);
/* ********************************************************************************** */
bool alloc1DRealArray(string ptrname,const int n,solreal* &thptr,const solreal inval);
/* ********************************************************************************** */
bool dealloc1DRealArray(solreal* &tp);
/* ********************************************************************************** */
bool alloc1DIntArray(string ptrname,const int n,int* &thptr,const int inval);
/* ********************************************************************************** */
bool alloc1DIntArray(string ptrname,const int n,int* &thptr);
/* ********************************************************************************** */
bool dealloc1DIntArray(int* &tp);
/* ********************************************************************************** */
bool alloc1DStringArray(string ptrname,const int n, string* &thptr);
/* ********************************************************************************** */
bool dealloc1DStringArray(string* &tp);
/* ********************************************************************************** */
bool alloc1DBoolArray(string ptrname,const int n,bool* &thptr);
/* ********************************************************************************** */
bool alloc1DBoolArray(string ptrname,const int n,bool* &thptr,const bool inval);
/* ********************************************************************************** */
bool dealloc1DBoolArray(bool* &tp);
/* ********************************************************************************** */
bool alloc2DRealArray(string ptrname,const int rows,const int cols,solreal** &thptr);
/* ********************************************************************************** */
bool alloc2DRealArray(string ptrname,const int rows,const int cols,solreal** &thptr,const solreal inval);
/* ********************************************************************************** */
bool dealloc2DRealArray(solreal** & tp,const int nc);
/* ********************************************************************************** */
bool alloc3DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,solreal*** &thptr);
/* ********************************************************************************** */
bool dealloc3DRealArray(solreal*** & tp,const int idx1,const int idx2);
/* ********************************************************************************** */
bool alloc4DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,
      const int idx4,solreal**** &thptr,const solreal val=0.0e0);
/* ********************************************************************************** */
bool dealloc4DRealArray(solreal**** & tp,const int idx1,const int idx2,const int idx3);
/* ********************************************************************************** */
bool alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr,const int val);
/* ********************************************************************************** */
bool alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr);
/* ********************************************************************************** */
bool dealloc2DIntArray(int** & tp,const int nc);
/* ********************************************************************************** */
bool alloc3DIntArray(string ptrname,const int idx1,const int idx2, const int idx3,\
      int*** &thptr,const int val=0);
/* ************************************************************************************ */
bool dealloc3DIntArray(int*** & tp,const int idx1,const int idx2);
/* ************************************************************************************ */
/** This function extends the array to which thptr points to, and assign
 *  theptr[n]=thenewval. @b Caution: this function @b do @b not modify n.  */
bool appendTo1DRealArray(string ptrname,const int n,solreal* &thptr,solreal thenewval);
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
/* ************************************************************************************ */
#endif//_SOLMEMHANDLE_H_



