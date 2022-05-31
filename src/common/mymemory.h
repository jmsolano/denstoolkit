/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
#ifndef _MYMEMORY_H_
#define _MYMEMORY_H_

#include <cstdlib>
#include <iostream>
#include <string>
using std::string;

/* ************************************************************************** */
class MyMemory {
/* ************************************************************************** */
public:
   /* ********************************************************************************** */
   MyMemory() {}
   /* ********************************************************************************** */
   static bool Alloc1DRealArray(string ptrname,const int n,double* &thptr);
   /* ********************************************************************************** */
   static bool Alloc1DRealArray(string ptrname,const int n,double* &thptr,const double inval);
   /* ********************************************************************************** */
   static bool Dealloc1DRealArray(double* &tp);
   /* ********************************************************************************** */
   /* ********************************************************************************** */
   /* ********************************************************************************** */
   static bool Alloc1DIntArray(string ptrname,const int n,int* &thptr,const int inval);
   /* ********************************************************************************** */
   static bool Alloc1DIntArray(string ptrname,const int n,int* &thptr);
   /* ********************************************************************************** */
   static bool Dealloc1DIntArray(int* &tp);
   /* ********************************************************************************** */
   static bool Alloc1DStringArray(string ptrname,const int n, string* &thptr);
   /* ********************************************************************************** */
   static bool Dealloc1DStringArray(string* &tp);
   /* ********************************************************************************** */
   static bool Alloc1DBoolArray(string ptrname,const int n,bool* &thptr);
   /* ********************************************************************************** */
   static bool Alloc1DBoolArray(string ptrname,const int n,bool* &thptr,const bool inval);
   /* ********************************************************************************** */
   static bool Dealloc1DBoolArray(bool* &tp);
   /* ********************************************************************************** */
   /** Allocates memory for 2D arrays. @param rows is the number of rows of the matrix
    * (in mathematical notation, i.e. in the c-style arrays is actually the first index).
    * @param: cols is the number of columns of the matrix (i.e. in the c-style is the
    * second index)  */
   static bool Alloc2DRealArray(string ptrname,const int rows,const int cols,double** &thptr);
   /* ********************************************************************************** */
   static bool Alloc2DRealArray(string ptrname,const int rows,const int cols,double** &thptr,const double inval);
   /* ********************************************************************************** */
   static bool Dealloc2DRealArray(double** & tp,const int nc);
   /* ********************************************************************************** */
   static bool Alloc3DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,double*** &thptr);
   /* ********************************************************************************** */
   static bool Dealloc3DRealArray(double*** & tp,const int idx1,const int idx2);
   /* ********************************************************************************** */
   static bool Alloc4DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,
      const int idx4,double**** &thptr,const double val=0.0e0);
   /* ********************************************************************************** */
   static bool Dealloc4DRealArray(double**** & tp,const int idx1,const int idx2,const int idx3);
   /* ********************************************************************************** */
   static bool Alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr,const int val);
   /* ********************************************************************************** */
   static bool Alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr);
   /* ********************************************************************************** */
   static bool Dealloc2DIntArray(int** & tp,const int nc);
   /* ********************************************************************************** */
   static bool Alloc3DIntArray(string ptrname,const int idx1,const int idx2, const int idx3,\
         int*** &thptr,const int val=0);
   /* ************************************************************************************ */
   static bool Dealloc3DIntArray(int*** & tp,const int idx1,const int idx2);
   /* ************************************************************************************ */
   /** This function extends the array to which thptr points to, and assign
    *  theptr[n]=thenewval. @b Caution: this function @b do @b not modify n.  */
   static bool AppendTo1DRealArray(string ptrname,const int n,double* &thptr,double thenewval);
   /* ************************************************************************** */
protected:
   /* ************************************************************************** */
};
/* ************************************************************************** */
#endif//_MYMEMORY_H_



