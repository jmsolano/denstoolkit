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



