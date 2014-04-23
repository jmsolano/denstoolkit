/*
 *  solmemhand.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 25/03/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _SOLMEMHANDLE_H_
#define _SOLMEMHANDLE_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
#endif

#include <iostream>

#include <string>
using std::string;

//**************************************************************************************************
//**************************************************************************************************
bool alloc1DRealArray(string ptrname,const int n,solreal* &thptr);
//**************************************************************************************************
bool alloc1DRealArray(string ptrname,const int n,solreal* &thptr,const solreal inval);
//**************************************************************************************************
bool dealloc1DRealArray(solreal* &tp);
//**************************************************************************************************
bool alloc1DIntArray(string ptrname,const int n,int* &thptr,const int inval);
//**************************************************************************************************
bool alloc1DIntArray(string ptrname,const int n,int* &thptr);
//**************************************************************************************************
bool dealloc1DIntArray(int* &tp);
//**************************************************************************************************
bool alloc1DStringArray(string ptrname,const int n, string* &thptr);
//**************************************************************************************************
bool dealloc1DStringArray(string* &tp);
//**************************************************************************************************
bool alloc1DBoolArray(string ptrname,const int n,bool* &thptr);
//**************************************************************************************************
bool alloc1DBoolArray(string ptrname,const int n,bool* &thptr,const bool inval);
//**************************************************************************************************
bool dealloc1DBoolArray(bool* &tp);
//**************************************************************************************************
bool alloc2DRealArray(string ptrname,const int rows,const int cols,solreal** &thptr);
//**************************************************************************************************
bool alloc2DRealArray(string ptrname,const int rows,const int cols,solreal** &thptr,const solreal inval);
//**************************************************************************************************
bool dealloc2DRealArray(solreal** & tp,const int nc);
//**************************************************************************************************
bool alloc3DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,solreal*** &thptr);
//**************************************************************************************************
bool dealloc3DRealArray(solreal*** & tp,const int idx1,const int idx2);
//**************************************************************************************************
bool alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr,const int val);
//**************************************************************************************************
bool alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr);
//**************************************************************************************************
bool dealloc2DIntArray(int** & tp,const int nc);
//**************************************************************************************************
#endif//_SOLMEMHANDLE_H_



