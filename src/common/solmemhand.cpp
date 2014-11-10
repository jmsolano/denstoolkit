/*
 *  solmemhand.cpp
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 25/03/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _SOLMEMHANDLE_CPP_
#define _SOLMEMHANDLE_CPP_

#include "solmemhand.h"
/* ************************************************************************************** */
/* ************************************************************************************** */
/* ************************************************************************************** */
bool alloc1DRealArray(string ptrname,const int n,solreal* &thptr)
{
   if (!(thptr=new solreal[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc1DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=0.0;}
      return true;
   }
}

/* ************************************************************************************** */
bool alloc1DRealArray(string ptrname,const int n,solreal* &thptr,const solreal inval)
{
   if (!(thptr=new solreal[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc1DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
/* ************************************************************************************** */
bool alloc1DIntArray(string ptrname,const int n,int* &thptr)
{
   if (!(thptr=new int[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc1DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=0;}
      return true;
   }
}
/* ************************************************************************************** */
bool alloc1DIntArray(string ptrname,const int n,int* &thptr,const int inval)
{
   if (!(thptr=new int[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc1DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
/* ************************************************************************************** */
bool alloc1DBoolArray(string ptrname,const int n,bool* &thptr,const bool inval)
{
   if (!(thptr=new bool[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc1DBoolArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
/* ************************************************************************************** */
bool alloc1DBoolArray(string ptrname,const int n,bool* &thptr)
{
   return alloc1DBoolArray(ptrname,n,thptr,false);
}
/* ************************************************************************************** */
bool dealloc1DBoolArray(bool* &tp)
{
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
/* ************************************************************************************** */
bool dealloc1DRealArray(solreal* &tp)
{
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
/* ************************************************************************************** */
bool dealloc1DIntArray(int* &tp)
{
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
/* ************************************************************************************** */
bool alloc1DStringArray(string ptrname,const int n, string* &thptr)
{
   if (!(thptr=new string[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc1DStringArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=' ';}
      return true;
   }
}
/* ************************************************************************************** */
bool dealloc1DStringArray(string* &tp)
{
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }

}
/* ************************************************************************************** */
bool alloc2DRealArray(string ptrname,const int rows,const int cols,solreal** &thptr)
{
   if (!(thptr=new solreal*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc2DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new solreal[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc2DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
      }
      for (int i=0; i<rows; i++) {
         for (int j=0; j<cols; j++) {
            thptr[i][j]=0.000000e0;
         }
      }
      return true;
   }
}
/* ************************************************************************************** */
bool alloc2DRealArray(string ptrname,const int rows,const int cols,solreal** &thptr,const solreal inval)
{
   if (!(thptr=new solreal*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc2DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new solreal[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc2DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
      }
      for (int i=0; i<rows; i++) {
         for (int j=0; j<cols; j++) {
            thptr[i][j]=inval;
         }
      }
      return true;
   }
}
/* ************************************************************************************** */
bool dealloc2DRealArray(solreal** & tp,const int nc)
{
   if (tp!=NULL) {
      for (int i=0; i<nc; i++) {
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
/* ************************************************************************************** */
bool alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr)
{
   return alloc2DIntArray(ptrname,rows,cols,thptr,0);
}
/* ************************************************************************************** */
bool alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr,const int val)
{
   if (!(thptr=new int*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc2DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new int[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc2DIntArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
      }
      for (int i=0; i<rows; i++) {
         for (int j=0; j<cols; j++) {
            thptr[i][j]=val;
         }
      }
      return true;
   }
}
/* ************************************************************************************** */
bool dealloc2DIntArray(int** & tp,const int nc)
{
   if (tp!=NULL) {
      for (int i=0; i<nc; i++) {
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
   }
   return true;
}
/* ************************************************************************************** */
bool alloc3DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,solreal*** &thptr)
{
   if (!(thptr=new solreal**[idx1])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc3DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new solreal*[idx2])) {
            std::cout << "Warning: cannot allocate "<< ptrname
            <<", in alloc3DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
         else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new solreal[idx3])) {
                  std::cout << "Warning: cannot allocate "<< ptrname
                  <<", in alloc3DRealArray(...) function.\n";
                  std::cout << __FILE__ << "" << __LINE__ << std::endl;
               }
            }
         }
      }
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for (int k=0; k<idx3; k++) {
               thptr[i][j][k]=0.0e0;
            }
         }
      }
      return true;
   }
}
/* ************************************************************************************** */
bool dealloc3DRealArray(solreal*** &tp,const int idx1,const int idx2)
{
   if (tp!=NULL) {
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            delete[] tp[i][j];
            tp[i][j]=NULL;
         }
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
/* ************************************************************************************** */
bool alloc3DIntArray(string ptrname,const int idx1,const int idx2,const int idx3,\
      int*** &thptr,const int val)
{
   if (!(thptr=new int**[idx1])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in alloc3DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new int*[idx2])) {
            std::cout << "Warning: cannot allocate "<< ptrname
            <<", in alloc3DIntArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
         else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new int[idx3])) {
                  std::cout << "Warning: cannot allocate "<< ptrname
                  <<", in alloc3DIntArray(...) function.\n";
                  std::cout << __FILE__ << "" << __LINE__ << std::endl;
               }
            }
         }
      }
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for (int k=0; k<idx3; k++) {
               thptr[i][j][k]=val;
            }
         }
      }
      return true;
   }
}
/* ************************************************************************************** */
bool dealloc3DIntArray(int*** &tp,const int idx1,const int idx2)
{
   if (tp!=NULL) {
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            delete[] tp[i][j];
            tp[i][j]=NULL;
         }
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
/* ************************************************************************************** */
bool appendTo1DRealArray(string ptrname,const int n,solreal* &thptr,solreal thenewval)
{
   solreal *tmpptr;
   bool res=alloc1DRealArray("tmpptr",(n+1),tmpptr);
   if ( res ) {
      //for ( int i=0 ; i<n ; i++ ) {tmpptr[i]=thptr[i];}
      memcpy( tmpptr, thptr, n * sizeof(solreal) );
      tmpptr[n]=thenewval;
      dealloc1DRealArray(thptr);
      thptr=tmpptr;
   } else {
      std::cout << "Error: something went wrong while trying to append a new"
         << std::endl << "value in array " << ptrname << std::endl;
   }
   return res;
}
/* ************************************************************************************** */
#endif//_SOLMEMHANDLE_CPP_

