/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.3.0
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



#ifndef _SOLMEMHANDLE_CPP_
#define _SOLMEMHANDLE_CPP_

#include "solmemhand.h"
#include <cstring>
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
bool alloc4DRealArray(string ptrname,const int idx1,const int idx2,\
      const int idx3,const int idx4,solreal**** &thptr,solreal val)
{
   string errmsg=string("Warning: cannot allocate ")+ptrname\
                 +string(", in alloc4DRealArray(...) function.\n");
   if (!(thptr=new solreal***[idx1])) {
      std::cout << errmsg;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new solreal**[idx2])) {
            std::cout << errmsg;
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
         } else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new solreal*[idx3])) {
                  std::cout << errmsg;
                  std::cout << __FILE__ << " " << __LINE__ << std::endl;
               } else {
                  for ( int k=0 ; k<idx3 ; k++ ) {
                     if ( !(thptr[i][j][k]=new solreal[idx4]) ) {
                        std::cout << errmsg;
                        std::cout << __FILE__ << " " << __LINE__ << std::endl;
                     }
                  }
               }
            }
         }
      }
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for (int k=0; k<idx3; k++) {
               for ( int l=0 ; l<idx4 ; l++ ) {
                  thptr[i][j][k][l]=val;
               }
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
bool dealloc4DRealArray(solreal**** &tp,const int idx1,const int idx2,\
      const int idx3)
{
   if (tp!=NULL) {
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for ( int k=0 ; k<idx3 ; k++ ) {
               delete[] tp[i][j][k];
               tp[i][j][k]=NULL;
            }
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

