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
#include <cstdlib>
#include <cstring>
#include "mymemory.h"

bool MyMemory::Alloc1DRealArray(string ptrname,const int n,double* &thptr) {
   if (!(thptr=new double[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=0.0;}
      return true;
   }
}
bool MyMemory::Alloc1DRealArray(string ptrname,const int n,double* &thptr,const double inval) {
   if (!(thptr=new double[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
bool MyMemory::Alloc1DIntArray(string ptrname,const int n,int* &thptr) {
   if (!(thptr=new int[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=0;}
      return true;
   }
}
bool MyMemory::Alloc1DIntArray(string ptrname,const int n,int* &thptr,const int inval) {
   if (!(thptr=new int[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
bool MyMemory::Alloc1DBoolArray(string ptrname,const int n,bool* &thptr,const bool inval) {
   if (!(thptr=new bool[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DBoolArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
bool MyMemory::Alloc1DBoolArray(string ptrname,const int n,bool* &thptr) {
   return Alloc1DBoolArray(ptrname,n,thptr,false);
}
bool MyMemory::Dealloc1DBoolArray(bool* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Dealloc1DRealArray(double* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Dealloc1DIntArray(int* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Alloc1DStringArray(string ptrname,const int n, string* &thptr) {
   if (!(thptr=new string[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DStringArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=' ';}
      return true;
   }
}
bool MyMemory::Dealloc1DStringArray(string* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
bool MyMemory::Alloc2DRealArray(string ptrname,const int rows,const int cols,double** &thptr) {
   if (!(thptr=new double*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new double[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
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
bool MyMemory::Alloc2DRealArray(string ptrname,const int rows,const int cols,double** &thptr,const double inval) {
   if (!(thptr=new double*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new double[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
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
bool MyMemory::Dealloc2DRealArray(double** & tp,const int nr) {
   if (tp!=NULL) {
      for (int i=0; i<nr; i++) {
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
bool MyMemory::Alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr) {
   return Alloc2DIntArray(ptrname,rows,cols,thptr,0);
}
bool MyMemory::Alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr,const int val) {
   if (!(thptr=new int*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new int[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DIntArray(...) function.\n";
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
bool MyMemory::Dealloc2DIntArray(int** & tp,const int nr) {
   if (tp!=NULL) {
      for (int i=0; i<nr; i++) {
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Alloc3DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,double*** &thptr) {
   if (!(thptr=new double**[idx1])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc3DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new double*[idx2])) {
            std::cout << "Warning: cannot allocate "<< ptrname
            <<", in alloc3DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
         else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new double[idx3])) {
                  std::cout << "Warning: cannot allocate "<< ptrname
                  <<", in Alloc3DRealArray(...) function.\n";
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
bool MyMemory::Alloc4DRealArray(string ptrname,const int idx1,const int idx2,\
      const int idx3,const int idx4,double**** &thptr,double val) {
   string errmsg=string("Warning: cannot allocate ")+ptrname\
                 +string(", in Alloc4DRealArray(...) function.\n");
   if (!(thptr=new double***[idx1])) {
      std::cout << errmsg;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new double**[idx2])) {
            std::cout << errmsg;
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
         } else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new double*[idx3])) {
                  std::cout << errmsg;
                  std::cout << __FILE__ << " " << __LINE__ << std::endl;
               } else {
                  for ( int k=0 ; k<idx3 ; k++ ) {
                     if ( !(thptr[i][j][k]=new double[idx4]) ) {
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
bool MyMemory::Dealloc3DRealArray(double*** &tp,const int idx1,const int idx2) {
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
bool MyMemory::Dealloc4DRealArray(double**** &tp,const int idx1,const int idx2,\
      const int idx3) {
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
bool MyMemory::Alloc3DIntArray(string ptrname,const int idx1,const int idx2,const int idx3,\
      int*** &thptr,const int val) {
   if (!(thptr=new int**[idx1])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc3DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new int*[idx2])) {
            std::cout << "Warning: cannot allocate "<< ptrname
            <<", in Alloc3DIntArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
         else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new int[idx3])) {
                  std::cout << "Warning: cannot allocate "<< ptrname
                  <<", in Alloc3DIntArray(...) function.\n";
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
bool MyMemory::Dealloc3DIntArray(int*** &tp,const int idx1,const int idx2) {
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
bool MyMemory::AppendTo1DRealArray(string ptrname,const int n,double* &thptr,double thenewval) {
   double *tmpptr;
   bool res=Alloc1DRealArray("tmpptr",(n+1),tmpptr);
   if ( res ) {
      //for ( int i=0 ; i<n ; i++ ) {tmpptr[i]=thptr[i];}
      std::memcpy( tmpptr, thptr, n * sizeof(double) );
      tmpptr[n]=thenewval;
      Dealloc1DRealArray(thptr);
      thptr=tmpptr;
   } else {
      std::cout << "Error: something went wrong while trying to append a new"
         << std::endl << "value in array " << ptrname << std::endl;
   }
   return res;
}

