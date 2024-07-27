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
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <iomanip>
using std::setprecision;
#include <algorithm>
#include <cmath>
#include <cassert>
#include "screenutils.h"
#include "molecule.h"
#include "matrixvectoroperations3d.h"

Molecule::Molecule() {
   Init();
}
void Molecule::Init() {
   atom.clear();
   for ( size_t i=0 ; i<bond.size() ; ++i ) { bond[i].clear(); }
   bond.clear();
   for ( size_t i=0 ; i<bndDist.size() ; ++i ) { bndDist[i].clear(); }
   bndDist.clear();
   for ( size_t i=0 ; i<cell.size() ; ++i ) {
      for ( size_t j=0 ; j<cell[i].size() ; ++j ) {
         for ( size_t k=0 ; k<cell[i][j].size() ; ++k ) {
            cell[i][j][k].clear();
         }
         cell[i][j].clear();
      }
      cell[i].clear();
   }
   cell.clear();
   cm.resize(3);
   cd.resize(3);
   xmin.resize(3);
   xmax.resize(3);
   origCent.resize(3);
   imsetup=false;
}
Molecule::~Molecule() {
   atom.clear();
   for ( size_t i=0 ; i<bond.size() ; ++i ) { bond[i].clear(); }
   bond.clear();
   for ( size_t i=0 ; i<bndDist.size() ; ++i ) { bndDist[i].clear(); }
   bndDist.clear();
   for ( size_t i=0 ; i<cell.size() ; ++i ) {
      for ( size_t j=0 ; j<cell[i].size() ; ++j ) {
         for ( size_t k=0 ; k<cell[i][j].size() ; ++k ) {
            cell[i][j][k].clear();
         }
         cell[i][j].clear();
      }
      cell[i].clear();
   }
   cell.clear();
   cm.clear();
   cd.clear();
   origCent.clear();
   imsetup=false;
}
void Molecule::AddAtom(vector<double> &ux,int an) {
   atom.push_back(Atom(ux,an));
}
void Molecule::AddAtom(vector<double> &ux,string &usymb) {
   atom.push_back(Atom(ux,usymb));
}
void Molecule::DisplayAtomProperties() {
   size_t k=atom.size();
   for ( size_t i=0 ; i<k ; ++i ) { atom[i].DisplayProperties(); }
}
void Molecule::DisplayProperties() {
   DisplayAtomProperties();
   DisplayBondProperties();
}
void Molecule::DisplayBondProperties() {
   size_t nb=bond.size();
   cout << "There are " << nb << " bonds." << '\n';
   //*
   for ( size_t i=0 ; i<nb ; ++i ) {
      cout << "   Bonded: " << i << " (" << std::setw(2) << atom[i].symbol << "): ";
      for ( size_t j=0 ; j<bond[i].size() ; ++j ) {
         cout << bond[i][j] << '(' << std::setw(2) << atom[bond[i][j]].symbol << ") ";
      }
     // cout << setprecision(4);
     // cout << "\nDistances: " << i << " (" << std::setw(2) << atom[i].symbol << "): ";
     // for ( size_t j=0 ; j<bond[i].size() ; ++j ) {
     //    cout << std::setw(5) << bndDist[i][j] << ' ';
     // }
      cout << '\n';
   }
   // */
}
void Molecule::SortCoordinates() {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("I'm not setup! Nothing to sort.");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      return;
   }
   //Sorts first column
   QuickSort(0);
   //Sorts second column...
   size_t n=atom.size()-1;
   size_t count,currpos;
   for ( size_t i=0 ; i<n ; ++i ) {
      count=0;
      currpos=i;
      //cout << "currpos: " << currpos << "(" << i << ")";
      while ( currpos<n && (fabs(atom[currpos+1].x[0]-atom[currpos].x[0])<=SINGLECOORDEPS) ) { //see globaldefs.h
         ++count;
         ++currpos;
      }
      //cout << ", count: " << count << endl;
      if ( count!=0 ) {
         QuickSort(i+count,i,1);
         i+=count;
      }
   }
   //Sorts third column
   for ( size_t i=0 ; i<n ; ++i ) {
      count=0;
      currpos=i;
      //cout << "currpos: " << currpos << "(" << i << ")";
      while ( currpos<n &&\
              (fabs(atom[currpos+1].x[1]-atom[currpos].x[1])<=SINGLECOORDEPS) &&\
              (fabs(atom[currpos+1].x[0]-atom[currpos].x[0])<=SINGLECOORDEPS) ) { //see globaldefs.h
         ++count;
         ++currpos;
      }
      //cout << ", count: " << count << endl;
      if ( count!=0 ) {
         QuickSort(i+count,i,2);
         i+=count;
      }
   }
}
void Molecule::QuickSort(int high,int low,int srtIdx) {
   int i=low;
   int j=high;
   int hplo2=(high+low)/2;
   double mid[3];
   mid[0]=atom[hplo2].x[0];
   mid[1]=atom[hplo2].x[1];
   mid[2]=atom[hplo2].x[2];
   do {
      while (atom[i].x[srtIdx]<mid[srtIdx]) { ++i; }
      while (atom[j].x[srtIdx]>mid[srtIdx]) { --j; }
      if ( i<=j ) {
         std::swap(atom[i],atom[j]);
         ++i;
         --j;
      }
   } while ( i<j );
   if ( low<j ) { QuickSort(j,low,srtIdx); }
   if ( i<high ) { QuickSort(high,i,srtIdx); }
}
string Molecule::EmpiricalFormula() const {
   size_t n=MAXATNUMDEF; //see globaldefs.h
   vector<int> count(n);
   for ( size_t i=0 ; i<n ; ++i ) { count[i]=0; }
   for ( size_t i=0 ; i<atom.size() ; ++i ) { count[atom[i].num-1]++; }
   string res;
   for ( size_t i=0 ; i<n ; ++i ) {
      if ( count[i]>0 ) {
         res+=Atom::GetAtomicSymbol(i+1);
         res+=std::to_string(count[i]);
      }
   }
   return res;
}
int Molecule::CountAtomsOfType(string ss) {
   int atnum=Atom::GetAtomicNumberFromSymbol(ss);
   return CountAtomsOfType(atnum);
}
int Molecule::CountAtomsOfType(int nn) {
   int count=0;
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      if ( atom[i].num==nn ) { ++count; }
   }
   return count;
}
void Molecule::DetermineBoundingBox() {
   size_t nNuc=atom.size();
   double bmin[3]={1.0e+50,1.0e+50,1.0e+50};
   double bmax[3]={-1.0e+50,-1.0e+50,-1.0e+50};
   double xx[3];
   for ( size_t i=0 ; i<nNuc ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) { xx[j]=atom[i].x[j]; }
      if ( xx[0]<bmin[0] ) { bmin[0]=xx[0]; }
      if ( xx[1]<bmin[1] ) { bmin[1]=xx[1]; }
      if ( xx[2]<bmin[2] ) { bmin[2]=xx[2]; }
      if ( xx[0]>bmax[0] ) { bmax[0]=xx[0]; }
      if ( xx[1]>bmax[1] ) { bmax[1]=xx[1]; }
      if ( xx[2]>bmax[2] ) { bmax[2]=xx[2]; }
   }
   for ( size_t i=0 ; i<3 ; ++i ) { xmin[i]=bmin[i]; }
   for ( size_t i=0 ; i<3 ; ++i ) { xmax[i]=bmax[i]; }
   rmax=-1.0e+50;
   for ( size_t i=0 ; i<3 ; ++i ) { if ( fabs(bmin[i]) > rmax ) { rmax = fabs(bmin[i]); } }
   for ( size_t i=0 ; i<3 ; ++i ) { if ( fabs(bmax[i]) > rmax ) { rmax = fabs(bmax[i]); } }
}
void Molecule::ComputeCenterOfMass() {
   if ( !ImSetup() ) { cout << __FILE__ << ", line: " << __LINE__ << '\n'; return; }
   for ( size_t i=0 ; i<3 ; ++i ) { cm[i]=0.0e0; }
   double m,totmass=0.0e0;
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      m=atom[i].weight;
      totmass+=m;
      for ( size_t j=0 ; j<3 ; ++j ) { cm[j]+= (m*atom[i].x[j]); }
   }
   for ( size_t i=0 ; i<3 ; ++i ) { cm[i]/=totmass; }
}
void Molecule::ComputeCentroid() {
   if ( !ImSetup() ) { cout << __FILE__ << ", line: " << __LINE__ << '\n'; return; }
   for ( size_t i=0 ; i<3 ; ++i ) { cd[i]=0.0e0; }
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      for ( size_t j=0 ; j<3 ; ++j ) { cd[j]+=atom[i].x[j]; }
   }
   double f=1.0e0/double(atom.size());
   for ( size_t i=0 ; i<3 ; ++i ) { cd[i]*=f; }
}
bool Molecule::ImSetup() const {
   if ( !imsetup ) {
      ScreenUtils::DisplayErrorMessage("The molecule is not setup!");
   }
   return imsetup;
}
void Molecule::CenterAtCentroid() {
   ComputeCentroid();
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      atom[i].x[0]-=cd[0];
      atom[i].x[1]-=cd[1];
      atom[i].x[2]-=cd[2];
   }
   for ( size_t i=0 ; i<3 ; ++i ) { origCent[i]-=cd[i]; cd[i]=0.0e0; }
}
void Molecule::CenterAtCenterOfMass() {
   ComputeCenterOfMass();
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      atom[i].x[0]-=cm[0];
      atom[i].x[1]-=cm[1];
      atom[i].x[2]-=cm[2];
   }
   for ( size_t i=0 ; i<3 ; ++i ) { origCent[i]-=cm[i]; cm[i]=0.0e0; }
}
void Molecule::ResetOriginOfCoordinates() {
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      atom[i].x[0]-=origCent[0];
      atom[i].x[1]-=origCent[1];
      atom[i].x[2]-=origCent[2];
   }
   for ( size_t i=0 ; i<3 ; ++i ) { origCent[i]=0.0e0; }
}
void Molecule::SetupBonds() {
   SetupCells();
   int nNuc=int(atom.size());
   if ( int(bond.size()) != nNuc ) {
      for ( int i=0 ; i<int(bond.size()) ; ++i ) { bond[i].clear(); }
      bond.resize(nNuc);
      for ( int i=0 ; i<nNuc ; ++i ) {
         bond[i].reserve(8);
      }
   }
   if ( int(bndDist.size()) != nNuc ) {
      for ( int i=0 ; i<int(bndDist.size()) ; ++i ) { bndDist[i].clear(); }
      bndDist.resize(nNuc);
      for ( int i=0 ; i<nNuc ; ++i ) {
         bndDist[i].reserve(8);
      }
   }
   maxBondDist=-1.0e+50;
   double d,vdwd,clsstatd=1.0e+50;
   int al,am;
   for ( int i=1 ; i<int(cell.size()-1) ; ++i ) {
      for ( int j=1 ; j<int(cell[i].size()-1) ; ++j ) {
         for ( int k=1 ; k<int(cell[i][j].size()-1) ; ++k ) {
            for ( int l=0 ; l<int(cell[i][j][k].size()) ; ++l ) {
               al=cell[i][j][k][l];
               for ( int p=-1 ; p<2 ; ++p ) {
                  for ( int q=-1 ; q<2 ; ++q ) {
                     for ( int r=-1 ; r<2 ; ++r ) {
                        for ( int s=0 ; s<int(cell[i+p][j+q][k+r].size()) ; ++s ) {
                           am=cell[i+p][j+q][k+r][s];
                           if ( al>=am ) { continue; }
                           d=MatrixVectorOperations3D::Distance(atom[al].x,atom[am].x);
                           vdwd=(atom[al].GetVDWRadius()+atom[am].GetVDWRadius());
                           if ( d<=vdwd ) {
                              bond[al].push_back(am);
                              bndDist[al].push_back(d);
                              if ( d>maxBondDist ) { maxBondDist=d; }
                           }
                           if ( d<clsstatd ) {clsstatd=d;}
                        }
                     }
                  }
               }
            }
         }
      }
   }
   /*
   for ( int i=0 ; i<nNuc ; ++i ) {
      for ( int j=i+1 ; j<nNuc ; ++j ) {
         d=MatrixVectorOperations3D::Distance(atom[i].x,atom[j].x);
         vdwd=(atom[i].GetVDWRadius()+atom[j].GetVDWRadius());
         if ( d<=vdwd ) {
            bond[i].push_back(j);
            bndDist[i].push_back(d);
            if ( d>maxBondDist ) { maxBondDist=d; }
         }
         if ( d<clsstatd ) {clsstatd=d;}
      }
   }
   // */
   if (nNuc==1) {maxBondDist=3.0e0;}
   if (maxBondDist<0.0e0) {
      if ( nNuc==2 ) {
         maxBondDist=MatrixVectorOperations3D::Distance(atom[0].x,atom[1].x);
      } else {
         maxBondDist=clsstatd;
      }
   }
}
void Molecule::SetupCells() {
   DetermineBoundingBox();
   double len[3],lmin[3],lmax[3];
   for ( size_t i=0 ; i<3 ; ++i ) { lmin[i]=xmin[i]-cellLen; }
   for ( size_t i=0 ; i<3 ; ++i ) { lmax[i]=xmax[i]+cellLen; }
   cout << "xmn: " << xmin[0] << ' ' << xmin[1] << ' ' << xmin[2] << '\n';
   cout << "lmn: " << lmin[0] << ' ' << lmin[1] << ' ' << lmin[2] << '\n';
   cout << "xmx: " << xmax[0] << ' ' << xmax[1] << ' ' << xmax[2] << '\n';
   cout << "lmx: " << lmax[0] << ' ' << lmax[1] << ' ' << lmax[2] << '\n';
   for ( size_t i=0 ; i<3 ; ++i ) { len[i]=xmax[i]-xmin[i]; }
   //cout << "actlen: " << len[0] << ' ' << len[1] << ' ' << len[2] << '\n';
   //for ( size_t i=0 ; i<3 ; ++i ) { len[i]=lmax[i]-lmin[i]; }
   cout << "len: " << len[0] << ' ' << len[1] << ' ' << len[2] << '\n';
   int nax=ceil(len[0]/cellLen), nx=nax+3;
   int nay=ceil(len[1]/cellLen), ny=nay+3;
   int naz=ceil(len[2]/cellLen), nz=naz+3;
   cout << "nc: " << nx << ' ' << ny << ' ' << nz << '\n';
   cell.resize(nx);
   for ( int i=0 ; i<nx ; ++i ) {
      cell[i].resize(ny);
      for ( int j=0 ; j<ny ; ++j ) {
         cell[i][j].resize(nz);
         for ( int k=0 ; k<nz ; ++k ) {
            cell[i][j][k].reserve(8);
         }
      }
   }
   int px,py,pz;
   //int mxx=-10,mxy=-10,mxz=-10;
   //int mnx=1000,mny=1000,mnz=1000;
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      px=floor(double(nax)*(atom[i].x[0]-xmin[0])/len[0])+1;
      py=floor(double(nay)*(atom[i].x[1]-xmin[1])/len[1])+1;
      pz=floor(double(naz)*(atom[i].x[2]-xmin[2])/len[2])+1;
      if ( px>=(nx-1) || 1>nx ) { cout << "Error!" << px << '\n'; continue; }
      if ( py>=(ny-1) || 1>ny ) { cout << "Error!" << py << '\n'; continue; }
      if ( pz>=(nz-1) || 1>nz ) { cout << "Error!" << pz << '\n'; continue; }
      cell[px][py][pz].push_back(int(i));
      //if ( px > mxx ) { mxx=px; }
      //if ( py > mxy ) { mxy=py; }
      //if ( pz > mxz ) { mxz=pz; }
      //if ( px < mnx ) { mnx=px; }
      //if ( py < mny ) { mny=py; }
      //if ( pz < mnz ) { mnz=pz; }
   }
   //cout << "minpos: " << mnx << ' ' << mny << ' ' << mnz << '\n';
   //cout << "maxpos: " << mxx << ' ' << mxy << ' ' << mxz << '\n';
}
vector<size_t> Molecule::ListOfNeighbours(const size_t atpos) {
   vector<size_t> res(0);
   res.reserve(4);
   if ( atpos==string::npos || atpos>=atom.size() ) {
      ScreenUtils::DisplayErrorMessage("Non existent atom!");
      cout << __FILE__ << ", fnc: " << __FUNCTION__ << ", line: " << __LINE__ << '\n';
      res.clear();
      return res;
   }
   res.push_back(atpos);
   double dx2;
   double refVdWr=atom[atpos].GetVDWRadius();
   double x0[3],x1[3];
   x0[0]=atom[atpos].x[0]; x0[1]=atom[atpos].x[1]; x0[2]=atom[atpos].x[2];
   double dVdWr;
   for ( size_t i=0 ; i<atom.size() ; ++i ) {
      if ( atpos==i ) { continue; }
      dVdWr=refVdWr+atom[i].GetVDWRadius();
      dVdWr*=dVdWr;
      x1[0]=atom[i].x[0]; x1[1]=atom[i].x[1]; x1[2]=atom[i].x[2];
      dx2=(x1[0]-x0[0])*(x1[0]-x0[0])+(x1[1]-x0[1])*(x1[1]-x0[1])+(x1[2]-x0[2])*(x1[2]-x0[2]);
      if ( dx2<=dVdWr ) { res.push_back(i); }
   }
   return res;
}
std::ostream &operator<<(std::ostream &out,const Molecule (&mol)) {
   if ( !mol.ImSetup() ) { return out; }
   out << "The molecule has " << mol.Size() << " atoms." << endl;
   if ( mol.Size()<1 ) { return out; }
   for ( size_t i=0 ; i<(mol.Size()-1) ; ++i ) {
      out << mol.atom[i] << endl;
   }
   out << mol.atom[mol.Size()-1] << endl;
   out << "Empirical formula: " << mol.EmpiricalFormula();
   return out;
}
std::ostream &operator<<(std::ostream &out,const Molecule* mol) {
   out << *mol;
   return out;
}
bool operator== (const Molecule &m1, const Molecule &m2) {
   size_t n1=m1.Size();
   if ( n1!=m2.Size() ) { return false; }
   bool tmpbool;
   vector<double> tmpvec;
   double d;
   int nat;
   for ( size_t i=0 ; i<n1 ; ++i ) {
      tmpvec=m1.atom[i].x;
      tmpbool=false;
      nat=m1.atom[i].num;
      for ( size_t j=0 ; j<n1 ; ++j ) {
         d=MatrixVectorOperations3D::Distance(tmpvec,m2.atom[j].x);
         if ( d<COORDSEPSILON && (nat==m2.atom[j].num) ) {
            tmpbool=true;
            break;
         }
      }
      if ( !tmpbool ) { return false; }
   }
   return true;
}
void DuplicateAtomOrder(Molecule &m1,Molecule &m2) {
   size_t n1=m1.Size();
   vector<double> tmpvec;
   double d;
   int nat;
   for ( size_t i=0 ; i<n1 ; ++i ) {
      tmpvec=m1.atom[i].x;
      nat=m1.atom[i].num;
      for ( size_t j=0 ; j<n1 ; ++j ) {
         d=MatrixVectorOperations3D::Distance(tmpvec,m2.atom[j].x);
         if ( d<COORDSEPSILON && (nat==m2.atom[j].num) ) {
            std::swap(m2.atom[i],m2.atom[j]);
            break;
         }
      }
   }
}

