/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.0.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2024, Juan Manuel Solano-Altamirano
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
#ifndef _GAUSSIANCUBE_H_
#define _GAUSSIANCUBE_H_
#include <string>
using std::string;
#include "inputmolecule_cub.h"

/* ************************************************************************** */
class GaussianCube {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   GaussianCube();
   GaussianCube(const string &cubname);
   ~GaussianCube();
   Molecule& GetMolecule(void) {return mol;}
   size_t Nx() const {return nx;}
   size_t Ny() const {return ny;}
   size_t Nz() const {return nz;}
   double X0(int idx) const {return x0[idx];}
   double DX(int idx1,int idx2) const {return dx[idx1][idx2];};
   double DX(int idx) const {return dx[idx][idx];};
   vector<double> DX() const {vector<double> r{dx[0][0],dx[1][1],dx[2][2]}; return r;}
   vector<double> X0() const {vector<double> r{x0[0],x0[1],x0[2]}; return r;}
   double Data(size_t xIdx,size_t yIdx,size_t zIdx) const {return data[xIdx*nyz+yIdx*nz+zIdx];}
   inline vector<double> GetX(int xIdx,int yIdx,int zIdx) const {
      vector<double> s=x0;
      for ( size_t i=0 ; i<3 ; ++i ) {
         s[i]+=(double(xIdx)*dx[i][0]+double(yIdx)*dx[i][1]+double(zIdx)*dx[i][2]);
      }
      return s;}
   bool CubeLoaded() { return cubeLoaded; }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   InputMoleculeCub mol;
   vector<double> x0;
   vector<vector<double> > dx;
   size_t nx,ny,nz;
   size_t nyz;
   vector<double> data;
   bool useBohrUnits;
   bool cubeLoaded;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _GAUSSIANCUBE_H_ */

