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

