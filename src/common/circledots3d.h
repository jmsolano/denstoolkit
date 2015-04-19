/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.0.1
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

#ifndef _CIRCLEDOTS3D_H_
#define _CIRCLEDOTS3D_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
#endif
#include <string>
using std::string;

/* ************************************************************************** */
class CircleDots3D {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   CircleDots3D();
   ~CircleDots3D();
/* ************************************************************************** */
   solreal getCartCoord(const int i,const int j);
/* ************************************************************************** */
   solreal getPhi(const int i);
/* ************************************************************************** */
   void getCartCoords(const int i,solreal (&xx)[3]);
/* ************************************************************************** */
   void setE1(const solreal x,const solreal y,const solreal z);
/* ************************************************************************** */
   void setE2(const solreal x,const solreal y,const solreal z);
/* ************************************************************************** */
   void setE1AndE2(const solreal (&ee1)[3],const solreal (&ee2)[3]);
/* ************************************************************************** */
   void setOrigin(const solreal x,const solreal y,const solreal z);
/* ************************************************************************** */
   void computeUE3(void);
/* ************************************************************************** */
   void setNPts(const int nn) {npts_=nn;}
/* ************************************************************************** */
   void setRadius(const solreal rr) {radius_=rr;}
/* ************************************************************************** */
   void setupCircle(void);
/* ************************************************************************** */
   void displayCoordinates(void);
/* ************************************************************************** */
   void writeCoordinates(const string &oname,bool wrtoo=false);
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void init();
/* ************************************************************************** */
   int npts_; /*!< Number of points of the circle discretization  */
/* ************************************************************************** */
   solreal radius_; /*!< Radius of the circle  */
/* ************************************************************************** */
   solreal e1_[3]; /*!< First unit vector. The circle will be drawn in 
                        the plane spanned by \f$\{\vec{e}_1\vec{e}_2\}\f$.  */
   solreal e2_[3]; /*!< Second unit vector. The circle will be drawn in 
                        the plane spanned by \f$\{\vec{e}_1\vec{e}_2\}\f$.  */
   solreal e3_[3]; /*!< This will contain the vector 
                        \f$\vec{e}_3\equiv\vec{e_1}\times\vec{e}_2/
                        |\vec{e_1}\times\vec{e}_2|\f$  */
   solreal ue1_[3]; /*!< First user unit vector. The circle will be drawn in 
                        the plane spanned by \f$\{\vec{e}_1\vec{e}_2\}\f$.  */
   solreal ue2_[3]; /*!< Second user unit vector. The circle will be drawn in 
                        the plane spanned by \f$\{\vec{e}_1\vec{e}_2\}\f$.  */
   solreal ue3_[3]; /*!< This will contain the user vector 
                        \f$\vec{e}_3\equiv\vec{e_1}\times\vec{e}_2/
                        |\vec{e_1}\times\vec{e}_2|\f$  */
   solreal oo_[3]; /*!< This will contain the origin vector (center of
                        coordinates).  */
/* ************************************************************************** */
   /** This array contains the Cartesian coordinates of the points within
    * the circle's boundary. In \f$xx\_[i][j]\f$, the jth-Cartesian coordinate
    * of the i-th point is stored. When j==3 the stored value is the 
    * angle measured counter clock-wise w.r.t. the vector \f$\hat e_1\f$. */
   solreal **xx_;
/* ************************************************************************** */
   bool havee1,havee2,havee3;
   bool imsetup;
/* ************************************************************************** */
   solreal dphi_;
/* ************************************************************************** */
   static const solreal twoPi;
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _CIRCLEDOTS3D_H_ */

