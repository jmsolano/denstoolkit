/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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
#ifndef _PALETTE_H_
#define _PALETTE_H_
#include <vector>
using std::vector;
#include <string>
using std::string;

/* ************************************************************************** */
class Palette {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   void GetRGB(const size_t pos,double &r,double  &g,double  &b) const;
   void GetRGB(const size_t pos,int &r,int  &g,int  &b) const;
   inline void GetHSL(const size_t pos,double &h,double  &s,double  &l) const
      { h=((*pal)[pos][0]); s=((*pal)[pos][1]); l=((*pal)[pos][2]); }
/* ************************************************************************** */
   Palette();
   void SelectPalette(const string &unm);
   /* cat the contents of ../expandpalette/SETPALFNCTS (see ../expandpalette/README  */
   void SetBentcoolwarm() { pal=&bentcoolwarm; }
   void SetBlues() { pal=&blues; }
   void SetBugn() { pal=&bugn; }
   void SetGnbu() { pal=&gnbu; }
   void SetGreens() { pal=&greens; }
   void SetGreys() { pal=&greys; }
   void SetInferno() { pal=&inferno; }
   void SetMagma() { pal=&magma; }
   void SetMoreland() { pal=&moreland; }
   void SetOranges() { pal=&oranges; }
   void SetOrrd() { pal=&orrd; }
   void SetPlasma() { pal=&plasma; }
   void SetPubu() { pal=&pubu; }
   void SetPurples() { pal=&purples; }
   void SetRdbu() { pal=&rdbu; }
   void SetRdylbu() { pal=&rdylbu; }
   void SetRdylgn() { pal=&rdylgn; }
   void SetReds() { pal=&reds; }
   void SetSpectral() { pal=&spectral; }
   void SetViridis() { pal=&viridis; }
   void SetYlgn() { pal=&ylgn; }
   void SetYlgnbu() { pal=&ylgnbu; }
   void SetYlorbr() { pal=&ylorbr; }
   void SetYlorrd() { pal=&ylorrd; }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   /* cat the contents of ../expandpalette/VECPALDEFS (see ../expandpalette/README  */
   static vector<vector<double> > bentcoolwarm;
   static vector<vector<double> > blues;
   static vector<vector<double> > bugn;
   static vector<vector<double> > gnbu;
   static vector<vector<double> > greens;
   static vector<vector<double> > greys;
   static vector<vector<double> > inferno;
   static vector<vector<double> > magma;
   static vector<vector<double> > moreland;
   static vector<vector<double> > oranges;
   static vector<vector<double> > orrd;
   static vector<vector<double> > plasma;
   static vector<vector<double> > pubu;
   static vector<vector<double> > purples;
   static vector<vector<double> > rdbu;
   static vector<vector<double> > rdylbu;
   static vector<vector<double> > rdylgn;
   static vector<vector<double> > reds;
   static vector<vector<double> > spectral;
   static vector<vector<double> > viridis;
   static vector<vector<double> > ylgn;
   static vector<vector<double> > ylgnbu;
   static vector<vector<double> > ylorbr;
   static vector<vector<double> > ylorrd;
/* ************************************************************************** */
   vector<vector<double> > *pal;
   static constexpr double oo255=1.0e0/255.0e0;
};
/* ************************************************************************** */


#endif  /* _PALETTE_H_ */

