/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.1
  
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
      { h=(pal[pos*3+0]); s=(pal[pos*3+1]); l=(pal[pos*3+2]); }
/* ************************************************************************** */
   Palette();
   void SelectPalette(const string &unm);
   /* cat the contents of ../expandpalette/SETPALFNCTS (see ../expandpalette/README  */
   void SetBentcoolwarm() { pal=bentcoolwarm; }
   void SetBlues() { pal=blues; }
   void SetBugn() { pal=bugn; }
   void SetGnbu() { pal=gnbu; }
   void SetGreens() { pal=greens; }
   void SetGreys() { pal=greys; }
   void SetInferno() { pal=inferno; }
   void SetMagma() { pal=magma; }
   void SetMoreland() { pal=moreland; }
   void SetOranges() { pal=oranges; }
   void SetOrrd() { pal=orrd; }
   void SetPlasma() { pal=plasma; }
   void SetPubu() { pal=pubu; }
   void SetPurples() { pal=purples; }
   void SetRdbu() { pal=rdbu; }
   void SetRdylbu() { pal=rdylbu; }
   void SetRdylgn() { pal=rdylgn; }
   void SetReds() { pal=reds; }
   void SetSpectral() { pal=spectral; }
   void SetViridis() { pal=viridis; }
   void SetYlgn() { pal=ylgn; }
   void SetYlgnbu() { pal=ylgnbu; }
   void SetYlorbr() { pal=ylorbr; }
   void SetYlorrd() { pal=ylorrd; }
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   /* cat the contents of ../expandpalette/VECPALDEFS (see ../expandpalette/README  */
   const static double bentcoolwarm[768];
   const static double blues[768];
   const static double bugn[768];
   const static double gnbu[768];
   const static double greens[768];
   const static double greys[768];
   const static double inferno[768];
   const static double magma[768];
   const static double moreland[768];
   const static double oranges[768];
   const static double orrd[768];
   const static double plasma[768];
   const static double pubu[768];
   const static double purples[768];
   const static double rdbu[768];
   const static double rdylbu[768];
   const static double rdylgn[768];
   const static double reds[768];
   const static double spectral[768];
   const static double viridis[768];
   const static double ylgn[768];
   const static double ylgnbu[768];
   const static double ylorbr[768];
   const static double ylorrd[768];
/* ************************************************************************** */
   const double *pal;
   static constexpr double oo255=1.0e0/255.0e0;
};
/* ************************************************************************** */


#endif  /* _PALETTE_H_ */

