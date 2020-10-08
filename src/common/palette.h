#ifndef _PALETTE_H_
#define _PALETTE_H_
#include <vector>
using std::vector;

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

