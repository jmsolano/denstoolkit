#ifndef _PALETTE_H_
#define _PALETTE_H_
#include <vector>
using std::vector;

/* ************************************************************************** */
class Palette {
/* ************************************************************************** */
public:
/* ************************************************************************** */
   void SetBlues() { pal=&blues; }
   void SetBugn() { pal=&bugn; }
   void SetGnbu() { pal=&gnbu; }
   void SetGreens() { pal=&greens; }
   void SetGreys() { pal=&greys; }
   void SetOranges() { pal=&oranges; }
   void SetOrrd() { pal=&orrd; }
   void SetPubu() { pal=&pubu; }
   void SetPurples() { pal=&purples; }
   void SetRdbu() { pal=&rdbu; }
   void SetRdylbu() { pal=&rdylbu; }
   void SetRdylgn() { pal=&rdylgn; }
   void SetReds() { pal=&reds; }
   void SetSpectral() { pal=&spectral; }
   void SetYlgn() { pal=&ylgn; }
   void SetYlgnbu() { pal=&ylgnbu; }
   void SetYlorbr() { pal=&ylorbr; }
   void SetYlorrd() { pal=&ylorrd; }
/* ************************************************************************** */
   void GetRGB(const size_t pos,double &r,double  &g,double  &b) const;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   /* cat the contents of ../expandpalette/VECPALDEFS (see ../expandpalette/README  */
   static vector<vector<double> > blues;
   static vector<vector<double> > bugn;
   static vector<vector<double> > gnbu;
   static vector<vector<double> > greens;
   static vector<vector<double> > greys;
   static vector<vector<double> > oranges;
   static vector<vector<double> > orrd;
   static vector<vector<double> > pubu;
   static vector<vector<double> > purples;
   static vector<vector<double> > rdbu;
   static vector<vector<double> > rdylbu;
   static vector<vector<double> > rdylgn;
   static vector<vector<double> > reds;
   static vector<vector<double> > spectral;
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

