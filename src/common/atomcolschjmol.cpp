/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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
#ifndef _COL_SCHEME_JMOL_H_
#define _COL_SCHEME_JMOL_H_

#define _HAVE_SELECTED_ATOM_PALETTE_ 1

#define MAXDEFINEDATOMICCOLORS 109

int getAtomicRGBColorInt(int nat,int rgb) {
   static const int atomicColorInt[MAXDEFINEDATOMICCOLORS][3]={ //109
      {255, 255, 255},  //  H  1 --
      {217, 255, 255},  //  He 2 --
      {204, 128, 255},  //  Li 3 --
      {194, 255, 0},  //  Be 4 --
      {255, 181, 181},  //  B  5 --
      {144, 144, 144},  //  C  6 - changed from ghemical --
      {48, 80, 248},  //  N  7 - changed from ghemical --
      {255, 13, 13},  //  O  8 --
      {144, 224, 80},  //  F  9 - changed from ghemical --
      {179, 227, 245},  //  Ne 10 --
      {171, 92, 242},  //  Na 11 --
      {138, 255, 0},  //  Mg 12 --
      {191, 166, 166},  //  Al 13 --
      {240, 200, 160},  //  Si 14 - changed from ghemical --
      {255, 128, 0},  //  P  15 --
      {255, 255, 48},  //  S  16 --
      {31, 240, 31},  //  Cl 17 --
      {128, 209, 227},  //  Ar 18 --
      {143, 64, 212},  //  K  19 --
      {61, 255, 0},  //  Ca 20 --
      {230, 230, 230},  //  Sc 21 --
      {191, 194, 199},  //  Ti 22 --
      {166, 166, 171},  //  V  23 --
      {138, 153, 199},  //  Cr 24 --
      {156, 122, 199},  //  Mn 25 --
      {224, 102, 51},  //  Fe 26 - changed from ghemical --
      {240, 144, 160},  //  Co 27 - changed from ghemical --
      {80, 208, 80},  //  Ni 28 - changed from ghemical --
      {200, 128, 51},  //  Cu 29 - changed from ghemical --
      {125, 128, 176},  //  Zn 30 --
      {194, 143, 143},  //  Ga 31 --
      {102, 143, 143},  //  Ge 32 --
      {189, 128, 227},  //  As 33 --
      {255, 161, 0},  //  Se 34 --
      {166, 41, 41},  //  Br 35 --
      {92, 184, 209},  //  Kr 36 --
      {112, 46, 176},  //  Rb 37 --
      {0, 255, 0},  //  Sr 38 --
      {148, 255, 255},  //  Y  39 --
      {148, 224, 224},  //  Zr 40 --
      {115, 194, 201},  //  Nb 41 --
      {84, 181, 181},  //  Mo 42 --
      {59, 158, 158},  //  Tc 43 --
      {36, 143, 143},  //  Ru 44 --
      {10, 125, 140},  //  Rh 45 --
      {0, 105, 133},  //  Pd 46 --
      {192, 192, 192},  //  Ag 47 - changed from ghemical --
      {255, 217, 143},  //  Cd 48 --
      {166, 117, 115},  //  In 49 --
      {102, 128, 128},  //  Sn 50 --
      {158, 99, 181},  //  Sb 51 --
      {212, 122, 0},  //  Te 52 --
      {148, 0, 148},  //  I  53 --
      {66, 158, 176},  //  Xe 54 --
      {87, 23, 143},  //  Cs 55 --
      {0, 201, 0},  //  Ba 56 --
      {112, 212, 255},  //  La 57 --
      {255, 255, 199},  //  Ce 58 --
      {217, 255, 199},  //  Pr 59 --
      {199, 255, 199},  //  Nd 60 --
      {163, 255, 199},  //  Pm 61 --
      {143, 255, 199},  //  Sm 62 --
      {97, 255, 199},  //  Eu 63 --
      {69, 255, 199},  //  Gd 64 --
      {48, 255, 199},  //  Tb 65 --
      {31, 255, 199},  //  Dy 66 --
      {0, 255, 156},  //  Ho 67 --
      {0, 230, 117},  //  Er 68 --
      {0, 212, 82},  //  Tm 69 --
      {0, 191, 56},  //  Yb 70 --
      {0, 171, 36},  //  Lu 71 --
      {77, 194, 255},  //  Hf 72 --
      {77, 166, 255},  //  Ta 73 --
      {33, 148, 214},  //  W  74 --
      {38, 125, 171},  //  Re 75 --
      {38, 102, 150},  //  Os 76 --
      {23, 84, 135},  //  Ir 77 --
      {208, 208, 224},  //  Pt 78 - changed from ghemical --
      {255, 209, 35},  //  Au 79 - changed from ghemical --
      {184, 184, 208},  //  Hg 80 - changed from ghemical --
      {166, 84, 77},  //  Tl 81 --
      {87, 89, 97},  //  Pb 82 --
      {158, 79, 181},  //  Bi 83 --
      {171, 92, 0},  //  Po 84 --
      {117, 79, 69},  //  At 85 --
      {66, 130, 150},  //  Rn 86 --
      {66, 0, 102},  //  Fr 87 --
      {0, 125, 0},  //  Ra 88 --
      {112, 171, 250},  //  Ac 89 --
      {0, 186, 255},  //  Th 90 --
      {0, 161, 255},  //  Pa 91 --
      {0, 143, 255},  //  U  92 --
      {0, 128, 255},  //  Np 93 --
      {0, 107, 255},  //  Pu 94 --
      {84, 92, 242},  //  Am 95 --
      {120, 92, 227},  //  Cm 96 --
      {138, 79, 227},  //  Bk 97 --
      {161, 54, 212},  //  Cf 98 --
      {179, 31, 212},  //  Es 99 --
      {179, 31, 186},  //  Fm 100 --
      {179, 13, 166},  //  Md 101 --
      {189, 13, 135},  //  No 102 --
      {199, 0, 102},  //  Lr 103 --
      {204, 0, 89},  //  Rf 104 --
      {209, 0, 79},  //  Db 105 --
      {217, 0, 69},  //  Sg 106 --
      {224, 0, 56},  //  Bh 107 --
      {230, 0, 46},  //  Hs 108 --
      {235, 0, 38} //  Mt 109 --
   };
   return atomicColorInt[nat][rgb];
}
int getAtomicRColorInt(int nat) {
   return getAtomicRGBColorInt(nat,0);
}
int getAtomicGColorInt(int nat) {
   return getAtomicRGBColorInt(nat,1);
}
int getAtomicBColorInt(int nat) {
   return getAtomicRGBColorInt(nat,2);
}
solreal getAtomicRGBColorReal(int nat,int rgb) {
   static const solreal atomicColor[MAXDEFINEDATOMICCOLORS][3]={
      {1.000000,1.000000,1.000000}, //  H  1 --
      {0.8509804,1.000000,1.000000}, //  He 2 --
      {0.8000000,0.5019608,1.000000}, //  Li 3 --
      {0.7607843,1.000000,0.000000}, //  Be 4 --
      {1.000000,0.7098039,0.7098039}, //  B  5 --
      {0.5647059,0.5647059,0.5647059}, //  C  6 - changed from ghemical --
      {0.1882353,0.3137255,0.9725490}, //  N  7 - changed from ghemical --
      {1.000000,0.05098039,0.05098039}, //  O  8 --
      {0.5647059,0.8784314,0.3137255}, //  F  9 - changed from ghemical --
      {0.7019608,0.8901961,0.9607843}, //  Ne 10 --
      {0.6705882,0.3607843,0.9490196}, //  Na 11 --
      {0.5411765,1.000000,0.000000}, //  Mg 12 --
      {0.7490196,0.6509804,0.6509804}, //  Al 13 --
      {0.9411765,0.7843137,0.6274510}, //  Si 14 - changed from ghemical --
      {1.000000,0.5019608,0.000000}, //  P  15 --
      {1.000000,1.000000,0.1882353}, //  S  16 --
      {0.1215686,0.9411765,0.1215686}, //  Cl 17 --
      {0.5019608,0.8196078,0.8901961}, //  Ar 18 --
      {0.5607843,0.2509804,0.8313725}, //  K  19 --
      {0.2392157,1.000000,0.000000}, //  Ca 20 --
      {0.9019608,0.9019608,0.9019608}, //  Sc 21 --
      {0.7490196,0.7607843,0.7803922}, //  Ti 22 --
      {0.6509804,0.6509804,0.6705882}, //  V  23 --
      {0.5411765,0.6000000,0.7803922}, //  Cr 24 --
      {0.6117647,0.4784314,0.7803922}, //  Mn 25 --
      {0.8784314,0.4000000,0.2000000}, //  Fe 26 - changed from ghemical --
      {0.9411765,0.5647059,0.6274510}, //  Co 27 - changed from ghemical --
      {0.3137255,0.8156863,0.3137255}, //  Ni 28 - changed from ghemical --
      {0.7843137,0.5019608,0.2000000}, //  Cu 29 - changed from ghemical --
      {0.4901961,0.5019608,0.6901961}, //  Zn 30 --
      {0.7607843,0.5607843,0.5607843}, //  Ga 31 --
      {0.4000000,0.5607843,0.5607843}, //  Ge 32 --
      {0.7411765,0.5019608,0.8901961}, //  As 33 --
      {1.000000,0.6313725,0.000000}, //  Se 34 --
      {0.6509804,0.1607843,0.1607843}, //  Br 35 --
      {0.3607843,0.7215686,0.8196078}, //  Kr 36 --
      {0.4392157,0.1803922,0.6901961}, //  Rb 37 --
      {0.000000,1.000000,0.000000}, //  Sr 38 --
      {0.5803922,1.000000,1.000000}, //  Y  39 --
      {0.5803922,0.8784314,0.8784314}, //  Zr 40 --
      {0.4509804,0.7607843,0.7882353}, //  Nb 41 --
      {0.3294118,0.7098039,0.7098039}, //  Mo 42 --
      {0.2313725,0.6196078,0.6196078}, //  Tc 43 --
      {0.1411765,0.5607843,0.5607843}, //  Ru 44 --
      {0.03921569,0.4901961,0.5490196}, //  Rh 45 --
      {0.000000,0.4117647,0.5215686}, //  Pd 46 --
      {0.7529412,0.7529412,0.7529412}, //  Ag 47 - changed from ghemical --
      {1.000000,0.8509804,0.5607843}, //  Cd 48 --
      {0.6509804,0.4588235,0.4509804}, //  In 49 --
      {0.4000000,0.5019608,0.5019608}, //  Sn 50 --
      {0.6196078,0.3882353,0.7098039}, //  Sb 51 --
      {0.8313725,0.4784314,0.000000}, //  Te 52 --
      {0.5803922,0.000000,0.5803922}, //  I  53 --
      {0.2588235,0.6196078,0.6901961}, //  Xe 54 --
      {0.3411765,0.09019608,0.5607843}, //  Cs 55 --
      {0.000000,0.7882353,0.000000}, //  Ba 56 --
      {0.4392157,0.8313725,1.000000}, //  La 57 --
      {1.000000,1.000000,0.7803922}, //  Ce 58 --
      {0.8509804,1.000000,0.7803922}, //  Pr 59 --
      {0.7803922,1.000000,0.7803922}, //  Nd 60 --
      {0.6392157,1.000000,0.7803922}, //  Pm 61 --
      {0.5607843,1.000000,0.7803922}, //  Sm 62 --
      {0.3803922,1.000000,0.7803922}, //  Eu 63 --
      {0.2705882,1.000000,0.7803922}, //  Gd 64 --
      {0.1882353,1.000000,0.7803922}, //  Tb 65 --
      {0.1215686,1.000000,0.7803922}, //  Dy 66 --
      {0.000000,1.000000,0.6117647}, //  Ho 67 --
      {0.000000,0.9019608,0.4588235}, //  Er 68 --
      {0.000000,0.8313725,0.3215686}, //  Tm 69 --
      {0.000000,0.7490196,0.2196078}, //  Yb 70 --
      {0.000000,0.6705882,0.1411765}, //  Lu 71 --
      {0.3019608,0.7607843,1.000000}, //  Hf 72 --
      {0.3019608,0.6509804,1.000000}, //  Ta 73 --
      {0.1294118,0.5803922,0.8392157}, //  W  74 --
      {0.1490196,0.4901961,0.6705882}, //  Re 75 --
      {0.1490196,0.4000000,0.5882353}, //  Os 76 --
      {0.09019608,0.3294118,0.5294118}, //  Ir 77 --
      {0.8156863,0.8156863,0.8784314}, //  Pt 78 - changed from ghemical --
      {1.000000,0.8196078,0.1372549}, //  Au 79 - changed from ghemical --
      {0.7215686,0.7215686,0.8156863}, //  Hg 80 - changed from ghemical --
      {0.6509804,0.3294118,0.3019608}, //  Tl 81 --
      {0.3411765,0.3490196,0.3803922}, //  Pb 82 --
      {0.6196078,0.3098039,0.7098039}, //  Bi 83 --
      {0.6705882,0.3607843,0.000000}, //  Po 84 --
      {0.4588235,0.3098039,0.2705882}, //  At 85 --
      {0.2588235,0.5098039,0.5882353}, //  Rn 86 --
      {0.2588235,0.000000,0.4000000}, //  Fr 87 --
      {0.000000,0.4901961,0.000000}, //  Ra 88 --
      {0.4392157,0.6705882,0.9803922}, //  Ac 89 --
      {0.000000,0.7294118,1.000000}, //  Th 90 --
      {0.000000,0.6313725,1.000000}, //  Pa 91 --
      {0.000000,0.5607843,1.000000}, //  U  92 --
      {0.000000,0.5019608,1.000000}, //  Np 93 --
      {0.000000,0.4196078,1.000000}, //  Pu 94 --
      {0.3294118,0.3607843,0.9490196}, //  Am 95 --
      {0.4705882,0.3607843,0.8901961}, //  Cm 96 --
      {0.5411765,0.3098039,0.8901961}, //  Bk 97 --
      {0.6313725,0.2117647,0.8313725}, //  Cf 98 --
      {0.7019608,0.1215686,0.8313725}, //  Es 99 --
      {0.7019608,0.1215686,0.7294118}, //  Fm 100 --
      {0.7019608,0.05098039,0.6509804}, //  Md 101 --
      {0.7411765,0.05098039,0.5294118}, //  No 102 --
      {0.7803922,0.000000,0.4000000}, //  Lr 103 --
      {0.8000000,0.000000,0.3490196}, //  Rf 104 --
      {0.8196078,0.000000,0.3098039}, //  Db 105 --
      {0.8509804,0.000000,0.2705882}, //  Sg 106 --
      {0.8784314,0.000000,0.2196078}, //  Bh 107 --
      {0.9019608,0.000000,0.1803922}, //  Hs 108 --
      {0.9215686,0.000000,0.1490196} //  Mt 109 --
   };
   return atomicColor[nat][rgb];
}
solreal getAtomicRColorReal(int nat) {
   return getAtomicRGBColorReal(nat,0);
}
solreal getAtomicGColorReal(int nat) {
   return getAtomicRGBColorReal(nat,1);
}
solreal getAtomicBColorReal(int nat) {
   return getAtomicRGBColorReal(nat,2);
}
void getAtomicRGBColorsReal(int nat,solreal &rr,solreal &gg,solreal &bb) {
   rr=getAtomicRGBColorReal(nat,0);
   gg=getAtomicRGBColorReal(nat,1);
   bb=getAtomicRGBColorReal(nat,2);
   return;
}
void getAtomicRGBColorsInt(int nat,int &rr,int &gg,int &bb) {
   rr=getAtomicRGBColorInt(nat,0);
   gg=getAtomicRGBColorInt(nat,1);
   bb=getAtomicRGBColorInt(nat,2);
   return;
}

#endif /* defined(_COL_SCHEME_JMOL_H_) */
