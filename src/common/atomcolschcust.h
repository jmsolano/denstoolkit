

#ifndef _COL_SCHEME_JMHP_H_
#define _COL_SCHEME_JMHP_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#define _HAVE_SELECTED_ATOM_PALETTE_ 1

#define MAXDEFINEDATOMICCOLORS 19


/** This table contains the colors for the atoms. The numbers are RGB integers (\f$0-255\f$)
 */
static int atomicColorInt[MAXDEFINEDATOMICCOLORS][3]={
   {255, 255, 255},
   {255, 192, 203},
   {178, 34, 34},
   {0, 255, 0},
   {211, 211, 211},
   {143, 143, 255},
   {255, 0, 0},
   {218, 165, 32},
   {0, 0, 255},
   {34, 139, 34},
   {0, 0, 255},
   {34, 139, 34},
   {128, 128, 144},
   {218, 165, 32},
   {255, 165, 0},
   {255, 200, 50},
   {0, 255, 0},
   {255, 20, 147},
   {255, 20, 147}
};

/** This table contains the colors for the atoms. The numbers are RGB floats (\f$0.0-1.0\f$) */
static solreal atomicColor[MAXDEFINEDATOMICCOLORS][3]={
   {1.00000, 1.00000, 1.00000},
   {1.00000, 0.75294, 0.79608},
   {0.69804, 0.13333, 0.13333},
   {0.00000, 1.00000, 0.00000},
   {0.82745, 0.82745, 0.82745},
   {0.56078, 0.56078, 1.00000},
   {1.00000, 0.00000, 0.00000},
   {0.85490, 0.64706, 0.12549},
   {0.00000, 0.00000, 1.00000},
   {0.13333, 0.54510, 0.13333},
   {0.00000, 0.00000, 1.00000},
   {0.13333, 0.54510, 0.13333},
   {0.50196, 0.50196, 0.56471},
   {0.85490, 0.64706, 0.12549},
   {1.00000, 0.64706, 0.00000},
   {1.00000, 0.78431, 0.19608},
   {0.00000, 1.00000, 0.00000},
   {1.00000, 0.07843, 0.57647},
   {1.00000, 0.07843, 0.57647}
};

#endif /* defined(_COL_SCHEME_JMHP_H_) */
