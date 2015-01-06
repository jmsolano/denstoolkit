
#ifndef _SOLDEFINES_H_
#define _SOLDEFINES_H_


#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef PARALLELISEDTK
#define PARALLELISEDTK 0
#endif

#define USEPROGRESSBAR 1
#define CURRENTVERSION "0.1.5a"
#define PROGRAMCONTRIBUTORS "JMSA/JMHP"
#define EPSFORELFVALUE (2.871e-05)
#define DEFAULTPOINTSPERDIRECTION (80)

#define DISPLAYDEBUGINFOFILELINE (std::cout << __FILE__ << ", line: " << __LINE__ << std::endl)

#if DEBUG
#define _SOL_USE_SAFE_CHECKS_ 1
#else
#define _SOL_USE_SAFE_CHECKS_ 0
#endif

#define _HAVE_GNUPLOT_ 1
#define _HAVE_GRAPHICSMAGICK_ 1
#define _HAVE_IMAGEMAGICK_ 0
#define _HAVE_POVRAY_ 1
#define CHOOSEPOVVERSION36 1

#define _MAX_MEM_ALLOWANCE_ (1024*1024*1024)

#ifndef CMD_POVRAY
#define CMD_POVRAY "povray"
#endif

#endif /* defined(_SOLDEFINES_H_) */

