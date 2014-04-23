
#ifndef _SOLDEFINES_H_
#define _SOLDEFINES_H_


#ifndef _HAVE_DEF_REAL_TYPE_
#define _HAVE_DEF_REAL_TYPE_
typedef double real;
typedef double solreal;
//typedef float real;
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef PARALLELISEDTK
#define PARALLELISEDTK 0
#endif

#define USEPROGRESSBAR 1
#define CURRENTVERSION "1.0.1"
#define PROGRAMCONTRIBUTORS "JMSA/JMHP"
#define EPSFORELFVALUE (2.871e-05)
#define DEFAULTPOINTSPERDIRECTION (200)
#define DEFAULTNPOINTSFORCUBE (80)
#define DEFAULTMAXVALUEOFP (5.0e0)

#define DISPLAYDEBUGINFOFILELINE (std::cout << __FILE__ << ", line: " << __LINE__ << std::endl)

#if DEBUG
#define _SOL_USE_SAFE_CHECKS_ 1
#else
#define _SOL_USE_SAFE_CHECKS_ 0
#endif

#define _HAVE_GNUPLOT_ 1

#define _HAVE_EPSTOOL_ 1

#define _HAVE_EPSTOPDF_ 1

#define _MAX_MEM_ALLOWANCE_ (1024*1024*1024)

#endif /* defined(_SOLDEFINES_H_) */