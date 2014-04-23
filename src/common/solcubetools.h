/*
 *  solcubetools.h
 *  
 *
 *  Created by Juan Manuel Solano Altamirano on 07/05/13.
 *  Copyright 2013. All rights reserved.
 *
 */

#ifndef _SOLCUBETOOLS_H_
#define _SOLCUBETOOLS_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <fstream>
using std::fstream;
using std::ofstream;
#include <cstdlib>
using std::exit;
#include <string>
using namespace std;
#include <iomanip>
using std::setprecision;
using std::scientific;

//*************************************************************************************************
//*************************************************************************************************
void writeCubeHeader(ofstream &ofil,string &t1,string &t2,int (&bdim)[3],
                solreal (&x0)[3],solreal (&dx)[3][3],int nat,solreal* (&atchrg),solreal* (&x));
//**********************************************************************************************
void writeCubeProp(ofstream &ofil,int dim,solreal* (&prop));
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
//*************************************************************************************************
#endif//_SOLCUBETOOLS_H_

