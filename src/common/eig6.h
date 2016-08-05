

/*
 
 Eigen-decomposition code for symmetric 6x6 matrices.
 
 The code below is a modification from the code found in 
 
 http://barnesc.blogspot.ca/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
 
 In turn, the above modification has the following text:
 
 Eigen decomposition code for symmetric 6x6 matrices, copied from the public
 domain Java Matrix library JAMA.
 
 */

#ifndef _EIG6_H_
#define _EIG6_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
 eigenvalues in d. */
void eigen_decomposition6(solreal A[6][6], solreal V[6][6], solreal d[6]);

#endif//_EIG6_H_

