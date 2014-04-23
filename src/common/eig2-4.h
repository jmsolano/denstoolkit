/* 
 
 Eigen-decomposition code for symmetric 3x3 and 4x4 matrices.
 
 The code below is a modification from the code found in 
 
 http://barnesc.blogspot.ca/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
 
 In turn, the above modification has the following text:
 
 Eigen decomposition code for symmetric 3x3 matrices, copied from the public
 domain Java Matrix library JAMA.
 
 */

#ifndef _EIG2_4_H_
#define _EIG2_4_H_

#ifndef _HAVE_DEF_SOLREAL_TYPE_
#define _HAVE_DEF_SOLREAL_TYPE_
typedef double solreal;
//typedef float solreal;
#endif

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
 eigenvalues in d. */
void eigen_decomposition4(solreal A[4][4], solreal V[4][4], solreal d[4]);

void eigen_decomposition3(solreal A[3][3], solreal V[3][3], solreal d[3]);

void eigen_decomposition2(solreal A[2][2], solreal V[2][2], solreal d[2]);

#endif//_EIG2_4_H_

