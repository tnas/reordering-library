#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef struct
{
	double*     AA;
	double*      D;
	int*        JA;
	int*        IA;
	int     m,n,nz;
} MAT;

typedef struct 
{
	double arr1;
	int    arr2;
	int    arr3;
} ARRAY;

/*----------------------------------------------------------------------------
 * MATRIX HEADER FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern void     MATRIX_readCSR           (MAT* A, FILE* f);
extern double   MATRIX_aij               (MAT* A, int i, int j);
extern void     MATRIX_printCSR          (MAT* A);
extern void     MATRIX_printFULL         (MAT* A);
extern long int MATRIX_envelope          (MAT* A);
extern long int MATRIX_bandwidth         (MAT* A);
extern void     MATRIX_clean             (MAT* A);
extern void     MATRIX_matvec            (MAT* A, double* x, double* b);
extern void     MATRIX_forward           (MAT* L, double* b, double* y);
extern void     MATRIX_backward          (MAT* U, double* y, double* x);
extern void     MATRIX_permutation       (MAT* A, int* p);
extern void     MATRIX_writeCSR          (MAT* A, double* f, int* s, int nP, int bandwidth);
extern void 	MATRIX_readCSR_SymmUpper (MAT* A, FILE* f);
extern long int MATRIX_wavefront	 (MAT* A);

extern int      COMPARE_array            (const void * a, const void * b);
extern int      COMPARE_eig              (const void * a, const void * b);
extern 	void 	write_output_after	 (const MAT *A);
extern 	void 	write_output_before	 (const MAT *A);

#endif