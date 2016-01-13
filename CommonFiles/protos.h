#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "heads.h"

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

typedef enum LABEL { UNREACHED, LABELED } LABEL;

/*----------------------------------------------------------------------------
 * ALGEBRA LINEAR FUNCTIONS PROTOTYPE 
 *--------------------------------------------------------------------------*/
extern int      daxpy                    (int n, double  a, double *x, double *y);
extern int      dcopy                    (int n, double *x, double *y);
extern double   ddot                     (int n, double *x, double *y);
extern int      dscal                    (int n, double  a, double *x, double *y);
extern int      ddiff                    (int n, double *a, double *x, double *y);
extern int      izero                    (int n, int *v);
extern int      iszero                   (double x);

/*----------------------------------------------------------------------------
 * COMPARE FUNCTIONS FOR QSORT 
 *--------------------------------------------------------------------------*/
extern int      COMPARE_array            (const void * a, const void * b);
extern int      COMPARE_eig              (const void * a, const void * b);
extern int      COMPARE_degr_ASC         (const void * a, const void * b);
extern int      COMPARE_dist_degr_DES    (const void * a, const void * b);
extern int      COMPARE_dist_degr_ASC    (const void * a, const void * b);

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

/*----------------------------------------------------------------------------
 * GRAPH FUNCTIONS PROTOTYPE IN CSR FORMAT
 *--------------------------------------------------------------------------*/
extern int      GRAPH_degree             (MAT* A, int x);
extern int*     GRAPH_adjacent           (MAT* A, int x);
extern int 	GRAPH_degree_per_level   (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors);
extern GRAPH* 	GRAPH_adjacent_per_level (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors);
extern void     GRAPH_bfs                (MAT* A, int x, int* dist);
extern int*     GRAPH_bfs_RCM            (MAT* A, int x, int* dist);
extern int      GRAPH_LS_depth           (int* LS, int n);
extern int      GRAPH_LS_width           (int* LS, int n);
extern LIST*    GRAPH_LS_last_level      (MAT* A, int* LS, int n);
extern int*     GRAPH_LS_peripheral      (MAT* A, int *node_s, int* node_e);
extern int* 	GRAPH_fixedpoint_bfs	 (MAT* adjacency, int root, int* levels);


/*----------------------------------------------------------------------------
 * LINKED LIST FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern LIST*    LIST_insert_IF_NOT_EXIST (LIST* L, int x);
extern LIST*    LIST_remove              (LIST* L, int x);
extern LIST*    LIST_remove_first        (LIST* L);
extern void     LIST_print               (LIST* L);
extern int      LIST_first               (LIST* L);
extern void     LIST_destroy             (LIST* L);
extern LIST*    LIST_insert		 (LIST* L, int x);

/*----------------------------------------------------------------------------
 * SOLVERS FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern double   GMRES_algorithm          (MAT* A, double *x, double *b, int n, int k, int *iter, double eps);
extern void     GMRES                    (MAT* A, double *x, double *b, double tol, int restart, int maxiter);
extern double   PGMRES_algorithm         (MAT* A, MAT* L, MAT* U, double *x, double *b, int n, int k, int *iter, double eps);
extern void     PGMRES                   (MAT* A, MAT* L, MAT* U, double *x, double *b, double tol, int restart, int maxiter);

/*----------------------------------------------------------------------------
 * PRECONDITIONERS FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern void     SPARMAT_setup            (SparMAT* mat, int n);
extern void     SPARILU_setup            (SparILU* lu, int n);
extern void     SPARILU_row              (SparILU* lu, int nrow);
extern void     CSR_to_SPARMAT           (MAT* A, SparMAT* mat);
extern void     SPARILU_to_CSR           (SparILU* lu, MAT* L, MAT* U);
extern int      LEVEL_OF_FILL            (SparMAT* csmat, SparILU* lu, int p);
extern void     ILUP                     (SparMAT* csmat, SparILU* lu, int p);
extern void     SPARMAT_clean            (SparMAT* mat);
extern void     SPARILU_clean            (SparILU* lu);
extern void     SPARILU_print            (SparILU* lu);
extern void     QSPLIT                   (double *a, int *ind, int n, int Ncut);
extern void     ILUT                     (SparMAT* csmat, SparILU* lu, int lfil, double tol);

/*----------------------------------------------------------------------------
 * REORDERING FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern void     REORDERING_RCM_opt       (MAT* A, int** p, int s);
extern void     REORDERING_RCM           (MAT* A, int** p);
extern void     REORDERING_SLOAN         (MAT* A, int** p);
extern void     REORDERING_SPECTRAL      (MAT* A, int** p);
extern void     REORDERING_SPECTRAL_WGT  (MAT* A, int** p);