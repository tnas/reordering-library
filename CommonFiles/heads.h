/*----------------------------------------------------------------------------
 * MATRIX STRUCTURE
 *--------------------------------------------------------------------------*/
#ifndef MATRIX_H
#define MATRIX_H

typedef struct
{
	double*     AA;
	double*      D;
	int*        JA;
	int*        IA;
	int     m,n,nz;
} MAT;

#endif /* MATRIX_H */

/*----------------------------------------------------------------------------
 * GRAPH STRUCTURE IN CSR FORMAT
 *--------------------------------------------------------------------------*/
#ifndef GRAPH_H
#define GRAPH_H

typedef struct 
{
	int label;
	int degree;
	int distance;
	int parent;
	int status;
	int chnum;
	int index;
} GRAPH;

#endif /* GRAPH_H */

/*----------------------------------------------------------------------------
 * ARRAY STRUCTURE
 *--------------------------------------------------------------------------*/
#ifndef ARRAY_H
#define ARRAY_H

typedef struct 
{
	double arr1;
	int    arr2;
	int    arr3;
} ARRAY;

#endif /* ARRAY_H */

/*----------------------------------------------------------------------------
 * LINKED LIST STRUCTURE
 *--------------------------------------------------------------------------*/
#ifndef LINKED_LIST_H
#define LINKED_LIST_H

typedef struct node
{
	int          data;
	int	     size;
	struct node* next;
} LIST;

#endif /* LINKED_LIST_H */

/*----------------------------------------------------------------------------
 * HSL LIBRARY STRUCTURES
 *--------------------------------------------------------------------------*/
#ifndef HSL_LIBRARY_H
#define HSL_LIBRARY_H

typedef enum MC60_ALGORITHM { SLOAN, RCM } MC60_ALGORITHM;

typedef enum MC60_CONTROL { 
	AUTOMATIC_PERIPHERAL, ESPECIFIED_PERIPHERAL, GLOBAL_PRIORITY_VECTOR } MC60_CONTROL;

#endif

/*----------------------------------------------------------------------------
 * PRECONDITIONERS STRUCTURES
 *--------------------------------------------------------------------------*/
#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

typedef struct
{
	int           n;
	int*    nzcount;  /* length of each row                          */
	int**        ja;  /* pointer-to-pointer to store column indices  */
	double**     ma;  /* pointer-to-pointer to store nonzero entries */
} SparMAT;

typedef struct
{
	int         n;
	SparMAT*    L;   /* L part elements   */
	double*     D;   /* diagonal elements */
	SparMAT*    U;   /* U part elements   */
	int*     work;   /* working buffer    */
} SparILU;

#endif /* PRECONDITIONERS_H */