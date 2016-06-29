#include "spectral.h"

/*----------------------------------------------------------------------------
 * SPECTRAL reordering
 *--------------------------------------------------------------------------*/
void REORDERING_HSL_SPECTRAL (MAT* A, int** Fp)
{
	int i;
	int n = A->n;
	
	double *fvector = calloc (n ,sizeof(double));
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int    *p       = calloc (n ,sizeof(int));
	int    *ia = A->IA;
	int    *ja = A->JA;
	int   lirn = A->nz;
	int    nnz = A->nz;
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,NULL);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	
	ARRAY* R = malloc(n*sizeof(ARRAY));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}
	
	
	qsort (R,n,sizeof(ARRAY),COMPARE_eig); 
	
	for (i = 0; i < n; ++i) 
		p[i] = R[i].arr2; 

	free(R);
	free(info);
	free(list);
	free(fvector);
	
	(*Fp) = p;
}

void REORDERING_HSL_SPECTRAL_WGT (MAT* A, int** Fp)
{
	int i;
	int n = A->n;
	
	double *fvector = calloc (n ,sizeof(double));
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int    *p       = calloc (n ,sizeof(int));
	int    *ia = A->IA;
	int    *ja = A->JA;
	double *a  = A->AA;
	int   lirn = A->nz;
	int    nnz = A->nz;
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,a);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	
	ARRAY* R = malloc(n*sizeof(ARRAY));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}
	
	
	qsort (R,n,sizeof(ARRAY),COMPARE_eig); 
	
	for (i = 0; i < n; ++i) 
		p[i] = R[i].arr2; 

	free(R);
	free(info);
	free(list);
	free(fvector);
	
	(*Fp) = p;
}
