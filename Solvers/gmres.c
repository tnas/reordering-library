/*----------------------------------------------------------------------------
 * GMRES SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

/*----------------------------------------------------------------------------
 * GMRES algorithm
 *--------------------------------------------------------------------------*/
double GMRES_algorithm (MAT* A, double *x, double *b, int n, int k, int *iter, double eps)
{ 
	int i, j, hz = k + 1;  

	double * ax = calloc(n,sizeof(double));
	double * r  = calloc(n,sizeof(double));
  
	double e0;
	double beta;
	double ret;

	MATRIX_matvec(A,x,ax);
	ddiff(n,b,ax,r);
  
	e0 = sqrt(ddot(n,r,r));

	if (e0 <= eps)
	{
		ret = e0;
		free(ax);
		free(r); 
	} 
 
	double** q = malloc(hz*sizeof(double*));
	for (i = 0; i < hz; ++i) 
		q[i] = malloc(n*sizeof(double));
  
	double** h = calloc(hz,sizeof(double*));
	for (i = 0; i < hz; ++i) 
		h[i] = calloc(hz,sizeof(double));
  
	double* e = malloc(hz*sizeof(double)); 
	double* s = malloc(hz*sizeof(double));
	double* c = malloc(hz*sizeof(double));
  
	dscal(n,1.0/e0,r,q[0]); 
  
	e[0]  = e0;
	*iter = 0;
  
	for (j = 0; e0 > eps && j < k; ++j)
	{
		*iter= *iter +1;

		MATRIX_matvec (A,q[j],ax);
    
		double nr1 = sqrt(ddot(n,ax,ax));

		for (i = 0; i <= j; ++i)
		{
			h[i][j] = ddot(n,q[i],ax);
			daxpy(n,-h[i][j],q[i],ax);
		}
    
		h[j+1][j] = sqrt(ddot(n,ax,ax));
		double nr2 = 0.001 * h[j+1][j] + nr1;
    
		if (fabs(nr2  - nr1) < eps)
		{
			for (i = 0; i <= j; ++i)
			{
				double hr = ddot(n,q[i],ax);
				h[i][j] += hr;
				daxpy(n,-hr,q[i],ax);
			}
			h[j+1][j] = sqrt(ddot(n,ax,ax));
		}

		for (i = 0; i <= j - 1; ++i)
		{
			double x  = h[i][j];
			double y  = h[i+1][j];
			h[i][j]   = x * c[i+1] + y * s[i+1];
			h[i+1][j] = x * s[i+1] - y * c[i+1];
		}

		beta    = sqrt(h[j][j]*h[j][j] + h[j+1][j]*h[j+1][j]);
		s[j+1]  = h[j+1][j] / beta;
		c[j+1]  = h[j][j]   / beta;
		h[j][j] = beta;
		e[j+1]  = s[j+1] * e[j];
		e[j]    = c[j+1] * e[j];
		e0      = e[j+1];
   
		dscal(n,1.0 / h[j+1][j],ax,q[j+1]); 
	}
	--j;
  
	ret = e[j+1];
	double * y = malloc(hz * sizeof(double));
	for (i = j; i >= 0; --i) 
	{
		double sum = 0.0;
		for (k = i + 1; k <= j; ++k)
			sum += h[i][k] * y[k];
		y[i] = (e[i] - sum) / h[i][i];
	}
	for (i = 0; i <= j; ++i) 
		daxpy(n,y[i],q[i],x);
	free(y);
  
	/* Free up memory */
	for (i=0;i<hz;++i) 
		free(q[i]); free(q); 
	for (i=0;i<hz;++i) 
		free(h[i]); free(h);
	free(ax);
	free(r); 
	free(e);
	free(s); 
	free(c);
	return ret;
}

/*----------------------------------------------------------------------------
 * GMRES algorithm call
 *--------------------------------------------------------------------------*/
void GMRES (MAT* A, double *x, double *b, double tol, int restart, int maxiter)
{  
	int    i, it;  
	int    n  = A->n;
	double bn = sqrt(ddot(n,b,b));
	double e;

	fprintf(stderr, "\n ||b||   : %le\n", bn);

	if (bn < tol) 
		return;

	for (i = 0; i < maxiter; i++)
	{
		e  = GMRES_algorithm (A, x, b, n, restart, &it, tol * bn);
		e /= bn;
		fprintf(stderr, " ITERS   : %d - %d, eps = %le\n", i+1,i*restart + it, e);
		if (e < tol)
			return;
	}      
}

/*----------------------------------------------------------------------------
 * Preconditioned GMRES algorithm
 *--------------------------------------------------------------------------*/
double PGMRES_algorithm (MAT* A, MAT* L, MAT* U, double *x, double *b, int n, int k, int *iter, double eps)
{ 
	int i, j, hz = k + 1;  

	double * ax = calloc(n,sizeof(double));
	double * r  = calloc(n,sizeof(double));
	double * p1 = calloc(n,sizeof(double));
	double * p2 = calloc(n,sizeof(double));
  
	double e0;
	double beta;
	double ret;

	MATRIX_matvec  (A,x,p1);
	/* Preconditioned Operations in ax */
	MATRIX_forward  (L,p1,p2);  /* Lp2 = p1 */
	MATRIX_backward (U,p2,ax);  /* Uax = p2 */

	ddiff (n,b,ax,r);
  
	e0 = sqrt (ddot(n,r,r));

	if (e0 <= eps)
	{
		free(ax);
		free(r);
		free(p1);
		free(p2);
		return e0;
	} 
 
	double** q = malloc(hz*sizeof(double*));
	for (i = 0; i < hz; ++i) 
		q[i] = malloc(n*sizeof(double));
  
	double** h = calloc(hz,sizeof(double*));
	for (i = 0; i < hz; ++i) 
		h[i] = calloc(hz,sizeof(double));
  
	double* e = malloc(hz*sizeof(double)); 
	double* s = malloc(hz*sizeof(double));
	double* c = malloc(hz*sizeof(double));
  
	dscal (n,1.0/e0,r,q[0]); 
  
	e[0]  = e0;
	*iter = 0;
  
	for (j = 0; e0 > eps && j < k; ++j)
	{
		*iter= *iter +1;

		MATRIX_matvec  (A,q[j],p1);
		/* Preconditioned Operations in ax */
		MATRIX_forward  (L,p1,p2); /* Lp2 = p1 */
		MATRIX_backward (U,p2,ax); /* Uax = p2 */
    
		double nr1 = sqrt (ddot(n,ax,ax));

		for (i = 0; i <= j; ++i)
		{
			h[i][j] = ddot (n,q[i],ax);
			daxpy (n,-h[i][j],q[i],ax);
		}
    
		h[j+1][j] = sqrt(ddot (n,ax,ax));
		double nr2 = 0.001 * h[j+1][j] + nr1;
    
		if (fabs(nr2  - nr1) < eps)
		{
			for (i = 0; i <= j; ++i)
			{
				double hr = ddot (n,q[i],ax);
				h[i][j] += hr;
				daxpy (n,-hr,q[i],ax);
			}
			h[j+1][j] = sqrt(ddot (n,ax,ax));
		}

		for (i = 0; i <= j - 1; ++i)
		{
			double x  = h[i][j];
			double y  = h[i+1][j];
			h[i][j]   = x * c[i+1] + y * s[i+1];
			h[i+1][j] = x * s[i+1] - y * c[i+1];
		}

		beta    = sqrt(h[j][j]*h[j][j] + h[j+1][j]*h[j+1][j]);
		s[j+1]  = h[j+1][j] / beta;
		c[j+1]  = h[j][j]   / beta;
		h[j][j] = beta;
		e[j+1]  = s[j+1] * e[j];
		e[j]    = c[j+1] * e[j];
		e0      = e[j+1];
   
		dscal(n,1.0 / h[j+1][j],ax,q[j+1]); 
	}
	--j;
  
	ret = e[j+1];
	double * y = malloc(hz * sizeof(double));
	for (i = j; i >= 0; --i) 
	{
		double sum = 0.0;
		for (k = i + 1; k <= j; ++k)
			sum += h[i][k] * y[k];
		y[i] = (e[i] - sum) / h[i][i];
	}
	for (i = 0; i <= j; ++i) 
		daxpy (n,y[i],q[i],x);
	free(y);
  
	/* Free up memory */
	for (i=0;i<hz;++i) 
		free(q[i]); free(q); 
	for (i=0;i<hz;++i) 
		free(h[i]); free(h);
	free(ax);
	free(r); 
	free(e);
	free(s); 
	free(c);
	free(p1);
	free(p2);
	return ret;
}

/*----------------------------------------------------------------------------
 * Preconditioned GMRES algorithm call
 *--------------------------------------------------------------------------*/
void PGMRES (MAT* A, MAT* L, MAT* U, double *x, double *b, double tol, int restart, int maxiter)
{  
	int     i, it;  
	int     n = A->n;
	double* c = calloc (n,sizeof(double));

	/* Preconditioned Operations in b */
	MATRIX_forward  (L,b,c); /* Lc = b */
	MATRIX_backward (U,c,b); /* Ub = c */
	free(c);
	
	double bn = sqrt(ddot(n,b,b));
	double e;

	fprintf(stderr, "\n ||b||   : %le\n", bn);

	if (bn < tol) 
		return;

	for (i = 0; i < maxiter; i++)
	{
		e  = PGMRES_algorithm (A, L, U, x, b, n, restart, &it, tol * bn);
		e /= bn;
		fprintf(stderr, " ITERS   : %d - %d, eps = %le\n", i+1,i*restart + it, e);
		if (e < tol)
			return;
	}      
}