/*----------------------------------------------------------------------------
 * ALGEBRA LINEAR OPERATIONS
 *--------------------------------------------------------------------------*/
#include "protos.h"

/*----------------------------------------------------------------------------
 * Constant times a vector plus a vector
 *--------------------------------------------------------------------------*/
inline int daxpy (int n, double a, double *x, double *y)
{
	long int i, m;
	register double sa;

	sa = a;
	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   += sa * x[i];
		y[i+1] += sa * x[i+1];
		y[i+2] += sa * x[i+2];
		y[i+3] += sa * x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i]   += sa * x[i];
	
	return 0;
}

/*----------------------------------------------------------------------------
 * Copies a vector x to a vector y
 *--------------------------------------------------------------------------*/
inline int dcopy (int n, double *x, double *y)
{
	long int i, m;

	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   = x[i];
		y[i+1] = x[i+1];
		y[i+2] = x[i+2];
		y[i+3] = x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i]   = x[i];
	
	return 0;
}

/*----------------------------------------------------------------------------
 * Forms the dot product of two vectors
 *--------------------------------------------------------------------------*/
inline double ddot (int n, double *x, double *y)
{
	long int i, m;
	double stemp;

	stemp = 0.0;
	m = n-4;

	for (i = 0; i < m; i += 5)
		stemp += x[i] * y[i] + x[i+1] * y[i+1] + x[i+2] * y[i+2] + x[i+3] * y[i+3] + x[i+4] * y[i+4];

	for ( ; i < n; i++)        /* clean-up loop */
		stemp += x[i] * y[i];

	return stemp;
} 

/*----------------------------------------------------------------------------
 * Scales a vector by a constant
 *--------------------------------------------------------------------------*/
inline int dscal (int n, double a, double *x, double* y)
{
	long int i, m;
	register double sa;

	sa = a;
	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   = sa * x[i];
		y[i+1] = sa * x[i+1];
		y[i+2] = sa * x[i+2];
		y[i+3] = sa * x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i]   = sa * x[i];

	return 0;
}

/*----------------------------------------------------------------------------
 * Erase a vector of doubles
 *--------------------------------------------------------------------------*/
inline int ddiff (int n, double *a, double *x, double *y)
{
	long int i, m;

	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   = a[i]   - x[i];
		y[i+1] = a[i+1] - x[i+1];
		y[i+2] = a[i+2] - x[i+2];
		y[i+3] = a[i+3] - x[i+3];
	}
	for ( ; i < n; ++i) 
		y[i]   = a[i]   - x[i];;

	return 0;
}

/*----------------------------------------------------------------------------
 * Erase a vector of integers
 *--------------------------------------------------------------------------*/
inline int izero (int n, int *v)
{
	long int i, m;

	m = n-3;
	for (i = 0; i < m; i += 4){
		v[i]   = 0;
		v[i+1] = 0;
		v[i+2] = 0;
		v[i+3] = 0;
	}
	for ( ; i < n; ++i) 
		v[i]   = 0;

	return 0;
}

/*----------------------------------------------------------------------------
 * Test if an element is zero
 *--------------------------------------------------------------------------*/
int iszero (double x)
{
	return (x < 1e-15 && x > -(1e-15));
}