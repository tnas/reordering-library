/*----------------------------------------------------------------------------
 * ILUP PRECONDITIONER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

/*----------------------------------------------------------------------
|     does a quick-sort split of a real array.
|     on input a[0 : (n-1)] is a real array
|     on output is permuted such that its elements satisfy:
|
|     abs(a[i]) >= abs(a[Ncut-1]) for i < Ncut-1 and
|     abs(a[i]) <= abs(a[Ncut-1]) for i > Ncut-1
|
|     ind[0 : (n-1)] is an integer array permuted in the same way as a.
|---------------------------------------------------------------------*/
void QSPLIT (double *a, int *ind, int n, int Ncut)
{
	
	double tmp, abskey;
	int j, itmp, first, mid, last, ncut;
	ncut = Ncut - 1;

	first = 0;
	last = n-1;
	if (ncut<first || ncut>last)
		return;
	
	/* outer loop -- while mid != ncut */
	do
	{
		mid    = first;
		abskey = fabs(a[mid]);
		for (j=first+1; j<=last; j++)
		{
			if (fabs(a[j]) > abskey)
			{
				mid      = mid+1;
				tmp      = a[mid];
				itmp     = ind[mid];
				a[mid]   = a[j];
				ind[mid] = ind[j];
				a[j]     = tmp;
				ind[j]   = itmp;
			}
		}
		/*-------------------- interchange */
		tmp        = a[mid];
		a[mid]     = a[first];
		a[first]   = tmp;
		itmp       = ind[mid];
		ind[mid]   = ind[first];
		ind[first] = itmp;
		/*-------------------- test for while loop */
		if (mid == ncut) break;
		if (mid > ncut) 
			last = mid-1;
		else
			first = mid+1;
	}while (mid != ncut);
	
	return;
}
	
/*----------------------------------------------------------------------------
 * ILUT preconditioner
 * incomplete LU factorization with dual truncation mechanism
 *--------------------------------------------------------------------------*/
void ILUT (SparMAT* csmat, SparILU* lu, int lfil, double tol)
{
	int n = csmat->n; 
	int len, lenu, lenl;
	int nzcount, *ja, *jbuf, *iw, i, j, k;
	int col, jpos, jrow, upos;
	double t, tnorm, tolnorm, fact, lxu, *wn, *ma, *w;
	SparMAT* L;
	SparMAT* U;
	double *D;

	if (lfil < 0)
	{
		printf("ILUT: Illegal value for lfil.\n");
		return;
	}    

	SPARILU_setup (lu,n);
	L = lu->L;
	U = lu->U;
	D = lu->D;

	iw   = (int*)    malloc(n*sizeof(int));
	jbuf = (int*)    malloc(n*sizeof(int));
	wn   = (double*) malloc(n*sizeof(double));
	w    = (double*) malloc(n*sizeof(double));  

	/* set indicator array jw to -1 */
	for (i = 0; i < n; i++)
		iw[i] = -1;

	/* beginning of main loop */
	for (i = 0; i < n; i++)
	{
		nzcount = csmat->nzcount[i];
		ja      = csmat->ja[i];
		ma      = csmat->ma[i];
		tnorm   = 0;
		
		for (j = 0; j < nzcount; j++)
			tnorm += fabs( ma[j] );

		if (tnorm == 0.0)
		{
			printf ("ILUT: zero row encountered.\n");
			return;
		}
		
		tnorm  /= (double)nzcount;
		tolnorm = tol * tnorm;

		/* unpack L-part and U-part of column of A in arrays w */
		lenu    = 0;
		lenl    = 0;
		jbuf[i] = i;
		w[i]    = 0;
		iw[i]   = i;
		for (j = 0; j < nzcount; j++)
		{
			col = ja[j];
			t   = ma[j];
			if (col < i)
			{
				iw[col]    = lenl;
				jbuf[lenl] = col;
				w[lenl]    = t;
				lenl++;
			} else if (col == i)
			{
				w[i] = t;
			} else
			{
				lenu++;
				jpos       = i + lenu;
				iw[col]    = jpos;
				jbuf[jpos] = col;
				w[jpos]    = t;
			}
		}

		j   = -1;
		len = 0;
		/* eliminate previous rows */
		while (++j < lenl)
		{
			/*----------------------------------------------------------------------------
			*  in order to do the elimination in the correct order we must select the
			*  smallest column index among jbuf[k], k = j+1, ..., lenl
			*--------------------------------------------------------------------------*/
			jrow = jbuf[j];
			jpos = j;
			/* determine smallest column index */
			for(k = j + 1; k < lenl; k++)
			{
				if (jbuf[k] < jrow)
				{
					jrow = jbuf[k];
					jpos = k;
				}
			}
			if (jpos != j)
			{
				col        = jbuf[j];
				jbuf[j]    = jbuf[jpos];
				jbuf[jpos] = col;
				iw[jrow]   = j;
				iw[col]    = jpos;
				t          = w[j];
				w[j]       = w[jpos];
				w[jpos]    = t;
			}

			/* get the multiplier */
			fact = w[j] * D[jrow];
			w[j] = fact;
			/* zero out element in row by resetting iw(n+jrow) to -1 */
			iw[jrow] = -1;

			/* combine current row and row jrow */
			nzcount = U->nzcount[jrow];
			ja = U->ja[jrow];
			ma = U->ma[jrow];
			for (k = 0; k < nzcount; k++)
			{
				col  = ja[k];
				jpos = iw[col];
				lxu  = - fact * ma[k];
				/* if fill-in element is small then disregard */
				if (fabs(lxu) < tolnorm && jpos == -1) 
					continue;

				if (col < i) 
				{
				/* dealing with lower part */
					if (jpos == -1)
					{
						/* this is a fill-in element */
						jbuf[lenl] = col;
						iw[col]    = lenl;
						w[lenl]    = lxu;
						lenl++;
					} else
					{
						w[jpos] += lxu;
					}
				} else
				{
					/* dealing with upper part */
					if (jpos == -1 && fabs(lxu) > tolnorm)
					{
					/* this is a fill-in element */
						lenu++;
						upos = i + lenu;
						jbuf[upos] = col;
						iw[col] = upos;
						w[upos] = lxu;
					} else
					{
						w[jpos] += lxu;
					}
				}
			}
		}

		/* restore iw */
		iw[i] = -1;
		for (j = 0; j < lenu; j++)
			iw[jbuf[i+j+1]] = -1;

		/*---------- case when diagonal is zero */
		if (w[i] == 0.0)
		{
			printf("zero diagonal encountered.\n");
			for (j = i; j < n; j++)
			{
				L->ja[j] = NULL; 
				L->ma[j] = NULL;
				U->ja[j] = NULL; 
				U->ma[j] = NULL;
			}
			return;
		}
		/*-----------Update diagonal */    
		D[i] = 1 / w[i];

		/* update L-matrix */
		len = lenl < lfil ? lenl : lfil;
		for (j = 0; j < lenl; j++)
		{
			wn[j] = fabs( w[j] );
			iw[j] = j;
		}
		QSPLIT (wn, iw, lenl, len);
		L->nzcount[i] = len;
		if(len > 0) 
		{
			ja = L->ja[i] = (int*)    malloc(len*sizeof(int));
			ma = L->ma[i] = (double*) malloc(len*sizeof(double));
		}
		for (j = 0; j < len; j++)
		{
			jpos  = iw[j];
			ja[j] = jbuf[jpos];
			ma[j] = w[jpos];
		}
		for(j = 0; j < lenl; j++) 
			iw[j] = -1;

		/* update U-matrix */
		len = lenu < lfil ? lenu : lfil;
		for (j = 0; j < lenu; j++)
		{
			wn[j] = fabs( w[i+j+1] );
			iw[j] = i+j+1;
		}
		QSPLIT (wn, iw, lenu, len);
		U->nzcount[i] = len;
		if (len > 0)
		{
			ja = U->ja[i] = (int*)    malloc(len*sizeof(int));
			ma = U->ma[i] = (double*) malloc(len*sizeof(double));
		}
		for (j = 0; j < len; j++)
		{
			jpos  = iw[j];
			ja[j] = jbuf[jpos];
			ma[j] = w[jpos];
		}
		for (j = 0; j < lenu; j++)
			iw[j] = -1;
	}

	free(iw);
	free(jbuf);
	free(wn);
	free(w);
	
	return;
}