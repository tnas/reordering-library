/*----------------------------------------------------------------------------
 * MATRIX FUNCTIONS: IN CSR FORMAT
 *--------------------------------------------------------------------------*/
#include "protos.h"


int COMPARE_array (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr3 <  ((ARRAY*)b)->arr3) return -1;
	if (((ARRAY*)a)->arr3 >  ((ARRAY*)b)->arr3) return  1;
	if (((ARRAY*)a)->arr3 == ((ARRAY*)b)->arr3)
	{
		if (((ARRAY*)a)->arr2 < ((ARRAY*)b)->arr2) return -1;
		if (((ARRAY*)a)->arr2 > ((ARRAY*)b)->arr2) return  1;
	}
	return 0;
}


/*----------------------------------------------------------------------------
 * Read matrix from file in MM format to a CSR structure
 *--------------------------------------------------------------------------*/
void MATRIX_readCSR (MAT* A, FILE* f)
{
	int M, N, nz;
	int i, j, k, I, J, elem = 0;
	double VAL;
	char line[1025], header[1025];
	char *symm;
	ARRAY* a;
	
	fgets(header,200,f);	
	strtok(header," \n");
	strtok(NULL," \n");
	strtok(NULL," \n");
	strtok(NULL," \n");
	
	symm  = strtok(NULL," \n");
	
	if ((strcmp(symm,"symmetric") != 0 && strcmp(symm,"general") != 0))
		fprintf (stderr,"\n Error. Matrix format not supported. Exiting.. [MATRIX_readCSR]\n\n");
	
	do 
	{
		if (fgets(line,1025,f) == NULL) 
			exit(0);
	}while (line[0] == '%');   
	sscanf(line,"%d %d %d", &N, &M, &nz);
	
	if (strcmp(symm,"symmetric") == 0)
	{
		a = malloc ((2*nz)*sizeof(ARRAY));
		for (i = 0, k = 0; i < nz; ++i)
		{
			fscanf(f,"%d%d%lf",&I,&J,&VAL);
			a[k].arr1 = VAL;
			a[k].arr2 = J - 1;
			a[k].arr3 = I - 1;
			k++;
			if (I != J)
			{
				a[k].arr1 = VAL;
				a[k].arr2 = I - 1;
				a[k].arr3 = J - 1;
				k++;
			}
		}
		nz = k;
		a  = realloc (a,nz*sizeof(ARRAY));
		qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	}
	
	if (strcmp(symm,"general") == 0)
	{
		a = malloc (nz*sizeof(ARRAY));
		for (i = 0; i < nz; ++i)
		{
			fscanf(f,"%d%d%lf",&I,&J,&VAL);
			a[i].arr1 = VAL;
			a[i].arr2 = J - 1;
			a[i].arr3 = I - 1;
		}
		qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	}

	/* reseve memory for matrices */
	A->m   = M;
	A->n   = N;
	A->nz  = nz;
	
	A->AA  = (double *) calloc(nz, sizeof (double));
	A->D   = (double *) calloc(N , sizeof (double));
	A->JA  = (int    *) calloc(nz, sizeof (int));
	A->IA  = (int    *) calloc(N+1,sizeof (int));
	
	for (i = 0; i < nz; ++i)
	{
		A->AA[i]   = a[i].arr1;
		A->JA[i]   = a[i].arr2;
		elem      += 1;
		A->IA[a[i].arr3+1] = elem;
	}

	free(a);
	
	/* Adjusting IA array */
	for (i = 1; i < N + 1; ++i)
		if (A->IA[i] == 0) 
			A->IA[i] = A->IA[i-1];

	/* Diagonal */
	if (M == N) /* square matrix */
	{
		for (i = 0; i < A->n; ++i)
		{
			int k1 = A->IA[i];
			int k2 = A->IA[i+1]-1;
			for (j = k1; j <= k2; ++j)
				if (i == A->JA[j]) 
					A->D[i] = A->AA[j];
		}
	}
}


/*----------------------------------------------------------------------------
 * Read upper triangular matrix from file in MM format to a CSR structure
 *--------------------------------------------------------------------------*/
void MATRIX_readCSR_SymmUpper (MAT* A, FILE* f)
{
	int M, N, nz, N_upper;
	int i, j, k, I, J, elem = 0;
	double VAL;
	char line[1025], header[1025];
	char *symm;
	ARRAY* a;
	
	fgets(header,200,f);	
	strtok(header," \n");
	strtok(NULL," \n");
	strtok(NULL," \n");
	strtok(NULL," \n");
	
	symm  = strtok(NULL," \n");
	
	if ((strcmp(symm,"symmetric") != 0 && strcmp(symm,"general") != 0))
		fprintf (stderr,"\n Error. Matrix format not supported. Exiting.. [MATRIX_readCSR]\n\n");
	
	do 
	{
		if (fgets(line,1025,f) == NULL) 
			exit(0);
	}while (line[0] == '%');   
	sscanf(line,"%d %d %d", &N, &M, &nz);
	
	if (strcmp(symm,"symmetric") == 0)
	{
		a = malloc (nz * sizeof(ARRAY));
		
		for (i = 0, k = 0; i < nz; ++i)
		{
			fscanf(f,"%d%d%lf",&I,&J,&VAL);
			
			if (I > J)
			{
				a[k].arr1 = VAL;
				a[k].arr2 = I - 1;
				a[k].arr3 = J - 1;
				k++;
			}
			else 
			{
				a[k].arr1 = VAL;
				a[k].arr2 = J - 1;
				a[k].arr3 = I - 1;
				k++;
			}
		}
		
		nz = k;
		a  = realloc (a,nz*sizeof(ARRAY));
		qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	}
	else if (strcmp(symm,"general") == 0)
	{
		N_upper = (nz + N)/2;
		a = malloc (N_upper * sizeof(ARRAY));
		
		for (i = 0, k = 0; i < nz; ++i)
		{
			fscanf(f, "%d%d%lf",&I, &J, &VAL);
			
			if (I <= J)
			{
				a[k].arr1 = VAL;
				a[k].arr2 = J - 1;
				a[k].arr3 = I - 1;
				k++;
			}
		}
		
		qsort(a, N_upper, sizeof(ARRAY), COMPARE_array);
		nz = N_upper;
	}

	/* reseve memory for matrices */
	A->m   = M;
	A->n   = N;
	A->nz  = nz;
	
	A->AA  = (double *) calloc(nz, sizeof (double));
	A->D   = (double *) calloc(N , sizeof (double));
	A->JA  = (int    *) calloc(nz, sizeof (int));
	A->IA  = (int    *) calloc(N+1,sizeof (int));
	
	for (i = 0; i < nz; ++i)
	{
		A->AA[i]   = a[i].arr1;
		A->JA[i]   = a[i].arr2;
		elem      += 1;
		A->IA[a[i].arr3+1] = elem;
	}

	free(a);
	
	/* Adjusting IA array */
	for (i = 1; i < N + 1; ++i)
		if (A->IA[i] == 0) 
			A->IA[i] = A->IA[i-1];

	/* Diagonal */
	if (M == N) /* square matrix */
	{
		for (i = 0; i < A->n; ++i)
		{
			int k1 = A->IA[i];
			int k2 = A->IA[i+1]-1;
			for (j = k1; j <= k2; ++j)
				if (i == A->JA[j]) 
					A->D[i] = A->AA[j];
		}
	}
}



/*----------------------------------------------------------------------------
 * Get the element ij stored as CSR
 *--------------------------------------------------------------------------*/
double MATRIX_aij (MAT* A, int i, int j)
{
	int k;
	int k1 = A->IA[i];
	int k2 = A->IA[i+1]-1;
	for (k = k1; k <= k2; ++k)
		if (A->JA[k] == j)
			return A->AA[k];
	return 0;
}

/*----------------------------------------------------------------------------
 * Print the matrix A in CSR format
 *--------------------------------------------------------------------------*/
void MATRIX_printCSR (MAT* A)
{
	int i, j, k1, k2;
	int n = A->n;
	for (i = 0; i < n; ++i)
	{
		k1 = A->IA[i];
		k2 = A->IA[i+1]-1;
		for (j = k1; j <= k2; j++)
			printf ("%d %d %le\n",A->JA[j],i,A->AA[j]); 
	}
}

/*----------------------------------------------------------------------------
 *  Print the matrix A in full format
 *--------------------------------------------------------------------------*/
void MATRIX_printFULL (MAT* A)
{
	int i, j;
	int n = A->n;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
			printf ("%6.2lf", MATRIX_aij(A,i,j));
		printf ("\n");
	}
}

/*----------------------------------------------------------------------------
 * Compute the matrix envelope/profile
 *--------------------------------------------------------------------------*/
long int MATRIX_envelope (MAT* A)
{
	int i;
	int n = A->n;
	unsigned long int envelope = 0;
	
	for (i = 0; i < n; ++i)
		envelope += i - A->JA[A->IA[i]];
	
	return envelope;
}

/*----------------------------------------------------------------------------
 * Compute the matrix bandwidth
 *--------------------------------------------------------------------------*/
long int MATRIX_bandwidth (MAT* A)
{
	int i, min;
	int n = A->n;
	unsigned long int bandwidth=0;
	
	for (i = 0; i < n; ++i)
	{
		min = (i - A->JA[A->IA[i]]);
		if (bandwidth < min) bandwidth = min; 
	}
	
	return bandwidth;
}



/*----------------------------------------------------------------------------
 * Compute the matrix wavefront
 *--------------------------------------------------------------------------*/
long int MATRIX_wavefront(MAT* A)
{
	
	return 0;
}

/*----------------------------------------------------------------------------
 * Compute the matrix-vector operation Ax = b
 *--------------------------------------------------------------------------*/
void MATRIX_matvec (MAT* A, double* x, double* b)
{
	int i, j, k1, k2;
	int n = A->n;
	double soma = 0;
	
	for (i = 0; i < n; ++i)
	{
		soma = 0;
		k1 = A->IA[i];
		k2 = A->IA[i+1]-1;
		for (j = k1; j <= k2; ++j)
			soma = soma + (A->AA[j]) * x[A->JA[j]];
		b[i] = soma;
	}  
}


/*----------------------------------------------------------------------------
 * Perform forward substituition to solve Ly = b
 *--------------------------------------------------------------------------*/
void MATRIX_forward (MAT* L, double* b, double* y)
{
	int i, j, k1, k2;
	int n = L->n;
	double sum;
	
	y[0] = b[0] / L->D[0];
  
	for (i = 1; i < n; ++i)
	{
		sum = 0.0;
		k1 = L->IA[i];
		k2 = L->IA[i+1]-1;
		for (j = k1; j <= k2; j++)
		{
			sum += L->AA[j]*y[L->JA[j]];
		}
		y[i] = (b[i] - sum) / L->D[i];
	}
}

/*----------------------------------------------------------------------------
 * Perform back substituition to solve Ux = y
 *--------------------------------------------------------------------------*/
void MATRIX_backward (MAT* U, double* y, double* x)
{
	int i, j, k1, k2;
	int n = U->n;
	double sum;
	
	x[n-1] = y[n-1] / U->D[n-1];
	
	for (i = n-2; i >= 0; --i)
	{
		sum = 0.0;
		k1 = U->IA[i];
		k2 = U->IA[i+1]-1;
		for (j = k1; j <= k2; j++)
		{
			sum += U->AA[j]*x[U->JA[j]];
		}
		x[i] = (y[i] - sum) / U->D[i];
	}
}

/*----------------------------------------------------------------------------
 * Perform the operation P*A*P' in CSR format
 *--------------------------------------------------------------------------*/
void MATRIX_permutation (MAT* A, int* p)
{

	int i, j, k;
	int n   = A->n;
	int m   = A->m;
	int nz  = A->nz;  

	MAT*  B     = malloc(      sizeof(MAT));
	      B->AA = malloc( nz  *sizeof(double));
	      B->JA = malloc( nz  *sizeof(int));
	      B->IA = malloc((n+1)*sizeof(int));
	      B->D  = malloc( n   *sizeof(double));
  
	B->IA[0] = 0;
	B->n     = n;
	B->m     = m;
	B->nz    = nz;

	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = A->IA[p[i]]; j <= A->IA[p[i]+1] - 1; ++j)
		{
			B->AA[k] = A->AA[j];
			B->JA[k] = A->JA[j];
			k  = k + 1;
		}
		B->IA[i+1] = k;    
	}
	
	ARRAY* a = malloc (nz*sizeof(ARRAY));
	int*   q = malloc (n *sizeof(int));
	
	for (i = 0; i < n; ++i) 
		q[p[i]] = i; 
	
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = B->IA[i]; j <= B->IA[i+1] - 1; ++j)
		{
			a[k].arr1 = B->AA[j];
			a[k].arr2 = q[B->JA[j]];
			a[k].arr3 = i;
			k = k + 1;
		}
		A->IA[i+1] = k;    
	}
	
	qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	
	for (i = 0; i < nz; ++i)
	{
		A->AA[i] = a[i].arr1;
		A->JA[i] = a[i].arr2;
	}
	
	free(a);
	free(q);
	MATRIX_clean(B);
}


/*----------------------------------------------------------------------------
 * Free up memory allocated for MAT struct
 *--------------------------------------------------------------------------*/
void MATRIX_clean (MAT* A)
{
	free(A->AA);
	free(A->JA);
	free(A->IA);
	if (A->m == A->n) free(A->D);
	free(A);  
}

