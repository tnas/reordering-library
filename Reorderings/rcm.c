#include "rcm.h"

/*----------------------------------------------------------------------------
 * RCM reordering from the LEVEL STRUCTURE in PSEUDO-PERIPHERAL algorithm
 *--------------------------------------------------------------------------*/
void REORDERING_RCM_opt (MAT* A, int** Fp, int s)
{ 
	int i;
	int n = A->n;
	
	int* q = calloc (n,sizeof(int));
	int* p = calloc (n,sizeof(int));
	
	q = GRAPH_bfs_RCM (A,s,q);
	
	/* Reverse order */
	for (i = 0; i < n; ++i)
		p[n-1-i] = q[i]; 
	
	(*Fp) = p;
	
	free(q);
}

/*----------------------------------------------------------------------------
 * Original RCM reordering
 *--------------------------------------------------------------------------*/
void REORDERING_RCM (MAT* A, int** Fp)
{ 
	int i, j, k = 0, x, s, e, nadj;
	int n = A->n;
	
	int* g      = GRAPH_LS_peripheral (A,&s,&e);
	int* p      = malloc (n*sizeof(int));
	int* q      = malloc (n*sizeof(int));
	int* close  = malloc (n*sizeof(int));
	
	GRAPH* G;
	int* adj;
	
	LIST* queue = NULL; 
	      
	for (i = 0; i < n; ++i)
		close[i] = -1;
	
	/* Insert in queue all the adj(x) in ascending order of degree */
	adj  = GRAPH_adjacent (A,s);
	nadj = GRAPH_degree (A,s);
	G    = malloc (nadj*sizeof(GRAPH));
	
	for (i = 0; i < nadj; ++i)
	{
		G[i].label  = adj[i];
		G[i].degree = GRAPH_degree (A,adj[i]);
	}
	
	qsort(G,nadj,sizeof(GRAPH),COMPARE_degr_ASC);
	
	for (i = 0; i < nadj; ++i)
		queue = LIST_insert_IF_NOT_EXIST (queue,G[i].label);
	
	free(G);
	free(adj);
		
	 /* Close (x) */
	close[s] = k; ++k;
	
	for (j = 1; j < n; ++j)
	{
		/* Get x from queue */
		x = LIST_first (queue);
		
		adj  = GRAPH_adjacent (A,x);
		nadj = GRAPH_degree (A,x);		
		G    = malloc (nadj*sizeof(GRAPH));
	
		for (i = 0; i < nadj; ++i)
		{
			G[i].label  = adj[i];
			G[i].degree = GRAPH_degree (A,adj[i]);
		}
		
		/* Insert in queue all the adj(x) in ascending order of degree */
		qsort(G,nadj,sizeof(GRAPH),COMPARE_degr_ASC);
		
		for (i = 0; i < nadj ; ++i)
		{
			if (close[G[i].label] == -1)
				queue = LIST_insert_IF_NOT_EXIST (queue,G[i].label);  
		}
		free(G);
		free(adj);
		
		close[x] = k; ++k;
		
		queue = LIST_remove (queue,x);

	}
		
	for (i = 0; i < n; ++i)
		q[close[i]] = i;
	
	/* Reverse order */
	for (i = 0, j = n - 1; i < n; ++i, --j)
		p[j] = q[i]; 
	
	(*Fp) = p;
	
	free(g);
	free(q);
	free(close);
}


/*----------------------------------------------------------------------------
 * HSL_MC60 RCM Reordering
 *--------------------------------------------------------------------------*/
long int REORDERING_HSL_RCM (MAT* A)
{
	int i, nsup, pair_lenght;
	int n        = A->n; 
	int *irn     = A->JA;
	int *icptr   = A->IA;
	int lirn     = A->nz;
	int jcntl[2] = { RCM, AUTOMATIC_PERIPHERAL };
	double weight[2];
	int info[4];
	double rinfo[4];
	int* vars;
	int** pair;
	int* iw;
	double* w;
	int* permsv;
	int* svar;
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	
	for (i = 0; i < n; i++) 
	{
		++irn[i];
		++icptr[i];
	}
	
	for (i = n; i < lirn; i++) ++irn[i];
	
	svar   = calloc(n, sizeof(int));
	vars   = calloc(n, sizeof(int));
	iw     = calloc(2*n + 2, sizeof(int));
	
	// To find supervariables and compress pattern
	mc60bd_(&n, &lirn, irn, icptr, &nsup, svar, vars, iw);
	
	free(iw);
	pair_lenght = nsup/2;
	permsv = calloc(nsup, sizeof(int));
	iw     = calloc(3*nsup + 1, sizeof(int));
	w      = calloc(nsup, sizeof(double));
	pair   = (int**) calloc(pair_lenght, sizeof(int*));
	
	for (i = 0; i < pair_lenght; ++i)
		pair[i] = (int*) calloc(2, sizeof(int));
	
	// To find supervariable permutation
	mc60cd_(&n, &nsup, &lirn, irn, icptr, vars, jcntl, permsv, weight, pair, info, iw, w);
	
	free(iw);
	iw = calloc(2*nsup + 1, sizeof(int));
	
	// To compute the profile and wavefront for a supervariable permutation
	mc60fd_(&n, &nsup, &lirn, irn, icptr, vars, permsv, iw, rinfo);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < n; i++) 
	{
		--irn[i];
		--icptr[i];
	}
	
	for (i = n; i < lirn; i++) --irn[i];
	
	free(vars);
	free(iw);
	free(w);
	free(permsv);
	free(svar);
	free(pair);
	
	return rinfo[SEMI_BANDWIDTH];
}
