/*----------------------------------------------------------------------------
 * RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
// #include "../CommonFiles/protos.h"
#include "../CommonFiles/protos_parallel.h"

void mc60ad_(int* n, int* lirn, int* irn, int* icptr, int* icntl, int* iw, int* info);
void mc60cd_(int* n, int* nsup, int* lirn, int* irn, int* icptr, int* vars, int* jcntl, int* permsv, double* weight, int** pair, int* info, int* iw, double* w);

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
void REORDERING_RCM_HSL (MAT* A, int** perm, int root)
{
	int i;
	int n        = A->n; 
	int nsup     = n;
	int *irn     = A->JA;
	int *icptr   = A->IA;
	int lirn     = 2 * icptr[n-1];
	int vars[n];
	int jcntl[2] = { RCM, AUTOMATIC_PERIPHERAL };
	double weight[2];
	int pair[2][nsup/2];
	int info[4];
	int iw[3*nsup + 1];
	int iiw[n];
	double w[nsup];
	int icntl[2] = { 0, 6 };
	
	*perm = calloc(nsup, sizeof(int));
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < lirn; i++) ++irn[i];
	
	for (i = 0; i < n; i++) 
	{
		++icptr[i];
		vars[i] = 1;
	}
	
	mc60ad_(&n, &lirn, irn, icptr, icntl, iiw, info);
	
	mc60cd_(&n, &nsup, &lirn, irn, icptr, vars, jcntl, *perm, weight, (int**) pair, info, iw, w);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < lirn; i++) --irn[i];
	
	for (i = 0; i < n; i++) 
	{
		--icptr[i];
		--((*perm)[i]);
	}

}


