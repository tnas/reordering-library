/*----------------------------------------------------------------------------
 * RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

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
	
	printf("The vector permutation is: ");
	for (i = 0; i < n; ++i) printf("%d ", p[i]);
	printf("\n");fflush(stdout);
	
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
void REORDERING_HSL_RCM (MAT* A, int** p)
{
	int i, j;
	int n        = A->n; 
	int nsup     = n;
	int *irn     = A->JA;
	int *icptr   = A->IA;
	int lirn     = A->nz;
	int jcntl[2] = { RCM, AUTOMATIC_PERIPHERAL };
	double weight[2];
	int info[4];
	double rinfo[4];
	int pair_lenght;
	int* vars;
	int* pair;
	int* iw;
	double* w;
	int* permsv;
	int* svar;
	int* possv;
	int* perm;
// 	int* perminv;
	
	pair_lenght = nsup/2;
	vars   = calloc(n, sizeof(int));
	pair   = calloc(2 * pair_lenght, sizeof(int));
	iw     = calloc(3*n + 1, sizeof(int));
	w      = calloc(n, sizeof(double));
	permsv = calloc(nsup, sizeof(int));
	svar   = calloc(n, sizeof(int));
	possv  = calloc(n, sizeof(int));
	perm   = calloc(n, sizeof(int));
	*p     = calloc(nsup, sizeof(int));
// 	perminv = calloc(n, sizeof(int));
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	
	for (i = 0; i < n; i++) 
	{
		++irn[i];
		++icptr[i];
	}
	
	for (i = n; i < lirn; i++) ++irn[i];
	
	mc60bd_(&n, &lirn, irn, icptr, &nsup, svar, vars, iw);
	
	mc60cd_(&n, &nsup, &lirn, irn, icptr, vars, jcntl, permsv, weight, pair, info, iw, w);
	
	mc60fd_(&n, &nsup, &lirn, irn, icptr, vars, permsv, iw, rinfo);
	
	mc60dd_(&n, &nsup, svar, vars, permsv, perm, possv);

	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < n; i++) 
	{
		--irn[i];
		--icptr[i];
// 		(*p)[--(perm[i])] = i;
		(*p)[i] = --(perm[i]);
// 		perminv[i] = --(perm[i]);
	}
	
	for (i = n; i < lirn; i++) --irn[i];
	
// 	/* Reverse order */
// 	for (i = 0, j = n - 1; i < n; ++i, --j)
// 		(*p)[j] = perminv[i];
	
	printf("The bandwidth is %f\n", rinfo[2]);fflush(stdout);
	
	printf("The chosen permutation is: ");
	for (i = 0; i < n; ++i) printf("%d ", (*p)[i]);
	printf("\n");fflush(stdout);
	
	free(vars);
	free(iw);
	free(w);
	free(permsv);
	free(svar);
	free(possv);
	free(perm);
	free(pair);
// 	free(perminv);
}
