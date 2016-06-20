/*----------------------------------------------------------------------------
 * RCM REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos.h"

/*----------------------------------------------------------------------------
 * SLOAN reordering
 *--------------------------------------------------------------------------*/
void REORDERING_SLOAN (MAT* A, int** Fp, int node_s, int node_e)
{
	int i,j,k,I,J,K,n = A->n,max_P,nadj1,nadj2;
	int *adj1,*adj2;
	
	int W1 = 1;
	int W2 = 2;
	
	LIST* L = NULL;
	LIST* q;
	
	int* d = calloc (n,sizeof (int));
	int* p = calloc (n,sizeof (int));
	int* P = calloc (n,sizeof (int));
	int* status = calloc (n,sizeof (int));
	
	/* start vertex s and end vertex e */
	GRAPH_bfs (A,node_e,d);
	
	/* INACTIVE status and priority P(i) = W1*d(i,e) - W2*(degree(i)+1) */
	/* INACTIVE  = -2 | PREACTIVE  = -1 | ACTIVE  =  0 \ POSACTIVE  = 1 */
	
	for (i = 0; i < n; ++i)
	{
		status[i] = -2;
		P[i]      = W1*d[i] - W2*(GRAPH_degree(A,i) + 1);
	}
	
	L = LIST_insert_IF_NOT_EXIST (L,node_s);
	status[node_s] = -1;

	for (I = 0; I < n; ++I)
	{
		max_P = -9999999;
		for (q = L; q != NULL; q = q->next)
		{
			if (max_P < P[q->data])
			{
				max_P = P[q->data];
				i = q->data;
			}
		}

		L = LIST_remove (L,i);
		
		if (status[i] == -1) // PREACTIVE
		{
			adj1  = GRAPH_adjacent(A,i);
			nadj1 = GRAPH_degree  (A,i);
			for (J = 0; J < nadj1; ++J)
			{
				j = adj1[J];
				P[j] += W2;
				if (status[j] == -2) // INACTIVE
				{
					L = LIST_insert_IF_NOT_EXIST (L,j);
					status[j] = -1;
				}
			}
			free(adj1);
		}
		
		p[I] = i;	
		status[i] = 1;
	
		adj1  = GRAPH_adjacent(A,i);
		nadj1 = GRAPH_degree  (A,i);
		for (J = 0; J < nadj1; ++J)
		{
			j = adj1[J];
			if (status[j] == -1) // PREACTIVE
			{
				P[j] += W2;
				status[j] = 0; // ACTIVE
				adj2  = GRAPH_adjacent(A,j);
				nadj2 = GRAPH_degree  (A,j);
				for (K = 0; K < nadj2; ++K)
				{
					k = adj2[K];
					if (status[k] == -1 || status[k] == 0)
						P[k] += W2;
					if (status[k] == -2)
					{
						P[k] += W2;
						L = LIST_insert_IF_NOT_EXIST (L,k);
						status[k] = -1;
					}
				}
				free(adj2);
			}
		}
		free(adj1);
	}
	free(d);
	free(P);
	free(status);
	(*Fp) = p;
}

/*----------------------------------------------------------------------------
 * HSL Sloan Reordering
 * 
 * @author: Thiago Nascimento - nascimenthiago@gmail.com
 * @since: 15-06-2016
 *--------------------------------------------------------------------------*/
void REORDERING_SLOAN_HSL (MAT* A, int** p, int start_node, int end_node)
{
	int i;
	int n        = A->n; 
	int nsup     = n;
	int *irn     = A->JA;
	int *icptr   = A->IA;
	int lirn     = A->nz;
// 	int jcntl[2] = { SLOAN, AUTOMATIC_PERIPHERAL };
	int jcntl[2] = { SLOAN, ESPECIFIED_PERIPHERAL };
	double weight[2] = { SLOAN_W1, SLOAN_W2 };
	int info[4];
	int pair_lenght;
	int* vars;
	int** pair;
	int* iw;
	double* w;
	int* permsv;
	int* svar;
	int* possv;
	int* perm;
	double rinfo[4];
	
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
	
	for (i = 0; i < pair_lenght; i++)
	{
		pair[0][i] = end_node;
		pair[1][i] = start_node;
	}
	
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
		(*p)[--(perm[i])] = i;
	}
	
	for (i = n; i < lirn; i++) --irn[i];
	
	free(vars);
	free(iw);
	free(w);
	free(permsv);
	free(svar);
	free(possv);
	free(perm);
	free(pair);
}
