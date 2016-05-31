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
// 	int node_e, node_s;
	int *adj1,*adj2;
	
	int W1 = 1;
	int W2 = 2;
	
	LIST* L = NULL;
	LIST* q;
	
// 	int* g = GRAPH_LS_peripheral (A,&node_s,&node_e);
// 	free(g);
	
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