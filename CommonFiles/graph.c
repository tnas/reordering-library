#include "graph.h"

/*----------------------------------------------------------------------------
 * GRAPH FUNCTIONS
 *--------------------------------------------------------------------------*/

int COMPARE_degr_ASC (const void * a, const void * b)
{ 
	return ((GRAPH*)a)->degree - ((GRAPH*)b)->degree;
}

int COMPARE_dist_degr_DES (const void * a, const void * b)
{ 
	if (((GRAPH*)a)->distance > ((GRAPH*)b)->distance) return -1;
	if (((GRAPH*)a)->distance < ((GRAPH*)b)->distance) return 1;
	return ((GRAPH*)a)->degree - ((GRAPH*)b)->degree;
}

int COMPARE_dist_degr_ASC (const void * a, const void * b)
{ 
	if (((GRAPH*)a)->distance > ((GRAPH*)b)->distance) return 1;
	if (((GRAPH*)a)->distance < ((GRAPH*)b)->distance) return -1;
	return ((GRAPH*)a)->degree - ((GRAPH*)b)->degree;
}

/*----------------------------------------------------------------------------
 * Return the degree of the vertex x in graph (i.e. matrix) A
 *--------------------------------------------------------------------------*/
int GRAPH_degree (MAT* A, int x)
{
	if (A->n == 0)
	{
		printf("error: GRAPH does not exist. Exiting.. [GRAPH_degree]\n");
		exit(0);
	}
	if (x > A->n - 1 || x < 0)
	{
		printf("error: Vertex %d does not exist in this GRAPH. Exiting.. [GRAPH_degree]\n", x);
		exit(0);
	}
	
	return A->IA[x+1] - A->IA[x]  /*Diagonal*/ - 1;
}


/*----------------------------------------------------------------------------
 * Return all vertices adjacents to x in graph (i.e. matrix) A
 *--------------------------------------------------------------------------*/
int* GRAPH_adjacent (MAT* A, int x)
{
	if (A->n == 0)
	{
		printf("error: GRAPH does not exist. Exiting.. [GRAPH_degree]\n");
		exit(0);
	}
	if (x > A->n - 1 || x < 0)
	{
		printf("error: Vertex %d does not exist in this GRAPH. Exiting.. [GRAPH_adjacent]\n", x);
		exit(0);
	}
	
	int i, k = 0;
	
	int  n   = GRAPH_degree (A,x);
	int *adj = malloc(n*sizeof(int));

	for (i = A->IA[x]; i < A->IA[x+1]; ++i)
	{
		if (A->JA[i] != x)
		{
			adj[k] = A->JA[i];
			k++;
		}
	}
	
	return adj;
}


/*-------------------------------------------------------------------------------------
 * Return the degree of the vertex x in graph (i.e. matrix) A at level adjacency_level
 *-------------------------------------------------------------------------------------*/
int GRAPH_degree_per_level (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors)
{
	if (A->n == 0)
	{
		printf("error: GRAPH does not exist. Exiting.. [GRAPH_degree_per_level]\n");
		exit(0);
	}
	if (x > A->n - 1 || x < 0)
	{
		printf("error: Vertex %d does not exist in this GRAPH. Exiting.. [GRAPH_degree_per_level]\n", x);
		exit(0);
	}
	
	int i, k = 0;
	
	for (i = A->IA[x]; i < A->IA[x+1]; ++i)
	{
		if (A->JA[i] != x && levels[A->JA[i]] == adjacency_level && colors[A->JA[i]] == UNREACHED) ++k;
	}
	
	return k;
}


/*-----------------------------------------------------------------------------------
 * Return all vertices adjacents x in graph (i.e. matrix) A at level adjacency_level
 *-----------------------------------------------------------------------------------*/
GRAPH* GRAPH_adjacent_per_level (MAT* A, int x, const int* levels, const int adjacency_level, const int* colors)
{
	if (A->n == 0)
	{
		printf("error: GRAPH does not exist. Exiting.. [GRAPH_degree]\n");
		exit(0);
	}
	if (x > A->n - 1 || x < 0)
	{
		printf("error: Vertex %d does not exist in this GRAPH. Exiting.. [GRAPH_adjacent]\n", x);
		exit(0);
	}
	
	int i, k = 0;
	
	int  degree = GRAPH_degree_per_level(A, x, levels, adjacency_level, colors);
	GRAPH* adj  = malloc(degree * sizeof(GRAPH));

	for (i = A->IA[x]; i < A->IA[x+1]; ++i)
	{
		if (A->JA[i] != x && levels[A->JA[i]] == adjacency_level && colors[A->JA[i]] == UNREACHED)
		{
			adj[k].label = A->JA[i];
			adj[k].degree = GRAPH_degree(A, A->JA[i]);
			k++;
		}
	}
	
	return adj;
}




/*----------------------------------------------------------------------------
 * Perform the Breadth-First Search in the graph (i.e. matrix) A
 *--------------------------------------------------------------------------*/
void GRAPH_bfs (MAT* A, int x, int* dist)
{
	int i, u, v, dgr; 
	int n = A->n;

	int* color = (int*) calloc (n,sizeof (int));
	int* adj;

	LIST *L = NULL;

	/* not visited  = -1  */
	/* visited      =  0  */
	/* not reached  =  1  */

	for (i = 0; i < n; ++i)
	{
		color[i] = -1;
		dist[i]  = 0;
	}

	color[x] = 0;
	L = LIST_insert_IF_NOT_EXIST (L,x);

	while (L != NULL)
	{
		u   = LIST_first (L);
		adj = GRAPH_adjacent (A,u);
		dgr = GRAPH_degree (A,u);

		for (i = 0; i < dgr; ++i)
		{
			v = adj[i];
			if (color[v] == -1)
			{
				color[v] = 0;
				dist[v] = dist[u] + 1; 
				L = LIST_insert_IF_NOT_EXIST (L,v);
			}
		}

		L = LIST_remove (L,u);
		color[u] = 1;
		free(adj); 
	}

	for (i = 0; i < n; ++i)
	{
		if (color[i] == -1)
		{
			printf ("error: Nodes unreacheable. Exiting.. [GRAPH_bfs]\n");
			exit(0);
		}
	}

	free(L);
	free(color);
}


/*----------------------------------------------------------------------------
 * Perform a modified Breadth-First Search for the RCM reordering
 *--------------------------------------------------------------------------*/
int* GRAPH_bfs_RCM (MAT* A, int x, int* p)
{
	int i, u, v, k = 0, dgr; 
	int n = A->n;

	int* color = (int*) calloc (n,sizeof (int));
	int* dist  = (int*) calloc (n,sizeof (int));
	int* adj;

	LIST  *L = NULL;
	GRAPH *G;

	/* not visited  = -1  */
	/* visited      =  0  */
	/* not reached  =  1  */

	for (i = 0; i < n; ++i)
	{
		color[i] = -1;
		p[i]     = 0;
	}

	color[x] = 0;
	L = LIST_insert_IF_NOT_EXIST (L,x);

	while (L != NULL)
	{
		u    = LIST_first (L);
		p[k] = u; ++k;
		adj  = GRAPH_adjacent (A,u);
		dgr  = GRAPH_degree (A,u);
		G    = malloc (dgr*sizeof(GRAPH));
		
		for (i = 0; i < dgr; ++i)
		{
			G[i].label  = adj[i];
			G[i].degree = GRAPH_degree (A,adj[i]);
		}
		
		qsort(G,dgr,sizeof(GRAPH),COMPARE_degr_ASC);
		
		for (i = 0; i < dgr; ++i)
		{
			v = G[i].label;
			if (color[v] == -1)
			{
				color[v] = 0;
				dist[v] = dist[u] + 1; 
				L = LIST_insert_IF_NOT_EXIST (L,v);
			}
		}

		L = LIST_remove (L,u);
		color[u] = 1;
		free(adj);
		free(G);
	}

	for (i = 0; i < n; ++i)
	{
		if (color[i] == -1)
		{
			printf ("error: Nodes unreacheable. Exiting.. [GRAPH_bfs]\n");
			exit(0);
		}
	}
	
	
	free(L);
	free(color);
	free(dist);
	return p;
}

/*----------------------------------------------------------------------------
 * Find the depth of an LEVEL STRUCTURE
 * (i.e. the total number of levels - same as find the maximum of an array)
 *--------------------------------------------------------------------------*/
int GRAPH_LS_depth (int* LS, int n)
{
	if (LS == NULL || n <= 0)
	{
		printf ("error: Array is empty. Exiting.. [GRAPH_LS_depth]\n");
		exit(0);
	}

	int i, depth = 0;
	
	for (i = 0; i < n; ++i)
		if (depth < LS[i])
			depth = LS[i];

	return depth + 1;
}

/*----------------------------------------------------------------------------
 * Find the width of an LEVEL STRUCTURE
 * (i.e. maximum number of nodes which belong to a single level)
 *--------------------------------------------------------------------------*/
int GRAPH_LS_width (int* LS, int n)
{
	if (LS == NULL || n <= 0)
	{
		printf ("error: Array is empty. Exiting.. [GRAPH_LS_width]\n");
		exit(0);
	}
	
	int i, width = 0;	
	int size_L = GRAPH_LS_depth (LS,n);
	int* c     = calloc (size_L,sizeof(int));

        for (i = 0; i < n; ++i)
		c[LS[i]]++;
	
	for (i = 0; i < size_L; ++i)
		if (width < c[i])
			width = c[i];
	free(c);
	return width;
}

/*----------------------------------------------------------------------------
 * Return the Last level shrinked from a LEVEL STRUCTURE
 *--------------------------------------------------------------------------*/
LIST* GRAPH_LS_last_level (MAT* A, int* LS, int n)
{
	int i,k1,k2;
	int last = GRAPH_LS_depth (LS,n) - 1;
	
	GRAPH *G = (GRAPH*) malloc (n*sizeof(GRAPH));
	LIST  *L = NULL;
	
	for (i = 0; i < n; ++i)
	{
		G[i].label    = i;
		G[i].distance = LS[i];
		G[i].degree   = GRAPH_degree (A,i);			
	}

	qsort (G,n,sizeof(GRAPH),COMPARE_dist_degr_DES);

	k1 = k2 = 0;
	for (i = 0; i < n; ++i)
	{
		k1 = G[i].degree; 
		if (G[i].distance != last) 
			break;
		if (k1 != k2)
		{ 
			L  = LIST_insert_IF_NOT_EXIST (L,G[i].label);
			k2 = k1; 
		}
	}
	free(G);
	return L;
}

/*----------------------------------------------------------------------------
 * Find the pseudo-peripheral nodes in the graph (i.e. matrix) A 
 *--------------------------------------------------------------------------*/
int* GRAPH_LS_peripheral (MAT* A, int *node_s, int* node_e)
{
	int i, s, e, x, dgr, min_dgr, width;
	int n = width = A->n;
	
	int*  LS1 =  (int*) calloc (n,sizeof (int));;
	int*  LS2 =  (int*) calloc (n,sizeof (int));;
	LIST* L = NULL;
	
	/* Choose a vertex s with minimum degree */
	min_dgr = n;
	for (i = 0; i < n; ++i)
	{
		dgr = GRAPH_degree (A,i);
		if (min_dgr > dgr)
		{
			min_dgr = dgr;
			s = i;
		}
	}
	
	/* Construct the Level Strucure of s */
	GRAPH_bfs (A,s,LS1);
	L    = GRAPH_LS_last_level (A,LS1,n);
	
	while (L!=NULL)
	{
		x   = LIST_first (L);
		L   = LIST_remove (L,x);
		GRAPH_bfs (A,x,LS2);

		if (GRAPH_LS_depth (LS2,n) > GRAPH_LS_depth (LS1,n) 
		 && GRAPH_LS_width (LS2,n) < width)
		{
			s     = x; 
			GRAPH_bfs (A,s,LS1);
			L     = GRAPH_LS_last_level (A,LS1,n);
			width = n;
		}
		else if (GRAPH_LS_width (LS2,n) < width)
		{
			e     = x; 
			width = GRAPH_LS_width (LS2,n);
		}
	}

	free(LS2);
    
	*node_s = s;
	*node_e = e;
	
	return LS1;	
}