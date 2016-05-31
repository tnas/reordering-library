/*----------------------------------------------------------------------------
 * Sloan REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN reordering
 *--------------------------------------------------------------------------*/

void update_far_neighbors() 
{
}

void Parallel_Sloan (MAT* A, int** Fp, int node_s, int node_e)
{
	int vertex,neighbor,k,next_id,ngb,K,n = A->n;
	int max_P,vertex_degree,vertex_degree2;
	int *neighbors,*adjacents2;
	int num_threads;
	int W1 = 1;
	int W2 = 2;
	
	LIST* queue = NULL;
	LIST* q;
	omp_lock_t lock_queue, lock_priority;
	
	int* distance = calloc (n,sizeof (int));
	int* p        = calloc (n,sizeof (int));
	int* priority = calloc (n,sizeof (int));
	int* status   = calloc (n,sizeof (int));
	
	/* start vertex s and end vertex e */
	GRAPH_bfs(A, node_e, distance);
	
	/* INACTIVE status and priority P(i) = W1*d(i,e) - W2*(degree(i)+1) */
		
	#pragma omp parallel 
	{
		#pragma omp single
		num_threads = omp_get_num_threads();
		
		#pragma omp for private(vertex)
		for (vertex = 0; vertex < n; ++vertex)
		{
			status[vertex]   = INACTIVE;
			priority[vertex] = W1*distance[vertex] - W2*(GRAPH_degree(A, vertex) + 1);
		}
		
		#pragma omp single nowait
		queue = LIST_insert_IF_NOT_EXIST(queue, node_s);
		
		#pragma omp single nowait
		status[node_s] = PREACTIVE;
		
		#pragma omp single nowait
		omp_init_lock(&lock_queue);
		
		#pragma omp single nowait
		omp_init_lock(&lock_priority);
	}
	
	
	#pragma omp parallel num_threads(num_threads) private(max_P, q, J, K, j, k, vertex, adjacents1, vertex_degree1, adjacents2, vertex_degree2)
	{
		#pragma omp single
		next_id = 0;
		
		while (queue != NULL || next_id < n)
		{
			max_P = -9999999;
			vertex = UNDEF_NODE;
			
			omp_set_lock(&lock_queue);
			
			// Chosing vertex that maximizes the priority
			for (q = queue; q != NULL; q = q->next)
			{
				if (max_P < priority[q->data])
				{
					max_P = priority[q->data];
					vertex = q->data;
				}
			}
			
			if (vertex != UNDEF_NODE) 
				queue = LIST_remove (queue, vertex);
			
			omp_unset_lock(&lock_queue);
			
			if (vertex != UNDEF_NODE) 
			{
				neighbors     = GRAPH_adjacent(A, vertex);
				vertex_degree = GRAPH_degree  (A, vertex);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor = neighbors[ngb];
					
					if ((status[vertex] == PREACTIVE) &&
						(status[neighbor] == INACTIVE || status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += W2; 
						status[neighbor] = ACTIVE;
						update_far_neighbors();
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == ACTIVE))
					{
						priority[neighbor] += W2;
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += W2;
						status[neighbor] = ACTIVE;
						update_far_neighbors();
					}
				}
				
				#pragma omp critical
				{
					p[next_id++] = vertex;
					status[vertex] = NUMBERED;
				}
				
				free(neighbors);
			}
		}
	}
	
	#pragma omp parallel
	{
		#pragma omp single nowait
		omp_destroy_lock(&lock_queue);
		
		#pragma omp single nowait
		omp_destroy_lock(&lock_priority);
		
		#pragma omp single nowait
		free(distance);
		
		#pragma omp single nowait
		free(priority);
		
		#pragma omp single nowait
		free(status);
		
		#pragma omp single nowait
		(*Fp) = p;
	}
	
// 	printf("Permutation vector of size %d: ", n);fflush(stdout);
// 	int i;
// 	for (i = 0; i < n; ++i)
// 		printf("%d ", (*Fp)[i]);
// 	printf("\n");fflush(stdout);
}
