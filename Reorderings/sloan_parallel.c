/*----------------------------------------------------------------------------
 * Sloan REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN reordering
 *--------------------------------------------------------------------------*/
void Parallel_Sloan (MAT* A, int** Fp, int node_s, int node_e)
{
	int vertex,j,k,I,J,K,n = A->n,max_P,vertex_degree1,vertex_degree2;
	int *adjacents1,*adjacents2;
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
	GRAPH_bfs (A, node_e, distance);
	
	/* INACTIVE status and priority P(i) = W1*d(i,e) - W2*(degree(i)+1) */
	/* INACTIVE  = -2 | PREACTIVE  = -1 | ACTIVE  =  0 \ POSACTIVE  = 1 */
	
	#pragma omp parallel 
	{
		#pragma omp single
		num_threads = omp_get_num_threads();
		
		#pragma omp for private(vertex)
		for (vertex = 0; vertex < n; ++vertex)
		{
			status[vertex]   = -2;
			priority[vertex] = W1*distance[vertex] - W2*(GRAPH_degree(A, vertex) + 1);
		}
		
		#pragma omp single nowait
		queue = LIST_insert_IF_NOT_EXIST (queue, node_s);
		
		#pragma omp single nowait
		status[node_s] = -1;
		
		#pragma omp single nowait
		omp_init_lock(&lock_queue);
		
		#pragma omp single nowait
		omp_init_lock(&lock_priority);
	}
	
	
	#pragma omp parallel num_threads(num_threads) private(max_P, q, J, K, j, k, vertex, adjacents1, vertex_degree1, adjacents2, vertex_degree2)
	{
		#pragma omp single
		I = 0;
		
		while (queue != NULL || I < n)
		{
			max_P = -9999999;
			vertex = -1;
			
			omp_set_lock(&lock_queue);
			
			for (q = queue; q != NULL; q = q->next)
			{
				if (max_P < priority[q->data])
				{
					max_P = priority[q->data];
					vertex = q->data;
				}
			}
			
			if (vertex != -1) 
				queue = LIST_remove (queue, vertex);
			
			omp_unset_lock(&lock_queue);
			
			if (vertex != -1) 
			{
// 				printf("Thread %d get vertex %d\n", omp_get_thread_num(), vertex);fflush(stdout);
				
				if (status[vertex] == -1)
				{
					adjacents1     = GRAPH_adjacent(A, vertex);
					vertex_degree1 = GRAPH_degree  (A, vertex);
					
					for (J = 0; J < vertex_degree1; ++J)
					{
						j = adjacents1[J];
						
						omp_set_lock(&lock_priority);
						priority[j] += W2;
						omp_unset_lock(&lock_priority);
						
						omp_set_lock(&lock_queue);
						if (status[j] == -2)
						{
							queue = LIST_insert_IF_NOT_EXIST (queue, j);
							status[j] = -1;
						}
						omp_unset_lock(&lock_queue);
					}
					
					free(adjacents1);
				}
				
				#pragma omp critical
				{
// 					printf("Thread %d setting vertex %d at position %d.\n", omp_get_thread_num(), vertex, I);
					p[I++] = vertex;	
				}
				
				status[vertex] = 1;
			
				adjacents1     = GRAPH_adjacent(A,vertex);
				vertex_degree1 = GRAPH_degree  (A,vertex);
				
				for (J = 0; J < vertex_degree1; ++J)
				{
					j = adjacents1[J];
					
					if (status[j] == -1)
					{
						omp_set_lock(&lock_priority);
						priority[j] += W2;
						status[j] = 0;
						omp_unset_lock(&lock_priority);
						
						adjacents2  = GRAPH_adjacent(A,j);
						vertex_degree2 = GRAPH_degree  (A,j);
						
						for (K = 0; K < vertex_degree2; ++K)
						{
							k = adjacents2[K];
							
							if (status[k] == -1 || status[k] == 0)
							{
								omp_set_lock(&lock_priority);
								priority[k] += W2;
								omp_unset_lock(&lock_priority);
							}
							
							omp_set_lock(&lock_queue);
							if (status[k] == -2)
							{
								omp_set_lock(&lock_priority);
								priority[k] += W2;
								omp_unset_lock(&lock_priority);
								
								
								queue = LIST_insert_IF_NOT_EXIST (queue,k);
								status[k] = -1;
								
							}
							omp_unset_lock(&lock_queue);
						}
						
						free(adjacents2);
					}
				}
				
				free(adjacents1);
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