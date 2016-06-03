/*----------------------------------------------------------------------------
 * Sloan REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN reordering
 *--------------------------------------------------------------------------*/

void inline update_far_neighbors(MAT* A, int** status, int** priority, int node, LIST** worklist) 
{
	int node_degree, ngb, neighbor;
	int* neighbors;
	
	neighbors   = GRAPH_adjacent(A, node);
	node_degree = GRAPH_degree(A, node);
	
// 	int i;
// 	printf("Processing FAR neighbors: ");fflush(stdout);
// 	for (i = 0; i< node_degree; i++) printf("%d ", neighbors[i]); fflush(stdout);
// 	printf("\n");fflush(stdout);
	
	for (ngb = 0; ngb < node_degree; ++ngb) 
	{
		neighbor = neighbors[ngb];
		
		if ((*status)[neighbor] == INACTIVE)
		{
			(*status)[neighbor] = PREACTIVE;
			*worklist = LIST_insert_IF_NOT_EXIST(*worklist, neighbor);
		}
		
		(*priority)[node] += SLOAN_W2;
	}
	
	free(neighbors);
}


void Parallel_Sloan (MAT* adjacency, int** Fp, int start_node, int end_node)
{
	int node, next_id, num_nodes;
	LIST* queue;
	omp_lock_t lock_queue;
	
	num_nodes        = adjacency->n;
	int* distance    = calloc (num_nodes,sizeof (int));
	int* permutation = calloc (num_nodes,sizeof (int));
	int* priority    = calloc (num_nodes,sizeof (int));
	int* status      = calloc (num_nodes,sizeof (int));
	
	GRAPH_bfs(adjacency, end_node, distance);
	queue = NULL;
	next_id = 0;
	
	int min = 999999;
	int max = -999999;
	
	#pragma omp parallel 
	{
		#pragma omp for private(node)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(GRAPH_degree(adjacency, node) + 1);
			if (priority[node] < min) min = priority[node];
			if (priority[node] > max) max = priority[node];
// 			printf("priority of node %d: %d\n", node, priority[node]);fflush(stdout);
		}
		printf("minimum priority: %d\n", min);fflush(stdout);
		printf("max priority: %d\n", max);fflush(stdout);
		min *= -1;
		
		#pragma omp for private(node)
		for (node = 0; node < num_nodes; ++node)
			priority[node] += min;
		
		#pragma omp single nowait
		queue = LIST_insert_IF_NOT_EXIST(queue, start_node);
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		omp_init_lock(&lock_queue);
	}
	
	
	#pragma omp parallel 
	{
		int vertex, max_priority, vertex_degree, neighbor, ngb;
		int* neighbors;
		LIST* q;
		
		while (queue != NULL || next_id < num_nodes)
		{
			max_priority = INIT_PRIORITY;
			vertex       = UNDEF_NODE;
			
// 			printf("worklist: ");LIST_print(queue);fflush(stdout);
			
			omp_set_lock(&lock_queue);
			
			// Chosing vertex that maximizes the priority
			for (q = queue; q != NULL; q = q->next)
			{
				if (max_priority < priority[q->data])
				{
					max_priority = priority[q->data];
					vertex = q->data;
				}
			}
			
			if (vertex != UNDEF_NODE) 
				queue = LIST_remove(queue, vertex);
			
			omp_unset_lock(&lock_queue);
			
			if (vertex != UNDEF_NODE) 
			{
// 				printf("Processing vertex %d\n", vertex);fflush(stdout);
				
				neighbors     = GRAPH_adjacent(adjacency, vertex);
				vertex_degree = GRAPH_degree  (adjacency, vertex);
				
// 				int i;
// 				printf("Processing neighbors: ");fflush(stdout);
// 				for (i = 0; i< vertex_degree; i++) printf("%d ", neighbors[i]); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
// 				printf("Status of nodes: ");fflush(stdout);
// 				for (i = 0; i< num_nodes; i++) printf("%d ", status[i]); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor = neighbors[ngb];
					
					if ((status[vertex] == PREACTIVE) && (status[neighbor] == INACTIVE))
					{
						omp_set_lock(&lock_queue);
						if (status[neighbor] == INACTIVE)
							queue = LIST_insert_IF_NOT_EXIST(queue, neighbor);
						
						update_far_neighbors(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
						
						priority[neighbor] += SLOAN_W2; 
						status[neighbor] = ACTIVE;
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += SLOAN_W2; 
						status[neighbor] = ACTIVE;
						
						omp_set_lock(&lock_queue);
						update_far_neighbors(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == ACTIVE))
					{
						priority[neighbor] += SLOAN_W2;
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += SLOAN_W2;
						status[neighbor] = ACTIVE;
						
						omp_set_lock(&lock_queue);
						update_far_neighbors(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
					}
				}
				
				#pragma omp critical
				{
// 					printf("Setting node %d at position %d of permutation\n", vertex, next_id);fflush(stdout);
					permutation[next_id++] = vertex;
					status[vertex] = NUMBERED;
				}
				
				free(neighbors);
			}
		}
	}
	
	min = 999999;
	max = -999999;
	
	#pragma omp single
	for (node = 0; node < num_nodes; ++node)
	{
		if (priority[node] < min) min = priority[node];
		if (priority[node] > max) max = priority[node];
	}
	
	printf("minimum priority: %d\n", min);fflush(stdout);
	printf("max priority: %d\n", max);fflush(stdout);
	printf("start/end node: %d/%d\n", start_node, end_node);fflush(stdout);
	
	
	#pragma omp parallel
	{
		#pragma omp single nowait
		omp_destroy_lock(&lock_queue);
		
		#pragma omp single nowait
		free(distance);
		
		#pragma omp single nowait
		free(priority);
		
		#pragma omp single nowait
		free(status);
		
		#pragma omp single nowait
		*Fp = permutation;
	}
	
// 	printf("Permutation vector of size %d: ", n);fflush(stdout);
// 	int i;
// 	for (i = 0; i < n; ++i)
// 		printf("%d ", (*Fp)[i]);
// 	printf("\n");fflush(stdout);
}





void Parallel_Sloan_v1(MAT* adjacency, int** Fp, int start_node, int end_node)
{
	int node, next_id, num_nodes;
	LIST* queue;
	omp_lock_t lock_queue;
	
	num_nodes        = adjacency->n;
	int* distance    = calloc (num_nodes,sizeof (int));
	int* permutation = calloc (num_nodes,sizeof (int));
	int* priority    = calloc (num_nodes,sizeof (int));
	int* status      = calloc (num_nodes,sizeof (int));
	
	GRAPH_bfs(adjacency, end_node, distance);
	queue = NULL;
	next_id = 0;
	
	#pragma omp parallel 
	{
		#pragma omp for private(node)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(GRAPH_degree(adjacency, node) + 1);
		}
		
		#pragma omp single nowait
		queue = LIST_insert_IF_NOT_EXIST(queue, start_node);
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		omp_init_lock(&lock_queue);
	}
	
	
	#pragma omp parallel 
	{
		int vertex, max_priority, vertex_degree, neighbor, ngb;
		int* neighbors;
		LIST* q;
		
		while (queue != NULL || next_id < num_nodes)
		{
			max_priority = INIT_PRIORITY;
			vertex       = UNDEF_NODE;
			
// 			printf("worklist: ");LIST_print(queue);fflush(stdout);
			
			omp_set_lock(&lock_queue);
			
			// Chosing vertex that maximizes the priority
			for (q = queue; q != NULL; q = q->next)
			{
				if (max_priority < priority[q->data])
				{
					max_priority = priority[q->data];
					vertex = q->data;
				}
			}
			
			if (vertex != UNDEF_NODE) 
				queue = LIST_remove(queue, vertex);
			
			omp_unset_lock(&lock_queue);
			
			if (vertex != UNDEF_NODE) 
			{
// 				printf("Processing vertex %d\n", vertex);fflush(stdout);
				
				neighbors     = GRAPH_adjacent(adjacency, vertex);
				vertex_degree = GRAPH_degree  (adjacency, vertex);
				
// 				int i;
// 				printf("Processing neighbors: ");fflush(stdout);
// 				for (i = 0; i< vertex_degree; i++) printf("%d ", neighbors[i]); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
// 				printf("Status of nodes: ");fflush(stdout);
// 				for (i = 0; i< num_nodes; i++) printf("%d ", status[i]); fflush(stdout);
// 				printf("\n");fflush(stdout);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor = neighbors[ngb];
					
					if ((status[vertex] == PREACTIVE) && (status[neighbor] == INACTIVE))
					{
						omp_set_lock(&lock_queue);
						if (status[neighbor] == INACTIVE)
							queue = LIST_insert_IF_NOT_EXIST(queue, neighbor);
						
						update_far_neighbors(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
						
						priority[neighbor] += SLOAN_W2; 
						status[neighbor] = ACTIVE;
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += SLOAN_W2; 
						status[neighbor] = ACTIVE;
						
						omp_set_lock(&lock_queue);
						update_far_neighbors(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == ACTIVE))
					{
						priority[neighbor] += SLOAN_W2;
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += SLOAN_W2;
						status[neighbor] = ACTIVE;
						
						omp_set_lock(&lock_queue);
						update_far_neighbors(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
					}
				}
				
				#pragma omp critical
				{
// 					printf("Setting node %d at position %d of permutation\n", vertex, next_id);fflush(stdout);
					permutation[next_id++] = vertex;
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
		free(distance);
		
		#pragma omp single nowait
		free(priority);
		
		#pragma omp single nowait
		free(status);
		
		#pragma omp single nowait
		*Fp = permutation;
	}
	
// 	printf("Permutation vector of size %d: ", n);fflush(stdout);
// 	int i;
// 	for (i = 0; i < n; ++i)
// 		printf("%d ", (*Fp)[i]);
// 	printf("\n");fflush(stdout);
}
