/*----------------------------------------------------------------------------
 * Sloan REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN reordering
 *--------------------------------------------------------------------------*/

void inline update_far_neighbors_v1(MAT* A, int** status, int** priority, int node, LIST** worklist) 
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

void inline update_far_neighbors(MAT* A, int** status, int** priority, int node, LIST*** worklist) 
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
			(*worklist)[(*priority)[neighbor]] = 
				LIST_insert_IF_NOT_EXIST((*worklist)[(*priority)[neighbor]], neighbor);
		}
		printf("updating priority of node %d\n", node);fflush(stdout);
		if ((*worklist)[(*priority)[node]] != NULL)
			(*worklist)[(*priority)[node]] = 
					LIST_remove((*worklist)[(*priority)[node]], node);
		(*priority)[node] += SLOAN_W2;
		(*worklist)[(*priority)[node]] = 
				LIST_insert_IF_NOT_EXIST((*worklist)[(*priority)[node]], node);
	}
	
	free(neighbors);
}


void Parallel_Sloan (MAT* adjacency, int** Fp, int start_node, int end_node)
{
	int node, next_id, num_nodes, min_priority, size_prior_bags;
	LIST** priority_bags;
	omp_lock_t lock_queue;
	
	num_nodes        = adjacency->n;
	int* distance    = calloc (num_nodes,sizeof (int));
	int* permutation = calloc (num_nodes,sizeof (int));
	int* priority    = calloc (num_nodes,sizeof (int));
	int* status      = calloc (num_nodes,sizeof (int));
	
	GRAPH_bfs(adjacency, end_node, distance);
	next_id = 0;
	min_priority = INFINITY_LEVEL;
	printf("**************Sloan start*****************\n");fflush(stdout);	
	#pragma omp parallel 
	{
		#pragma omp for private(node)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(GRAPH_degree(adjacency, node) + 1);
			if (priority[node] < min_priority) min_priority = priority[node];
		}

// 		int i;
// 		printf("Priorities before: ");fflush(stdout);
// 		for (i = 0; i< num_nodes; i++) printf("%d ", priority[i]); fflush(stdout);
// 		printf("\n");fflush(stdout);
		
		#pragma omp single
		if (min_priority > 0) 
		{
			printf("*** [Error] Sloan minimum initial priority is higher than 0 ***\n");
			exit(1);
		}
		
		#pragma omp for private(node)
		for (node = 0; node < num_nodes; ++node) 
			priority[node] -= min_priority;
		
		int j;
		printf("Priorities After: ");fflush(stdout);
		for (j = 0; j< num_nodes; j++) printf("%d ", priority[j]); fflush(stdout);
		printf("\n");fflush(stdout);
		
		printf("Priority start node %d: %d\n", start_node, priority[start_node]);fflush(stdout);
		
		#pragma omp single
		{
			size_prior_bags = 2 * priority[start_node];
			priority_bags = calloc(size_prior_bags, sizeof(LIST));
		}
		
		#pragma omp single nowait
		priority_bags[priority[start_node]] = 
				LIST_insert_IF_NOT_EXIST(priority_bags[priority[start_node]], start_node);
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		omp_init_lock(&lock_queue);
	}
	
	#pragma omp parallel 
	{
		int vertex, vertex_degree, neighbor, ngb, prior;
		int* neighbors;
		
		while (next_id < num_nodes)
		{
			vertex = UNDEF_NODE;
			
// 			printf("worklist: ");LIST_print(queue);fflush(stdout);
			
			omp_set_lock(&lock_queue);
			
			// Chosing vertex that maximizes the priority
			for (prior = size_prior_bags-1; prior >=0; --prior)
			{
// 				printf("Searching vertex at position %d\n", prior);
				
				if (priority_bags[prior] != NULL && priority_bags[prior]->size > 0)
				{
					vertex = priority_bags[prior]->data;
					priority_bags[prior] = LIST_remove(priority_bags[prior], vertex);
					prior = 0;
				}
			}
			
			omp_unset_lock(&lock_queue);
			
// 			printf("Selecteds vertex %d\n", vertex);fflush(stdout);
			
			if (vertex != UNDEF_NODE) 
			{
				printf("Processing vertex %d\n", vertex);fflush(stdout);
				
				neighbors     = GRAPH_adjacent(adjacency, vertex);
				vertex_degree = GRAPH_degree  (adjacency, vertex);
				
				int i;
				printf("Processing neighbors: ");fflush(stdout);
				for (i = 0; i< vertex_degree; i++) printf("%d ", neighbors[i]); fflush(stdout);
				printf("\n");fflush(stdout);
				
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
						{
							priority[neighbor] += SLOAN_W2; 
							status[neighbor]    = ACTIVE;
							priority_bags[priority[neighbor]] = 
								LIST_insert_IF_NOT_EXIST(priority_bags[priority[neighbor]], neighbor);
						}
						
						update_far_neighbors(adjacency, &status, &priority, neighbor, &priority_bags);
						omp_unset_lock(&lock_queue);
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == PREACTIVE))
					{
						omp_set_lock(&lock_queue);
						
						priority_bags[priority[neighbor]] = 
							LIST_remove(priority_bags[priority[neighbor]], neighbor);
						priority[neighbor] += SLOAN_W2;
						priority_bags[priority[neighbor]] = 
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[neighbor]], neighbor);
						status[neighbor] = ACTIVE;
						
						update_far_neighbors(adjacency, &status, &priority, neighbor, &priority_bags);
						omp_unset_lock(&lock_queue);
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == ACTIVE))
					{
						omp_set_lock(&lock_queue);
						priority_bags[priority[neighbor]] = 
							LIST_remove(priority_bags[priority[neighbor]], neighbor);
						priority[neighbor] += SLOAN_W2;
						priority_bags[priority[neighbor]] = 
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[neighbor]], neighbor);
						omp_unset_lock(&lock_queue);
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						omp_set_lock(&lock_queue);
						priority_bags[priority[neighbor]] = 
							LIST_remove(priority_bags[priority[neighbor]], neighbor);
						priority[neighbor] += SLOAN_W2;
						priority_bags[priority[neighbor]] = 
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[neighbor]], neighbor);
							
						status[neighbor] = ACTIVE;
						
						
						update_far_neighbors(adjacency, &status, &priority, neighbor, &priority_bags);
						omp_unset_lock(&lock_queue);
					}
				}
				
				#pragma omp critical
				{
					printf("Setting node %d at position %d of permutation\n", vertex, next_id);fflush(stdout);
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
			max_priority = MIN_PRIORITY;
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
						
						update_far_neighbors_v1(adjacency, &status, &priority, neighbor, &queue);
						omp_unset_lock(&lock_queue);
						
						priority[neighbor] += SLOAN_W2; 
						status[neighbor] = ACTIVE;
					}
					else if ((status[vertex] == PREACTIVE) && (status[neighbor] == PREACTIVE))
					{
						priority[neighbor] += SLOAN_W2; 
						status[neighbor] = ACTIVE;
						
						omp_set_lock(&lock_queue);
						update_far_neighbors_v1(adjacency, &status, &priority, neighbor, &queue);
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
						update_far_neighbors_v1(adjacency, &status, &priority, neighbor, &queue);
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
