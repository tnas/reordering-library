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
	
	for (ngb = 0; ngb < node_degree; ++ngb) 
	{
		neighbor = neighbors[ngb];
		
		if ((*status)[neighbor] == INACTIVE)
		{
			(*status)[neighbor] = PREACTIVE;
			(*worklist)[(*priority)[neighbor]] = 
				LIST_insert_IF_NOT_EXIST((*worklist)[(*priority)[neighbor]], neighbor);
		}
		
		if ((*worklist)[(*priority)[node]] != NULL)
			(*worklist)[(*priority)[node]] = 
					LIST_remove((*worklist)[(*priority)[node]], node);
		(*priority)[node] += SLOAN_W2;
		(*worklist)[(*priority)[node]] = 
				LIST_insert_IF_NOT_EXIST((*worklist)[(*priority)[node]], node);
	}
	
	free(neighbors);
}


void inline update_far_log_neighbors(MAT* A, int** status, int node, LIST** local_log) 
{
	int neighbor_degree, far_ngb, far_neighbor;
	int* far_neighbors;
	
	far_neighbors   = GRAPH_adjacent(A, node);
	neighbor_degree = GRAPH_degree(A, node);
	
	for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
	{
		far_neighbor = far_neighbors[far_ngb];
		
		if ((*status)[far_neighbor] == INACTIVE)
		{
			#pragma omp critical
			(*status)[far_neighbor] = PREACTIVE;
			
			*local_log = LIST_add_IF_NOT_EXIST(*local_log, far_neighbor, SLOAN_W2, omp_get_thread_num());
		}
		
		*local_log = LIST_add_IF_NOT_EXIST(*local_log, node, SLOAN_W2, omp_get_thread_num());
	}
	
	free(far_neighbors);
}



void Parallel_Sloan (MAT* adjacency, int** Fp, int start_node, int end_node)
{
	int num_nodes, next_id, min_priority, max_priority, size_prior_bags, 
	    threads_working, num_threads;
	LIST** priority_bags;
	int* distance;
	int* permutation;
	int* priority;
	int* status;
	int* status_threads;
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor, ngb, prior,
		    neighbor_degree, far_ngb, far_neighbor, update_far, thread, id;
		int* neighbors;
		int* far_neighbors;
		LIST* log;
		LIST* p_log;
	
		#pragma omp single 
		{
			num_threads    = omp_get_num_threads();
			status_threads = calloc(num_threads, sizeof(int));
		}
		
		#pragma omp single nowait
		num_nodes = adjacency->n;
		
		#pragma omp single nowait
		permutation = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		priority = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		status = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		min_priority = INFINITY_LEVEL;
		
		#pragma omp single
		{
			distance = calloc (num_nodes, sizeof (int));
			// TODO: Usar a BFS paralela
			GRAPH_bfs(adjacency, end_node, distance);
		}
		
		#pragma omp for 
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(GRAPH_degree(adjacency, node) + 1);
			if (priority[node] < min_priority) min_priority = priority[node];
		}

		#pragma omp single
		if (min_priority > 0) 
		{
			printf("*** [Error] Sloan minimum initial priority is higher than 0 ***\n");
			exit(1);
		}
		
		#pragma omp for 
		for (node = 0; node < num_nodes; ++node) 
			priority[node] -= min_priority;
		
		#pragma omp single
		{
			size_prior_bags = 2 * priority[start_node];
			priority_bags = calloc(size_prior_bags, sizeof(LIST*));
			printf(">>Thread %d - priority start node %d: %d\n", omp_get_thread_num(), start_node, priority[start_node]);fflush(stdout);
		}
		
		#pragma omp single nowait
		priority_bags[priority[start_node]] = 
				LIST_insert_IF_NOT_EXIST(priority_bags[priority[start_node]], start_node);
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		#pragma omp single nowait
		next_id = 0;
		
		#pragma omp single nowait
		threads_working = 0;
		
		#pragma omp for
		for (thread = 0; thread < num_threads; ++thread)
			status_threads[thread] = THREAD_OFF;
		
		#pragma omp barrier
		
		thread = omp_get_thread_num();
		
		while (next_id < num_nodes)
		{
// 			printf(">>Thread %d working over the next_id = %d\n", omp_get_thread_num(), next_id);fflush(stdout);
			
			if (!threads_working)
			{
				// Chosing maximum priority
				#pragma omp single
				{
					for (prior = size_prior_bags - 1; prior >= 0; --prior)
					{
						if (priority_bags[prior] != NULL)
						{
							max_priority = prior;
							prior = 0;
						}
					}
					
					printf(">>Thread %d set maximum priority: %d\n", omp_get_thread_num(), max_priority);fflush(stdout);
				}
				
				log = NULL;
				
				while (priority_bags[max_priority] != NULL)
				{
// 					#pragma omp single
// 					{
// 						printf("Priority Bag: ");LIST_print(priority_bags[max_priority]);fflush(stdout);
// 					}
					
					printf(">>Thread %d start processing\n", omp_get_thread_num());fflush(stdout);
					vertex = UNDEF_NODE;
					
					#pragma omp critical
					{
						if (priority_bags[max_priority] != NULL)
						{
							vertex = priority_bags[max_priority]->data;
							priority_bags[max_priority] = LIST_remove(priority_bags[max_priority], vertex);
							status_threads[thread] = THREAD_ON;
							printf(">>Thread %d getting node %d - STATUS %d\n", omp_get_thread_num(), vertex, status_threads[thread]);fflush(stdout);
						}
					}
					
					if (vertex != UNDEF_NODE)
					{
						neighbors     = GRAPH_adjacent(adjacency, vertex);
						vertex_degree = GRAPH_degree  (adjacency, vertex);
						
						for (ngb = 0; ngb < vertex_degree; ++ngb)
						{
							update_far = 0;
							neighbor   = neighbors[ngb];
							
							if (status[vertex] == PREACTIVE)
							{
								if (status[neighbor] == ACTIVE)
								{
									log = LIST_add_IF_NOT_EXIST(log, neighbor, SLOAN_W2, omp_get_thread_num());
								}
								else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
								{
									log = LIST_add_IF_NOT_EXIST(log, neighbor, SLOAN_W2, omp_get_thread_num());
									
									#pragma omp critical
									status[neighbor] = ACTIVE;
									
									update_far = 1;
								}
							}
							else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
							{
								log = LIST_add_IF_NOT_EXIST(log, neighbor, SLOAN_W2, omp_get_thread_num());
									
								#pragma omp critical
								status[neighbor] = ACTIVE;
								
								update_far = 1;
							}
					
							if (update_far)
							{
								/* ***********************
								* Updating far neighbors
								* ***********************
								*/
								
								far_neighbors   = GRAPH_adjacent(adjacency, neighbor);
								neighbor_degree = GRAPH_degree(adjacency, neighbor);
								
								for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
								{
									far_neighbor = far_neighbors[far_ngb];
									
									if (status[far_neighbor] == INACTIVE)
									{
										#pragma omp critical
										status[far_neighbor] = PREACTIVE;
										
										log = LIST_add_IF_NOT_EXIST(log, far_neighbor, SLOAN_W2, omp_get_thread_num());
									}
									
									log = LIST_add_IF_NOT_EXIST(log, neighbor, SLOAN_W2, omp_get_thread_num());
								}
								
								free(far_neighbors);
							}
						}
						
						free(neighbors);
						
						#pragma omp critical
						{
							printf("Thread %d setting node %d at position %d of permutation\n", omp_get_thread_num(), vertex, next_id);fflush(stdout);
							permutation[next_id++] = vertex;
							status[vertex] = NUMBERED;
						}
					}
				}
				
				#pragma omp critical
				{
					while (log != NULL)
					{
						node = log->data;
						printf(">>Thread %d updating node %d\n", omp_get_thread_num(), node);fflush(stdout);
						if (priority_bags[priority[node]] != NULL)
							priority_bags[priority[node]] = 
								LIST_remove(priority_bags[priority[node]], node);
						priority[node] += log->value;
						priority_bags[priority[node]] =	
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[node]], node);
							
						p_log = log->next;
						free(log);
						log = p_log;
					}
					
					if (status_threads[thread] == THREAD_ON)
						status_threads[thread] = THREAD_OFF;
					
					printf(">>Thread %d - STATUS %d\n", omp_get_thread_num(), status_threads[thread]);fflush(stdout);
					
				}
				
				#pragma omp critical
				{
					printf("status threads array: ");
					for (id = 0; id < num_threads; ++id) printf("");
					printf("\n");
					
					threads_working = 0;
					for (id = 0; id < num_threads; ++id)
						threads_working |= status_threads[id]; 
					printf("**Thread %d - threads working = %d\n", omp_get_thread_num(), threads_working);fflush(stdout);
				}
				
				
			}
		}
		
		#pragma omp single nowait
		free(distance);
		
		#pragma omp single nowait
		free(priority);
		
		#pragma omp single nowait
		free(status);
		
		#pragma omp single nowait
		free(priority_bags);
		
		#pragma omp single nowait
		free(status_threads);
		
		#pragma omp single nowait
		*Fp = permutation;
	}
	
	printf("Permutation vector of size %d: ", num_nodes);fflush(stdout);
	int i;
	for (i = 0; i < num_nodes; ++i)
		printf("%d ", (*Fp)[i]);
	printf("\n");fflush(stdout);
}




void Parallel_Sloan_v2 (MAT* adjacency, int** Fp, int start_node, int end_node)
{
	int node, next_id, num_nodes, min_priority, size_prior_bags;
	LIST** priority_bags;
	omp_lock_t lock_queue;
	
	// TODO: paralelizar as alocações
	num_nodes        = adjacency->n;
	int* distance    = calloc (num_nodes,sizeof (int));
	int* permutation = calloc (num_nodes,sizeof (int));
	int* priority    = calloc (num_nodes,sizeof (int));
	int* status      = calloc (num_nodes,sizeof (int));
	
	// TODO: Usar a BFS paralela
	GRAPH_bfs(adjacency, end_node, distance);
	next_id = 0;
	min_priority = INFINITY_LEVEL;
	
// 	printf("**************Sloan start*****************\n");fflush(stdout);	
	
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
		
// 		printf("Priority start node %d: %d\n", start_node, priority[start_node]);fflush(stdout);
		
		#pragma omp single
		{
			size_prior_bags = 2 * priority[start_node];
			priority_bags = calloc(size_prior_bags, sizeof(LIST*));
		}
		
		#pragma omp single nowait
		priority_bags[priority[start_node]] = 
				LIST_insert_IF_NOT_EXIST(priority_bags[priority[start_node]], start_node);
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		omp_init_lock(&lock_queue);
	}
	
	// TODO: eliminar essa divisão e fazer um único bloco paralelo
	
	#pragma omp parallel 
	{
		int vertex, vertex_degree, neighbor, ngb, prior;
		int* neighbors;
		
		while (next_id < num_nodes)
		{
			
// 			int j;
// 			printf("Priorities : ");fflush(stdout);
// 			for (j = 0; j< num_nodes; j++) printf("%02d ", priority[j]); fflush(stdout);
// 			printf("\n");fflush(stdout);
// 			printf("Status     : ");fflush(stdout);
// 			for (j = 0; j< num_nodes; j++) printf("%02d ", status[j]); fflush(stdout);
// 			printf("\n");fflush(stdout);
			
			
			
// 			printf("worklist: ");LIST_print(queue);fflush(stdout);
			
			omp_set_lock(&lock_queue);
			
// 			for (prior = size_prior_bags-1; prior >=0; --prior)
// 			{
// 				if (priority_bags[prior] != NULL)
// 				{
// 					printf("[Thread %d] Priority bag %d: ", omp_get_thread_num(), prior);
// 					LIST_print(priority_bags[prior]);
// 					fflush(stdout);
// 				}
// 			}

			vertex = UNDEF_NODE;
			
			
			// Chosing vertex that maximizes the priority
			
			for (prior = size_prior_bags-1; prior >=0; --prior)
			{
				if (priority_bags[prior] != NULL)
				{
// 					printf("Thread %d getting vertex %d at priority %d\n", omp_get_thread_num(), priority_bags[prior]->data, prior); fflush(stdout);
					vertex = priority_bags[prior]->data;
					priority_bags[prior] = LIST_remove(priority_bags[prior], vertex);
					prior = 0;
				}
			}
			
// 			// Chosing vertex that maximizes the priority
// 			for (prior = size_prior_bags-1; prior >=0; --prior)
// 			{
// 				if (priority_bags[prior] != NULL)
// 				{
// 					printf("Priority bag AFTER %d: ", prior);
// 					LIST_print(priority_bags[prior]);
// 					fflush(stdout);
// 				}
// 			}
			
			omp_unset_lock(&lock_queue);
			
// 			printf("Selecteds vertex %d\n", vertex);fflush(stdout);
			
			if (vertex != UNDEF_NODE) 
			{
// 				printf("---------------------------------------\n");
// 				printf("Processing vertex %d\n", vertex);
// 				printf("---------------------------------------\n");fflush(stdout);
				
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
						priority[neighbor] += SLOAN_W2; 
						status[neighbor]    = ACTIVE;
						priority_bags[priority[neighbor]] = 
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[neighbor]], neighbor);
						
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
					printf("Thread %d setting node %d at position %d of permutation\n", omp_get_thread_num(), vertex, next_id);fflush(stdout);
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
		free(priority_bags);
		
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
