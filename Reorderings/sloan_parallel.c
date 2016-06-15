/*----------------------------------------------------------------------------
 * Sloan REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../CommonFiles/protos_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN Reordering
 *--------------------------------------------------------------------------*/

void Parallel_Sloan (MAT* adjacency, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, min_priority, max_priority, num_prior_bags, 
	    num_threads, chunk_size;
	int* distance;
	int** priority;
	int* status;
	int* size_bags;
	int* degree;
	
	num_nodes = adjacency->n;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_bfs(adjacency, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor, ngb, prior_bag, n_node,
		    neighbor_degree, far_ngb, far_neighbor, update_far;
		int* neighbors;
		int* far_neighbors;
		
		#pragma omp single 
		{
			num_threads = omp_get_num_threads();
			chunk_size = num_nodes / num_threads;
			chunk_size = chunk_size > 0 ? chunk_size : 1;
		}
		
		#pragma omp single nowait
		*permutation = calloc (num_nodes, sizeof (int));
		
		#pragma omp single nowait
		priority = calloc (num_nodes, sizeof (int*));
		
		#pragma omp single nowait
		status = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		degree = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		min_priority = INFINITY_LEVEL;
		
		#pragma omp barrier
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			degree[node]   = GRAPH_degree(adjacency, node);
			priority[node] = calloc(2, sizeof(int));
			priority[node][SLOAN_CURR_PRIOR] = SLOAN_W1*distance[node] - SLOAN_W2*(degree[node] + 1);
				
			#pragma omp critical
			if (priority[node][SLOAN_CURR_PRIOR] < min_priority) 
				min_priority = priority[node][SLOAN_CURR_PRIOR];
		}
		
		#pragma omp single nowait
		if (min_priority > 0) 
		{
			printf("*** [Error] Sloan minimum initial priority is higher than 0 ***\n");
			exit(1);
		}
		
		#pragma omp single
		{
			num_prior_bags = SLOAN_PRIORITY_FACTOR * 
				(priority[start_node][SLOAN_CURR_PRIOR] - min_priority);
			size_bags = calloc(num_prior_bags, sizeof(LIST*));
		}
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node) 
		{
			priority[node][SLOAN_CURR_PRIOR] -= min_priority;
			priority[node][SLOAN_NEW_PRIOR]   = priority[node][SLOAN_CURR_PRIOR];
			
			#pragma omp atomic
			++size_bags[priority[node][SLOAN_CURR_PRIOR]];
		}
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		next_id = 0;
		
		#pragma omp barrier
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while (next_id < num_nodes)
		{
			// Chosing maximum priority
			#pragma omp single
			{ 	
// 				for (node = 0; node < num_nodes; ++node) 
// 				{		
// 					if (status[node] == NUMBERED)
// 					{
// 						printf("Node %d => priority %d - N\n", node, priority[node][SLOAN_CURR_PRIOR]);fflush(stdout);
// 					}
// 					else 
// 					{
// 						printf("Node %d => priority %d\n", node, priority[node][SLOAN_CURR_PRIOR]);fflush(stdout);
// 					}
// 				}
// 				
// 				printf("----------------------------------\n");fflush(stdout);
// 				
// 				for (prior_bag = num_prior_bags - 1; prior_bag >= 0; --prior_bag)
// 				{
// 					if (size_bags[prior_bag] > 0)
// 					{
// 						printf("Priority bag: %d => size %d\n", prior_bag, size_bags[prior_bag]);fflush(stdout);
// 					}
// 				}
				
				
				for (prior_bag = num_prior_bags - 1; prior_bag >= 0; --prior_bag)
				{
					if (size_bags[prior_bag] > 0)
					{
						max_priority = prior_bag;
						prior_bag = 0;
					}
				}
			}
			
// 			printf("Max priority: %d\n", max_priority);fflush(stdout);
			
			// Processing the Logic Bag
			#pragma omp for schedule(static, chunk_size)
			for (n_node = 0; n_node < num_nodes; ++n_node)
			{
				if (status[n_node] != NUMBERED)
				{
					if (priority[n_node][SLOAN_CURR_PRIOR] == max_priority)
					{
// 						printf("Getting node %d from priority bag %d\n", n_node, max_priority);fflush(stdout);
						vertex        = n_node;
						neighbors     = GRAPH_adjacent(adjacency, vertex);
						vertex_degree = degree[vertex];
						
						for (ngb = 0; ngb < vertex_degree; ++ngb)
						{
							update_far = OFF;
							neighbor   = neighbors[ngb];
							
							if (status[vertex] == PREACTIVE)
							{
								if (status[neighbor] == ACTIVE)
								{
									#pragma omp atomic
									priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
								}
								else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
								{
									#pragma omp atomic
									priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
									
									#pragma omp critical
									status[neighbor] = ACTIVE;
									
									update_far = ON;
								}
								
								
// 								switch (status[neighbor])
// 								{
// 									case ACTIVE :
// 										
// 										#pragma omp atomic
// 										priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
// 										
// 										break;
// 										
// 									case INACTIVE :
// 										
// 										// Setting status to PREACTIVE
// 										#pragma omp atomic
// 										++status[neighbor];
// 									
// 									case PREACTIVE :
// 										
// 										// Setting status to ACTIVE
// 										#pragma omp atomic
// 										++status[neighbor];
// 										
// 										#pragma omp atomic
// 										priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
// 										
// 										update_far = ON;
// 										
// 										break;
// 								}
								
							}
							else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
							{
								#pragma omp atomic
								priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
								
// 								// Setting status to ACTIVE
// 								#pragma omp atomic
// 								++status[neighbor];
								
								#pragma omp critical
								status[neighbor] = ACTIVE;
								
								update_far = ON;
							}
					
							if (update_far)
							{
								/* ***********************
								* Updating far neighbors
								* ***********************
								*/
								
								far_neighbors   = GRAPH_adjacent(adjacency, neighbor);
								neighbor_degree = degree[neighbor];
								
								for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
								{
									far_neighbor = far_neighbors[far_ngb];
									
									if (status[far_neighbor] == INACTIVE)
									{
										// Setting status to PREACTIVE
// 										#pragma omp atomic
// 										++status[far_neighbor];
										
										#pragma omp critical
										status[far_neighbor] = PREACTIVE;
										
										#pragma omp atomic
										priority[far_neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
									}
									
									#pragma omp atomic
									priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
								}
								
								free(far_neighbors);
							}
						}
						
						free(neighbors);
						
						#pragma omp critical
						{
							size_bags[priority[vertex][SLOAN_CURR_PRIOR]]--;
// 							printf("Setting node %d at position %d => Size Bag %d  = %d\n", vertex, next_id, priority[vertex][SLOAN_CURR_PRIOR], size_bags[priority[vertex][SLOAN_CURR_PRIOR]]);fflush(stdout);
							(*permutation)[next_id++] = vertex;
							status[vertex] = NUMBERED;
							
						}
					}
				}
			}
			
			#pragma omp for schedule(static, chunk_size)
			for (node = 0; node < num_nodes; ++node)
			{
				if ((status[node] != NUMBERED) && 
					(priority[node][SLOAN_CURR_PRIOR] != priority[node][SLOAN_NEW_PRIOR]))
				{
// 					printf("Removing node %d from priority bag %d to %d\n", node, priority[node][SLOAN_CURR_PRIOR], priority[node][SLOAN_NEW_PRIOR]);fflush(stdout);
					
					#pragma omp atomic
					size_bags[priority[node][SLOAN_CURR_PRIOR]]--;
						
					#pragma omp atomic
					size_bags[priority[node][SLOAN_NEW_PRIOR]]++;
					
					priority[node][SLOAN_CURR_PRIOR] = priority[node][SLOAN_NEW_PRIOR];
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
		free(size_bags);
		
		#pragma omp single nowait
		free(degree);
	}
}


void Parallel_Sloan_second (MAT* adjacency, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, min_priority, max_priority, num_prior_bags, 
	    size_max_bag;
	LIST** priority_bags;
	LIST* p_bag;
	int* distance;
	int** priority;
	int* status;
	int* size_bags;
	int* current_bag;
	int* degree;
	
	num_nodes = adjacency->n;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_bfs(adjacency, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor, ngb, prior, n_node,
		    neighbor_degree, far_ngb, far_neighbor, update_far;
		int* neighbors;
		int* far_neighbors;
		
		#pragma omp single nowait
		*permutation = calloc (num_nodes, sizeof (int));
		
		#pragma omp single nowait
		priority = calloc (num_nodes, sizeof (int*));
		
		#pragma omp single nowait
		status = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		degree = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		min_priority = INFINITY_LEVEL;
		
		#pragma omp barrier
		
		#pragma omp for 
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			degree[node]   = GRAPH_degree(adjacency, node);
			priority[node] = calloc(2, sizeof(int));
			priority[node][SLOAN_CURR_PRIOR] = SLOAN_W1*distance[node] - SLOAN_W2*(degree[node] + 1);
				
			#pragma omp critical
			if (priority[node][SLOAN_CURR_PRIOR] < min_priority) 
				min_priority = priority[node][SLOAN_CURR_PRIOR];
		}
		
		#pragma omp single
		if (min_priority > 0) 
		{
			printf("*** [Error] Sloan minimum initial priority is higher than 0 ***\n");
			exit(1);
		}
		
		#pragma omp for 
		for (node = 0; node < num_nodes; ++node) 
		{
			priority[node][SLOAN_CURR_PRIOR] -= min_priority;
			priority[node][SLOAN_NEW_PRIOR] = priority[node][SLOAN_CURR_PRIOR];
		}
		
		#pragma omp single
		num_prior_bags = SLOAN_PRIORITY_FACTOR * priority[start_node][SLOAN_CURR_PRIOR];
		
		#pragma omp single
		{
			priority_bags = calloc(num_prior_bags, sizeof(LIST*));
			priority_bags[priority[start_node][SLOAN_CURR_PRIOR]] = 
				LIST_insert_IF_NOT_EXIST(priority_bags[priority[start_node][SLOAN_CURR_PRIOR]], start_node);
		}
		
		#pragma omp single nowait
		{
			size_bags = calloc(num_prior_bags, sizeof(LIST*));
			size_bags[priority[start_node][SLOAN_CURR_PRIOR]]++;
		}
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		next_id = 0;
		
		#pragma omp barrier
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while (next_id < num_nodes)
		{
			// Chosing maximum priority
			#pragma omp single
			{ 				
				for (prior = num_prior_bags - 1; prior >= 0; --prior)
				{
					if (priority_bags[prior] != NULL)
					{
						max_priority = prior;
						prior = 0;
					}
				}
				
				// Cloning bag for processing
				current_bag  = calloc(size_bags[max_priority], sizeof(int));
				size_max_bag = size_bags[max_priority];
				p_bag        = priority_bags[max_priority];
				n_node       = 0;
				
				while (p_bag != NULL)
				{
					current_bag[n_node++] = p_bag->data;
					p_bag = LIST_remove(p_bag, p_bag->data);
					size_bags[max_priority]--;
				}
				
				priority_bags[max_priority] = p_bag;
			}
			
			#pragma omp for
			for (n_node = 0; n_node < size_max_bag; ++n_node)
			{
				vertex        = current_bag[n_node];
				neighbors     = GRAPH_adjacent(adjacency, vertex);
				vertex_degree = degree[vertex];
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					update_far = 0;
					neighbor   = neighbors[ngb];
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							
							#pragma omp critical
							status[neighbor] = ACTIVE;
							
							update_far = 1;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						#pragma omp atomic
						priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
						
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
						neighbor_degree = degree[neighbor];
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (status[far_neighbor] == INACTIVE)
							{
								#pragma omp critical
								status[far_neighbor] = PREACTIVE;
								
								#pragma omp atomic
								priority[far_neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							}
							
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
						}
						
						free(far_neighbors);
					}
				}
				
				free(neighbors);
				
				#pragma omp critical
				{
					(*permutation)[next_id++] = vertex;
					status[vertex] = NUMBERED;
				}
			}
			
			#pragma omp single
			free(current_bag);
			
			#pragma omp for 
			for (node = 0; node < num_nodes; ++node)
			{
				if (priority[node][SLOAN_CURR_PRIOR] != priority[node][SLOAN_NEW_PRIOR])
				{
					#pragma omp critical
					{
						if (LIST_contains(priority_bags[priority[node][SLOAN_CURR_PRIOR]], node))
						{
							size_bags[priority[node][SLOAN_CURR_PRIOR]]--;
							priority_bags[priority[node][SLOAN_CURR_PRIOR]] = 
								LIST_remove(priority_bags[priority[node][SLOAN_CURR_PRIOR]], node);
						}
						
						if (priority[node][SLOAN_NEW_PRIOR] > num_prior_bags)
						{
							printf("*** [Error] Priority is higher than maximum number of bags ***\n");
							exit(1);
						}
						
						size_bags[priority[node][SLOAN_NEW_PRIOR]]++;
						priority_bags[priority[node][SLOAN_NEW_PRIOR]] =	
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[node][SLOAN_NEW_PRIOR]], node);
							
					}
					
					priority[node][SLOAN_CURR_PRIOR] = priority[node][SLOAN_NEW_PRIOR];
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
		free(size_bags);
		
		#pragma omp single nowait
		free(degree);
	}
}


void Parallel_Sloan_first (MAT* adjacency, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, min_priority, max_priority, num_prior_bags, 
	    num_threads, size_max_bag;
	LIST** priority_bags;
	LIST** logs;
	LIST* p_bag;
	int* distance;
	int* priority;
	int* status;
	int* size_bags;
	int* current_bag;
	int* degree;
	
	num_nodes = adjacency->n;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_bfs(adjacency, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor, ngb, prior, n_node,
		    neighbor_degree, far_ngb, far_neighbor, update_far, thread, id;
		int* neighbors;
		int* far_neighbors;
		LIST* p_log;
		
		#pragma omp single 
		num_threads = omp_get_num_threads();
		
		#pragma omp single nowait
		*permutation = calloc (num_nodes, sizeof (int));
		
		#pragma omp single nowait
		priority = calloc (num_nodes, sizeof (int));
		
		#pragma omp single nowait
		status = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		degree = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		logs = (LIST**) malloc(num_threads * sizeof(LIST*));
		
		#pragma omp single nowait
		min_priority = INFINITY_LEVEL;
		
		#pragma omp barrier
		
		#pragma omp for 
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			degree[node]   = GRAPH_degree(adjacency, node);
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(degree[node] + 1);
			#pragma omp critical
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
		num_prior_bags = SLOAN_PRIORITY_FACTOR * priority[start_node];
		
		#pragma omp single
		{
			priority_bags = calloc(num_prior_bags, sizeof(LIST*));
			priority_bags[priority[start_node]] = 
					LIST_insert_IF_NOT_EXIST(priority_bags[priority[start_node]], start_node);
		}
		
		#pragma omp single nowait
		{
			size_bags = calloc(num_prior_bags, sizeof(LIST*));
			size_bags[priority[start_node]]++;
		}
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		next_id = 0;
		
		thread = omp_get_thread_num();
		
		#pragma omp barrier
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while (next_id < num_nodes)
		{
			// Chosing maximum priority
			#pragma omp single
			{ 				
				for (prior = num_prior_bags - 1; prior >= 0; --prior)
				{
					if (priority_bags[prior] != NULL)
					{
						max_priority = prior;
						prior = 0;
					}
				}
				
				// Cloning bag for processing
				current_bag  = calloc(size_bags[max_priority], sizeof(int));
				size_max_bag = size_bags[max_priority];
				p_bag        = priority_bags[max_priority];
				n_node       = 0;
				
				while (p_bag != NULL)
				{
					current_bag[n_node++] = p_bag->data;
					p_bag = LIST_remove(p_bag, p_bag->data);
					size_bags[max_priority]--;
				}
				
				priority_bags[max_priority] = p_bag;
			}
			
			logs[thread] = NULL;
			
			#pragma omp for
			for (n_node = 0; n_node < size_max_bag; ++n_node)
			{
				vertex        = current_bag[n_node];
				neighbors     = GRAPH_adjacent(adjacency, vertex);
				vertex_degree = degree[vertex];
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					update_far = 0;
					neighbor   = neighbors[ngb];
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							logs[thread] = LIST_add_IF_NOT_EXIST(logs[thread], neighbor, SLOAN_W2);
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							logs[thread] = LIST_add_IF_NOT_EXIST(logs[thread], neighbor, SLOAN_W2);
							
							#pragma omp critical
							status[neighbor] = ACTIVE;
							
							update_far = 1;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						logs[thread] = LIST_add_IF_NOT_EXIST(logs[thread], neighbor, SLOAN_W2);
						
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
						neighbor_degree = degree[neighbor];
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (status[far_neighbor] == INACTIVE)
							{
								#pragma omp critical
								status[far_neighbor] = PREACTIVE;
								
								logs[thread] = LIST_add_IF_NOT_EXIST(logs[thread], far_neighbor, SLOAN_W2);
								
							}
							
							logs[thread] = LIST_add_IF_NOT_EXIST(logs[thread], neighbor, SLOAN_W2);
						}
						
						free(far_neighbors);
					}
				}
				
				free(neighbors);
				
				#pragma omp critical
				{
					(*permutation)[next_id++] = vertex;
					status[vertex] = NUMBERED;
				}
			}
			
			#pragma omp single
			free(current_bag);
			
			#pragma omp for schedule(static, 1)
			for (id = 0; id < num_threads; ++id)
			{
				while (logs[id] != NULL)
				{
					node = logs[id]->data;
					
					#pragma omp critical
					{
						if (LIST_contains(priority_bags[priority[node]], node))
						{
							size_bags[priority[node]]--;
							priority_bags[priority[node]] = 
								LIST_remove(priority_bags[priority[node]], node);
						}
						
						priority[node] += logs[id]->value;
						
						if (priority[node] > num_prior_bags)
						{
							printf("*** [Error] Priority is higher than maximum number of bags ***\n");
							exit(1);
						}
						
						size_bags[priority[node]]++;
						priority_bags[priority[node]] =	
							LIST_insert_IF_NOT_EXIST(priority_bags[priority[node]], node);
							
						p_log = logs[id]->next;
						free(logs[id]);
						logs[id] = p_log;
					}
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
		free(size_bags);
		
		#pragma omp single nowait
		free(logs);
		
		#pragma omp single nowait
		free(degree);
	}
}
