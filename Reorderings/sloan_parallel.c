/*
 * Copyright 2016 Thiago Nascimento nascimenthiago@gmail.com
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include "sloan_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN Reordering
 *--------------------------------------------------------------------------*/

void Parallel_Sloan_original(const METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, min_priority, max_priority, num_prior_bags, 
	    num_threads, chunk_size, count_threads_on;
	int* distance;
	int** priority;
	int* status;
	int* size_bags;
	int* degree;
	
	num_nodes = mgraph->size;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_static_BFS(mgraph, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor, ngb, prior_bag, 
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
			degree[node]   = mgraph->graph[node].degree;
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
			#pragma omp sections
			{
				// Chosing maximum priority == Defining the Logical Bag
				#pragma omp section
				{ 	
					for (prior_bag = num_prior_bags - 1; prior_bag >= 0; --prior_bag)
					{
						if (size_bags[prior_bag] > 0)
						{
							max_priority = prior_bag;
							prior_bag = 0; // stop seaching
						}
					}
				}
				
				#pragma omp section
				count_threads_on = 0;
			}
			
			#pragma omp for schedule(static, chunk_size)
			for (vertex = 0; vertex < num_nodes; ++vertex)
			{
				// Processing Logical Bag
				if ((status[vertex] != NUMBERED) && (priority[vertex][SLOAN_CURR_PRIOR] == max_priority))
				{
					neighbors     = mgraph->graph[vertex].neighboors;
					vertex_degree = mgraph->graph[vertex].degree;
					
					#pragma omp atomic
					count_threads_on++;
					
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
						}
						else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							
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
							
							far_neighbors   = mgraph->graph[neighbor].neighboors;
							neighbor_degree = mgraph->graph[neighbor].degree;
							
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
						}
					}
					
					#pragma omp atomic
					--count_threads_on;
					
					while (count_threads_on > 0); // Barrier
					
					#pragma omp critical
					{
						size_bags[priority[vertex][SLOAN_CURR_PRIOR]]--;
						(*permutation)[next_id++] = vertex;
						status[vertex] = NUMBERED;
					}
				}
			}
			
			#pragma omp for schedule(static, chunk_size)
			for (node = 0; node < num_nodes; ++node)
			{
				if ((status[node] != NUMBERED) && 
					(priority[node][SLOAN_CURR_PRIOR] != priority[node][SLOAN_NEW_PRIOR]))
				{
					#pragma omp atomic
					size_bags[priority[node][SLOAN_CURR_PRIOR]]--;
						
					#pragma omp atomic
					size_bags[priority[node][SLOAN_NEW_PRIOR]]++;
					
					priority[node][SLOAN_CURR_PRIOR] = priority[node][SLOAN_NEW_PRIOR];
				}
			}
		}
		
		#pragma omp barrier
		
		#pragma omp single nowait
		free(distance);
		
		#pragma omp single nowait
		free(status);
		
		#pragma omp single nowait
		free(size_bags);
		
		#pragma omp single nowait
		free(degree);
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
			free(priority[node]);
		
		#pragma omp single nowait
		free(priority);
	}
}


inline static void prioritize_node(const int* priority, const int node, const int num_nodes, priority_set* th_max_priority_set)
{
	if (priority[node] > th_max_priority_set->priority)
	{
		th_max_priority_set->priority = priority[node];
		th_max_priority_set->head = th_max_priority_set->tail = 0;
		QUEUE_enque(&th_max_priority_set->elements, num_nodes, &(th_max_priority_set->tail), node);
	}
	else if (priority[node] == th_max_priority_set->priority)
	{
		QUEUE_enque(&th_max_priority_set->elements, num_nodes, &(th_max_priority_set->tail), node);
	}
}


void Parallel_Sloan(const METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, num_threads, chunk_size;
	int* distance;
	int* priority;
	int* status;
	int* degree;
	omp_lock_t* lock_priority;
	priority_set* max_priority_set;
	
	num_nodes = mgraph->size;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_static_BFS(mgraph, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, index_vertex, vertex, vertex_degree, neighbor, ngb, neighbor_degree, 
		    far_ngb, far_neighbor, update_far;
		int* neighbors;
		int* far_neighbors;
		priority_set* th_max_priority_set;
		
		#pragma omp sections
		{
			#pragma omp section
			{
				num_threads = omp_get_num_threads();
				chunk_size = num_nodes / num_threads;
				chunk_size = chunk_size > 0 ? chunk_size : 1;
			}
			
			#pragma omp section
			*permutation = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			priority = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			lock_priority = malloc(num_nodes * sizeof(omp_lock_t));
			
			#pragma omp section
			status = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			degree = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			{
				max_priority_set           = malloc(sizeof(priority_set));
				max_priority_set->priority = -INFINITY_LEVEL;
				max_priority_set->head     = 0;
				max_priority_set->tail     = 0;
				max_priority_set->elements = calloc(num_nodes, sizeof(int));
			}
		}
		
		// Thread-local maximum priority set
		th_max_priority_set           = malloc(sizeof(priority_set));
		th_max_priority_set->priority = -INFINITY_LEVEL;
		th_max_priority_set->head     = 0;
		th_max_priority_set->tail     = 0;
		th_max_priority_set->elements = calloc(num_nodes, sizeof(int));
		
		/* ********************************
		 * ****** Pre-processing **********
		 * ********************************
		 */
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			omp_init_lock(&(lock_priority[node]));
			status[node]   = INACTIVE;
			degree[node]   = mgraph->graph[node].degree;
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(degree[node] + 1);
			printf("node: %d, priority: %d\n", node, priority[node]);fflush(stdout);
			
			// Thread-local maximmum priority
			if (priority[node] >= th_max_priority_set->priority)
				prioritize_node(priority, node, num_nodes, th_max_priority_set);
		}
		
		// Defining initial maximum priority 
		#pragma omp critical
		{
			if (th_max_priority_set->priority > max_priority_set->priority)
				max_priority_set->priority = th_max_priority_set->priority;
		}
		
		#pragma omp barrier
		
		// Defining initial maximum priority set
		#pragma omp critical
		{
			if (th_max_priority_set->priority == max_priority_set->priority)
			{
				for (node = th_max_priority_set->head; node < th_max_priority_set->tail; ++node)
					QUEUE_enque(&max_priority_set->elements, num_nodes, 
							&(max_priority_set->tail), th_max_priority_set->elements[node]);
			}
		}
		
		// Reseting each thread-local max priority set
		th_max_priority_set->priority = -INFINITY_LEVEL;
		th_max_priority_set->head = th_max_priority_set->tail = 0;
		
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
			// Processing maximum priority set of nodes
			#pragma omp for schedule(static, chunk_size)
			for (index_vertex = max_priority_set->head; 
			     index_vertex < max_priority_set->tail - max_priority_set->head; ++index_vertex)
			{
				vertex        = max_priority_set->elements[index_vertex];
				neighbors     = mgraph->graph[vertex].neighboors;
				vertex_degree = mgraph->graph[vertex].degree;
				printf("processing node %d with priority %d\n", vertex, priority[vertex]);fflush(stdout);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					update_far = OFF;
					neighbor   = neighbors[ngb];
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							omp_set_lock(&lock_priority[neighbor]);
							priority[neighbor] += SLOAN_W2;
							prioritize_node(priority, neighbor, num_nodes, th_max_priority_set);
							omp_unset_lock(&lock_priority[neighbor]);
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							omp_set_lock(&lock_priority[neighbor]);
							priority[neighbor] += SLOAN_W2;
							prioritize_node(priority, neighbor, num_nodes, th_max_priority_set);
							omp_unset_lock(&lock_priority[neighbor]);
							
							#pragma omp critical
							status[neighbor] = ACTIVE;
							
							update_far = ON;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						omp_set_lock(&lock_priority[neighbor]);
						priority[neighbor] += SLOAN_W2;
						prioritize_node(priority, neighbor, num_nodes, th_max_priority_set);
						omp_unset_lock(&lock_priority[neighbor]);
						
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
						
						far_neighbors   = mgraph->graph[neighbor].neighboors;
						neighbor_degree = mgraph->graph[neighbor].degree;
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (status[far_neighbor] == INACTIVE)
							{
								#pragma omp critical
								status[far_neighbor] = PREACTIVE;
								
								omp_set_lock(&lock_priority[far_neighbor]);
								priority[far_neighbor] += SLOAN_W2;
								prioritize_node(priority, far_neighbor, num_nodes, th_max_priority_set);
								omp_unset_lock(&lock_priority[far_neighbor]);
							}
							
							omp_set_lock(&lock_priority[far_neighbor]);
							priority[neighbor] += SLOAN_W2;
							prioritize_node(priority, far_neighbor, num_nodes, th_max_priority_set);
							omp_unset_lock(&lock_priority[far_neighbor]);
						}
					}
				}
				
				// Placing vertex in permutation array
				#pragma omp critical
				{
					(*permutation)[next_id++] = vertex;
					status[vertex] = NUMBERED;
				}
			}
			
			// Reseting max priority set
			#pragma omp single
			{
				max_priority_set->priority = -INFINITY_LEVEL;
				max_priority_set->head = max_priority_set->tail = 0;
			}
			
			// Defining next max priority set
			#pragma omp critical
			{
				if (th_max_priority_set->priority > max_priority_set->priority)
					max_priority_set->priority = th_max_priority_set->priority;
			}
			
			#pragma omp barrier
			
			#pragma omp critical
			{
				if (th_max_priority_set->priority == max_priority_set->priority)
				{
					for (node = th_max_priority_set->head; node < th_max_priority_set->tail; ++node)
						QUEUE_enque(&max_priority_set->elements, num_nodes, 
							    &(max_priority_set->tail), th_max_priority_set->elements[node]);
				}
			}
			
			// Reseting each thread-local max priority set
			th_max_priority_set->priority = -INFINITY_LEVEL;
			th_max_priority_set->head = th_max_priority_set->tail = 0;
		}
		
		free(th_max_priority_set->elements);
		free(th_max_priority_set);
		
		#pragma omp barrier
			
		#pragma omp for schedule(static, chunk_size) nowait
		for (node = 0; node < num_nodes; ++node)
			omp_destroy_lock(&(lock_priority[node]));
		
		#pragma omp sections
		{
			#pragma omp section
			free(distance);
			
			#pragma omp section
			free(status);
			
			#pragma omp section
			free(degree);
			
			#pragma omp section
			{
				free(max_priority_set->elements);
				free(max_priority_set);
			}
			
			#pragma omp section
			free(priority);
		}
	}
}

