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

void Parallel_Logical_Bag_Sloan(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
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
	
// 	start_node = 2; end_node = 8;
	num_nodes = mgraph->size;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_static_BFS(mgraph, end_node, &distance, BFS_PERCENT_CHUNK);
	
// 	printf(">>start_node, end_node: %d(%d), %d(%d)\n", start_node, distance[start_node], end_node, distance[end_node]);fflush(stdout);

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
// 			printf("-------->>>>>>>>>>>>min_priority: %d\n", min_priority);fflush(stdout);
		}
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node) 
		{
			priority[node][SLOAN_CURR_PRIOR] -= min_priority;
			priority[node][SLOAN_NEW_PRIOR]   = priority[node][SLOAN_CURR_PRIOR];
		}
		
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		++size_bags[priority[start_node][SLOAN_CURR_PRIOR]];
		
		#pragma omp single nowait
		next_id = 0;
		
		#pragma omp barrier
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while (next_id < num_nodes)
		{
// 			printf("searching max priority...\n");fflush(stdout);
			
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
// 							printf("setting max priority as %d\n", max_priority);fflush(stdout);
							prior_bag = 0; // stop seaching
						}
					}
				}
				
				#pragma omp section
				count_threads_on = 0;
			}
			
// 			printf("Priorities: ");
// 			for (vertex = 0; vertex < num_nodes; ++vertex)
// 			{
// 				if ((status[vertex] == ACTIVE) || (status[vertex] == PREACTIVE))
// 				{
// 					printf("%d(*%d), ", vertex, priority[vertex][SLOAN_CURR_PRIOR]);fflush(stdout);
// 				}
// 				else 
// 				{
// 					printf("%d(%d), ", vertex, priority[vertex][SLOAN_CURR_PRIOR]);fflush(stdout);
// 				}
// 			}
// 			printf("\n");fflush(stdout);
// 			
// 			printf("max_priority: %d\n", max_priority);fflush(stdout);
			
			#pragma omp for schedule(static, chunk_size) 
			for (vertex = 0; vertex < num_nodes; ++vertex)
			{
				// Processing Logical Bag
				if ((status[vertex] == ACTIVE || status[vertex] == PREACTIVE) && (priority[vertex][SLOAN_CURR_PRIOR] == max_priority))
				{
					vertex_degree = mgraph->graph[vertex].degree;
					neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
					
					#pragma omp atomic
					count_threads_on++;
					
// 					printf("numbered vertex/priority: %d/%d\n", vertex, priority[vertex][SLOAN_CURR_PRIOR]);fflush(stdout);
					
					for (ngb = 0; ngb < vertex_degree; ++ngb)
					{
						update_far = OFF;
						neighbor   = neighbors[ngb];
						
// 						if (status[vertex] == PREACTIVE)
// 						{
// 							if (status[neighbor] == ACTIVE)
// 							{
// 								#pragma omp atomic
// 								priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
// 							}
// 							else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
// 							{
// 								#pragma omp atomic
// 								priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
// 								
// 								#pragma omp critical
// 								status[neighbor] = ACTIVE;
// 								
// 								update_far = ON;
// 							}
// 						}
// 						else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
// 						{
// 							#pragma omp atomic
// 							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
// 							
// 							#pragma omp critical
// 							status[neighbor] = ACTIVE;
// 							
// 							update_far = ON;
// 						}
						
						if (status[vertex] == PREACTIVE && (status[neighbor] == INACTIVE || status[neighbor] == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							
							#pragma omp critical
							status[neighbor] = ACTIVE;
							
							update_far = ON;
						}
						else if (status[vertex] == PREACTIVE && status[neighbor] == ACTIVE)
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
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
							
							neighbor_degree = mgraph->graph[neighbor].degree;
							far_neighbors   = GRAPH_neighboors(mgraph->mat, neighbor, neighbor_degree);
							
							for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
							{
								far_neighbor = far_neighbors[far_ngb];
								
								if (far_neighbor != vertex)
								{
									if (status[far_neighbor] == INACTIVE)
									{
										#pragma omp critical
										status[far_neighbor] = PREACTIVE;
									}
									
									#pragma omp atomic
									priority[far_neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
								}
							}
							
							free(far_neighbors);
						}
					}
					
					free(neighbors);
					
					#pragma omp atomic
					--count_threads_on;
					
					while (count_threads_on > 0); // Barrier
					
					#pragma omp critical
					{
						if (size_bags[priority[vertex][SLOAN_CURR_PRIOR]] > 0)
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
					if (size_bags[priority[node][SLOAN_CURR_PRIOR]] > 0)
					{
						#pragma omp critical
						{
							if (size_bags[priority[node][SLOAN_CURR_PRIOR]] > 0)
								size_bags[priority[node][SLOAN_CURR_PRIOR]]--;
						}
					}
					
					#pragma omp atomic
					size_bags[priority[node][SLOAN_NEW_PRIOR]]++;
					
					priority[node][SLOAN_CURR_PRIOR] = priority[node][SLOAN_NEW_PRIOR];
				}
			}
			
// 			for (prior_bag = num_prior_bags - 1; prior_bag >= 0; --prior_bag)
// 			{
// 				if (size_bags[prior_bag] < 0)
// 				{
// 					printf("Negative Bag: %d, Size: %d\n", prior_bag, size_bags[prior_bag]);fflush(stdout);
// 				}
// 			}
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


void Parallel_Bag_Sloan(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, num_threads, chunk_size, min_priority, max_priority, 
	    num_prior_bags, size_max_bag, max_processed_priority;
	bag* bags_priority;
	int* status;
	int* max_bag;
	
	num_nodes = mgraph->size;
	GRAPH_parallel_fixedpoint_sloan_BFS(mgraph, end_node, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, index_vertex, vertex, vertex_degree, neighbor, ngb, neighbor_degree, 
		    far_ngb, far_neighbor, update_far, n_bag, p_bag, tail_bag_change, c_node;
		int* neighbors;
		int* far_neighbors;
		GRAPH* bag_change;
		
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
			status = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			max_bag = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			min_priority = mgraph->min_sloan_priority;
			
			#pragma omp section
			max_priority = mgraph->max_sloan_priority;
		}
		
		/* ********************************
		 * ****** Pre-processing **********
		 * ** Generating Priority Bags **** 
		 * ********************************
		 */
		#pragma omp single
		{
			// Oversizing estimate
			num_prior_bags = SLOAN_PRIORITY_FACTOR * (mgraph->sloan_priority[start_node] - min_priority); 
			bags_priority = malloc(num_prior_bags * sizeof(bag));
		}
		
		// Initializing priority bags - So expensive loop!!!
		#pragma omp for 
		for (n_bag = 0; n_bag < num_prior_bags; ++n_bag)
		{
			ARRAY_LIST_init(&bags_priority[n_bag].list_elements);
			omp_init_lock(&bags_priority[n_bag].lock);
		}
		
		// Loading priority bags
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			mgraph->sloan_priority[node] -= min_priority;
			omp_set_lock(&(bags_priority[mgraph->sloan_priority[node]].lock));
			ARRAY_LIST_insert(&bags_priority[mgraph->sloan_priority[node]].list_elements, node);
			omp_unset_lock(&(bags_priority[mgraph->sloan_priority[node]].lock));
		}
		
		#pragma omp sections
		{
			#pragma omp section
			max_priority -= min_priority;
			
			#pragma omp section
			status[start_node] = PREACTIVE; 
			
			#pragma omp section
			next_id = 0;
		}
		
		// Bag that stores each node change per thread
		bag_change = malloc((num_nodes/2) * sizeof(GRAPH));
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while (next_id < num_nodes)
		{
			// Cloning current max priority bag
			#pragma omp single	
			{
				size_max_bag = bags_priority[max_priority].list_elements->size;
				c_node       = 0;
				
				while (bags_priority[max_priority].list_elements->size > 0)
				{
					max_bag[c_node++] = 
						ARRAY_LIST_remove_first(&bags_priority[max_priority].list_elements);
				}
			}
			
			tail_bag_change = 0;
			
			#pragma omp single
			max_processed_priority = -INFINITY_LEVEL;
			
			// Processing maximum priority bag of nodes
			#pragma omp for 
			for (index_vertex = 0; index_vertex < size_max_bag; ++index_vertex)
			{
				vertex = max_bag[index_vertex];
				
				if (status[vertex] == NUMBERED) 
					continue;
				
				vertex_degree = mgraph->graph[vertex].degree;
				neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor   = neighbors[ngb];
					update_far = OFF;
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							bag_change[tail_bag_change].label    = neighbor;
							bag_change[tail_bag_change].status   = ACTIVE;
							bag_change[tail_bag_change].distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
							++tail_bag_change;
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							bag_change[tail_bag_change].label    = neighbor;
							bag_change[tail_bag_change].status   = ACTIVE;
							bag_change[tail_bag_change].distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
							++tail_bag_change;
							
							update_far = ON;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						bag_change[tail_bag_change].label    = neighbor;
						bag_change[tail_bag_change].status   = ACTIVE;
						bag_change[tail_bag_change].distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
						++tail_bag_change;
						
						update_far = ON;
					}
			
					if (update_far)
					{
						/* ***********************
						* Updating far neighbors
						* ***********************
						*/
						
						far_neighbors   = GRAPH_adjacent(mgraph->mat, neighbor);
						neighbor_degree = mgraph->graph[neighbor].degree;
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (far_neighbor == vertex) continue;
							
							bag_change[tail_bag_change].label = far_neighbor;
							
							if (status[far_neighbor] == INACTIVE)
							{
								bag_change[tail_bag_change].status   = PREACTIVE;
								bag_change[tail_bag_change].distance = mgraph->sloan_priority[far_neighbor] + SLOAN_W2;
							}
							else
							{
								bag_change[tail_bag_change].status   = status[far_neighbor];
								bag_change[tail_bag_change].distance = mgraph->sloan_priority[far_neighbor];
							}
							
							bag_change[tail_bag_change].distance += SLOAN_W2;
							++tail_bag_change;
						}
						
						free(far_neighbors);
					}
				}
				
				free(neighbors);
				
				// Placing vertex in permutation array
				#pragma omp critical
				{
					if (status[vertex] != NUMBERED)
					{
						(*permutation)[next_id++] = vertex;
						status[vertex] = NUMBERED;
					}
				}
			} // priority bag loop
			
			// Unloading each bag change
			if (tail_bag_change > 0)
			{
				#pragma omp critical
				{
					for (index_vertex = 0; index_vertex < tail_bag_change; ++index_vertex)
					{
						vertex = bag_change[index_vertex].label;
						
						if (bag_change[index_vertex].status > status[vertex])
							status[vertex] = bag_change[index_vertex].status;
						
  						// Defining max processed priority
						if (bag_change[index_vertex].distance > max_processed_priority) 
							max_processed_priority = bag_change[index_vertex].distance;
						
						mgraph->sloan_priority[vertex] = bag_change[index_vertex].distance;
						
						ARRAY_LIST_insert(&bags_priority[mgraph->sloan_priority[vertex]].list_elements, vertex);
					}
				}
			}
			
			// Defining next max priority bag
			#pragma omp single
			{
				if (max_processed_priority >= max_priority)
				{
					max_priority = max_processed_priority;
				}
				else
				{
					// Searching next max priority bag
					for (p_bag = max_priority; p_bag >= 0; --p_bag)
					{
						if (bags_priority[p_bag].list_elements->size > 0)
						{
							// Find out next max priority bag
							max_priority = p_bag;
							p_bag = 0;
						}
					}
				}
			}
		} // main loop - while
		
		free(bag_change);
		
		#pragma omp barrier
		
		/* ********************************
		 * ****** Post-processing *********
		 * ********************************
		 */
		
		#pragma omp for 
		for (n_bag = 0; n_bag < num_prior_bags; ++n_bag)
		{
			ARRAY_LIST_destroy(&bags_priority[n_bag].list_elements);
			omp_destroy_lock(&bags_priority[n_bag].lock);
		}
		
		#pragma omp sections
		{
			#pragma omp section
			free(status);
			
			#pragma omp section
			free(bags_priority);
			
			#pragma omp section
			free(max_bag);
		}
	}
}


int COMPARE_priority_DESC(const void* st, const void* nd)
{ 
	volatile SLOAN_GRAPH* g_st = (SLOAN_GRAPH*) st;
	volatile SLOAN_GRAPH* g_nd = (SLOAN_GRAPH*) nd;
		       
	if (g_st->priorities[g_st->label] > g_nd->priorities[g_nd->label]) 
	{
		return -1;
	}
	if (g_st->priorities[g_st->label] < g_nd->priorities[g_nd->label])
	{
		return 1;
	}
	
	return 0;
}



void Parallel_Sloan(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, prior_head, prior_tail, pqueue_size;
	int* status;
	SLOAN_GRAPH* priority_queue;
	omp_lock_t queue_lock;
	
	num_nodes = mgraph->size;
	pqueue_size = (omp_get_max_threads() + 1) * num_nodes; // oversizing estimate
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor_degree, ngb, far_ngb, neighbor, 
		    far_neighbor, update_far, th_head, th_tail, size_chunk, 
		    count_chunk, snap_head, snap_tail, dirty_head, dirty_tail;
		SLOAN_GRAPH active_node, dirty_node;
		int* neighbors;
		int* far_neighbors;
		int* priority_snapshot;
		SLOAN_GRAPH* dirty_priority;
		
		/* ********************************
		 * ****** Pre-processing **********
		 * ********************************
		 */
		
		#pragma omp sections
		{
			#pragma omp section
			*permutation = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			status = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			priority_queue = calloc(pqueue_size, sizeof(SLOAN_GRAPH));
			
			#pragma omp section
			next_id = prior_head = prior_tail = 0;
			
			#pragma omp section
			omp_init_lock(&queue_lock);
		}
		
		#pragma omp for private(node)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]                    = INACTIVE;
			mgraph->sloan_priority[node]   -= mgraph->min_sloan_priority;  // Normalizing priorities at 0
		}
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single
		{
			dirty_node.label      = start_node;
			dirty_node.status     = PREACTIVE;
			dirty_node.distance   = mgraph->sloan_priority[start_node];
			dirty_node.priorities = mgraph->sloan_priority;
			GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
		}
		
		priority_snapshot = calloc(num_nodes, sizeof(int));
		snap_head = snap_tail = 0;
		
		dirty_priority = calloc(num_nodes, sizeof(SLOAN_GRAPH));
		dirty_head = dirty_tail = 0;
		dirty_node.priorities = mgraph->sloan_priority;
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while ((!QUEUE_empty(priority_queue, prior_head, prior_tail)) && next_id < num_nodes)
		{
			
			if (!QUEUE_empty(priority_queue, prior_head, prior_tail)) 
			{
				omp_set_lock(&queue_lock);
				
				if (prior_head < prior_tail)
					qsort(&priority_queue[prior_head], prior_tail-prior_head, sizeof(SLOAN_GRAPH), COMPARE_priority_DESC);
				
				th_tail = prior_tail;
				th_head = prior_head;	
				size_chunk = ceil(BFS_PERCENT_CHUNK * (th_tail - th_head));
				prior_head += size_chunk;
				
				omp_unset_lock(&queue_lock);
					
				for (count_chunk = 0; count_chunk < size_chunk; ++count_chunk)
				{
					active_node = GRAPH_deque(&priority_queue, pqueue_size, &th_head);
					QUEUE_enque(&priority_snapshot, num_nodes, &snap_tail, active_node.label);
				}
			}
			
			// Processing snapshot of maximum priority nodes
			while (!QUEUE_empty(priority_snapshot, snap_head, snap_tail))
			{
				vertex        = QUEUE_deque(&priority_snapshot, num_nodes, &snap_head);
				vertex_degree = mgraph->graph[vertex].degree;
				neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor   = neighbors[ngb];
					update_far = OFF;
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							dirty_node.label    = neighbor;
							dirty_node.status   = ACTIVE;
							dirty_node.distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
							GRAPH_enque(&dirty_priority, num_nodes, &dirty_tail, dirty_node);
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							dirty_node.label    = neighbor;
							dirty_node.status   = ACTIVE;
							dirty_node.distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
							GRAPH_enque(&dirty_priority, num_nodes, &dirty_tail, dirty_node);
							
							update_far = ON;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						dirty_node.label    = neighbor;
						dirty_node.status   = ACTIVE;
						dirty_node.distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
						GRAPH_enque(&dirty_priority, num_nodes, &dirty_tail, dirty_node);
						
						update_far = ON;
					}
			
					if (update_far)
					{
						/* ***********************
						* Updating far neighbors
						* ***********************
						*/
						
						far_neighbors   = GRAPH_adjacent(mgraph->mat, neighbor);
						neighbor_degree = mgraph->graph[neighbor].degree;
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (far_neighbor == vertex) continue;
							
							dirty_node.label = far_neighbor;
							
							if (status[far_neighbor] == INACTIVE)
							{
								dirty_node.status   = PREACTIVE;
								dirty_node.distance = mgraph->sloan_priority[far_neighbor] + SLOAN_W2;
							}
							else
							{
								dirty_node.status   = status[far_neighbor];
								dirty_node.distance = mgraph->sloan_priority[far_neighbor];
							}
							
							dirty_node.distance += SLOAN_W2;
							
							GRAPH_enque(&dirty_priority, num_nodes, &dirty_tail, dirty_node);
						}
						
						free(far_neighbors);
					}
				}
				
				free(neighbors);
				
				// Placing vertex in permutation array
				#pragma omp critical
				{
					if (status[vertex] != NUMBERED)
					{
						(*permutation)[next_id++] = vertex;
						status[vertex] = NUMBERED;
					}
				}
			} // priority bag loop
			
			snap_head = snap_tail = 0;
			
			// Updating priorities and status from dirty nodes
			if (!QUEUE_empty(dirty_priority, dirty_head, dirty_tail))
			{
				while (!QUEUE_empty(dirty_priority, dirty_head, dirty_tail))
				{
					dirty_node = GRAPH_deque(&dirty_priority, num_nodes, &dirty_head);
					vertex     = dirty_node.label;
					
					if (status[vertex] == NUMBERED) continue;
					
					// Updating status
					if (status[vertex] == INACTIVE)
					{
						#pragma omp critical
						{
							if (status[vertex] == INACTIVE)
							{
								status[vertex] = dirty_node.status;
								GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
							}
						}
					}
					else if (dirty_node.status > status[vertex]) 
					{
						#pragma omp critical
						if (dirty_node.status > status[vertex])
							status[vertex] = dirty_node.status;
					}
					
					// Updating priority
					mgraph->sloan_priority[vertex] = dirty_node.distance;
				}
				
				dirty_head = dirty_tail = 0;
			}
			
		} // main loop - while
		
		/* ********************************
		 * ****** Post-processing *********
		 * ********************************
		 */
		
		free(priority_snapshot);
		free(dirty_priority);
		
		#pragma omp barrier
		
		#pragma omp sections
		{
			#pragma omp section
			free(priority_queue);
			
			#pragma omp section
			free(status);
			
			#pragma omp section
			omp_destroy_lock(&queue_lock);
		}
	}
}
