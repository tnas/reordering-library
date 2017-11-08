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

void Parallel_Logical_Bag_Sloan_new(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, max_priority, num_prior_bags, 
	    num_threads, chunk_size, count_threads_on;
	int** priority;
	int* size_bags;
	
	num_nodes = mgraph->size;
	GRAPH_parallel_fixedpoint_sloan_BFS(mgraph, start_node, BFS_PERCENT_CHUNK);
	
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
		
		#pragma omp barrier
		
		#pragma omp single
		{
			num_prior_bags = SLOAN_PRIORITY_FACTOR * 
				(mgraph->graph[start_node].priority - mgraph->min_sloan_priority);
			size_bags = calloc(num_prior_bags, sizeof(LIST*));
		}
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node) 
		{
			priority[node] = calloc(2, sizeof(int));
			priority[node][SLOAN_CURR_PRIOR] = mgraph->graph[node].priority - mgraph->min_sloan_priority;
			priority[node][SLOAN_NEW_PRIOR]  = priority[node][SLOAN_CURR_PRIOR];
			printf("node %d priority %d\n", node, priority[node][SLOAN_CURR_PRIOR]);fflush(stdout);
		}
		
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
				if ((mgraph->graph[vertex].status == ACTIVE || mgraph->graph[vertex].status == PREACTIVE) && (priority[vertex][SLOAN_CURR_PRIOR] == max_priority))
				{
					vertex_degree = mgraph->graph[vertex].degree;
					neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
					
					#pragma omp atomic
					count_threads_on++;
					
					for (ngb = 0; ngb < vertex_degree; ++ngb)
					{
						update_far = OFF;
						neighbor   = neighbors[ngb];
						
						if (mgraph->graph[vertex].status == PREACTIVE && 
							(mgraph->graph[neighbor].status == INACTIVE || mgraph->graph[neighbor].status == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
							
							#pragma omp critical
							mgraph->graph[neighbor].status = ACTIVE;
							
							update_far = ON;
						}
						else if (mgraph->graph[vertex].status == PREACTIVE && mgraph->graph[neighbor].status == ACTIVE)
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
						}
						else if ((mgraph->graph[vertex].status == ACTIVE) && (mgraph->graph[neighbor].status == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
							
							#pragma omp critical
							mgraph->graph[neighbor].status = ACTIVE;
							
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
									if (mgraph->graph[far_neighbor].status == INACTIVE)
									{
										#pragma omp critical
										mgraph->graph[far_neighbor].status = PREACTIVE;
									}
									
									#pragma omp atomic
									priority[far_neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
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
						mgraph->graph[vertex].status = NUMBERED;
					}
				}
			}
			
			#pragma omp for schedule(static, chunk_size)
			for (node = 0; node < num_nodes; ++node)
			{
				if ((mgraph->graph[node].status != NUMBERED) && 
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
		}
		
		#pragma omp barrier
		
		#pragma omp single nowait
		free(size_bags);
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
			free(priority[node]);
		
		#pragma omp single nowait
		free(priority);
	}
}


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
			priority[node][SLOAN_CURR_PRIOR] = BAG_SLOAN_W1*distance[node] - BAG_SLOAN_W2*(degree[node] + 1);
				
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
				if ((status[vertex] == ACTIVE || status[vertex] == PREACTIVE) && (priority[vertex][SLOAN_CURR_PRIOR] == max_priority))
				{
					vertex_degree = mgraph->graph[vertex].degree;
					neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
					
					#pragma omp atomic
					count_threads_on++;
					
					for (ngb = 0; ngb < vertex_degree; ++ngb)
					{
						update_far = OFF;
						neighbor   = neighbors[ngb];
						
						if (status[vertex] == PREACTIVE && (status[neighbor] == INACTIVE || status[neighbor] == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
							
							#pragma omp critical
							status[neighbor] = ACTIVE;
							
							update_far = ON;
						}
						else if (status[vertex] == PREACTIVE && status[neighbor] == ACTIVE)
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
						}
						else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
							
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
									priority[far_neighbor][SLOAN_NEW_PRIOR] += BAG_SLOAN_W2;
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


int COMPARE_priority_DESC(const void* st, const void* nd)
{ 
	volatile SLOAN_GRAPH* g_st = (SLOAN_GRAPH*) st;
	volatile SLOAN_GRAPH* g_nd = (SLOAN_GRAPH*) nd;
		       
	if (g_st->priority > g_nd->priority) 
	{
		return -1;
	}
	if (g_st->priority < g_nd->priority)
	{
		return 1;
	}
	
	return 0;
}


void Parallel_Relaxed_Order_Sloan_current(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, prior_head, prior_tail, pqueue_size, norm;
	SLOAN_GRAPH* priority_queue;
	
	num_nodes   = mgraph->size;
	pqueue_size = (omp_get_max_threads() + 10) * num_nodes; // oversizing estimate
	GRAPH_parallel_fixedpoint_sloan_BFS(mgraph, start_node, BFS_PERCENT_CHUNK);
	
	norm = floor((mgraph->graph[end_node].distance - mgraph->graph[start_node].distance) / 
		mgraph->max_degree);
	
	#pragma omp parallel 
	{
		int vertex, vertex_degree, neighbor_degree, ngb, far_ngb, neighbor, 
		    far_neighbor, update_far, th_head, th_tail, size_chunk, 
		    count_chunk, dirty_head, dirty_tail;
		SLOAN_GRAPH dirty_node;
		int* neighbors;
		int* far_neighbors;
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
			priority_queue = calloc(pqueue_size, sizeof(SLOAN_GRAPH));
			
			#pragma omp section
			next_id = prior_head = prior_tail = 0;
		}
		
		#pragma omp for
		for (vertex = 0; vertex < num_nodes; ++vertex)
		{
			mgraph->graph[vertex].status = INACTIVE;
			mgraph->graph[vertex].priority = 
				(-norm)*SLOAN_W1*(mgraph->graph[end_node].distance-mgraph->graph[vertex].distance) +
				SLOAN_W2*(mgraph->graph[vertex].degree + 1);
		}
		
		#pragma omp single nowait
		mgraph->graph[start_node].status = PREACTIVE;
		
		#pragma omp single
		{
			dirty_node.label      = start_node;
			dirty_node.status     = PREACTIVE;
			GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
		}
		
		dirty_priority = calloc(pqueue_size, sizeof(SLOAN_GRAPH));
		dirty_head = dirty_tail = 0;
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while ((!QUEUE_empty(priority_queue, prior_head, prior_tail)) || next_id < num_nodes)
		{
			if (!QUEUE_empty(priority_queue, prior_head, prior_tail)) 
			{
				#pragma omp critical
				{
					th_tail = prior_tail;
					th_head = prior_head;
					qsort(&priority_queue[th_head], th_tail-th_head, sizeof(SLOAN_GRAPH), COMPARE_priority_DESC);
					size_chunk = ceil(BFS_PERCENT_CHUNK * (th_tail - th_head));
					prior_head += size_chunk;
				}
			}
			
			// Processing snapshot of maximum priority nodes
			for (count_chunk = 0; count_chunk < size_chunk; ++count_chunk)
			{
				vertex        = GRAPH_deque(&priority_queue, pqueue_size, &th_head).label;
				vertex_degree = mgraph->graph[vertex].degree;
				neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor   = neighbors[ngb];
					update_far = OFF;
					
					if (mgraph->graph[vertex].status == PREACTIVE && 
						(mgraph->graph[neighbor].status == INACTIVE || mgraph->graph[neighbor].status == PREACTIVE))
					{
						dirty_node.label    = neighbor;
						dirty_node.status   = ACTIVE;
						GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
						
						update_far = ON;
					}
					else if (mgraph->graph[vertex].status == PREACTIVE && mgraph->graph[neighbor].status == ACTIVE)
					{
						dirty_node.label    = neighbor;
						dirty_node.status   = ACTIVE;
						GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
					}
					else if (mgraph->graph[vertex].status == ACTIVE && mgraph->graph[neighbor].status == PREACTIVE)
					{
						dirty_node.label    = neighbor;
						dirty_node.status   = ACTIVE;
						GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
						
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
							
							if (far_neighbor == vertex) continue;
							
							dirty_node.label = far_neighbor;
							dirty_node.status = mgraph->graph[far_neighbor].status == INACTIVE ? 
								PREACTIVE : mgraph->graph[far_neighbor].status;
							GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
						}
						
						free(far_neighbors);
					}
				}
				
				free(neighbors);
				
				// Placing vertex in permutation array
				if (mgraph->graph[vertex].status != NUMBERED)
				{
					#pragma omp critical
					{
						if (mgraph->graph[vertex].status != NUMBERED)
						{
							(*permutation)[next_id++] = vertex;
							mgraph->graph[vertex].status = NUMBERED;
						}
					}
				}
				
			} // priority snapshot loop
			
			// Updating priorities and status from dirty nodes
			if (!QUEUE_empty(dirty_priority, dirty_head, dirty_tail))
			{
				while (!QUEUE_empty(dirty_priority, dirty_head, dirty_tail))
				{
					dirty_node = GRAPH_deque(&dirty_priority, pqueue_size, &dirty_head);
					vertex     = dirty_node.label;
					
					if (mgraph->graph[vertex].status == NUMBERED) continue;
					
					// Updating priority
					mgraph->graph[vertex].priority += SLOAN_W1;
					
					// Updating status
					if (dirty_node.status > mgraph->graph[vertex].status) 
					{
						#pragma omp critical
						{
							if (dirty_node.status > mgraph->graph[vertex].status)
							{
								mgraph->graph[vertex].status = dirty_node.status;
								GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
							}
						}
					}
// 					if (dirty_node.status > mgraph->graph[vertex].status)
// 					{
// 						mgraph->graph[vertex].status = dirty_node.status;
// 						
// 						#pragma omp critical
// 						GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
// 					}
				}
				
				dirty_tail = dirty_head = 0;
			}
			
		} // main loop - while
		
		/* ********************************
		 * ****** Post-processing *********
		 * ********************************
		 */
		free(dirty_priority);
		
		#pragma omp barrier
		
		#pragma omp single
		free(priority_queue);
	}
}



int get_current_degree(int node, int* adjacents, METAGRAPH* mgraph)
{
	int n_adj, curr_dgr, const_k;
	curr_dgr = 0;
	const_k  = 1;
	
	for (n_adj = 0; n_adj < mgraph->graph[node].degree; n_adj++)
		if (mgraph->graph[adjacents[n_adj]].status == NUMBERED ||
			mgraph->graph[adjacents[n_adj]].status == ACTIVE) ++curr_dgr;
	
	if (mgraph->graph[node].status == NUMBERED || mgraph->graph[node].status == ACTIVE)	
		const_k = 0;
		
	return mgraph->graph[node].degree - curr_dgr + const_k;
}


void Parallel_Relaxed_Order_Sloan(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, prior_head, prior_tail, pqueue_size;
	SLOAN_GRAPH* priority_queue;
	
	num_nodes   = mgraph->size;
	pqueue_size = (omp_get_max_threads() + 10) * num_nodes; // oversizing estimate
	GRAPH_parallel_fixedpoint_sloan_BFS(mgraph, end_node, BFS_PERCENT_CHUNK);
	
	
// 	int x;
// 	for (x=0; x<mgraph->size; ++x)
// 	{
// 		printf("d[%d] = %d, ", x, mgraph->graph[x].distance);
// 	}
// 	printf("\n");fflush(stdout);
// 	
// 	printf("start_node: %d, end_node: %d\n", start_node, end_node);fflush(stdout);
	
// 	printf("end_distance: %d, start_distance: %d, max_degree: %d\n", 
// 	       mgraph->graph[end_node].distance, mgraph->graph[start_node].distance, mgraph->max_degree);fflush(stdout);
// 	
// 	norm = floor((mgraph->graph[end_node].distance - mgraph->graph[start_node].distance) / 
// 		mgraph->max_degree);
// 	printf("start_node: %d, end_node: %d, norm: %d\n", start_node, end_node, norm);fflush(stdout);
	
	#pragma omp parallel 
	{
		int vertex, vertex_degree, neighbor_degree, ngb, far_ngb, neighbor, 
		    far_neighbor, th_head, th_tail, size_chunk, 
		    count_chunk, dirty_head, dirty_tail;
		SLOAN_GRAPH dirty_node;
		int* neighbors;
		int* far_neighbors;
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
			priority_queue = calloc(pqueue_size, sizeof(SLOAN_GRAPH));
			
			#pragma omp section
			next_id = prior_head = prior_tail = 0;
		}
		
		#pragma omp for
		for (vertex = 0; vertex < num_nodes; ++vertex)
		{
			mgraph->graph[vertex].status = INACTIVE;
			mgraph->graph[vertex].priority = 
			
				SLOAN_W2*mgraph->graph[vertex].distance - SLOAN_W1*(mgraph->graph[vertex].degree + 1);
			
// 				SLOAN_W1*(mgraph->size - get_current_degree(vertex, adjacents, mgraph)) + mgraph->graph[vertex].distance * SLOAN_W2; 
				
// 				SLOAN_W1*(mgraph->size - mgraph->graph[vertex].degree) + mgraph->graph[vertex].distance * SLOAN_W2; 
			
// 				(-norm)*SLOAN_W1*(mgraph->graph[end_node].distance-mgraph->graph[vertex].distance) +
// 				SLOAN_W2*(mgraph->graph[vertex].degree + 1);

// 			printf("node: %d, distance: %d, priority:%d\n", vertex, mgraph->graph[vertex].distance, mgraph->graph[vertex].priority);fflush(stdout);
		}
		
		#pragma omp single nowait
		mgraph->graph[start_node].status = PREACTIVE;
		
		#pragma omp single nowait
		{
			dirty_node.label      = start_node;
			dirty_node.status     = PREACTIVE;
			dirty_node.priority   = mgraph->graph[start_node].priority;
			GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
		}
		
		dirty_priority = calloc(pqueue_size, sizeof(SLOAN_GRAPH));
		dirty_head = dirty_tail = 0;
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while ((!QUEUE_empty(priority_queue, prior_head, prior_tail)) || next_id < num_nodes)
		{
			if (!QUEUE_empty(priority_queue, prior_head, prior_tail)) 
			{
				#pragma omp critical
				{
					th_tail = prior_tail;
					th_head = prior_head;
					
					size_chunk = ceil(BFS_PERCENT_CHUNK * (th_tail - th_head));
					
// 					printf("th_tail: %d, th_head: %d, size_chunk: %d\n", th_tail, th_head, size_chunk);fflush(stdout);
					
// 					printf("Before sorting: ");
// 					for (ngb = th_head; ngb < th_tail; ++ngb)
// 						printf("%d(%d) - ", priority_queue[ngb].label, priority_queue[ngb].priority);
// 					printf("\n");fflush(stdout);
					
					// Nodes sorting does not impact on the reordering quality???
// 					qsort(&priority_queue[th_head], th_tail-th_head, sizeof(SLOAN_GRAPH), COMPARE_priority_DESC); 
					
					prior_head += size_chunk;
					
// 					printf("After sorting: ");
// 					for (ngb = th_head; ngb < th_tail; ++ngb)
// 						printf("%d(p:%d/s:%d) - ", priority_queue[ngb].label, priority_queue[ngb].priority, mgraph->graph[priority_queue[ngb].label].status);
// 					printf("\n-----------------------------------\n");fflush(stdout);
				}
			}
			
			// Processing snapshot of maximum priority nodes
			for (count_chunk = 0; count_chunk < size_chunk; ++count_chunk)
			{
				vertex = GRAPH_deque(&priority_queue, pqueue_size, &th_head).label;
				
				if (mgraph->graph[vertex].status == PREACTIVE)
				{
					vertex_degree = mgraph->graph[vertex].degree;
					neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
					
					for (ngb = 0; ngb < vertex_degree; ++ngb)
					{
						neighbor = neighbors[ngb];
						
						dirty_node.label = neighbor;
						
						if (mgraph->graph[neighbor].status == INACTIVE)
						{
							dirty_node.status = PREACTIVE;
// 							GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
							
						}	
						else if (mgraph->graph[neighbor].status == PREACTIVE) 
						{
							dirty_node.status = ACTIVE;
// 							GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
						}
						
						GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
						
						if (dirty_node.status == PREACTIVE)
						{
							/* ***********************
							* Updating far neighbors
							* ***********************
							*/
							dirty_node.status = ACTIVE;
							
							neighbor_degree = mgraph->graph[neighbor].degree;
							far_neighbors   = GRAPH_neighboors(mgraph->mat, neighbor, neighbor_degree);
							
							for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
							{
								far_neighbor = far_neighbors[far_ngb];
								
								if (far_neighbor == vertex) continue;
								
								dirty_node.label = far_neighbor;
								
								if (mgraph->graph[far_neighbor].status == PREACTIVE || 
									mgraph->graph[far_neighbor].status == ACTIVE)
								{
									GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
								}
								else if (mgraph->graph[far_neighbor].status == INACTIVE)
								{
									dirty_node.status = PREACTIVE;
									GRAPH_enque(&dirty_priority, pqueue_size, &dirty_tail, dirty_node);
								}
							}
							
							free(far_neighbors);
						}
					}
					
					free(neighbors);
				}
				
				// Placing vertex in permutation array
				if (mgraph->graph[vertex].status != NUMBERED)
				{
					#pragma omp critical
					{
						if (mgraph->graph[vertex].status != NUMBERED)
						{
							(*permutation)[next_id++] = vertex;
							mgraph->graph[vertex].status = NUMBERED;
// 							printf("Node numbered: %d\n", vertex);fflush(stdout);
						}
					}
				}
				
			} // priority snapshot loop
			
			// Updating priorities and status from dirty nodes
			while (!QUEUE_empty(dirty_priority, dirty_head, dirty_tail))
			{
				dirty_node = GRAPH_deque(&dirty_priority, pqueue_size, &dirty_head);
				vertex     = dirty_node.label;
				
				if (mgraph->graph[vertex].status == NUMBERED) continue;
				
				// Updating priority
				mgraph->graph[vertex].priority += SLOAN_W1;
				dirty_node.priority = mgraph->graph[vertex].priority;
				
				// Updating status and adding to priority queue
				if (dirty_node.status == PREACTIVE) 
				{
					#pragma omp critical
					{
						if (dirty_node.status > mgraph->graph[vertex].status)
						{
							mgraph->graph[vertex].status = dirty_node.status;
							GRAPH_enque(&priority_queue, pqueue_size, &prior_tail, dirty_node);
						}
					}
				}
			}
			
			dirty_tail = dirty_head = 0;
			
		} // main loop - while
		
		/* ********************************
		 * ****** Post-processing *********
		 * ********************************
		 */
		free(dirty_priority);
		
		#pragma omp barrier
		
		#pragma omp single
		free(priority_queue);
	}
}


