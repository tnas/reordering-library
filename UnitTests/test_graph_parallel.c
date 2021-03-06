/*
 * Copyright 2016 Thiago Nascimento <nascimenthiago@gmail.com>
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

#include "test_graph_parallel.h"


void test_GRAPH_parallel_build_METAGRAPH()
{
	METAGRAPH* mgraph;
	MAT* matrix;
	char* matrix_path = "../Matrices/sample.mtx";
	int expected_graph_size = 10;
	int expected_node_min_degree = 1;
	
	MATRIX_read_from_path(matrix_path, &matrix);
	
	mgraph = GRAPH_parallel_build_METAGRAPH(matrix);
	
	assert(mgraph->size == expected_graph_size);
	assert(mgraph->vertex_min_degree == expected_node_min_degree);
	
	assert(mgraph->graph[0].degree == 2);
	assert(mgraph->graph[1].degree == 1);
	assert(mgraph->graph[2].degree == 2);
	assert(mgraph->graph[3].degree == 2);
	assert(mgraph->graph[4].degree == 2);
	assert(mgraph->graph[5].degree == 2);
	assert(mgraph->graph[6].degree == 2);
	assert(mgraph->graph[7].degree == 4);
	assert(mgraph->graph[8].degree == 3);
	assert(mgraph->graph[9].degree == 2);
	
	GRAPH_parallel_destroy_METAGRAPH(mgraph);
	
	printf("test_GRAPH_parallel_build_METAGRAPH -- SUCCESS\n");fflush(stdout);
}


void test_GRAPH_shrinking_strategy_half_sorted()
{
	int length = 5;
	GRAPH* nodes;
	GRAPH* shrinked_nodes;
	
	nodes = calloc(length, sizeof(GRAPH));
	
	nodes[0].label = 0;
	nodes[0].degree = 5;
	
	nodes[1].label = 1;
	nodes[1].degree = 4;
	
	nodes[2].label = 2;
	nodes[2].degree = 3;
	
	nodes[3].label = 3;
	nodes[3].degree = 2;
	
	nodes[4].label = 4;
	nodes[4].degree = 1;
	
	shrinked_nodes = GRAPH_shrinking_strategy_half_sorted(nodes, &length);
	
	assert(length == 2);
	
	assert(shrinked_nodes[0].label = 4);
	assert(shrinked_nodes[1].label = 3);
	
	assert(shrinked_nodes[0].degree = 1);
	assert(shrinked_nodes[1].degree = 2);
	
	free(nodes);
	free(shrinked_nodes);
	
	printf("test_GRAPH_shrinking_strategy_half_sorted -- SUCCESS\n");fflush(stdout);
}


void test_GRAPH_shrinking_strategy_half_sorted_can24()
{
	int node, length;
	int half_length = 12;
	GRAPH* nodes;
	GRAPH* shrinked_nodes;
	GRAPH expected_shrink_nodes[half_length];
	MAT* matrix;
	char* matrix_path = "../Matrices/can24.mtx";
	
	expected_shrink_nodes[0].label  = 8;    expected_shrink_nodes[6].label  = 3;
        expected_shrink_nodes[0].degree = 3;    expected_shrink_nodes[6].degree = 5;
	
	expected_shrink_nodes[1].label  = 16;   expected_shrink_nodes[7].label  = 4;
	expected_shrink_nodes[1].degree = 3;	expected_shrink_nodes[7].degree = 5;
	
	expected_shrink_nodes[2].label  = 22;   expected_shrink_nodes[8].label  = 5;
	expected_shrink_nodes[2].degree = 3;	expected_shrink_nodes[8].degree = 5;
	
	expected_shrink_nodes[3].label  = 23;   expected_shrink_nodes[9].label  = 10;
	expected_shrink_nodes[3].degree = 3;	expected_shrink_nodes[9].degree = 5;
	
	expected_shrink_nodes[4].label  = 1;    expected_shrink_nodes[10].label  = 11;
	expected_shrink_nodes[4].degree = 5;	expected_shrink_nodes[10].degree = 5;
	
	expected_shrink_nodes[5].label  = 2;    expected_shrink_nodes[11].label  = 12;
	expected_shrink_nodes[5].degree = 5;	expected_shrink_nodes[11].degree = 5;
	
	MATRIX_read_from_path(matrix_path, &matrix);
	length = matrix->n;
	nodes = calloc(length, sizeof(GRAPH));
	
	for (node = 0; node < length; node++)
	{
		nodes[node].label = node;
		nodes[node].degree = GRAPH_degree(matrix, node);
	}
	
	shrinked_nodes = GRAPH_shrinking_strategy_half_sorted(nodes, &length);
	
	assert(length == half_length);
	
	for (node = 0; node < half_length; node++)
	{
		assert(shrinked_nodes[node].label == expected_shrink_nodes[node].label);
		assert(shrinked_nodes[node].degree == expected_shrink_nodes[node].degree);
	}
	
	free(nodes);
	free(shrinked_nodes);
	MATRIX_clean(matrix);
	
	printf("test_GRAPH_shrinking_strategy_half_sorted_can24 -- SUCCESS\n");fflush(stdout);
}


void test_GRAPH_shrinking_strategy_vertex_by_degree()
{
	int length = 5;
	GRAPH* nodes;
	GRAPH* shrinked_nodes;
	
	nodes = calloc(length, sizeof(GRAPH));
	
	nodes[0].label = 0;
	nodes[0].degree = 5;
	
	nodes[1].label = 1;
	nodes[1].degree = 3;
	
	nodes[2].label = 2;
	nodes[2].degree = 5;
	
	nodes[3].label = 3;
	nodes[3].degree = 2;
	
	nodes[4].label = 4;
	nodes[4].degree = 3;
	
	shrinked_nodes = GRAPH_shrinking_strategy_vertex_by_degree(nodes, &length);
	
	assert(length == 3);
	
	assert(shrinked_nodes[0].label = 3);
	assert(shrinked_nodes[0].degree = 2);
	
	assert(shrinked_nodes[1].label = 1);
	assert(shrinked_nodes[1].degree = 3);
	
	assert(shrinked_nodes[2].label = 2);
	assert(shrinked_nodes[2].degree = 5);
	
	free(nodes);
	free(shrinked_nodes);
	
	printf("test_GRAPH_shrinking_strategy_vertex_by_degree -- SUCCESS\n");fflush(stdout);
}


void test_GRAPH_shrinking_strategy_five_non_adjacent()
{
	int length = 7;
	GRAPH* nodes;
	GRAPH* shrinked_nodes;
	
	nodes = calloc(length, sizeof(GRAPH));
	
	nodes[0].label = 0;
	nodes[0].degree = 5;
	int neighs0[] = {51, 52, 53, 54, 55};
	nodes[0].neighboors = neighs0;
	
	nodes[1].label = 1;
	nodes[1].degree = 4;
	int neighs1[] = {41, 42, 43, 44};
	nodes[1].neighboors = neighs1;
	
	// Should be discarded - adjacent with label 3
	nodes[2].label = 2;
	nodes[2].degree = 3;
	int neighs2[] = {31, 22, 32}; 
	nodes[2].neighboors = neighs2;
	// -------------------------------
	
	nodes[3].label = 3;
	nodes[3].degree = 2;
	int neighs3[] = {21, 22};
	nodes[3].neighboors = neighs3;
	
	nodes[4].label = 4;
	nodes[4].degree = 1;
	int neighs4[] = {11};
	nodes[4].neighboors = neighs4;
	
	nodes[5].label = 5;
	nodes[5].degree = 6;
	int neighs5[] = {61, 62, 63, 64, 65, 66};
	nodes[5].neighboors = neighs5;
	
	nodes[6].label = 6;
	nodes[6].degree = 7;
	int neighs6[] = {71, 72, 73, 74, 75, 76, 77};
	nodes[6].neighboors = neighs6;
	
	shrinked_nodes = GRAPH_shrinking_strategy_five_non_adjacent(nodes, &length);
	
	assert(length == 5);
	
	assert(shrinked_nodes[0].label = 4);
	assert(shrinked_nodes[0].degree = 1);
	
	assert(shrinked_nodes[1].label = 3);
	assert(shrinked_nodes[1].degree = 2);
	
	assert(shrinked_nodes[2].label = 1);
	assert(shrinked_nodes[2].degree = 4);
	
	assert(shrinked_nodes[3].label = 5);
	assert(shrinked_nodes[3].degree = 6);
	
	assert(shrinked_nodes[4].label = 6);
	assert(shrinked_nodes[4].degree = 7);
	
	free(nodes);
	free(shrinked_nodes);
	
	printf("test_GRAPH_shrinking_strategy_five_non_adjacent -- SUCCESS\n");fflush(stdout);
}


void test_GRAPH_parallel_build_BFS()
{
	BFS* bfs;
	METAGRAPH* mgraph;
	MAT* matrix;
	char* matrix_path = "../Matrices/sample.mtx";
	int root = 1;
	int expected_height = 5;
	
	MATRIX_read_from_path(matrix_path, &matrix);
	mgraph = GRAPH_parallel_build_METAGRAPH(matrix);
	
	bfs = GRAPH_parallel_build_BFS(mgraph, root);
	
	assert(bfs->height == expected_height);
	
	GRAPH_parallel_destroy_BFS(bfs);
	GRAPH_parallel_destroy_METAGRAPH(mgraph);
	
	printf("test_GRAPH_parallel_build_BFS -- SUCCESS\n");fflush(stdout);
}


void test_GRAPH_parallel_pseudodiameter_sample()
{
	METAGRAPH* mgraph;
	graph_diameter* diameter;
	MAT* matrix;
	char* matrix_path = "../Matrices/sample.mtx";
	
	MATRIX_read_from_path(matrix_path, &matrix);
	mgraph = GRAPH_parallel_build_METAGRAPH(matrix);
	
	diameter = GRAPH_parallel_pseudodiameter(mgraph, HALF_SORTED);
	
	GRAPH_parallel_destroy_METAGRAPH(mgraph);
	free(diameter);
	
	printf("test_GRAPH_parallel_pseudodiameter_sample -- SUCCESS\n");fflush(stdout);
}


void run_all_test_GRAPH_parallel()
{
	printf("[UNIT TESTS]\n");fflush(stdout);
	
	omp_set_num_threads(1);
	printf("Tests with a single thread\n");fflush(stdout);
	test_GRAPH_parallel_build_METAGRAPH();
	test_GRAPH_shrinking_strategy_half_sorted();
	test_GRAPH_shrinking_strategy_half_sorted_can24();
	test_GRAPH_shrinking_strategy_vertex_by_degree();
	test_GRAPH_shrinking_strategy_five_non_adjacent();
	test_GRAPH_parallel_build_BFS();
	test_GRAPH_parallel_pseudodiameter_sample();
	
	printf("Tests with %d threads\n", NUM_THREADS);fflush(stdout);
	omp_set_num_threads(NUM_THREADS);
	test_GRAPH_parallel_build_METAGRAPH();
	test_GRAPH_shrinking_strategy_half_sorted();
	test_GRAPH_shrinking_strategy_half_sorted_can24();
	test_GRAPH_shrinking_strategy_vertex_by_degree();
	test_GRAPH_shrinking_strategy_five_non_adjacent();
	test_GRAPH_parallel_build_BFS();
}


void compare_hsl_pseudodiameter()
{
	METAGRAPH* mgraph;
	graph_diameter* diameter;
	MAT* matrix;
	int* peripherals;
	double start_time_hsl, end_time_hsl, start_time_parallel, end_time_parallel;
	int mat, strategy;
	
	char* matrix_path[] = {
		"../Matrices/sample.mtx",
		"../Matrices/can24.mtx",
		"../Matrices/rail_5177.mtx",
		"../Matrices/FEM_3D_thermal1.mtx",
	};
	
	int num_matrices = sizeof(matrix_path)/sizeof(matrix_path[0]);
	
	for (mat = 0; mat < num_matrices; mat++)
	{
		// Executing HSL pseudo-diameter
		MATRIX_read_from_path(matrix_path[mat], &matrix);
		start_time_hsl = omp_get_wtime();
		peripherals = get_pseudo_diameter_hsl(matrix);
		end_time_hsl = omp_get_wtime();
		MATRIX_clean(matrix);
		
		printf("[%s] Pseudo-diameter HSL (%d, %d): %.5f\n", matrix_path[mat],
			peripherals[0], peripherals[1], (end_time_hsl - start_time_hsl)/100);fflush(stdout);
		
		// Executing Parallel pseudo-diameter
		MATRIX_read_from_path(matrix_path[mat], &matrix);
		mgraph   = GRAPH_parallel_build_METAGRAPH(matrix);
		
		for (strategy = 0; strategy < 3; strategy++)
		{
			start_time_parallel = omp_get_wtime();
			diameter = GRAPH_parallel_pseudodiameter(mgraph, strategy);
			end_time_parallel = omp_get_wtime();
			
			printf("[%s] Pseudo-diameter Parallel - Strategy %d (%d, %d): %.5f\n", 
			matrix_path[mat], strategy, diameter->start, diameter->end, 
			(end_time_parallel - start_time_parallel)/100);fflush(stdout);
		}
		
		MATRIX_clean(matrix);
		
		
	}
}


void run_all_comparison_GRAPH_parallel()
{
	printf("[PERFORMANCE TESTS]\n");fflush(stdout);
	
	omp_set_num_threads(1);
	printf("Tests with a single thread\n");fflush(stdout);
	compare_hsl_pseudodiameter();
	
	omp_set_num_threads(NUM_THREADS);
	printf("Tests with %d threads\n", NUM_THREADS);fflush(stdout);
	compare_hsl_pseudodiameter();
}