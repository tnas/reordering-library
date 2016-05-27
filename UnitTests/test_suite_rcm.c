#include "test_suite_rcm.h"


void test_prefix_sum()
{
	int* sums;
	int max_level, i;
	
	max_level = 15;
	int counts[15] =       { 9, 8, 3, 2, 7, 1, 6, 4, 5, 5, 8, 3, 4, 9, 1 };
	int correct_sums[15] = { 9, 17, 20, 22, 29, 30, 36, 40, 45, 50, 58, 61, 65, 74, 75 };
	
	prefix_sum(counts, &sums, max_level);
	
	for (i = 0; i < max_level; ++i)
		assert(sums[i] == correct_sums[i]);
}


void test_prefix_sum_parallel_serial(int* counts, int max_level)
{
	int i;
	int* p_sums;
	int* s_sums;
	
	// Parallel prefix sum
	prefix_sum(counts, &p_sums, max_level);
	
	// Serial prefix sum
	s_sums = calloc(max_level, sizeof(int));
	s_sums[0] = counts[0];
	for (i = 1; i < max_level; ++i)
		s_sums[i] = s_sums[i-1] + counts[i];
	
	// Comparing parallel and serial prefix sum
	for (i = 0; i < max_level; ++i) 
	{
		if (s_sums[i] != p_sums[i]) 
		{
			printf("Serial prefix sum not equals Parallel prefix sum: position %d: s = %d and p = %d\n",
				i, s_sums[i], p_sums[i]);
			exit(1);
		}
	}
		
	free(s_sums);
	free(p_sums);
}


int get_node_peripheral(const char* path_matrix_file) {
	
	int root, e;
	int* level_structure;
	MAT* matrix;
	FILE* matrix_file;
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	
	level_structure = GRAPH_LS_peripheral(matrix, &root, &e);
	
	MATRIX_clean(matrix);
	free(level_structure);
	
	return root;
}


test_def test_reorder_algorithm(test_def defs)
{
	long int bandwidth, envelope, bandwidth_after, envelope_after;
	int* permutation;
	double time;
	MAT* matrix;
	FILE* matrix_file;
	int* g;
	int node_s, node_e;
	
	if ((matrix_file = fopen(defs.path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	
	write_output_before(matrix);
	
	bandwidth_after = MATRIX_bandwidth(matrix);
	envelope_after  = MATRIX_envelope(matrix);
	
	if (defs.strategy == parallel)
		omp_set_num_threads(defs.num_threads);
	
	switch (defs.algorithm)
	{
		case serial_rcm :
			time = get_time(); 
			REORDERING_RCM_opt(matrix, &permutation, defs.root);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case hsl_rcm :
			time = get_time(); 
			REORDERING_HSL_RCM(matrix, &permutation);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case hsl_spectral :
			time = get_time(); 
			REORDERING_HSL_SPECTRAL(matrix, &permutation);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case serial_sloan :
			g = GRAPH_LS_peripheral (matrix, &node_s, &node_e);
			free(g);
			
			time = get_time();
			REORDERING_SLOAN(matrix, &permutation, node_s, node_e);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case unordered_rcm :
			time = get_time(); 
			Unordered_RCM(matrix, &permutation, defs.root, defs.percent_chunk);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case leveled_rcm :
			time = get_time(); 
			Leveled_RCM(matrix, &permutation, defs.root);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case bucket_rcm :
			time = get_time(); 
			Bucket_RCM(matrix, &permutation, defs.root);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case parallel_sloan :
			g = GRAPH_LS_peripheral (matrix, &node_s, &node_e);
			free(g);
			
			time = get_time();
			Parallel_Sloan(matrix, &permutation, node_s, node_e);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
	}
	
	MATRIX_permutation(matrix, permutation);
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	defs.bandwidth = bandwidth;
	
	write_output_after(matrix);
	
	free(permutation);
	MATRIX_clean(matrix);
	
	printf("%s: Band/Env [ %ld / %ld => %ld / %ld ] Time [ %.6f ]\n",
		defs.algorithm_name, bandwidth_after, envelope_after, bandwidth, envelope, time); 
	fflush(stdout);
	
	return defs;
}


/* ******************************************************
 * **************** Serial Strategies *******************
 * ******************************************************
 */


test_def test_serial_rcm(const char* path_matrix_file, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.root = root;
	defs.algorithm_name = "Serial RCM";
	defs.algorithm = serial_rcm;
	defs.strategy = serial;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}



test_def test_hsl_spectral(const char* path_matrix_file)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "HSL Spectral";
	defs.algorithm = hsl_spectral;
	defs.strategy = serial;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


test_def test_serial_sloan(const char* path_matrix_file)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Serial Sloan";
	defs.algorithm = serial_sloan;
	defs.strategy = serial;
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


// test_def test_hsl_rcm(const char* path_matrix_file)
// {
// 	long int bandwidth, envelope, bandwidth_after, envelope_after;
// 	int* permutation;
// 	double time;
// 	MAT* matrix;
// 	FILE* matrix_file;
// 	test_def defs;
// 	
// 	defs.algorithm_name = "HSL RCM";
// 	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
// 		exit(1);
// 	
// 	matrix = (MAT*) malloc (sizeof(MAT));
// 	MATRIX_readCSR_SymmUpper (matrix, matrix_file);
// 	fclose(matrix_file);
// 	
// 	MATRIX_printFULL(matrix);
// 	write_output_before(matrix);
// 	
// 	time = get_time(); 
// 	REORDERING_HSL_RCM(matrix, &permutation);
// 	time = (get_time() - time)/100.0;
// 	defs.time = time;
// 	
// 	MATRIX_clean(matrix);
// 	
// 	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
// 		exit(1);
// 	
// 	matrix = (MAT*) malloc (sizeof(MAT));
// 	MATRIX_readCSR (matrix, matrix_file);
// 	fclose(matrix_file);
// 		
// 	bandwidth_after = MATRIX_bandwidth(matrix);
// 	envelope_after  = MATRIX_envelope(matrix);
// 	
// 	MATRIX_permutation(matrix, permutation);
// 	bandwidth = MATRIX_bandwidth(matrix);
// 	envelope  = MATRIX_envelope(matrix);	
// 	defs.bandwidth = bandwidth;
// 	
// 	write_output_after(matrix);
// 	
// 	free(permutation);
// 	MATRIX_clean(matrix);
// 	
// 	printf("%s: Band/Env [ %ld / %ld => %ld / %ld ] Time [ %.6f ]\n",
// 		defs.algorithm_name, bandwidth_after, envelope_after, bandwidth, envelope, time); 
// 	fflush(stdout);
// 	
// 	return defs;
// }


test_def test_hsl_rcm(const char* path_matrix_file)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "HSL RCM";
	defs.algorithm = hsl_rcm;
	defs.strategy = serial;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


/* ******************************************************
 * ************** Parallel Strategies *******************
 * ******************************************************
 */

test_def test_unordered_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_percent, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Unordered RCM";
	defs.algorithm = unordered_rcm;
	defs.root = root;
	defs.percent_chunk = bfs_chunk_percent;
	defs.num_threads = num_threads;
	defs.strategy = parallel;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


test_def test_leveled_rcm(const char* path_matrix_file, const int num_threads, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Leveled RCM";
	defs.algorithm = leveled_rcm;
	defs.root = root;
	defs.num_threads = num_threads;
	defs.strategy = parallel;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


test_def test_bucket_rcm(const char* path_matrix_file, const int num_threads, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Bucket RCM";
	defs.algorithm = bucket_rcm;
	defs.root = root;
	defs.num_threads = num_threads;
	defs.strategy = parallel;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


test_def test_parallel_sloan(const char* path_matrix_file, const int num_threads)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Parallel Sloan";
	defs.algorithm = parallel_sloan;
	defs.strategy = parallel;
	defs.num_threads = num_threads;
	defs = test_reorder_algorithm(defs);
	
	return defs;
}



void normalize_tests(const test_def* results, test_def* result)
{
	int pos_min_band, pos_max_band, pos_min_time, pos_max_time, times;
	double sum_time, max_time, min_time;
	long int sum_band, max_band, min_band;
	
	max_band = max_time = 0;
	min_band = INT_MAX;
	min_time = DBL_MAX;
	sum_band = sum_time = 0;
	
	// Finding out max/min bandwidth and time execution
	for (times = 0; times < TEST_EXEC_TIMES; ++times)
	{
		if (results[times].bandwidth > max_band)
		{
			pos_max_band = times;
			max_band = results[times].bandwidth;
		}
		else if (results[times].bandwidth < min_band)
		{
			pos_min_band = times;
			min_band = results[times].bandwidth;
		}
		
		if (results[times].time > max_time)
		{
			pos_max_time = times;
			max_time = results[times].time;
		}
		else if (results[times].time < min_time)
		{
			pos_min_time = times;
			min_time = results[times].time;
		}
	}
	
	for (times = 0; times < TEST_EXEC_TIMES; ++times)
	{
		if (times != pos_min_band && times != pos_max_band)
			sum_band += results[times].bandwidth;
		
		if (times != pos_max_time && times != pos_min_time)
			sum_time += results[times].time;
	}
	
	result->bandwidth = sum_band / (TEST_EXEC_TIMES - 2);
	result->time      = sum_time / (TEST_EXEC_TIMES - 2);
}


int is_serial_algorithm(reorder_algorithm algorithm)
{
	if (algorithm == serial_rcm || algorithm == serial_sloan ||
		algorithm == hsl_rcm || algorithm == hsl_spectral)
		return 1;
	
	return 0;
}


void run_all_tests()
{
	int root, count_matrix, count_alg, count_exec, count_nthreads, num_matrices, num_nthreads, num_algorithms;
	FILE* out_file;
	test_def* test_results;
	test_def result;
	
	/* *******************************
	 * Definition of tests parameters
	 * *******************************
	 */
	float bfs_chunk_percent = .5;
	
	char* matrices[] = {
// 		"../Big-Matrices/dw8192.mtx",
// 		"../Big-Matrices/rail_79841.mtx",
// 		"../Big-Matrices/Dubcova3.mtx",
// 		"../Big-Matrices/inline_1.mtx",
// 		"../Big-Matrices/audikw_1.mtx",
// 		"../Big-Matrices/dielFilterV3real.mtx",
// 		"../Big-Matrices/atmosmodj.mtx",
// 		"../Big-Matrices/G3_circuit.mtx"
// 		"../Matrices/rail_5177.mtx",
// 		"../Matrices/bcspwr01.mtx",
// 		"../Matrices/bcspwr02.mtx",
// 		"../Matrices/rail_5177.mtx"
		"../Matrices/FEM_3D_thermal1.mtx",
// 		"../Matrices/Dubcova2.mtx"
	};
	
	int nthreads[] = { 4, 8, 16, 32, 64, 128 };
	
	reorder_algorithm algorithm[] = { serial_rcm };
	
	/* *****************
	 * Tests execution
	 * *****************
	 */
	
	if ((out_file = fopen("run_all_tests_normalized_output.txt", "w")) == NULL) 
		exit(1);
	
	num_matrices   = sizeof(matrices)/sizeof(matrices[0]);
	num_nthreads   = sizeof(nthreads)/sizeof(nthreads[0]);
	num_algorithms = sizeof(algorithm)/sizeof(algorithm[0]);

	for (count_matrix = 0; count_matrix < num_matrices; ++count_matrix)
	{
		root = get_node_peripheral(matrices[count_matrix]);
		
		fprintf(out_file, "-----------------------------------------------------------------------\n");
		fprintf(out_file, "Tests Execution - Matrix: %s\n", matrices[count_matrix]);
		fflush(out_file);
		
		for (count_alg = 0; count_alg < num_algorithms; ++count_alg)
		{
			fprintf(out_file, "-----------------------------------------------------------------------\n");
			fprintf(out_file, "Algorithm: %d\n", algorithm[count_alg]);
			fprintf(out_file, "-----------------------------------------------------------------------\n");
			fflush(out_file);
			
// 			count_nthreads = is_serial_algorithm(algorithm[count_alg]) ?
// 				num_nthreads - 1 : 0;
			
			for (count_nthreads = 0; count_nthreads < num_nthreads; ++count_nthreads)
			{
				test_results = calloc(TEST_EXEC_TIMES, sizeof(test_def));
				
				for (count_exec = 0; count_exec < TEST_EXEC_TIMES; ++count_exec)
				{
					switch (algorithm[count_alg])
					{
						case unordered_rcm :
							test_results[count_exec] = test_unordered_rcm(matrices[count_matrix], nthreads[count_nthreads], bfs_chunk_percent, root);
							break;
							
						case serial_rcm :
							test_results[count_exec] = test_serial_rcm(matrices[count_matrix], root);
							break;
							
						case leveled_rcm :
							test_results[count_exec] = test_leveled_rcm(matrices[count_matrix], nthreads[count_nthreads], root);
							break;
							
						case bucket_rcm :
							test_results[count_exec] = test_bucket_rcm(matrices[count_matrix], nthreads[count_nthreads], root);
							break;
							
						case serial_sloan :
							test_results[count_exec] = test_serial_sloan(matrices[count_matrix]);
							break;
							
						case hsl_spectral :
							test_results[count_exec] = test_hsl_spectral(matrices[count_matrix]);
							break;
							
						case parallel_sloan :
							test_results[count_exec] = test_parallel_sloan(matrices[count_matrix], nthreads[count_nthreads]);
							break;
							
						default :
							break;
						
					}
				}
				
				normalize_tests(test_results, &result);
				fprintf(out_file, "[%s] Threads: %d -- Bandwidth: %ld -- Time: %.6f\n", 
					test_results[0].algorithm_name, nthreads[count_nthreads], result.bandwidth, result.time);
				fflush(out_file);
				
				free(test_results);
			}
		}
		
	}
	
	fclose(out_file);
}

