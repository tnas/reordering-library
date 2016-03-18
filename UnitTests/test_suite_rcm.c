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
		case rcm :
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
			
		case leveled_rcm :
			time = get_time(); 
			Leveled_RCM(matrix, &permutation, defs.root);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case leveled_rcm_v1 :
			time = get_time(); 
			Leveled_RCM_v1(matrix, &permutation, defs.root);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case leveled_rcm_v2 :
			time = get_time(); 
			Leveled_RCM_v2(matrix, &permutation, defs.root);
			time = (get_time() - time)/100.0;
			defs.time = time;
			break;
			
		case unordered_rcm :
			time = get_time(); 
			Unordered_RCM(matrix, &permutation, defs.root, defs.percent_chunk);
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



test_def test_serial_rcm(const char* path_matrix_file, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.root = root;
	defs.algorithm_name = "Serial RCM";
	defs.algorithm = rcm;
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




test_def test_leveled_rcm_v1(const char* path_matrix_file, const int num_threads, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Leveled RCM v1";
	defs.algorithm = leveled_rcm_v1;
	defs.root = root;
	defs.num_threads = num_threads;
	defs.strategy = parallel;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}



test_def test_leveled_rcm_v2(const char* path_matrix_file, const int num_threads, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Leveled RCM v2";
	defs.algorithm = leveled_rcm_v2;
	defs.root = root;
	defs.num_threads = num_threads;
	defs.strategy = parallel;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}


test_def test_unordered_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_size, int root)
{
	test_def defs;
	
	defs.path_matrix_file = path_matrix_file;
	defs.algorithm_name = "Unordered RCM";
	defs.algorithm = unordered_rcm;
	defs.root = root;
	defs.percent_chunk = .5;
	defs.num_threads = num_threads;
	defs.strategy = parallel;
	
	defs = test_reorder_algorithm(defs);
	
	return defs;
}



void print_head_table_result(FILE* out_file)
{
	int count_exec;
	
	fprintf(out_file, "|Reordering  | Threads |");
	for (count_exec = 1; count_exec <= TEST_EXEC_TIMES; ++count_exec)
		fprintf(out_file, " Exec %d  |", count_exec);
	fprintf(out_file, "\n");
	
	for (count_exec = 1; count_exec <= TEST_EXEC_TIMES + 2; ++count_exec)
		fprintf(out_file, "---------|");
	fprintf(out_file, "\n");
}




void run_all_tests_unordered_rcm(FILE* out_file, int* nthreads, int num_nthreads, char* matrix, int root)
{
	int count_percents, count_nthreads, count_exec;
	test_def tresult;
	
	int num_bfs_percents = 2;
	float bfs_chunk_percent[] = { .5, .8 };
	
	for (count_percents = 0; count_percents < num_bfs_percents; ++count_percents)
	{
		fprintf(out_file, "Dinamyc Chunk: %f\n\n", bfs_chunk_percent[count_percents]);
		printf("\nDinamyc Chunk: %f\n", bfs_chunk_percent[count_percents]);
		
		print_head_table_result(out_file);
		
		for (count_nthreads = 0; count_nthreads < num_nthreads; ++count_nthreads)
		{
			fprintf(out_file, "|Unordered      |   %d    |", nthreads[count_nthreads]);
			
			for (count_exec = 1; count_exec <= TEST_EXEC_TIMES; ++count_exec)
			{
				tresult = test_unordered_rcm(matrix, 
					nthreads[count_nthreads], bfs_chunk_percent[count_percents], root);
				fprintf(out_file, "%.6f |", tresult.time);
			}
			fprintf(out_file, "\n\n");
			fflush(out_file);
		}
	}
	
}


void run_all_tests_serial_rcm(FILE* out_file, char* matrix, int root)
{
	int count_exec;
	test_def tresult;
	
	fprintf(out_file, "|Serial      |    1    |");
			
	for (count_exec = 1; count_exec <= TEST_EXEC_TIMES; ++count_exec)
	{
		tresult = test_serial_rcm(matrix, root);
		fprintf(out_file, "%.6f |", tresult.time);
	}
	
	fprintf(out_file, "\n\n");
	fflush(out_file);
	
}



void run_all_tests_leveled_rcm(FILE* out_file, int* nthreads, int num_nthreads, char* matrix, int root, reorder_algorithm alg)
{
	int count_nthreads, count_exec;
	test_def tresult;
	char* alg_name;
	
	print_head_table_result(out_file);
	
	switch (alg)
	{
		case (leveled_rcm) :
			alg_name = "Leveled";
			break;
			
		case (leveled_rcm_v1) :
			alg_name = "Leveled v1";
			break;
			
		case (leveled_rcm_v2) :
			alg_name = "Leveled v2";
			break;
			
		default :
			break;
		
	}
	
	for (count_nthreads = 0; count_nthreads < num_nthreads; ++count_nthreads)
	{
		fprintf(out_file, "|%s      |   %d    |", alg_name, nthreads[count_nthreads]);
		
		long int vec_band[TEST_EXEC_TIMES];
		
		for (count_exec = 1; count_exec <= TEST_EXEC_TIMES; ++count_exec)
		{
			switch (alg)
			{
				case (leveled_rcm) :
					tresult = test_leveled_rcm(matrix, nthreads[count_nthreads], root);
					
					break;
					
				case (leveled_rcm_v1) :
					tresult = test_leveled_rcm_v1(matrix, nthreads[count_nthreads], root);
					break;
					
				case (leveled_rcm_v2) :
					tresult = test_leveled_rcm_v2(matrix, nthreads[count_nthreads], root);
					break;
					
				default :
					break;
				
			}
			
			vec_band[count_exec-1] = tresult.bandwidth;
			fprintf(out_file, "%.6f |", tresult.time);
			
			
		}
		
		fprintf(out_file, "\n");
		fprintf(out_file, "|%s      |   %d    |", alg_name, nthreads[count_nthreads]);
		
		for (count_exec = 0; count_exec < TEST_EXEC_TIMES; ++count_exec)
			fprintf(out_file, "%ld    |", vec_band[count_exec]);
		
		fprintf(out_file, "\n");
		fflush(out_file);
	}
	
}



void run_all_tests()
{
	int root, count_matrix, count_alg;
	FILE* out_file;
	
	int num_matrices = 6;
	char* matrices[] = {
		"../Big-Matrices/atmosmodj.mtx",
		"../Big-Matrices/Dubcova3.mtx",
		"../Big-Matrices/dw8192.mtx",
		"../Big-Matrices/inline_1.mtx"
		"../Big-Matrices/nlpkkt120.mtx",
		"../Big-Matrices/nlpkkt240.mtx"
// 		"../Matrices/aft01.mtx",
// 		"../Matrices/bcspwr01.mtx",
// 		"../Matrices/bcspwr02.mtx",
// 		"../Matrices/rail_5177.mtx"
// 		"../Matrices/FEM_3D_thermal1.mtx",
// 		"../Matrices/Dubcova2.mtx"
	};
	
	int num_nthreads = 6;
	int nthreads[] = { 4, 8, 16, 32, 64, 128 };
	
	int num_algorithms = 1;
	reorder_algorithm algorithm[] = { leveled_rcm, leveled_rcm_v1, leveled_rcm_v2 };
	
	if ((out_file = fopen("run_all_tests_output.txt", "w")) == NULL) 
		exit(1);

	
	for (count_matrix = 0; count_matrix < num_matrices; ++count_matrix)
	{
		fprintf(out_file, "-----------------------------------------------------------------------\n");
		fprintf(out_file, "Tests Execution - Matrix: %s\n", matrices[count_matrix]);
		fprintf(out_file, "-----------------------------------------------------------------------\n");
		fflush(out_file);
		printf("-------------- Matrix: %s\n", matrices[count_matrix]);
		
		root = get_node_peripheral(matrices[count_matrix]);
		
		print_head_table_result(out_file);
		
		for (count_alg = 0; count_alg < num_algorithms; ++count_alg)
		{
			switch (algorithm[count_alg])
			{
				case unordered_rcm :
					run_all_tests_unordered_rcm(out_file, nthreads, num_nthreads, matrices[count_matrix], root);
					break;
					
				case rcm :
					run_all_tests_serial_rcm(out_file, matrices[count_matrix], root);
					break;
					
				case leveled_rcm :
				case leveled_rcm_v1 :
				case leveled_rcm_v2 :
					run_all_tests_leveled_rcm(out_file, nthreads, num_nthreads, matrices[count_matrix], root, algorithm[count_alg]);
					break;
					
				default :
					break;
				
			}
			
		}
		
	}
	
	fclose(out_file);
}
