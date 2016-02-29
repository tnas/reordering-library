#include "rcm_parallel_test.h"

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


test_result run_test_serial_rcm(const char* path_matrix_file, int root)
{
	long int bandwidth, envelope, bandwidth_after, envelope_after;
	int* permutation;
	double time;
	MAT* matrix;
	FILE* matrix_file;
	test_result result;
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	
	write_output_before(matrix);
	
	bandwidth_after = MATRIX_bandwidth(matrix);
	envelope_after  = MATRIX_envelope(matrix);
	
	time = get_time(); 
	REORDERING_RCM_opt(matrix, &permutation, root);
	time = (get_time() - time)/100.0;
	result.time = time;
	
	MATRIX_permutation(matrix, permutation);
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	result.bandwidth = bandwidth;
	
	write_output_after(matrix);
	
	free(permutation);
	MATRIX_clean(matrix);
	
	printf("Serial RCM: Band/Env [ %ld / %ld => %ld / %ld ] Time [ %.6f ]\n",
		bandwidth_after, envelope_after, bandwidth, envelope, time); fflush(stdout);
	
	return result;
}


test_result run_test_leveled_rcm(const char* path_matrix_file, int root)
{
	long int bandwidth, envelope, bandwidth_after, envelope_after;
	int* permutation;
	double time;
	MAT* matrix;
	FILE* matrix_file;
	test_result result;
	
	permutation = NULL;
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	
	write_output_before(matrix);
	
	bandwidth_after = MATRIX_bandwidth(matrix);
	envelope_after  = MATRIX_envelope(matrix);
	
	time = get_time(); 
	Leveled_RCM(matrix, permutation, root);
	time = (get_time() - time)/100.0;
	result.time = time;
	
	MATRIX_permutation(matrix, permutation);
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	result.bandwidth = bandwidth;
	
	write_output_after(matrix);
	
	free(permutation);
	MATRIX_clean(matrix);
	
	printf("Leveled RCM: Band/Env [ %ld / %ld => %ld / %ld ] Time [ %.6f ]\n",
		bandwidth_after, envelope_after, bandwidth, envelope, time); fflush(stdout);
	
	return result;
}





test_result run_test_unordered_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_size, int root)
{
	long int bandwidth, envelope;
	int* permutation;
	double time;
	MAT* matrix;
	FILE* matrix_file;
	test_result result;
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	
	omp_set_num_threads(num_threads);
	
	time = get_time(); 
	Unordered_RCM(matrix, &permutation, root, bfs_chunk_size);
	time = (get_time() - time)/100.0;
	result.time = time;
	
	MATRIX_permutation(matrix, permutation); 
	
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	result.bandwidth = bandwidth;
	
	free(permutation);
	MATRIX_clean(matrix);
	
	printf("Unordered RCM (Threads %d, Chunk %f): Band/Env [ %ld / %ld ] Time [ %.6f ]\n",
		num_threads, bfs_chunk_size, bandwidth, envelope, time); fflush(stdout);
	
	return result;
}


int get_node_peripheral(const char* path_matrix_file) {
	
	int root, e;
	MAT* matrix;
	FILE* matrix_file;
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	
	GRAPH_LS_peripheral(matrix, &root, &e);
	
	return root;
}


void run_test_serial_parallel_rcm(const char* path_matrix_file, const int num_threads, const float bfs_chunk_size)
{
	
	int root;
	
	root = get_node_peripheral(path_matrix_file);
	
	run_test_serial_rcm(path_matrix_file, root);
	run_test_unordered_rcm(path_matrix_file, num_threads, bfs_chunk_size, root);
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


void run_all_tests()
{
	int root, count_matrix, count_percents, count_exec, count_nthreads;
	FILE* out_file;
	test_result tresult;
	
	int num_matrices = 1;
	char* matrices[] = {
// 		"../Big-Matrices/atmosmodj.mtx",
// 		"../Big-Matrices/Dubcova3.mtx",
// 		"../Big-Matrices/dw8192.mtx",
		"../Big-Matrices/inline_1.mtx"
// 		"../Big-Matrices/nlpkkt120.mtx",
// 		"../Big-Matrices/nlpkkt240.mtx"
	};	
																																
	int num_bfs_percents = 2;
	float bfs_chunk_percent[] = { .5, .8 };
	
	int num_nthreads = 3;
	int nthreads[] = { 32, 64, 128 };
	
	
	if ((out_file = fopen("run_tests_output.txt", "w")) == NULL) 
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
		
		fprintf(out_file, "|Serial      |    1    |");
		for (count_exec = 1; count_exec <= TEST_EXEC_TIMES; ++count_exec)
		{
			tresult = run_test_serial_rcm(matrices[count_matrix], root);
			fprintf(out_file, "%.6f |", tresult.time);
		}
		fprintf(out_file, "\n\n");
		fflush(out_file);
		
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
					tresult = run_test_unordered_rcm(matrices[count_matrix], 
						nthreads[count_nthreads], bfs_chunk_percent[count_percents], root);
					fprintf(out_file, "%.6f |", tresult.time);
				}
				fprintf(out_file, "\n\n");
				fflush(out_file);
			}
		}
	}
	
	fclose(out_file);
}
