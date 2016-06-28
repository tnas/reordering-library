#include "test_suite_reordering.h"

void mc60hd_(int* n, int* nsup, int* lirn, int* irn, int* icptr, int* vars, int* mask, int* ls, int* xls, int* list, int* info);


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


int* get_node_peripheral_hsl(const char* path_matrix_file) {
	
	MAT* matrix;
	FILE* matrix_file;
	int* peripheral_nodes;
	
	peripheral_nodes = calloc(2, sizeof(int));
	
	if ((matrix_file = fopen(path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);

	int i;
	int n        = matrix->n; 
	int nsup     = n;
	int *irn     = matrix->JA;
	int *icptr   = matrix->IA;
	int lirn     = matrix->nz;
	int *vars;
	int *mask;
	int *ls;
	int *xls;
	int *list;
	int info[6];
	
	vars = calloc(n, sizeof(int));
	mask = calloc(n, sizeof(int));
	ls   = calloc(n, sizeof(int));
	xls  = calloc(n, sizeof(int));
	list = calloc(n, sizeof(int));
	
	for (i = 0; i < n; i++) 
	{
		++irn[i];
		++icptr[i];
		vars[i] = mask[i] = 1;
	}
	
	for (i = n; i < lirn; i++) ++irn[i];
	
	mc60hd_(&n, &nsup, &lirn, irn, icptr, vars, mask, ls, xls, list, info);
	
	for (i = 0; i < n; i++) 
	{
		++irn[i];
		++icptr[i];
	}
	
	for (i = n; i < lirn; i++) ++irn[i];
	
	peripheral_nodes[0] = info[0] - 1;
	peripheral_nodes[1] = info[1] - 1;
	
	MATRIX_clean(matrix);
	free(vars);
	free(mask);
	free(ls);
	free(xls);
	free(list);
	
	return peripheral_nodes;
}

int inline is_hsl_algorithm(reorder_algorithm algorithm)
{
	if (algorithm == hsl_rcm || 
	    algorithm == hsl_spectral || 
	    algorithm == hsl_sloan)
		return 1;
	
	return 0;
}

int inline is_serial_algorithm(reorder_algorithm algorithm)
{
	if (algorithm == serial_rcm || 
	    algorithm == serial_sloan ||
	    is_hsl_algorithm(algorithm))
		return 1;
	
	return 0;
}

int inline is_parallel_algorithm(reorder_algorithm algorithm)
{
	if (algorithm == unordered_rcm || algorithm == leveled_rcm ||
	    algorithm == bucket_rcm    || algorithm == parallel_sloan)
		return 1;
	
	return 0;
}


int inline is_sloan_algorithm(reorder_algorithm algorithm)
{
	if (algorithm == hsl_sloan || algorithm == parallel_sloan)
		return 1;
	
	return 0;
}

double inline get_total_time(test defs)
{
	return defs.time_peripheral + defs.time_reordering + defs.time_permutation;
}


test test_reorder_algorithm(test defs)
{
	int* permutation;
	double time;
	MAT* matrix;
	FILE* matrix_file;
	int* g;
	int* peripheral_nodes;
	
	if ((matrix_file = fopen(defs.path_matrix_file, "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, matrix_file);
	fclose(matrix_file);
	write_output_before(matrix);
	defs.original_band = MATRIX_bandwidth(matrix);
	
	// Setting threads for parallel algorithms
	if (is_parallel_algorithm(defs.algorithm)) 
	{
		omp_set_num_threads(defs.num_threads);
		
		// Getting pseudo peripheral nodes
		time = omp_get_wtime();
		peripheral_nodes = get_node_peripheral_hsl(defs.path_matrix_file);
		defs.time_peripheral = (omp_get_wtime() - time)/100.0;
	}
	else
	{
		// Getting pseudo peripheral nodes
		time = get_time();
		peripheral_nodes = get_node_peripheral_hsl(defs.path_matrix_file);
		defs.time_peripheral = (get_time() - time)/100.0;
	}
	
	defs.root       = peripheral_nodes[START];
	defs.start_node = peripheral_nodes[START];
	defs.end_node   = peripheral_nodes[END];
	
	switch (defs.algorithm)
	{
		case serial_rcm : // t = 0
			defs.algorithm = serial_rcm;
			defs.algorithm_name = "Serial RCM";
			time = get_time(); 
			REORDERING_RCM_opt(matrix, &permutation, defs.root);
			defs.time_reordering = (get_time() - time)/100.0;
			break;
			
		case serial_sloan : // t = 1
			defs.algorithm = serial_sloan;
			defs.algorithm_name = "Serial Sloan";
			g = GRAPH_LS_peripheral (matrix, &(defs.start_node), &(defs.end_node));
			free(g);
			
			time = get_time();
			REORDERING_SLOAN(matrix, &permutation, defs.start_node, defs.end_node);
			defs.time_reordering = (get_time() - time)/100.0;
			break;
			
		case hsl_rcm : // t = 2
			defs.algorithm = hsl_rcm;
			defs.algorithm_name = "HSL RCM";
			time = get_time(); 
			defs.reorder_band = REORDERING_HSL_RCM(matrix);
			defs.time_reordering = (get_time() - time)/100.0;
			break;
			
		case hsl_spectral : // t = 3
			defs.algorithm      = hsl_spectral;
			defs.algorithm_name = "HSL Spectral";
			time = get_time(); 
			REORDERING_HSL_SPECTRAL(matrix, &permutation);
			defs.time_reordering = (get_time() - time)/100.0;
			break;
			
		case hsl_sloan : // t = 4
			defs.algorithm_name = "HSL Sloan";
			defs.algorithm      = hsl_sloan;
			time = get_time();
			defs.wavefront = REORDERING_SLOAN_HSL(matrix);
			defs.time_reordering = (get_time() - time)/100.0;
			break;
			
		case unordered_rcm : // t = 5
			defs.algorithm_name = "Unordered RCM";
			defs.algorithm      = unordered_rcm;
			time = omp_get_wtime();
			Unordered_RCM(matrix, &permutation, defs.root, defs.percent_chunk);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case leveled_rcm : // t = 6
			defs.algorithm_name = "Leveled RCM";
			defs.algorithm      = leveled_rcm;
			time = omp_get_wtime();
			Leveled_RCM(matrix, &permutation, defs.root);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case bucket_rcm : // t = 7
			defs.algorithm_name = "Bucket RCM";
			defs.algorithm      = bucket_rcm;
			time = omp_get_wtime(); 
			Bucket_RCM(matrix, &permutation, defs.root);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case parallel_sloan : // t = 8
			defs.algorithm_name = "Parallel Sloan";
			defs.algorithm      = parallel_sloan;
			time = omp_get_wtime();
			Parallel_Sloan(matrix, &permutation, defs.start_node, defs.end_node);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		default :
			printf("*** [Error] Algorithm must be between 0 and 8 ***\n");
			exit(1);
	}
	
	// Setting threads for parallel algorithms
	if (is_parallel_algorithm(defs.algorithm)) omp_set_num_threads(defs.num_threads);
	
	if (is_sloan_algorithm(defs.algorithm))
	{
		if (is_hsl_algorithm(defs.algorithm))
		{
			printf("%s: Wavefront [ %ld ] Time [ %.6f ]\n", defs.algorithm_name, defs.wavefront, defs.time_reordering); 
			fflush(stdout);
		}
		else
		{
			if (is_parallel_algorithm(defs.algorithm))
			{
				time = omp_get_wtime();
				MATRIX_PARALLEL_permutation(matrix, permutation);
				defs.wavefront = MATRIX_PARALLEL_wavefront(matrix);
				defs.time_permutation = (omp_get_wtime() - time)/100.0;
			}
			else
			{
				time = get_time();
				MATRIX_PARALLEL_permutation(matrix, permutation);
				defs.wavefront = MATRIX_PARALLEL_wavefront(matrix);
				defs.time_permutation = (get_time() - time)/100.0;
			}
			
			printf("%s: Wavefront [ %ld ] => Time (Periph/Reorder/Permut/Total) [ %.6f || %.6f || %.6f || %.6f ]\n",
				defs.algorithm_name, defs.wavefront, 
				defs.time_peripheral, defs.time_reordering, defs.time_permutation, get_total_time(defs)); 
			fflush(stdout);
		}
	}
	else
	{
		if (is_hsl_algorithm(defs.algorithm))
		{
			printf("%s: (Before/After) [ %ld/%ld ] Time [ %.6f ]\n", defs.algorithm_name, 
			       defs.original_band, defs.reorder_band, defs.time_reordering); 
			fflush(stdout);
		}
		else 
		{
			if (is_parallel_algorithm(defs.algorithm))
			{
				time = omp_get_wtime();
				MATRIX_PARALLEL_permutation(matrix, permutation);
				defs.reorder_band = MATRIX_PARALLEL_wavefront(matrix);
				defs.time_permutation = (omp_get_wtime() - time)/100.0;
			}
			else
			{
				time = get_time();
				MATRIX_PARALLEL_permutation(matrix, permutation);
				defs.reorder_band = MATRIX_PARALLEL_wavefront(matrix);
				defs.time_permutation = (get_time() - time)/100.0;
			}
		
			printf("%s: Bandwidth (Before/After) [ %ld/%ld ] => Time (Periph/Reorder/Permut/Total) [ %.6f || %.6f || %.6f || %.6f ]\n",
				defs.algorithm_name, defs.original_band, defs.reorder_band, 
				defs.time_peripheral, defs.time_reordering, defs.time_permutation, get_total_time(defs)); 
			fflush(stdout);
			
			free(permutation);
		}
		
	}
	
	write_output_after(matrix);
	MATRIX_clean(matrix);
	free(peripheral_nodes);
	
	return defs;
}


long int normalize_int_results(const long int* results)
{
	int exec_min_val, exec_max_val, exec;
	long int sum_val, max_val, min_val;
	
	max_val = sum_val = 0;
	min_val = INT_MAX;
	
	if (TEST_EXEC_TIMES == 1) return results[0];
	
	// Finding out max/min value
	for (exec = 0; exec < TEST_EXEC_TIMES; ++exec)
	{
		if (results[exec] > max_val)
		{
			exec_max_val = exec;
			max_val = results[exec];
		}
		else if (results[exec] < min_val)
		{
			exec_min_val = exec;
			min_val = results[exec];
		}
	}
	
	// Discarding max and min results value
	for (exec = 0; exec < TEST_EXEC_TIMES; ++exec)
	{
		if (exec != exec_min_val && exec != exec_max_val)
			sum_val += results[exec];
	}
	
	return sum_val / (TEST_EXEC_TIMES - 2);
}


double normalize_double_results(const double* results)
{
	int exec_min_val, exec_max_val, exec;
	double sum_val, max_val, min_val;
	
	max_val = sum_val = 0;
	min_val = INT_MAX;
	
	if (TEST_EXEC_TIMES == 1) return results[0];
	
	// Finding out max/min value
	for (exec = 0; exec < TEST_EXEC_TIMES; ++exec)
	{
		if (results[exec] > max_val)
		{
			exec_max_val = exec;
			max_val = results[exec];
		}
		else if (results[exec] < min_val)
		{
			exec_min_val = exec;
			min_val = results[exec];
		}
	}
	
	// Discarding max and min results value
	for (exec = 0; exec < TEST_EXEC_TIMES; ++exec)
	{
		if (exec != exec_min_val && exec != exec_max_val)
			sum_val += results[exec];
	}
	
	return sum_val / (TEST_EXEC_TIMES - 2);
}


void run_all_tests()
{
	int count_matrix, count_alg, exec, num_threads,
	    count_nthreads, num_matrices, size_set_nthreads, num_algorithms;
	FILE* out_file;
	test result;
	double time_reorderings[TEST_EXEC_TIMES];
	double time_peripherals[TEST_EXEC_TIMES];
	double time_permutations[TEST_EXEC_TIMES];
	long int bandwidths[TEST_EXEC_TIMES];
	long int wavefronts[TEST_EXEC_TIMES];
	
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
// 		"../Matrices/FEM_3D_thermal1.mtx",
// 		"../Matrices/Dubcova2.mtx"
	};
	
	int nthreads[] = { 4, 8, 16, 32, 64, 128 };
	
	reorder_algorithm algorithm[] = { bucket_rcm };
	
	/* *****************
	 * Tests execution
	 * *****************
	 */
	
	if ((out_file = fopen("run_all_tests_normalized_output.txt", "w")) == NULL) 
		exit(1);
	
	num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	num_algorithms    = sizeof(algorithm)/sizeof(algorithm[0]);

	for (count_matrix = 0; count_matrix < num_matrices; ++count_matrix)
	{
		fprintf(out_file, "-----------------------------------------------------------------------\n");
		fprintf(out_file, "Tests Execution - Matrix: %s\n", matrices[count_matrix]); fflush(out_file);
		
		for (count_alg = 0; count_alg < num_algorithms; ++count_alg)
		{
			fprintf(out_file, "-----------------------------------------------------------------------\n");
			fprintf(out_file, "Algorithm: %d\n", algorithm[count_alg]);
			fprintf(out_file, "-----------------------------------------------------------------------\n");
			fflush(out_file);
			
			num_threads = is_serial_algorithm(algorithm[count_alg]) ? 1 :
				size_set_nthreads;
			
			for (count_nthreads = 0; count_nthreads < num_threads; ++count_nthreads)
			{
				result.algorithm        = algorithm[count_alg];
				result.path_matrix_file = matrices[count_matrix];
				result.num_threads      = nthreads[count_nthreads];
				result.percent_chunk    = bfs_chunk_percent;
					
				for (exec = 0; exec < TEST_EXEC_TIMES; ++exec)
				{
					result = test_reorder_algorithm(result);
					time_reorderings[exec]  = result.time_reordering;
					time_peripherals[exec]  = result.time_peripheral;
					time_permutations[exec] = result.time_permutation;
					bandwidths[exec]        = result.reorder_band;
					wavefronts[exec]        = result.wavefront;
				}
				
				result.time_peripheral  = normalize_double_results(time_peripherals); 
				result.time_reordering  = normalize_double_results(time_reorderings);
				result.time_permutation = normalize_double_results(time_permutations);
				result.reorder_band     = normalize_int_results(bandwidths);
				result.wavefront        = normalize_int_results(wavefronts);
				
				if (is_sloan_algorithm(result.algorithm))
				{
					if (is_hsl_algorithm(result.algorithm))
					{
						fprintf(out_file, "%s: Wavefront [ %ld ] Time [ %.6f ]\n", 
						       result.algorithm_name, result.wavefront, result.time_reordering); 
						fflush(out_file);
					}
					else
					{
						fprintf(out_file, "[%s] Threads: %d -- Wavefront [ %ld ] -- Time (Periph/Reorder/Permut/Total) [ %.6f || %.6f || %.6f || %.6f ]\n",
							result.algorithm_name, is_serial_algorithm(algorithm[count_alg]) ? 1 : nthreads[count_nthreads],
							result.wavefront, result.time_peripheral, result.time_reordering, result.time_permutation, get_total_time(result)); 
						fflush(out_file);
					}
				}
				else
				{
					if (is_hsl_algorithm(result.algorithm))
					{
						fprintf(out_file, "%s: (Before/After) [ %ld/%ld ] Time [ %.6f ]\n", result.algorithm_name, 
						result.original_band, result.reorder_band, result.time_reordering); 
						fflush(out_file);
					}
					else 
					{
						fprintf(out_file, "[%s] Threads: %d -- Bandwidth (Before/After) [ %ld/%ld ] -- Time (Periph/Reorder/Permut/Total) [ %.6f || %.6f || %.6f || %.6f ]\n",
							result.algorithm_name, is_serial_algorithm(algorithm[count_alg]) ? 1 : nthreads[count_nthreads],
							result.original_band, result.reorder_band, result.time_peripheral, 
							result.time_reordering, result.time_permutation, get_total_time(result)); 
						fflush(out_file);
					}
					
				}
			}
		}
	}
	
	fclose(out_file);
}
