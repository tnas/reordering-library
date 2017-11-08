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

#include "test_suite_reordering.h"


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
	if (algorithm == serial_rcm || algorithm == serial_sloan || 
	    algorithm == boost_rcm ||  algorithm == boost_sloan ||
	    is_hsl_algorithm(algorithm))
		return 1;
	
	return 0;
}

int inline is_parallel_algorithm(reorder_algorithm algorithm)
{
	if (!is_serial_algorithm(algorithm)) return 1;
	return 0;
}


int inline is_sloan_algorithm(reorder_algorithm algorithm)
{
	if (algorithm == hsl_sloan || algorithm == relaxed_order_sloan ||
	    algorithm == serial_sloan || algorithm == boost_sloan ||
		algorithm == logbag_sloan)
		return 1;
	
	return 0;
}

double inline get_total_time(test defs)
{
	return defs.time_peripheral + defs.time_reordering;
}


test test_reorder_algorithm(test defs)
{
	int* permutation;
	double time;
	MAT* matrix;
	int* peripheral_nodes;
	graph_diameter* diameter;
	METAGRAPH* mgraph;

	diameter         = NULL;
	peripheral_nodes = NULL;
	
	// Setting threads for parallel algorithms
	if (is_parallel_algorithm(defs.algorithm)) omp_set_num_threads(defs.num_threads);
	
	MATRIX_read_from_path(defs.path_matrix_file, &matrix);
	mgraph = GRAPH_parallel_build_METAGRAPH(matrix);
	
	MATRIX_write_gnuplot(matrix, "original_matrix");
	
	if (is_sloan_algorithm(defs.algorithm))
	{
		defs.original_wavefront = MATRIX_envelope(matrix);
// 		defs.original_wavefront = MATRIX_PARALLEL_max_wavefront(matrix);
	}
	else
	{
		defs.original_band = MATRIX_bandwidth(matrix);
	}
		
	// Getting pseudo peripheral nodes
	if (is_hsl_algorithm(defs.algorithm))
	{
		time = omp_get_wtime();
		peripheral_nodes = get_pseudo_diameter_hsl(matrix);
		defs.time_peripheral = (omp_get_wtime() - time)/100.0;
		
		defs.root       = peripheral_nodes[START];
		defs.start_node = peripheral_nodes[START];
		defs.end_node   = peripheral_nodes[END];
	}
	else
	{
		time = omp_get_wtime();
		diameter = GRAPH_parallel_pseudodiameter(mgraph, VERTEX_BY_DEGREE);
// 		diameter = GRAPH_parallel_pseudodiameter(mgraph, HALF_SORTED);
// 		diameter = GRAPH_parallel_pseudodiameter(mgraph, FIVE_NON_ADJACENT);
		
		defs.time_peripheral = (omp_get_wtime() - time)/100.0;
		
		defs.root       = diameter->start;
		defs.start_node = diameter->start;
		defs.end_node   = diameter->end;
	}
	
	switch (defs.algorithm)
	{
		case serial_rcm : // t = 0
			defs.algorithm = serial_rcm;
			defs.algorithm_name = "Serial RCM";
			time = omp_get_wtime(); 
			REORDERING_RCM_opt(matrix, &permutation, defs.root);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case serial_sloan : // t = 1
			defs.algorithm = serial_sloan;
			defs.algorithm_name = "Serial Sloan";
			time = omp_get_wtime();
			REORDERING_SLOAN(matrix, &permutation, defs.start_node, defs.end_node);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case hsl_rcm : // t = 2
			defs.algorithm = hsl_rcm;
			defs.algorithm_name = "HSL RCM";
			time = omp_get_wtime();
			defs.reorder_band = Reordering_RCM_pseudodiameter_HSL(matrix, defs.start_node, defs.end_node);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			defs.time_permutation = 0;
			break;
			
		case hsl_spectral : // t = 3
			defs.algorithm      = hsl_spectral;
			defs.algorithm_name = "HSL Spectral";
			time = omp_get_wtime(); 
			REORDERING_HSL_SPECTRAL(matrix, &permutation);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			defs.time_permutation = 0;
			break;
			
		case hsl_sloan : // t = 4
			defs.algorithm_name = "HSL Sloan";
			defs.algorithm      = hsl_sloan;
			time = omp_get_wtime();
			defs.reorder_wavefront = REORDERING_SLOAN_pseudodiameter_HSL(matrix, defs.start_node, defs.end_node);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			defs.time_permutation = 0;
			break;
			
		case unordered_rcm : // t = 5
			defs.algorithm_name = "Unordered RCM";
			defs.algorithm      = unordered_rcm;
			time = omp_get_wtime();
			Unordered_RCM(mgraph, &permutation, defs.root, defs.percent_chunk);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case leveled_rcm : // t = 6
			defs.algorithm_name = "Leveled RCM";
			defs.algorithm      = leveled_rcm;
			time = omp_get_wtime();
			Leveled_RCM(mgraph, &permutation, defs.root);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case bucket_rcm : // t = 7
			defs.algorithm_name = "Bucket RCM";
			defs.algorithm      = bucket_rcm;
			time = omp_get_wtime(); 
			Bucket_RCM(mgraph, &permutation, defs.root);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case relaxed_order_sloan : // t = 8
			defs.algorithm_name = "Relaxed Order Sloan";
			defs.algorithm      = relaxed_order_sloan;
			time = omp_get_wtime();
			Parallel_Relaxed_Order_Sloan(mgraph, &permutation, defs.start_node, defs.end_node);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		case boost_rcm : // t = 9
			defs.algorithm_name  = "Boost RCM";
			defs.algorithm       = boost_rcm;
			defs.time_reordering = Boost_RCM(mgraph, &permutation, defs.root);
			break;
			
		case boost_sloan : // t = 10
			defs.algorithm_name  = "Boost Sloan";
			defs.algorithm       = boost_sloan;
			defs.time_reordering = Boost_Sloan(mgraph, &permutation, defs.start_node, defs.end_node);
			break;
		
		case logbag_sloan : // t = 11
			defs.algorithm_name = "Parallel Logical Bag Sloan";
			defs.algorithm      = logbag_sloan;
			time = omp_get_wtime();
			Parallel_Logical_Bag_Sloan(mgraph, &permutation, defs.start_node, defs.end_node);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
		
		case shrinked_rcm : // t = 12
			defs.algorithm_name = "Shrinked RCM";
			defs.algorithm      = shrinked_rcm;
			time = omp_get_wtime(); 
			Shrinked_RCM(mgraph, &permutation, defs.root);
			defs.time_reordering = (omp_get_wtime() - time)/100.0;
			break;
			
		default :
			printf("*** [Error] Algorithm must be between 0 and 13 ***\n");
			exit(1);
	}
	
	// Setting threads for parallel algorithms
	if (is_parallel_algorithm(defs.algorithm)) omp_set_num_threads(defs.num_threads);
	
	if (is_sloan_algorithm(defs.algorithm))
	{
		if (is_hsl_algorithm(defs.algorithm))
		{
			printf("%s: (Before/After) [ %ld/%ld ] Time => (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n", defs.algorithm_name, 
			       defs.original_wavefront, defs.reorder_wavefront, defs.time_reordering, defs.time_reordering, get_total_time(defs)); 
			fflush(stdout);
		}
		else
		{
			time = omp_get_wtime();
			MATRIX_permutation(matrix, permutation);
			defs.time_permutation = (omp_get_wtime() - time)/100.0;
			defs.reorder_wavefront = MATRIX_envelope(matrix);
// 			defs.reorder_wavefront = MATRIX_PARALLEL_max_wavefront(matrix);
			
			printf("%s: (Before/After) [ %ld/%ld ] => Time (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n",
				defs.algorithm_name, defs.original_wavefront, defs.reorder_wavefront, 
				defs.time_peripheral, defs.time_reordering, get_total_time(defs)); 
			fflush(stdout);
			
			free(permutation);
		}
	}
	else // RCM Algorithms
	{ 
		if (is_hsl_algorithm(defs.algorithm))
		{
			printf("%s: (Before/After) [ %ld/%ld ] => Time (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n", defs.algorithm_name, 
			       defs.original_band, defs.reorder_band, defs.time_peripheral, defs.time_reordering, get_total_time(defs)); 
			fflush(stdout);
		}
		else 
		{
			if (is_parallel_algorithm(defs.algorithm))
			{
				time = omp_get_wtime();
				MATRIX_PARALLEL_permutation(matrix, permutation);
				defs.time_permutation = (omp_get_wtime() - time)/100.0;
				defs.reorder_band = MATRIX_PARALLEL_bandwidth(matrix);
			}
			else
			{
				time = omp_get_wtime();
				MATRIX_permutation(matrix, permutation);
				defs.time_permutation = (omp_get_wtime() - time)/100.0;
				defs.reorder_band = MATRIX_bandwidth(matrix);
			}
		
			printf("%s: Bandwidth (Before/After) [ %ld/%ld ] => Time (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n",
				defs.algorithm_name, defs.original_band, defs.reorder_band, 
				defs.time_peripheral, defs.time_reordering, get_total_time(defs)); 
			fflush(stdout);
			
			free(permutation);
		}
		
	}
	
	MATRIX_write_gnuplot(matrix, "reordered_matrix");
	GRAPH_parallel_destroy_METAGRAPH(mgraph);
	if (peripheral_nodes != NULL) free(peripheral_nodes);
	if (diameter != NULL) free(diameter);
	
	return defs;
}



void execute_test_suite(const char** matrices, const int* nthreads, const reorder_algorithm* algorithm, const int exec_times,
	const int num_matrices, const int size_set_nthreads, const int num_algorithms)
{
	int count_matrix, count_alg, exec, num_threads, count_nthreads;
	FILE* out_file;
	test result;
	statistic norm_values;
	double time_reorderings[exec_times];
	double time_peripherals[exec_times];
	double time_permutations[exec_times];
	long int bandwidths[exec_times];
	long int wavefronts[exec_times];
	
	float bfs_chunk_percent = .5;
	
	if ((out_file = fopen("reordering_library_output.txt", "w")) == NULL) 
		exit(1);
	
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
					
				for (exec = 0; exec < exec_times; ++exec)
				{
					result = test_reorder_algorithm(result);
					time_reorderings[exec]  = result.time_reordering;
					time_peripherals[exec]  = result.time_peripheral;
					time_permutations[exec] = result.time_permutation;
					bandwidths[exec]        = result.reorder_band;
					wavefronts[exec]        = result.reorder_wavefront;
				}
				
				normalize_results(time_peripherals, exec_times, &norm_values);
				result.time_peripheral = norm_values.average_value; 
				
				normalize_results(time_reorderings, exec_times, &norm_values);
				result.time_reordering = norm_values.average_value; 
				
				normalize_results(time_permutations, exec_times, &norm_values);
				result.time_permutation = norm_values.average_value; 
				
				normalize_int_results(bandwidths, exec_times, &norm_values);
				result.reorder_band = norm_values.average_value; 
				
				normalize_int_results(wavefronts, exec_times, &norm_values);
				result.reorder_wavefront = norm_values.average_value; 
				
				if (is_sloan_algorithm(result.algorithm))
				{
					if (is_hsl_algorithm(result.algorithm))
					{
						fprintf(out_file, "%s: Wavefront [ %ld/%ld ] Time => (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n", result.algorithm_name, 
							result.original_wavefront, result.reorder_wavefront, result.time_peripheral, result.time_reordering, get_total_time(result)); 
						fflush(out_file);
					}
					else
					{
						fprintf(out_file, "[%s] Threads: %d -- Wavefront [ %ld/%ld ] -- Time (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n",
							result.algorithm_name, is_serial_algorithm(algorithm[count_alg]) ? 1 : nthreads[count_nthreads],
							result.original_wavefront, result.reorder_wavefront, result.time_peripheral, result.time_reordering,  
							get_total_time(result)); 
						fflush(out_file);
					}
				}
				else
				{
					if (is_hsl_algorithm(result.algorithm))
					{
						fprintf(out_file, "%s: (Before/After) [ %ld/%ld ] => Time (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n", result.algorithm_name, 
							result.original_band, result.reorder_band, result.time_peripheral, result.time_reordering, get_total_time(result)); 
						fflush(out_file);
					}
					else 
					{
						fprintf(out_file, "[%s] Threads: %d -- Bandwidth (Before/After) [ %ld/%ld ] -- Time (Periph/Reorder/Total) [ %.6f || %.6f || %.6f ]\n",
							result.algorithm_name, is_serial_algorithm(algorithm[count_alg]) ? 1 : nthreads[count_nthreads],
							result.original_band, result.reorder_band, result.time_peripheral, 
							result.time_reordering, get_total_time(result)); 
						fflush(out_file);
					}
				}
			}
		}
	}
	
	fclose(out_file);
}


// t = 1
void run_reordering_tests()
{
	int num_executions  = 1;
	
	const char* matrices[] = {
// 		"../Matrices/apothen.mtx",
// 		"../Matrices/fidap001.mtx",
// 		"../Matrices/fidapm08.mtx",
// 		"../Matrices/aft01.mtx",
// 		"../Matrices/hsl.mtx",
// 		"../Matrices/sample.mtx",
// 		"../Matrices/bcspwr01.mtx",
// 		"../Matrices/can24.mtx",
// 		"./Matrices/bcspwr02.mtx",
// 		"../Matrices/rail_5177.mtx",
// 		"./Matrices/FEM_3D_thermal1.mtx",
// 		"../Matrices/Dubcova2.mtx",
		"../Matrices/thermomech_TC.mtx",
	};
	
	int nthreads[] = { 2 };
	
	reorder_algorithm algorithms[] = {hsl_sloan, relaxed_order_sloan};
	
	int num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	int size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	int num_algorithms    = sizeof(algorithms)/sizeof(algorithms[0]);
	
	// Running tests
	execute_test_suite(matrices, nthreads, algorithms, num_executions, 
		      num_matrices, size_set_nthreads, num_algorithms);
}

// t = 7
void run_tema_journal_tests()
{
	int num_executions  = 5;
	
	const char* matrices[] = {
		"../Matrices/benzene.mtx",
		"../Matrices/FEM_3D_thermal1.mtx",
		"../Matrices/rail_79841.mtx",
		"../Matrices/thermomech_TC.mtx",
		"../Matrices/Dubcova3.mtx",
		"../Matrices/Ga41As41H72.mtx",
		"../Matrices/F1.mtx",
		"../Matrices/helm2d03.mtx",
		"../Matrices/msdoor.mtx",
		"../Matrices/inline_1.mtx",		
		"../Matrices/gsm_106857.mtx",
		"../Matrices/Fault_639.mtx",
		"../Matrices/tmt_sym.mtx",
		"../Matrices/tmt_unsym.mtx",
		"../Matrices/audikw_1.mtx",
		"../Matrices/nlpkkt80.mtx",
		"../Matrices/dielFilterV3real.mtx",
		"../Matrices/dielFilterV2real.mtx",
		"../Matrices/Serena.mtx",
		"../Matrices/G3_circuit.mtx",
	};
	
	int nthreads[] = { 1, 2, 4, 6, 8, 10, 12 };
	
	reorder_algorithm algorithms[] = { boost_rcm, unordered_rcm, bucket_rcm };
	
	int num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	int size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	int num_algorithms    = sizeof(algorithms)/sizeof(algorithms[0]);
	
	// Running tests
	execute_test_suite(matrices, nthreads, algorithms, num_executions, 
		      num_matrices, size_set_nthreads, num_algorithms);
}

// t = 8
void run_dissertation_largest_matrices()
{
	int num_executions  = 5;
	
	const char* matrices[] = {
		"../Big-Matrices/10-inline_1.mtx",
		"../Big-Matrices/09-gsm_106857.mtx",
		"../Big-Matrices/08-Fault_639.mtx",
		"../Big-Matrices/07-tmt_sym.mtx",
		"../Big-Matrices/06-boneS10.mtx",
		"../Big-Matrices/05-audikw_1.mtx",
		"../Big-Matrices/04-nlpkkt80.mtx",
		"../Big-Matrices/03-dielFilterV2real.mtx",
		"../Big-Matrices/02-Serena.mtx",
		"../Big-Matrices/01-G3_circuit.mtx",		
	};
	
	int nthreads[] = { 1, 2, 4, 6, 8, 10, 12 };
	
	reorder_algorithm algorithms[] = { hsl_rcm, boost_rcm, unordered_rcm, shrinked_rcm, bucket_rcm, hsl_sloan, boost_sloan, logbag_sloan, relaxed_order_sloan };
	
	int num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	int size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	int num_algorithms    = sizeof(algorithms)/sizeof(algorithms[0]);
	
	// Running tests
	execute_test_suite(matrices, nthreads, algorithms, num_executions, 
		      num_matrices, size_set_nthreads, num_algorithms);
}

// t = 9
void run_dissertation_smallest_matrices()
{
	int num_executions  = 5;
	
	const char* matrices[] = {
		"../Big-Matrices/20-m_t1.mtx",
		"../Big-Matrices/19-filter3D.mtx",
		"../Big-Matrices/18-SiO2.mtx",
		"../Big-Matrices/17-d_pretok.mtx",
		"../Big-Matrices/16-CO.mtx",
		"../Big-Matrices/15-offshore.mtx",
		"../Big-Matrices/14-Ga41As41H72.mtx",
		"../Big-Matrices/13-F1.mtx",
		"../Big-Matrices/12-mario002.mtx",
		"../Big-Matrices/11-msdoor.mtx",
	};
	
	int nthreads[] = { 1, 2, 4, 6, 8, 10, 12 };
	reorder_algorithm algorithms[] = { hsl_rcm, boost_rcm, unordered_rcm, shrinked_rcm, bucket_rcm, hsl_sloan, boost_sloan, logbag_sloan, relaxed_order_sloan };
	
	int num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	int size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	int num_algorithms    = sizeof(algorithms)/sizeof(algorithms[0]);
	
	// Running tests
	execute_test_suite(matrices, nthreads, algorithms, num_executions, 
		      num_matrices, size_set_nthreads, num_algorithms);
}

// t = 10
void run_sbpo2017_tests()
{
	int num_executions  = 5;
	
	const char* matrices[] = {
		"../Big-Matrices/20-m_t1.mtx",
		"../Big-Matrices/19-filter3D.mtx",
		"../Big-Matrices/18-SiO2.mtx",
		"../Big-Matrices/17-d_pretok.mtx",
		"../Big-Matrices/16-CO.mtx",
		"../Big-Matrices/15-offshore.mtx",
		"../Big-Matrices/14-Ga41As41H72.mtx",
		"../Big-Matrices/13-F1.mtx",
		"../Big-Matrices/12-mario002.mtx",
		"../Big-Matrices/11-msdoor.mtx",
		"../Big-Matrices/10-inline_1.mtx",
		"../Big-Matrices/09-gsm_106857.mtx",
		"../Big-Matrices/08-Fault_639.mtx",
		"../Big-Matrices/07-tmt_sym.mtx",
		"../Big-Matrices/06-boneS10.mtx",
		"../Big-Matrices/05-audikw_1.mtx",
		"../Big-Matrices/04-nlpkkt80.mtx",
		"../Big-Matrices/03-dielFilterV2real.mtx",
		"../Big-Matrices/02-Serena.mtx",
		"../Big-Matrices/01-G3_circuit.mtx",
	};
	
	int nthreads[] = { 1, 2, 4, 6, 8, 10, 12 };
	
	reorder_algorithm algorithms[] = { boost_sloan, logbag_sloan };
	
	int num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	int size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	int num_algorithms    = sizeof(algorithms)/sizeof(algorithms[0]);
	
	// Running tests
	execute_test_suite(matrices, nthreads, algorithms, num_executions, 
		      num_matrices, size_set_nthreads, num_algorithms);
}


// t = 11
void run_ijca2017_tests()
{
	int num_executions  = 5;
	
	const char* matrices[] = {
		"../Big-Matrices/inline_1.mtx",
		"../Big-Matrices/audikw_1.mtx",
		"../Big-Matrices/dielFilterV3real.mtx",
		"../Big-Matrices/G3_circuit.mtx",
		"../Big-Matrices/333SP.mtx",
	};
	
	int nthreads[] = { 1, 2, 4, 6, 8, 10, 12 };
	
	reorder_algorithm algorithms[] = { hsl_sloan, relaxed_order_sloan };
	
	int num_matrices      = sizeof(matrices)/sizeof(matrices[0]);
	int size_set_nthreads = sizeof(nthreads)/sizeof(nthreads[0]);
	int num_algorithms    = sizeof(algorithms)/sizeof(algorithms[0]);
	
	// Running tests
	execute_test_suite(matrices, nthreads, algorithms, num_executions, 
		      num_matrices, size_set_nthreads, num_algorithms);
}

