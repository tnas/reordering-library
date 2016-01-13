#include "./CommonFiles/protos_parallel.h"



void run_tests(char* matrix_file)
{
	double time;
	FILE* out_file;
	FILE* f;
	int exe_times, exe_count;
	MAT* matrix;
	int root, e;
	long int bandwidth, envelope;
	int* permutation;
	
	if ((out_file = fopen("run_tests_output.txt", "w")) == NULL) 
		exit(1);
	
	exe_times = 10;
	root = -1;
	
	fprintf(out_file, "Tests Execution - Matrix: %s - Threads: %d\n", matrix_file, NUM_THREADS);
	printf("Tests Execution - Matrix: %s - Threads: %d\n", matrix_file, NUM_THREADS); fflush(stdout);
	
	for (exe_count = 1; exe_count <= exe_times; ++exe_count)
	{
		fprintf(out_file, "*---------------------------------------------*\n");
		fprintf(out_file, "Execution: %d\n", exe_count);
		fprintf(out_file, "*---------------------------------------------*\n");
		
		printf("*---------------------------------------------*\n");fflush(stdout);
		printf("Execution: %d\n", exe_count);fflush(stdout);
		printf("*---------------------------------------------*\n");fflush(stdout);
		
		if ((f = fopen(matrix_file, "r")) == NULL) 
			exit(1);
	
		matrix = (MAT*) malloc (sizeof(MAT));
		MATRIX_readCSR (matrix, f);
		fclose(f);
		
		if (root == -1) GRAPH_LS_peripheral(matrix, &root, &e);
		
		bandwidth = MATRIX_bandwidth(matrix);
		envelope  = MATRIX_envelope(matrix);
		
		time = get_time(); 
		REORDERING_RCM_opt(matrix, &permutation, root);
		time = (get_time() - time)/100.0;
		
		MATRIX_permutation(matrix, permutation);
		
		bandwidth = MATRIX_bandwidth(matrix);
		envelope  = MATRIX_envelope(matrix);	
		
		free(permutation);
		MATRIX_clean(matrix);
		
		fprintf(out_file, "REORDERING SERIAL RCM\n");
		fprintf(out_file, "  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope); 
		fprintf(out_file, "  - Elapsed time: %.6f sec\n", time);
		
		printf("REORDERING SERIAL RCM\n");fflush(stdout);
		printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope); fflush(stdout);
		printf("  - Elapsed time: %.6f sec\n", time);fflush(stdout);
		
		if ((f = fopen(matrix_file, "r")) == NULL) 
			exit(1);
	
		matrix = (MAT*) malloc (sizeof(MAT));
		MATRIX_readCSR (matrix, f);
		fclose(f);
		
		bandwidth = MATRIX_bandwidth(matrix);
		envelope  = MATRIX_envelope(matrix);
		
		time = get_time(); 
		Unordered_RCM(matrix, &permutation, root);
		time = (get_time() - time)/100.0;
		
		MATRIX_permutation(matrix, permutation);
		
		bandwidth = MATRIX_bandwidth(matrix);
		envelope  = MATRIX_envelope(matrix);	
		
		free(permutation);
		MATRIX_clean(matrix);
		
		fprintf(out_file, "REORDERING PARALLEL UNORDERED RCM\n");
		fprintf(out_file, "  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope);
		fprintf(out_file, "  - Elapsed time: %.6f sec\n\n", time);
		
		printf("REORDERING PARALLEL UNORDERED RCM\n");fflush(stdout);
		printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope);fflush(stdout);
		printf("  - Elapsed time: %.6f sec\n\n", time);fflush(stdout);
	}
	
	fclose(out_file);
}



int main (int argc, char* argv[]){
  
	if (argc != 2)
	{
		printf("\n Error! No matrix file (.mtx) specified."); 
		printf("\n Usage: ./program <filename>. Exiting... [main]\n\n");
		return 0;
	}
	
	if (RUN_AUTOMATIC_TESTS)
	{
		run_tests(argv[1]);
		exit(0);
	}
	
	int root, e;
	long int bandwidth,envelope;
	int* permutation;
	double time, time_peripheral;
	MAT* matrix;
	FILE* f;
	
	if ((f = fopen(argv[1], "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, f);
	fclose(f);
	
	/*---------------------------------------------*/
	/*------------PSEUDO PERIPHERAL----------------*/
	/*---------------------------------------------*/
	time_peripheral = get_time(); 
	GRAPH_LS_peripheral(matrix, &root, &e);
	time_peripheral = (get_time() - time_peripheral)/100.0;
	
	/*---------------------------------------------*/
	/*---REORDERING--------------------------------*/
	/*---------------------------------------------*/
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);
	printf("\n  [ REORDERING ]\n"); 
	printf("  - Matrix : %s\n", argv[1]); 
	printf("  - Bandwidth/Envelope before = %ld    / %ld\n\n", bandwidth,envelope); 
	
	/*--------------SERIAL RCM-------------------------*/
	
	printf("REORDERING SERIAL RCM\n"); fflush(stdout);
	time = get_time(); 
	REORDERING_RCM_opt(matrix, &permutation, root);
	time = (get_time() - time)/100.0;
	
	MATRIX_permutation(matrix, permutation);
	
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope); 
	printf("  - Elapsed time: %.6f sec\n\n", time); 
	
	free(permutation);
	MATRIX_clean(matrix);
	
	/*--------------UNORDERED RCM-------------------------*/
	
	if ((f = fopen(argv[1], "r")) == NULL) 
		exit(1);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, f);
	fclose(f);
	
	printf("REORDERING PARALLEL UNORDERED RCM\n"); fflush(stdout);
	time = get_time(); 
	Unordered_RCM(matrix, &permutation, root);
	time = (get_time() - time)/100.0;
	
	MATRIX_permutation(matrix, permutation); 
	
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope);
	printf("  - Elapsed time: %.6f sec\n\n", time);
	
	free(permutation);
	MATRIX_clean(matrix);    
	
	return 0;
}