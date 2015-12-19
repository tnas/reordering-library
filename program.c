#include "./CommonFiles/protos_parallel.h"




int main (int argc, char* argv[]){
  
	if (argc != 2)
	{
		printf("\n Error! No matrix file (.mtx) specified."); 
		printf("\n Usage: ./program <filename>. Exiting... [main]\n\n");
		return 0;
	}
	
	long int bandwidth,envelope;
	int* permutation;
	double time;
	MAT* matrix;
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, argv[1]);
	
	/*---------------------------------------------*/
	/*---REORDERING--------------------------------*/
	/*---------------------------------------------*/
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);
	printf("\n  [ REORDERING ]\n"); 
	printf("  - Matrix : %s\n", argv[1]); 
	printf("  - Bandwidth/Envelope before = %ld    / %ld\n\n", bandwidth,envelope); 
	
	
	printf("REORDERING SERIAL RCM\n"); fflush(stdout);
	/*---START TIME---------------> */ time = get_time(); 
	REORDERING_RCM_opt(matrix, &permutation);
	MATRIX_permutation(matrix, permutation); 
	/*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;
	
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope); 
	printf("  - Elapsed time: %.6f sec\n\n", time); 
	
	free(permutation);
	MATRIX_clean(matrix);
	
	matrix = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (matrix, argv[1]);
	
	printf("REORDERING PARALLEL UNORDERED RCM\n"); fflush(stdout);
	/*---START TIME---------------> */ time = get_time(); 
	Unordered_RCM(matrix, &permutation);
	MATRIX_permutation(matrix, permutation); 
	/*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;
	
	bandwidth = MATRIX_bandwidth(matrix);
	envelope  = MATRIX_envelope(matrix);	
	printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope);
	printf("  - Elapsed time: %.6f sec\n\n", time);
	
	
	free(permutation);
	MATRIX_clean(matrix);    
	
	return 0;
}