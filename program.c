#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "./CommonFiles/protos.h"


double get_time ()
{
	struct timeval tv; gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}

int main (int argc, char* argv[]){
  
	if (argc != 2)
	{
		printf("\n Error! No matrix file (.mtx) and/or number of machines specified."); 
		printf("\n Usage: ./program <filename>. Exiting... [main]\n\n");
		return 0;
	}
	
	long int bandwidth,envelope;
	int    *p;
	double time;
	
	MAT *A = (MAT*) malloc (sizeof(MAT));
	MATRIX_readCSR (A,argv[1]);
	
	/*---------------------------------------------*/
	/*---REORDERING--------------------------------*/
	/*---------------------------------------------*/
	bandwidth = MATRIX_bandwidth(A);
	envelope  = MATRIX_envelope(A);
	printf("\n  [ REORDERING ]\n");
	printf("  - Bandwidth/Envelope before = %ld    / %ld\n", bandwidth,envelope);
	
	/*---START TIME---------------> */ time = get_time(); 
// 	UNORDERED_RCM(A, &p);
// 	REORDERING_SPECTRAL(A,&p);
	REORDERING_RCM_opt(A,&p);
// 	REORDERING_SLOAN(A,&p);
	MATRIX_permutation(A,p); 
	/*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;
	
	bandwidth = MATRIX_bandwidth(A);
	envelope  = MATRIX_envelope(A);	
	printf("  - Bandwidth/Envelope after  = %ld    / %ld\n", bandwidth,envelope);
	printf("  - Elapsed time: %.6f sec\n\n", time);

	/*---------------------------------------------*/
	/*---WRITING OUTPUT----------------------------*/
	/*---------------------------------------------*/
	int i,j;
	
	FILE *f = fopen("output.mtx", "w");
	fprintf(f,"%d %d %d\n", A->n, A->m, A->nz);
	for (i = 0; i < A->n; ++i)
	{
		for (j = A->IA[i]; j <= A->IA[i+1]-1; j++)
			fprintf (f,"%d %d %le\n",A->JA[j],i,A->AA[j]); 
	}
	fclose(f);
	FILE *g = fopen("plot.plt", "w");
	fprintf(g,"#!/bin/sh\n");
	fprintf(g,"set xrange [0:%d]\n", A->n);
	fprintf(g,"set yrange [0:%d] reverse\n", A->n);
	fprintf(g,"set size square\n");
	fprintf(g,"set terminal png\n");
	fprintf(g,"set output \"output.png\"\n");
	fprintf(g,"plot \"output.mtx\" using 1:2 notitle with dots linecolor rgb \"#000000\"\n");
	fclose(g);
	
	free(p);
	MATRIX_clean(A);    
	return 0;
}