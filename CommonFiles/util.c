#include "util.h"

int get_random_integer(int max)
{
	srand((unsigned int) time(NULL));
	return rand() % (max + 1);
}


int pow_uint(int base, const int exp)
{
	int count, val_pow;
	if (exp == 0) return 1;
	val_pow = base;
	for (count = 1; count < exp; ++count)
		val_pow *= base;
	return val_pow;
}

double get_time ()
{
	struct timeval tv; gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}


/*---------------------------------------------*/
/*---WRITING OUTPUT----------------------------*/
/*---------------------------------------------*/
void write_output(const MAT *A)
{
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
}