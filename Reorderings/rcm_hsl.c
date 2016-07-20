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

#include "rcm_hsl.h"

long int Reordering_RCM_HSL(MAT* A)
{
	int i, nsup, pair_lenght;
	int n        = A->n; 
	int *irn     = A->JA;
	int *icptr   = A->IA;
	int lirn     = A->nz;
	int jcntl[2] = { RCM, AUTOMATIC_PERIPHERAL };
	double weight[2];
	int info[4];
	double rinfo[4];
	int* vars;
	int** pair;
	int* iw;
	double* w;
	int* permsv;
	int* svar;
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	
	for (i = 0; i < n; i++) 
	{
		++irn[i];
		++icptr[i];
	}
	
	for (i = n; i < lirn; i++) ++irn[i];
	
	svar   = calloc(n, sizeof(int));
	vars   = calloc(n, sizeof(int));
	iw     = calloc(2*n + 2, sizeof(int));
	
	// To find supervariables and compress pattern
	mc60bd_(&n, &lirn, irn, icptr, &nsup, svar, vars, iw);
	
	free(iw);
	pair_lenght = nsup/2;
	permsv = calloc(nsup, sizeof(int));
	iw     = calloc(3*nsup + 1, sizeof(int));
	w      = calloc(nsup, sizeof(double));
	pair   = (int**) calloc(pair_lenght, sizeof(int*));
	
	for (i = 0; i < pair_lenght; ++i)
		pair[i] = (int*) calloc(2, sizeof(int));
	
	// To find supervariable permutation
	mc60cd_(&n, &nsup, &lirn, irn, icptr, vars, jcntl, permsv, weight, pair, info, iw, w);
	
	free(iw);
	iw = calloc(2*nsup + 1, sizeof(int));
	
	// To compute the profile and wavefront for a supervariable permutation
	mc60fd_(&n, &nsup, &lirn, irn, icptr, vars, permsv, iw, rinfo);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < n; i++) 
	{
		--irn[i];
		--icptr[i];
	}
	
	for (i = n; i < lirn; i++) --irn[i];
	
	free(vars);
	free(iw);
	free(w);
	free(permsv);
	free(svar);
	free(pair);
	
	return rinfo[SEMI_BANDWIDTH];
}
