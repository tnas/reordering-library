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

#include "graph_hsl.h"

int* get_pseudo_diameter_hsl(const MAT* matrix) 
{
	int* peripheral_nodes;
	
	peripheral_nodes = calloc(2, sizeof(int));

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
	
	free(vars);
	free(mask);
	free(ls);
	free(xls);
	free(list);
	
	return peripheral_nodes;
}