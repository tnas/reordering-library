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
#include <omp.h>
#include <math.h>
#include "util.h"

#define iseven(n) ((n)%(2)==(0)?(1):(0))
#define isdivisor(d, n) ((n)%(d)==(0)?(1):(0))

typedef struct 
{
	int initial_prefix_sum;
	int curr_prefix_sum;
	int curr_total_sum;
	int last_prefix_sum;
	int last_total_sum;
} status_prefix_sum;


void  prefix_sum(const int* counts, int** sums, const int max_level);


