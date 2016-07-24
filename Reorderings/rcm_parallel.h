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

#include "../CommonFiles/graph_parallel.h"
#include "../CommonFiles/util_parallel.h"
#include "reordering.h"


void Unordered_RCM(MAT* A, int** perm, int root, const float percent_chunk);
void Unordered_RCM_METAGRAPH(const METAGRAPH* mgraph, int** perm, int root, const float percent_chunk);
void Leveled_RCM(MAT* mat, int** perm, int root);
void Leveled_RCM_METAGRAPH(METAGRAPH* mgraph, int** perm, int root);
void Bucket_RCM(MAT* mat, int** perm, int root);
void Bucket_RCM_METAGRAPH(const METAGRAPH* mgraph, int** perm, int root);