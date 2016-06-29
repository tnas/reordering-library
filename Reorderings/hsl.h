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

#ifndef __HSL_H__
#define __HSL_H__

typedef enum MC60_ALGORITHM { 
	SLOAN, RCM 
} MC60_ALGORITHM;

typedef enum MC60_CONTROL { 
	AUTOMATIC_PERIPHERAL, 
	ESPECIFIED_PERIPHERAL, 
	GLOBAL_PRIORITY_VECTOR 
} MC60_CONTROL;

typedef enum MATRIX_PROPERTY { 
	PROFILE, 
	MAX_WAVEFRONT, 
	SEMI_BANDWIDTH, 
	RMS_WAVEFRONT 
} MATRIX_PROPERTY;


void mc60bd_     (int* n, int* lirn, int* irn, int* icptr, int* nsup, int* svar, int* vars, int* iw);
void mc60cd_     (int* n, int* nsup, int* lirn, int* irn, int* icptr, int* vars, int* jcntl, int* permsv, double* weight, int** pair, int* info, int* iw, double* w);
void mc60dd_     (int* n, int* nsup, int* svar, int* vars, int* permsv, int* perm, int* possv);
void mc60fd_     (int* n, int* nsup, int* lirn, int* irn, int* icptr, int* vars, int* permsv, int* iw, double* rinfo);
void mc60hd_	 (int* n, int* nsup, int* lirn, int* irn, int* icptr, int* vars, int* mask, int* ls, int* xls, int* list, int* info);
void mc73_fiedler(int* n, int* lirn, int* irn, int* ip, int* list, double* fvector, int* info, double* a);

#endif