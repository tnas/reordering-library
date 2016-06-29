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