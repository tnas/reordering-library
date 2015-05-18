#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>

typedef struct {
  int     coord1;
  int     coord2;
  double  coord3;
} Vet;

typedef struct {
  double*     AA;
  int*        JA;
  int*        IA;
  int       n,nz;
} Mat;


/* Leitura matrix market */
int     mm_read_mtx_crd_size      (FILE *f, int *M, int *N, int *nz ){
  char line[1025];
  int num_items_read;

  /* set return null parameter values, in case we exit with errors */
  *M = *N = *nz = 0;

  /* now continue scanning until you reach the end-of-comments */
  do 
    {
      if (fgets(line,1025,f) == NULL) return 12;
    } 
    while (line[0] == '%');

  /* line[] is either blank or has M,N, nz */
  if (sscanf(line, "%d %d %d", M, N, nz) == 3) return 0;
  else
  do
    { 
      num_items_read = fscanf(f, "%d %d %d", M, N, nz); 
      if (num_items_read == EOF) return 12;
    }
    while (num_items_read != 3);
  return 0;
}

void    mat_read                  (Mat* A, char* p){

  int     code,M,N,nz,temp2=0;
  int     i,j,k,temp1=0;
  FILE    *f;
  
  if ((f = fopen(p, "r")) == NULL) 
  exit(1);
  
  /* encontra o tamanho da matriz esparsa .... */
  if ((code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
  exit(1);
  
  A->n  = N;
  A->nz = nz;

  /* reservando memoria para as matrizes */ 
  A->AA  =    (double *) malloc((nz+1)*sizeof (double));
  A->JA  =    (int    *) malloc((nz+1)*sizeof (int));
  A->IA  =    (int    *) malloc((N+2)*sizeof (int));

  /* criando AA, JA e IA */ 
  j=0;k=1;
  for (i=1; i<=nz; i++)
  {
    fscanf(f, "%d %d %lf\n", &A->JA[i], &temp1, &A->AA[i]);

    if (i==1)
    { 
      A->IA[k] = temp1; k++;
    }
    else if (temp2!=temp1)
    { 
      A->IA[k] = A->IA[k-1]+j; k++; j=0;
    }
    temp2 = temp1; j++;       
  }
  A->IA[k] = nz+1;
  /* fechando o arquivo */
  fclose(f);
}

int     vet_compare               (const void *a, const void *b){ 
  if (((Vet*)a)->coord1 <  ((Vet*)b)->coord1) return -1;
  if (((Vet*)a)->coord1 >  ((Vet*)b)->coord1) return  1;
  if (((Vet*)a)->coord1 == ((Vet*)b)->coord1)
  {
    if (((Vet*)a)->coord2 < ((Vet*)b)->coord2) return -1;
    if (((Vet*)a)->coord2 > ((Vet*)b)->coord2) return  1;
  }
  return 0;
}

int main (int argc, char *argv[]){
   
  if (argc < 3)
  {
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
    exit(1);
  }
  
  int     code,M,N,nz;
  int     i,j;
  FILE    *f;

  int choose = atoi(argv[2]);
  
  if (choose == 1)
  {
    
    if ((f = fopen(argv[1], "r")) == NULL) 
    exit(1);
    
    /* encontra o tamanho da matriz esparsa .... */
    if ((code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    exit(1);

    Vet* vet = (Vet*) malloc((2*nz-N)*sizeof(Vet));
    
    for (i=0; i<=(2*nz-N-1); i++)
    {
      fscanf(f,"%d %d %le\n", &vet[i].coord1,&vet[i].coord2,&vet[i].coord3);
      if (vet[i].coord1!=vet[i].coord2) 
      {
	i++;
	vet[i].coord1 = vet[i-1].coord2;
	vet[i].coord2 = vet[i-1].coord1;
	vet[i].coord3 = vet[i-1].coord3;
      }
    }
    /* ordering result */
    qsort(vet, (2*nz-N-1), sizeof(Vet), vet_compare);
    FILE *p = fopen ("new.mtx", "w");
    fprintf (p,"%c%cMatrixMarket matrix coordinate real general\n",37,37);
    fprintf (p,"%d %d %d\n",M,N,2*nz-N);
  
    for (i=0;i<=(2*nz-N-1);i++)
    {
      fprintf(p,"%d %d %.16le\n",vet[i].coord2,vet[i].coord1,vet[i].coord3);
    }
    fclose (p);
  }
  
  if (choose == 2)
  {
    
    if ((f = fopen(argv[1], "r")) == NULL) 
    exit(1);
    
    /* encontra o tamanho da matriz esparsa .... */
    if ((code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    exit(1);
    
    Vet* vet = (Vet*) malloc((2*nz-N)*sizeof(Vet));
    
    for (i=0; i<=(2*nz-N-1); i++)
    {
      fscanf(f,"%d %d\n", &vet[i].coord1,&vet[i].coord2);
      if (vet[i].coord1!=vet[i].coord2) 
      {
	i++;
	vet[i].coord1 = vet[i-1].coord2;
	vet[i].coord2 = vet[i-1].coord1;     
      }
    }
    /* ordering result */
    qsort(vet, (2*nz-N-1), sizeof(Vet), vet_compare);
    FILE *p = fopen ("new.mtx", "w");
    fprintf (p,"%c%cMatrixMarket matrix coordinate real general\n",37,37);
    fprintf (p,"%d %d %d\n",M,N,2*nz-N);
  
    for (i=0;i<=(2*nz-N-1);i++)
    {
      if (vet[i].coord2!=vet[i].coord1)
      fprintf(p,"%d %d 1\n",vet[i].coord2,vet[i].coord1);
      else
      fprintf(p,"%d %d 10\n",vet[i].coord2,vet[i].coord1);
    }
    fclose (p);    
  }
  if (choose == 3)
  {
    
    struct timeb t_start, t_current;
    ftime(&t_start);
    
    Mat* a = (Mat*) malloc(sizeof (Mat));
    /* lendo e criando os vetores AA, JA e IA */
       
    mat_read      (a,argv[1]);
    int n = a->n;
    int nz = a->nz;
    
    FILE *p = fopen ("new.graph", "w");
    fprintf (p," %d %d\n",n,(nz-n)/2);
    
    for (i=1;i<=n;i++){
      for (j=a->IA[i];j<=a->IA[i+1]-1;j++){
	if (i != a->JA[j])
	fprintf (p," %d",a->JA[j]);
      }
      fprintf (p,"\n");
    }
    fclose(p);
   
    ftime(&t_current);
    printf ("TIME: %lf\n", (double) (1000.0 * (t_current.time - t_start.time)
        + (t_current.millitm - t_start.millitm)));
   }
   return 0;
}
  