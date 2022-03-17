#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "funcs.h"
#include <omp.h>
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int main(int argc, char *argv[]){
    double start, end;
    int i;
    
    int N=atoi(argv[1]);
    int nThreads=atoi(argv[2]);
    omp_set_nested(1);
    double min=1;
    double max=2;
    double ** G1 = malloc(N*sizeof(double *));
    double ** G2 = malloc(N*sizeof(double *));
    double ** G3;
    for(i=0; i< N; i++){
        G1[i] = malloc(N*sizeof(double));
        G2[i] = malloc(N*sizeof(double));
    }
    rand_init(N,G1,min,max);
    rand_init(N,G2,min-3.0,max);

    /* start clock */
    start = get_wall_seconds();
    G3 = strassen(G1,G2,N,nThreads);
    end = get_wall_seconds()-start;
    printf("Strassen: Nthreads: %d, N: %d, Time: %lf\n", nThreads, N, end);
/*
    double **ref;
   int j;
    start = get_wall_seconds();
    ref=matrix_mult(G1,G2,N);
    end = get_wall_seconds()-start;
    printf("Normal: N: %d, Time: %lf\n", N,end);
    //printf("Reference:\n");
    double eps=0.01;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            if(ref[i][j]>G3[i][j]+eps){
                if(ref[i][j]<G3[i][j]-eps)
                printf("Matrices do not match\n");
            };
            //printf("%lf ",ref[i][j]);
        }
        //printf("\n");
    }
    */
    //Free the original matrices after running
    
    for(i=0;i<N;i++){
        free(G1[i]);
        free(G2[i]);
        free(G3[i]);
        //free(ref[i]);
    }
    free(G1);
    free(G2);
    free(G3);
    //free(ref);
    



    return 0;
}