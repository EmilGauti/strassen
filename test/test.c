#include <stdlib.h>
#include <stdio.h>

double **strassen(int N){
    
    double ** G = malloc(N*sizeof(double *));
    for(int i=0; i< N; i++){
        G[i] = malloc(N*sizeof(double));
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            G[i][j]=i+j;
        }
    }
    return G;
}
double** mat_add(double **G1, double **G2, int N){
    double ** sum = malloc(N*sizeof(double *));
    for(int i=0; i< N; i++){
        sum[i] = malloc(N*sizeof(double));
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sum[i][j]=G1[i][j]+G2[i][j];
        }
    }
    return sum;
}
double** mat_sub(double **G1, double **G2, int N){
    double ** sum = malloc(N*sizeof(double *));
    for(int i=0; i< N; i++){
        sum[i] = malloc(N*sizeof(double));
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sum[i][j]=G1[i][j]-G2[i][j];
        }
    }
    return sum;
}

int main(){
    double **G3, **G3_add, **G3_sub;
    int N=4;
    G3=strassen(N);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%lf ",G3[i][j]);
        }
        printf("\n");
    }
    G3_add=mat_add(G3,G3,N);
    printf("add\n");
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%lf ",G3_add[i][j]);
        }
        printf("\n");
    }    
    G3_sub=mat_sub(G3_add,G3,N);
    printf("sub\n");
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%lf ",G3_sub[i][j]);
        }
        printf("\n");
    }
}