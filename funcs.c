#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "funcs.h"
#include <omp.h>
double **matrix_mult(double **G1, double **G2, int N){
    int i,j,k;
    double sum;
    double ** prod = malloc(N*sizeof(double *));
    for(i=0; i< N; i++){
        prod[i] = malloc(N*sizeof(double));
    }
    for(i=0; i<N;i++){
        for(j=0;j<N;j++){
            sum=0;
            for(k=0;k<N;k=k+4){
                sum=sum+G1[i][k]*G2[k][j];
                sum=sum+G1[i][k+1]*G2[k+1][j];
                sum=sum+G1[i][k+2]*G2[k+2][j];
                sum=sum+G1[i][k+3]*G2[k+3][j];

            }
            prod[i][j]=sum;
        }
    }
    return prod;
}

void split(mSplit_t *subMat, double** G, int N){
    /*
    In order to not have to create new submatrices I opted for a structure
    that stores the borders of each and every submatrix. This way we can operate
    on the original matrix.
    */
    //subMat->N=N;
    //subMat->mat=G;
    int i,j;
    double ** nw = malloc(N/2*sizeof(double *));
    double ** ne = malloc(N/2*sizeof(double *));
    double ** sw = malloc(N/2*sizeof(double *));
    double ** se = malloc(N/2*sizeof(double *));
    for(i=0; i< N/2; i++){
        nw[i] = malloc(N/2*sizeof(double));
        ne[i] = malloc(N/2*sizeof(double));
        sw[i] = malloc(N/2*sizeof(double));
        se[i] = malloc(N/2*sizeof(double));
    }

    for(i=0;i<N/2;i++){
        for(j=0;j<N/2;j++){
            nw[i][j]=G[i][j];
            ne[i][j]=G[i][j+N/2];
            sw[i][j]=G[i+N/2][j];
            se[i][j]=G[i+N/2][j+N/2];
        }
    }
    subMat->nw=nw;
    subMat->ne=ne;   
    subMat->sw=sw;
    subMat->se=se;
}
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
void rand_init(int N, double **G,double min, double max){
    srand(time(NULL));
    int i,j;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            //if(i==j)
            G[i][j]=randfrom(min,max);
        }
    }
}
double** mat_add(double **G1, double **G2, int N){
    int i,j;
    double ** sum = malloc(N*sizeof(double *));
    for(i=0; i< N; i++){
        sum[i] = malloc(N*sizeof(double));
    }
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            sum[i][j]=G1[i][j]+G2[i][j];
        }
    }
    return sum;
}
double** mat_sub(double **G1, double **G2, int N){
    int i,j;
    double ** sum = malloc(N*sizeof(double *));
    for(i=0; i< N; i++){
        sum[i] = malloc(N*sizeof(double));
    }
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            sum[i][j]=G1[i][j]-G2[i][j];
        }
    }
    return sum;
}
double** stackP(double** p1, double** p2, double** p3, double** p4, double** p5, double** p6, double** p7, int N){
    int i,j;
    double ** c = malloc(N*sizeof(double *));
    for(i=0; i< N; i++){
        c[i] = malloc(N*sizeof(double));
    }
    double **c11_1, **c11_2;
    double **c22_1, **c22_2;
    double **c11, **c12, **c21, **c22;
    //c11 = mat_sub(mat_add(p5, p4,N/2), mat_add(p2, p6, N/2),N/2); //p5+p4-p2+p6 
    c11_1=mat_sub(p4,p2,N/2);
    c11_2=mat_add(p5,c11_1,N/2);
    c11= mat_add(c11_2,p6,N/2);
    //c11 = mat_add(mat_add(p5,mat_sub(p4,p2,N/2),N/2),p6,N/2);//p5+p4-p2+p6 =((p5+p4)-(p2-p6)) = ((p5+(p4-p2))+p6)
    c12 = mat_add(p1, p2,N/2);//p1+p2
    c21 = mat_add(p3, p4,N/2);//p3+p4
    //c22 = mat_sub(mat_add(p1, p5, N/2), mat_sub(p3, p7,N/2),N/2);//p1+p5-p3-p7
    c22_1=mat_add(p1,p5,N/2);
    c22_2=mat_add(p3,p7,N/2);
    c22 = mat_sub(c22_1,c22_2,N/2);//p1+p5-p3-p7=((p1+p5)-(p3+p7))
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            if(i<N/2){
                if(j<N/2){
                    c[i][j]=c11[i][j];
                }
                else{
                    c[i][j]=c12[i][j-N/2];
                }
            }
            else{
                if(j<N/2){
                    c[i][j]=c21[i-N/2][j];
                }
                else{
                    c[i][j]=c22[i-N/2][j-N/2];
                }
            }
        }
    }
    free_mat(c11_1,N/2);
    free_mat(c11_2,N/2);
    free_mat(c22_1,N/2);
    free_mat(c22_2,N/2);

    free_mat(c11,N/2);
    free_mat(c12,N/2);
    free_mat(c21,N/2);
    free_mat(c22,N/2);
    return c;
}
void free_p(double** p1, double** p2, double** p3, double** p4, double** p5, double** p6, double** p7, int N){
    int i;
    for(i=0;i<N;i++){
        free(p1[i]);
        free(p2[i]);
        free(p3[i]);
        free(p4[i]);
        free(p5[i]);
        free(p6[i]);
        free(p7[i]);
    }
    free(p1);
    free(p2);
    free(p3);
    free(p4);
    free(p5);
    free(p6);
    free(p7);
}
void free_s(mSplit_t *s, int N){
    int i;
    for(i=0; i< N/2; i++){
        //printf("free_s: i=%d , N=%d\n",i,N/2);
        free(s->nw[i]);
        free(s->ne[i]);
        free(s->sw[i]);
        free(s->se[i]);
    }
    //printf("free_s loop done\n");
    free(s->nw);
    free(s->ne);
    free(s->sw);
    free(s->se);
    free(s);
}
void free_mat(double **G, int N){
    int i;
    for(i=0; i<N;i++){
        //printf("free_mat: i=%d , N=%d\n",i,N);
        free(G[i]);
    }
    //printf("Free_mat loop done\n");
    free(G);
    //printf("Free_mat done\n");
}
void free_sub_p(double **p1_2, double **p2_1, double **p3_1, double **p4_2, double **p5_1, double **p5_2, double **p6_1, double **p6_2, double **p7_1, double **p7_2, int N){
    free_mat(p1_2,N);
    free_mat(p2_1,N);
    free_mat(p3_1,N);
    free_mat(p4_2,N);
    free_mat(p5_1,N);
    free_mat(p5_2,N);
    free_mat(p6_1,N);
    free_mat(p6_2,N);
    free_mat(p7_1,N);
    free_mat(p7_2,N);

}

double **strassen(double** G1, double** G2, int N, int nThreads){

    /*
    if(N==1){
        //This is never used if the below is used
        //Testing showed it to be faster
        //show in report....
        double ** sum = malloc(sizeof(double *));
        sum[0] = malloc(sizeof(double));
        sum[0][0]=G1[0][0]*G2[0][0];
        return sum;
    }
    */
    if(N<=64){
        //N=64 is optimal
        //more testing required
        //printf("nThreads=%d\n",nThreads);
        return matrix_mult(G1,G2,N);
    }
    double **p1,**p2,**p3,**p4,**p5,**p6,**p7;
    double **c;
    mSplit_t *s1;
    mSplit_t *s2;
    s1=(mSplit_t *)malloc(sizeof(mSplit_t));
    s2=(mSplit_t *)malloc(sizeof(mSplit_t));
    split(s1,G1,N);
    split(s2,G2,N);
    
    
    
    //               a b 
    //               c d
    //               e f 
    //               g h
    double **p1_2, **p2_1, **p3_1, **p4_2, **p5_1, **p5_2, **p6_1, **p6_2, **p7_1, **p7_2;
    
    
    
    
    
    
    int nr_threads_send=nThreads/7;
    int N_2=N/2; 
    if(nThreads<=0){
        nThreads=1;
    }
    #pragma omp parallel num_threads(nThreads)
    {
        #pragma omp single nowait
        {
            #pragma omp task
            {
                p1_2=mat_sub(s2->ne,s2->se,N_2);
		        p1 = strassen(s1->nw,p1_2,N_2,nr_threads_send);//(a, f - h); x
            }
            #pragma omp task
            {
                p2_1=mat_add(s1->nw,s1->ne,N_2);
                p2 = strassen(p2_1,s2->se,N_2,nr_threads_send);//(a + b, h); x
            }
            #pragma omp task
            {
                p3_1=mat_add(s1->sw,s1->se,N_2);
                p3 = strassen(p3_1,s2->nw,N_2,nr_threads_send);//(c + d, e); x
            }
            #pragma omp task
            {
                p4_2=mat_sub(s2->sw,s2->nw,N_2);
                p4 = strassen(s1->se,p4_2,N_2,nr_threads_send);//(d, g - e); x
            }
            #pragma omp task
            {
                p5_1=mat_add(s1->nw,s1->se,N_2);
                p5_2=mat_add(s2->nw,s2->se,N_2);
                p5 = strassen(p5_1,p5_2,N_2,nr_threads_send);//(a + d, e + h); x x
            }
            #pragma omp task
            {
                p6_1=mat_sub(s1->ne,s1->se,N_2);
                p6_2=mat_add(s2->sw,s2->se,N_2);
                p6 = strassen(p6_1,p6_2,N_2,nr_threads_send);//(b - d, g + h); x x
            }
            #pragma omp task
            {
                p7_1=mat_sub(s1->nw,s1->sw,N_2);
                p7_2=mat_add(s2->nw,s2->ne,N_2);
                p7 = strassen(p7_1,p7_2,N_2,nr_threads_send);//(a - c, e + f); x x
            }
            
        }
        

    }

    
    free_sub_p(p1_2, p2_1, p3_1, p4_2, p5_1, p5_2, p6_1, p6_2, p7_1, p7_2,N_2);
    free_s(s1,N);
    free_s(s2,N);
    c=stackP(p1,p2,p3,p4,p5,p6,p7,N);
    free_p(p1,p2,p3,p4,p5,p6,p7,N_2);
    
    
    
    return c;
}
