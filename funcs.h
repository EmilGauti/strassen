typedef struct mSplit{
    //int N;
    //double **mat;
    double **nw;
    double **ne;
    double **sw;
    double **se;
}mSplit_t;

double **matrix_mult(double **G1, double **G2, int N);
void split(mSplit_t *subMat, double** G, int N);
double randfrom(double min, double max);
void rand_init(int N, double **G,double min, double max);
double** mat_add(double **G1, double **G2, int N);
double** mat_sub(double **G1, double **G2, int N);
double** stackP(double** p1, double** p2, double** p3, double** p4, double** p5, double** p6, double** p7, int N);
double **strassen(double** G1, double** G2, int N);

void free_p(double** p1, double** p2, double** p3, double** p4, double** p5, double** p6, double** p7, int N);
void free_s(mSplit_t *s, int N);
void free_mat(double **G, int N);
void free_sub_p(double **p1_2, double **p2_1, double **p3_1, double **p4_2, double **p5_1, double **p5_2, double **p6_1, double **p6_2, double **p7_1, double **p7_2, int N);