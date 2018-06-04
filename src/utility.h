
void reorg(double *v, double ***m, int nrow, int ncol);
void reorg_int(int *v, int ***m, int nrow, int ncol);
double mean(double *v, int size);
double mean_j(double **v, int j, int size);
double var(double *v, int size);
void rsample(int* per, int n);
double rinvGauss1(double mu, double lambda);
void rinvGauss(double *x, int* n, double* mu, double* lambda);
void dinvGauss(double* y, double* x, int* n, double* Rmu, double* Rlambda);
void readtext(double **matrix, char *str, int nrows, int ncols, int offsetrow, int offsetcol, int transpose);
void determinant(double** A, int n, double* det);
void detR(double* RA, int* n, double* det);

void print_v(double* v, int nrl, int nrh);
void Rprint_v(double* v, int nrl, int nrh);
void Rprint_ve(double* v, int nrl, int nrh);
void Rprint_vi(int* v, int nrl, int nrh);

void print_me(double** m, long nrl, long nrh, long ncl, long nch);
void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch);
void Rprint_mi(int** m, long nrl, long nrh, long ncl, long nch);

void print_mf(double** m, long nrl, long nrh, long ncl, long nch);
void Rprint_mf(double** m, long nrl, long nrh, long ncl, long nch);
