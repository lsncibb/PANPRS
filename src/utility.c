/*
 *  utility.c
 *
 *  Created by Wei Sun on Fri Feb 15 2008.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"

/**********************************************************************
 * 
 * reorg
 *
 * Reorganize a vector to a matrix of given size. 
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/

void reorg(double *v, double ***m, int nrow, int ncol)
{
    int i;
    
    *m = (double **)R_alloc(nrow, sizeof(double*));
    
    (*m)[0] = v;
    if(nrow>1){
        for(i=1; i<nrow; i++){
            (*m)[i] = (*m)[i-1] + ncol;
        }
    }
}

void reorg_int(int *v, int ***m, int nrow, int ncol)
{
    int i;
    
    *m = (int **)R_alloc(nrow, sizeof(int*));
    
    (*m)[0] = v;
    if(nrow>1){
        for(i=1; i<nrow; i++){
            (*m)[i] = (*m)[i-1] + ncol;
        }
    }
}


/**********************************************************************
 * 
 * mean
 *
 * calculate the mean value of a vector
 *
 **********************************************************************/
double mean(double *v, int size){
  int i;
  double ave = 0.0;
  for(i=0; i<size; i++){ ave += v[i]; }
  ave /= size;
  return(ave);
}

double mean_j(double **v, int j, int size){
  int i;
  double ave = 0.0;
  for(i=0; i<size; i++){ ave += v[i][j]; }
  ave /= size;
  return(ave);
}

/**********************************************************************
 * 
 * var
 *
 * calculate the variance of a vector
 *
 **********************************************************************/
double var(double *v, int size){
  int i;
  double ave, ave2;
  ave = ave2 = 0.0;
  for(i=0; i<size; i++){ 
    ave += v[i]; 
    ave2 += v[i]*v[i];
  }
  ave /= size;
  ave2 /= size;
  return(ave2 - ave*ave);
}

/**********************************************************************
 * 
 * rsample
 *
 * generate a permutaion from 0 to n-1 using Knuth shuffle
 * http://en.wikipedia.org/wiki/Knuth_shuffle 
 *
 **********************************************************************/

void rsample(int* per, int n) {
  int i, j, v;
  double rU;
  
  for (i=0; i<n; i++) {
    per[i] = i;
  }
 
  for (i=n-1; i>0; i--) {
    rU = runif(0,1);
    j = floor(rU*(i+1)); /* a integer from 0 to i */
    v = per[i];
    per[i] = per[j];
    per[j] = v;
  }
}

/**********************************************************************
 * 
 * rinvGauss
 *
 * generate random numbers from invers Gaussian distribution
 *
 **********************************************************************/
 
double rinvGauss1(double mu, double lambda){
 	double x1, y, z;
 	y  = rchisq(1);
 	x1 = (mu/(2*lambda))*(2*lambda + mu*y - sqrt(4*lambda*mu*y + mu*mu*y*y));
 	z  = runif(0,1);
 	// Rprintf("y=%f, x1=%f, z=%f\n", y, x1, z);
 	if(z < mu/(mu + x1)){
 		return(x1);
 	}else{
 	  return(mu*mu/x1);
  }
}

void rinvGauss(double *x, int* nR, double* mu, double* lambda){
  int i, n;
  n = *nR;
  
  GetRNGstate();

  for(i=0; i < n; i++){
    x[i] = rinvGauss1(*mu, *lambda);
  }
  PutRNGstate();
}

/**********************************************************************
 * 
 * dinvGauss
 *
 * density function for invers Gaussian distribution
 *
 **********************************************************************/
 
void dinvGauss(double* y, double* x, int* n, double* Rmu, double* Rlambda){
  int i;
  double mu, lambda;
  mu = *Rmu;
  lambda = *Rlambda;
  
  for(i=0; i<*n; i++){
  	y[i] = exp(-lambda*(x[i] - mu)*(x[i] - mu)/(2*mu*mu*x[i]));
  	y[i] *= sqrt(lambda/6.283187)*pow(x[i], -1.5);
  }
}


/**********************************************************************
 *
 * readtext (double**, char*, int, int, int, int, int) 
 * 
 * Read the data from the text file
 * The text file must be separated by "whitespace" 
 * The number of rows and columns have to be provided.
 * The element cannot be empty
 * The elements are checked to see whether it is numerical or not
 * Each elements has the maximum 255 characters
 **********************************************************************/

void readtext(double **matrix, char *str, int nrows, int ncols, 
  int offsetrow, int offsetcol, int transpose)
{
  FILE *file;
  int i,j;
  char temp[255];
  file = fopen(str,"r+t");
  if(file == NULL){
    error("fail to open file %s\n", str);
  }
  Rprintf("start reading file\n");
  
  if (transpose==0) {
    for (i=-offsetrow; i<nrows; i++) {
      for (j=-offsetcol; j<ncols; j++) {
        fscanf(file,"%s",temp);
        if(i>=0 && j>=0){
          matrix[i][j] = atof(temp);
        }
      }
    }
  } else {
    // Need to transpose the matrix when loading the file.
    for (j=-offsetcol; j<ncols; j++) {
      for (i=-offsetrow; i<nrows; i++) {
        fscanf(file,"%s",temp);
        if(i>=0 && j>=0){
          matrix[i][j] = atof(temp);
        }
      }
    }
  }
  fclose(file);
}


/**********************************************************************
 *
 * determinant (double** A, int n, double* det) 
 *
 * clacluate the determinant of square matrix A by LU decomposition
 * n is the dimension of the matrix, and det is the determinant
 **********************************************************************/

void determinant(double** A, int n, double* det){
  
  int i, j, k, kb, imax, s1, s2;
  double tmp, pivot;

  /**
   * use det to record the number of row exchange. 
   * det = 1.0 if there are even number of exchagnes and 
   * det = -1.0 if there are odd number of exchanges.
   */
  *det = 1.0;
    
  for(j=0; j<n; j++){
  
    pivot = 0.0;
    imax  = j;
    
    for(i=0; i<n; i++){
      tmp = A[i][j];
      kb  = i; if(i >= j) kb = j;
      for(k=0; k<kb; k++) tmp -= A[i][k]*A[k][j];
      A[i][j] = tmp;
      
      /**
       * If i >= j, consider to exchange rows from here. 
       */
      if(i>=j){
        tmp = fabs(tmp);
        if(tmp >= pivot){ pivot = tmp; imax  = i; }
      }
    }
    
    /**
     * if the pivot is not the current diagnal one, exchange rows
     * well, one easier way is to exchange row ponters
     * however, this will chagne the matrix A in the upper level function
     * if A[0] is switched to somewhere else, it may cause problem
     * when you free the matrix
     */
    if(imax != j){
      for(k=0; k<n; k++){
        tmp = A[j][k];
        A[j][k] = A[imax][k];
        A[imax][k] = tmp;
      }
      (*det)  = -(*det);
    }
    
    if(fabs(A[j][j]) < 1e-16){
      Rprintf("j=%d\n", j);
      for(s1=0; s1<n; s1++){
        for(s2=0; s2<n; s2++){
          printf("%f ", A[s1][s2]);
        }
        printf("\n");
      }
      error("singular matrix in determinant calculation\n");
    }
    
    tmp = A[j][j];
    for(i=j+1; i<n; i++) A[i][j] /= tmp;
  }
  
  for(j=0; j<n; j++) (*det) *= A[j][j];
}

void detR(double* RA, int* n, double* det){
  double **A;
  reorg(RA, &A, *n, *n);
  determinant(A, *n, det);
}

/**********************************************************************
 * 
 * print_v
 *
 * print out a vector of type double
 *
 **********************************************************************/

void print_v(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		printf ("%f\t", v[i]);
	}
	printf ("%f\n", v[i]);
}

void Rprint_v(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%f\t", v[i]);
	}
	Rprintf ("%f\n", v[i]);
}

void Rprint_ve(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%.2e\t", v[i]);
	}
	Rprintf ("%.2e\n", v[i]);
}

void Rprint_vi(int* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%d\t", v[i]);
	}
	Rprintf ("%d\n", v[i]);
}

/**********************************************************************
 * 
 * print_me
 *
 * print out a matrix of format %e
 *
 **********************************************************************/

void print_me(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
    for(j = ncl; j <= nch; j++){
      printf ("%.2e\t", m[i][j]);
    }
    printf("\n");
	}
}

void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
    for(j = ncl; j <= nch; j++){
      Rprintf ("%.2e\t", m[i][j]);
    }
    Rprintf("\n");
	}
}

void Rprint_mi(int** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
    for(j = ncl; j <= nch; j++){
      Rprintf ("%i\t", m[i][j]);
    }
    Rprintf("\n");
	}
}

/**********************************************************************
 * 
 * print_mf
 *
 * print out a matrix of format %f
 *
 **********************************************************************/

void print_mf(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
    for(j = ncl; j <= nch; j++){
      printf ("%f\t", m[i][j]);
    }
    printf("\n");
	}
}

void Rprint_mf(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
    for(j = ncl; j <= nch; j++){
      Rprintf ("%f\t", m[i][j]);
    }
    Rprintf("\n");
	}
}

