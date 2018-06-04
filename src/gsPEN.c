/*
 *  IterativeSNPC.c
 *
 *  Created by Ting-Huei Chen on Feb 06, 2015.
 *
 *  Last updated by Ting-Huei Chen on Feb 25, 2018.
 */
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "gsPEN.h"



/*********************************************************************
 *
 * 
 *
 *********************************************************************/
 
void gsPEN(double* RsummaryBetas,int* ldJ, int* dims, int* Numitervec, 
     int* RIndexMatrix, int* IndJ, 
     double* ldvec, double* inv_summaryBetas, int*ChrIndexBeta, double*RupperVal,
     double* Init_summaryBetas, double* RSDvec, 
     double* RtuningMatrix, double* RBetaMatrix, char **penalty_)
{
    int NumTuning = dims[0], P = dims[1];
    int nrow_IndexMatrix = dims[4], ncol_IndexMatrix = dims[5];
    int **IndexMatrix;
    int ncol_tuningMatrix = dims[9], Q = dims[10], ncolBetaMatrix = dims[11];
    double **BetaMatrix, **tuningMatrix, **summaryBetas, **SDvec;
    
    /* reorganize vector into matrix */
    reorg(RsummaryBetas, &summaryBetas, P, Q);
    reorg(RSDvec, &SDvec, P, Q);
    reorg(RBetaMatrix, &BetaMatrix, NumTuning, ncolBetaMatrix);
    
    reorg(RtuningMatrix, &tuningMatrix, NumTuning, ncol_tuningMatrix);
    reorg_int(RIndexMatrix, &IndexMatrix, nrow_IndexMatrix, ncol_IndexMatrix);
        
    gsPENC(summaryBetas, ldJ, dims, Numitervec,
     IndexMatrix, IndJ,
     ldvec, inv_summaryBetas, ChrIndexBeta, RupperVal,
     Init_summaryBetas, SDvec, 
     tuningMatrix, BetaMatrix,penalty_);
} 


void gsPENC(double** summaryBetas, int* ldJ, int* dims, int* Numitervec, 
     int** IndexMatrix, int* IndJ,
     double* ldvec, double* inv_summaryBetas, int*ChrIndexBeta, double*RupperVal, 
     double* Init_summaryBetas, double** SDvec, 
     double**tuningMatrix, double** BetaMatrix,
     char **penalty_)
{
  
    int NumTuning = dims[0], P = dims[1], NumIter = dims[2], breaking = dims[3];
    int niter, t1, i, j,j1, j2, found, **skipb;
    int NumSNP = dims[4], NumInd = dims[6];
    int Zscale = dims[7];

    int Q = dims[10], k1, cindex;
    int OutputForPath = dims[12], dfMax = dims[13], df_q1;
    char *penalty=penalty_[0];

    
    double tmp0, *tmpbetas, threshold, bj_bar, epsilon = 0.0001;
    double *sumBetas;
    double lambda1, lambda2, tau2, upperVal = RupperVal[0];
    
    
    /* allocate memory */
    tmpbetas = (double *)calloc(P, sizeof(double));
    sumBetas = (double *)calloc(P, sizeof(double));
    
    GetRNGstate();
    
    double **jointBmatrix, **tempBmatrix;
    
    jointBmatrix = (double **) malloc(P * sizeof(double*));
    tempBmatrix = (double **) malloc(P * sizeof(double*));
    skipb = (int **) malloc(P * sizeof(int*));
  
    jointBmatrix[0] = (double *) calloc(P*Q, sizeof(double));
    if(jointBmatrix[0] == NULL){ error("fail to allocate memory of jointBmatrix"); }
    for(j=0; j<P; j++){  
      jointBmatrix[j] = jointBmatrix[0] + j*Q; 
    }
    
    tempBmatrix[0] = (double *) calloc(P*Q, sizeof(double));
    if(tempBmatrix[0] == NULL){ error("fail to allocate memory of tempBmatrix"); }
    for(j=0; j<P; j++){  
      tempBmatrix[j] = tempBmatrix[0] + j*Q; 
    }
    
    skipb[0] = (int *) calloc(P*Q, sizeof(int));
    if(skipb[0] == NULL){ error("fail to allocate memory of skipb"); }
    for(j=0; j<P; j++){  
      skipb[j] = skipb[0] + j*Q; 
    }

    for(i=0; i<P; i++){  
      sumBetas[i] = 0.0;
      for(j=0; j<Q; j++){  
        jointBmatrix[i][j] = 0.0;
        tempBmatrix[i][j] = 0.0;
        skipb[i][j] = 0;
      } 
    }
    
    for(j1=0; j1<NumSNP; j1++){
      j = IndexMatrix[j1][0];
      
      for(k1=0; k1<Q; k1++){
        jointBmatrix[j][k1] = 0.0;
      }  
    }
    

    for(t1=0; t1<NumTuning; t1++){
      Numitervec[t1] = 0;

      if (strcmp(penalty,"mixLOG")==0){
        lambda1 = tuningMatrix[t1][0];
        lambda2 = tuningMatrix[t1][2];
        tau2 = tuningMatrix[t1][3];
      }
      
      
      for(i=0; i<P; i++){  
        tmp0 = 0.0;
        for(j=0; j<Q; j++){  
          jointBmatrix[i][j] = 0.0;
          tempBmatrix[i][j] = 0.0;
          skipb[i][j] = 0;

        }

        sumBetas[i] = tmp0;
      }

      
      cindex = 0;
      for(niter=0; niter<NumIter; niter++){
          
        if(NumInd!=0){
          for(j1=0; j1<NumInd; j1++){
            j = IndJ[j1];
            
            for(k1=0; k1<Q; k1++){
              bj_bar = summaryBetas[j][k1];
              
              if(bj_bar!=0.0){
                if (strcmp(penalty,"mixLOG")==0){
                  threshold = lambda1 + lambda2*(1/(sumBetas[j] + tau2));
                }  
                if(Zscale==1){
                  threshold = threshold*SDvec[j][k1];
                }

                if(bj_bar > threshold){
                  jointBmatrix[j][k1] = bj_bar - threshold;
                }else if(bj_bar < -threshold){
                  jointBmatrix[j][k1] = bj_bar + threshold;
                }else{
                  jointBmatrix[j][k1] = 0.0;
                }
              }
              if(summaryBetas[j][k1]*jointBmatrix[j][k1]<0){
                Rprintf("summaryBetas[j]=%d\n",j);
                Rprintf("summaryBetas[k1]=%d\n",k1);
                Rprintf("summaryBetas[j][k1]=%e\n",summaryBetas[j][k1]);
                Rprintf("jointBmatrix[j][k1]=%e\n",jointBmatrix[j][k1]);
                error("sign inverse");
              }
            }
          } 
        }

        for(j1=0; j1<NumSNP; j1++){
          j = IndexMatrix[j1][0];
          
          for(k1=0; k1<Q; k1++){
          
            if(skipb[j][k1]==0){
              if(summaryBetas[j][k1]!=0.0){
                tmp0 = 0.0;
      
                for(j2=IndexMatrix[j1][1]; j2<(IndexMatrix[j1][2]+1); j2++){
                  tmp0 = tmp0 + ldvec[j2]*jointBmatrix[ldJ[j2]][k1];
                }
          
                bj_bar = (summaryBetas[j][k1] - tmp0);
          
                if (strcmp(penalty,"mixLOG")==0){
                  threshold = lambda1 + lambda2*(1/(sumBetas[j] + tau2));
                }

              
                if(fabs(bj_bar)>upperVal){
                  if(breaking==1){
                    Numitervec[t1] = -1;
                    break;
                  }else{
                    bj_bar = 0.0;
                    skipb[j][k1] = 1;
                  }
                }
            
                if(Zscale==1){
                  threshold = threshold*SDvec[j][k1];
                }

              
                if(bj_bar > threshold){
                  jointBmatrix[j][k1] = bj_bar - threshold;
          
                }else if(bj_bar < -threshold){
                  jointBmatrix[j][k1] = bj_bar + threshold;
                }else{
                  jointBmatrix[j][k1] = 0.0;
                }
              
              }else{
                jointBmatrix[j][k1] = 0.0;
              }  
            }else{
              jointBmatrix[j][k1] = 0.0;
            }
          }
        }
        
       /**
        * check convergence
        */
              
        found = 0;
        df_q1 = 0;
        for(j=0; j<P; j++){
          for(k1=0; k1<Q; k1++){
            if(fabs(jointBmatrix[j][k1])>upperVal){
              skipb[j][k1] = 1;
            }
            if(k1==0){
              if(fabs(jointBmatrix[j][k1])!=0){
                df_q1 = df_q1 + 1;
              }
              if(df_q1 > dfMax){
                cindex = 1;
                break;
              }
            }
            if(fabs(tempBmatrix[j][k1] - jointBmatrix[j][k1]) > epsilon){
              found = 1;
              break;
            }
          } 
        }
        if(cindex == 1){
          Numitervec[t1] = -2;
          break;
        }
        
        if(found==0){
          k1=0;
          for(i=0; i<Q; i++){
            for(j=0; j<P; j++){
              BetaMatrix[t1][k1] = jointBmatrix[j][i];
              k1 = k1 + 1;
            } 
          } 
          Numitervec[t1] = (niter+1);
          break;
        }
        for(j1=0; j1<NumSNP; j1++){
          j = IndexMatrix[j1][0];
          tmp0 = 0.0;
          for(k1=0; k1<Q; k1++){
            tempBmatrix[j][k1] = jointBmatrix[j][k1];
            if(jointBmatrix[j][k1]!=0){
              tmp0 = tmp0 + fabs(jointBmatrix[j][k1]);
            }
          }

          sumBetas[j] = tmp0;
        }

        if(NumInd!=0){
          for(j1=0; j1<NumInd; j1++){
            j = IndJ[j1];
            tmp0 = 0.0;
            for(k1=0; k1<Q; k1++){
              tempBmatrix[j][k1] = jointBmatrix[j][k1];
              if(jointBmatrix[j][k1]!=0){
                tmp0 = tmp0 + fabs(jointBmatrix[j][k1]);
              }
            }
            sumBetas[j] = tmp0;
          }
        }

        
        if(OutputForPath==1 && niter==(NumIter-1)){
          k1=0;
          for(i=0; i<Q; i++){
            for(j=0; j<P; j++){
              BetaMatrix[t1][k1] = jointBmatrix[j][i];
              k1 = k1 + 1;
            } 
          } 
        }
      }
    } 
    PutRNGstate();
    
    Free(jointBmatrix[0]);
    Free(tempBmatrix[0]);
    Free(skipb[0]);
    Free(jointBmatrix);
    Free(tempBmatrix);
    Free(skipb);

    
 }

 
