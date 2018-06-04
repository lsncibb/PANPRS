/*
 *  IterativeSNPC.c
 *
 *  Created by Ting-Huei Chen on Feb 06, 2015.
 *
 *  Last updated by Ting-Huei Chen on Feb 24, 2018.
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
#include "gsfPEN.h"



/*********************************************************************
 *
 *  gsfFunc function
 *
 *********************************************************************/
 


void gsfFunc(double* RsummaryBetas,int* ldJ, int* dims, int* Numitervec, 
     int* RIndexMatrix, int* IndJ, 
     double* ldvec, int*ChrIndexBeta, double*RupperVal,
     double* Init_summaryBetas, double* RSDvec, 
     double* RtuningMatrix, double* RBetaMatrix, char **penalty_,
     double* lambda0vec, double* RZmatrix, int* dims2, 
     double* RAllTuningMatrix, double* lambdavec, int* RfuncLambda, 
     int* IfuncSNP)
{
    int P = dims[0];
    int nrow_IndexMatrix = dims[3], ncol_IndexMatrix = dims[4];
    int **IndexMatrix, nrow_tuningMatrix= dims[8];
    int ncol_tuningMatrix = dims[9], Q = dims[10], ncolBetaMatrix = dims[11];
    int NumTuning = dims2[0], nrow_funcLambda = dims2[1];
    int ncol_funcLambda = dims2[2], nrow_Zmatrix = dims2[3], ncol_Zmatrix = dims2[4];
    int nrow_AllTuningMatrix = dims2[6], ncol_AllTuningMatrix = dims2[7];
    
    double **BetaMatrix, **tuningMatrix, **summaryBetas, **SDvec;

    
    /* reorganize vector into matrix */
    reorg(RsummaryBetas, &summaryBetas, P, Q);
    reorg(RSDvec, &SDvec, P, Q);
    reorg(RBetaMatrix, &BetaMatrix, NumTuning, ncolBetaMatrix);
    
    
    reorg(RtuningMatrix, &tuningMatrix, nrow_tuningMatrix, ncol_tuningMatrix);
    reorg_int(RIndexMatrix, &IndexMatrix, nrow_IndexMatrix, ncol_IndexMatrix);
    
    double **Zmatrix, **AllTuningMatrix;
    int **funcLambda;
    
    /* reorganize vector into matrix */
    reorg(RAllTuningMatrix, &AllTuningMatrix, nrow_AllTuningMatrix, ncol_AllTuningMatrix);
    reorg(RZmatrix, &Zmatrix, nrow_Zmatrix, ncol_Zmatrix);
    reorg_int(RfuncLambda, &funcLambda, nrow_funcLambda, ncol_funcLambda);


    
    gsfFuncC(summaryBetas, ldJ, dims, Numitervec, IndexMatrix, IndJ,
     ldvec, ChrIndexBeta, RupperVal, Init_summaryBetas, SDvec, 
     tuningMatrix, BetaMatrix,penalty_, lambda0vec, Zmatrix, dims2, AllTuningMatrix, 
     lambdavec, funcLambda, IfuncSNP);
} 



void gsfFuncC(double** summaryBetas, int* ldJ, int* dims, int* Numitervec, 
     int** IndexMatrix, int* IndJ,
     double* ldvec, int*ChrIndexBeta, double*RupperVal, 
     double* Init_summaryBetas, double** SDvec, 
     double**tuningMatrix, double** BetaMatrix,
     char **penalty_,      double* lambda0vec, double** Zmatrix,
     int* dims2, double** AllTuningMatrix, double* lambdavec, int **funcLambda,
     int* IfuncSNP)
{
  
    int **skipb, niter, tuningIndex = -1, i, j,j1, j2, found, stuning = 0, cindex;
    int P = dims[0], NumIter = dims[1], breaking = dims[2], NumSNP = dims[3];
    int NumInd = dims[5], Zscale = dims[7];
    int nrow_tuningMatrix = dims[8];
    int Q = dims[10], k1, dfMax = dims[12];
    
    int nrow_funcLambda = dims2[1];
    int ncol_funcLambda = dims2[2];
    int leng_p_Threshold = dims2[5], df_q1, ncol_AllTuningMatrix = dims2[7];
    
    
    
    double tmp0, *tmpbetas, threshold, bj_bar, epsilon = 0.0001;
    double *sumBetas;
    double lambda1, lambda2, tau2, upperVal = RupperVal[0];
    
    
    double lambda0;
    double *lamtemp;
    int len_lamtemp = ncol_funcLambda, numfunc = ncol_funcLambda;
    lamtemp = (double *)calloc(len_lamtemp,sizeof(double));
    
  
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
    

    int kp, tui, jf, tui2, jf2; 
    
    for(kp=0; kp<leng_p_Threshold; kp++){
      
      lambda0 = lambda0vec[kp];
      for(tui=0; tui<nrow_funcLambda; tui++){
        for(tui2=0; tui2<nrow_tuningMatrix; tui2++){
          tuningIndex = tuningIndex + 1;
          Numitervec[tuningIndex] = 0;
          AllTuningMatrix[tuningIndex][0] = lambda0;
          
          for(jf=0; jf<numfunc; jf++){
            lamtemp[jf] = lambdavec[funcLambda[tui][jf]];
            AllTuningMatrix[tuningIndex][(jf+1)] = lamtemp[jf];
          }   
          lambda2 = tuningMatrix[tui2][2];
          tau2 = tuningMatrix[tui2][3];
          
          AllTuningMatrix[tuningIndex][(ncol_AllTuningMatrix - 2)] = lambda2;
          AllTuningMatrix[tuningIndex][(ncol_AllTuningMatrix - 1)] = tau2;
          

          /* Initiation */
      
          for(i=0; i<P; i++){  
            sumBetas[i] = 0.0;
            for(j=0; j<Q; j++){  
              jointBmatrix[i][j] = 0.0;
              tempBmatrix[i][j] = 0.0;
              skipb[i][j] = 0;
            } 
          }

          cindex = 0;      
          for(niter=0; niter<NumIter; niter++){
          
            if(NumInd!=0){
          
              for(j1=0; j1<NumInd; j1++){
                j = IndJ[j1];
                
                lambda1 = lambda0;
                
                if(IfuncSNP[j]==1){
                  for(jf2=0; jf2<numfunc; jf2++){
                    lambda1 = lambda1 + Zmatrix[j][jf2]*lamtemp[jf2];
                  }
                }
            
                for(k1=0; k1<Q; k1++){
                  bj_bar = summaryBetas[j][k1];
                  
                  if(bj_bar!=0.0){
                  
                    if(Q==1){
                      threshold = lambda1;
                    }else{
                      threshold = lambda1 + lambda2/(sumBetas[j] + tau2);
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
              
              lambda1 = lambda0;
                
              if(IfuncSNP[j]==1){
                for(jf2=0; jf2<numfunc; jf2++){
                  lambda1 = lambda1 + Zmatrix[j][jf2]*lamtemp[jf2];
                }
              }
              
              for(k1=0; k1<Q; k1++){
                if(skipb[j][k1]==0){
                  if(summaryBetas[j][k1]!=0.0){
                    tmp0 = 0.0;
                    for(j2=IndexMatrix[j1][1]; j2<(IndexMatrix[j1][2]+1); j2++){
                      tmp0 = tmp0 + ldvec[j2]*jointBmatrix[ldJ[j2]][k1];
                    }
                    bj_bar = (summaryBetas[j][k1] - tmp0);
            
                    if(Q==1){
                      threshold = lambda1;
                    }else{
                      threshold = lambda1 + lambda2/(sumBetas[j] + tau2);
                    }
                
                    //Rprintf("threshold=%e\n",threshold);
              
                    if(fabs(bj_bar)>upperVal){
                      if(breaking==1){
                        Numitervec[tuningIndex] = -1;
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
               Numitervec[tuningIndex] = -2;
               break;
            }


            if(found==0){            
              k1=0;
              for(i=0; i<Q; i++){
                for(j=0; j<P; j++){
                  BetaMatrix[tuningIndex][k1] = jointBmatrix[j][i];
                  k1 = k1 + 1;
                } 
              } 
            
              Numitervec[tuningIndex] = (niter+1);
              stuning = tuningIndex;
              break;
            }
            
            for(j1=0; j1<NumSNP; j1++){
              j = IndexMatrix[j1][0];
              tmp0 = 0.0;
              
              for(k1=0; k1<Q; k1++){
                tempBmatrix[j][k1] = jointBmatrix[j][k1];
                tmp0 = tmp0 + fabs(jointBmatrix[j][k1]);
              }
              sumBetas[j] = tmp0;
            }

            if(NumInd!=0){
              for(j1=0; j1<NumInd; j1++){
                j = IndJ[j1];
                tmp0 = 0.0;
                for(k1=0; k1<Q; k1++){
                  tempBmatrix[j][k1] = jointBmatrix[j][k1];
                  tmp0 = tmp0 + fabs(jointBmatrix[j][k1]);
                }
                sumBetas[j] = tmp0;
              }
            }
            Numitervec[tuningIndex] = (niter+1);
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


 
