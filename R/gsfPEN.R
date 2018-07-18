`gsfPEN` <- function(summaryZ, Nvec, plinkLD, NumIter = 1000, RupperVal = NULL, breaking = 1, 
 numChrs=22, ChrIndexBeta = 0, Init_summaryBetas= 0, Zscale = 1, 
 tuningMatrix = NULL, penalty="mixLOG",
 funcIndex, numfunc, p.Threshold = NULL, p.Thresholdpara = c(0.5, 10^-4, 4),
 taufactor = c(1/25, 1, 3), llim_length = 4, subtuning = 4, Lambda_limit = c(0.5,0.9), Lenlam = 4,
 lambdavec_func = NULL, lambdavec_func_limit_len = c(1.5, 3), dfMax = NULL, outputAll = 0){
 

 
  if(Zscale!=1){error("Tuning values set-up for multiple traits analysis requires Zscale=1.")}
  
  
  Nq = length(Nvec) 
  

  summaryBetas = matrix(0, nrow(summaryZ), Nq)
  SDvec = matrix(0, nrow(summaryZ), Nq)
    
  for(ii in 1:Nq){
    summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
    SDvec[,ii] = 1/sqrt(Nvec[ii])
  }
    
  rownames(summaryBetas) = rownames(summaryZ)
  
  if(is.null(dfMax)){
    dfMax = ceiling(0.7*nrow(summaryZ))
  }

  
  if(is.null(p.Threshold)){
    p.Threshold = seq(p.Thresholdpara[1], p.Thresholdpara[2], length.out=p.Thresholdpara[3])
  }

  
  if(any(c(is.null(tuningMatrix), is.null(lambdavec_func)))){
    medianval = median(apply(abs(summaryBetas),1,sum),na.rm=T)
    tauvec = sort(medianval*taufactor)
    Lambda_limit = quantile(abs(summaryZ[,1]), Lambda_limit)
    
    output = Tuning_setup_group_func(lambdavec_func, lambdavec_func_limit_len, p.Threshold, numfunc,  tauvec, subtuning, Lambda_limit,
      Lenlam, llim_length, medianval)
      
    funcLambda = output[[1]]
    lambdavec_func = output[[2]]
    tuningMatrix = output[[3]]
  }else{
    funcLambda0 = permutations(length(lambdavec_func), numfunc, repeats.allowed=T)
    funcLambda = funcLambda0 - 1
  }
  

  
  Betaindex = c(1:nrow(summaryBetas))-1
  SNPnames = rownames(summaryBetas)
  
  ldJ = PlinkLD_transform(plinkLD, SNPnames)


  JidMatrix = matrix(,nrow(ldJ),2)
  mat1 = match(ldJ[,1],SNPnames)
  mat2 = match(ldJ[,2],SNPnames)

  JidMatrix[,1] = Betaindex[mat1]
  JidMatrix[,2] = Betaindex[mat2]

  ldJ[,1] = JidMatrix[,1]
  ldJ[,2] = JidMatrix[,2]

  od = order(JidMatrix[,1],JidMatrix[,2], decreasing=F)
  ldJ = ldJ[od,]
  
  wind = which(! Betaindex %in% ldJ[,1])

  IndJ = -1

  if(length(wind) > 0){
    IndJ = Betaindex[wind]
  }

  Counts = table(ldJ[,1])
  NumSNP = length(Counts)

  
  IndexS = c(0,cumsum(Counts)[-NumSNP])
  IndexE = cumsum(Counts)-1

  IndexMatrix = matrix(,NumSNP,3)

  IndexMatrix[,1] = as.numeric(names(Counts))

  nrow_IndexMatrix = NumSNP
  ncol_IndexMatrix = ncol(IndexMatrix)


  IndexMatrix[,2] = IndexS
  IndexMatrix[,3] = IndexE


  ldvec = ldJ[,3]

  ldJ = ldJ[,2]

  length_ldJ = length(ldJ)

  if(is.null(RupperVal)){
    RupperVal = ceiling(max(abs(summaryBetas),na.rm=T)*50)
  }


  P = nrow(summaryBetas)
  Q = ncol(summaryBetas)
  

 
  if(nrow(funcIndex)!=nrow(summaryBetas)){
    stop("nrow of summaryBetas and row of funcIndex do not match.")
  }
    
  Zmatrix = 1 - funcIndex
    
  sumfuncIndex = apply(Zmatrix,1,sum)
  IfuncSNP = rep(0, P)
  IfuncSNP[which(sumfuncIndex!=0)] = 1

  
  nrow_funcLambda = nrow(funcLambda)
  ncol_funcLambda = ncol(funcLambda)
  leng.p.Threshold = length(p.Threshold)
  
  nrow_Zmatrix = nrow(Zmatrix)
  ncol_Zmatrix = ncol(Zmatrix)
    
  NumTuning = nrow(tuningMatrix)*leng.p.Threshold*nrow_funcLambda
  print(paste0("Number of total tuning combinations = ",NumTuning))
    
  
  nrow_AllTuningMatrix = NumTuning
  ncol_AllTuningMatrix = (numfunc+1) + 2
  AllTuningMatrix = matrix(0, nrow_AllTuningMatrix, ncol_AllTuningMatrix)
  
  dims2 = rep(0,8)
  dims2[1] = NumTuning
  dims2[2] = nrow_funcLambda
  dims2[3] = ncol_funcLambda
  dims2[4] = nrow_Zmatrix 
  dims2[5] = ncol_Zmatrix
  dims2[6] = leng.p.Threshold
  dims2[7] = nrow_AllTuningMatrix
  dims2[8] = ncol_AllTuningMatrix
  
  
  lambda0vec = abs(-qnorm(p.Threshold/2))

  nrow_BetaMatrix = NumTuning
  ncolBetaMatrix = P*Q
  dims = rep(0,13)
  dims[1] = P
  dims[2] = NumIter
  dims[3] = breaking
  dims[4] = nrow_IndexMatrix
  dims[5] = ncol_IndexMatrix
  dims[6] = length(wind)
  dims[7] = numChrs
  dims[8] = Zscale 
  dims[9] = nrow(tuningMatrix)
  dims[10] = ncol(tuningMatrix)
  dims[11] = Q  
  dims[12] = ncolBetaMatrix
  dims[13] = dfMax
  
  Numitervec = rep(0, NumTuning)
  BetaMatrix = matrix(0, nrow_BetaMatrix, ncolBetaMatrix)
  

  Z = .C("gsfFunc", as.double(t(summaryBetas)),as.integer(ldJ), as.integer(dims), Numitervec=as.integer(Numitervec),
    as.integer(t(IndexMatrix)), as.integer(IndJ),
    as.double(ldvec), as.integer(ChrIndexBeta), as.double(RupperVal),
    as.double(Init_summaryBetas), as.double(t(SDvec)), 
    as.double(t(tuningMatrix)), BetaMatrix = as.double(t(BetaMatrix)),
    penalty, as.double(lambda0vec), as.double(t(Zmatrix)),
    as.integer(dims2), AllTuningMatrix = as.double(t(AllTuningMatrix)), as.double(lambdavec_func),
    as.integer(t(funcLambda)), as.integer(IfuncSNP),
    PACKAGE="SummaryLasso")

  BetaMatrix = matrix(Z$BetaMatrix, nrow = NumTuning, ncol=ncolBetaMatrix, byrow=TRUE)
  colnames(BetaMatrix) = paste0(rep(SNPnames, times = Q),".trait", rep(c(1:Q), each = P))

  AllTuningMatrix = matrix(Z$AllTuningMatrix, nrow = nrow_AllTuningMatrix, ncol = ncol_AllTuningMatrix, byrow = TRUE)
  Numitervec = Z$Numitervec
  
  colnames(AllTuningMatrix) = c("lam0",paste0("lam.f",c(1:numfunc)),"lam2","tau")
  
  output = Cleaning(BetaMatrix=BetaMatrix, Numitervec=Numitervec, AllTuningMatrix=AllTuningMatrix)
  Numitervec = output[[1]]
  BetaMatrix = output[[2]]
  AllTuningMatrix = output[[3]]
  rm(output)
  if(outputAll==0){
    convergeIndex = which(Numitervec > 0)
    Numitervec = Numitervec[convergeIndex]
    BetaMatrix = BetaMatrix[convergeIndex,]
    AllTuningMatrix = AllTuningMatrix[convergeIndex,]
  }
       
  ll = list(BetaMatrix = BetaMatrix, Numitervec = Numitervec, AllTuningMatrix = AllTuningMatrix)
  
  return(ll)
   
}

