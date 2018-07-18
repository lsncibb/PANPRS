`gsPEN` <- function(summaryZ, Nvec, plinkLD, NumIter = 100, breaking = 1, 
  numChrs=22, ChrIndexBeta = 0, Init_summaryBetas= 0, Zscale = 1, RupperVal = NULL, 
  tuningMatrix = NULL, penalty=c("mixLOG"), OutputForPath=0,
  taufactor = c(1/25, 1, 10), llim_length = 10, subtuning = 50, Lambda_limit = c(0.5,0.9), 
  Lenlam_singleTrait = 200, dfMax = NULL, outputAll = 0 ){
 
 
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
    

  if(is.null(tuningMatrix)){
    medianval = median(apply(abs(summaryBetas),1,sum),na.rm=T)
    tauvec = sort(medianval*taufactor)
    
    Lambda_limit = quantile(abs(summaryZ[,1]), Lambda_limit)

    tuningMatrix = Tuning_setup_group_only(tauvec, subtuning, Lambda_limit, Lenlam_singleTrait, llim_length, medianval)
  }            
     

  inv_summaryBetas = 0
  count_nonzero = function(xx){
    length(which(xx!=0))
  }
  counts = apply(summaryBetas,1,count_nonzero)


  
  Betaindex = c(1:nrow(summaryBetas))-1
  SNPnames = rownames(summaryBetas)
  
  ldJ = PlinkLD_transform(plinkLD, SNPnames)
  rm(plinkLD)

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
  
  if(nrow(tuningMatrix)>1){
    NumTuning = nrow(tuningMatrix)
  }

  ncolBetaMatrix = P*Q
  dims = rep(0,14)
  dims[1] = NumTuning
  dims[2] = P
  dims[3] = NumIter
  dims[4] = breaking
  dims[5] = nrow_IndexMatrix
  dims[6] = ncol_IndexMatrix
  dims[7] = length(wind)  
  dims[8] = Zscale
  dims[9] = nrow(tuningMatrix)
  dims[10] = ncol(tuningMatrix)
  dims[11] = Q
  dims[12] = ncolBetaMatrix
  dims[13] = OutputForPath 
  dims[14] = dfMax 
  
  #dims[14] = length_ldJ
  #dims[15] = numChrs


  
  Numitervec = rep(0, NumTuning)
  BetaMatrix = matrix(0, NumTuning, ncolBetaMatrix)
  

  Z = .C("gsPEN", as.double(t(summaryBetas)),as.integer(ldJ), as.integer(dims), Numitervec=as.integer(Numitervec),
   as.integer(t(IndexMatrix)), as.integer(IndJ),
   as.double(ldvec), as.double(inv_summaryBetas), as.integer(ChrIndexBeta), as.double(RupperVal),
   as.double(Init_summaryBetas), as.double(t(SDvec)),
   as.double(t(tuningMatrix)), BetaMatrix = as.double(t(BetaMatrix)),
   penalty, PACKAGE="SummaryLasso")

  BetaMatrix = matrix(Z$BetaMatrix, nrow = NumTuning, ncol = ncolBetaMatrix, byrow = TRUE)
   
  colnames(tuningMatrix) = c("lam1","null","lam2","tau")
  tuningMatrix = tuningMatrix[,c(1,3,4)]
  tuningMatrix[c(1:Lenlam_singleTrait), c(2:3)] = NA
   
  colnames(BetaMatrix) = paste0(rep(SNPnames, times = Q),".trait", rep(c(1:Q), each = P))
   
  Numitervec = Z$Numitervec
  output = Cleaning(BetaMatrix=BetaMatrix, Numitervec=Numitervec, AllTuningMatrix=tuningMatrix)
  Numitervec = output[[1]]
  BetaMatrix = output[[2]]
  tuningMatrix = output[[3]]
  rm(output)
  if(outputAll==0){
    convergeIndex = which(Numitervec > 0)
    Numitervec = Numitervec[convergeIndex]
    BetaMatrix = BetaMatrix[convergeIndex,]
    tuningMatrix = tuningMatrix[convergeIndex,]
  }
   
  ll = list(BetaMatrix = BetaMatrix, Numitervec = Numitervec, tuningMatrix = tuningMatrix)
  return(ll)
   
}

