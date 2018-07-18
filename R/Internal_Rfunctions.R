
Tuning_setup_group_only = function(tauvec, subtuning, Lambda_limit, Lenlam, llim_length, medianval){
   
   lambdavec = seq((min(Lambda_limit)),(max(Lambda_limit)),len=Lenlam)
   tuningMatrix = cbind(lambdavec,1,0,1)
   
   topvec = seq(from = min(Lambda_limit), to = max(Lambda_limit),len=llim_length)
   
   for(i1 in 1:length(tauvec)){
     tau = tauvec[i1]
     
     for(t1 in 2:length(topvec)){
       top = topvec[t1]        
       lambdavec = seq(min(Lambda_limit),(top-0.05),len=subtuning)

       temp = (top - lambdavec)*(medianval + tau)
       tuningMatrix = rbind(tuningMatrix,cbind(lambdavec,tau,temp,tau))
     }
   }
   return(tuningMatrix)
}




Tuning_setup_group_func = function(lambdavec_func, lambdavec_func_limit_len, p.Threshold, numfunc,  tauvec, subtuning, Lambda_limit, Lenlam, llim_length, medianval){

  lambdavec = seq((min(Lambda_limit)),(max(Lambda_limit)),len=Lenlam)
   
  topvec = seq(from = min(Lambda_limit), to = max(Lambda_limit),len=llim_length)
  
  Start = 1  
  for(i1 in 1:length(tauvec)){
    tau = tauvec[i1]
     
    for(t1 in 2:length(topvec)){
      top = topvec[t1]        
      lambdavec = seq(min(Lambda_limit),(top-0.05),len=subtuning)

      temp = (top - lambdavec)*(medianval + tau)
      if(Start==1){
        tuningMatrix = cbind(lambdavec,tau,temp,tau)
        Start = 0
      }else{
        tuningMatrix = rbind(tuningMatrix,cbind(lambdavec,tau,temp,tau))
      }
    }
  }
  if(is.null(lambdavec_func)){
    lambdavec_func = seq(0, lambdavec_func_limit_len[1], length.out=lambdavec_func_limit_len[2])
  }
  n.threshold = length(p.Threshold)
  leng.lambda = length(lambdavec_func)
  funcLambda0 = permutations(length(lambdavec_func), numfunc, repeats.allowed=T)
  funcLambda = funcLambda0 - 1
  
  
  
  output = list()
  output[[1]] = funcLambda
  output[[2]] = lambdavec_func
  output[[3]] = tuningMatrix

  return(output)
}

nonzero = function(xx){
  return(length(which(xx!=0)))
}

Cleaning = function(BetaMatrix, Numitervec, AllTuningMatrix){
   tuningvec = apply(AllTuningMatrix,1,paste0,collapse=":")
   uniqvec = unique(tuningvec)
   mat = match(uniqvec,tuningvec)
   Numitervec = Numitervec[mat]
   BetaMatrix = BetaMatrix[mat,]
   AllTuningMatrix = AllTuningMatrix[mat,]
   NumCounts = apply(BetaMatrix,1,nonzero)
   od = order(NumCounts)
   Numitervec = Numitervec[od]
   BetaMatrix = BetaMatrix[od,]
   AllTuningMatrix = AllTuningMatrix[od,]
   output = list()
   output[[1]] = Numitervec
   output[[2]] = BetaMatrix
   output[[3]] = AllTuningMatrix
   return(output)
}


