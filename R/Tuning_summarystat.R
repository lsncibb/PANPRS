Tuning_summarystat <-
function(betavec,family="gaussian",Penalty="LOG", n.tau=6, n.lambdas=100, lambda.min=NULL, pfactor=0.1, minTomax = FALSE){
  
  output = list()
  output[["lambda"]] = NULL
  output[["tau"]] = NULL

  if (family=="gaussian"){
    
    if(Penalty == "LOG"){
    
   	  l1.max = max(abs(betavec))
      if(is.null(lambda.min)){
        minLambda = l1.max*pfactor}else{
        minLambda = lambda.min
        }
	  maxTau = max(abs(betavec))
      if(minTomax){
        Thresholds <- c(exp(seq(log(minLambda),log(l1.max),len=n.lambdas)))
      }else{
        Thresholds <- c(exp(seq(log(l1.max),log(minLambda),len=n.lambdas)))
      }
      tauset = c(exp(seq(log(1e-6),log(maxTau),len=n.tau)))

      tau = lambda = c()

      for(t in 1:length(Thresholds)){
        
        thres = Thresholds[t]  
        slambda = tauset*thres
  
        tau = c(tau, tauset)
        lambda = c(lambda, slambda)
      }

      if(length(tau)!=length(lambda)){
        stop("error")}
      
      output[["lambda"]] = lambda
      output[["tau"]] = tau
    }
  }
  return(output)
}
