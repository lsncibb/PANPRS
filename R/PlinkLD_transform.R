PlinkLD_transform <- function(plinkLD, Allkeepsnps){

  wkeep = which( plinkLD$SNP_A  %in% Allkeepsnps & plinkLD$SNP_B %in% Allkeepsnps)
  plinkLD = plinkLD[wkeep,]

  if(length(which(is.na(plinkLD)))>0){stop("")}
	
  ldJ = plinkLD[,c("SNP_B","SNP_A","R")]
  names(ldJ) = c("SNP_A","SNP_B","R")
	  
  ldJ = rbind(plinkLD[,c("SNP_A","SNP_B","R")],ldJ)

  return(ldJ)
}
