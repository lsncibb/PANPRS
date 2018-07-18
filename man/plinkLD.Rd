\name{plinkLD}
\alias{plinkLD}
\title{The LD info from output of the software (plink)}
\description{ The LD information is crucial for the analysis by SummaryLasso. 
	The reference alleles used to obtained for the Z statsitics or the 
	regression coefficients have to be the sames as those used for the LD calculation. 
	This file can be obtained directly from the output of the LD calculation by the 
	software (plink); for example the output can be like plink.ld.
    On the other hand, the user can calcuate the LD based on their prefered tools.
    The variables are as
  follows: 

 \itemize{
   \item CHR_A: The chromosome of SNP_A 
   \item BP_A: The positions of SNP_A  
   \item SNP_A: The names of SNP_A 
   \item CHR_B: The chromosome of SNP_B 
   \item BP_B: The positions of SNP_B 
   \item SNP_B: The names of SNP_B    
   \item R: The correlation between SNP_A and SNP_B
 }}
\usage{data(plinkLD)}
\format{A data frame with 205959 rows and 7 columns}
\references{
  \itemize{
    \item Purcell S, et al. (2007) 
        PLINK: a toolset for whole-genome association and population-based 
	linkage analysis. \emph{American Journal of Human Genetics}, \strong{81}.
  }
}
\keyword{datasets}
