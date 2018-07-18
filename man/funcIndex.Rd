\name{funcIndex}
\alias{funcIndex}
\title{Inputs for the functional annotations of SNPs.}
\description{A 3614 x 3 matrix with (0,1) entry with 3614 SNPs and 
3 functional annotations. For the element at i-th row, j-th column, the entry 0 means SNP i without j-th functional annotation; entry 1 means otherwise.
  follows: 

 \itemize{
   \item f1: The binary index for functional annotation 1.
   \item f2: The binary index for functional annotation 2.
   \item f3: The binary index for functional annotation 3.
 }}
\usage{data(summaryZ)}
\format{A matrix with 3614 rows for the 3614 SNPs and 3 columns for functional annotations.}
\keyword{datasets}
